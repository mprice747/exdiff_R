# Given a vector of alpha parameters, sample one point from a dirichlet distribution
rdirichlet <- function(alpha_vec) {
  k <- length(alpha_vec)
  gamma_samp <- rgamma(k, alpha_vec)
  return(gamma_samp/sum(gamma_samp))
}

# Sample from a multivariate normal where ||x[1:dim]|| < pi. Uses accept/reject strategy.
# We call this a circular mulivariate normal distribution
rcnorm_accept_reject <- function(n, mu, Sigma, r, dim = length(mu)) {
  
  # sampled_mat will contain samples
  num_accepted <- 0
  p <- nrow(Sigma)
  sampled_mat <- matrix(rep(0, p * n), nrow = n)
  
  # Sample from multivariate normal and see if norm is less than r. If yes, add to 
  # sample_mat
  while (num_accepted < n) {
    
    one_sample <- rmvnorm(1, mu, Sigma)
    
    if (sum(one_sample[1:dim]^2) < r^2) {
      num_accepted <- num_accepted + 1
      sampled_mat[num_accepted, ] <- one_sample
    }
    
  }
  
  return(sampled_mat)
}



# Obtain the normalizing constant (1/ P(||x[1:dim]|| < r)) for a circular normal distribution
# Needed for log likelihood calculations 
get_normalize_const_cnorm <- function(mu, Sigma, r, dim) {
  
  eigen_decomp <- eigen(Sigma[1:dim, 1:dim])
  
  B_matrix <- eigen_decomp$vectors %*% sqrt(diag(eigen_decomp$values))
  
  a_vector <- solve(B_matrix, mu[1:dim])
  
  prob <- 1 - imhof(r^2, eigen_decomp$values, h = rep(1, dim),
                    delta = a_vector^2,
                    epsabs = 10^(-6), epsrel = 10^(-6), limit = 10000)$Qq
  
  return(1/prob)
}

# Sample from a mixed multivariate normal distribution where ||x[1:dim]|| < pi. 
# Uses accept/reject strategy. Call this a mixture of circular normal distributions
mixed_rcnorm_accept_reject <- function(mean_cov_list, probs, 
                                       norm_consts, r, dim) {
  
  # Sample a component based on probabilities 
  comp_samp <- sample.int(length(probs), 1,
                          prob = probs)
  
  # Get sampled component's mean, covariance matrix and normalizing constant 
  mu_value <- mean_cov_list[[comp_samp]]$mean
  Sigma_value <- mean_cov_list[[comp_samp]]$cov
  
  # Normalizing constant will be number of points to sample at a time 
  # (Serves as expected number of samples until one is within the radius)
  norm_const_value <- max(2, round(norm_consts[comp_samp]))
  
  found_result <- FALSE
  
  while(!found_result){
    
    # Sample and see if they are within boundaries
    pot_samples <- rmvnorm(norm_const_value, mu_value, Sigma_value)
    
    in_circle <- which(rowSums(pot_samples[, 1:dim]^2) < r^2)
    
    # If yes, return result. If not, try again
    if (length(in_circle) > 0){
      final_sample <- pot_samples[in_circle[1], ]
      found_result <- TRUE
    }
  }
  return(final_sample)
}

# Sample n points from a mixture of mixture of circular normal distributions. Uses 
# parallelism to speed up computations
r_mixed_cnorm <- function(n, mean_cov_list, probs, 
                          norm_consts, r, dim) {
  
  p <- length(mean_cov_list[[1]]$mean)
  
  # Intialize progress bar
  print("Importance Sampling Progress Bar:")
  
  cl <- makeCluster(detectCores())
  registerDoSNOW(cl)
  funcs <- ls(envir = .GlobalEnv)
  
  pb <- txtProgressBar(max = n, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Sample and use parallelism to speed up computations 
  mixed_cnorm_samples <- foreach(i = 1:n, .combine = rbind,
                                 .packages = c('stats', 'MASS', 'mvtnorm' ), 
                                 .export = funcs, .options.snow = opts) %dopar% {
                                   
                                   
             mixed_cnorm_sample <- mixed_rcnorm_accept_reject(mean_cov_list, 
                                                              probs, 
                                                              norm_consts, r, dim)
             
             mixed_cnorm_sample
                                 }
  
  close(pb)
  stopCluster(cl)
  
  # Return result
  row.names(mixed_cnorm_samples) <- seq(1, n)
  
  return(mixed_cnorm_samples)
}


# Using WLB Samples and the prior distribution, create a mixture proposal density
# We fit a mixed normal with the Weighted Likelihood Bootstrap samples and then add the prior
# distribution as a component
get_proposal_cluster <- function(wlb_list, prior_mean, prior_sd, 
                                 prop_mix, prop_prior, num_betas) {
  
  
  normalize_props <- c(prop_mix, prop_prior)/sum(c(prop_mix, prop_prior))
  
  wlb_samples <- do.call(rbind, 
                         lapply(wlb_list$real_wlb_samples, function(x){x$min_param}))
  
  # Use Mclust to fit a mixed normal 
  mclust_wlb <- Mclust(wlb_samples, verbose = FALSE)
  num_comps <- length(mclust_wlb$parameters$pro)

  # Get the fitted probabilities, means and covariances for each compoent 
  mclust_probs <- c(normalize_props[1] * mclust_wlb$parameters$pro, 
                    normalize_props[2])
  mclust_means <- unname(cbind(mclust_wlb$parameters$mean, prior_mean))
  orig_covs <- mclust_wlb$parameters$variance$sigma
  mclust_covs <- abind(orig_covs, diag(prior_sd^2), along = 3)
  

  # Proposal density will be contained in a list of lists, with each list 
  # containing a vector representing the mean of the component and a matrix 
  # containing the covariance of the component 
  mixed_normal_prop <- list()
  norm_consts <- rep(0, length(mclust_probs))
  
  for(i in 1:length(mclust_probs)){
    
    # Add mean and covariance to a list
    mixed_comp <- list()
    mixed_comp$mean <- mclust_means[, i]
    mixed_comp$cov <- mclust_covs[, , i]
    mixed_normal_prop[[i]] <- mixed_comp
    
    norm_consts[i] <- get_normalize_const_cnorm(mixed_comp$mean,  mixed_comp$cov, 
                                             pi, num_betas)
  }
  
  # Return list of means and covairance, covariance probs (weights) and normalizing
  # constant because ||betas|| < pi
  return(list(mixed_normal_prop = mixed_normal_prop,
              weights = mclust_probs,
              norm_consts = norm_consts))
}

# Obtain a WLB sample by optimizing a weighted log likelihood. Because of the non-concavity
# of the log likelihood, we use differential evolution (a genetic approach) to find global 
# optima

# Inputs:
# input_X - length n vector, transformed X points to be on 0-1 scale
# input_Y - length n vector, Y points
# num_betas - integer, number of parameters for the diffeomorphism 
# num_stationary - integer, number of stationary points 
# first direction - 1 or -1, 1 if the difference between the function value the first 
# stationary point and f(0) is positve, -1 if it is negative
# zero_is_zero - boolean, TRUE if f(0) = 0, FALSE if f(0) needs to be estimated
# de_CR - number in [0, 1], crossover probability for differential evolution. Defaults to 0.75.
# de_F - number in [0, 2], differential weighting factor for differential evolution. Defaults to 0.8.
# de_itermax - positive integer, maximum iterations for differential evolution. Defaults to 200.

# Ouputs: 
# min_param - WLB sample, global optimum of weighted log likelihood found with differntial evolution
# min_cost - cost value of WLB sample

get_wlb_sample <- function(input_X, input_Y, num_betas, num_stationary,
                           first_direction, zero_is_zero,
                           de_CR = 0.75, de_F = 0.8, de_itermax = 200){
  
  # Simulate weight vector from rdirichlet
  weight_vector <- rdirichlet(rep(1, length(input_X)))
  
  
  num_changes <- num_stationary + 1
  max_change <- max(abs(input_Y)) - min(min(input_Y), 0)
  
  # Number of population members in differential evolution
  NP_param <- 12 * length(prior_mean)
  
  # Define search bounds (min and max) for differential evolution
  if (zero_is_zero){
    
    min_bounds <- c(rep(-pi, num_betas), 
                    rep(log(min(abs(diff(input_Y[order(input_X)])))/2), num_changes), 
                    log(1e-6))
    
    max_bounds <- c(rep(pi, num_betas), rep(log(max_change), num_changes + 1))
    
  } else {
    
    min_bounds <- c(rep(-pi, num_betas), min(input_Y) - sd(input_Y), 
                    rep(log(min(abs(diff(input_Y[order(input_X)])))/2), num_changes), 
                    log(1e-6))
    
    max_bounds <- c(rep(pi, num_betas), max(input_Y) + sd(input_Y), 
                    rep(log(max_change), num_changes + 1))
  }
  
  
  # Optimize using differential evolution
  first_generation <- rcnorm_accept_reject(NP_param, prior_mean, diag(prior_sd^2), 
                                           pi,num_betas)
  
  de_res <- DEoptim(neg_bayes_value, min_bounds, max_bounds, 
                    DEoptim.control(NP = NP_param, itermax = de_itermax,
                                    CR = de_CR, 
                                    F = de_F, trace = FALSE, 
                                    initialpop = first_generation), 
                    num_betas = num_betas,
                    input_X = input_X, input_Y = input_Y,
                    first_direction = first_direction,
                    prior_mean = NULL,
                    prior_sd = NULL, 
                    zero_is_zero = zero_is_zero,
                    weight_vector = weight_vector)
  
  
  final_res <- tryCatch({
    
    # Take best result from differential evolutiona and use as starting point for
    # BFGS. Theoretically should be global optima 
    bfgs_res <- optim(as.vector(de_res$optim$bestmem), 
                      fn = neg_bayes_value, 
                      method = 'BFGS', control = list(maxit = 10000),
                      num_betas = num_betas,
                      input_X = input_X, input_Y = input_Y,
                      first_direction = first_direction,
                      prior_mean = NULL,
                      prior_sd = NULL, 
                      zero_is_zero = zero_is_zero,
                      weight_vector = weight_vector, 
                      hessian = TRUE)
    
    final_res <- list(min_param = bfgs_res$par,
                      min_cost = bfgs_res$value)
    
    final_res
    
  }, error = function(e){
    
    # If BFGS produced error, just take the differential evolution result
    min_param <- as.vector(de_res$optim$bestmem)
    min_cost <- de_res$optim$bestval
    
    final_res <- list(min_param = min_param, 
                      min_cost = min_cost)
    final_res
    
  })
  
  return(final_res)
}

# Obtain weighted likelihood bootstrap samples by optimizinf a weighted log likelihood with
# differential evolution

# Inputs:
# input_X - length n vector, transformed X points to be on 0-1 scale
# input_Y - length n vector, Y points
# num_betas - integer, number of parameters for the diffeomorphism 
# num_stationary - integer, number of stationary points 
# first direction - 1 or -1, 1 if the difference between the function value the first 
# prior_mean - p dimensional vector of prior mean (Joint Independent Normal)
# prior_sd - p dimension vector of prior sds (Joint Independent Normal)
# zero_is_zero - boolean, TRUE if f(0) = 0, FALSE if f(0) needs to be estimated
# num_wlb_samples - positive integer, number of WLB samples to obtain
# de_CR - number in [0, 1], crossover probability for differential evolution. Defaults to 0.75.
# de_F - number in [0, 2], differential weighting factor for differential evolution. Defaults to 0.8.
# de_itermax - positive integer, maximum iterations for differential evolution. Defaults to 200.

# Outputs: 
# real_wlb_samples - list containing the WLB samples and its optimization cost
# wlb_weights - weight given to each sample proportional to its posterior value

get_wlb_samples <- function(input_X, input_Y, num_betas, num_stationary,
                            first_direction, prior_mean, prior_sd, zero_is_zero, 
                            num_wlb_samples,  de_CR = 0.75,
                            de_F = 0.8, de_itermax = 200) {
  
  # Make a progress bar for the Weighted Likelihood Bootstrap 
  print("WLB Progress Bar:")
  cl <- makeCluster(detectCores())
  registerDoSNOW(cl)
  funcs <- ls(envir = .GlobalEnv)
  
  pb <- txtProgressBar(max = num_wlb_samples, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Obtain WLB samples with parallelization 
  real_wlb_samples <- foreach(i = 1:num_wlb_samples, .packages = c('pracma', 'stats', 'MASS', 
                                                              'mvtnorm', 'DEoptim'), 
                         .export = funcs, .options.snow = opts) %dopar% {
                           
                           
                           wlb_result <- get_wlb_sample(input_X, input_Y, num_betas, num_stationary,
                                                        first_direction, zero_is_zero, 
                                                        de_CR, de_F, de_itermax)
                           
                           wlb_result
                           
                           
                         }
  
  close(pb)
  stopCluster(cl)
  
  # Transform 
  num_wlbs <- length(real_wlb_samples)
  
  sampled_params <-matrix(rep(0, num_wlbs * length(real_wlb_samples[[1]]$min_param)), 
                          nrow = num_wlbs)
  
  for (i in 1:num_wlbs){
    sampled_params[i, ] <- real_wlb_samples[[i]]$min_param
  }
  
  # Get the unnormalized posterior value for each of the WLB samples
  posterior_likelihood <- apply(sampled_params, 1, bayes_value, num_betas = num_betas,
                                input_X = input_X, input_Y = input_Y, first_direction = first_direction,
                                prior_mean = prior_mean, prior_sd = prior_sd, 
                                zero_is_zero = zero_is_zero,
                                weight_vector = NULL)
  
  # Normalize the weighted likelihood bootstrap
  wlb_weights <- softmax(posterior_likelihood)
  
  return(list(real_wlb_samples = real_wlb_samples, 
              wlb_weights = wlb_weights))
}

# Obtain posterior samples via importance sampling. The proposal distribution is created with 
# the WLB samples and consists of three components: 
# 1) A categorical distribution with the WLB samples, with weights proportional the samples' posterior values
# 2) A mixed normal fitted with the WLB samples
# 3) The prior distribution

# Inputs:
# n - positive integer, number of importance samples
# wlb_list - list, output from get_wlb_samples
# prior_mean - p dimensional vector of prior mean (Joint Independent Normal)
# prior_sd - p dimension vector of prior sds (Joint Independent Normal)
# prop_imp_sampler - 3 dimensional vector adding to 1, weights to be given to each of the 
# three components in the proposal density
# num_betas - integer, number of parameters for the diffeomorphism 

# Outputs: 
# importance_samples - n x p matrix containing importance samples
# likelihood - n dimensional vector contaning log likelihood of proposal distribution

wlb_importance_sampler <- function(n, wlb_list, prior_mean, prior_sd, 
                                   prop_imp_sampler, num_betas) {

  # Get mixed normal proposal density (includes prior distribution)
  prop_dist <- get_proposal_cluster(wlb_list, prior_mean, prior_sd, 
                                      prop_imp_sampler[2], prop_imp_sampler[3], 
                                      num_betas)
    
  wlb_samples <- do.call(rbind, 
                         lapply(wlb_list$real_wlb_samples, function(x){x$min_param}))
  
  
  # Determine how many samples from WLB samples and mixture density
  wlb_prop <- prop_imp_sampler[1]
  component_sample <- sample.int(2, n, replace = TRUE,  
                                 prob = c(wlb_prop, 1 - wlb_prop))
  
  num_from_wlb <- sum(component_sample == 1)
  num_from_mix <- n - num_from_wlb
  
  # Resample WLB samples
  resampled_wlb <- sample.int(nrow(wlb_samples), 
                              num_from_wlb, replace = TRUE, 
                              prob = wlb_list$wlb_weights)
  wlb_imp_samples <- wlb_samples[resampled_wlb, ]
  
  # Sample from Mixture density and combine results
  mixed_cnorm_samples <- r_mixed_cnorm(num_from_mix, prop_dist$mixed_normal_prop, 
                                       prop_dist$weights, 
                                       prop_dist$norm_consts, pi, num_betas)
  
  all_imp_samps <- rbind(wlb_imp_samples, mixed_cnorm_samples)
  
  num_comps <- length(prop_dist$weights)
  
  # Now will et log likelihood of the proposal distribution
  log_like_mat <- matrix(rep(0, n * (num_comps + 1)), nrow = n)
  
  if (num_from_wlb > 0){
    
  # First add categorical probabilities for the WLB Samples
  log_like_mat[1:num_from_wlb, 1] <- log(wlb_list$wlb_weights * wlb_prop)[resampled_wlb]}
  
  # For each component, get log likelihood for each important sampler
  for (comp in 1:num_comps) {
    
    mu_value <- prop_dist$mixed_normal_prop[[comp]]$mean
    Sigma_value <- prop_dist$mixed_normal_prop[[comp]]$cov
    
    dnorm_value <- dmvnorm(all_imp_samps, mu_value, Sigma_value,
                           log = TRUE)
    
    log_like_mat[, comp + 1] <- dnorm_value + log(prop_dist$weights[comp]) + 
      log(prop_dist$norm_consts[comp])
  }
  
  full_mixed_likelihood <- apply(log_like_mat, 1, log_sum_exp)
  
  return(list(importance_samples = all_imp_samps, 
              likelihood = full_mixed_likelihood))
  
}


# Get posterior samples using importance sampling with proposal fitted with WLB samples

# Inputs:
# num_samples - integer, number of importance samples
# apply_psis - boolean, whether to apply Pareto Smoothed Importance Sampling
# input_X - length n vector, transformed X points to be on 0-1 scale
# input_Y - length n vector, Y points
# num_betas - integer, number of parameters for the diffeomorphism 
# num_stationary - integer, number of stationary points 
# first direction - 1 or -1, 1 if the difference between the function value the first 
# prior_mean - p dimensional vector of prior mean (Joint Independent Normal)
# prior_sd - p dimension vector of prior sds (Joint Independent Normal)
# stationary point and f(0) is positve, -1 if it is negative
# zero_is_zero - boolean, TRUE if f(0) = 0, FALSE if f(0) needs to be estimated
# prop_imp_sampler - 3 dimensional vector adding to 1, weights to be given to each of the 
# three components in the proposal density (WLB categorical, fitted mixed normal, prior)
# de_CR - number in [0, 1], crossover probability for differential evolution. Defaults to 0.75.
# de_F - number in [0, 2], differential weighting factor for differential evolution. Defaults to 0.8.
# de_itermax - positive integer, maximum iterations for differential evolution. Defaults to 200.
# apply_psis - boolean, whether to apply Pareto Smoothed Importance Sampling to importance weights
wlb_posterior_samples <- function(num_samples, num_wlb_samples, 
                                  input_X, input_Y, num_betas,
                               num_stationary, first_direction, prior_mean, 
                               prior_sd, zero_is_zero, prop_imp_sampler,
                               de_CR = 0.75, de_F = 0.8, de_itermax = 200, 
                               apply_psis = TRUE) {
  
  # Get WLB samples
  wlb_list <- get_wlb_samples(input_X, input_Y, num_betas, num_stationary,
                                first_direction, prior_mean, prior_sd, zero_is_zero, 
                                num_wlb_samples, de_CR = de_CR, de_F = de_F, 
                              de_itermax = de_itermax)
  
  # Get Importance Samples
  importance_samples <- wlb_importance_sampler(num_samples, wlb_list, prior_mean, 
                                               prior_sd, prop_imp_sampler, num_betas)
  # Get posterior likelihood of Importance Samples
  posterior_likelihood_imp <- apply(importance_samples$importance_samples, 1, 
                                bayes_value, num_betas = num_betas,
                                input_X = input_X, input_Y = input_Y, 
                                first_direction = first_direction,
                                prior_mean = prior_mean, prior_sd = prior_sd, 
                                zero_is_zero = zero_is_zero,
                                weight_vector = NULL)
  
  log_importance_weights <- posterior_likelihood_imp - importance_samples$likelihood
  
  # Apply PSIS to importance weights
  if (apply_psis){
    log_importance_weights<- psis(log_importance_weights)$log_weights
  }
  
  # Apply softmax to get weights 
  posterior_imp_weights_mix <- softmax(log_importance_weights)
  
  return(list(samples = importance_samples$importance_samples, 
              weights = posterior_imp_weights_mix))

}


