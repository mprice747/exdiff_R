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



# Runs a single chain of Metropolis Hastings with proposal N(x_{t -1}, \Sigma) 
# for the Diffeomorphism Model

# n - Positive integer, Number of MCMC samples
# proposal_mean - vector, Mean for initial proposal distribution
# proposal_cov - p x p positive definite matrix, Covariance for proposal distribution
# num_betas - integer, number of parameters for the diffeomorphism 
# input_X - length n vector, transformed X points to be on 0-1 scale
# input_Y - length n vector, Y points
# first direction - 1 or -1, 1 if the difference between the function value the first 
# stationary point and f(0) is positve, -1 if it is negative
# prior_mean - p dimensional vector of prior mean (Joint Independent Normal)
# prior_sd - p dimension vector of prior sds (Joint Independent Normal)
# zero_is_zero - boolean, TRUE if f(0) = 0, FALSE if f(0) needs to be estimated
# interpolation - 'cubic' or 'linear', referring to interpolation type
# dim_correction - boolean, whether to multiply dim_parameter * (2.4^2)/p to 
# the proposal covairance (generally recommended)
# dim_parameter - positive real number, Generally needs to be tuned, smaller number lead 
# to higher acceptance rates

# Outputs: 
# A list with the following elements: 
# mcmc_samples - n x p matrix containing MCMC samples
# acp_prob - MCMC acceptance rate (Should be around 15-30%)
mcmc_one_chain <- function(n, proposal_mean, proposal_cov, num_betas,
                           input_X, input_Y, first_direction,
                           prior_mean, prior_sd, zero_is_zero,
                           interpolation, dim_correction = TRUE, 
                           dim_paramater = 1) {
  
  # Apply correction to proposal covariance
  if(dim_correction) {
    proposal_cov <- dim_parameter * (2.4^2)/nrow(proposal_cov) * proposal_cov
  }
  
  mcmc_samples <- matrix(rep(0, n *nrow(proposal_cov) ), nrow = n)
  
  # First sample
  old_samp <- rcnorm_accept_reject(1, proposal_mean, proposal_cov, pi, num_betas)
  
  # Log Posterior of old sample
  old_post_like <- bayes_value(old_samp[1, ], num_betas,
                               input_X, input_Y, first_direction,
                               prior_mean, prior_sd, 
                               zero_is_zero,
                               interpolation)
  
  num_accepted <- 0
  for (i in 1:n){
    
    # Get new sample
    new_samp <- rcnorm_accept_reject(1, old_samp, proposal_cov, pi, num_betas)
    
    # Log Posterior of new sample
    new_post_like <- bayes_value(new_samp[1, ], num_betas,
                                 input_X, input_Y, first_direction,
                                 prior_mean, prior_sd, 
                                 zero_is_zero,
                                 interpolation)
    
    new_given_old <- log(get_normalize_const_cnorm(old_samp, proposal_cov, pi, num_betas))
    
    old_given_new <- log(get_normalize_const_cnorm(new_samp, proposal_cov, pi, num_betas))
    
    # Calculate Metropolis Hastings ratio
    mh_ratio <- exp(new_post_like + old_given_new - old_post_like - new_given_old)
    
    if(runif(1) <= mh_ratio){
      
      # If accepted, update sample
      old_samp <- new_samp
      old_post_like <- new_post_like
      
      num_accepted <-  num_accepted + 1
    }
    
    mcmc_samples[i, ] <- old_samp[1, ]
    
  }
  return(list(mcmc_samples = mcmc_samples, acp_prob = num_accepted/n))
}

# Run multiple parallel chains of Metropolis Hastings for the Diffeomorphism Model


# posterior_modes - length m list of output from find_modes, will contain proposal_mean 
# proposal_cov's for each MCMC chain
# n_per_chain - Positive integer, Number of samples per chain
# num_betas - integer, number of parameters for the diffeomorphism 
# input_X - length n vector, transformed X points to be on 0-1 scale
# input_Y - length n vector, Y points
# first direction - 1 or -1, 1 if the difference between the function value the first 
# stationary point and f(0) is positve, -1 if it is negative
# prior_mean - p dimensional vector of prior mean (Joint Independent Normal)
# prior_sd - p dimension vector of prior sds (Joint Independent Normal)
# zero_is_zero - boolean, TRUE if f(0) = 0, FALSE if f(0) needs to be estimated
# interpolation - 'cubic' or 'linear', referring to interpolation type
# dim_parameter - positive real number, Generally needs to be tuned, smaller number lead 
# to higher acceptance rates

# Outputs
# List with following elements
# mcmc_samples - List of matrices of size n_per_chain x p, all MCMC Samples
# acceptance_probs - vector of acceptance rate for each MCMC Chain
mcmc_parallel <- function(posterior_modes, n_per_chain, num_betas,
                          input_X, input_Y, first_direction,
                          prior_mean, prior_sd, zero_is_zero,
                          interpolation, dim_parameter) {
  
  # Make a progress bar for MCMC CHain
  print("MCMC Progress:")
  cl <- makeCluster(detectCores())
  registerDoSNOW(cl)
  funcs <- ls(envir = .GlobalEnv)
  
  pb <- txtProgressBar(max = length(posterior_modes), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Obtain MCMC Chains with Parallelizations
  mcmc_samples_all <- foreach(i = 1:length(posterior_modes), .packages = c('mvtnorm', 'pracma', 
                                                                       'CompQuadForm', 'MASS', 'stats'), 
                              .export = funcs, .options.snow = opts) %dopar% {
                              
                                mcmc_samples <- mcmc_one_chain(n_per_chain, posterior_modes[[i]]$min_par,
                                                               posterior_modes[[i]]$laplace_approx, 
                                                               num_betas,
                                                               input_X, input_Y, first_direction,
                                                               prior_mean, prior_sd, 
                                                               zero_is_zero,
                                                               interpolation, 
                                                               dim_correction = TRUE, 
                                                               dim_paramater = dim_parameter)
                                
                                mcmc_samples
                                
                              }
  
  close(pb)
  stopCluster(cl)
  
  mcmc_samples <- lapply(mcmc_samples_all, function(x){x$mcmc_samples})
  acceptance_probs <- do.call(rbind, lapply(mcmc_samples_all, function(x){x$acp_prob}))
  
  return(list(mcmc_samples = mcmc_samples, acceptance_probs = acceptance_probs))
}












