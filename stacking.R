
# Using one MCMC chain get Leave One Out Posterior Predictive Probability for each data point
# using PSIS smoothing

# Inputs:
# one_chain - m x p matrix representing one MCMC chain
# num_betas - integer, number of parameters for the diffeomorphism 
# input_X - length n vector, transformed X points to be on 0-1 scale
# input_Y - length n vector, Y points
# first direction - 1 or -1, 1 if the difference between the function value the first 
# stationary point and f(0) is positve, -1 if it is negative
# zero_is_zero - boolean, TRUE if f(0) = 0, FALSE if f(0) needs to be estimated
# interpolation - 'cubic' or 'linear', referring to interpolation type

# Outputs:
# loo_estimate - length n vector, estimate of p(y_{i} | y_{-i})
get_loo_estimate <- function(one_chain, num_betas,
                             input_X, input_Y, first_direction,
                             zero_is_zero,
                             interpolation){
  
  # Get log likelihood for each data point for each MCMC chain
  log_like_mat <- apply(one_chain, 1, like_value_w_processing, num_betas = num_betas,
                        input_X = input_X, input_Y = input_Y, first_direction = first_direction,
                        zero_is_zero = zero_is_zero,
                        interpolation = interpolation,
                        weight_vector = NULL, 
                        return_all = TRUE)
  
  # Smooth negative log likelihood using PSIS
  psis_mat <- psis(-1 * log_like_mat)$log_weights
  
  # Importance Sampling estimate of p(y_{i} | y_{-i})
  imp_ratio_num <- apply(psis_mat + log_like_mat, 1, log_sum_exp)
  imp_ratio_dem <- apply(psis_mat, 1, log_sum_exp)
  
  loo_estimate <- exp(imp_ratio_num - imp_ratio_dem)
  
  return(loo_estimate)
}

# Get Leave One Out Posterior Predictive Probability for each data point using PSIS smoothing
# for each MCMC Chain

# Inputs: 
# mcmc_samples - length r list of m x p matrices, each containing samples of a single MCMC Chain
# num_betas - integer, number of parameters for the diffeomorphism 
# input_X - length n vector, transformed X points to be on 0-1 scale
# input_Y - length n vector, Y points
# first direction - 1 or -1, 1 if the difference between the function value the first 
# stationary point and f(0) is positve, -1 if it is negative
# zero_is_zero - boolean, TRUE if f(0) = 0, FALSE if f(0) needs to be estimated
# interpolation - 'cubic' or 'linear', referring to interpolation type

# Outputs: 
# r x n matrix representing each chain's estimate of the Leave One Out Posterior 
# Predictive Probability for each data point
get_loo_estimates <- function(mcmc_samples, num_betas,
                              input_X, input_Y, first_direction,
                              zero_is_zero,
                              interpolation) {
  # LOO Progress PBar
  print("LOO Progress:")
  cl <- makeCluster(detectCores())
  registerDoSNOW(cl)
  funcs <- ls(envir = .GlobalEnv)
  
  pb <- txtProgressBar(max = length(mcmc_samples), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Obtain LOO with parallelization 
  loo_estimates <- foreach(i = 1:length(mcmc_samples), .packages = c('pracma', 'loo', 'stats'), 
                          .export = funcs, .options.snow = opts) %dopar% {
                            
                            
                            
                            loo_estimate <- get_loo_estimate(mcmc_samples[[i]], num_betas,
                                                                input_X, input_Y, first_direction,
                                                                zero_is_zero,
                                                                interpolation)
                            
                            loo_estimate
                            
                          }
  
  close(pb)
  stopCluster(cl)
  
  return(do.call(rbind, loo_estimates))
  
}

# Resample from one MCMC chain
sample_from_mcmc_chain <- function(mcmc_chain, num_samples){
  
  sampled_inds <- sample.int(nrow(mcmc_chain), num_samples, TRUE)
  return(mcmc_chain[sampled_inds, ])
  
}

# Resample from Ensemble of MCMC Chains using Stacking Method

# num_resamples - positive integer, number of total resamples
# mcmc_samples - length r list of m x p matrices, each containing samples of a single MCMC Chain
# prior_lambda - positive real number greater than 1, lambda value for solving for optimal simplex
# for stacking procedure
# num_betas - integer, number of parameters for the diffeomorphism 
# input_X - length n vector, transformed X points to be on 0-1 scale
# input_Y - length n vector, Y points
# first direction - 1 or -1, 1 if the difference between the function value the first 
# stationary point and f(0) is positve, -1 if it is negative
# zero_is_zero - boolean, TRUE if f(0) = 0, FALSE if f(0) needs to be estimated
# interpolation - 'cubic' or 'linear', referring to interpolation type

# Outputs: 
# List with the following elements:
# resampled_mcmc - num_resamples x p matrix, resampled MCMC samples
# mcmc_weights - length r vector, calculated resampled probabilities
resample_via_stacking <- function(num_resamples, mcmc_samples, prior_lambda, num_betas,
                              input_X, input_Y, first_direction, zero_is_zero,
                              interpolation) {
  
  num_chains <- length(mcmc_samples)
  
  # Get LOO estimates for each data point for each chain
  loo_estimates <- t(get_loo_estimates(mcmc_samples, num_betas,
                                       input_X, input_Y, first_direction,
                                       zero_is_zero,
                                       interpolation))
  
  prior_alpha_vec <- rep(prior_lambda, num_chains)
  
  # Optimization function for finding stacking weights for each chain
  # Chains with higher LOO estimates will be weighted higher as that chain can better predict
  # future data points
  stacking_opt <- function(weights){
    
    mixed_log <- sum(log(loo_estimates %*% weights))
    
    log_prior <- ddirichlet(weights, prior_alpha_vec, TRUE)
    
    return(-1 * (mixed_log + log_prior))
  }
  
  # Find optimal simplex
  opt_weights <- solnp(pars = rep(1/num_chains, num_chains), fun = stacking_opt,
    eqfun = function(x){sum(x)}, eqB = 1, ineqfun = function(x){x},
    ineqLB = rep(1e-8, num_chains), ineqUB = rep(1, num_chains), 
    control = list(trace = 0))$pars
  
  # Determine how many samples to take per chain
  sampled_mcmc_chains <- table(sample.int(num_chains, size = num_resamples, 
                                          TRUE, prob = opt_weights))
  
  # Resample MCMC samples
  resampled_mcmc <- lapply(names(sampled_mcmc_chains),function(x){
    sample_from_mcmc_chain(mcmc_samples[[as.integer(x)]], sampled_mcmc_chains[x])})
  
  resampled_mcmc <- do.call(rbind, resampled_mcmc)
  
  return(list(resampled_mcmc = resampled_mcmc, mcmc_weights = opt_weights))
}



