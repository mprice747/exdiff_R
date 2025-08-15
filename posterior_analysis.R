# Calculate HPD Intervals 
compute_HPD <- function(samples, prob = 0.95) {
  
  n <- length(samples)
  sorted <- sort(samples)
  m <- floor(prob * n)
  
  intervals <- sapply(1:(n - m), function(i) sorted[i + m] - sorted[i])
  min_idx <- which.min(intervals)
  
  return(c(sorted[min_idx], sorted[min_idx + m]))
}


# Make stationary points HPD Intervals (90%, 95%, 99% and 95%
# Bonferonni Corrected Intervals)
make_SP_hpd_intervals <- function(stat_points){
  
  dimnames_list <- list(unlist(lapply(seq(1, ncol(stat_points)), 
                                      function(x){paste0("S", x)})),
                        c("lower", "upper"))
  
  
  intervals_90 <- t(apply(stat_points, 2, compute_HPD, prob = 0.90))
  dimnames(intervals_90) <- dimnames_list
  
  intervals_95 <- t(apply(stat_points, 2, compute_HPD, prob = 0.95))
  dimnames(intervals_95) <- dimnames_list
  
  intervals_99 <- t(apply(stat_points, 2, compute_HPD, prob = 0.99))
  dimnames(intervals_99) <- dimnames_list
  
  bonf_corr <- 1 - 0.05/ncol(stat_points)
  
  intervals_bonf <- t(apply(stat_points, 2, compute_HPD, prob = bonf_corr))
  dimnames(intervals_bonf) <- dimnames_list
  
  intervals_list <- list(intervals_90 = intervals_90, 
                         intervals_95 = intervals_95, 
                         intervals_99 = intervals_99, 
                         intervals_bonf = intervals_bonf)
  
  return(intervals_list)
  
}

# Given posterior samples, return stationary points and posterior estimates for E(Y|X_{i})

# Inputs: 
# posterior samples - m x p matrix of posterior samples
# num_betas - integer, number of parameters for the diffeomorphism 
# first direction - 1 or -1, 1 if the difference between the function value the first 
# zero_is_zero - boolean, TRUE if f(0) = 0, FALSE if f(0) needs to be estimated

# Outputs: 
# stationary_points - m x num_stationary matrix of stationary point estimates
# posterior_predictions - n x m matrix of posterior estimates of E(Y | X_{i})

get_stat_points_predictions <- function(posterior_samples, num_betas, 
                                        first_direction, 
                                        zero_is_zero, interpolation){
  
  
  # Transform posterior samples
  post_process_posterior <- t(apply(posterior_samples, 1,
                                    process_post_samples_as_vec, 
                                    num_betas = num_betas,
                                    first_direction = first_direction, 
                                    zero_is_zero = zero_is_zero))
  
  total_samples <- nrow(post_process_posterior)
  ncol_post <- ncol(post_process_posterior)
  
  # Get Betas, Lambdas and B points
  final_Betas <- post_process_posterior[, 1:num_betas]
  final_Lambdas <- post_process_posterior[, (num_betas + 1):(ncol_post - 1)]
  
  final_BPoints <- matrix(rep(seq(0, 1, length.out = ncol(final_Lambdas)), total_samples), 
                          nrow = total_samples, byrow = TRUE)
  
  # Obtain Posterior predictions for a grid of points
  posterior_predictions <- diffeo_predictions(seq(0, 1, length.out = 500), 
                                              final_Betas, final_BPoints,
                                              final_Lambdas, 
                                              interpolation = interpolation,
                                              transform_X = FALSE, 
                                              return_diffeo = TRUE)
  # Get mean posterior predictive values
  b_vals <- final_BPoints[1, 2:(ncol(final_BPoints) - 1)]
  
  # Get new stationary points via interpolation
  stat_points_fun <- function(diffeo_res){
    return(interp1(x = diffeo_res, y = seq(0, 1, length.out = 500), 
                   xi = b_vals, method = 'linear'))
  }
  
  stationary_points <- t(apply(posterior_predictions$diffeomorphism, 
                               2, stat_points_fun))
  
  # Return stationary points and posterior predictions 
  return(stationary_points)
}


# Sample from posterior of Diffeomorphism Model using multiple chains of MCMC Model, with 
# Laplace approximations at the posterior modes as the covariance matrices. We resample
# using stacking, reducing the effect chains with low posterior probability. These resamples
# are used to create HPD intervals of the stationary points

# X - vector of length n, sample X points
# Y - vector of length n, sample Y points
# num_betas - integer, number of weight parameters for the diffeomorphism 
# num_stationary - integer, number of stationary points 
# first direction - 1 or -1, 1 if the difference between the function value the first 
# stationary point and f(0) is positve, -1 if it is negative
# zero_is_zero - boolean, TRUE if f(0) = 0, FALSE if f(0) needs to be estimated
# interpolation - 'cubic' or 'linear', referring to interpolation type
# prior_mean - p dimensional vector of prior mean (Joint Independent Normal)
# prior_sd - p dimension vector of prior sds (Joint Independent Normal)
# n_chains - positive integer, number of MCMC Chains
# n_per_chain - positive integer, number of samples per chain
# cov_scale - positive real number, multiple proposal covariance matrices by this number
# Generally needs to be tuned, smaller number leads to higher acceptance rate
# stacking_param - positive real number greater than 1, lambda value for solving for optimal simplex
# for stacking procedure

# Outputs: 
# A list with the following elements: 
# diff_stat_points - A length num_stationary list containing the sampled stationary points
# intervals_list_diff - A list of matrices containing 90%, 95%, 99% and 95% Bonferonni Corrected
# HPD Intervals for the Stationary Points
# acceptance_rates - vector of acceptance rates for MCMC Chains (Should be 15-30%)
# chain_weights - vector of chain weights used for resampling 
posterior_stat_points_diff <- function(X, Y, num_betas, num_stationary, first_direction, 
                                       zero_is_zero, interpolation, prior_mean, prior_sd, 
                                       n_chains, n_per_chain, cov_scale,
                                       stacking_param, total_resampled) {
  
  # Transform Data to be [0, 1]
  X_trans <- min_max_transform_X(X)
  input_X <- X_trans$new_X
  input_Y <- Y
  
  # Find modes and inverse Hessians, to be used as proposals for MCMC Chain
  posterior_modes <- find_modes(n_chains, num_betas,
                                input_X, input_Y, first_direction,
                                prior_mean, prior_sd, zero_is_zero,
                                interpolation, 1)
  # Sample from posterior using standard Metropolis Hastings
  mcmc_samples <- mcmc_parallel(posterior_modes, n_per_chain, num_betas,
                                input_X, input_Y, first_direction,
                                prior_mean, prior_sd, zero_is_zero,
                                interpolation, cov_scale)
  
  # Resample MCMC Samples using Stacking
  resampled <- resample_via_stacking(total_resampled, mcmc_samples$mcmc_samples, 
                                     stacking_param, num_betas, input_X, input_Y, 
                                     first_direction, zero_is_zero,
                                     interpolation)
  
  print("Calculating Stationary Points")
  
  # Get stationary points
  stat_points_mcmc <- get_stat_points_predictions(resampled$resampled_mcmc, num_betas, 
                                                  first_direction, zero_is_zero, 
                                                  interpolation)
  
  diff_stat_points_mat <- as.matrix(apply(stat_points_mcmc, 2, reverse_min_max_transform_X, 
                            lower_bound = X_trans$lower_bound, 
                            upper_bound = X_trans$upper_bound))
  
  intervals_list_diff <- make_SP_hpd_intervals(diff_stat_points_mat)
  
  diff_stat_points <- list()
  
  for (i in 1:num_stationary){
    diff_stat_points[[i]] <- diff_stat_points_mat[, i]
  }
  
  # Returns sampled stationary points
  return(list(diff_stat_points = diff_stat_points,
              intervals_list_diff = intervals_list_diff, 
              acceptance_rates = mcmc_samples$acceptance_probs[, 1],
              chain_weights = resampled$mcmc_weights))
}





