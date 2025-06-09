# Apply LogSumExp function to a vector
log_sum_exp <- function(vec){
  
  m_val <- max(vec)
  exp_vec <- exp(vec - m_val)
  
  return(m_val + log(sum(exp_vec)))
}

# Apply SoftMax function to a vector 
softmax <- function(vec) {
  
  m_val <- max(vec)
  exp_vec <- exp(vec - m_val)
  
  return(exp_vec/sum(exp_vec))
  
}

# Check if lambdas height vector is valid (peak, valley, peak, valley, etc)
check_fluctuating_lambdas <- function(lambdas, first_direction) {
  
  lam_len <- length(lambdas)
  
  # True Sign of the difference between the consecutive lambdas
  if (first_direction == 1) {
    true_diff_vec <- rep(c(1, -1), length.out = lam_len - 1)
  } else{
    true_diff_vec <- rep(c(-1, 1), length.out = lam_len - 1)
  }
  
  # See if sign of difference between lambdas matches what it should be
  diff_vec <- sign(lambdas[2:lam_len] - lambdas[1:(lam_len - 1)])
  
  if(sum(diff_vec == true_diff_vec) == lam_len - 1) {
    return(TRUE)
  } else{
    return(FALSE)
  }
}

# Return likelihood value of applying diffeomorphism
# Also includes weight vector for WLB
like_value <- function(betas, lambdas, sigma, input_X, 
                       input_Y, first_direction, 
                       weight_vector = NULL, return_all = FALSE) {
  
  # If norm > pi, return infinity (invalid point)
  if (sum(betas^2) >(pi^2) + 1e-05) {
    return(-Inf)
  }
  
  # Check if lambdas are fluctuating
  if(!check_fluctuating_lambdas(lambdas, first_direction)){
    return(-Inf)
  }
  
  # Checks is sigma is positive
  if (sigma <= 0) {
    return(-Inf)
  }
  
  # Use an equally spaced b_vec
  b_vec <- seq(0, 1, length.out = length(lambdas))
  
  # Return diffeomorphism likelihood values for all (X, Y)
  all_log <- diffeo_all(input_X, input_Y, betas, b_vec, 
                        lambdas, c(sigma),  interpolation = 'cubic', 
                        transform_X = FALSE)
  
  # Return all likelihood values, the sum or 
  # weighted sum for WLB
  if(return_all){
    return(as.vector(all_log$log_likelihood))
  } else{
    if (is.null(weight_vector)){
      return(sum(all_log$log_likelihood))
    } else {
      return(sum(weight_vector * as.vector(all_log$log_likelihood)))
    }
  }
  
}

# Post process samples so that everything is in (-infinty, infinty scale)
# Except for beta_vec which will remain on continuous scale where ||beta|| < pi 
process_post_samples <- function(param_vec, num_betas, first_direction, 
                                 zero_is_zero = FALSE){
  
  total_len <- length(param_vec)
  
  # Get betas
  betas <- param_vec[1:num_betas]
  
  # Get the untransformed lambda variables
  lambdas_raw <- param_vec[(num_betas + 1):(total_len - 1)] 
  
  # Depending on if the first lambda value is 0
  if (zero_is_zero){ 
    correction <- 0
  } else{correction <- -1}
  
  if(first_direction == 1){
    switch_vec <- rep(c(1, -1), 
                      length.out = length(lambdas_raw) + correction)
  } else {
    switch_vec <- rep(c(-1, 1), 
                      length.out = length(lambdas_raw) + correction)
  }
  
  # Transform lambda values into alternating heights
  if(zero_is_zero){
    lambdas <- c(0, cumsum(switch_vec * exp(lambdas_raw)))
  } else{
    lambda_1 <- lambdas_raw[1]
    adding_to <- switch_vec * exp(lambdas_raw[2:length(lambdas_raw)])
    lambdas <- cumsum(c(lambda_1, adding_to))
    
  }
  
  # Turn sigma positive
  sigma <- exp(param_vec[total_len])
  
  return(list(betas = betas, lambdas = lambdas, 
              sigma = sigma))
  
}

# Transform sampled values back into original scale
process_post_samples_as_vec <- function(param_vec, num_betas, first_direction,
                                        zero_is_zero = FALSE) {
  
  process_list <- process_post_samples(param_vec, num_betas, first_direction, 
   zero_is_zero)
  
  return(unlist(process_list, use.names = FALSE))
  
}


# Given unprocessed samples (||betas|| < pi and (-infty, infty)) find log likelihood value
# for entire dataset

# Inputs: 

# param_vec - vector of unprocessed samples
# num_betas - integer, number of parameters for the diffeomorphism 
# input_X - length n vector, transformed X points to be on 0-1 scale
# input_Y - length n vector, Y points
# zero_is_zero - boolean, TRUE if f(0) = 0, FALSE if f(0) needs to be estimated
# weight_vector - NULL or n dimensional vector sum to 1, for WLB, calculates the weighted log likelihood 

# Outputs:
# log_like - Raw or Weighted Log Likelihood for (X, Y) data
like_value_w_processing <- function(param_vec, num_betas,
                                    input_X, input_Y, first_direction,
                                    zero_is_zero = FALSE,
                                     weight_vector = NULL) {
  
  process_values <- process_post_samples(param_vec, num_betas, 
                                         first_direction, 
                                         zero_is_zero)
 
  
  log_like <- like_value(process_values$betas, process_values$lambdas, 
                         process_values$sigma, 
                         input_X, input_Y, first_direction, 
                         weight_vector)
  
  return(log_like)
  
}


# Given unprocessed samples (||betas|| < pi and (-infty, infty)) find log likelihood +
# log prior for entire dataset

# Inputs: 

# param_vec - vector of unprocessed samples
# num_betas - integer, number of parameters for the diffeomorphism 
# input_X - length n vector, transformed X points to be on 0-1 scale
# input_Y - length n vector, Y points
# zero_is_zero - boolean, TRUE if f(0) = 0, FALSE if f(0) needs to be estimated
# weight_vector - NULL or n dimensional vector sum to 1, for WLB, calculates the weighted log likelihood 
# prior_mean - p dimensional vector of prior mean (Joint Independent Normal)
# prior_sd - p dimension vector of prior sds (Joint Independent Normal)
# zero_is_zero - boolean, TRUE if f(0) = 0, FALSE if f(0) needs to be estimated
# weight_vector - NULL or n dimensional vector sum to 1, for WLB, calculates the weighted log likelihood

# Outputs:
# log_like - Raw or Weighted Log Likelihood + Log Prior for (X, Y) data
bayes_value <- function(param_vec, num_betas,
                                 input_X, input_Y, first_direction,
                                prior_mean, prior_sd, 
                        zero_is_zero = FALSE,
                        weight_vector = NULL
                       ) {
  
  # Get log likelihood
  log_like <- like_value_w_processing(param_vec, num_betas,
                                      input_X, input_Y, first_direction,
                                      zero_is_zero = zero_is_zero, 
                                      weight_vector = weight_vector)
 
  
  # Get prior likelihood 
  if ((is.null(prior_sd)) | (is.null(prior_mean))){
    return(log_like)
  } else{

    log_prior <- sum(dnorm(param_vec, mean = prior_mean, 
                           sd = prior_sd, log = TRUE))
    return(log_like + log_prior)
  }
}

# Get negative of bayes_value, used for minimization
neg_bayes_value <- function(param_vec, num_betas,
                            input_X, input_Y, first_direction,
                            prior_mean, prior_sd, zero_is_zero = FALSE,
                            weight_vector = NULL){
  
  return(-1 * bayes_value(param_vec, num_betas,
                          input_X, input_Y, first_direction,
                          prior_mean, prior_sd, 
                          zero_is_zero,
                          weight_vector))
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
                                  zero_is_zero){
  
  
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
                                              interpolation = 'cubic',
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
  return(list(stationary_points = stationary_points, 
              posterior_predictions = posterior_predictions))
}

# Given data and vector of weights, estimate a weighted kernel density and calculated a
# HPD interval
kde_hpd_interval <- function(data_vec, weights, conf) {
  
  # Estimated weighted kde
  kde_density <- density(data_vec, bw = "nrd", weights = weights, 
                         n = 10000, from = 0, to = 1)
  
  # Estimate cdf via trapezoid area
  heights <- diff(kde_density$x)
  trapz_areas <- (kde_density$y[1:9999] + kde_density$y[2:10000])/2 * heights
  cdf_estimate <- c(0, cumsum(trapz_areas))
  
  # Sample via inverse cdf method and calculate highest posterior density
  icdf_sample <- approx(cdf_estimate, kde_density$x, runif(10000))$y

  hdi_int <- hdi(icdf_sample, ci = conf)
  
  return(c(hdi_int$CI_low, hdi_int$CI_high))
}

# Create a credible interval for the endpoints from the weights 
cred_interval <- function(data_vec, weights, endpoints){
  
  order_data_vec <- order(data_vec)
  
  # Create CDF
  x_axis <- cumsum(weights[order_data_vec])
  y_axis <- data_vec[order_data_vec]
  
  # Use linear interpolation to estimate CDF values 
  return(interp1(x_axis, y_axis, endpoints, 'linear'))
  
  
}

# Create credible intervals for each posterior dimension from probability weights 

# Inputs: 
# posterior_samples - m x p dimensional matrix representing posterior samples
# weights - p dimensional probability vector 
# conf - number between 0 and 1, confidence level for credible interval 

# Outputs:
# cred_intervals - p posterior credible intervals at specified confidence level 
weighted_cred_intervals <- function(posterior_samples, weights, conf) {
  
  alpha_val <- (1 - conf)/2
  
  endpoints <- c(alpha_val, 1 - alpha_val)
  
  cred_intervals <- apply(posterior_samples, 2, cred_interval, 
                          weights = weights, endpoints = endpoints)
  
  return(cred_intervals)
  
}

# Create HPD intervals for each posterior dimension from probability weights 

# Inputs: 
# posterior_samples - m x p dimensional matrix representing posterior samples
# weights - p dimensional probability vector 
# conf - number between 0 and 1, confidence level for HPD interval 

# Outputs:
# hpd_intervals - p posterior HPD intervals at specified confidence level 
weighted_hpd_intervals <- function(posterior_samples, weights, conf) {
  
  hpd_intervals <- apply(posterior_samples, 2, kde_hpd_interval, 
                         weights = weights, conf = conf)
  
  return(hpd_intervals)
  
}



