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
                       input_Y, first_direction, interpolation,  
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
                        lambdas, c(sigma),  interpolation = interpolation, 
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
# first direction - 1 or -1, 1 if the difference between the function value the first 
# stationary point and f(0) is positve, -1 if it is negative
# zero_is_zero - boolean, TRUE if f(0) = 0, FALSE if f(0) needs to be estimated
# interpolation - 'cubic' or 'linear', referring to interpolation type
# weight_vector - NULL or n dimensional vector sum to 1, for WLB, calculates the weighted log likelihood 

# Outputs:
# log_like - Raw or Weighted Log Likelihood for (X, Y) data
like_value_w_processing <- function(param_vec, num_betas,
                                    input_X, input_Y, first_direction,
                                    zero_is_zero,
                                    interpolation,
                                     weight_vector = NULL, 
                                    return_all = FALSE) {
  
  process_values <- process_post_samples(param_vec, num_betas, 
                                         first_direction, 
                                         zero_is_zero)
 
  
  log_like <- like_value(process_values$betas, process_values$lambdas, 
                         process_values$sigma, 
                         input_X, input_Y, first_direction, interpolation,
                         weight_vector, return_all)
  
  return(log_like)
}


# Given unprocessed samples (||betas|| < pi and (-infty, infty)) find log likelihood +
# log prior for entire dataset

# Inputs: 

# param_vec - vector of unprocessed samples
# num_betas - integer, number of parameters for the diffeomorphism 
# input_X - length n vector, transformed X points to be on 0-1 scale
# input_Y - length n vector, Y points
# first direction - 1 or -1, 1 if the difference between the function value the first 
# stationary point and f(0) is positve, -1 if it is negative
# prior_mean - p dimensional vector of prior mean (Joint Independent Normal)
# prior_sd - p dimension vector of prior sds (Joint Independent Normal)
# zero_is_zero - boolean, TRUE if f(0) = 0, FALSE if f(0) needs to be estimated
# interpolation - 'cubic' or 'linear', referring to interpolation type
# weight_vector - NULL or n dimensional vector sum to 1, for WLB, calculates the weighted log likelihood
# w_0 - real number > 0, weight to add to prior 

# Outputs:
# log_like - Raw or Weighted Log Likelihood + Log Prior for (X, Y) data
bayes_value <- function(param_vec, num_betas,
                                 input_X, input_Y, first_direction,
                                prior_mean, prior_sd, 
                        zero_is_zero,
                        interpolation,
                        weight_vector = NULL, 
                        w_0 = 1
                       ) {
  
  # Get log likelihood
  log_like <- like_value_w_processing(param_vec, num_betas,
                                      input_X, input_Y, first_direction,
                                      zero_is_zero = zero_is_zero, 
                                      interpolation = interpolation,
                                      weight_vector = weight_vector)
 
  
  # Get prior likelihood 
  if ((is.null(prior_sd)) | (is.null(prior_mean))){
    return(log_like)
  } else{

    log_prior <- sum(w_0 * dnorm(param_vec, mean = prior_mean, 
                           sd = prior_sd, log = TRUE))
    return(log_like + log_prior)
  }
}

# Get negative of bayes_value, used for minimization
neg_bayes_value <- function(param_vec, num_betas,
                            input_X, input_Y, first_direction,
                            prior_mean, prior_sd, zero_is_zero,
                            interpolation,
                            weight_vector = NULL, 
                            w_0 = 1){
  
  return(-1 * bayes_value(param_vec, num_betas,
                          input_X, input_Y, first_direction,
                          prior_mean, prior_sd, 
                          zero_is_zero, interpolation,
                          weight_vector, 
                          w_0))
}

