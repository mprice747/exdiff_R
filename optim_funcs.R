# Sample from a multivariate normal where ||x[1:dim]|| < pi. Uses accept/reject strategy.
# We call this a circular multivariate normal distribution
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


# Uses BFGS to find local minimum of negative log likelihood or negative log likelihood
# + w_0 *log prior

# Inputs: 
# num_betas - integer, number of parameters for the diffeomorphism 
# input_X - length n vector, transformed X points to be on 0-1 scale
# input_Y - length n vector, Y points
# first direction - 1 or -1, 1 if the difference between the function value the first 
# stationary point and f(0) is positve, -1 if it is negative
# prior_mean - p dimensional vector of prior mean (Joint Independent Normal)
# prior_sd - p dimension vector of prior sds (Joint Independent Normal)
# zero_is_zero - boolean, TRUE if f(0) = 0, FALSE if f(0) needs to be estimated
# interpolation - 'cubic' or 'linear', referring to interpolation type
# w_0 - positive real number, weight to be added to prior 

# Outputs: 
# A list with the following elements
# min_par - estimated local minimum
# min_cost - cost function value of min_value
# laplace_approx - inverse of Hessian at min_par
find_local_minimum <- function(num_betas,
                              input_X, input_Y, first_direction,
                              prior_mean, prior_sd, zero_is_zero,
                              interpolation, w_0) {
  
  found_result <- FALSE
  
  while(!found_result){
    
    bfgs_result <- tryCatch({
      
      initial_point <- rcnorm_accept_reject(1, prior_mean, diag(prior_sd^2), 
                                            pi, num_betas)
      
      # BFGS Optimization with initial starting point
      final_res <- optim(as.vector(initial_point), 
                         fn = neg_bayes_value, 
                         method = 'BFGS', control = list(maxit = 10000),
                         num_betas = num_betas,
                         input_X = input_X, input_Y = input_Y,
                         first_direction = first_direction,
                         prior_mean = prior_mean,
                         prior_sd = prior_sd, 
                         zero_is_zero = zero_is_zero,
                         interpolation = interpolation,
                         weight_vector = NULL, 
                         w_0 = w_0, 
                         hessian = TRUE)
      
      
      if (sum(eigen(final_res$hessian)$values <= 0) == 0){
        found_result <- TRUE
      }
      
      # Return parameter, the optimization cost, and Inverse of Hessian
      # (Laplace Approximation)
      list(min_par = final_res$par,
           min_cost = final_res$value, 
           laplace_approx = solve(final_res$hessian))
      
    },  error = function(e){
      
      
    })
    
  }
  
  return(bfgs_result)
}

# Runs BFGS multiple times to find modes of likelihood or posterior distribution. 

# Inputs: 
# num_searches - positive integer, Number of runs of BFGS to find modes
# num_betas - integer, number of parameters for the diffeomorphism 
# input_X - length n vector, transformed X points to be on 0-1 scale
# input_Y - length n vector, Y points
# first direction - 1 or -1, 1 if the difference between the function value the first 
# stationary point and f(0) is positve, -1 if it is negative
# prior_mean - p dimensional vector of prior mean (Joint Independent Normal)
# prior_sd - p dimension vector of prior sds (Joint Independent Normal)
# zero_is_zero - boolean, TRUE if f(0) = 0, FALSE if f(0) needs to be estimated
# interpolation - 'cubic' or 'linear', referring to interpolation type
# w_0 - positive real number, weight to be added to prior 

# Outputs: 
# A list of lists with each list containing the following elements: 
# min_par - estimated mode
# min_cost - cost function value of min_par
# laplace_approx - inverse of Hessian at min_par
find_modes <- function(num_searches, num_betas,
                       input_X, input_Y, first_direction,
                       prior_mean, prior_sd, zero_is_zero,
                       interpolation, w_0){
  
  print("Finding Modes:")
  cl <- makeCluster(detectCores())
  registerDoSNOW(cl)
  funcs <- ls(envir = .GlobalEnv)
  
  pb <- txtProgressBar(max = num_searches, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Obtain WLB samples with parallelization 
  mode_samples <- foreach(i = 1:num_searches, .packages = c('pracma', 'stats', 'MASS', 
                                                                   'mvtnorm', 'DEoptim'), 
                              .export = funcs, .options.snow = opts) %dopar% {
                                
                                
                                mode_result <- find_local_minimum(num_betas,
                                                                  input_X, input_Y, first_direction,
                                                                  prior_mean, prior_sd, zero_is_zero,
                                                                  interpolation, w_0)
                                
                                mode_result
                                
                                
                              }
  
  close(pb)
  stopCluster(cl)
  
  return(mode_samples)
  
}

# Finds MAP of diffeomorphism model. Runs BFGS multiple times and takes value
# with minimum cost

# Inputs: 
# num_searches - positive integer, Number of runs of BFGS to find modes
# num_betas - integer, number of parameters for the diffeomorphism 
# input_X - length n vector, transformed X points to be on 0-1 scale
# input_Y - length n vector, Y points
# first direction - 1 or -1, 1 if the difference between the function value the first 
# stationary point and f(0) is positve, -1 if it is negative
# prior_mean - p dimensional vector of prior mean (Joint Independent Normal)
# prior_sd - p dimension vector of prior sds (Joint Independent Normal)
# zero_is_zero - boolean, TRUE if f(0) = 0, FALSE if f(0) needs to be estimated
# interpolation - 'cubic' or 'linear', referring to interpolation type
# w_0 - positive real number, weight to be added to prior 

# Outputs: 
# A list with the following elements
# min_par - estimated MLE or MAP
# min_cost - cost function value of min_value
# laplace_approx - inverse of Hessian at min_par
find_map <- function(num_searches, num_betas,
                         input_X, input_Y, first_direction,
                         prior_mean, prior_sd, zero_is_zero,
                         interpolation) {
    
  modes_found <- find_modes(num_searches, num_betas,input_X, input_Y, 
                              first_direction, prior_mean, prior_sd, zero_is_zero,
                              interpolation, 1)

  
  min_list <- modes_found[[which.min(sapply(modes_found, function(x) x$min_cost))]]
  
  return(min_list)
  
}


