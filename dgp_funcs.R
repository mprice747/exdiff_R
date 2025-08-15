# Obtain intervals of stationary points using DGP method 

# Inputs: 
# sample_t - sampled univariate stationary points from DGP stationary point posterior
# num_stationary - positive integer, number of stationary points

# Outputs: 
# A list with the following elements
# dgp_stat_points - list containing posterior samples of each stationary point 
# intervals_list - list containing 90, 95, 99 and Bonferonni corrrected intervals for each 
# stationary point
dgp_post_processing <- function(sample_t, num_stationary){
  
  # Run K means on sample t's
  dgp_kmeans <- kmeans(sample_t, num_stationary)
  centers_order <- order(dgp_kmeans$centers[, 1])
  
  dgp_stat_points <- list()
  
  # Get sampled points for each cluster
  for (i in 1:num_stationary){
    dgp_stat_points[[i]] <- sample_t[dgp_kmeans$cluster == centers_order[i]]
  }
  
  dimnames_list <- list(unlist(lapply(seq(1, length(dgp_stat_points)), 
                                      function(x){paste0("S", x)})),
                        c("lower", "upper"))
  
  # Calculate 90, 95, 99 and Bonferonnu HPD Intervals 
  intervals_90 <- do.call(rbind, lapply(dgp_stat_points, compute_HPD, prob = 0.90))
  dimnames(intervals_90) <- dimnames_list
  
  intervals_95 <- do.call(rbind, lapply(dgp_stat_points, compute_HPD, prob = 0.95))
  dimnames(intervals_95) <- dimnames_list
  
  intervals_99 <- do.call(rbind, lapply(dgp_stat_points, compute_HPD, prob = 0.99))
  dimnames(intervals_99) <- dimnames_list
  
  bonf_corr <- 1 - 0.05/num_stationary
  
  intervals_bonf <- do.call(rbind, lapply(dgp_stat_points, compute_HPD, prob = bonf_corr))
  dimnames(intervals_bonf) <- dimnames_list
  
  intervals_list <- list(intervals_90 = intervals_90, 
                         intervals_95 = intervals_95, 
                         intervals_99 = intervals_99, 
                         intervals_bonf = intervals_bonf)
  
  return(list(dgp_stat_points = dgp_stat_points,
              intervals_list_dgp = intervals_list))
}

# Create difference matrix (Used for kernel of Gaussian process)
create_H0_matrix <- function(x_vec) {
  return(outer(x_vec, x_vec,
               FUN = function(x1, x2) (x1 - x2)))
}

# Sample from the stationary point's posterior for the derivative constrained 
# Gaussian process method

# X - length n vector, x points (input data)
# Y - length n yecotr, y points (output data)
# num_stationary - positive integer, number of stationary points
# total_samples - positive integer, number of posterior samples to obtain
# beta_shape1 - positive real number, First shape parmater for beta prior
# beta_shape2 - positive real number, Second shape parmater for beta prior
posterior_stat_points_dgp <- function(X, Y, num_stationary, total_samples,
                                      beta_shape1 = 1, beta_shape2 = 1) {
  
  # Obtain Gaussian process parameters
  n <- length(X)
  H0_mat <- create_H0_matrix(X)
  
  EB_gp <- Rsolnp::solnp(pars = c(1, 1, 1), fun = log_mar_lik_gp,
                         LB = c(0.0001, 0.0001, 0.0001),
                         UB = c(1 / 0.0001, 1 / 0.0001, 1 / 0.0001),
                         control = list(TOL = 1e-5, trace = 0),
                         y = Y, H0 = H0_mat)$par
  
  
  # Gaussian Process parameters
  lambda_gp <- EB_gp[1]^2 / (n * EB_gp[2]^2)
  Kff_gp <- se_ker(H0 = H0_mat, tau = 1, h = EB_gp[3])
  A_gp <- Kff_gp + diag((n * lambda_gp), n)
  
  # Obtain log Posterior of stationary point parameter
  grid_t <- seq(min(X), max(X), length.out = 5000)
  log_post_gp <- rep(0, 5000)
  
  for (j in 1:5000){
    
    log_post_gp[j] <- log_post_t_theory(t = grid_t[j], y = Y, x = X,
                                        Kff = Kff_gp, A = A_gp, lambda = lambda_gp,
                                        h = EB_gp[3], sig2 = EB_gp[1]^2,
                                        shape1 = beta_shape1, shape2 = beta_shape2, 
                                        a = min(X), b = max(X))
  }
  
  # Obtain unnormalized posterior values at grid_t
  post_prob_gp <- exp(log_post_gp - max(log_post_gp))
  
  # Estimate CDF values of posterior
  cdf_post_prob_gp <- cumtrapz(grid_t, post_prob_gp)
  cdf_post_prob_gp <- cdf_post_prob_gp/cdf_post_prob_gp[5000]
  
  # Obtain Samples by Inverse transform method
  sample_t <- interp1(cdf_post_prob_gp[, 1], grid_t, runif(total_samples))
  
  
  return(dgp_post_processing(sample_t, num_stationary))
  
}

