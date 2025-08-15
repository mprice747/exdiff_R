
# WIll containg results
n <- 50
num_sims <- 100

stat_points_1_90_diff  <- matrix(rep(0, num_sims * 2), ncol = 2)
stat_points_1_95_diff  <- matrix(rep(0, num_sims * 2), ncol = 2)
stat_points_1_99_diff  <- matrix(rep(0, num_sims * 2), ncol = 2)

stat_points_1_90_dgp  <- matrix(rep(0, num_sims * 2), ncol = 2)
stat_points_1_95_dgp  <- matrix(rep(0, num_sims * 2), ncol = 2)
stat_points_1_99_dgp  <- matrix(rep(0, num_sims * 2), ncol = 2)

rmse_diff_1 <- rep(0, num_sims)
rmse_dgp_1 <- rep(0, num_sims) 


# Parameters
num_betas <- 4
num_stationary <- 1
first_direction <- 1
zero_is_zero <- FALSE
interpolation <- "cubic"
n_chains <- 100
n_per_chain <- 2000
cov_scale <- 0.9
stacking_param <- 1.0025
total_resampled <- 100000

X_matrix <- matrix(rep(0, num_sims * n), nrow = num_sims)
Y_matrix <- matrix(rep(0, num_sims * n), nrow = num_sims)

for (i in 1:num_sims){
  print(i)

  #0.39973
  X <- runif(n, 0, 1)
  Y <- (-2*X^2) + (3*X) + sin(2*X) + cos(3 *X) + 1 + rnorm(n, sd = 0.2)
  
  #X_matrix[i, ] <- X
  #Y_matrix[i, ] <- Y

  
  prior_mean <- c(rep(0, num_betas), mean(Y),
                  rep(0, num_stationary + 1), 0)
  prior_sd <- c(rep(1, num_betas), 2 * sd(Y), 
                rep(max(Y) - min(Y), num_stationary + 1), 
                0.5)
  
  # HPD Intervals of Stationary Points using Diffeomorphism Method
  diff_method <- posterior_stat_points_diff(X, Y, num_betas, num_stationary, first_direction, 
                              zero_is_zero, interpolation, prior_mean, prior_sd, 
                              n_chains, n_per_chain, cov_scale,
                              stacking_param, total_resampled)
  
  # HPD Intervals of Stationary Points using DGP Method
  dgp_method <- posterior_stat_points_dgp(X, Y, num_stationary, total_resampled)
  
  # Get Stationary Point INtervals
  stat_points_1_90_diff[i, ] <- diff_method$intervals_list_diff$intervals_90[1, ]
  stat_points_1_95_diff[i, ] <- diff_method$intervals_list_diff$intervals_95[1, ]
  stat_points_1_99_diff[i, ] <- diff_method$intervals_list_diff$intervals_99[1, ]
  
  stat_points_1_90_dgp[i, ] <- dgp_method$intervals_list_dgp$intervals_90[1, ]
  stat_points_1_95_dgp[i, ] <- dgp_method$intervals_list_dgp$intervals_95[1, ]
  stat_points_1_99_dgp[i, ] <- dgp_method$intervals_list_dgp$intervals_99[1, ] 
  
  # Calculate RMSE of sampled points
  rmse_diff_1[i] <- sqrt(mean((diff_method$diff_stat_points[[1]] - 0.39973)^2))
  rmse_dgp_1[i] <- sqrt(mean((dgp_method$dgp_stat_points[[1]] - 0.39973)^2))
}


sum((stat_points_1_90_diff[, 1] <  0.39973 ) & (stat_points_1_90_diff[, 2] >  0.39973))
sum((stat_points_1_95_diff[, 1] <  0.39973 ) & (stat_points_1_95_diff[, 2] >  0.39973))
sum((stat_points_1_99_diff[, 1] <  0.39973 ) & (stat_points_1_99_diff[, 2] >  0.39973))

sum((stat_points_1_90_dgp[, 1] <  0.39973 ) & (stat_points_1_90_dgp[, 2] >  0.39973))
sum((stat_points_1_95_dgp[, 1] <  0.39973 ) & (stat_points_1_95_dgp[, 2] >  0.39973))
sum((stat_points_1_99_dgp[, 1] <  0.39973 ) & (stat_points_1_99_dgp[, 2] >  0.39973))
