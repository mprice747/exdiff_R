num_betas <- 4
num_stationary <- 2
first_direction <- 1
zero_is_zero <- FALSE

prop_imp_sampler <- c(0.05, 0.9, 0.05)
apply_psis <- TRUE
de_CR <- 0.75
de_F <- 0.8
de_itermax <- 200

prior_mean <- c(rep(0, num_betas), mean(input_Y),
                rep(0, num_stationary + 1), 0)
prior_sd <- c(rep(pi, num_betas), 2 * sd(input_Y), 
              rep(max(input_Y) - min(input_Y), num_stationary + 1), 
              0.5)

set.seed(1247)
X <- runif(50, 0, 2)
Y <- 0.3 + 0.4 * X + 0.5 * sin(3.2* X) + 1.1/(1 + X^2) + rnorm(50, sd = 0.25)

X_trans <- min_max_transform_X(X)
input_X <- X_trans$new_X
input_Y <- Y



num_samples <- 250000
num_wlb_samples <- 500


post_samples <- wlb_posterior_samples(num_samples, num_wlb_samples, 
                                      input_X, input_Y, num_betas,
                                      num_stationary, first_direction, prior_mean, 
                                      prior_sd, zero_is_zero, prop_imp_sampler,
                                      de_CR = 0.75, de_F = 0.8, de_itermax = 200, 
                                      apply_psis = TRUE)

posterior_sps_predictions <- get_stat_points_predictions(post_samples$samples, num_betas, 
                                                         first_direction, zero_is_zero)

stat_point_intervals <- weighted_cred_intervals(posterior_sps_predictions$stationary_points,
                                               post_samples$weights, 0.95)

t(apply(stat_point_intervals, 2, reverse_min_max_transform_X, lower_bound = X_trans$lower_bound, 
        upper_bound = X_trans$upper_bound))

Y_predictive <- posterior_sps_predictions$posterior_predictions$predictions %*% post_samples$weights

plot(X, Y)
lines(seq(0, 2, length.out = 500), Y_predictive, col = 'blue')


stat_point_intervals
