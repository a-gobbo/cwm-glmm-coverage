# Data generator.
#
#   z_i ~ Cat(w)                       true cluster, w = (0.85, 0.15)
#   x_i | z_i = c ~ N(mu_c, 1)         clustering covariate, mu = (0, sep)
#   rf_i ~ Bern(0.3)                   binary regressor
#   logit P(y_i = 1) = alpha_c + beta_c * rf_i + b_{j[i],c},   b_{j,c} ~ N(0, sigma_b^2)
#
# With 70 patients per hospital, the large cluster has ~60 patients per hospital
# and the small one ~10. beta differs between clusters, otherwise misassigned
# patients would not bias the estimates.

sim_data <- function(J = 30, n_j = 70, w = c(0.85, 0.15), sep = 3.5,
                     alpha = c(-0.6, -0.4), beta = c(0.2, -1.0),
                     p_rf = 0.3, sigma_b = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  C <- length(w)
  mu <- c(0, sep)
  b <- matrix(rnorm(J * C, 0, sigma_b), J, C)       # b[j, c]: hospital j in cluster c
  N <- J * n_j
  hosp <- rep(1:J, each = n_j)
  z <- sample(1:C, N, replace = TRUE, prob = w)
  x <- rnorm(N, mu[z], 1)
  rf <- rbinom(N, 1, p_rf)
  y <- rbinom(N, 1, plogis(alpha[z] + beta[z] * rf + b[cbind(hosp, z)]))
  list(data = data.frame(y, x, rf, hosp = factor(hosp), z_true = z),
       b_true = b,
       truth = list(w = w, mu = mu, alpha = alpha, beta = beta, sigma_b = sigma_b))
}
