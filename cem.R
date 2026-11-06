# Classification EM for 2 clusters, as in the paper: hard assignment, E-step
# with the fitted probabilities of glmer, M-step with one regression per
# cluster (fit_cluster). Covariate model: x | c ~ N(mu_c, sd_c^2).
# sigma = NULL estimates sigma_b in the M-step; sigma = s fixes it.
# Initialisation: k-means on x. The best iteration is kept, since with glmer
# the algorithm is not an exact EM and the objective is not monotone.

source("fit.R")

# P(y = 1) of every patient under the regression of one cluster
# (hospitals with no patients in the cluster: b = 0)
predict_p <- function(f, d) {
  b <- f$b[as.character(d$hosp)]
  b[is.na(b)] <- 0
  plogis(f$beta[1] + f$beta[2] * d$rf + b)
}

cem <- function(d, C = 2, sigma = NULL, max_iter = 50, tol = 1e-6) {
  N <- nrow(d)
  z <- kmeans(d$x, C, nstart = 10)$cluster
  best <- list(ll = -Inf); ll_old <- -Inf; n_worse <- 0
  for (it in 1:max_iter) {
    # M-step
    nc <- tabulate(z, C)
    if (any(nc < 10)) stop("nearly empty cluster")
    w  <- nc / N
    mu <- sapply(1:C, function(c) mean(d$x[z == c]))
    sd <- sapply(1:C, function(c) sqrt(mean((d$x[z == c] - mu[c])^2)))
    fits <- lapply(1:C, function(c) fit_cluster(d[z == c, ], sigma))
    # E-step: L[i, c] = log w_c + log phi(x_i) + log p(y_i | c)
    L <- sapply(1:C, function(c)
      log(w[c]) + dnorm(d$x, mu[c], sd[c], log = TRUE) +
        dbinom(d$y, 1, predict_p(fits[[c]], d), log = TRUE))
    m <- apply(L, 1, max)
    ll <- sum(m + log(rowSums(exp(L - m))))    # mixture log-likelihood (log-sum-exp)
    # C-step: hard assignment
    z_new <- max.col(L)
    if (ll > best$ll) {
      best <- list(ll = ll, mu = mu, fits = fits, z = z, iter = it)
      n_worse <- 0
    } else n_worse <- n_worse + 1
    converged <- all(z_new == z) || abs(ll - ll_old) < tol * abs(ll) || n_worse >= 3
    z <- z_new; ll_old <- ll
    if (converged) break
  }
  # align labels with the true clusters: order by increasing mu (true mu = 0, sep)
  o <- order(best$mu)
  best$fits <- best$fits[o]
  best$z <- match(best$z, o)
  best
}
