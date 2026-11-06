# Mixed-effects logistic regression of one cluster: y ~ rf + (1 | hosp).
#
#   sigma = NULL   glmer, sigma_b estimated (as in the paper)
#   sigma = s      same lme4 Laplace likelihood with sigma_b fixed at s,
#                  minimised over the fixed effects only
#
# Both return: beta, beta_se (Wald), b (conditional modes, used by the CEM),
# sigma_b (estimated or fixed) and model (a merMod object, for REsim).

suppressPackageStartupMessages({library(lme4); library(merTools)})

FORMULA <- y ~ rf + (1 | hosp)

# Hospital intervals of Fig. 2 of the paper: merTools::REsim with its default
# 200 draws; returns the mean and sd of each b_j. The random number state is
# restored, so the draws do not change the rest of the simulation.
resim_b <- function(m, lev, n.sims = 200) {
  had <- exists(".Random.seed", envir = globalenv())
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()))
  r <- suppressWarnings(REsim(m, n.sims = n.sims))
  i <- match(lev, r$groupID)
  list(mean = r$mean[i], sd = r$sd[i])
}

fit_cluster <- function(d, sigma = NULL) {
  d$hosp <- droplevels(d$hosp)
  if (is.null(sigma)) {
    m <- suppressWarnings(suppressMessages(
      glmer(FORMULA, data = d, family = binomial,
            control = glmerControl(optimizer = "bobyqa"))))
    cf <- summary(m)$coefficients
    re <- ranef(m)$hosp
    return(list(beta = cf[, 1], beta_se = cf[, 2],
                b = setNames(re[, 1], rownames(re)),
                sigma_b = attr(VarCorr(m)$hosp, "stddev"), model = m))
  }

  # sigma_b fixed. devf(c(theta, beta)) runs PIRLS and returns the Laplace
  # deviance -2 l_LA; for a random intercept, theta = sigma_b.
  devf <- suppressWarnings(suppressMessages(
    glmer(FORMULA, data = d, family = binomial, devFunOnly = TRUE)))
  obj <- function(beta) devf(c(sigma, beta))
  start <- coef(glm(y ~ rf, data = d, family = binomial))
  opt <- optim(start, obj, method = "BFGS")
  V <- 2 * solve(optimHess(opt$par, obj))      # deviance = -2 loglik
  devf(c(sigma, opt$par))                      # leave the internal state at the optimum
  env <- environment(devf)
  lev <- levels(d$hosp)
  b <- setNames(as.vector(env$pp$b(1)), lev)   # conditional modes from the last PIRLS
  # merMod object with sigma_b fixed, built as glmer does, so that REsim applies
  lf <- glFormula(FORMULA, data = d, family = binomial)
  m <- mkMerMod(env, list(par = c(sigma, opt$par), fval = devf(c(sigma, opt$par)),
                          conv = 0, message = ""), lf$reTrms, fr = lf$fr)
  list(beta = setNames(opt$par, names(start)), beta_se = sqrt(diag(V)),
       b = b, sigma_b = sigma, model = m)
}
