################################################################################
##  simulation/06_bounded_design_constants.R
##  Computes, once, the fixed constants of the bounded validation design used by
##  03_oracle_study.R and 05_near_cancellation.R with DESIGN=bounded:
##    - covariates clipped to [-3, 3] (Assumption 4.4: bounded covariates and
##      bounded linear predictors);
##    - censoring hazards lc_k fixed in advance (30% censoring before the end of
##      follow-up), so that censoring is independent of the event times;
##    - an administrative end of follow-up tau_k, the 90% quantile of the follow-up
##      time, so that Pr(Y >= tau_k | x) > 0 for every x (Assumption 4.2).
##  Constants are computed from 100,000 participants per study (seed 20290001) and
##  printed; they are copied into BOUNDED in 00_setup.R.
################################################################################
suppressPackageStartupMessages(library(hiermetacox))
here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else "." })
XB <- 3; NPILOT <- 1e5
designs <- list(
  oracle = local({ K <- 5; alpha <- c(1, -0.8, 0.6, -0.5, 0.4, -0.3); eps <- matrix(0, K, 6)
    for (j in 1:4) eps[c(j, j %% K + 1), j] <- sign(alpha[j]) * c(0.5, 1.0)
    list(alpha = alpha, eps = eps) }),
  nearcancel = local({ K <- 5; alpha <- c(1, -0.8, 1.2, -1.2, 0.4, -0.3); eps <- matrix(0, K, 6)
    eps[1:2, 1] <- c(0.5, 1.0); eps[2:3, 2] <- c(-0.5, -1.0)
    eps[3:4, 3] <- c(-1.19, 1.2); eps[4:5, 4] <- c(1.19, 1.19)
    list(alpha = alpha, eps = eps) }))
scales <- c(0.05, 0.075, 0.1, 0.125, 0.15)
for (nm in names(designs)) {
  set.seed(20290001)
  s <- sim_ipd(K = 5, n = NPILOT, p = 100, alpha = designs[[nm]]$alpha, eps = designs[[nm]]$eps,
               scales = scales, rho = 0.5, cens = 0, xbound = XB)        # event times only
  lc <- tau <- numeric(5); B <- numeric(5)
  for (k in 1:5) {
    Tt <- s$data$time[s$data$study == k]
    lc[k] <- uniroot(function(l) mean(1 - exp(-l * Tt)) - 0.30, c(1e-10, 1e4), tol = 1e-12)$root
    Ct <- rexp(length(Tt), lc[k]); tau[k] <- unname(quantile(pmin(Tt, Ct), 0.90))
    B[k] <- XB * sum(abs(s$theta[k, ]))                                   # bound on |x'theta_0k|
  }
  cat(sprintf('%s: cens.rate = c(%s), tau = c(%s); max |x\'theta| <= %.1f\n', nm,
              paste(signif(lc, 4), collapse = ", "), paste(signif(tau, 4), collapse = ", "), max(B)))
}
