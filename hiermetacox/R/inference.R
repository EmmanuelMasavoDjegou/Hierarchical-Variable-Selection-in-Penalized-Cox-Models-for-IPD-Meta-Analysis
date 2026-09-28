#' Post-selection (oracle-type) Wald inference on the selected support
#'
#' Refits an unpenalized stratified Cox model (Breslow ties, one stratum per
#' study) on the support selected by the penalized fit: one column per selected
#' global effect and one column per nonzero study-specific deviation. By the
#' oracle property (Theorems 1-2 of the accompanying paper), the resulting Wald
#' standard errors are asymptotically valid on the event of correct selection;
#' in finite samples they are conditional on the selected model and should be
#' read together with the bootstrap selection frequencies from
#' \code{\link{hmcox_boot}}.
#'
#' @param object a \code{cv.hmcox} (or \code{hmcox}) fit.
#' @param rule tuning rule ("min" or "1se") of the \code{cv.hmcox} fit.
#' @param level confidence level.
#' @return list with \code{global} (one row per selected covariate: penalized
#'   estimate, refitted log-HR, SE, HR and CI, Wald p-value) and \code{study}
#'   (study-specific log-HR theta_kj with SEs for selected covariates).
#' @export
hmcox_inference <- function(object, rule = object$rule, level = 0.95) {
  D <- object$data
  cf <- if (inherits(object, "cv.hmcox")) coef.cv.hmcox(object, rule) else
    coef.hmcox(object, length(object$lambda))
  A <- which(cf$alpha != 0)
  if (length(A) == 0) return(list(global = NULL, study = NULL))
  z <- qnorm(1 - (1 - level) / 2)
  cols <- list(); lab <- character(0); map <- list()
  for (j in A) {
    cols[[length(cols) + 1]] <- D$X[, j]
    lab <- c(lab, paste0("a_", D$xnames[j])); map[[length(map) + 1]] <- c(j, 0)
  }
  for (j in A) for (k in seq_len(D$K)) if (cf$eps[k, j] != 0) {
    cols[[length(cols) + 1]] <- D$X[, j] * (D$sid == k)
    lab <- c(lab, paste0("e", k, "_", D$xnames[j])); map[[length(map) + 1]] <- c(j, k)
  }
  Z <- do.call(cbind, cols); colnames(Z) <- make.names(lab)
  map <- do.call(rbind, map)
  ## 'strata' must appear unqualified in the formula to be recognised as a
  ## special; it is imported from survival in NAMESPACE.
  df <- data.frame(time = D$time, status = D$status, sid = D$sid)
  fit <- survival::coxph(Surv(time, status) ~ Z + strata(sid), data = df, ties = "breslow")
  b <- stats::coef(fit); V <- stats::vcov(fit)
  if (anyNA(b)) warning("refit is singular for some columns (collinear support).")
  ga <- which(map[, 2] == 0)
  se <- sqrt(diag(V))[ga]
  glob <- data.frame(covariate = D$xnames[map[ga, 1]],
                     alpha_penalized = cf$alpha[map[ga, 1]],
                     logHR = b[ga], SE = se,
                     HR = exp(b[ga]), HR_lower = exp(b[ga] - z * se),
                     HR_upper = exp(b[ga] + z * se),
                     p_value = 2 * pnorm(-abs(b[ga] / se)),
                     n_deviating_studies = sapply(map[ga, 1], function(j) sum(cf$eps[, j] != 0)),
                     row.names = NULL, stringsAsFactors = FALSE)
  st <- list()
  for (j in A) for (k in seq_len(D$K)) {
    w <- numeric(length(b)); w[ga[map[ga, 1] == j]] <- 1
    ek <- which(map[, 1] == j & map[, 2] == k); if (length(ek)) w[ek] <- 1
    est <- sum(w * b); s <- sqrt(drop(t(w) %*% V %*% w))
    st[[length(st) + 1]] <- data.frame(covariate = D$xnames[j], study = D$levels[k],
                                       theta_penalized = cf$theta[k, j],
                                       deviation_selected = cf$eps[k, j] != 0,
                                       logHR = est, SE = s, HR = exp(est),
                                       HR_lower = exp(est - z * s), HR_upper = exp(est + z * s),
                                       stringsAsFactors = FALSE)
  }
  list(global = glob, study = do.call(rbind, st), refit = fit)
}

#' Stratified bootstrap of selection and estimation
#'
#' Resamples subjects with replacement within each study, refits the model at the
#' tuning parameters selected by cross-validation, and records how often each
#' covariate (global effect) and each deviation is selected, together with
#' percentile intervals of the penalized estimates.
#' @param object a \code{cv.hmcox} fit.
#' @param B number of bootstrap resamples.
#' @param seed optional seed.
#' @param level percentile interval level.
#' @param rule tuning rule whose (lambda, ratio) is held fixed.
#' @return list with \code{selection} (per covariate), \code{deviation_freq}
#'   (K x p) and the matrix of bootstrap global estimates.
#' @export
hmcox_boot <- function(object, B = 200, seed = NULL, level = 0.95, rule = object$rule) {
  if (!is.null(seed)) set.seed(seed)
  D <- object$data
  ix <- if (rule == "1se") object$sel.1se else object$sel.min
  lam <- object$lambda.grid[seq_len(ix[1])]
  ratio <- object$ratio.grid[ix[2]]
  p <- D$p; K <- D$K
  Ab <- matrix(NA_real_, B, p, dimnames = list(NULL, D$xnames))
  Ef <- matrix(0, K, p)
  ok <- 0
  for (b in seq_len(B)) {
    idx <- unlist(lapply(seq_len(K), function(k) {
      r <- (D$sstart[k] + 1):D$sstart[k + 1]; r[sample.int(length(r), replace = TRUE)] }))
    Db <- D
    Db$X <- D$X[idx, , drop = FALSE]; Db$time <- D$time[idx]; Db$status <- D$status[idx]
    Db$sid <- D$sid[idx]
    o <- order(Db$sid, Db$time, -Db$status)
    Db <- hmc_subset(within_order(Db, o), rep(TRUE, length(idx)))
    f <- try(hmc_fit_prepared(Db, object$penalty, object$gamma, ratio, lam,
                              object$heterogeneity, 1e-5, 50, 2000, p), silent = TRUE)
    if (inherits(f, "try-error") || length(f$lambda) < length(lam)) next
    ok <- ok + 1
    Ab[ok, ] <- f$alpha[, length(lam)]
    Ef <- Ef + (f$eps[, , length(lam)] != 0)
  }
  Ab <- Ab[seq_len(ok), , drop = FALSE]
  lo <- (1 - level) / 2
  sel <- data.frame(covariate = D$xnames,
                    selection_freq = colMeans(Ab != 0),
                    boot_lower = apply(Ab, 2, quantile, lo),
                    boot_upper = apply(Ab, 2, quantile, 1 - lo),
                    stringsAsFactors = FALSE)
  list(selection = sel, deviation_freq = Ef / max(ok, 1), alpha_boot = Ab, B_used = ok)
}

within_order <- function(D, o) {
  D$X <- D$X[o, , drop = FALSE]; D$time <- D$time[o]; D$status <- D$status[o]
  D$sid <- D$sid[o]; D
}

#' Simulate multi-study survival data under the hierarchical model
#'
#' Covariates are AR(1) Gaussian with correlation \code{rho}; study k has a
#' Weibull baseline hazard with shape \code{shape} and scale drawn from
#' \code{scale.range}; \eqn{\theta_k=\alpha+\varepsilon_k}. For each covariate in
#' \code{het.idx}, \code{n.dev} randomly chosen studies receive a deviation
#' \eqn{\varepsilon_{kj}\sim N(0,\sigma_\varepsilon^2)} (all K studies if
#' \code{n.dev = K}); all other deviations are zero. Independent exponential
#' censoring is calibrated within each study by root-finding so that the
#' expected censoring proportion equals \code{cens}.
#' @param K number of studies.
#' @param n study sample sizes (recycled to length K).
#' @param p number of covariates.
#' @param alpha nonzero global effects (the first \code{length(alpha)} covariates).
#' @param rho AR(1) correlation of the covariates.
#' @param cens target censoring proportion in each study.
#' @param het.idx covariates that can carry study-specific deviations.
#' @param n.dev number of studies (chosen at random) that deviate on each of them.
#' @param sd.eps standard deviation of the nonzero deviations.
#' @param shape Weibull shape of the baseline hazards.
#' @param scale.range range of the study-specific Weibull scales.
#' @param eps optional fixed K x p0 matrix of deviations (p0 <= p; remaining columns
#'   zero); when given, \code{het.idx}, \code{n.dev} and \code{sd.eps} are ignored.
#' @param scales optional fixed Weibull scales of the K baseline hazards.
#' @return list with \code{data} (data frame: study, time, status, X1..Xp),
#'   \code{alpha}, \code{eps} (K x p) and \code{theta} (K x p).
#' @export
sim_ipd <- function(K = 5, n = 200, p = 100, alpha = c(1, -0.8, 0.6, -0.5, 0.4, -0.3),
                    rho = 0.5, cens = 0.3, het.idx = 1:4, n.dev = 2, sd.eps = 0.4,
                    shape = 1.5, scale.range = c(0.05, 0.15), eps = NULL, scales = NULL) {
  n <- rep_len(n, K)
  a <- c(alpha, rep(0, p - length(alpha)))
  if (is.null(eps)) {
    eps <- matrix(0, K, p)
    if (sd.eps > 0 && length(het.idx) > 0)
      for (j in het.idx) {
        ks <- if (n.dev >= K) seq_len(K) else sample.int(K, n.dev)
        eps[ks, j] <- rnorm(length(ks), 0, sd.eps)
      }
  } else {                                  # fixed deviations supplied by the user
    stopifnot(nrow(eps) == K, ncol(eps) <= p)
    eps <- cbind(eps, matrix(0, K, p - ncol(eps)))
  }
  theta <- sweep(eps, 2, a, "+")
  Sig <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
  Rch <- chol(Sig)
  lam_k <- if (is.null(scales)) runif(K, scale.range[1], scale.range[2]) else rep_len(scales, K)
  out <- vector("list", K)
  for (k in seq_len(K)) {
    X <- matrix(rnorm(n[k] * p), n[k], p) %*% Rch
    eta <- drop(X %*% theta[k, ])
    Tt <- (-log(runif(n[k])) / (lam_k[k] * exp(eta)))^(1 / shape)
    if (cens > 0) {
      f <- function(lc) mean(1 - exp(-lc * Tt)) - cens
      lc <- uniroot(f, c(1e-10, 1e4), tol = 1e-10)$root
      Ct <- rexp(n[k], lc)
    } else Ct <- rep(Inf, n[k])
    out[[k]] <- data.frame(study = k, time = pmin(Tt, Ct), status = as.integer(Tt <= Ct),
                           X)
  }
  dat <- do.call(rbind, out)
  names(dat)[-(1:3)] <- paste0("X", seq_len(p))
  list(data = dat, alpha = a, eps = eps, theta = theta)
}

#' Selection and estimation metrics against the truth
#'
#' Metrics are reported at three levels: the global support
#' \eqn{\{j:\alpha_j\neq0\}} (TPR, FDR, F1, MCC, exact recovery, size), the
#' deviation support \eqn{\{(k,j):\varepsilon_{kj}\neq0\}} (TPR_dev, FDR_dev,
#' n_dev), and the coefficients themselves (MSE of \eqn{\alpha},
#' \eqn{\varepsilon} and \eqn{\theta}). \code{err_theta_nc} is the mean absolute
#' error of \eqn{\hat\theta_{kj}} on near-cancelled cells, i.e. true deviations
#' with \eqn{|\theta_{kj}|<0.15} although \eqn{\alpha_j\neq0} (NA if none).
#' @param alpha_hat estimated global effects (p).
#' @param eps_hat estimated deviations (K x p) or NULL for common-effect methods.
#' @param truth output of \code{\link{sim_ipd}}.
#' @return one-row data frame.
#' @export
hmc_metrics <- function(alpha_hat, eps_hat, truth) {
  p <- length(truth$alpha); K <- nrow(truth$eps)
  act <- truth$alpha != 0; sel <- alpha_hat != 0
  TP <- sum(sel & act); FP <- sum(sel & !act); FN <- sum(!sel & act); TN <- sum(!sel & !act)
  mccd <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  if (is.null(eps_hat)) eps_hat <- matrix(0, K, p)
  theta_hat <- sweep(eps_hat, 2, alpha_hat, "+")
  te <- truth$eps != 0; se <- eps_hat != 0
  nc <- te & abs(truth$theta) < 0.15
  data.frame(
    TPR = TP / max(sum(act), 1), FDR = if (sum(sel)) FP / sum(sel) else 0,
    F1 = 2 * TP / max(2 * TP + FP + FN, 1),
    MCC = if (mccd > 0) (TP * TN - FP * FN) / mccd else 0,
    exact = as.integer(all(sel == act)), size = sum(sel),
    alpha_bias = mean(abs(alpha_hat[act] - truth$alpha[act])),
    fp_bias = if (sum(!act)) mean(abs(alpha_hat[!act])) else 0,
    MSE_alpha = mean((alpha_hat - truth$alpha)^2),
    MSE_eps = mean((eps_hat - truth$eps)^2),
    MSE_theta = mean((theta_hat - truth$theta)^2),
    TPR_dev = if (sum(te)) sum(se & te) / sum(te) else NA_real_,
    FDR_dev = if (sum(se)) sum(se & !te) / sum(se) else 0,
    n_dev = sum(se),
    n_near_cancel = sum(nc),
    err_theta_nc = if (sum(nc)) mean(abs(theta_hat[nc] - truth$theta[nc])) else NA_real_,
    hier_violation = sum(colSums(se) > 0 & !sel))
}
