#' Study-stratified fold assignment
#'
#' Folds are assigned separately within each study and, within study,
#' separately for events and censored observations, so that every training
#' set contains every study and a similar proportion of events.
#' @param study study labels.
#' @param status event indicator (1 = event).
#' @param nfolds number of folds.
#' @param seed optional seed.
#' @return integer vector of fold labels.
#' @export
hmc_folds <- function(study, status, nfolds = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  f <- integer(length(study))
  for (s in unique(study)) for (e in c(0, 1)) {
    idx <- which(study == s & status == e)
    if (length(idx) == 0) next
    start <- sample.int(nfolds, 1)
    f[idx[sample.int(length(idx))]] <- ((seq_along(idx) + start - 2) %% nfolds) + 1
  }
  f
}

#' Joint cross-validation of the two tuning parameters
#'
#' Selects \eqn{(\lambda_\alpha,\lambda_\varepsilon)} over the two-dimensional grid
#' formed by a path of \code{nlambda} values of \eqn{\lambda_\alpha} and the ratios
#' \eqn{r=\lambda_\varepsilon/\lambda_\alpha} in \code{ratio}. The criterion is the
#' cross-validated linear-predictor deviance of the stratified partial likelihood
#' (Dai and Breheny, 2019): out-of-fold linear predictors are assembled for every
#' subject and the stratified partial likelihood deviance is evaluated on the
#' full data. Folds are stratified by study and event status.
#'
#' @inheritParams hmcox
#' @param ratio vector of ratios \eqn{\lambda_\varepsilon/\lambda_\alpha} (default 0.5, 1, 2).
#' @param nfolds number of folds V (default 5).
#' @param foldid optional fold labels (length N, original order).
#' @param seed optional seed for the fold assignment.
#' @param rule "min" (default) or "1se" (smallest model within one fold-level
#'   standard error of the minimum).
#' @return An object of class \code{c("cv.hmcox","hmcox")}: the full-data fit at
#'   the selected ratio plus \code{cvm} (nlambda x nratio matrix of CV deviance),
#'   \code{lambda.min}, \code{ratio.min} and \code{index.min}.
#' @export
cv.hmcox <- function(X, time, status, study, penalty = c("MCP", "SCAD", "lasso"),
                     gamma = switch(penalty, SCAD = 3.7, 3),
                     ratio = c(0.5, 1, 2), nfolds = 5, foldid = NULL, seed = NULL,
                     nlambda = 25, lambda.min.ratio = NULL, heterogeneity = TRUE,
                     standardize = TRUE, tol = 1e-5, max.outer = 50,
                     max.inner = 2000, dfmax = NULL, rule = c("min", "1se")) {
  penalty <- match.arg(penalty); rule <- match.arg(rule)
  if (!heterogeneity) ratio <- 1
  D <- hmc_prepare(X, time, status, study, standardize)
  if (is.null(dfmax)) dfmax <- min(D$p, max(5, floor(sum(D$status) / 3)))
  lambda <- hmc_lambda_seq(D, nlambda, lambda.min.ratio)
  L <- length(lambda); R <- length(ratio)
  if (is.null(foldid)) foldid <- hmc_folds(study, status, nfolds, seed)
  fs <- foldid[D$ord]                                   # sorted order
  V <- length(unique(fs))

  full <- vector("list", R)
  for (r in seq_len(R))
    full[[r]] <- hmc_fit_prepared(D, penalty, gamma, ratio[r], lambda,
                                  heterogeneity, tol, max.outer, max.inner, dfmax)

  eta_cv <- array(NA_real_, dim = c(D$N, L, R))
  for (v in sort(unique(fs))) {
    tr <- fs != v
    Dtr <- hmc_subset(D, tr)
    te <- which(!tr)
    Xte <- D$X[te, , drop = FALSE]; ste <- D$sid[te]
    for (r in seq_len(R)) {
      f <- hmc_fit_prepared(Dtr, penalty, gamma, ratio[r], lambda,
                            heterogeneity, tol, max.outer, max.inner, dfmax)
      for (l in seq_along(f$lambda)) {
        th <- sweep(f$eps[, , l, drop = TRUE], 2, f$alpha[, l], "+")
        if (is.null(dim(th))) th <- matrix(th, nrow = D$K)
        eta_cv[te, l, r] <- rowSums(Xte * th[ste, , drop = FALSE])
      }
    }
  }
  cvm <- matrix(NA_real_, L, R, dimnames = list(NULL, paste0("ratio=", ratio)))
  cvsd <- cvm
  fl <- sort(unique(fs))
  Dfold <- lapply(fl, function(v) hmc_subset(D, fs == v))
  for (r in seq_len(R)) for (l in seq_len(L)) {
    if (anyNA(eta_cv[, l, r])) next
    cvm[l, r] <- 2 * D$N * hmc_loss(D, eta_cv[, l, r])
    ## fold-level contributions for a 1-SE rule (basic approach approximation)
    dv <- vapply(seq_along(fl), function(i)
      2 * Dfold[[i]]$N * hmc_loss(Dfold[[i]], eta_cv[fs == fl[i], l, r]), 0)
    cvsd[l, r] <- sd(dv) * sqrt(V)
  }
  dfm <- matrix(NA_real_, L, R)
  for (r in seq_len(R)) {
    nf <- length(full[[r]]$lambda)
    dfm[seq_len(nf), r] <- colSums(full[[r]]$alpha != 0) + apply(full[[r]]$eps != 0, 3, sum)
  }
  sel <- hmc_rule_index(cvm, cvsd, dfm)
  l.min <- sel$min[1]; r.min <- sel$min[2]
  if (rule == "1se") { l.min <- sel$se1[1]; r.min <- sel$se1[2] }
  fit <- full[[r.min]]
  structure(c(fit, list(penalty = penalty, gamma = gamma, ratio = ratio[r.min],
                        heterogeneity = heterogeneity, data = D,
                        cvm = cvm, cvsd = cvsd, lambda.grid = lambda,
                        ratio.grid = ratio, lambda.min = lambda[l.min],
                        ratio.min = ratio[r.min], index.min = unname(l.min),
                        sel.min = sel$min, sel.1se = sel$se1, df = dfm,
                        foldid = foldid, nfolds = V, rule = rule,
                        all.fits = full, call = match.call())),
            class = c("cv.hmcox", "hmcox"))
}

## indices (lambda, ratio) chosen by the minimum rule and by the one-standard-error
## rule: among grid points whose CV deviance is within one fold-level standard
## error of the minimum, the model with the fewest nonzero parameters
## (global effects + deviations); ties broken by larger lambda_alpha, then
## larger ratio.
hmc_rule_index <- function(cvm, cvsd, df) {
  best <- which(cvm == min(cvm, na.rm = TRUE), arr.ind = TRUE)[1, ]
  thr <- cvm[best[1], best[2]] + cvsd[best[1], best[2]]
  ok <- which(!is.na(cvm) & cvm <= thr, arr.ind = TRUE)
  o <- order(df[ok], ok[, 1], -ok[, 2])
  list(min = unname(best), se1 = unname(ok[o[1], ]))
}

#' Coefficients of a cross-validated fit
#' @param object a \code{cv.hmcox} object.
#' @param rule "min" or "1se"; defaults to the rule used when fitting.
#' @param ... unused.
#' @export
coef.cv.hmcox <- function(object, rule = object$rule, ...) {
  ix <- if (rule == "1se") object$sel.1se else object$sel.min
  f <- object$all.fits[[ix[2]]]
  a <- f$alpha[, ix[1]]; e <- f$eps[, , ix[1]]
  if (is.null(dim(e))) e <- matrix(e, nrow = object$data$K,
                                   dimnames = list(object$data$levels, object$data$xnames))
  list(alpha = a, eps = e, theta = sweep(e, 2, a, "+"),
       lambda.alpha = object$lambda.grid[ix[1]], ratio = object$ratio.grid[ix[2]])
}

#' @export
predict.cv.hmcox <- function(object, newX, study, rule = object$rule, global.only = FALSE, ...) {
  cf <- coef.cv.hmcox(object, rule)
  obj <- object
  obj$alpha <- matrix(cf$alpha, ncol = 1); obj$eps <- array(cf$eps, dim = c(dim(cf$eps), 1))
  predict.hmcox(obj, newX, study, index = 1, global.only = global.only)
}

#' @export
print.cv.hmcox <- function(x, ...) {
  cf <- coef(x)
  cat("Cross-validated hierarchical penalized stratified Cox model\n")
  cat("penalty =", x$penalty, "| folds =", x$nfolds, "| rule =", x$rule, "\n")
  cat(sprintf("lambda_alpha = %.4g, lambda_eps/lambda_alpha = %g (lambda_eps = %.4g)\n",
              x$lambda.min, x$ratio.min, x$lambda.min * x$ratio.min))
  cat("global effects selected:", sum(cf$alpha != 0), "of", x$data$p,
      "| nonzero deviations:", sum(cf$eps != 0), "\n")
  invisible(x)
}
