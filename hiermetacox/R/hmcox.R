#' @useDynLib hiermetacox, .registration = TRUE
#' @importFrom stats coef sd quantile qnorm pnorm uniroot rnorm runif rexp median setNames
NULL

## ---------------------------------------------------------------------------
## Internal: sort, standardise within study, build tie / stratum structure
## ---------------------------------------------------------------------------
hmc_prepare <- function(X, time, status, study, standardize = TRUE,
                        center = NULL, scale = NULL) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  n <- nrow(X)
  if (length(time) != n || length(status) != n || length(study) != n)
    stop("X, time, status and study must have the same number of rows.")
  if (any(!is.finite(time)) || any(time <= 0)) stop("time must be positive and finite.")
  if (!all(status %in% c(0, 1))) stop("status must be 0/1.")
  if (anyNA(X)) stop("X contains missing values; impute before fitting.")
  study <- as.factor(study)
  lev <- levels(droplevels(study))
  sid <- as.integer(factor(study, levels = lev))
  K <- length(lev)
  p <- ncol(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(p))

  if (standardize) {
    if (is.null(center)) {
      center <- matrix(0, K, p); scale <- matrix(1, K, p)
      for (k in seq_len(K)) {
        rows <- sid == k
        mu <- colMeans(X[rows, , drop = FALSE])
        s <- sqrt(colMeans(sweep(X[rows, , drop = FALSE], 2, mu)^2))
        s[!is.finite(s) | s < 1e-12] <- Inf          # constant column -> 0
        center[k, ] <- mu; scale[k, ] <- s
      }
    }
    for (k in seq_len(K)) {
      rows <- sid == k
      X[rows, ] <- sweep(sweep(X[rows, , drop = FALSE], 2, center[k, ]), 2, scale[k, ], "/")
    }
  }
  ord <- order(sid, time, -status)
  Xs <- X[ord, , drop = FALSE]
  ts <- as.double(time[ord]); ds <- as.double(status[ord]); ss <- sid[ord]
  nk <- tabulate(ss, K)
  if (any(nk == 0)) stop("empty study level.")
  sstart <- c(0L, cumsum(nk))
  ## tie groups (same study and same time), 0-based indices
  grp <- paste(ss, ts)
  rl <- rle(grp)
  ends <- cumsum(rl$lengths); starts <- ends - rl$lengths + 1
  tfirst <- rep(starts, rl$lengths) - 1L
  tlast <- rep(ends, rl$lengths) - 1L
  structure(list(X = Xs, time = ts, status = ds, sid = ss, ord = ord,
                 sstart = as.integer(sstart), tfirst = as.integer(tfirst),
                 tlast = as.integer(tlast), omega = sqrt(nk / sum(nk)),
                 nk = nk, K = K, p = p, N = n, levels = lev,
                 center = center, scale = scale, standardize = standardize,
                 xnames = colnames(X)),
            class = "hmc_data")
}

hmc_subset <- function(D, rows_sorted) {
  ## rows_sorted: logical/indices in the *sorted* order of D
  idx <- if (is.logical(rows_sorted)) which(rows_sorted) else sort(rows_sorted)
  X <- D$X[idx, , drop = FALSE]
  ss <- D$sid[idx]; ts <- D$time[idx]; ds <- D$status[idx]
  nk <- tabulate(ss, D$K)
  if (any(nk == 0)) stop("a training fold contains no subject from some study.")
  sstart <- c(0L, cumsum(nk))
  grp <- paste(ss, ts); rl <- rle(grp)
  ends <- cumsum(rl$lengths); starts <- ends - rl$lengths + 1
  out <- D
  out$X <- X; out$time <- ts; out$status <- ds; out$sid <- ss
  out$sstart <- as.integer(sstart)
  out$tfirst <- as.integer(rep(starts, rl$lengths) - 1L)
  out$tlast <- as.integer(rep(ends, rl$lengths) - 1L)
  out$omega <- sqrt(nk / sum(nk)); out$nk <- nk; out$N <- length(idx)
  out$ord <- NULL
  out
}

## stratified negative log partial likelihood / N (Breslow) and score
hmc_loss <- function(D, eta, score = FALSE) {
  res <- .Call("hmc_loss_score", if (score) D$X else double(0), D$status,
               D$sstart, D$tfirst, D$tlast, as.double(eta), as.integer(score),
               PACKAGE = "hiermetacox")
  if (score) list(loss = res[[1]], score = res[[2]]) else res[[1]]
}

## per-observation contributions to the stratified partial-likelihood deviance,
## -2 * delta_i * (eta_i - log sum_{l in R_k(Y_i)} exp(eta_l)), in the sorted order
## of D (Breslow ties: the risk set of a tied time includes all tied subjects).
## Their sum is 2 * N * hmc_loss(D, eta).
hmc_dev_terms <- function(D, eta) {
  out <- numeric(D$N)
  for (k in seq_len(D$K)) {
    rows <- (D$sstart[k] + 1):D$sstart[k + 1]
    e <- eta[rows]; m <- max(e); hz <- exp(e - m)
    R <- rev(cumsum(rev(hz)))
    rsk <- R[D$tfirst[rows] - D$sstart[k] + 1]
    out[rows] <- -2 * D$status[rows] * (e - m - log(rsk))
  }
  out
}

pen_code <- function(penalty) switch(penalty, lasso = 0L, MCP = 1L, SCAD = 2L, enet = 3L,
                                     stop("penalty must be MCP, SCAD, lasso or enet"))

#' Lambda sequence for the global effects
#' @keywords internal
hmc_lambda_seq <- function(D, nlambda = 25, lambda.min.ratio = NULL, l1 = 1) {
  sc <- hmc_loss(D, rep(0, D$N), score = TRUE)$score
  lmax <- max(abs(sc)) / l1 * 1.0001   # l1 = elastic-net mixing (1 otherwise)
  if (is.null(lambda.min.ratio)) lambda.min.ratio <- if (D$N > D$p) 0.05 else 0.10
  exp(seq(log(lmax), log(lmax * lambda.min.ratio), length.out = nlambda))
}

hmc_fit_prepared <- function(D, penalty, gamma, ratio, lambda, heterogeneity,
                             tol, max.outer, max.inner, dfmax,
                             alpha.init = NULL, eps.init = NULL) {
  p <- D$p; K <- D$K
  a0 <- if (is.null(alpha.init)) rep(0, p) else alpha.init
  e0 <- if (is.null(eps.init)) rep(0, K * p) else as.double(eps.init)
  res <- .Call("hmc_path", D$X, D$status, D$sstart, D$tfirst, D$tlast,
               as.double(D$omega), pen_code(penalty), as.double(gamma),
               as.double(lambda), as.double(ratio), as.integer(heterogeneity),
               as.double(tol), as.integer(max.outer), as.integer(max.inner),
               as.integer(dfmax), as.double(a0), as.double(e0),
               PACKAGE = "hiermetacox")
  L <- length(lambda)
  alpha <- res[[1]]
  dimnames(alpha) <- list(D$xnames, NULL)
  eps <- array(res[[2]], dim = c(K, p, L), dimnames = list(D$levels, D$xnames, NULL))
  keep <- !is.na(alpha[1, ])
  list(alpha = alpha[, keep, drop = FALSE], eps = eps[, , keep, drop = FALSE],
       lambda = lambda[keep], loss = res[[3]][keep], objective = res[[4]][keep],
       iter = res[[5]][keep], converged = as.logical(res[[6]][keep]),
       halvings = res[[7]][keep])
}

#' Fit a hierarchical penalized stratified Cox model for IPD meta-analysis
#'
#' Fits the model \eqn{h_{ki}(t)=h_{0k}(t)\exp\{x_{ki}^\top(\alpha+\varepsilon_k)\}}
#' along a path of \eqn{\lambda_\alpha} values, with
#' \eqn{\lambda_\varepsilon = \code{ratio}\times\lambda_\alpha}. The objective is
#' \deqn{-\frac{1}{N}\sum_k \ell_k(\alpha+\varepsilon_k) + \sum_j P(|\alpha_j|;\lambda_\alpha)
#'  + \sum_{k,j} P(\omega_k|\varepsilon_{kj}|;\lambda_\varepsilon)}
#' subject to the hierarchy constraint \eqn{\alpha_j=0 \Rightarrow \varepsilon_{kj}=0},
#' where \eqn{\ell_k} is the Breslow partial log-likelihood of study \eqn{k}
#' (each study has its own baseline hazard) and \eqn{\omega_k=(n_k/N)^{1/2}}.
#'
#' @param X numeric covariate matrix (N x p).
#' @param time,status follow-up time and event indicator (1 = event).
#' @param study study (stratum) membership, length N.
#' @param penalty "MCP" (default), "SCAD", "lasso" or "enet" (elastic net,
#'   penalty \eqn{\lambda\{a|b|+(1-a)b^2/2\}} as in glmnet).
#' @param gamma concavity parameter (default 3 for MCP, 3.7 for SCAD); for
#'   \code{penalty = "enet"} it is the mixing parameter \eqn{a} (default 0.5).
#' @param ratio \eqn{\lambda_\varepsilon/\lambda_\alpha}.
#' @param lambda optional decreasing sequence of \eqn{\lambda_\alpha}.
#' @param nlambda,lambda.min.ratio path length and smallest/largest ratio.
#' @param heterogeneity if FALSE, all deviations are fixed at zero
#'   (penalized stratified Cox model with common coefficients).
#' @param standardize standardise each covariate within each study (default TRUE);
#'   coefficients are then per within-study standard deviation.
#' @param tol,max.outer,max.inner convergence controls.
#' @param dfmax stop the path once more than dfmax global effects are active.
#' @return An object of class \code{"hmcox"}.
#' @export
hmcox <- function(X, time, status, study, penalty = c("MCP", "SCAD", "lasso", "enet"),
                  gamma = switch(penalty, SCAD = 3.7, enet = 0.5, 3), ratio = 1,
                  lambda = NULL, nlambda = 25, lambda.min.ratio = NULL,
                  heterogeneity = TRUE, standardize = TRUE, tol = 1e-5,
                  max.outer = 50, max.inner = 2000, dfmax = NULL) {
  penalty <- match.arg(penalty)
  D <- hmc_prepare(X, time, status, study, standardize)
  if (is.null(dfmax)) dfmax <- min(D$p, max(5, floor(sum(D$status) / 3)))
  if (penalty == "enet" && !(gamma > 0 && gamma <= 1)) stop("for enet, gamma is the mixing parameter in (0, 1]")
  if (is.null(lambda)) lambda <- hmc_lambda_seq(D, nlambda, lambda.min.ratio,
                                                if (penalty == "enet") gamma else 1)
  fit <- hmc_fit_prepared(D, penalty, gamma, ratio, lambda, heterogeneity,
                          tol, max.outer, max.inner, dfmax)
  structure(c(fit, list(penalty = penalty, gamma = gamma, ratio = ratio,
                        heterogeneity = heterogeneity, data = D,
                        call = match.call())),
            class = "hmcox")
}

#' @export
print.hmcox <- function(x, ...) {
  cat("Hierarchical penalized stratified Cox model (", x$penalty, ", gamma = ",
      x$gamma, ", lambda_eps/lambda_alpha = ", x$ratio, ")\n", sep = "")
  cat("K =", x$data$K, "studies | N =", x$data$N, "| p =", x$data$p,
      "| events =", sum(x$data$status), "\n")
  df <- colSums(x$alpha != 0)
  dfe <- apply(x$eps != 0, 3, sum)
  print(data.frame(lambda = signif(x$lambda, 4), n_global = df, n_deviation = dfe,
                   converged = x$converged), row.names = FALSE)
  invisible(x)
}

#' Coefficients of a fitted hierarchical model
#' @param object an \code{hmcox} or \code{cv.hmcox} object.
#' @param index path index (defaults to the CV-selected index for cv objects).
#' @param ... unused.
#' @return list with \code{alpha} (p), \code{eps} (K x p) and \code{theta} (K x p).
#' @export
coef.hmcox <- function(object, index = length(object$lambda), ...) {
  a <- object$alpha[, index]
  e <- object$eps[, , index]
  if (is.null(dim(e))) e <- matrix(e, nrow = object$data$K,
                                   dimnames = list(object$data$levels, object$data$xnames))
  list(alpha = a, eps = e, theta = sweep(e, 2, a, "+"))
}

#' Linear predictor for new data
#' @param object an \code{hmcox} or \code{cv.hmcox} object.
#' @param newX new covariates (raw scale).
#' @param study study labels of the new rows; rows from studies not seen in
#'   training use the global effects only.
#' @param index path index.
#' @param global.only use \eqn{\alpha} only, ignoring deviations.
#' @param ... unused.
#' @export
predict.hmcox <- function(object, newX, study, index = length(object$lambda),
                          global.only = FALSE, ...) {
  D <- object$data
  newX <- as.matrix(newX)
  cf <- coef.hmcox(object, index)
  lp <- numeric(nrow(newX))
  st <- as.character(study)
  for (s in unique(st)) {
    rows <- st == s
    Z <- newX[rows, , drop = FALSE]
    k <- match(s, D$levels)
    if (D$standardize) {
      if (!is.na(k)) {
        Z <- sweep(sweep(Z, 2, D$center[k, ]), 2, D$scale[k, ], "/")
      } else {                           # unseen study: standardise within itself
        mu <- colMeans(Z); sdv <- sqrt(colMeans(sweep(Z, 2, mu)^2)); sdv[sdv < 1e-12] <- Inf
        Z <- sweep(sweep(Z, 2, mu), 2, sdv, "/")
      }
    }
    b <- if (is.na(k) || global.only) cf$alpha else cf$theta[k, ]
    lp[rows] <- drop(Z %*% b)
  }
  lp
}
