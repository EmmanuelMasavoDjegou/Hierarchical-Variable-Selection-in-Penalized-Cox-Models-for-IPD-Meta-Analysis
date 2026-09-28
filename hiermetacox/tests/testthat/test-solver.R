## Correctness tests for the compiled solver (Section 3 of the manuscript).

make_data <- function(seed = 1, K = 3, n = c(80, 120, 100), p = 5, ties = TRUE) {
  set.seed(seed)
  s <- sim_ipd(K = K, n = n, p = p, alpha = c(0.8, -0.5), het.idx = 1, n.dev = 1,
               sd.eps = 0.6, cens = 0.3)
  if (ties) s$data$time <- round(s$data$time, 1) + 0.1
  s
}

test_that("unpenalized common-effect fit equals coxph with strata (Breslow ties)", {
  s <- make_data(); d <- s$data; X <- as.matrix(d[, -(1:3)])
  f <- hmcox(X, d$time, d$status, d$study, penalty = "lasso", lambda = 1e-10,
             heterogeneity = FALSE, tol = 1e-10)
  D <- f$data
  library(survival)
  cx <- coxph(Surv(D$time, D$status) ~ D$X + strata(D$sid), ties = "breslow")
  expect_equal(unname(f$alpha[, 1]), unname(coef(cx)), tolerance = 1e-6)
  expect_equal(f$loss, -cx$loglik[2] / D$N, tolerance = 1e-8)
})

test_that("unpenalized heterogeneous fit equals separate per-study Cox fits", {
  s <- make_data(); d <- s$data; X <- as.matrix(d[, -(1:3)])
  f <- hmcox(X, d$time, d$status, d$study, penalty = "lasso", lambda = 1e-7,
             tol = 1e-11, max.outer = 500)
  D <- f$data; th <- coef(f)$theta
  sep <- t(sapply(1:3, function(k) {
    r <- D$sid == k
    coef(survival::coxph(survival::Surv(D$time[r], D$status[r]) ~ D$X[r, ], ties = "breslow"))
  }))
  expect_equal(unname(th), unname(sep), tolerance = 1e-4)
})

test_that("lasso path matches glmnet's stratified Cox model", {
  skip_if_not_installed("glmnet")
  s <- make_data(2, n = 150, p = 20, ties = FALSE); d <- s$data; X <- as.matrix(d[, -(1:3)])
  f <- hmcox(X, d$time, d$status, d$study, penalty = "lasso", heterogeneity = FALSE,
             nlambda = 10, tol = 1e-10)
  D <- f$data
  y <- glmnet::stratifySurv(survival::Surv(D$time, D$status), D$sid)
  g <- glmnet::glmnet(D$X, y, family = "cox", lambda = f$lambda, standardize = FALSE,
                      thresh = 1e-13)
  expect_lt(max(abs(as.matrix(coef(g)) - f$alpha)), 1e-4)
})

test_that("elastic-net path matches glmnet's stratified Cox model", {
  skip_if_not_installed("glmnet")
  s <- make_data(2, n = 150, p = 20, ties = FALSE); d <- s$data; X <- as.matrix(d[, -(1:3)])
  f <- hmcox(X, d$time, d$status, d$study, penalty = "enet", gamma = 0.5,
             heterogeneity = FALSE, nlambda = 10, tol = 1e-10)
  D <- f$data
  y <- glmnet::stratifySurv(survival::Surv(D$time, D$status), D$sid)
  g <- glmnet::glmnet(D$X, y, family = "cox", alpha = 0.5, lambda = f$lambda,
                      standardize = FALSE, thresh = 1e-13)
  expect_lt(max(abs(as.matrix(coef(g)) - f$alpha)), 1e-4)
})

test_that("MCP and SCAD solutions satisfy the KKT conditions for global effects", {
  s <- make_data(2, n = 150, p = 20, ties = FALSE); d <- s$data; X <- as.matrix(d[, -(1:3)])
  for (pen in c("MCP", "SCAD")) {
    f <- hmcox(X, d$time, d$status, d$study, penalty = pen, nlambda = 15, tol = 1e-10)
    D <- f$data; l <- 10; a <- f$alpha[, l]; th <- coef(f, l)$theta
    eta <- rowSums(D$X * th[D$sid, ])
    sc <- hiermetacox:::hmc_loss(D, eta, score = TRUE)$score
    lam <- f$lambda[l]; g0 <- f$gamma
    dP <- function(b) { b <- abs(b)
      if (pen == "MCP") pmax(lam - b / g0, 0) else ifelse(b <= lam, lam, pmax(g0 * lam - b, 0) / (g0 - 1)) }
    ## unconstrained coordinates: alpha_j != 0 with no deviations
    nz <- a != 0 & colSums(coef(f, l)$eps != 0) == 0
    expect_lt(max(abs(sc[nz] - dP(a[nz]) * sign(a[nz]))), 1e-6)
    expect_true(all(abs(sc[!nz]) <= lam + 1e-8))
  }
})

test_that("hierarchy holds, a sharing study exists, and the path converges", {
  set.seed(3)
  s <- sim_ipd(K = 5, n = 150, p = 30, sd.eps = 0.8)
  d <- s$data; X <- as.matrix(d[, -(1:3)])
  for (pen in c("MCP", "SCAD", "lasso", "enet")) {
    f <- hmcox(X, d$time, d$status, d$study, penalty = pen, ratio = 0.5, nlambda = 20)
    for (l in seq_along(f$lambda)) {
      e <- f$eps[, , l]; a <- f$alpha[, l]
      expect_true(all(colSums(e != 0)[a == 0] == 0))      # alpha_j = 0 => eps_.j = 0
      expect_true(all(colSums(e != 0)[a != 0] < nrow(e)))  # some study shares alpha_j
    }
    expect_true(all(f$converged))
  }
})

test_that("cross-validation returns valid selections under both rules", {
  set.seed(4)
  s <- sim_ipd(K = 4, n = 120, p = 20)
  d <- s$data; X <- as.matrix(d[, -(1:3)])
  cvf <- cv.hmcox(X, d$time, d$status, d$study, penalty = "MCP", nfolds = 3, seed = 1, nlambda = 12)
  expect_equal(dim(cvf$cvm), c(12L, 3L))
  df <- function(r) { cf <- coef(cvf, rule = r); sum(cf$alpha != 0) + sum(cf$eps != 0) }
  expect_lte(df("1se"), df("min"))
  inf <- hmcox_inference(cvf, rule = "1se")
  expect_true(all(inf$global$HR_lower <= inf$global$HR & inf$global$HR <= inf$global$HR_upper))
})

test_that("every training fold contains every study", {
  set.seed(5)
  st <- rep(1:3, c(20, 30, 25)); ev <- rbinom(75, 1, 0.6)
  f <- hmc_folds(st, ev, nfolds = 5, seed = 1)
  for (v in 1:5) expect_equal(length(unique(st[f != v])), 3)
})

test_that("a single stratum (naive pooling) works and equals the unstratified Cox lasso", {
  skip_if_not_installed("glmnet")
  set.seed(6)
  s <- sim_ipd(K = 3, n = 100, p = 15); d <- s$data; X <- as.matrix(d[, -(1:3)])
  one <- rep(1L, nrow(d))
  f <- hmcox(X, d$time, d$status, one, penalty = "lasso", heterogeneity = FALSE,
             nlambda = 8, tol = 1e-10)
  g <- glmnet::glmnet(f$data$X, survival::Surv(f$data$time, f$data$status), family = "cox",
                      lambda = f$lambda, standardize = FALSE, thresh = 1e-13)
  expect_lt(max(abs(as.matrix(coef(g)) - f$alpha)), 1e-4)
  cvf <- cv.hmcox(X, d$time, d$status, one, penalty = "MCP", heterogeneity = FALSE,
                  nfolds = 3, seed = 1, nlambda = 10)
  expect_equal(nrow(coef(cvf)$eps), 1L)
})

test_that("the CV criterion is exactly the sum of its fold-level scores", {
  set.seed(7)
  s <- sim_ipd(K = 4, n = 120, p = 20); d <- s$data; X <- as.matrix(d[, -(1:3)])
  f <- hmcox(X, d$time, d$status, d$study, nlambda = 5)
  D <- f$data; th <- coef(f, 3)$theta; eta <- rowSums(D$X * th[D$sid, ])
  expect_equal(sum(hiermetacox:::hmc_dev_terms(D, eta)),
               2 * D$N * hiermetacox:::hmc_loss(D, eta), tolerance = 1e-10)
})
