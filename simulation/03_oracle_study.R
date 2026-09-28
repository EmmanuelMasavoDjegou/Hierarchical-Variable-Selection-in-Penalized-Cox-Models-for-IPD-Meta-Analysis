################################################################################
##  simulation/03_oracle_study.R
##  Empirical verification of the oracle results (Section 5.5): Theorem 4.1
##  (support recovery and exact equality with the oracle), Theorem 4.3
##  (asymptotic normality and efficiency) and Corollary 4.1 (Wald coverage),
##  for the proposed hierarchical MCP and SCAD estimators, with the hierarchical
##  lasso as a negative control (not oracle by theory).
##
##  Design  : K = 5, p = 100, fixed truth satisfying Assumption 1 and condition
##            (a) of Theorem 4.2; study sizes n_k in {200, 400, 800, 1600}.
##  Usage   : Rscript 03_oracle_study.R [R] ;  NCORES=8 Rscript 03_oracle_study.R 200
##  Output  : results/oracle/nNNNN_rYYY.csv   one file per (size, replicate)
##  Seeds   : replicate r at study size n uses set.seed(20270000 + 10 * n + r).
##  Next    : 04_oracle_tables_figures.R builds Table 4 and Figures 5-6.
################################################################################
args <- commandArgs(trailingOnly = TRUE)
here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
source(file.path(here, "00_setup.R"))
Rrep <- if (length(args) >= 1) as.integer(args[1]) else 200L
SIZES <- if (length(args) >= 2) as.integer(args[-1]) else c(200L, 400L, 800L, 1600L)
ncores <- max(1L, as.integer(Sys.getenv("NCORES", "1")))
odir <- file.path(here, "results", "oracle"); dir.create(odir, showWarnings = FALSE, recursive = TRUE)

## ---- fixed truth ------------------------------------------------------------------
ORACLE <- local({
  K <- 5; p <- 100
  alpha <- c(1, -0.8, 0.6, -0.5, 0.4, -0.3)
  eps <- matrix(0, K, 6)
  ## covariate j = 1..4 deviates in studies j and j+1, moving away from zero, so
  ## that every theta_kj is well separated from zero and from the shared value
  for (j in 1:4) {
    ks <- c(j, j %% K + 1)
    eps[ks, j] <- sign(alpha[j]) * c(0.5, 1.0)
  }
  list(K = K, p = p, alpha = alpha, eps = eps, scales = c(0.05, 0.075, 0.1, 0.125, 0.15),
       rho = 0.5, cens = 0.3)
})
ORACLE_METHODS <- c("MCP-H", "SCAD-H", "LASSO-H")
## Tuning strategies: "theory" uses the rate of Assumption 7,
## lambda_alpha = lambda_eps = KAPPA * sqrt(L_N / N) with L_N = log(max(p, N));
## "min" and "1se" are the cross-validation rules of Section 3.2.
KAPPA <- 1.1

## reduced parameter of the true support: alpha_j (j in A0) and, for the deviating
## cells, the study-specific coefficient theta_kj = alpha_j + eps_kj
true_params <- function() {
  A0 <- which(ORACLE$alpha != 0)
  D0 <- which(ORACLE$eps != 0, arr.ind = TRUE)
  rbind(data.frame(param = paste0("alpha_", A0), type = "alpha", k = NA, j = A0,
                   truth = ORACLE$alpha[A0]),
        data.frame(param = sprintf("theta_%d_%d", D0[, 1], D0[, 2]), type = "theta",
                   k = D0[, 1], j = D0[, 2],
                   truth = ORACLE$alpha[D0[, 2]] + ORACLE$eps[D0]))
}

## unpenalised stratified Cox fit (Breslow) on a given support; returns estimates
## and SEs of alpha_j and of theta_kj for the requested parameters
fit_support <- function(X, time, status, study, A, D, par) {
  cols <- list(); lab <- character(0)
  for (j in A) { cols[[length(cols) + 1]] <- X[, j]; lab <- c(lab, paste0("a", j)) }
  if (length(D)) for (r in seq_len(nrow(D))) {
    cols[[length(cols) + 1]] <- X[, D[r, 2]] * (study == D[r, 1])
    lab <- c(lab, sprintf("e%d_%d", D[r, 1], D[r, 2]))
  }
  Z <- do.call(cbind, cols); colnames(Z) <- lab
  df <- data.frame(time = time, status = status, sid = study)
  f <- coxph(Surv(time, status) ~ Z + strata(sid), data = df, ties = "breslow")
  b <- coef(f); V <- vcov(f); names(b) <- lab; dimnames(V) <- list(lab, lab)
  out <- t(sapply(seq_len(nrow(par)), function(i) {
    w <- setNames(numeric(length(b)), lab)
    ja <- paste0("a", par$j[i])
    if (!ja %in% lab) return(c(NA, NA))
    w[ja] <- 1
    if (par$type[i] == "theta") {
      je <- sprintf("e%d_%d", par$k[i], par$j[i]); if (je %in% lab) w[je] <- 1
    }
    c(sum(w * b), sqrt(drop(t(w) %*% V %*% w)))
  }))
  data.frame(est = out[, 1], se = out[, 2])
}

run_one <- function(n, r) {
  f_out <- file.path(odir, sprintf("n%04d_r%03d.csv", n, r))
  if (file.exists(f_out)) return(invisible(NULL))
  set.seed(20270000 + 10 * n + r)
  sim <- sim_ipd(K = ORACLE$K, n = n, p = ORACLE$p, alpha = ORACLE$alpha, rho = ORACLE$rho,
                 cens = ORACLE$cens, eps = ORACLE$eps, scales = ORACLE$scales)
  d <- sim$data; X <- as.matrix(d[, -(1:3)])       # unit-variance design, not re-standardised
  par <- true_params()
  A0 <- which(ORACLE$alpha != 0); D0 <- unname(which(sim$eps != 0, arr.ind = TRUE))
  orc <- fit_support(X, d$time, d$status, d$study, A0, D0, par)      # the oracle
  foldid <- hmc_folds(d$study, d$status, TUNE$nfolds)
  rows <- list()
  N <- nrow(d); lam_th <- KAPPA * sqrt(log(max(ORACLE$p, N)) / N)
  for (m in ORACLE_METHODS) {
    pen <- switch(sub("-H$", "", m), MCP = "MCP", SCAD = "SCAD", LASSO = "lasso")
    gam <- switch(pen, SCAD = 3.7, 3)
    cvf <- cv.hmcox(X, d$time, d$status, d$study, penalty = pen,
                    gamma = gam, ratio = TUNE$ratio, foldid = foldid,
                    nlambda = TUNE$nlambda, lambda.min.ratio = TUNE$lambda.min.ratio,
                    tol = 1e-7, dfmax = TUNE$dfmax, standardize = FALSE)
    fth <- hmcox(X, d$time, d$status, d$study, penalty = pen, gamma = gam, ratio = 1,
                 lambda = exp(seq(log(0.5), log(lam_th), length.out = 20)),
                 tol = 1e-7, standardize = FALSE)
    for (rule in c("theory", "min", "1se")) {
      cf <- if (rule == "theory") coef(fth, length(fth$lambda)) else coef(cvf, rule = rule)
      Ahat <- unname(which(cf$alpha != 0)); Dhat <- unname(which(cf$eps != 0, arr.ind = TRUE))
      correct <- setequal(Ahat, A0) && setequal(paste(Dhat[, 1], Dhat[, 2]),
                                                paste(D0[, 1], D0[, 2]))
      pen_est <- unname(ifelse(par$type == "alpha", cf$alpha[par$j],
                               cf$theta[cbind(pmax(par$k, 1), par$j)]))
      ref <- if (length(Ahat)) fit_support(X, d$time, d$status, d$study, Ahat,
                                           if (nrow(Dhat)) Dhat else NULL, par) else
        data.frame(est = rep(NA, nrow(par)), se = rep(NA, nrow(par)))
      rows[[paste(m, rule)]] <- data.frame(
        n = n, rep = r, method = m, rule = rule, par,
        correct_support = as.integer(correct),
        pen = pen_est, oracle = orc$est, oracle_se = orc$se,
        refit = ref$est, refit_se = ref$se,
        equal_oracle = as.integer(correct && max(abs(pen_est - orc$est)) < 1e-3),
        covered = as.integer(!is.na(ref$se) &
                               abs(ref$est - par$truth) <= qnorm(0.975) * ref$se))
    }
  }
  tmp <- paste0(f_out, ".tmp")
  write.csv(do.call(rbind, rows), tmp, row.names = FALSE); file.rename(tmp, f_out)
  invisible(NULL)
}

jobs <- expand.grid(r = seq_len(Rrep), n = SIZES)
jobs <- jobs[!file.exists(file.path(odir, sprintf("n%04d_r%03d.csv", jobs$n, jobs$r))), ]
cat(sprintf("[%s] %d oracle replicate(s) on %d core(s)\n", format(Sys.time(), "%H:%M:%S"),
            nrow(jobs), ncores))
if (nrow(jobs)) {
  runner <- function(i) run_one(jobs$n[i], jobs$r[i])
  if (ncores > 1 && .Platform$OS.type == "unix") {
    invisible(parallel::mclapply(seq_len(nrow(jobs)), runner, mc.cores = ncores,
                                 mc.preschedule = FALSE))
  } else for (i in seq_len(nrow(jobs))) runner(i)
}
cat("done\n")
