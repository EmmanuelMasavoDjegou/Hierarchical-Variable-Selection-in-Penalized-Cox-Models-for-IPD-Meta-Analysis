################################################################################
##  simulation/05_near_cancellation.R
##  Near-cancellation experiment (Section 5.5): tests condition (b) of Theorem 4.2
##  (representation consistency), including the multiplicity q(c) of alternative
##  coefficient values.
##
##  Design  : K = 5, p = 100, n_k = 800. Covariate 3 (alpha = 1.2) deviates in
##            two studies with theta = 0.01 (near zero, q = 1, delta = m - q = 2)
##            and theta = 2.4; covariate 4 (alpha = -1.2) deviates in two studies
##            sharing theta = -0.01 (q = 2, delta = 1). All deviations are at least
##            1.19 in absolute value, inside the flat region for every ratio used.
##            Condition (b), F(lambda_a) < delta F(lambda_e), requires
##            r = lambda_e / lambda_a > 1/sqrt(2) = 0.71 for covariate 3 and r > 1
##            for covariate 4 (MCP and SCAD: F proportional to lambda^2).
##  Tuning  : lambda_a = 1.1 sqrt(L_N / N) (Remark 4.2), r in {0.5, 0.85, 2}:
##            r = 0.85 separates the two thresholds.
##  Records : theorem_majority - Theorem 4.2 checked exactly: among all feasible
##            representations of the ORACLE fit, is the majority one cheapest?
##            majority - does the ALGORITHM's output use the majority representation?
##            alg_worse - is the algorithm's objective above that of the oracle fit
##            in its majority representation (an inferior local minimum)?
##  Usage   : NCORES=8 Rscript 05_near_cancellation.R [R]      (default R = 100)
##  Output  : results/nearcancel/rYYY.csv, one file per replicate
##  Seeds   : replicate r uses set.seed(20280000 + r).
################################################################################
args <- commandArgs(trailingOnly = TRUE)
here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
source(file.path(here, "00_setup.R"))
Rrep <- if (length(args) >= 1) as.integer(args[1]) else 100L
ncores <- max(1L, as.integer(Sys.getenv("NCORES", "1")))
ndir <- file.path(here, "results", "nearcancel"); dir.create(ndir, showWarnings = FALSE, recursive = TRUE)

NC <- local({
  K <- 5; alpha <- c(1, -0.8, 1.2, -1.2, 0.4, -0.3)
  eps <- matrix(0, K, 6)
  eps[1:2, 1] <- c(0.5, 1.0)            # well separated, as in the oracle design
  eps[2:3, 2] <- c(-0.5, -1.0)
  eps[3:4, 3] <- c(-1.19, 1.2)          # theta = 0.01 (near zero, q = 1) and 2.4
  eps[4:5, 4] <- c(1.19, 1.19)          # theta = -0.01 in two studies (q = 2)
  list(K = K, p = 100, n = 800, alpha = alpha, eps = eps,
       scales = c(0.05, 0.075, 0.1, 0.125, 0.15), ratios = c(0.5, 0.85, 2), kappa = 1.1)
})

## penalty value, as in the solver (MCP, SCAD)
pen_val <- function(t, lam, pen, gam = if (pen == "SCAD") 3.7 else 3) {
  t <- abs(t)
  if (pen == "MCP") ifelse(t <= gam * lam, lam * t - t^2 / (2 * gam), gam * lam^2 / 2)
  else ifelse(t <= lam, lam * t, ifelse(t <= gam * lam,
              (2 * gam * lam * t - t^2 - lam^2) / (2 * (gam - 1)), lam^2 * (gam + 1) / 2))
}
## oracle fit on the true support, returned as a K x p matrix of theta
oracle_theta <- function(X, time, status, study, A, D) {
  cols <- lapply(A, function(j) X[, j]); lab <- paste0("a", A)
  for (r in seq_len(nrow(D))) { cols[[length(cols) + 1]] <- X[, D[r, 2]] * (study == D[r, 1])
    lab <- c(lab, sprintf("e%d_%d", D[r, 1], D[r, 2])) }
  Z <- do.call(cbind, cols)
  df <- data.frame(time = time, status = status, sid = study)
  b <- coef(coxph(Surv(time, status) ~ Z + strata(sid), data = df, ties = "breslow"))
  th <- matrix(0, NC$K, NC$p); th[, A] <- matrix(b[seq_along(A)], NC$K, length(A), byrow = TRUE)
  for (r in seq_len(nrow(D))) th[D[r, 1], D[r, 2]] <- th[D[r, 1], D[r, 2]] + b[length(A) + r]
  th
}

run_one <- function(r) {
  f_out <- file.path(ndir, sprintf("r%03d.csv", r))
  if (file.exists(f_out)) return(invisible(NULL))
  set.seed(20280000 + r)
  sim <- sim_ipd(K = NC$K, n = NC$n, p = NC$p, alpha = NC$alpha, eps = NC$eps,
                 scales = NC$scales, rho = 0.5, cens = 0.3)
  d <- sim$data; X <- as.matrix(d[, -(1:3)]); N <- nrow(d)
  lam_a <- NC$kappa * sqrt(log(max(NC$p, N)) / N)
  D0 <- unname(which(sim$eps != 0, arr.ind = TRUE))
  th_or <- oracle_theta(X, d$time, d$status, d$study, 1:6, D0)
  rows <- list()
  for (pen in c("MCP", "SCAD")) for (ratio in NC$ratios) {
    f <- hmcox(X, d$time, d$status, d$study, penalty = pen, ratio = ratio,
               lambda = exp(seq(log(0.5), log(lam_a), length.out = 20)),
               tol = 1e-7, standardize = FALSE)
    cf <- coef(f, length(f$lambda))
    Ahat <- unname(which(cf$alpha != 0)); Dhat <- unname(which(cf$eps != 0, arr.ind = TRUE))
    support_theta <- setequal(Ahat, 1:6)       # the fitted theta_.j have the true support
    lam_e <- ratio * lam_a; om <- f$data$omega
    ## objective at the oracle fit in its majority representation (alpha = shared value)
    a_or <- sapply(1:NC$p, function(j) { tb <- table(round(th_or[, j], 10))
      as.numeric(names(tb)[which.max(tb)]) })
    e_or <- sweep(th_or, 2, a_or)
    eta_or <- rowSums(f$data$X * th_or[f$data$sid, ])
    Q_or <- hiermetacox:::hmc_loss(f$data, eta_or) + sum(pen_val(a_or, lam_a, pen)) +
      sum(pen_val(om * e_or, lam_e, pen))
    Q_alg <- f$objective[length(f$lambda)]
    for (j in 3:4) {
      dj_true <- D0[D0[, 2] == j, 1]; dj_hat <- Dhat[Dhat[, 2] == j, 1]
      ## Theorem 4.2 check: among ALL feasible representations of the oracle fit
      ## (anchor c = one of its theta values, c != 0), is the majority one cheapest?
      cand <- unique(round(th_or[, j], 10)); cand <- cand[cand != 0]
      g <- sapply(cand, function(c) pen_val(c, lam_a, pen) + sum(pen_val(om * (th_or[, j] - c), lam_e, pen)))
      thm_majority <- as.integer(abs(cand[which.min(g)] - a_or[j]) < 1e-8)
      rows[[paste(pen, ratio, j)]] <- data.frame(
        rep = r, penalty = pen, ratio = ratio, covariate = j,
        q = if (j == 3) 1L else 2L, alpha_true = NC$alpha[j], alpha_hat = unname(cf$alpha[j]),
        n_dev_hat = length(dj_hat), theorem_majority = thm_majority,
        Q_alg = Q_alg, Q_oracle = Q_or, alg_worse = as.integer(Q_alg > Q_or + 1e-10),
        majority = as.integer(support_theta && setequal(dj_hat, dj_true) &&
                                abs(cf$alpha[j] - NC$alpha[j]) < 0.25))
    }
  }
  tmp <- paste0(f_out, ".tmp")
  write.csv(do.call(rbind, rows), tmp, row.names = FALSE); file.rename(tmp, f_out)
  invisible(NULL)
}

todo <- which(!file.exists(file.path(ndir, sprintf("r%03d.csv", seq_len(Rrep)))))
cat(sprintf("[%s] %d near-cancellation replicate(s) on %d core(s)\n",
            format(Sys.time(), "%H:%M:%S"), length(todo), ncores))
if (length(todo)) {
  if (ncores > 1 && .Platform$OS.type == "unix") {
    invisible(parallel::mclapply(todo, run_one, mc.cores = ncores, mc.preschedule = FALSE))
  } else for (r in todo) run_one(r)
}
cat("done\n")
