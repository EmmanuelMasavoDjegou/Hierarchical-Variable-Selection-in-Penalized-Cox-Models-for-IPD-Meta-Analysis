################################################################################
##  simulation/07_sensitivity.R
##  Sensitivity experiments for the revision (Sections 3.2 and 5.4), on exactly the
##  data sets and folds of the main study (same seeds as 01_run_simulation.R).
##
##  Part A  cross-validation criterion: hierarchical MCP and SCAD in all 13
##          scenarios; supports selected by the linear-predictor criterion (Dai and
##          Breheny, 2024; used in the paper) and by the fold-wise criterion of
##          Verweij and van Houwelingen (1993), both rules, from the same CV run.
##  Part B  optimisation: hierarchical MCP in three scenarios (Base, Strong
##          heterogeneity, p = 500): (B1) five fold assignments, (B2) tolerance
##          1e-4 versus 1e-6, (B3) ratio grid {0.5, 1, 2} versus {1, 2, 4}, and
##          (B4) at the selected (lambda, r) of the reference fit, the objective and
##          support reached from a cold start, from the stratified common-effect
##          fit and from five random perturbations of the reference solution.
##
##  Usage   : NCORES=8 Rscript simulation/07_sensitivity.R [R]      (default R = 20)
##  Output  : results/sensitivity/{A,B}_sXX_rYYY.csv, one file per job (resumable)
##  Next    : 08_sensitivity_tables.R summarises them.
################################################################################
args <- commandArgs(trailingOnly = TRUE)
here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) {   # folder of this script
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))     # works under Rscript
  if (length(f)) dirname(normalizePath(f[1])) else "." })                          # from any directory
source(file.path(here, "00_setup.R"))
Rrep <- if (length(args) >= 1) as.integer(args[1]) else 20L
ncores <- max(1L, as.integer(Sys.getenv("NCORES", "1")))
sdir <- file.path(here, "results", "sensitivity"); dir.create(sdir, showWarnings = FALSE, recursive = TRUE)
SCEN_B <- c("Base", "Strong heterogeneity", "p = 500")

## the data set and folds of replicate r of scenario s, exactly as in 01_run_simulation.R
make_data <- function(s, r) {
  par <- scenario_params(names(SCENARIOS)[s])
  set.seed(20260000 + 1000 * s + r)
  sim <- do.call(sim_ipd, par); d <- sim$data
  X <- std_within(as.matrix(d[, -(1:3)]), d$study)
  list(sim = sim, d = d, X = X, foldid = hmc_folds(d$study, d$status, TUNE$nfolds))
}
cv_fit <- function(dd, pen, ratio = TUNE$ratio, tol = TUNE$tol, foldid = dd$foldid)
  cv.hmcox(dd$X, dd$d$time, dd$d$status, dd$d$study, penalty = pen,
           gamma = if (pen == "SCAD") 3.7 else 3, ratio = ratio, foldid = foldid,
           nlambda = TUNE$nlambda, lambda.min.ratio = TUNE$lambda.min.ratio, tol = tol,
           dfmax = TUNE$dfmax, standardize = FALSE)
## coefficients at grid index ix = c(lambda index, ratio index)
coef_at <- function(cvf, ix) {
  f <- cvf$all.fits[[ix[2]]]
  list(alpha = unname(f$alpha[, ix[1]]), eps = unname(matrix(f$eps[, , ix[1]], nrow = cvf$data$K)))
}
supp <- function(cf) c(paste0("a", which(cf$alpha != 0)),
                       { w <- which(cf$eps != 0, arr.ind = TRUE); if (length(w)) paste0("e", w[, 1], "_", w[, 2]) })
jacc <- function(a, b) if (!length(union(a, b))) 1 else length(intersect(a, b)) / length(union(a, b))
metr <- function(cf, sim) hmc_metrics(cf$alpha, cf$eps, sim)[, c("TPR", "FDR", "MCC", "exact", "MSE_theta", "TPR_dev", "FDR_dev")]

partA <- function(s, r) {
  dd <- make_data(s, r); out <- list()
  for (pen in c("MCP", "SCAD")) {
    cvf <- cv_fit(dd, pen)
    sel_vvh <- hiermetacox:::hmc_rule_index(cvf$cvm.vvh, cvf$cvsd.vvh, cvf$df)
    for (rule in c("min", "1se")) {
      ix_lp <- if (rule == "1se") cvf$sel.1se else cvf$sel.min
      ix_vv <- if (rule == "1se") sel_vvh$se1 else sel_vvh$min
      c_lp <- coef_at(cvf, ix_lp); c_vv <- coef_at(cvf, ix_vv)
      out[[paste(pen, rule)]] <- rbind(
        data.frame(scenario = names(SCENARIOS)[s], rep = r, method = paste0(pen, "-H"), rule = rule,
                   criterion = "lp", metr(c_lp, dd$sim), lambda_index = ix_lp[1], ratio = cvf$ratio.grid[ix_lp[2]],
                   same_as_lp = 1),
        data.frame(scenario = names(SCENARIOS)[s], rep = r, method = paste0(pen, "-H"), rule = rule,
                   criterion = "vvh", metr(c_vv, dd$sim), lambda_index = ix_vv[1], ratio = cvf$ratio.grid[ix_vv[2]],
                   same_as_lp = as.integer(setequal(supp(c_lp), supp(c_vv)))))
    }
  }
  do.call(rbind, out)
}

partB <- function(s, r) {
  dd <- make_data(s, r); scn <- names(SCENARIOS)[s]; rows <- list()
  ref <- cv_fit(dd, "MCP")
  ref_c <- lapply(c(min = "min", `1se` = "1se"), function(rule) coef_at(ref, if (rule == "1se") ref$sel.1se else ref$sel.min))
  add <- function(exper, variant, cvf) for (rule in c("min", "1se")) {
    cf <- coef_at(cvf, if (rule == "1se") cvf$sel.1se else cvf$sel.min)
    rows[[length(rows) + 1]] <<- data.frame(scenario = scn, rep = r, experiment = exper, variant = variant,
      rule = rule, metr(cf, dd$sim), jaccard_ref = jacc(supp(cf), supp(ref_c[[rule]])),
      same_support = as.integer(setequal(supp(cf), supp(ref_c[[rule]]))), objective_diff = NA_real_)
  }
  add("reference", "reference", ref)
  for (sd in 1:4) add("folds", paste0("fold seed ", sd),                              # B1
                      cv_fit(dd, "MCP", foldid = hmc_folds(dd$d$study, dd$d$status, TUNE$nfolds, seed = 1e6 + 100 * r + sd)))
  add("tolerance", "1e-6", cv_fit(dd, "MCP", tol = 1e-6))                            # B2
  add("ratio grid", "{1,2,4}", cv_fit(dd, "MCP", ratio = c(1, 2, 4)))                 # B3
  ## B4: starting values at the selected grid point of the reference fit (minimum rule)
  ix <- ref$sel.min; lam <- ref$lambda.grid[ix[1]]; rat <- ref$ratio.grid[ix[2]]
  D <- ref$data; om <- D$omega
  warm <- coef_at(ref, ix); Qw <- ref$all.fits[[ix[2]]]$objective[ix[1]]
  fit1 <- function(a0, e0) {
    f <- hiermetacox:::hmc_fit_prepared(D, "MCP", 3, rat, lam, TRUE, TUNE$tol, 50, 2000, D$p,
                                        alpha.init = a0, eps.init = sweep(e0, 1, om, "*"))
    list(cf = list(alpha = unname(f$alpha[, 1]), eps = unname(matrix(f$eps[, , 1], nrow = D$K))), Q = f$objective[1])
  }
  stc <- cv.hmcox(dd$X, dd$d$time, dd$d$status, dd$d$study, penalty = "MCP", heterogeneity = FALSE,
                  foldid = dd$foldid, nlambda = TUNE$nlambda, lambda.min.ratio = TUNE$lambda.min.ratio,
                  tol = TUNE$tol, dfmax = TUNE$dfmax, standardize = FALSE)
  starts <- list(cold = list(rep(0, D$p), matrix(0, D$K, D$p)),
                 common = list(unname(coef(stc, rule = "min")$alpha), matrix(0, D$K, D$p)))
  set.seed(30000000 + 1000 * s + r)
  for (q in 1:5) {                       # perturb the active global effects; deviations kept: feasible
    a0 <- warm$alpha; act <- which(a0 != 0); a0[act] <- a0[act] + rnorm(length(act), 0, 0.25)
    starts[[paste0("perturbed ", q)]] <- list(a0, warm$eps)
  }
  for (nm in names(starts)) {
    z <- fit1(starts[[nm]][[1]], starts[[nm]][[2]])
    rows[[length(rows) + 1]] <- data.frame(scenario = scn, rep = r, experiment = "starting value", variant = nm,
      rule = "min", metr(z$cf, dd$sim), jaccard_ref = jacc(supp(z$cf), supp(warm)),
      same_support = as.integer(setequal(supp(z$cf), supp(warm))), objective_diff = z$Q - Qw)
  }
  do.call(rbind, rows)
}

jobs <- rbind(expand.grid(part = "A", s = seq_along(SCENARIOS), r = seq_len(Rrep), stringsAsFactors = FALSE),
              expand.grid(part = "B", s = match(SCEN_B, names(SCENARIOS)), r = seq_len(Rrep), stringsAsFactors = FALSE))
jobs$file <- file.path(sdir, sprintf("%s_s%02d_r%03d.csv", jobs$part, jobs$s, jobs$r))
jobs <- jobs[!file.exists(jobs$file), ]
cat(sprintf("[%s] %d sensitivity job(s) on %d core(s)\n", format(Sys.time(), "%H:%M:%S"), nrow(jobs), ncores))
runner <- function(i) {
  res <- tryCatch(if (jobs$part[i] == "A") partA(jobs$s[i], jobs$r[i]) else partB(jobs$s[i], jobs$r[i]),
                  error = function(e) { message(jobs$file[i], ": ", conditionMessage(e)); NULL })
  if (!is.null(res)) { tmp <- paste0(jobs$file[i], ".tmp"); write.csv(res, tmp, row.names = FALSE); file.rename(tmp, jobs$file[i]) }
  invisible(NULL)
}
if (nrow(jobs)) {
  if (ncores > 1 && .Platform$OS.type == "unix")
    invisible(parallel::mclapply(seq_len(nrow(jobs)), runner, mc.cores = ncores, mc.preschedule = FALSE))
  else for (i in seq_len(nrow(jobs))) runner(i)
}
cat("done\n")
