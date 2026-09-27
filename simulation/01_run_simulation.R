################################################################################
##  simulation/01_run_simulation.R
##  Monte Carlo study of Section 5 (13 scenarios x 8 methods x R replicates).
##
##  Usage   : Rscript 01_run_simulation.R [R] [scenario numbers...]
##            e.g.  Rscript 01_run_simulation.R 25           (all scenarios; paper)
##                  Rscript 01_run_simulation.R 25 1 2 3     (subset)
##  Output  : results/sim_raw_<scenario>.csv      per-replicate metrics
##            results/sim_coverage_<scenario>.csv per-replicate Wald CI checks
##  Next    : 02_make_tables_figures.R builds Tables 1-3, A1 and Figures 1-4.
##  Seeds   : replicate r of scenario s uses set.seed(20260000 + 1000 * s + r).
##  Runtime : about 1.3 h on one core for R = 25.
################################################################################

args <- commandArgs(trailingOnly = TRUE)
here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
source(file.path(here, "00_setup.R"))
Rrep <- if (length(args) >= 1) as.integer(args[1]) else 100L
which_sc <- if (length(args) >= 2) as.integer(args[-1]) else seq_along(SCENARIOS)
outdir <- file.path(here, "results"); dir.create(outdir, showWarnings = FALSE)
slug <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

for (s in which_sc) {
  scn <- names(SCENARIOS)[s]; par <- scenario_params(scn)
  f_raw <- file.path(outdir, sprintf("sim_raw_%02d_%s.csv", s, slug(scn)))
  f_cov <- file.path(outdir, sprintf("sim_coverage_%02d_%s.csv", s, slug(scn)))
  done <- if (file.exists(f_raw)) unique(read.csv(f_raw)$rep) else integer(0)
  cat(sprintf("[%s] scenario %d: %s (%d done)\n", format(Sys.time(), "%H:%M:%S"), s, scn, length(done)))
  for (r in setdiff(seq_len(Rrep), done)) {
    set.seed(20260000 + 1000 * s + r)
    sim <- do.call(sim_ipd, par)
    d <- sim$data
    X <- std_within(as.matrix(d[, -(1:3)]), d$study)
    foldid <- hmc_folds(d$study, d$status, TUNE$nfolds)
    rows <- list(); covs <- list()
    for (m in METHODS) {
      t0 <- proc.time()[3]
      fits <- tryCatch(fit_method(m, X, d$time, d$status, d$study, foldid),
                       error = function(e) { message(m, ": ", conditionMessage(e)); NULL })
      tm <- proc.time()[3] - t0
      if (is.null(fits)) next
      for (rule in names(fits)) {
        fit <- fits[[rule]]
        met <- hmc_metrics(fit$alpha, fit$eps, sim)
        rows[[paste(m, rule)]] <- data.frame(scenario = scn, rep = r, method = m, rule = rule, met,
                                             t(fit$extra), seconds = tm,
                                             cens_obs = 1 - mean(d$status))
        if (m %in% c("MCP-H", "SCAD-H") && any(fit$alpha != 0)) {
          inf <- tryCatch(suppressWarnings(hmcox_inference(fit$cv, rule = rule)$global),
                          error = function(e) NULL)
          if (!is.null(inf)) {
            j <- as.integer(sub("X", "", inf$covariate))
            lo <- log(inf$HR_lower); hi <- log(inf$HR_upper)
            covs[[paste(m, rule)]] <- data.frame(scenario = scn, rep = r, method = m, rule = rule,
                                                 j = j, truth = sim$alpha[j], est = inf$logHR,
                                                 se = inf$SE,
                                                 covered = as.integer(lo <= sim$alpha[j] & sim$alpha[j] <= hi))
          }
        }
      }
    }
    write.table(do.call(rbind, rows), f_raw, sep = ",", row.names = FALSE,
                col.names = !file.exists(f_raw), append = file.exists(f_raw))
    if (length(covs))
      write.table(do.call(rbind, covs), f_cov, sep = ",", row.names = FALSE,
                  col.names = !file.exists(f_cov), append = file.exists(f_cov))
    if (r %% 10 == 0) cat(sprintf("   rep %d  [%s]\n", r, format(Sys.time(), "%H:%M:%S")))
  }
}
cat("done\n")
