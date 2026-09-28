################################################################################
##  simulation/01_run_simulation.R
##  Monte Carlo study of Section 5: 13 scenarios x 12 methods x R replicates.
##
##  Usage   : Rscript 01_run_simulation.R [R] [scenario numbers...]
##            e.g.  Rscript 01_run_simulation.R 50           (all scenarios; paper)
##                  Rscript 01_run_simulation.R 50 1 2 3     (subset)
##            NCORES=8 Rscript 01_run_simulation.R 50        (parallel, Unix/macOS)
##  Output  : results/raw/sXX_rYYY.csv        metrics, one file per replicate
##            results/raw/sXX_rYYY_cov.csv    Wald-interval checks (MCP-H, SCAD-H)
##  Resume  : a replicate whose file exists is skipped, so the run can be
##            interrupted and re-issued at any time; files are written
##            atomically (temporary name, then rename).
##  Seeds   : replicate r of scenario s uses set.seed(20260000 + 1000 * s + r),
##            so every replicate can be regenerated in isolation.
##  Next    : 02_make_tables_figures.R builds Tables 1-3, A1 and Figures 1-4.
##  Runtime : about 4 h on one core for R = 50; divide by NCORES.
################################################################################

args <- commandArgs(trailingOnly = TRUE)
here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
source(file.path(here, "00_setup.R"))
Rrep <- if (length(args) >= 1) as.integer(args[1]) else 50L
which_sc <- if (length(args) >= 2) as.integer(args[-1]) else seq_along(SCENARIOS)
ncores <- max(1L, as.integer(Sys.getenv("NCORES", "1")))
rawdir <- file.path(here, "results", "raw"); dir.create(rawdir, showWarnings = FALSE, recursive = TRUE)

run_one <- function(s, r) {
  scn <- names(SCENARIOS)[s]; par <- scenario_params(scn)
  f_raw <- file.path(rawdir, sprintf("s%02d_r%03d.csv", s, r))
  if (file.exists(f_raw)) return(invisible(NULL))
  set.seed(20260000 + 1000 * s + r)
  sim <- do.call(sim_ipd, par)
  d <- sim$data
  X <- std_within(as.matrix(d[, -(1:3)]), d$study)
  foldid <- hmc_folds(d$study, d$status, TUNE$nfolds)
  rows <- list(); covs <- list()
  for (m in METHODS) {
    t0 <- proc.time()[3]
    fits <- tryCatch(fit_method(m, X, d$time, d$status, d$study, foldid),
                     error = function(e) { message(scn, " r", r, " ", m, ": ", conditionMessage(e)); NULL })
    tm <- proc.time()[3] - t0
    if (is.null(fits)) next
    for (rule in names(fits)) {
      fit <- fits[[rule]]
      rows[[paste(m, rule)]] <- data.frame(scenario = scn, rep = r, method = m, rule = rule,
                                           hmc_metrics(fit$alpha, fit$eps, sim),
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
  if (length(covs)) {
    tmp <- paste0(sub("\\.csv$", "_cov.csv", f_raw), ".tmp")
    write.csv(do.call(rbind, covs), tmp, row.names = FALSE)
    file.rename(tmp, sub("\\.csv$", "_cov.csv", f_raw))
  }
  tmp <- paste0(f_raw, ".tmp")                       # metrics file last: marks completion
  write.csv(do.call(rbind, rows), tmp, row.names = FALSE)
  file.rename(tmp, f_raw)
  invisible(NULL)
}

jobs <- expand.grid(r = seq_len(Rrep), s = which_sc)
jobs <- jobs[!file.exists(file.path(rawdir, sprintf("s%02d_r%03d.csv", jobs$s, jobs$r))), ]
cat(sprintf("[%s] %d replicate(s) to run on %d core(s)\n", format(Sys.time(), "%H:%M:%S"),
            nrow(jobs), ncores))
if (nrow(jobs)) {
  runner <- function(i) { run_one(jobs$s[i], jobs$r[i])
    if (i %% 10 == 0) cat(sprintf("[%s] %d/%d done\n", format(Sys.time(), "%H:%M:%S"), i, nrow(jobs))) }
  if (ncores > 1 && .Platform$OS.type == "unix") {
    invisible(parallel::mclapply(seq_len(nrow(jobs)), runner, mc.cores = ncores,
                                 mc.preschedule = FALSE))
  } else for (i in seq_len(nrow(jobs))) runner(i)
}
cat("done\n")
