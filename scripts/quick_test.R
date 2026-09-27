################################################################################
##  scripts/quick_test.R  --  run this first (about one minute).
##  Runs two replicates of the base simulation scenario through exactly the code
##  used for the paper (00_setup.R / fit_method / hmc_metrics) and prints the
##  metrics, then fits and summarises one cross-validated model.
################################################################################
source("simulation/00_setup.R")
set.seed(1)
par <- scenario_params("Base")
out <- list()
for (r in 1:2) {
  sim <- do.call(sim_ipd, par); d <- sim$data
  X <- std_within(as.matrix(d[, -(1:3)]), d$study)
  fid <- hmc_folds(d$study, d$status, TUNE$nfolds)
  for (m in c("MCP-H", "MCP-S", "LASSO-P")) {
    f <- fit_method(m, X, d$time, d$status, d$study, fid)
    for (rule in names(f))
      out[[length(out) + 1]] <- data.frame(rep = r, method = m, rule = rule,
                                           hmc_metrics(f[[rule]]$alpha, f[[rule]]$eps, sim)[, c("TPR", "FDR", "MCC", "MSE_theta", "n_dev")])
  }
}
print(do.call(rbind, out), digits = 3)
fit <- cv.hmcox(X, d$time, d$status, d$study, penalty = "MCP", nfolds = 5, seed = 1, rule = "1se")
print(fit)
print(hmcox_inference(fit)$global[, c("covariate", "HR", "HR_lower", "HR_upper", "p_value")], digits = 3)
cat("quick test finished\n")
