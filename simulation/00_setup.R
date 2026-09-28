################################################################################
##  simulation/00_setup.R
##  Scenario definitions and competitor wrappers for the simulation study
##  (Section 5 of the manuscript).
##
##  Manuscript : Hierarchical Variable Selection for Nonconvex Penalized Cox
##               Models in Individual-Participant-Data Meta-Analysis
##  Package    : hiermetacox (this repository, folder hiermetacox/)
##  Sourced by : 01_run_simulation.R, 03_coverage_and_timing.R
################################################################################

suppressPackageStartupMessages({
  library(hiermetacox)
  library(survival)
})

## Common tuning settings, stated in Section 3.2 of the manuscript ------------
TUNE <- list(nfolds = 5, ratio = c(0.5, 1, 2), nlambda = 20,
             lambda.min.ratio = 0.05, tol = 1e-4, dfmax = 100)

## Scenarios (Table 1 of the manuscript). One factor changed at a time. -------
BASE <- list(K = 5, n = 200, p = 100,
             alpha = c(1, -0.8, 0.6, -0.5, 0.4, -0.3),
             rho = 0.5, cens = 0.30, het.idx = 1:4, n.dev = 2, sd.eps = 0.4)

SCENARIOS <- list(
  "Base"                    = list(),
  "Independent covariates"  = list(rho = 0),
  "Strong correlation"      = list(rho = 0.8),
  "No heterogeneity"        = list(sd.eps = 0),
  "Strong heterogeneity"    = list(sd.eps = 0.8),
  "Dense heterogeneity"     = list(n.dev = 5),
  "Heavy censoring"         = list(cens = 0.60),
  "Weak signals"            = list(alpha = c(0.5, -0.4, 0.3, -0.25, 0.2, -0.15)),
  "p = 200"                 = list(p = 200),
  "p = 500"                 = list(p = 500),
  "Ten studies"             = list(K = 10, n = 100),
  "Unequal study sizes"     = list(n = c(100, 150, 200, 250, 300)),
  "Small studies, p > n_k"  = list(n = 100, p = 200)
)
scenario_params <- function(name) modifyList(BASE, SCENARIOS[[name]])

## Full 3 x 4 design: model structure x penalty.
##   -H hierarchical (proposed structure), -S stratified common effects,
##   -P naive pooling (unstratified).  Elastic net uses mixing 0.5 throughout.
METHODS <- c("MCP-H", "SCAD-H", "LASSO-H", "ENet-H",
             "MCP-S", "SCAD-S", "LASSO-S", "ENet-S",
             "MCP-P", "SCAD-P", "LASSO-P", "ENet-P")
METHOD_LABELS <- c(
  "MCP-H"   = "MCP, hierarchical (proposed)",
  "SCAD-H"  = "SCAD, hierarchical (proposed)",
  "LASSO-H" = "Lasso, hierarchical",
  "ENet-H"  = "Elastic net, hierarchical",
  "MCP-S"   = "MCP, stratified common effects",
  "SCAD-S"  = "SCAD, stratified common effects",
  "LASSO-S" = "Lasso, stratified common effects",
  "ENet-S"  = "Elastic net, stratified common effects",
  "MCP-P"   = "MCP, naive pooling",
  "SCAD-P"  = "SCAD, naive pooling",
  "LASSO-P" = "Lasso, naive pooling",
  "ENet-P"  = "Elastic net, naive pooling")
ENET_MIX <- 0.5

## within-study standardisation (applied once, used by every method) ----------
std_within <- function(X, study) {
  for (s in unique(study)) {
    r <- study == s
    mu <- colMeans(X[r, , drop = FALSE])
    sdv <- sqrt(colMeans(sweep(X[r, , drop = FALSE], 2, mu)^2)); sdv[sdv < 1e-12] <- Inf
    X[r, ] <- sweep(sweep(X[r, , drop = FALSE], 2, mu), 2, sdv, "/")
  }
  X
}

## one method on one data set; returns a list with one element per tuning rule
## ("min", "1se"), each list(alpha, eps, extra). The two rules are read off the
## same cross-validation run, so they cost nothing extra.
fit_method <- function(method, X, time, status, study, foldid) {
  ## All twelve methods are fitted by the same solver, grid, CV criterion and
  ## folds; they differ only in structure and penalty. Naive pooling is the
  ## common-effect model with a single stratum (one baseline hazard for all).
  pen <- switch(sub("-(H|S|P)$", "", method),
                MCP = "MCP", SCAD = "SCAD", LASSO = "lasso", ENet = "enet")
  structure_ <- sub("^.*-", "", method)
  het <- structure_ == "H"
  strata_ <- if (structure_ == "P") rep(1L, length(time)) else study
  gam <- switch(pen, SCAD = 3.7, enet = ENET_MIX, 3)
  cvf <- cv.hmcox(X, time, status, strata_, penalty = pen, gamma = gam, heterogeneity = het,
                    ratio = TUNE$ratio, foldid = foldid, nlambda = TUNE$nlambda,
                    lambda.min.ratio = TUNE$lambda.min.ratio, tol = TUNE$tol,
                    dfmax = TUNE$dfmax, standardize = FALSE)
    out <- list()
    for (rule in c("min", "1se")) {
      cf <- coef(cvf, rule = rule)
      ix <- if (rule == "1se") cvf$sel.1se else cvf$sel.min
      out[[rule]] <- list(alpha = unname(cf$alpha), eps = if (het) unname(cf$eps) else NULL,
                          cv = cvf,
                          extra = c(ratio = cf$ratio, index = ix[1],
                                    boundary = as.integer(ix[1] == length(cvf$lambda.grid)),
                                    converged = as.integer(all(sapply(cvf$all.fits, function(f) all(f$converged)))),
                                    halvings = sum(sapply(cvf$all.fits, function(f) sum(f$halvings)))))
    }
  out
}
