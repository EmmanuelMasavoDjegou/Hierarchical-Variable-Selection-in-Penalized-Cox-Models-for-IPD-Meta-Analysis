################################################################################
##  simulation/08_sensitivity_tables.R
##  Summarises results/sensitivity/ (07_sensitivity.R); nothing is refitted.
##  Writes:
##    output/tables/tab_cv_criterion.tex   linear-predictor vs Verweij-van Houwelingen
##    output/tables/tab_sensitivity.tex    folds, tolerance, ratio grid, starting values
##    results/sensitivity_cv_by_scenario.csv, results/sensitivity_opt_by_scenario.csv
##  All comparisons are paired (same data and folds); Monte Carlo standard errors are
##  those of the mean paired difference.
################################################################################
here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else "." })
sdir <- file.path(here, "results", "sensitivity")
tab_dir <- file.path(here, "..", "output", "tables"); dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
rd <- function(pat) do.call(rbind, lapply(list.files(sdir, pat, full.names = TRUE), read.csv))
A <- rd("^A_s[0-9]{2}_r[0-9]{3}\\.csv$"); B <- rd("^B_s[0-9]{2}_r[0-9]{3}\\.csv$")
f2 <- function(v) formatC(v, format = "f", digits = 2); f3 <- function(v) formatC(v, format = "f", digits = 3)
se <- function(v) sd(v) / sqrt(length(v))
cat("Part A replicates per scenario:", range(table(A$scenario) / 8), "| Part B:", range(table(B$scenario[B$experiment == "reference"]) / 2), "\n")

## ---- Part A: cross-validation criterion --------------------------------------------
w <- merge(A[A$criterion == "lp", ], A[A$criterion == "vvh", ],
           by = c("scenario", "rep", "method", "rule"), suffixes = c(".lp", ".vvh"))
w$dMCC <- w$MCC.vvh - w$MCC.lp
byscen <- aggregate(cbind(same = same_as_lp.vvh, MCC.lp, MCC.vvh, dMCC) ~ scenario + method + rule, w, mean)
write.csv(byscen, file.path(here, "results", "sensitivity_cv_by_scenario.csv"), row.names = FALSE)
LAB <- c("MCP-H" = "MCP, hierarchical", "SCAD-H" = "SCAD, hierarchical")
RUL <- c("1se" = "one-SE", min = "minimum")
rows <- c()
for (m in names(LAB)) for (rl in names(RUL)) {
  z <- w[w$method == m & w$rule == rl, ]
  sc <- aggregate(dMCC ~ scenario, z, mean)
  rows <- c(rows, sprintf("%s, %s rule & %s & %s & %s & %s (%s) & %s to %s & %s / %s & %s / %s \\\\",
    LAB[[m]], RUL[[rl]], f2(mean(z$same_as_lp.vvh)), f3(mean(z$MCC.lp)), f3(mean(z$MCC.vvh)),
    f3(mean(z$dMCC)), f3(se(z$dMCC)), f3(min(sc$dMCC)), f3(max(sc$dMCC)),
    f2(mean(z$TPR_dev.lp, na.rm = TRUE)), f2(mean(z$TPR_dev.vvh, na.rm = TRUE)),     # undefined without
    f2(mean(z$FDR_dev.lp, na.rm = TRUE)), f2(mean(z$FDR_dev.vvh, na.rm = TRUE))))    # true deviations
}
writeLines(c("\\begin{tabular}{@{}lrrrrrrr@{}}", "\\toprule",
  "& Same & \\multicolumn{2}{c}{MCC} & $\\Delta$MCC & $\\Delta$MCC range & Dev.\\ TPR & Dev.\\ FDR \\\\",
  "\\cmidrule(lr){3-4}",
  "Method and rule & support & LP & VVH & (MC SE) & over scenarios & LP / VVH & LP / VVH \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), file.path(tab_dir, "tab_cv_criterion.tex"))
cat("Part A written:", nrow(w), "paired comparisons\n")

## ---- Part B: optimisation sensitivity -----------------------------------------------
ref <- B[B$experiment == "reference", c("scenario", "rep", "rule", "MCC")]
Bv <- merge(B[B$experiment != "reference", ], ref, by = c("scenario", "rep", "rule"), suffixes = c("", ".ref"))
Bv$dMCC <- Bv$MCC - Bv$MCC.ref
Bv$group <- ifelse(Bv$experiment == "starting value",
                   ifelse(grepl("^perturbed", Bv$variant), "perturbed (5 per replicate)", Bv$variant),
                   ifelse(Bv$experiment == "folds", "four other fold assignments", Bv$variant))
Bv$exper <- c(folds = "Fold assignment", tolerance = "Tolerance", `ratio grid` = "Ratio grid",
              `starting value` = "Starting value")[Bv$experiment]
opt <- aggregate(cbind(same_support, jaccard_ref, dMCC) ~ scenario + exper + group + rule, Bv, mean)
write.csv(opt, file.path(here, "results", "sensitivity_opt_by_scenario.csv"), row.names = FALSE)
SCN <- c("Base", "Strong heterogeneity", "p = 500")
GR <- list(c("Fold assignment", "four other fold assignments", "min"),
           c("Fold assignment", "four other fold assignments", "1se"),
           c("Tolerance", "1e-6", "min"), c("Tolerance", "1e-6", "1se"),
           c("Ratio grid", "{1,2,4}", "min"), c("Ratio grid", "{1,2,4}", "1se"),
           c("Starting value", "cold", "min"), c("Starting value", "common", "min"),
           c("Starting value", "perturbed (5 per replicate)", "min"))
lab_g <- function(g) sprintf("%s: %s, %s rule", g[1],
  c("four other fold assignments" = "4 other fold assignments", "1e-6" = "$\\tau=10^{-6}$",
    "{1,2,4}" = "$r\\in\\{1,2,4\\}$", cold = "cold start", common = "common-effect start",
    "perturbed (5 per replicate)" = "5 perturbed starts")[[g[2]]], RUL[[g[3]]])
rows <- c()
for (g in GR) {
  cells <- sapply(SCN, function(sc) {
    z <- Bv[Bv$scenario == sc & Bv$exper == g[1] & Bv$group == g[2] & Bv$rule == g[3], ]
    sprintf("%s & %s", f2(mean(z$same_support)), f3(mean(z$dMCC)))
  })
  rows <- c(rows, sprintf("%s & %s \\\\", lab_g(g), paste(cells, collapse = " & ")))
}
st <- Bv[Bv$experiment == "starting value", ]
lower <- aggregate(objective_diff ~ scenario + group, st, function(v) mean(v < -1e-8))
mag <- aggregate(objective_diff ~ scenario + group, st, function(v) median(v[v < -1e-8]))
writeLines(c("\\begin{tabular}{@{}lrrrrrr@{}}", "\\toprule",
  sprintf("& %s \\\\", paste(sprintf("\\multicolumn{2}{c}{%s}", SCN), collapse = " & ")),
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(l){6-7}",
  "Variation & Same & $\\Delta$MCC & Same & $\\Delta$MCC & Same & $\\Delta$MCC \\\\",
  "\\midrule", rows, "\\midrule",
  sprintf("Lower objective than warm path, cold start & %s \\\\", paste(sapply(SCN, function(sc)
    sprintf("\\multicolumn{2}{c}{%s}", f2(lower$objective_diff[lower$scenario == sc & lower$group == "cold"]))), collapse = " & ")),
  sprintf("Lower objective than warm path, perturbed starts & %s \\\\", paste(sapply(SCN, function(sc)
    sprintf("\\multicolumn{2}{c}{%s}", f2(lower$objective_diff[lower$scenario == sc & lower$group == "perturbed (5 per replicate)"]))), collapse = " & ")),
  "\\bottomrule", "\\end{tabular}"), file.path(tab_dir, "tab_sensitivity.tex"))
cat("objective decrease from alternative starts (median among decreases):\n"); print(mag, digits = 3)
cat("Part B written\n")
