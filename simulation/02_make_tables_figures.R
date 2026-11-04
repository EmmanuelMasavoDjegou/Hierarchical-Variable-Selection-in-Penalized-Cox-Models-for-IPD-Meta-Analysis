################################################################################
##  simulation/02_make_tables_figures.R
##  Aggregates results/raw/sXX_rYYY*.csv (no model is
##  refitted) and writes every simulation table and figure of the manuscript.
##
##  Output (copy output/ to the Overleaf project; figures/ and tables/ keep these names):
##    output/tables/tab_scenarios.tex      Table 1
##    output/tables/tab_sim_base.tex       Table 2
##    output/tables/tab_sim_coverage.tex   Table 3
##    output/tables/tab_sim_all.tex        Table A1 (appendix)
##    output/figures/sim_*_{1se,min}.png   Figures 1-10 (drawn by 09_simulation_figures.R)
##    simulation/results/sim_summary.csv       long-format summary with MC SEs
################################################################################
suppressPackageStartupMessages({ library(ggplot2) })
## ggplot2 >= 4 needs S7's `@` method, which S7 registers only if it was installed under
## R >= 4.3.0; a copy built under an older R fails later with "no applicable method for `@`".
if (packageVersion("ggplot2") >= "4.0.0" && getRversion() >= "4.3.0" &&
    !isTRUE(tryCatch(is.logical(theme_bw()@complete), error = function(e) FALSE)))
  stop("package S7 was installed under an R version older than 4.3.0 and must be reinstalled ",
       "under this R (", R.version.string, "). In a fresh session run:\n",
       "  install.packages(c(\"S7\", \"ggplot2\"))\nsee the README, section 'Installation'.", call. = FALSE)
here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) {   # folder of this script
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))     # works under Rscript
  if (length(f)) dirname(normalizePath(f[1])) else "." })                          # from any directory
source(file.path(here, "00_setup.R"))
res_dir <- file.path(here, "results")
tab_dir <- file.path(here, "..", "output", "tables")
fig_dir <- file.path(here, "..", "output", "figures")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

rawdir <- file.path(res_dir, "raw")
raw <- do.call(rbind, lapply(list.files(rawdir, "^s[0-9]{2}_r[0-9]{3}\\.csv$", full.names = TRUE), read.csv))
cov <- do.call(rbind, lapply(list.files(rawdir, "^s[0-9]{2}_r[0-9]{3}_cov\\.csv$", full.names = TRUE), read.csv))
raw$scenario <- factor(raw$scenario, levels = names(SCENARIOS))
raw$method <- factor(raw$method, levels = METHODS)
R_used <- tapply(raw$rep, raw$scenario, function(x) length(unique(x)))
cat("replicates per scenario:\n"); print(R_used)

metrics <- c("TPR", "FDR", "F1", "MCC", "exact", "size", "alpha_bias", "fp_bias",
             "MSE_alpha", "MSE_eps", "MSE_theta", "TPR_dev", "FDR_dev", "n_dev",
             "n_near_cancel", "err_theta_nc", "seconds")
S <- do.call(rbind, lapply(split(raw, list(raw$scenario, raw$method, raw$rule), drop = TRUE), function(d) {
  data.frame(scenario = d$scenario[1], method = d$method[1], rule = d$rule[1], R = nrow(d),
             t(sapply(metrics, function(m) mean(d[[m]], na.rm = TRUE))),
             t(setNames(sapply(metrics, function(m) sd(d[[m]], na.rm = TRUE) / sqrt(sum(!is.na(d[[m]])))),
                        paste0(metrics, "_se"))))
}))
write.csv(S, file.path(res_dir, "sim_summary.csv"), row.names = FALSE)

f3 <- function(x) ifelse(is.na(x), "--", formatC(x, format = "f", digits = 3))
f2 <- function(x) ifelse(is.na(x), "--", formatC(x, format = "f", digits = 2))
bold_best <- function(v, fmt, best = c("max", "min")) {
  s <- fmt(v); ok <- !is.na(v)
  b <- if (match.arg(best) == "max") max(v[ok]) else min(v[ok])
  s[ok & abs(v - b) < 1e-12] <- paste0("\\textbf{", s[ok & abs(v - b) < 1e-12], "}")
  s
}

## ---- Table 1: scenarios -----------------------------------------------------
desc <- c("Reference configuration",
          "$\\rho=0$", "$\\rho=0.8$", "$\\sigma_\\varepsilon=0$ (no deviations)",
          "$\\sigma_\\varepsilon=0.8$", "all $K$ studies deviate",
          "censoring $60\\%$", "$\\boldsymbol\\alpha_{\\mathcal A_0}$ halved",
          "$p=200$", "$p=500$", "$K=10$, $n_k=100$",
          "$n_k=100,150,200,250,300$", "$n_k=100$, $p=200$ ($p>n_k$)")
lines <- c("\\begin{tabular}{@{}rll@{}}", "\\toprule",
           "\\# & Scenario & Change from the base configuration \\\\", "\\midrule",
           sprintf("%d & %s & %s \\\\", seq_along(SCENARIOS),
                   ifelse(names(SCENARIOS) == "Dense heterogeneity",
                          "Dense heterogeneity$^{\\dagger}$", names(SCENARIOS)), desc),
           "\\bottomrule", "\\end{tabular}")
writeLines(gsub("p > n_k", "$p>n_k$", lines), file.path(tab_dir, "tab_scenarios.tex"))

## ---- Table 2: base scenario, both rules ----------------------------------------
mk_rows <- function(sc, rule) {
  d <- S[S$scenario == sc & S$rule == rule, ]; d <- d[match(METHODS, d$method), ]
  cbind(METHOD_LABELS[as.character(d$method)],
        bold_best(d$TPR, f3), bold_best(d$FDR, f3, "min"), bold_best(d$MCC, f3),
        bold_best(d$exact, f2), f2(d$size), bold_best(d$MSE_alpha * 100, f3, "min"),
        bold_best(d$MSE_eps * 100, f3, "min"), bold_best(d$MSE_theta * 100, f3, "min"),
        f3(ifelse(grepl("-H$", d$method), d$TPR_dev, NA)),
        f3(ifelse(grepl("-H$", d$method), d$FDR_dev, NA)), f2(d$seconds))
}
hdr <- c("& \\multicolumn{5}{c}{Global support} & \\multicolumn{3}{c}{MSE $\\times100$} & \\multicolumn{2}{c}{Deviation support} & \\\\",
         "\\cmidrule(lr){2-6}\\cmidrule(lr){7-9}\\cmidrule(lr){10-11}",
         "Method & TPR & FDR & MCC & Exact & $|\\hat{\\mathcal A}|$ & $\\alpha$ & $\\varepsilon$ & $\\theta$ & TPR & FDR & Time (s) \\\\")
body <- c()
for (rule in c("1se", "min")) {
  body <- c(body, sprintf("\\multicolumn{12}{@{}l}{\\emph{%s}} \\\\",
                          if (rule == "1se") "One-standard-error rule" else "Minimum rule"),
            apply(mk_rows("Base", rule), 1, paste, collapse = " & "))
  body <- paste0(body, ifelse(grepl("\\\\\\\\$", body) | grepl("^\\\\midrule$", body), "", " \\\\"))
  if (rule == "1se") body <- c(body, "\\midrule")
}
writeLines(c("\\begin{tabular}{@{}lrrrrrrrrrrr@{}}", "\\toprule", hdr, "\\midrule", body,
             "\\bottomrule", "\\end{tabular}"), file.path(tab_dir, "tab_sim_base.tex"))

## ---- Table 3: coverage of post-selection Wald intervals ------------------------
cov$scenario <- factor(cov$scenario, levels = names(SCENARIOS))
ca <- cov[cov$truth != 0, ]
cc <- aggregate(cbind(covered, n = 1) ~ scenario + method + rule, ca, sum)
cc$coverage <- cc$covered / cc$n
cstr <- aggregate(covered ~ scenario + method + rule, ca[abs(ca$truth) >= 0.5, ], mean)
names(cstr)[4] <- "cov_strong"
cc <- merge(cc, cstr, all.x = TRUE)
sel_rate <- aggregate(TPR ~ scenario + method + rule, raw[raw$method %in% c("MCP-H", "SCAD-H"), ], mean)
cc <- merge(cc, sel_rate)
cc <- cc[order(cc$scenario, cc$method, cc$rule), ]
write.csv(cc, file.path(res_dir, "sim_coverage_summary.csv"), row.names = FALSE)
show <- c("Base", "No heterogeneity", "Strong heterogeneity", "Heavy censoring",
          "Weak signals", "p = 500", "Small studies, p > n_k")
cm <- cc[cc$scenario %in% show & cc$method == "MCP-H", ]
tl <- c()
for (sc in show) {
  a <- cm[cm$scenario == sc & cm$rule == "1se", ]; b <- cm[cm$scenario == sc & cm$rule == "min", ]
  if (!nrow(a) || !nrow(b)) next
  tl <- c(tl, sprintf("%s & %s & %s & %s & %s & %s & %s \\\\", gsub("p > n_k", "$p>n_k$", sc),
                      f3(a$TPR), f3(a$coverage), f3(a$cov_strong), f3(b$TPR), f3(b$coverage), f3(b$cov_strong)))
}
writeLines(c("\\begin{tabular}{@{}lrrrrrr@{}}", "\\toprule",
             "& \\multicolumn{3}{c}{One-SE rule} & \\multicolumn{3}{c}{Minimum rule} \\\\",
             "\\cmidrule(lr){2-4}\\cmidrule(l){5-7}",
             "Scenario & TPR & Cov. (all) & Cov. ($|\\alpha_j|\\ge 0.5$) & TPR & Cov. (all) & Cov. ($|\\alpha_j|\\ge 0.5$) \\\\",
             "\\midrule", tl, "\\bottomrule", "\\end{tabular}"),
           file.path(tab_dir, "tab_sim_coverage.tex"))

## ---- Table A1: all scenarios --------------------------------------------------
al <- c()
for (sc in names(SCENARIOS)) {
  a <- S[S$scenario == sc & S$rule == "1se", ]; a <- a[match(METHODS, a$method), ]
  b <- S[S$scenario == sc & S$rule == "min", ]; b <- b[match(METHODS, b$method), ]
  if (!nrow(a) || anyNA(a$method)) next
  rows <- cbind(METHOD_LABELS[METHODS], bold_best(a$MCC, f3), bold_best(a$FDR, f3, "min"), bold_best(a$exact, f2),
                bold_best(a$MSE_theta * 100, f3, "min"),
                bold_best(b$MCC, f3), bold_best(b$FDR, f3, "min"), bold_best(b$exact, f2),
                bold_best(b$MSE_theta * 100, f3, "min"))
  al <- c(al, sprintf("\\multicolumn{9}{@{}l}{\\textbf{%s} ($R=%d$)} \\\\",
                      gsub("p > n_k", "$p>n_k$", sc), R_used[[sc]]),
          paste(apply(rows, 1, paste, collapse = " & "), "\\\\"), "\\midrule")
}
al <- al[-length(al)]
writeLines(c("\\begin{longtable}{@{}lrrrrrrrr@{}}",
             "\\caption{Monte Carlo means for all scenarios. MSE$(\\theta)$ is multiplied by 100. Bold: best value within scenario and rule. Monte Carlo standard errors are written to \\texttt{simulation/results/sim\\_summary.csv} by \\texttt{simulation/02\\_make\\_tables\\_figures.R}.}\\label{tab:sim_all}\\\\",
             "\\toprule",
             "& \\multicolumn{4}{c}{One-SE rule} & \\multicolumn{4}{c}{Minimum rule} \\\\",
             "\\cmidrule(lr){2-5}\\cmidrule(l){6-9}",
             "Method & MCC & FDR & Exact & MSE$(\\theta)$ & MCC & FDR & Exact & MSE$(\\theta)$ \\\\",
             "\\midrule", "\\endfirsthead",
             "\\toprule Method & MCC & FDR & Exact & MSE$(\\theta)$ & MCC & FDR & Exact & MSE$(\\theta)$ \\\\ \\midrule",
             "\\endhead", "\\bottomrule", "\\endlastfoot", al, "\\end{longtable}"),
           file.path(tab_dir, "tab_sim_all.tex"))

## ---- Figures ----------------------------------------------------------------------
## Drawn by 09_simulation_figures.R from results/sim_summary.csv (written above).
source(file.path(here, "09_simulation_figures.R"), local = new.env())

## key numbers used in the text (printed so they can be checked against main.tex)
key <- S[S$scenario == "Base", c("method", "rule", "TPR", "FDR", "MCC", "exact", "MSE_theta", "TPR_dev", "FDR_dev", "seconds")]
print(key[order(key$rule, key$method), ], digits = 3)
cat("hierarchy violations (must be 0):", sum(raw$hier_violation), "\n")

## Termination of the solver (Section 3.4). Each CV run consists of (V + 1) x |ratios|
## path fits of nlambda points each (18 x 20 = 360 path points); every path point ends
## with the tolerance met, with no accepted halved step (stalled) or at the iteration
## limit. converged, stalled and iter_limit are properties of the CV run and are the
## same for both rules, so each fit is counted once (rows of the minimum rule);
## boundary (selected lambda = smallest grid value) depends on the rule.
n_pts <- (TUNE$nfolds + 1) * length(TUNE$ratio) * TUNE$nlambda
fits <- raw[raw$rule == "min", ]
conv <- do.call(rbind, lapply(split(fits, list(fits$scenario, fits$method), drop = TRUE), function(z)
  data.frame(scenario = z$scenario[1], method = z$method[1], fits = nrow(z),
             fits_with_iter_limit = sum(z$iter_limit > 0),
             mean_iter_limit_points = mean(z$iter_limit), mean_stalled_points = mean(z$stalled),
             path_points_per_fit = n_pts,
             boundary_min = mean(raw$boundary[raw$rule == "min" & raw$scenario == z$scenario[1] & raw$method == z$method[1]]),
             boundary_1se = mean(raw$boundary[raw$rule == "1se" & raw$scenario == z$scenario[1] & raw$method == z$method[1]]))))
conv <- conv[order(conv$scenario, conv$method), ]
write.csv(conv, file.path(res_dir, "sim_convergence_summary.csv"), row.names = FALSE)
cat(sprintf("fits (CV runs) with >= 1 path point at the iteration limit: %d of %d (%.1f%%); path points at the limit: %d of %d (%.3f%%); stalled path points: %d\n",
            sum(fits$iter_limit > 0), nrow(fits), 100 * mean(fits$iter_limit > 0),
            sum(fits$iter_limit), nrow(fits) * n_pts, 100 * sum(fits$iter_limit) / (nrow(fits) * n_pts),
            sum(fits$stalled)))
cat("  by structure (fits with >= 1 iteration-limit point):\n")
print(tapply(fits$iter_limit > 0, sub("^.*-", "", fits$method), function(v) sprintf("%d of %d", sum(v), length(v))))
cat(sprintf("selected lambda at the smallest grid value: %.1f%% of fits (minimum rule), %.1f%% (one-SE rule)\n",
            100 * mean(raw$boundary[raw$rule == "min"]), 100 * mean(raw$boundary[raw$rule == "1se"])))
cat("details by scenario and method: simulation/results/sim_convergence_summary.csv\n")
