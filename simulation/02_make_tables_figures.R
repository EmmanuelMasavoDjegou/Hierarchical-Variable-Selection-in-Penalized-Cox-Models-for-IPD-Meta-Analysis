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
##    output/figures/sim_mcc.png           Figure 1
##    output/figures/sim_fdr_tpr.png       Figure 2
##    output/figures/sim_mse_theta.png     Figure 3
##    output/figures/sim_deviations.png    Figure 4
##    simulation/results/sim_summary.csv       long-format summary with MC SEs
################################################################################
suppressPackageStartupMessages({ library(ggplot2) })
here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
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
                          if (rule == "1se") "One-standard-error rule" else "Minimum-deviance rule"),
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
             "\\caption{Monte Carlo means for all scenarios. MSE$(\\theta)$ is multiplied by 100. Bold: best value within scenario and rule. Monte Carlo standard errors are in \\texttt{simulation/results/sim\\_summary.csv}.}\\label{tab:sim_all}\\\\",
             "\\toprule",
             "& \\multicolumn{4}{c}{One-SE rule} & \\multicolumn{4}{c}{Minimum rule} \\\\",
             "\\cmidrule(lr){2-5}\\cmidrule(l){6-9}",
             "Method & MCC & FDR & Exact & MSE$(\\theta)$ & MCC & FDR & Exact & MSE$(\\theta)$ \\\\",
             "\\midrule", "\\endfirsthead",
             "\\toprule Method & MCC & FDR & Exact & MSE$(\\theta)$ & MCC & FDR & Exact & MSE$(\\theta)$ \\\\ \\midrule",
             "\\endhead", "\\bottomrule", "\\endlastfoot", al, "\\end{longtable}"),
           file.path(tab_dir, "tab_sim_all.tex"))

## ---- Figures ----------------------------------------------------------------------
## Colour encodes the penalty and marker shape the model structure, so that the
## 3 x 4 design can be read at a glance.
S$penalty <- factor(sub("-(H|S|P)$", "", S$method), levels = c("MCP", "SCAD", "LASSO", "ENet"),
                    labels = c("MCP", "SCAD", "Lasso", "Elastic net"))
S$structure <- factor(sub("^.*-", "", S$method), levels = c("H", "S", "P"),
                      labels = c("Hierarchical", "Stratified common effects", "Naive pooling"))
pen_cols <- c("MCP" = "#B2182B", "SCAD" = "#EF8A62", "Lasso" = "#2166AC", "Elastic net" = "#7F7F7F")
str_shp <- c("Hierarchical" = 16, "Stratified common effects" = 15, "Naive pooling" = 2)
S$scen <- factor(S$scenario, levels = rev(names(SCENARIOS)))
S$rule_lab <- factor(ifelse(S$rule == "1se", "One-standard-error rule", "Minimum-deviance rule"),
                     levels = c("One-standard-error rule", "Minimum-deviance rule"))
th <- theme_bw(base_size = 10) + theme(legend.position = "bottom", legend.box = "vertical",
                                       legend.margin = margin(0, 0, 0, 0),
                                       panel.grid.minor = element_blank(),
                                       strip.background = element_rect(fill = "grey93"))
dodge <- position_dodge(width = 0.8)
dotplot <- function(dat, facet, xlab, log = FALSE) {
  g <- ggplot(dat, aes(x = y, y = scen, colour = penalty, shape = structure, group = method)) +
    geom_errorbarh(aes(xmin = y - 1.96 * se, xmax = y + 1.96 * se), height = 0,
                   position = dodge, linewidth = 0.3) +
    geom_point(position = dodge, size = 1.4) +
    scale_colour_manual(values = pen_cols, name = "Penalty") +
    scale_shape_manual(values = str_shp, name = "Structure") +
    facet_wrap(facet, scales = "free_x") + labs(x = xlab, y = NULL) + th +
    guides(colour = guide_legend(nrow = 1, order = 1), shape = guide_legend(nrow = 1, order = 2))
  if (log) g <- g + scale_x_log10()
  g
}
dotfig <- function(v, lab, file, h = 9, w = 8.5, log = FALSE) {
  d <- S; d$y <- d[[v]]; d$se <- d[[paste0(v, "_se")]]
  ggsave(file.path(fig_dir, file), dotplot(d, ~rule_lab, lab, log), width = w, height = h, dpi = 300)
}
dotfig("MCC", "Matthews correlation coefficient (global selection)", "sim_mcc.png")
dotfig("MSE_theta", "MSE of study-specific log-hazard ratios (log scale)", "sim_mse_theta.png", log = TRUE)

L <- rbind(transform(S, metric = "False discovery rate", y = FDR, se = FDR_se),
           transform(S, metric = "True positive rate", y = TPR, se = TPR_se))
L <- L[L$rule == "1se", ]
ggsave(file.path(fig_dir, "sim_fdr_tpr.png"),
       dotplot(L, ~metric, "One-standard-error rule"), width = 8.5, height = 9, dpi = 300)

H <- S[S$structure == "Hierarchical" & !S$scenario %in% c("No heterogeneity"), ]
H2 <- rbind(transform(H, metric = "Deviations: true positive rate", y = TPR_dev, se = TPR_dev_se),
            transform(H, metric = "Deviations: false discovery rate", y = FDR_dev, se = FDR_dev_se))
H2$metric <- factor(H2$metric, levels = c("Deviations: true positive rate", "Deviations: false discovery rate"))
g <- ggplot(H2, aes(x = y, y = scen, colour = penalty, shape = rule_lab,
                    group = interaction(method, rule_lab))) +
  geom_errorbarh(aes(xmin = y - 1.96 * se, xmax = y + 1.96 * se), height = 0,
                 position = dodge, linewidth = 0.3) +
  geom_point(position = dodge, size = 1.4) +
  scale_colour_manual(values = pen_cols, name = "Penalty (hierarchical)") +
  scale_shape_manual(values = c(16, 1), name = NULL) +
  facet_wrap(~metric, scales = "free_x") + labs(x = NULL, y = NULL) + th +
  guides(colour = guide_legend(nrow = 1, order = 1), shape = guide_legend(nrow = 1, order = 2))
ggsave(file.path(fig_dir, "sim_deviations.png"), g, width = 8.5, height = 8, dpi = 300)

## key numbers used in the text (printed so they can be checked against main.tex)
key <- S[S$scenario == "Base", c("method", "rule", "TPR", "FDR", "MCC", "exact", "MSE_theta", "TPR_dev", "FDR_dev", "seconds")]
print(key[order(key$rule, key$method), ], digits = 3)
cat("hierarchy violations (must be 0):", sum(raw$hier_violation), "\n")
cat("non-converged hierarchical CV runs:", sum(raw$converged == 0, na.rm = TRUE), "of",
    sum(!is.na(raw$converged)), "| boundary lambda:", mean(raw$boundary, na.rm = TRUE), "\n")
