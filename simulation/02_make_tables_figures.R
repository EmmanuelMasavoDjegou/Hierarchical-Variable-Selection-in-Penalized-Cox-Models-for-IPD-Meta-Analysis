################################################################################
##  simulation/02_make_tables_figures.R
##  Aggregates results/sim_raw_*.csv and results/sim_coverage_*.csv (no model is
##  refitted) and writes every simulation table and figure of the manuscript.
##
##  Output (paths fixed by \input / \includegraphics in manuscript/main.tex):
##    manuscript/tables/tab_scenarios.tex      Table 1
##    manuscript/tables/tab_sim_base.tex       Table 2
##    manuscript/tables/tab_sim_coverage.tex   Table 3
##    manuscript/tables/tab_sim_all.tex        Table A1 (appendix)
##    manuscript/figures/sim_mcc.pdf           Figure 1
##    manuscript/figures/sim_fdr_tpr.pdf       Figure 2
##    manuscript/figures/sim_mse_theta.pdf     Figure 3
##    manuscript/figures/sim_deviations.pdf    Figure 4
##    simulation/results/sim_summary.csv       long-format summary with MC SEs
################################################################################
suppressPackageStartupMessages({ library(ggplot2) })
here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
source(file.path(here, "00_setup.R"))
res_dir <- file.path(here, "results")
tab_dir <- file.path(here, "..", "manuscript", "tables")
fig_dir <- file.path(here, "..", "manuscript", "figures")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

raw <- do.call(rbind, lapply(list.files(res_dir, "^sim_raw_.*csv$", full.names = TRUE), read.csv))
cov <- do.call(rbind, lapply(list.files(res_dir, "^sim_coverage_.*csv$", full.names = TRUE), read.csv))
raw$scenario <- factor(raw$scenario, levels = names(SCENARIOS))
raw$method <- factor(raw$method, levels = METHODS)
R_used <- tapply(raw$rep, raw$scenario, function(x) length(unique(x)))
cat("replicates per scenario:\n"); print(R_used)

metrics <- c("TPR", "FDR", "F1", "MCC", "exact", "size", "alpha_bias", "fp_bias",
             "MSE_theta", "TPR_dev", "FDR_dev", "n_dev", "seconds")
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
           sprintf("%d & %s & %s \\\\", seq_along(SCENARIOS), names(SCENARIOS), desc),
           "\\bottomrule", "\\end{tabular}")
writeLines(gsub("p > n_k", "$p>n_k$", lines), file.path(tab_dir, "tab_scenarios.tex"))

## ---- Table 2: base scenario, both rules ----------------------------------------
mk_rows <- function(sc, rule) {
  d <- S[S$scenario == sc & S$rule == rule, ]; d <- d[match(METHODS, d$method), ]
  cbind(METHOD_LABELS[as.character(d$method)],
        bold_best(d$TPR, f3), bold_best(d$FDR, f3, "min"), bold_best(d$MCC, f3),
        bold_best(d$exact, f2), f2(d$size), bold_best(d$MSE_theta * 100, f3, "min"),
        f3(ifelse(grepl("-H$", d$method), d$TPR_dev, NA)),
        f3(ifelse(grepl("-H$", d$method), d$FDR_dev, NA)), f2(d$seconds))
}
hdr <- "Method & TPR & FDR & MCC & Exact & $|\\hat{\\mathcal A}|$ & MSE$(\\theta)\\times 100$ & TPR$_{\\rm dev}$ & FDR$_{\\rm dev}$ & Time (s) \\\\"
body <- c()
for (rule in c("1se", "min")) {
  body <- c(body, sprintf("\\multicolumn{10}{@{}l}{\\emph{%s}} \\\\",
                          if (rule == "1se") "One-standard-error rule" else "Minimum-deviance rule"),
            apply(mk_rows("Base", rule), 1, paste, collapse = " & "))
  body <- paste0(body, ifelse(grepl("\\\\\\\\$", body), "", " \\\\"))
  if (rule == "1se") body <- c(body, "\\midrule")
}
writeLines(c("\\begin{tabular}{@{}lrrrrrrrrr@{}}", "\\toprule", hdr, "\\midrule", body,
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
  rows <- cbind(METHODS, bold_best(a$MCC, f3), bold_best(a$FDR, f3, "min"), bold_best(a$exact, f2),
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
cols <- c("MCP-H" = "#B2182B", "SCAD-H" = "#EF8A62", "LASSO-H" = "#FDB863",
          "MCP-S" = "#5E3C99", "MCP-P" = "#2166AC", "SCAD-P" = "#67A9CF",
          "LASSO-P" = "#878787", "ENet-P" = "#1A1A1A")
shp <- c("MCP-H" = 16, "SCAD-H" = 16, "LASSO-H" = 16, "MCP-S" = 15,
         "MCP-P" = 17, "SCAD-P" = 17, "LASSO-P" = 17, "ENet-P" = 17)
S$scen <- factor(S$scenario, levels = rev(names(SCENARIOS)))
S$rule_lab <- factor(ifelse(S$rule == "1se", "One-standard-error rule", "Minimum-deviance rule"),
                     levels = c("One-standard-error rule", "Minimum-deviance rule"))
th <- theme_bw(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
                                       strip.background = element_rect(fill = "grey93"))
dodge <- position_dodge(width = 0.7)
dotfig <- function(v, lab, file, h = 7.2, w = 8, log = FALSE, sub = S) {
  sub$y <- sub[[v]]; sub$se <- sub[[paste0(v, "_se")]]
  g <- ggplot(sub, aes(x = y, y = scen, colour = method, shape = method)) +
    geom_errorbarh(aes(xmin = y - 1.96 * se, xmax = y + 1.96 * se), height = 0, position = dodge,
                   linewidth = 0.4) +
    geom_point(position = dodge, size = 1.7) +
    scale_colour_manual(values = cols, labels = METHOD_LABELS, name = NULL) +
    scale_shape_manual(values = shp, labels = METHOD_LABELS, name = NULL) +
    facet_wrap(~rule_lab) + labs(x = lab, y = NULL) + th +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
  if (log) g <- g + scale_x_log10()
  ggsave(file.path(fig_dir, file), g, width = w, height = h)
}
dotfig("MCC", "Matthews correlation coefficient (global selection)", "sim_mcc.pdf")
dotfig("MSE_theta", "MSE of study-specific log-hazard ratios (log scale)", "sim_mse_theta.pdf", log = TRUE)

L <- rbind(transform(S, metric = "False discovery rate", y = FDR, se = FDR_se),
           transform(S, metric = "True positive rate", y = TPR, se = TPR_se))
L <- L[L$rule == "1se", ]
g <- ggplot(L, aes(x = y, y = scen, colour = method, shape = method)) +
  geom_errorbarh(aes(xmin = y - 1.96 * se, xmax = y + 1.96 * se), height = 0, position = dodge, linewidth = 0.4) +
  geom_point(position = dodge, size = 1.7) +
  scale_colour_manual(values = cols, labels = METHOD_LABELS, name = NULL) +
  scale_shape_manual(values = shp, labels = METHOD_LABELS, name = NULL) +
  facet_wrap(~metric, scales = "free_x") + labs(x = "One-standard-error rule", y = NULL) + th +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
ggsave(file.path(fig_dir, "sim_fdr_tpr.pdf"), g, width = 8, height = 7.2)

H <- S[S$method %in% c("MCP-H", "SCAD-H", "LASSO-H") & !S$scenario %in% c("No heterogeneity"), ]
H2 <- rbind(transform(H, metric = "Deviations: true positive rate", y = TPR_dev, se = TPR_dev_se),
            transform(H, metric = "Deviations: false discovery rate", y = FDR_dev, se = FDR_dev_se))
H2$metric <- factor(H2$metric, levels = c("Deviations: true positive rate", "Deviations: false discovery rate"))
g <- ggplot(H2, aes(x = y, y = scen, colour = method, shape = rule_lab)) +
  geom_errorbarh(aes(xmin = y - 1.96 * se, xmax = y + 1.96 * se), height = 0,
                 position = position_dodge(width = 0.7), linewidth = 0.4) +
  geom_point(position = position_dodge(width = 0.7), size = 1.7) +
  scale_colour_manual(values = cols[c("MCP-H", "SCAD-H", "LASSO-H")], labels = METHOD_LABELS, name = NULL) +
  scale_shape_manual(values = c(16, 1), name = NULL) +
  facet_wrap(~metric, scales = "free_x") + labs(x = NULL, y = NULL) + th +
  theme(legend.box = "vertical") +
  guides(colour = guide_legend(nrow = 1, order = 1), shape = guide_legend(nrow = 1, order = 2))
ggsave(file.path(fig_dir, "sim_deviations.pdf"), g, width = 8, height = 6.8)

## key numbers used in the text (printed so they can be checked against main.tex)
key <- S[S$scenario == "Base", c("method", "rule", "TPR", "FDR", "MCC", "exact", "MSE_theta", "TPR_dev", "FDR_dev", "seconds")]
print(key[order(key$rule, key$method), ], digits = 3)
cat("hierarchy violations (must be 0):", sum(raw$hier_violation), "\n")
cat("non-converged hierarchical CV runs:", sum(raw$converged == 0, na.rm = TRUE), "of",
    sum(!is.na(raw$converged)), "| boundary lambda:", mean(raw$boundary, na.rm = TRUE), "\n")
