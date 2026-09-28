################################################################################
##  simulation/04_oracle_tables_figures.R
##  Aggregates results/oracle/*.csv (nothing is refitted) and writes the tables and
##  figures of Section 5.5, which verify Theorems 4.1 and 4.3 and Corollary 4.1
##  empirically for the proposed hierarchical MCP and SCAD estimators.
##
##  Output (copy output/ to the Overleaf project):
##    output/tables/tab_oracle_selection.tex   Table 4  (support recovery, equality)
##    output/tables/tab_oracle_inference.tex   Table 5  (bias, SD/SE, coverage, ...)
##    output/figures/oracle_selection.png      Figure 5
##    output/figures/oracle_qq.png             Figure 6
##    output/tables/tab_nearcancel.tex         Table 6  (near-cancellation, Theorem 4.2)
##    simulation/results/oracle_summary.csv    per-parameter summaries
################################################################################
suppressPackageStartupMessages({ library(ggplot2) })
here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
odir <- file.path(here, "results", "oracle")
tab_dir <- file.path(here, "..", "output", "tables"); fig_dir <- file.path(here, "..", "output", "figures")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
x <- do.call(rbind, lapply(list.files(odir, "^n[0-9]{4}_r[0-9]{3}\\.csv$", full.names = TRUE), read.csv))
LAB <- c("MCP-H" = "MCP, hierarchical (proposed)", "SCAD-H" = "SCAD, hierarchical (proposed)",
         "LASSO-H" = "Lasso, hierarchical (control)")
RULES <- c(theory = "Theory rate", "1se" = "One-SE rule", min = "Minimum rule")
x$method <- factor(x$method, levels = names(LAB)); x$rule <- factor(x$rule, levels = names(RULES))
sizes <- sort(unique(x$n))
cat("replicates per size:", tapply(x$rep, x$n, function(r) length(unique(r))), "\n")
f2 <- function(v) formatC(v, format = "f", digits = 2); f3 <- function(v) formatC(v, format = "f", digits = 3)

## ---- per-replicate indicators (one row per method x rule x n x rep) ---------------
R1 <- aggregate(cbind(correct_support, equal_oracle) ~ method + rule + n + rep, x, max)
sel <- aggregate(cbind(correct_support, equal_oracle) ~ method + rule + n, R1, mean)

## ---- Table 4: support recovery and exact equality with the oracle -------------------
rows <- c()
for (rl in names(RULES)) {
  rows <- c(rows, sprintf("\\multicolumn{%d}{@{}l}{\\emph{%s}} \\\\", length(sizes) + 1, RULES[[rl]]))
  for (m in names(LAB)) {
    cells <- sapply(sizes, function(n) { s <- sel[sel$method == m & sel$rule == rl & sel$n == n, ]
      if (!nrow(s)) "--" else sprintf("%s / %s", f2(s$correct_support), f2(s$equal_oracle)) })
    rows <- c(rows, sprintf("%s & %s \\\\", LAB[[m]], paste(cells, collapse = " & ")))
  }
  if (rl != "min") rows <- c(rows, "\\midrule")
}
writeLines(c(sprintf("\\begin{tabular}{@{}l%s@{}}", strrep("c", length(sizes))), "\\toprule",
             sprintf("& \\multicolumn{%d}{c}{Study size $n_k$} \\\\", length(sizes)),
             sprintf("\\cmidrule(l){2-%d}", length(sizes) + 1),
             sprintf("Method & %s \\\\", paste(sizes, collapse = " & ")), "\\midrule", rows,
             "\\bottomrule", "\\end{tabular}"), file.path(tab_dir, "tab_oracle_selection.tex"))

## ---- per-parameter inference summaries -------------------------------------------------
S <- do.call(rbind, lapply(split(x, list(x$method, x$rule, x$n, x$param), drop = TRUE), function(d) {
  z <- (d$pen - d$truth) / d$oracle_se
  data.frame(method = d$method[1], rule = d$rule[1], n = d$n[1], param = d$param[1],
             type = d$type[1], truth = d$truth[1], R = nrow(d),
             bias = mean(d$pen - d$truth), sd_emp = sd(d$pen),
             se_mean = mean(d$refit_se, na.rm = TRUE),
             coverage = mean(d$covered, na.rm = TRUE),
             eff = mean((d$oracle - d$truth)^2) / mean((d$pen - d$truth)^2),
             sw_p = if (nrow(d) >= 8 && sd(z) > 0) shapiro.test(z)$p.value else NA_real_)
}))
write.csv(S, file.path(here, "results", "oracle_summary.csv"), row.names = FALSE)

## ---- Table 5: inference under the theory-rate tuning -------------------------------------
T5 <- S[S$rule == "theory", ]
rows <- c()
for (m in names(LAB)) {
  rows <- c(rows, sprintf("\\multicolumn{7}{@{}l}{\\emph{%s}} \\\\", LAB[[m]]))
  for (n in sizes) {
    d <- T5[T5$method == m & T5$n == n, ]
    if (!nrow(d)) next
    rows <- c(rows, sprintf("$n_k=%d$ & %s & %s & %s & %s & %s & %d/%d \\\\", n,
                            f3(mean(abs(d$bias))), f2(mean(d$sd_emp / d$se_mean)),
                            f3(mean(d$coverage)), f2(mean(d$eff)),
                            f2(stats::median(d$sw_p, na.rm = TRUE)),
                            sum(d$sw_p < 0.05, na.rm = TRUE), sum(!is.na(d$sw_p))))
  }
  if (m != "LASSO-H") rows <- c(rows, "\\midrule")
}
writeLines(c("\\begin{tabular}{@{}lrrrrrr@{}}", "\\toprule",
             "& Mean $|\\mathrm{bias}|$ & SD/SE & Coverage & Efficiency & Median SW $p$ & SW rejections \\\\",
             "\\midrule", rows, "\\bottomrule", "\\end{tabular}"),
           file.path(tab_dir, "tab_oracle_inference.tex"))

## ---- Figure 5: P(correct support) and P(equal to oracle) versus n -------------------------
L <- rbind(transform(sel, what = "Correct support", y = correct_support),
           transform(sel, what = "Equal to the oracle", y = equal_oracle))
L$rule_lab <- factor(RULES[as.character(L$rule)], levels = RULES)
g <- ggplot(L, aes(n, y, colour = method, shape = method)) +
  geom_line(position = position_dodge(width = 0.06)) +
  geom_point(size = 2, position = position_dodge(width = 0.06)) +
  scale_x_log10(breaks = sizes) +
  scale_colour_manual(values = c("MCP-H" = "#B2182B", "SCAD-H" = "#EF8A62", "LASSO-H" = "#2166AC"),
                      labels = LAB, name = NULL) +
  scale_shape_manual(values = c("MCP-H" = 16, "SCAD-H" = 17, "LASSO-H" = 1), labels = LAB, name = NULL) +
  coord_cartesian(ylim = c(0, 1)) + facet_grid(what ~ rule_lab) +
  labs(x = expression("Study size " * n[k] * " (log scale)"), y = "Proportion of replicates") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", legend.box = "vertical",
                                   panel.grid.minor = element_blank(),
                                   strip.background = element_rect(fill = "grey93"))
ggsave(file.path(fig_dir, "oracle_selection.png"), g, width = 8.5, height = 5.6, dpi = 300)

## ---- Figure 6: normal QQ plots of the standardised estimates -------------------------
## Under Theorem 4.3 each standardised estimate (hat - true) / SE_oracle is
## asymptotically N(0, 1); the 14 reduced parameters are pooled within each panel.
Q <- x[x$rule == "theory", ]
Q$z <- (Q$pen - Q$truth) / Q$oracle_se
Q$method_lab <- factor(LAB[as.character(Q$method)], levels = LAB)
Q$n_lab <- factor(paste0("n[k] == ", Q$n), levels = paste0("n[k] == ", sizes))
g <- ggplot(Q, aes(sample = z, colour = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3) +
  stat_qq(size = 0.5, alpha = 0.6) +
  facet_grid(method_lab ~ n_lab, scales = "free_y",
             labeller = labeller(n_lab = label_parsed, method_lab = label_wrap_gen(18))) +
  scale_colour_manual(values = c("MCP-H" = "#B2182B", "SCAD-H" = "#EF8A62", "LASSO-H" = "#2166AC"),
                      guide = "none") +
  labs(x = "Standard normal quantiles", y = expression((hat(vartheta) - vartheta[0]) / SE[oracle])) +
  theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(),
                                  strip.background = element_rect(fill = "grey93"))
ggsave(file.path(fig_dir, "oracle_qq.png"), g, width = 8.5, height = 5.8, dpi = 300)

## ---- Table 6: near-cancellation experiment (05_near_cancellation.R) --------------------
ndir <- file.path(here, "results", "nearcancel")
nf <- list.files(ndir, "^r[0-9]{3}\\.csv$", full.names = TRUE)
if (length(nf)) {
  y <- do.call(rbind, lapply(nf, read.csv))
  a <- aggregate(cbind(theorem_majority, majority, alg_worse) ~ covariate + q + ratio + penalty, y, mean)
  cond_b <- function(q, r) (3 - q) * r^2 > 1          # delta F(l_e) > F(l_a), m = 3
  rows <- c()
  for (cv in c(3, 4)) {
    qv <- a$q[a$covariate == cv][1]
    rows <- c(rows, sprintf("\\multicolumn{8}{@{}l}{\\emph{Covariate %d: $q=%d$, $\\delta=%d$}} \\\\", cv, qv, 3 - qv))
    for (r in sort(unique(a$ratio))) {
      m <- a[a$covariate == cv & a$ratio == r & a$penalty == "MCP", ]
      sc <- a[a$covariate == cv & a$ratio == r & a$penalty == "SCAD", ]
      rows <- c(rows, sprintf("$r=%s$ & %s & %s & %s & %s & %s & %s & %s \\\\", format(r),
        if (cond_b(qv, r)) "yes" else "no",
        f2(m$theorem_majority), f2(m$majority), f2(m$alg_worse),
        f2(sc$theorem_majority), f2(sc$majority), f2(sc$alg_worse)))
    }
    if (cv == 3) rows <- c(rows, "\\midrule")
  }
  writeLines(c("\\begin{tabular}{@{}lccccccc@{}}", "\\toprule",
    "& & \\multicolumn{3}{c}{MCP, hierarchical} & \\multicolumn{3}{c}{SCAD, hierarchical} \\\\",
    "\\cmidrule(lr){3-5}\\cmidrule(l){6-8}",
    "Ratio & (b) holds & Theorem & Algorithm & Inferior & Theorem & Algorithm & Inferior \\\\",
    "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), file.path(tab_dir, "tab_nearcancel.tex"))
  cat("near-cancellation replicates:", length(nf), "\n"); print(a, digits = 2)
}

print(sel, digits = 2)
print(aggregate(cbind(abs(bias), coverage, eff) ~ method + n, T5, mean), digits = 3)
