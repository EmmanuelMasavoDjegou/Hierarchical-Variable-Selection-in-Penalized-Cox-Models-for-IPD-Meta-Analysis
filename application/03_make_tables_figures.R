################################################################################
##  application/03_make_tables_figures.R
##  Builds every table and figure of Section 6 from the stored outputs of
##  01_prepare_data.R and 02_fit_and_infer.R (nothing is refitted).
##
##  Output (paths fixed by \input / \includegraphics in manuscript/main.tex):
##    manuscript/tables/tab_cohorts.tex       Table 4
##    manuscript/tables/tab_app_methods.tex   Table 5
##    manuscript/tables/tab_app_genes.tex     Table 6
##    manuscript/tables/tab_app_frailty.tex   Table 7
##    manuscript/tables/tab_app_loso.tex      Table 8
##    manuscript/figures/app_km.pdf           Figure 5
##    manuscript/figures/app_cv_surface.pdf   Figure 6
##    manuscript/figures/app_forest.pdf       Figure 7
##    manuscript/figures/app_stability.pdf    Figure 8
################################################################################
suppressPackageStartupMessages({ library(ggplot2); library(survival) })
here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
source(file.path(here, "..", "simulation", "00_setup.R"))
res <- file.path(here, "results")
tab_dir <- file.path(here, "..", "manuscript", "tables")
fig_dir <- file.path(here, "..", "manuscript", "figures")
A <- readRDS(file.path(res, "analysis_data.rds"))
F <- readRDS(file.path(res, "fits.rds")); fits <- F$fits
I <- readRDS(file.path(res, "inference.rds"))
coh <- read.csv(file.path(res, "table_cohorts.csv"))
gen <- read.csv(file.path(res, "table_genes.csv"))
sth <- read.csv(file.path(res, "table_study_hr.csv"))
fr <- read.csv(file.path(res, "table_frailty.csv"))
lo <- read.csv(file.path(res, "table_loso.csv"))
tm <- read.csv(file.path(res, "table_methods.csv"))
fs <- read.csv(file.path(res, "sensitivity_folds.csv"))
th <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(),
                                       strip.background = element_rect(fill = "grey93"))
f2 <- function(x) formatC(x, format = "f", digits = 2)
f3 <- function(x) formatC(x, format = "f", digits = 3)
pfmt <- function(p) ifelse(p < 1e-4, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 4))
it <- function(g) paste0("\\textit{", g, "}")

## Table 4: cohorts
writeLines(c("\\begin{tabular}{@{}llrrrrr@{}}", "\\toprule",
  "Study & GEO/TCGA series & $n_k$ & Deaths & Censored (\\%) & Median OS (yr) & Median follow-up (yr) \\\\",
  "\\midrule",
  sprintf("%s & %s & %d & %d & %d & %s & %s \\\\", coh$Study, gsub("_", "\\\\_", coh$Dataset), coh$N, coh$Events,
          coh$Censoring, f2(coh$MedianOS), f2(coh$MedianFU)),
  "\\midrule",
  sprintf("Total & & %s & %d & %d & & \\\\", format(sum(coh$N), big.mark = ","), sum(coh$Events),
          round(100 * (1 - sum(coh$Events) / sum(coh$N)))),
  "\\bottomrule", "\\end{tabular}"), file.path(tab_dir, "tab_cohorts.tex"))

## Table 5: methods
dv <- function(x) ifelse(is.na(x), "--", as.character(x))
writeLines(c("\\begin{tabular}{@{}lrrrrr@{}}", "\\toprule",
  "& \\multicolumn{2}{c}{Genes selected} & \\multicolumn{2}{c}{Deviations} & \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "Method & min & 1-SE & min & 1-SE & CV time (s) \\\\", "\\midrule",
  sprintf("%s & %d & %d & %s & %s & %s \\\\", tm$label, tm$genes_min, tm$genes_1se, dv(tm$dev_min),
          dv(tm$dev_1se), formatC(tm$seconds, format = "f", digits = 1)),
  "\\bottomrule", "\\end{tabular}"), file.path(tab_dir, "tab_app_methods.tex"))

## Table 6: genes (main MCP-H model)
gen <- gen[order(gen$p_value), ]
writeLines(c("\\begin{tabular}{@{}lrrrrrrcc@{}}", "\\toprule",
  "Gene & $\\exp(\\hat\\alpha_j)$ & HR & 95\\% CI & $p$ & $p_{\\rm Holm}$ & Dev. & Boot. freq. & Methods \\\\",
  "\\midrule",
  sprintf("%s%s & %s & %s & (%s, %s) & %s & %s & %d & %s & %d/8 \\\\", it(gen$covariate),
          ifelse(gen$in_1se_model, "$^{\\ast}$", ""), f2(exp(gen$alpha_penalized)), f2(gen$HR),
          f2(gen$HR_lower), f2(gen$HR_upper), pfmt(gen$p_value), pfmt(gen$p_holm),
          gen$n_deviating_studies, f2(gen$selection_freq), gen$selected_by),
  "\\bottomrule", "\\end{tabular}"), file.path(tab_dir, "tab_app_genes.tex"))

## Table 7: frailty sensitivity
writeLines(c("\\begin{tabular}{@{}llll@{}}", "\\toprule",
  "Gene & Stratified Cox & Shared gamma frailty & Naive pooling \\\\", "\\midrule",
  sprintf("%s & %s & %s & %s \\\\", it(fr$gene), fr$stratified, fr$frailty, fr$pooled),
  "\\bottomrule", "\\end{tabular}"), file.path(tab_dir, "tab_app_frailty.tex"))

## Table 8: leave-one-study-out C-index
lm_ <- aggregate(cbind(C, n_genes) ~ method + rule, lo, mean)
wide <- function(rule) {
  d <- lo[lo$rule == rule, ]
  sapply(unique(lo$method), function(m) sapply(sort(unique(lo$held_out)), function(k)
    d$C[d$method == m & d$held_out == k]))
}
ms <- unique(lo$method)
rows <- c()
for (rule in c("min", "1se")) {
  W <- wide(rule)
  rows <- c(rows, sprintf("\\multicolumn{%d}{@{}l}{\\emph{%s}} \\\\", length(ms) + 1,
                          if (rule == "min") "Minimum rule" else "One-standard-error rule"),
            sprintf("Study %s held out & %s \\\\", sort(unique(lo$held_out)),
                    apply(W, 1, function(r) paste(f3(r), collapse = " & "))),
            sprintf("Mean C & %s \\\\", paste(f3(colMeans(W)), collapse = " & ")),
            sprintf("Mean no.\\ of genes & %s \\\\", paste(formatC(sapply(ms, function(m)
              lm_$n_genes[lm_$method == m & lm_$rule == rule]), format = "f", digits = 1), collapse = " & ")))
  if (rule == "min") rows <- c(rows, "\\midrule")
}
writeLines(c(sprintf("\\begin{tabular}{@{}l%s@{}}", strrep("r", length(ms))), "\\toprule",
  sprintf(" & %s \\\\", paste(ms, collapse = " & ")), "\\midrule", rows,
  "\\bottomrule", "\\end{tabular}"), file.path(tab_dir, "tab_app_loso.tex"))

## Figure 5: Kaplan-Meier by cohort
d <- A$dat
km <- survfit(Surv(time, status) ~ study, data = d)
kd <- data.frame(time = km$time, surv = km$surv,
                 study = rep(sub("study=", "Study ", names(km$strata)), km$strata))
kd <- rbind(data.frame(time = 0, surv = 1, study = unique(kd$study)), kd)
g <- ggplot(kd, aes(time, surv, colour = study)) + geom_step(linewidth = 0.6) +
  labs(x = "Years since diagnosis", y = "Overall survival", colour = NULL) +
  scale_colour_brewer(palette = "Dark2") + coord_cartesian(xlim = c(0, 12)) + th +
  theme(legend.position = c(0.85, 0.72))
ggsave(file.path(fig_dir, "app_km.pdf"), g, width = 5.5, height = 3.6)

## Figure 6: CV surface of MCP-H
cvf <- fits[["MCP-H"]][["min"]]$cv
cs <- do.call(rbind, lapply(seq_along(cvf$ratio.grid), function(r)
  data.frame(lambda = cvf$lambda.grid, cv = cvf$cvm[, r], se = cvf$cvsd[, r],
             ratio = paste0("r = ", cvf$ratio.grid[r]),
             df = cvf$df[, r])))
mn <- cvf$sel.min; s1 <- cvf$sel.1se
pts <- data.frame(lambda = cvf$lambda.grid[c(mn[1], s1[1])],
                  cv = c(cvf$cvm[mn[1], mn[2]], cvf$cvm[s1[1], s1[2]]),
                  ratio = paste0("r = ", cvf$ratio.grid[c(mn[2], s1[2])]),
                  rule = c("minimum", "one-SE"))
g <- ggplot(cs, aes(log(lambda), cv, colour = ratio)) + geom_line() + geom_point(size = 0.9) +
  geom_hline(yintercept = min(cvf$cvm, na.rm = TRUE) + cvf$cvsd[mn[1], mn[2]], linetype = 3) +
  geom_point(data = pts, aes(shape = rule), size = 3, colour = "black") +
  scale_shape_manual(values = c(minimum = 4, "one-SE" = 1), name = NULL) +
  scale_colour_brewer(palette = "Set1", name = expression(lambda[epsilon] / lambda[alpha])) +
  labs(x = expression(log(lambda[alpha])), y = "Cross-validated deviance (5 folds)") + th +
  theme(legend.position = "right")
ggsave(file.path(fig_dir, "app_cv_surface.pdf"), g, width = 6, height = 3.6)

## Figure 7: global and study-specific hazard ratios
sth$study <- paste("Study", sth$study)
gl <- gen[, c("covariate", "HR", "HR_lower", "HR_upper")]; gl$study <- "Global"
fd <- rbind(gl, sth[, c("covariate", "HR", "HR_lower", "HR_upper", "study")])
fd$type <- ifelse(fd$study == "Global", "Global effect", "Study-specific")
fd$dev <- ifelse(fd$study == "Global", FALSE,
                 sth$deviation_selected[match(paste(fd$covariate, fd$study), paste(sth$covariate, sth$study))])
fd$study <- factor(fd$study, levels = rev(c("Global", paste("Study", LETTERS[1:6]))))
fd$covariate <- factor(fd$covariate, levels = gen$covariate)
g <- ggplot(fd, aes(HR, study, colour = interaction(type, dev))) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.25) +
  geom_point(aes(shape = type), size = 2) +
  scale_shape_manual(values = c("Global effect" = 18, "Study-specific" = 16), name = NULL) +
  scale_colour_manual(values = c("Global effect.FALSE" = "#B2182B", "Study-specific.FALSE" = "grey45",
                                 "Study-specific.TRUE" = "#2166AC"),
                      labels = c("Global effect", "Study shares global effect", "Deviation selected"),
                      name = NULL) +
  scale_x_log10() + facet_wrap(~covariate, ncol = 3, scales = "free_x") +
  labs(x = "Hazard ratio per within-study SD (log scale), 95% Wald CI", y = NULL) + th +
  theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "app_forest.pdf"), g, width = 7.5, height = 2.3 * ceiling(nrow(gen) / 3) + 1)

## Figure 8: bootstrap selection frequencies
sf <- I$boot$selection; sf <- sf[order(-sf$selection_freq), ][1:20, ]
sf$covariate <- factor(sf$covariate, levels = rev(sf$covariate))
sf$sel <- ifelse(sf$covariate %in% gen$covariate, "Selected on full data", "Not selected")
g <- ggplot(sf, aes(selection_freq, covariate, fill = sel)) + geom_col(width = 0.7) +
  scale_fill_manual(values = c("Selected on full data" = "#B2182B", "Not selected" = "grey65"), name = NULL) +
  labs(x = sprintf("Selection frequency in %d within-study bootstrap resamples", I$boot$B_used), y = NULL) +
  th + theme(legend.position = "bottom", axis.text.y = element_text(face = "italic"))
ggsave(file.path(fig_dir, "app_stability.pdf"), g, width = 5.5, height = 4.8)

cat("done\n"); print(fs); print(lm_)
