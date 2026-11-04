################################################################################
##  simulation/09_simulation_figures.R
##  Figures of Sections 5.2-5.4, drawn from results/sim_summary.csv only (no raw
##  replicate files needed), so they can be redrawn after 02_make_tables_figures.R,
##  which also calls this script at its end.
##
##  Two designs, each drawn separately for the one-standard-error rule (_1se) and the
##  minimum rule (_min):
##   * value tables: one row per scenario, one column per method, the Monte Carlo
##     mean printed in every cell; shading runs from lighter (worse) to darker
##     (better) on a fixed scale per metric, so the two rules are comparable.
##       output/figures/sim_mcc_{1se,min}.png          MCC of global selection
##       output/figures/sim_fdr_tpr_{1se,min}.png      FDR and TPR of global selection
##       output/figures/sim_mse_theta_{1se,min}.png    MSE of study-specific effects (x100)
##       output/figures/sim_deviations_{1se,min}.png   TPR and FDR of the deviations
##   * paired differences of MCC with +-2 Monte Carlo standard errors:
##       output/figures/sim_mcc_diff_{1se,min}.png     hierarchical - stratified and
##                                                     stratified - naive pooling
##  Usage: Rscript simulation/09_simulation_figures.R   (from the repository root)
################################################################################
suppressMessages(library(ggplot2))
here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) {   # folder of this script
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))     # works under Rscript
  if (length(f)) dirname(normalizePath(f[1])) else "." })                          # from any directory
fig_dir <- file.path(here, "..", "output", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
S <- read.csv(file.path(here, "results", "sim_summary.csv"))

SCEN <- c("Base", "Independent covariates", "Strong correlation", "No heterogeneity",
          "Strong heterogeneity", "Dense heterogeneity", "Heavy censoring", "Weak signals",
          "p = 200", "p = 500", "Ten studies", "Unequal study sizes", "Small studies, p > n_k")
S$scen <- factor(S$scenario, levels = rev(SCEN))
S$penalty <- factor(sub("-(H|S|P)$", "", S$method), levels = c("MCP", "SCAD", "LASSO", "ENet"),
                    labels = c("MCP", "SCAD", "Lasso", "Elastic net"))
S$structure <- factor(sub("^.*-", "", S$method), levels = c("H", "S", "P"),
                      labels = c("Hier.", "Strat.", "Pooled"))
RULES <- c("1se" = "one-standard-error rule", min = "minimum rule")

th <- theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey95", colour = NA),
        strip.text = element_text(face = "bold"), axis.ticks = element_blank(),
        axis.title = element_blank(), legend.position = "bottom",
        legend.title = element_text(size = 8), legend.text = element_text(size = 8),
        panel.spacing.x = unit(6, "pt"), strip.placement = "outside")

## Fill is a "goodness" score in [0, 1] on a fixed range per metric (reversed when lower
## is better), so lighter always means worse and darker better; the printed number is
## the Monte Carlo mean itself.
good <- function(v, lo, hi, higher_better, log = FALSE) {
  if (log) { v <- log10(v); lo <- log10(lo); hi <- log10(hi) }
  g <- pmin(1, pmax(0, (v - lo) / (hi - lo)))
  if (higher_better) g else 1 - g
}
fill_scale <- scale_fill_gradient(low = "#eef3fb", high = "#0d3d78", limits = c(0, 1),
                                  breaks = c(0, 1), labels = c("worse", "better"), name = NULL,
                                  guide = guide_colourbar(barwidth = unit(6, "cm"), barheight = unit(0.25, "cm")))
tile_layers <- function()
  list(geom_tile(colour = "white", linewidth = 0.6),
       geom_text(aes(label = lab, colour = g > 0.55), size = 2.5),
       scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"), guide = "none"),
       fill_scale, scale_x_discrete(position = "top"), th)

## one metric, all 12 methods
value_table <- function(rule, v, lo, hi, higher, digits = 2, mult = 1, log = FALSE) {
  d <- S[S$rule == rule, ]
  d$g <- good(d[[v]], lo, hi, higher, log)
  d$lab <- formatC(d[[v]] * mult, format = "f", digits = digits)
  ggplot(d, aes(structure, scen, fill = g)) + tile_layers() + facet_wrap(~ penalty, nrow = 1)
}
## two metrics stacked (rows), all 12 methods
value_table2 <- function(rule, specs) {
  d <- do.call(rbind, lapply(specs, function(sp) {
    z <- S[S$rule == rule, ]
    data.frame(z[, c("scen", "penalty", "structure")], metric = sp$name,
               g = good(z[[sp$v]], sp$lo, sp$hi, sp$higher),
               lab = formatC(z[[sp$v]], format = "f", digits = 2))
  }))
  d$metric <- factor(d$metric, levels = sapply(specs, `[[`, "name"))
  ggplot(d, aes(structure, scen, fill = g)) + tile_layers() +
    facet_grid(metric ~ penalty, switch = "y") +
    theme(strip.text.y.left = element_text(angle = 90), strip.background.y = element_blank())
}

for (rule in names(RULES)) {
  ggsave(file.path(fig_dir, sprintf("sim_mcc_%s.png", rule)),
         value_table(rule, "MCC", 0.25, 1, TRUE), width = 7.5, height = 4.6, dpi = 300)
  ggsave(file.path(fig_dir, sprintf("sim_fdr_tpr_%s.png", rule)),
         value_table2(rule, list(list(name = "False discovery rate", v = "FDR", lo = 0, hi = 0.9, higher = FALSE),
                                 list(name = "True positive rate", v = "TPR", lo = 0.5, hi = 1, higher = TRUE))),
         width = 7.5, height = 8.4, dpi = 300)
  ggsave(file.path(fig_dir, sprintf("sim_mse_theta_%s.png", rule)),
         value_table(rule, "MSE_theta", 0.00015, 0.02, FALSE, digits = 2, mult = 100, log = TRUE),
         width = 7.5, height = 4.6, dpi = 300)

  ## deviations: hierarchical methods only, no "No heterogeneity" (it has no true deviations)
  H <- S[S$rule == rule & S$structure == "Hier." & S$scenario != "No heterogeneity", ]
  D <- rbind(data.frame(H[, c("scen", "penalty")], metric = "True positive rate",
                        g = good(H$TPR_dev, 0, 0.75, TRUE), lab = formatC(H$TPR_dev, format = "f", digits = 2)),
             data.frame(H[, c("scen", "penalty")], metric = "False discovery rate",
                        g = good(H$FDR_dev, 0, 1, FALSE), lab = formatC(H$FDR_dev, format = "f", digits = 2)))
  D$metric <- factor(D$metric, levels = c("True positive rate", "False discovery rate"))
  D$scen <- droplevels(D$scen)
  ggsave(file.path(fig_dir, sprintf("sim_deviations_%s.png", rule)),
         ggplot(D, aes(penalty, scen, fill = g)) + tile_layers() + facet_wrap(~ metric, nrow = 1),
         width = 6.5, height = 4.4, dpi = 300)

  ## paired differences of MCC; SE computed as for independent estimators (conservative)
  z <- S[S$rule == rule, c("scen", "penalty", "structure", "MCC", "MCC_se")]
  pick <- function(s) setNames(z[z$structure == s, c("scen", "penalty", "MCC", "MCC_se")],
                               c("scen", "penalty", paste0("m", s), paste0("s", s)))
  w <- merge(merge(pick("Hier."), pick("Strat.")), pick("Pooled"))
  C <- rbind(data.frame(w[, 1:2], comparison = "Effect of the hierarchy\n(hierarchical minus stratified)",
                        diff = w$mHier. - w$mStrat., se = sqrt(w$sHier.^2 + w$sStrat.^2)),
             data.frame(w[, 1:2], comparison = "Effect of stratification\n(stratified minus naive pooling)",
                        diff = w$mStrat. - w$mPooled, se = sqrt(w$sStrat.^2 + w$sPooled^2)))
  C$comparison <- factor(C$comparison, levels = unique(C$comparison))
  C$sig <- factor(ifelse(abs(C$diff) > 2 * C$se, "Difference > 2 Monte Carlo SE", "Within 2 Monte Carlo SE"),
                  levels = c("Difference > 2 Monte Carlo SE", "Within 2 Monte Carlo SE"))
  pC <- ggplot(C, aes(diff, scen)) + geom_vline(xintercept = 0, colour = "grey40") +
    geom_errorbarh(aes(xmin = diff - 2 * se, xmax = diff + 2 * se), height = 0, colour = "grey55") +
    geom_point(aes(shape = sig), size = 1.9, colour = "#2a78d6", fill = "white", stroke = 0.8) +
    scale_shape_manual(values = c(16, 21), name = NULL, drop = FALSE) +
    facet_grid(comparison ~ penalty) +
    labs(x = sprintf("Difference in MCC, %s (positive favours the first model)", RULES[[rule]])) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey93"),
          strip.background.x = element_rect(fill = "grey95", colour = NA), strip.background.y = element_blank(),
          strip.text = element_text(face = "bold"), strip.text.y = element_text(angle = 0, hjust = 0),
          axis.title.y = element_blank(), legend.position = "bottom")
  ggsave(file.path(fig_dir, sprintf("sim_mcc_diff_%s.png", rule)), pC, width = 7.5, height = 6.6, dpi = 300)
}
cat("simulation figures written to output/figures/\n")
