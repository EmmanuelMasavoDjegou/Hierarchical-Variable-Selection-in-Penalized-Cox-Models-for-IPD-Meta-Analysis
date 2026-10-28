################################################################################
##  application/04_split_inference.R
##  Split-sample inference for the application (Section 6), added in the revision
##  because the Wald p-values of Table 11 are computed after selecting the genes on
##  the same data and their Holm adjustment covers only the selected genes.
##
##  Method  : multi-split p-values (Meinshausen, Meier and Buehlmann, 2009). For
##            b = 1..B, every cohort is split at random into two halves (within
##            cohort and event status). Hierarchical MCP is fitted on the first half
##            with the cross-validation of Section 3.2 (minimum rule); on the second
##            half the unpenalised stratified Cox model with the selected global
##            effects and deviations is refitted and the Wald p-value of each
##            selected global effect is multiplied by the number of selected genes
##            (genes not selected get 1). The adjusted p-value of gene j is
##            P_j = min{1, (1 - log g_min) inf_{g in (g_min, 1)} q_g(P_j^(b) / g)},
##            g_min = 0.05, which controls the family-wise error rate over all 500
##            genes asymptotically under the screening and sparsity conditions of
##            that paper; there the per-split p-values come from least squares, here
##            from the partial likelihood, so the guarantee is asymptotic only.
##            Covariates are standardised within cohort separately in each half.
##  Input   : results/analysis_data.rds (01_prepare_data.R)
##  Output  : results/table_split.csv, results/split_raw.rds,
##            ../output/tables/tab_app_split.tex
##  Seeds   : split b uses set.seed(4000 + b); CV folds within the selection half use
##            seed 5000 + b.     Runtime: about 20 minutes on one core for B = 50.
################################################################################
here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) {   # folder of this script
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))     # works under Rscript
  if (length(f)) dirname(normalizePath(f[1])) else "." })                          # from any directory
source(file.path(here, "..", "simulation", "00_setup.R"))
args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args)) as.integer(args[1]) else 50L
GMIN <- 0.05
out <- file.path(here, "results")
tab_dir <- file.path(here, "..", "output", "tables"); dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
A <- readRDS(file.path(out, "analysis_data.rds"))
dat <- A$dat; genes <- A$genes; Xraw <- A$X; colnames(Xraw) <- genes; p <- length(genes)

one_split <- function(b) {
  set.seed(4000 + b)
  half <- integer(nrow(dat))
  for (s in unique(dat$study)) for (e in 0:1) {
    idx <- which(dat$study == s & dat$status == e)
    half[idx] <- sample(rep_len(1:2, length(idx)))
  }
  i1 <- which(half == 1); i2 <- which(half == 2)
  X1 <- std_within(Xraw[i1, ], dat$study[i1]); X2 <- std_within(Xraw[i2, ], dat$study[i2])
  fid <- hmc_folds(dat$study[i1], dat$status[i1], TUNE$nfolds, seed = 5000 + b)
  cvf <- cv.hmcox(X1, dat$time[i1], dat$status[i1], dat$study[i1], penalty = "MCP",
                  ratio = TUNE$ratio, foldid = fid, nlambda = TUNE$nlambda,
                  lambda.min.ratio = TUNE$lambda.min.ratio, tol = TUNE$tol,
                  dfmax = TUNE$dfmax, standardize = FALSE, rule = "min")
  cf <- coef(cvf, rule = "min")
  S <- which(cf$alpha != 0); Dv <- which(cf$eps != 0, arr.ind = TRUE)
  praw <- rep(1, p); hr <- rep(NA_real_, p)
  if (length(S)) {
    lev <- rownames(cf$eps)
    cols <- lapply(S, function(j) X2[, j]); lab <- paste0("a", S)
    if (length(Dv)) for (r in seq_len(nrow(Dv))) {
      cols[[length(cols) + 1]] <- X2[, Dv[r, 2]] * (dat$study[i2] == lev[Dv[r, 1]])
      lab <- c(lab, sprintf("e%d_%d", Dv[r, 1], Dv[r, 2]))
    }
    Z <- do.call(cbind, cols); colnames(Z) <- lab                  # named columns: a single
    df2 <- data.frame(time = dat$time[i2], status = dat$status[i2], sid = dat$study[i2], Z)  # column keeps its name
    fml <- as.formula(paste("Surv(time, status) ~", paste(lab, collapse = " + "), "+ strata(sid)"))
    fit <- tryCatch(coxph(fml, data = df2, ties = "breslow"), error = function(e) NULL)
    if (!is.null(fit)) {
      sm <- summary(fit)$coefficients
      la <- lab[seq_along(S)]; ok <- la %in% rownames(sm)      # aliased terms count as p = 1
      pa <- rep(1, length(S)); pa[ok] <- sm[la[ok], "Pr(>|z|)"]; pa[is.na(pa)] <- 1
      praw[S] <- pmin(1, pa * length(S)); hr[S[ok]] <- exp(sm[la[ok], "coef"])
    }
  }
  list(b = b, selected = S, n_dev = nrow(Dv), p_adj_split = praw, hr = hr)
}

f_raw <- file.path(out, "split_raw.rds")
res <- if (file.exists(f_raw)) readRDS(f_raw) else list()
todo <- setdiff(seq_len(B), vapply(res, `[[`, 0, "b"))
for (b in todo) { res[[length(res) + 1]] <- one_split(b); saveRDS(res, f_raw)
  cat(sprintf("split %d: %d genes selected\n", b, length(res[[length(res)]]$selected))) }
res <- res[vapply(res, `[[`, 0, "b") %in% seq_len(B)]
P <- sapply(res, `[[`, "p_adj_split")                       # p x B
gam <- seq(GMIN, 1, by = 0.01)
Qg <- sapply(gam, function(g) pmin(1, apply(P / g, 1, quantile, probs = g, type = 7)))
P_ms <- pmin(1, (1 - log(GMIN)) * apply(Qg, 1, min))
selfreq <- rowMeans(sapply(res, function(z) seq_len(p) %in% z$selected))
tab <- data.frame(gene = genes, split_selection_freq = selfreq, multisplit_p = P_ms,
                  median_split_p = apply(P, 1, median))
tab <- tab[order(tab$multisplit_p, -tab$split_selection_freq), ]
write.csv(tab, file.path(out, "table_split.csv"), row.names = FALSE)

main6 <- c("GMPR", "CXCL9", "CRYAB", "OSR2", "GAS1", "IGF2")
show <- unique(c(main6, head(tab$gene[tab$split_selection_freq > 0], 10)))
show <- tab[tab$gene %in% show, ]
show <- show[order(show$multisplit_p, -show$split_selection_freq), ]
pf <- function(v) ifelse(v >= 0.9995, "1", ifelse(v < 0.001, formatC(v, format = "e", digits = 1), formatC(v, format = "f", digits = 3)))
rows <- sprintf("\\textit{%s}%s & %s & %s & %s \\\\", show$gene, ifelse(show$gene %in% main6, "$^{\\dagger}$", ""),
                formatC(show$split_selection_freq, format = "f", digits = 2), pf(show$median_split_p), pf(show$multisplit_p))
writeLines(c("\\begin{tabular}{@{}lrrr@{}}", "\\toprule",
             "Gene & Selected in first half & Median split $p$ & Multi-split $p$ \\\\", "\\midrule",
             rows, "\\bottomrule", "\\end{tabular}"), file.path(tab_dir, "tab_app_split.tex"))
cat(sprintf("B = %d splits; median genes selected per split = %g; genes with multi-split p < 0.05: %s\n",
            length(res), median(sapply(res, function(z) length(z$selected))),
            paste(tab$gene[tab$multisplit_p < 0.05], collapse = ", ")))
