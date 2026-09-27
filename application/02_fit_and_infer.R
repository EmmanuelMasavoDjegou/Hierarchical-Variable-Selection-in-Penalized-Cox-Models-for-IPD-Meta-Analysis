################################################################################
##  application/02_fit_and_infer.R
##  Real-data analysis of Section 6: six ovarian cancer cohorts, p = 500 genes.
##
##  Input   : results/analysis_data.rds  (from 01_prepare_data.R)
##  Output  : results/fits.rds                     all cross-validated fits
##            results/table_methods.csv            Table 5  (selected-set sizes)
##            results/table_genes.csv              Table 6  (HR, 95% CI, p, stability)
##            results/table_study_hr.csv           study-specific HRs (Figure 7)
##            results/table_frailty.csv            Table 7  (frailty sensitivity)
##            results/table_loso.csv               Table 8  (cross-study validation)
##            results/sensitivity_folds.csv        V = 5 vs V = 10 selection
##  Seeds   : fold assignment seed 2026; bootstrap seed 7; LOSO seed 11.
##  Runtime : about 25 minutes on one core.
################################################################################
here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) ".")
source(file.path(here, "..", "simulation", "00_setup.R"))       # TUNE, fit_method, std_within
out <- file.path(here, "results")
A <- readRDS(file.path(out, "analysis_data.rds"))
dat <- A$dat; genes <- A$genes
X <- std_within(A$X, dat$study); colnames(X) <- genes
RULE <- "min"   # main model: estimates study-specific effects (Section 3.4); 1-SE reported alongside

## 1. all eight methods, same folds --------------------------------------------
fid <- hmc_folds(dat$study, dat$status, TUNE$nfolds, seed = 2026)
fits <- list(); secs <- c()
for (m in METHODS) {
  t0 <- proc.time()[3]
  fits[[m]] <- fit_method(m, X, dat$time, dat$status, dat$study, fid)
  secs[m] <- proc.time()[3] - t0
  cat(m, round(secs[m], 1), "s\n")
}
saveRDS(list(fits = fits, secs = secs, foldid = fid), file.path(out, "fits.rds"))

sel <- function(m, rule) which(fits[[m]][[rule]]$alpha != 0)
tab_m <- do.call(rbind, lapply(METHODS, function(m) {
  data.frame(method = m, label = METHOD_LABELS[[m]],
             genes_min = length(sel(m, "min")), genes_1se = length(sel(m, "1se")),
             dev_min = if (is.null(fits[[m]]$min$eps)) NA else sum(fits[[m]]$min$eps != 0),
             dev_1se = if (is.null(fits[[m]][["1se"]]$eps)) NA else sum(fits[[m]][["1se"]]$eps != 0),
             ratio_1se = fits[[m]][["1se"]]$extra[["ratio"]],
             overlap_with_MCPH = length(intersect(sel(m, RULE), sel("MCP-H", RULE))),
             seconds = round(secs[[m]], 1))
}))
write.csv(tab_m, file.path(out, "table_methods.csv"), row.names = FALSE)
print(tab_m)

## 2. inference for the proposed MCP-H fit -------------------------------------
cvf <- fits[["MCP-H"]][[RULE]]$cv
inf <- hmcox_inference(cvf, rule = RULE)
bt <- hmcox_boot(cvf, B = 100, seed = 7, rule = RULE)
g <- merge(inf$global, bt$selection[, c("covariate", "selection_freq", "boot_lower", "boot_upper")],
           by = "covariate", sort = FALSE)
g$in_1se_model <- g$covariate %in% genes[sel("MCP-H", "1se")]
g$selected_by <- sapply(g$covariate, function(z)
  sum(sapply(METHODS, function(m) z %in% genes[sel(m, RULE)])))
g$p_holm <- p.adjust(g$p_value, "holm")
g <- g[order(g$p_value), ]
write.csv(g, file.path(out, "table_genes.csv"), row.names = FALSE)
st <- inf$study
write.csv(st, file.path(out, "table_study_hr.csv"), row.names = FALSE)
dev_freq <- bt$deviation_freq; colnames(dev_freq) <- genes; rownames(dev_freq) <- cvf$data$levels
saveRDS(list(inf = inf, boot = bt, dev_freq = dev_freq), file.path(out, "inference.rds"))
print(g)

## 3. frailty sensitivity: selected genes, stratified vs shared gamma frailty ---
Xs <- X[, g$covariate, drop = FALSE]
f_str <- coxph(Surv(dat$time, dat$status) ~ Xs + strata(dat$study), ties = "breslow")
f_fr  <- coxph(Surv(dat$time, dat$status) ~ Xs + frailty(dat$study, distribution = "gamma"))
f_pool <- coxph(Surv(dat$time, dat$status) ~ Xs, ties = "breslow")
ci <- function(f) { b <- coef(f)[seq_len(ncol(Xs))]; s <- sqrt(diag(vcov(f)))[seq_len(ncol(Xs))]
  sprintf("%.2f (%.2f, %.2f)", exp(b), exp(b - 1.96 * s), exp(b + 1.96 * s)) }
tab_fr <- data.frame(gene = g$covariate, stratified = ci(f_str), frailty = ci(f_fr), pooled = ci(f_pool))
theta_fr <- f_fr$history[[1]]$theta
write.csv(tab_fr, file.path(out, "table_frailty.csv"), row.names = FALSE)
writeLines(sprintf("%.4f", theta_fr), file.path(out, "frailty_variance.txt"))
print(tab_fr); cat("frailty variance:", theta_fr, "\n")

## 4. number of folds: V = 5 vs V = 10 ----------------------------------------
cv10 <- cv.hmcox(X, dat$time, dat$status, dat$study, penalty = "MCP", ratio = TUNE$ratio,
                 nfolds = 10, seed = 2026, nlambda = TUNE$nlambda,
                 lambda.min.ratio = TUNE$lambda.min.ratio, tol = TUNE$tol,
                 dfmax = TUNE$dfmax, standardize = FALSE)
fold_sens <- do.call(rbind, lapply(c("min", "1se"), function(rule) {
  s5 <- genes[sel("MCP-H", rule)]; s10 <- genes[coef(cv10, rule = rule)$alpha != 0]
  data.frame(rule = rule, V = c(5, 10), n_selected = c(length(s5), length(s10)),
             overlap = length(intersect(s5, s10)),
             n_dev = c(sum(fits[["MCP-H"]][[rule]]$eps != 0), sum(coef(cv10, rule = rule)$eps != 0)),
             genes = c(paste(s5, collapse = " "), paste(s10, collapse = " ")))
}))
write.csv(fold_sens, file.path(out, "sensitivity_folds.csv"), row.names = FALSE)
print(fold_sens)

## 5. leave-one-study-out validation of the global risk score -------------------
loso <- list()
for (k in unique(dat$study)) {
  tr <- dat$study != k; te <- !tr
  ftr <- hmc_folds(dat$study[tr], dat$status[tr], TUNE$nfolds, seed = 11)
  for (m in c("MCP-H", "SCAD-H", "MCP-S", "MCP-P", "LASSO-P", "ENet-P")) {
    ff <- fit_method(m, X[tr, ], dat$time[tr], dat$status[tr], dat$study[tr], ftr)
    for (rule in c("min", "1se")) {
      lp <- drop(X[te, , drop = FALSE] %*% ff[[rule]]$alpha)
      cc <- if (all(lp == 0)) 0.5 else
        concordance(Surv(dat$time[te], dat$status[te]) ~ lp, reverse = TRUE)$concordance
      loso[[length(loso) + 1]] <- data.frame(held_out = k, method = m, rule = rule,
                                             n_genes = sum(ff[[rule]]$alpha != 0), C = cc)
    }
  }
  cat("LOSO", k, "done\n")
}
loso <- do.call(rbind, loso)
write.csv(loso, file.path(out, "table_loso.csv"), row.names = FALSE)
print(aggregate(cbind(C, n_genes) ~ method + rule, loso, mean))
