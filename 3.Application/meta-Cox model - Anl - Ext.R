################################################################################
##                                                                            ##
##   SECTION 7: REAL DATA APPLICATION  (FINAL VERSION WITH BOOTSTRAP CIs)    ##
##   "Additive Hierarchical Variable Selection for Nonconvex Penalized        ##
##    Cox Models in Individual-Participant-Data Meta-Analysis"                ##
##                                                                            ##
##   Dataset:  curatedOvarianData                                             ##
##     Multi-cohort high-grade serous ovarian carcinoma microarray data.      ##
##     Outcome: overall survival (days -> years).                             ##
##     Ganzfried et al. (2013).                                               ##
##                                                                            ##
##   NEW in this version:                                                     ##
##     - bootstrap_mcp_hier()  : stratified bootstrap (B resamples)          ##
##       returning SD, 95% CI for every alpha_hat_j and eps_hat_kj           ##
##     - Table7_2 extended with SE and 95% CI columns for alpha_hat          ##
##     - Table7_4 (new): study-specific deviation estimates with 95% CIs     ##
##     - Figure 7.3 (forest plot) updated with CI error bars                 ##
##                                                                            ##
##   OUTPUT FILES (written to Section7_outputs/):                             ##
##     Figures : Fig7_1_cohort_bar.png      Fig7_2_km_curves.png             ##
##               Fig7_3_forest_CI.png       Fig7_4_deviation_tile.png        ##
##               Fig7_5_mcp_vs_scad.png     Fig7_6_selection_sizes.png       ##
##     Tables  : Table7_1_cohort_chars.csv  Table7_2_biomarkers_CI.csv       ##
##               Table7_3_method_summary.csv Table7_4_deviations_CI.csv      ##
##                                                                            ##
################################################################################


################################################################################
##  0.  PACKAGES                                                              ##
################################################################################

## Install once (uncomment on first run):
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("curatedOvarianData", "Biobase"))
# install.packages(c("survival","ncvreg","glmnet","MASS","ggplot2","tidyr",
#                    "dplyr","patchwork","ggrepel","RColorBrewer","scales",
#                    "survminer","stringr","gridExtra","parallel","pbapply"))

suppressPackageStartupMessages({
  library(survival)
  library(ncvreg)
  library(glmnet)
  library(MASS)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(patchwork)
  library(ggrepel)
  library(RColorBrewer)
  library(scales)
  library(survminer)
  library(stringr)
  library(gridExtra)
  library(parallel)
  library(Biobase)
  library(curatedOvarianData)
})

## pbapply is optional (progress bar for bootstrap); graceful fallback
has_pbapply <- requireNamespace("pbapply", quietly = TRUE)

out_dir <- "Section7_outputs"
dir.create(out_dir, showWarnings = FALSE)
fig <- function(name) file.path(out_dir, name)

theme_paper <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1),
      plot.subtitle    = element_text(colour = "grey40", size = base_size - 1),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text       = element_text(face = "bold"),
      legend.position  = "right",
      panel.grid.minor = element_blank()
    )
}

METHOD_COLORS <- c(
  "MCP (Hier.)"    = "#E6851E",
  "SCAD (Hier.)"   = "#C0392B",
  "MCP (Pooled)"   = "#7F7F7F",
  "SCAD (Pooled)"  = "#BDC3C7",
  "LASSO (Pooled)" = "#2980B9",
  "ENet (Pooled)"  = "#8E44AD"
)


################################################################################
##  1.  ESTIMATION HELPERS                                                    ##
################################################################################

build_pooled_data <- function(data_list, offset_list = NULL) {
  K <- length(data_list)
  X_list <- y_list <- off_list <- vector("list", K)
  for (k in seq_len(K)) {
    dk          <- data_list[[k]]
    Xk          <- as.matrix(dk[, -(1:3), drop = FALSE])
    X_list[[k]] <- Xk
    y_list[[k]] <- cbind(time = dk$time, status = dk$status)
    off_list[[k]] <- if (is.null(offset_list)) rep(0, nrow(dk))
    else offset_list[[k]]
  }
  list(
    X      = do.call(rbind, X_list),
    y      = Surv(do.call(c, lapply(y_list, `[`, , "time")),
                  do.call(c, lapply(y_list, `[`, , "status"))),
    offset = do.call(c, off_list)
  )
}

cv_select_lambda <- function(X, y, offset, pen, alpha_g = 1, nfolds = 5L) {
  if (pen %in% c("SCAD", "MCP")) {
    cv <- tryCatch(
      cv.ncvsurv(X, y, penalty = pen, offset = offset,
                 nfolds = nfolds, warn = FALSE),
      error = function(e) NULL)
    if (!is.null(cv)) return(cv$lambda.min)
    fit <- ncvsurv(X, y, penalty = pen, offset = offset)
    nnz <- colSums(fit$beta != 0)
    return(fit$lambda[which.min(abs(nnz - 5L))])
  } else {
    cv <- tryCatch(
      cv.glmnet(X, y, family = "cox", alpha = alpha_g, offset = offset,
                nfolds = nfolds, type.measure = "C"),
      error = function(e) NULL)
    if (!is.null(cv)) return(cv$lambda.min)
    fit <- glmnet(X, y, family = "cox", alpha = alpha_g, offset = offset)
    return(fit$lambda[max(1L, length(fit$lambda) %/% 2L)])
  }
}

#' Fit the hierarchical penalised IPD meta-Cox model (Algorithm 1).
#'
#' @param data_list    List of K study data frames (study_id, time, status, X..)
#' @param penalty      "mcp" | "scad" | "lasso" | "enet"
#' @param lambda_alpha Fixed lambda for global alpha (NULL = CV each iteration)
#' @param lambda_eps   Fixed lambda for deviations  (NULL = CV each iteration)
#' @param max_iter     Maximum outer alternating iterations
#' @param conv_tol     Convergence tolerance on max|delta theta|
#' @param active_tol   Hard-zero threshold after fitting
#' @param nfolds_cv    CV folds
#' @param verbose      Print progress
#' @return List: alpha_hat (p), eps_hat (K x p), theta_hat (K x p),
#'         active_set, lambda_alpha_used, lambda_eps_used, n_iter
fit_hierarchical_meta_cox <- function(
    data_list,
    penalty      = "mcp",
    lambda_alpha = NULL,
    lambda_eps   = NULL,
    max_iter     = 40L,
    conv_tol     = 1e-4,
    active_tol   = 1e-6,
    nfolds_cv    = 5L,
    verbose      = FALSE
) {
  K         <- length(data_list)
  p         <- ncol(data_list[[1]]) - 3L
  pen_upper <- toupper(penalty)
  use_ncv   <- pen_upper %in% c("SCAD", "MCP")
  alpha_g   <- ifelse(pen_upper == "ENET", 0.5, 1.0)
  
  alpha_hat  <- rep(0, p)
  eps_hat    <- matrix(0, K, p)
  theta_hat  <- matrix(0, K, p)
  lam_a_used <- lam_e_used <- NA_real_
  
  if (verbose)
    cat(sprintf("  [%s]  lambda_alpha=%s  lambda_eps=%s\n",
                pen_upper,
                ifelse(is.null(lambda_alpha), "CV", sprintf("%.5f", lambda_alpha)),
                ifelse(is.null(lambda_eps),   "CV", sprintf("%.5f", lambda_eps))))
  
  for (iter in seq_len(max_iter)) {
    theta_old <- theta_hat
    
    ## Stage 1: update global alpha
    off_eps <- lapply(seq_len(K), function(k)
      as.numeric(as.matrix(data_list[[k]][, -(1:3)]) %*% eps_hat[k, ]))
    pd    <- build_pooled_data(data_list, off_eps)
    lam_a <- if (is.null(lambda_alpha))
      cv_select_lambda(pd$X, pd$y, pd$offset, pen_upper, alpha_g, nfolds_cv)
    else lambda_alpha
    lam_a_used <- lam_a
    
    if (use_ncv) {
      fa        <- ncvsurv(pd$X, pd$y, penalty = pen_upper,
                           offset = pd$offset, lambda = lam_a, warn = FALSE)
      alpha_hat <- as.numeric(coef(fa))
    } else {
      fa        <- glmnet(pd$X, pd$y, family = "cox", alpha = alpha_g,
                          offset = pd$offset, lambda = lam_a)
      alpha_hat <- as.numeric(coef(fa))
    }
    alpha_hat[abs(alpha_hat) < active_tol] <- 0
    active_set <- which(alpha_hat != 0)
    
    ## Stage 2: update study-specific epsilon (restricted to active_set)
    eps_hat <- matrix(0, K, p)
    if (length(active_set) > 0L) {
      for (k in seq_len(K)) {
        dk     <- data_list[[k]]
        Xk     <- as.matrix(dk[, -(1:3)])
        nk     <- nrow(dk)
        surv_k <- Surv(dk$time, dk$status)
        off_a  <- as.numeric(Xk %*% alpha_hat)
        Xk_act <- Xk[, active_set, drop = FALSE]
        
        if (ncol(Xk_act) == 1L) {
          df_tmp <- data.frame(time = dk$time, status = dk$status,
                               v1 = Xk_act[, 1])
          fm <- tryCatch(
            coxph(Surv(time, status) ~ v1 + offset(off_a), data = df_tmp,
                  control = coxph.control(iter.max = 50L)),
            error = function(e) NULL)
          if (!is.null(fm)) {
            est <- coef(fm)["v1"]
            eps_hat[k, active_set] <- ifelse(abs(est) < active_tol, 0, est)
          }
        } else if (use_ncv) {
          lam_e <- if (is.null(lambda_eps))
            tryCatch(
              cv_select_lambda(Xk_act, surv_k, off_a, pen_upper,
                               nfolds = min(nfolds_cv, max(3L, nk %/% 20L))),
              error = function(e) lam_a * 0.5)
          else lambda_eps
          lam_e_used <- lam_e
          fe <- tryCatch(
            ncvsurv(Xk_act, surv_k, penalty = pen_upper,
                    offset = off_a, lambda = lam_e, warn = FALSE),
            error = function(e) NULL)
          if (!is.null(fe)) {
            ce <- as.numeric(coef(fe))
            ce[abs(ce) < active_tol] <- 0
            eps_hat[k, active_set] <- ce
          }
        } else {
          lam_e <- if (is.null(lambda_eps))
            tryCatch(
              cv_select_lambda(Xk_act, surv_k, off_a, pen_upper, alpha_g,
                               min(nfolds_cv, max(3L, nk %/% 20L))),
              error = function(e) lam_a * 0.5)
          else lambda_eps
          lam_e_used <- lam_e
          fe <- tryCatch(
            glmnet(Xk_act, surv_k, family = "cox", alpha = alpha_g,
                   offset = off_a, lambda = lam_e),
            error = function(e) NULL)
          if (!is.null(fe)) {
            ce <- as.numeric(coef(fe))
            ce[abs(ce) < active_tol] <- 0
            eps_hat[k, active_set] <- ce
          }
        }
      }
    }
    
    theta_hat  <- sweep(eps_hat, 2L, alpha_hat, "+")
    max_change <- max(abs(theta_hat - theta_old))
    if (verbose && (iter == 1L || iter %% 5L == 0L))
      cat(sprintf("    iter %2d | active=%3d | delta=%.2e\n",
                  iter, length(active_set), max_change))
    if (max_change < conv_tol) break
  }
  
  list(alpha_hat = alpha_hat, eps_hat = eps_hat, theta_hat = theta_hat,
       active_set        = which(alpha_hat != 0),
       lambda_alpha_used = lam_a_used,
       lambda_eps_used   = lam_e_used,
       n_iter            = iter)
}

fit_pooled_penalized <- function(data_list, method = "lasso",
                                 lambda = NULL, nfolds = 5L) {
  combined  <- do.call(rbind, data_list)
  X_pool    <- as.matrix(combined[, -(1:3)])
  surv_pool <- Surv(combined$time, combined$status)
  K         <- length(data_list)
  alpha_g   <- ifelse(method == "enet", 0.5, 1.0)
  if (is.null(lambda)) {
    cv     <- cv.glmnet(X_pool, surv_pool, family = "cox", alpha = alpha_g,
                        nfolds = nfolds, type.measure = "C")
    lambda <- cv$lambda.min
  }
  fit      <- glmnet(X_pool, surv_pool, family = "cox",
                     alpha = alpha_g, lambda = lambda)
  beta_hat <- as.numeric(coef(fit))
  list(alpha_hat  = beta_hat,
       eps_hat    = matrix(0, K, length(beta_hat)),
       theta_hat  = matrix(rep(beta_hat, K), nrow = K, byrow = TRUE),
       active_set = which(abs(beta_hat) > 1e-8),
       lambda     = lambda)
}

fit_ncvreg_pooled <- function(data_list, penalty = "MCP",
                              lambda = NULL, nfolds = 5L) {
  combined  <- do.call(rbind, data_list)
  X_pool    <- as.matrix(combined[, -(1:3)])
  surv_pool <- Surv(combined$time, combined$status)
  K         <- length(data_list)
  if (is.null(lambda)) {
    cv <- tryCatch(
      cv.ncvsurv(X_pool, surv_pool, penalty = penalty, nfolds = nfolds),
      error = function(e) NULL)
    lambda <- if (!is.null(cv)) cv$lambda.min else 0.05
  }
  fit      <- ncvsurv(X_pool, surv_pool, penalty = penalty,
                      lambda = lambda, warn = FALSE)
  beta_hat <- as.numeric(coef(fit))
  list(alpha_hat  = beta_hat,
       eps_hat    = matrix(0, K, length(beta_hat)),
       theta_hat  = matrix(rep(beta_hat, K), nrow = K, byrow = TRUE),
       active_set = which(abs(beta_hat) > 1e-8),
       lambda     = lambda)
}


################################################################################
##  2.  BOOTSTRAP UNCERTAINTY QUANTIFICATION                                  ##
################################################################################

#' Stratified bootstrap for the hierarchical MCP meta-Cox model.
#'
#' Within each bootstrap resample, patients are resampled with replacement
#' WITHIN each study (stratified bootstrap), preserving the K-study structure
#' and the per-study event rates. The tuning parameters lambda_alpha and
#' lambda_eps are fixed to the values from the original fit so that each
#' bootstrap replicate uses identical regularization strength, ensuring that
#' the bootstrap variance reflects estimation uncertainty rather than
#' tuning variability. This is the standard approach for bootstrap inference
#' in penalized survival models (Efron & Tibshirani 1993; Tibshirani 1997).
#'
#' @param data_list      Original K-study data list (study_id, time, status, X..)
#' @param fit_original   Output of fit_hierarchical_meta_cox() on the full data
#' @param B              Number of bootstrap resamples (500 recommended for paper)
#' @param gene_names     Character vector of gene names (length p)
#' @param conf_level     Nominal coverage for CIs (default 0.95)
#' @param n_cores        Cores for parallel bootstrap (1 = sequential)
#' @param seed           RNG seed for reproducibility
#'
#' @return List:
#'   alpha_boot   : B x p matrix of bootstrap alpha estimates
#'   eps_boot     : B x K x p array of bootstrap epsilon estimates
#'   alpha_se     : p-vector of bootstrap standard errors for alpha
#'   alpha_ci_lo  : p-vector of 2.5th percentile bootstrap CIs
#'   alpha_ci_hi  : p-vector of 97.5th percentile bootstrap CIs
#'   eps_se       : K x p matrix of bootstrap SEs for epsilon
#'   eps_ci_lo    : K x p matrix of 2.5th percentile CIs for epsilon
#'   eps_ci_hi    : K x p matrix of 97.5th percentile CIs for epsilon
#'   B_used       : number of successful resamples (< B if any failed)
bootstrap_mcp_hier <- function(
    data_list,
    fit_original,
    B            = 500L,
    gene_names   = NULL,
    conf_level   = 0.95,
    n_cores      = 1L,
    seed         = 42L
) {
  K   <- length(data_list)
  p   <- length(fit_original$alpha_hat)
  lam_a <- fit_original$lambda_alpha_used
  lam_e <- fit_original$lambda_eps_used
  
  alpha_lo <- (1 - conf_level) / 2
  alpha_hi <- 1 - alpha_lo
  
  cat(sprintf(
    "\n  Bootstrap: B=%d | K=%d studies | p=%d genes | cores=%d\n",
    B, K, p, n_cores))
  cat(sprintf(
    "  Fixed lambda_alpha=%.5f  lambda_eps=%.5f\n", lam_a, lam_e))
  cat(sprintf("  Percentile CIs at %.1f%%\n\n", conf_level * 100))
  
  ## ── Single bootstrap worker ──────────────────────────────────────────────
  one_boot <- function(b) {
    set.seed(seed + b)
    # Stratified resample: sample within each study independently
    boot_data <- lapply(data_list, function(dk) {
      idx <- sample(nrow(dk), replace = TRUE)
      dk[idx, ]
    })
    tryCatch(
      fit_hierarchical_meta_cox(
        boot_data,
        penalty      = "mcp",
        lambda_alpha = lam_a,
        lambda_eps   = lam_e,
        max_iter     = 40L,
        conv_tol     = 1e-4,
        active_tol   = 1e-6,
        nfolds_cv    = 5L,
        verbose      = FALSE),
      error = function(e) NULL)
  }
  
  ## ── Run bootstrap (parallel or sequential) ──────────────────────────────
  if (n_cores > 1L) {
    cl <- makeCluster(n_cores)
    ## Export local variables (data_list, lam_a, lam_e, seed live in the
    ## function's local environment, so envir=environment() is correct here).
    clusterExport(cl, c("data_list", "lam_a", "lam_e", "seed"),
                  envir = environment())
    ## Export global helper functions separately using globalenv() — using
    ## environment() for these would fail because they are not local objects.
    clusterExport(cl, c("fit_hierarchical_meta_cox", "build_pooled_data",
                        "cv_select_lambda"),
                  envir = globalenv())
    clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(survival); library(ncvreg); library(glmnet)
      })
    })
    boot_fits <- parLapply(cl, seq_len(B), one_boot)
    stopCluster(cl)
  } else if (has_pbapply) {
    boot_fits <- pbapply::pblapply(seq_len(B), one_boot)
  } else {
    boot_fits <- vector("list", B)
    for (b in seq_len(B)) {
      if (b %% 50L == 0L) cat(sprintf("    bootstrap %d / %d\n", b, B))
      boot_fits[[b]] <- one_boot(b)
    }
  }
  
  ## ── Discard failed resamples ─────────────────────────────────────────────
  ok        <- !sapply(boot_fits, is.null)
  B_used    <- sum(ok)
  boot_fits <- boot_fits[ok]
  if (B_used < B)
    warning(sprintf("  %d / %d bootstrap resamples failed and were dropped.",
                    B - B_used, B))
  cat(sprintf("  %d successful resamples.\n", B_used))
  
  ## ── Collect alpha estimates (B_used x p) ─────────────────────────────────
  alpha_boot <- do.call(rbind, lapply(boot_fits, `[[`, "alpha_hat"))
  
  ## ── Collect epsilon estimates (B_used x K x p array) ────────────────────
  eps_boot <- array(
    unlist(lapply(boot_fits, `[[`, "eps_hat")),
    dim = c(K, p, B_used))           # dim: K x p x B_used
  eps_boot <- aperm(eps_boot, c(3, 1, 2))  # -> B_used x K x p
  
  ## ── Compute SE and percentile CIs for alpha ──────────────────────────────
  alpha_se    <- apply(alpha_boot, 2, sd,
                       na.rm = TRUE)
  alpha_ci_lo <- apply(alpha_boot, 2, quantile,
                       probs = alpha_lo, na.rm = TRUE)
  alpha_ci_hi <- apply(alpha_boot, 2, quantile,
                       probs = alpha_hi, na.rm = TRUE)
  
  ## ── Compute theta_boot: B_used x K x p ──────────────────────────────────
  ## theta_kj = alpha_j + epsilon_kj, so the correct bootstrap distribution
  ## of theta_kj propagates uncertainty from BOTH alpha_j and epsilon_kj.
  ## We compute this directly from the bootstrap replicates rather than
  ## approximating with eps_se alone (which would ignore the alpha uncertainty).
  theta_boot <- sweep(eps_boot, c(1, 3), alpha_boot, "+")
  # alpha_boot is B_used x p; eps_boot is B_used x K x p.
  # sweep over dims 1 (B) and 3 (p): adds alpha_boot[b,j] to eps_boot[b,k,j].
  
  theta_se    <- apply(theta_boot, c(2, 3), sd,      na.rm = TRUE)  # K x p
  theta_ci_lo <- apply(theta_boot, c(2, 3), quantile,
                       probs = alpha_lo, na.rm = TRUE)               # K x p
  theta_ci_hi <- apply(theta_boot, c(2, 3), quantile,
                       probs = alpha_hi, na.rm = TRUE)               # K x p
  
  ## ── Compute SE and percentile CIs for epsilon ────────────────────────────
  # eps_boot is B_used x K x p; apply over dim 1 (bootstrap replicates)
  eps_se    <- apply(eps_boot, c(2, 3), sd,      na.rm = TRUE)  # K x p
  eps_ci_lo <- apply(eps_boot, c(2, 3), quantile,
                     probs = alpha_lo, na.rm = TRUE)             # K x p
  eps_ci_hi <- apply(eps_boot, c(2, 3), quantile,
                     probs = alpha_hi, na.rm = TRUE)             # K x p
  
  ## ── Attach gene / study names ─────────────────────────────────────────────
  ## Guard: gene_names length must match p (= length(fit_original$alpha_hat)).
  ## A mismatch arises when common_genes was computed before one gene was
  ## filtered out during preprocessing, making length(gene_names) != p.
  if (!is.null(gene_names)) {
    if (length(gene_names) != p) {
      warning(sprintf(
        paste0("bootstrap_mcp_hier: length(gene_names)=%d does not match ",
               "p=%d (derived from fit_original$alpha_hat). ",
               "Gene name assignment skipped to avoid a dimnames error. ",
               "Pass the gene names vector whose length equals p."),
        length(gene_names), p))
    } else {
      colnames(alpha_boot) <- gene_names
      names(alpha_se)      <- gene_names
      names(alpha_ci_lo)   <- gene_names
      names(alpha_ci_hi)   <- gene_names
      colnames(eps_se)     <- gene_names
      colnames(eps_ci_lo)  <- gene_names
      colnames(eps_ci_hi)  <- gene_names
      colnames(theta_se)    <- gene_names
      colnames(theta_ci_lo) <- gene_names
      colnames(theta_ci_hi) <- gene_names
    }
  }
  
  list(
    alpha_boot  = alpha_boot,
    eps_boot    = eps_boot,
    theta_boot  = theta_boot,
    alpha_se    = alpha_se,
    alpha_ci_lo = alpha_ci_lo,
    alpha_ci_hi = alpha_ci_hi,
    eps_se      = eps_se,
    eps_ci_lo   = eps_ci_lo,
    eps_ci_hi   = eps_ci_hi,
    theta_se    = theta_se,
    theta_ci_lo = theta_ci_lo,
    theta_ci_hi = theta_ci_hi,
    B_used      = B_used,
    conf_level  = conf_level
  )
}


################################################################################
##  3.  DATA LOADING AND PREPROCESSING                                        ##
################################################################################

cat("\n================================================================\n")
cat("  7.1  Loading curatedOvarianData\n")
cat("================================================================\n\n")

candidate_eset_names <- c(
  "TCGA_eset",
  "GSE32062.GPL6480_eset",
  "GSE9891_eset",
  "GSE26712_eset",
  "GSE17260_eset",
  "GSE30161_eset"
)

extract_cohort <- function(eset_name, top_var_genes = 2000L) {
  tryCatch({
    data(list = eset_name, package = "curatedOvarianData", envir = environment())
    eset <- get(eset_name, envir = environment())
    pd   <- pData(eset)
    
    time_col   <- intersect(c("days_to_death","os_time","survival_time",
                              "t.os","T.OS"), colnames(pd))[1]
    status_col <- intersect(c("vital_status","os_status","censored",
                              "e.os","E.OS"), colnames(pd))[1]
    if (is.na(time_col) || is.na(status_col)) return(NULL)
    
    time_raw   <- as.numeric(pd[[time_col]])
    status_raw <- pd[[status_col]]
    status_bin <- if (is.character(status_raw) || is.factor(status_raw))
      as.integer(tolower(as.character(status_raw)) %in%
                   c("deceased","dead","1","died","yes"))
    else as.integer(as.numeric(status_raw))
    
    keep <- !is.na(time_raw) & !is.na(status_bin) & time_raw > 0
    if (sum(keep) < 30 || sum(status_bin[keep]) < 10) return(NULL)
    
    time_use <- time_raw[keep] / 365.25
    expr_mat <- t(exprs(eset)[, keep])
    iqr_vals  <- apply(expr_mat, 2, IQR)
    top_genes <- order(iqr_vals, decreasing = TRUE)[
      seq_len(min(top_var_genes, ncol(expr_mat)))]
    expr_use  <- scale(expr_mat[, top_genes, drop = FALSE])
    
    df <- data.frame(study_id = eset_name, time = time_use,
                     status   = status_bin[keep],
                     expr_use, stringsAsFactors = FALSE)
    list(data       = df,
         n          = nrow(df),
         n_events   = sum(status_bin[keep]),
         cens_rate  = 1 - mean(status_bin[keep]),
         gene_names = colnames(expr_use),
         eset_name  = eset_name)
  }, error = function(e) {
    message(sprintf("  [SKIP] %s : %s", eset_name, conditionMessage(e)))
    NULL
  })
}

cat("Loading cohorts (top_var_genes = 2000)...\n")
cohort_list <- Filter(Negate(is.null),
                      lapply(candidate_eset_names, extract_cohort,
                             top_var_genes = 2000L))
cat(sprintf("  Loaded %d / %d cohorts\n",
            length(cohort_list), length(candidate_eset_names)))

gene_sets    <- lapply(cohort_list, `[[`, "gene_names")
common_genes <- Reduce(intersect, gene_sets)
cat(sprintf("  Common genes after intersection: p = %d\n", length(common_genes)))

study_data_list <- lapply(cohort_list, function(ch) {
  df        <- ch$data
  keep_cols <- c("study_id","time","status",
                 intersect(common_genes, colnames(df)))
  df_sub    <- df[, keep_cols, drop = FALSE]
  ecols     <- keep_cols[-(1:3)]
  scaled    <- scale(df_sub[, ecols])
  scaled[is.nan(scaled)] <- 0
  df_sub[, ecols] <- scaled
  df_sub
})

K       <- length(study_data_list)
p       <- length(common_genes)
N_total <- sum(sapply(study_data_list, nrow))
n_events <- sum(sapply(study_data_list, function(d) sum(d$status)))
cat(sprintf("\n  Final: K=%d | N=%d | p=%d | events=%d (%.1f%%)\n",
            K, N_total, p, n_events, 100 * n_events / N_total))


################################################################################
##  TABLE 7.1  Cohort Characteristics                                         ##
################################################################################

table_cohort <- do.call(rbind, lapply(seq_len(K), function(k) {
  ch  <- cohort_list[[k]]
  d   <- study_data_list[[k]]
  sf  <- survfit(Surv(time, status) ~ 1, data = d)
  med <- unname(summary(sf)$table["median"])
  data.frame(
    Study         = LETTERS[k],
    Dataset       = gsub("_eset$", "", ch$eset_name),
    N             = ch$n,
    Events        = ch$n_events,
    Censoring_pct = sprintf("%.0f%%", ch$cens_rate * 100),
    Median_OS_yr  = sprintf("%.1f", med),
    Study_Label   = paste0("Study ", LETTERS[k]),
    stringsAsFactors = FALSE
  )
}))
write.csv(table_cohort[, -7], fig("Table7_1_cohort_chars.csv"),
          row.names = FALSE)
cat("\nTable 7.1 saved.\n")


################################################################################
##  FIGURE 7.1  Cohort Composition Bar Chart                                  ##
################################################################################

cohort_bar_df <- table_cohort %>%
  mutate(Censored = N - Events) %>%
  tidyr::pivot_longer(c(Events, Censored),
                      names_to = "Type", values_to = "Count") %>%
  mutate(Type = factor(Type, levels = c("Censored", "Events")))

p_cohort_bar <- ggplot(cohort_bar_df,
                       aes(x = Study_Label, y = Count, fill = Type)) +
  geom_col(width = 0.65, colour = "white") +
  geom_text(
    data = table_cohort %>% mutate(ylab = N + max(N) * 0.02),
    aes(x = Study_Label, y = ylab, label = paste0("N=", N)),
    inherit.aes = FALSE, size = 3.2, fontface = "bold", vjust = 0) +
  scale_fill_manual(
    values = c(Events = "#2980B9", Censored = "#BDC3C7"), name = NULL) +
  labs(title    = "Cohort Composition",
       subtitle = sprintf("%d independent ovarian cancer cohorts (curatedOvarianData)", K),
       x = NULL, y = "Number of patients") +
  theme_paper() + theme(legend.position = "top")
ggsave(fig("Fig7_1_cohort_bar.png"), p_cohort_bar,
       width = 8, height = 5, dpi = 300)
cat("Figure 7.1 saved.\n")


################################################################################
##  FIGURE 7.2  Kaplan-Meier Curves                                           ##
################################################################################

km_data <- do.call(rbind, lapply(seq_len(K), function(k) {
  d <- study_data_list[[k]]
  d$study_label <- paste0("Study ", LETTERS[k])
  d
}))
km_fit <- survfit(Surv(time, status) ~ study_label, data = km_data)
p_km <- ggsurvplot(
  km_fit, data = km_data,
  palette           = c("#E41A1C","#377EB8","#4DAF4A",
                        "#984EA3","#FF7F00","#A65628"),
  conf.int          = FALSE,
  risk.table        = TRUE,
  risk.table.height = 0.28,
  xlab              = "Time (years)",
  ylab              = "Overall survival probability",
  title             = "Kaplan-Meier Curves by Cohort",
  legend.title      = "",
  ggtheme           = theme_paper(),
  fontsize          = 3.5)
ggsave(fig("Fig7_2_km_curves.png"), print(p_km), width = 9, height = 7)
cat("Figure 7.2 saved.\n")


################################################################################
##  7.3  METHOD FITTING                                                       ##
################################################################################

cat("\n================================================================\n")
cat("  7.3  Fitting all six methods\n")
cat("================================================================\n\n")

set.seed(2024L)

cat("Fitting MCP (Hierarchical)...\n")
fit_mcp_hier <- fit_hierarchical_meta_cox(
  study_data_list, penalty = "mcp",
  lambda_alpha = NULL, lambda_eps = NULL,
  max_iter = 40L, verbose = TRUE)

cat("\nFitting SCAD (Hierarchical)...\n")
fit_scad_hier <- fit_hierarchical_meta_cox(
  study_data_list, penalty = "scad",
  lambda_alpha = NULL, lambda_eps = NULL,
  max_iter = 40L, verbose = TRUE)

cat("\nFitting MCP (Pooled)...\n")
fit_mcp_pool <- fit_ncvreg_pooled(study_data_list, penalty = "MCP")

cat("\nFitting SCAD (Pooled)...\n")
fit_scad_pool <- fit_ncvreg_pooled(study_data_list, penalty = "SCAD")

cat("\nFitting LASSO (Pooled)...\n")
fit_lasso <- fit_pooled_penalized(study_data_list, method = "lasso")

cat("\nFitting Elastic Net (Pooled)...\n")
fit_enet <- fit_pooled_penalized(study_data_list, method = "enet")

all_fits <- list(
  "MCP (Hier.)"    = fit_mcp_hier,
  "SCAD (Hier.)"   = fit_scad_hier,
  "MCP (Pooled)"   = fit_mcp_pool,
  "SCAD (Pooled)"  = fit_scad_pool,
  "LASSO (Pooled)" = fit_lasso,
  "ENet (Pooled)"  = fit_enet
)

cat("\n  Selected gene set sizes:\n")
for (nm in names(all_fits))
  cat(sprintf("    %-18s : %d genes\n", nm, length(all_fits[[nm]]$active_set)))


################################################################################
##  7.3b  BOOTSTRAP UNCERTAINTY QUANTIFICATION FOR MCP (Hierarchical)        ##
################################################################################

## ── Configuration ────────────────────────────────────────────────────────────
## B = 500 resamples gives stable SEs; increase to 1000 for publication.
## Set n_cores > 1 to parallelise (e.g. detectCores() - 1).
BOOTSTRAP_B      <- 500L
BOOTSTRAP_CORES  <- 1L          # change to detectCores()-1 for speed
BOOTSTRAP_SEED   <- 2024L
CONF_LEVEL       <- 0.95

cat("\n================================================================\n")
cat("  7.3b  Bootstrap uncertainty quantification\n")
cat(sprintf("        B=%d  cores=%d  seed=%d\n",
            BOOTSTRAP_B, BOOTSTRAP_CORES, BOOTSTRAP_SEED))
cat("================================================================\n")

## Derive p from the fit (not from common_genes) to avoid a length mismatch
## that occurs when common_genes has a different length than fit_mcp_hier$alpha_hat
## (e.g. if a zero-variance gene was silently dropped during scale()).
p_fit        <- length(fit_mcp_hier$alpha_hat)
gene_names_p <- if (length(common_genes) == p_fit) {
  common_genes
} else {
  warning(sprintf(
    paste0("length(common_genes)=%d != p=%d from fit. ",
           "Using first %d common_genes for bootstrap naming."),
    length(common_genes), p_fit, p_fit))
  common_genes[seq_len(p_fit)]
}

boot_results <- bootstrap_mcp_hier(
  data_list    = study_data_list,
  fit_original = fit_mcp_hier,
  B            = BOOTSTRAP_B,
  gene_names   = gene_names_p,
  conf_level   = CONF_LEVEL,
  n_cores      = BOOTSTRAP_CORES,
  seed         = BOOTSTRAP_SEED
)

## Convenience extractors
alpha_se    <- boot_results$alpha_se
alpha_ci_lo <- boot_results$alpha_ci_lo
alpha_ci_hi <- boot_results$alpha_ci_hi
eps_se      <- boot_results$eps_se       # K x p
eps_ci_lo   <- boot_results$eps_ci_lo   # K x p
eps_ci_hi   <- boot_results$eps_ci_hi   # K x p
theta_se    <- boot_results$theta_se    # K x p  (propagates alpha + eps uncertainty)
theta_ci_lo <- boot_results$theta_ci_lo # K x p
theta_ci_hi <- boot_results$theta_ci_hi # K x p

cat(sprintf("\n  Bootstrap complete: %d / %d resamples used.\n",
            boot_results$B_used, BOOTSTRAP_B))


################################################################################
##  7.4  RESULTS: BIOMARKER SELECTION AND REPRODUCIBILITY                    ##
################################################################################

cat("\n================================================================\n")
cat("  7.4  Biomarker selection results\n")
cat("================================================================\n\n")

mcp_active <- fit_mcp_hier$active_set

if (length(mcp_active) > 0) {
  
  ord       <- order(abs(fit_mcp_hier$alpha_hat[mcp_active]),
                     decreasing = TRUE)
  top_idx   <- mcp_active[ord]
  top_genes <- common_genes[top_idx]
  
  ## ── TABLE 7.2  Extended biomarker table with SE and 95% CI ───────────────
  ## Columns: Gene | MCP(Hier.) | SE | CI_lo | CI_hi | SCAD(Hier.) |
  ##          MCP(Pool.) | SCAD(Pool.) | LASSO | ENet | N_methods
  tbl2 <- data.frame(Gene = top_genes, stringsAsFactors = FALSE)
  
  ## Point estimates from all six methods
  for (nm in names(all_fits)) {
    vals <- round(all_fits[[nm]]$alpha_hat[top_idx], 4)
    vals[abs(vals) < 1e-8] <- 0
    tbl2[[nm]] <- vals
  }
  
  ## Bootstrap SE and 95% CI for MCP (Hier.) — primary method
  tbl2[["SE_MCP_Hier"]]    <- round(alpha_se[top_idx],    4)
  tbl2[["CI95_lo_MCP"]]    <- round(alpha_ci_lo[top_idx], 4)
  tbl2[["CI95_hi_MCP"]]    <- round(alpha_ci_hi[top_idx], 4)
  
  ## Significance flag: CI excludes zero
  tbl2[["Sig_95"]] <- ifelse(
    tbl2[["CI95_lo_MCP"]] > 0 | tbl2[["CI95_hi_MCP"]] < 0, "*", "")
  
  ## Reproducibility score
  tbl2$N_methods <- rowSums(
    tbl2[, names(all_fits), drop = FALSE] != 0)
  
  ## Reorder: reproducibility desc, then |alpha| desc
  tbl2 <- tbl2[order(-tbl2$N_methods,
                     -abs(tbl2[["MCP (Hier.)"]])), ]
  
  write.csv(tbl2, fig("Table7_2_biomarkers_CI.csv"), row.names = FALSE)
  cat("Table 7.2 (with bootstrap CI) saved.\n\n")
  
  ## Pretty print to console
  cat(sprintf("  %-10s  %8s  %6s  %8s  %8s  %5s  %s\n",
              "Gene", "alpha_hat", "SE", "CI_lo", "CI_hi", "Nmeth", "Sig"))
  cat(strrep("-", 62), "\n")
  for (i in seq_len(nrow(tbl2))) {
    r <- tbl2[i, ]
    cat(sprintf("  %-10s  %8.4f  %6.4f  %8.4f  %8.4f  %5d  %s\n",
                r$Gene,
                r[["MCP (Hier.)"]],
                r$SE_MCP_Hier,
                r$CI95_lo_MCP,
                r$CI95_hi_MCP,
                r$N_methods,
                r$Sig_95))
  }
  cat("\n")
  
} else {
  cat("  No genes selected by MCP (Hier.).\n")
  tbl2 <- NULL
}


################################################################################
##  FIGURE 7.3  Forest Plot WITH Bootstrap CI Error Bars                     ##
################################################################################

if (!is.null(tbl2) && length(mcp_active) > 0) {
  
  alpha_active <- fit_mcp_hier$alpha_hat[mcp_active]
  n_show       <- min(12, length(mcp_active))
  show_idx     <- mcp_active[order(abs(alpha_active),
                                   decreasing = TRUE)][seq_len(n_show)]
  show_genes_f <- common_genes[show_idx]
  
  ## Build data frame: study-specific rows + global row per gene
  forest_df <- do.call(rbind, lapply(seq_along(show_idx), function(gi) {
    j    <- show_idx[gi]
    gene <- common_genes[j]
    rbind(
      do.call(rbind, lapply(seq_len(K), function(k)
        data.frame(
          gene     = gene,
          label    = paste0("Study ", LETTERS[k]),
          estimate = fit_mcp_hier$theta_hat[k, j],
          ci_lo    = theta_ci_lo[k, j],   # bootstrap percentile CI for theta_kj
          ci_hi    = theta_ci_hi[k, j],   # propagates alpha + epsilon uncertainty
          type     = "Study-specific",
          stringsAsFactors = FALSE))),
      data.frame(
        gene     = gene,
        label    = "Global (meta)",
        estimate = fit_mcp_hier$alpha_hat[j],
        ci_lo    = alpha_ci_lo[j],             # bootstrap percentile CI
        ci_hi    = alpha_ci_hi[j],
        type     = "Global",
        stringsAsFactors = FALSE))
  }))
  
  forest_df$gene  <- factor(forest_df$gene,  levels = rev(show_genes_f))
  forest_df$label <- factor(forest_df$label,
                            levels = c(paste0("Study ", LETTERS[seq_len(K)]),
                                       "Global (meta)"))
  
  p_forest <- ggplot(forest_df,
                     aes(x = estimate, y = label,
                         colour = type, size = type, shape = type)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                   height = 0.3, linewidth = 0.5, alpha = 0.7) +
    geom_point(alpha = 0.9) +
    scale_colour_manual(
      values = c("Global" = "#C0392B", "Study-specific" = "#2980B9"),
      name = NULL) +
    scale_size_manual(
      values = c("Global" = 4, "Study-specific" = 2.5), name = NULL) +
    scale_shape_manual(
      values = c("Global" = 18, "Study-specific" = 16), name = NULL) +
    facet_wrap(~gene, scales = "free_x", ncol = 3) +
    labs(
      title    = "Global and Study-Specific Log-Hazard Ratios",
      subtitle = paste0("Red diamond = global alpha_j (bootstrap 95% CI); ",
                        "blue circles = study-specific theta_kj (bootstrap 95% CI)"),
      x = "Log-hazard ratio", y = NULL) +
    theme_paper(base_size = 10) +
    theme(legend.position = "top")
  
  ggsave(fig("Fig7_3_forest_CI.png"), p_forest,
         width = 11,
         height = max(6, ceiling(n_show / 3) * 3),
         dpi = 300)
  cat("Figure 7.3 (with CIs) saved.\n")
}


################################################################################
##  7.5  CROSS-STUDY HETEROGENEITY                                           ##
################################################################################

cat("\n================================================================\n")
cat("  7.5  Cross-study heterogeneity\n")
cat("================================================================\n\n")

if (length(mcp_active) > 0) {
  
  alpha_active <- fit_mcp_hier$alpha_hat[mcp_active]
  active_genes <- common_genes[mcp_active]
  eps_active   <- fit_mcp_hier$eps_hat[, mcp_active, drop = FALSE]
  tau2         <- apply(eps_active, 2, var)
  I2           <- tau2 / (tau2 + 1e-6)
  
  het_summary <- data.frame(
    gene      = active_genes,
    alpha_hat = alpha_active,
    tau2      = tau2,
    I2        = I2,
    het_level = cut(I2, c(-Inf, 0.25, 0.50, 0.75, Inf),
                    labels = c("Low","Moderate","Substantial","High"))
  ) %>% arrange(desc(abs(alpha_hat)))
  cat("Heterogeneity summary (MCP Hierarchical):\n")
  print(het_summary)
  
  ## ── TABLE 7.4  Study-specific deviations with bootstrap 95% CIs ──────────
  ## Long-format table: one row per (gene, study) pair.
  ## Columns: Gene | Study | epsilon_hat | SE | CI_lo | CI_hi | Sig_95
  tbl4_rows <- do.call(rbind, lapply(seq_along(mcp_active), function(gi) {
    j    <- mcp_active[gi]
    gene <- common_genes[j]
    do.call(rbind, lapply(seq_len(K), function(k) {
      est  <- fit_mcp_hier$eps_hat[k, j]
      se_v <- eps_se[k, j]
      lo   <- eps_ci_lo[k, j]
      hi   <- eps_ci_hi[k, j]
      data.frame(
        Gene      = gene,
        Study     = paste0("Study ", LETTERS[k]),
        epsilon_hat = round(est,  4),
        SE          = round(se_v, 4),
        CI95_lo     = round(lo,   4),
        CI95_hi     = round(hi,   4),
        Sig_95      = ifelse(lo > 0 | hi < 0, "*", ""),
        stringsAsFactors = FALSE)
    }))
  }))
  ## Order by gene (reproducibility) then study
  tbl4_rows <- tbl4_rows[order(tbl4_rows$Gene, tbl4_rows$Study), ]
  write.csv(tbl4_rows, fig("Table7_4_deviations_CI.csv"), row.names = FALSE)
  cat("\nTable 7.4 (deviation estimates with 95% CI) saved.\n")
  cat("  Preview (first 18 rows):\n")
  print(head(tbl4_rows, 18))
  
  ## ── FIGURE 7.4  Deviation Heatmap ────────────────────────────────────────
  n_tile     <- min(15, ncol(eps_active))
  tile_ord   <- order(abs(alpha_active), decreasing = TRUE)[seq_len(n_tile)]
  tile_genes <- active_genes[tile_ord]
  eps_tile   <- eps_active[, tile_ord, drop = FALSE]
  
  eps_df           <- as.data.frame(eps_tile)
  colnames(eps_df) <- tile_genes
  eps_df$Study     <- paste0("Study ", LETTERS[seq_len(K)])
  eps_long         <- tidyr::pivot_longer(eps_df, -Study,
                                          names_to  = "Gene",
                                          values_to = "epsilon")
  eps_long$Gene  <- factor(eps_long$Gene, levels = rev(tile_genes))
  max_eps        <- max(abs(eps_long$epsilon), na.rm = TRUE)
  
  p_tile <- ggplot(eps_long, aes(x = Study, y = Gene, fill = epsilon)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", epsilon)),
              size = 2.6, colour = "grey20") +
    scale_fill_gradient2(
      low = "#2980B9", mid = "white", high = "#C0392B",
      midpoint = 0, limits = c(-max_eps, max_eps),
      name = expression(hat(epsilon)[kj])) +
    labs(
      title    = "Study-Specific Deviations from Global Effects",
      subtitle = expression(hat(epsilon)[kj] == hat(theta)[kj] - hat(alpha)[j] ~
                              "  |  MCP (Hierarchical)"),
      x = NULL, y = NULL) +
    theme_paper(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(fig("Fig7_4_deviation_tile.png"), p_tile,
         width  = max(6, K * 1.1 + 2),
         height = max(5, n_tile * 0.45 + 2),
         dpi    = 300)
  cat("Figure 7.4 saved.\n")
  
  ## ── FIGURE 7.5  MCP vs SCAD Study-Specific Estimates ────────────────────
  n_genes   <- ncol(fit_mcp_hier$theta_hat)
  theta_cmp <- do.call(rbind, lapply(seq_len(K), function(k)
    data.frame(
      study      = paste0("Study ", LETTERS[k]),
      theta_MCP  = fit_mcp_hier$theta_hat[k, ],
      theta_SCAD = fit_scad_hier$theta_hat[k, ],
      active     = ifelse(seq_len(n_genes) %in% fit_mcp_hier$active_set,
                          "Selected by MCP", "Not selected"),
      stringsAsFactors = FALSE)))
  
  p_cmp <- ggplot(theta_cmp,
                  aes(x = theta_MCP, y = theta_SCAD,
                      colour = active, alpha = active)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey50") +
    geom_point(size = 1.2) +
    scale_colour_manual(
      values = c("Selected by MCP" = "#C0392B", "Not selected" = "grey70"),
      name = NULL) +
    scale_alpha_manual(
      values = c("Selected by MCP" = 0.85, "Not selected" = 0.25),
      name = NULL) +
    facet_wrap(~study, ncol = min(K, 3)) +
    labs(
      title    = "MCP vs SCAD Hierarchical Study-Specific Estimates",
      subtitle = expression(hat(theta)[kj] ~
                              "comparison across all studies and p genes"),
      x = expression("MCP (Hier.)  " * hat(theta)[kj]),
      y = expression("SCAD (Hier.) " * hat(theta)[kj])) +
    theme_paper(base_size = 10) + theme(legend.position = "top")
  ggsave(fig("Fig7_5_mcp_vs_scad.png"), p_cmp,
         width  = max(7, ceiling(K / 3) * 3.5 + 1),
         height = 6, dpi = 300)
  cat("Figure 7.5 saved.\n")
}


################################################################################
##  FIGURE 7.6  Selected Set Sizes                                            ##
################################################################################

sel_sizes <- data.frame(
  Method   = names(all_fits),
  Selected = sapply(all_fits, function(f) length(f$active_set)),
  Type     = c("Proposed","Proposed","Competitor",
               "Competitor","Competitor","Competitor"),
  stringsAsFactors = FALSE
)
sel_sizes$Method <- factor(sel_sizes$Method, levels = names(all_fits))

p_sel <- ggplot(sel_sizes, aes(x = Method, y = Selected, fill = Type)) +
  geom_col(width = 0.65, colour = "white") +
  geom_text(aes(label = Selected), vjust = -0.35, size = 3.5,
            fontface = "bold") +
  scale_fill_manual(
    values = c(Proposed = "#C0392B", Competitor = "grey60"), name = NULL) +
  labs(
    title    = "Selected Gene Set Sizes by Method",
    subtitle = sprintf("p = %d genes | K = %d cohorts | N = %d patients",
                       p, K, N_total),
    x = NULL, y = "Genes selected") +
  theme_paper() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "top")
ggsave(fig("Fig7_6_selection_sizes.png"), p_sel,
       width = 7, height = 5, dpi = 300)
cat("Figure 7.6 saved.\n")


################################################################################
##  TABLE 7.3  Method Summary                                                 ##
################################################################################

sel_bin <- sapply(all_fits, function(fit)
  as.integer(common_genes %in% common_genes[fit$active_set]))
rownames(sel_bin) <- common_genes
repro_genes <- rownames(sel_bin)[rowSums(sel_bin) >= 4]

tbl3 <- data.frame(
  Method = names(all_fits),
  p_total           = p,
  Selected_genes    = sapply(all_fits, function(f) length(f$active_set)),
  Reproducible_core = sapply(all_fits, function(f)
    length(intersect(common_genes[f$active_set], repro_genes))),
  Lambda_alpha = sapply(all_fits, function(f)
    round(if (!is.null(f$lambda_alpha_used) && !is.na(f$lambda_alpha_used))
      f$lambda_alpha_used
      else if (!is.null(f$lambda)) f$lambda else NA_real_, 5)),
  stringsAsFactors = FALSE
)
write.csv(tbl3, fig("Table7_3_method_summary.csv"), row.names = FALSE)
cat("\nTable 7.3 saved.\n")
print(tbl3)


################################################################################
##  FINAL SUMMARY                                                             ##
################################################################################

cat("\n================================================================\n")
cat(sprintf("  Section 7 complete.  Outputs in: %s/\n", out_dir))
cat(sprintf("  Dataset: curatedOvarianData | K=%d | N=%d | p=%d\n",
            K, N_total, p))
cat(sprintf("  Bootstrap: B=%d resamples (%d used) | %.0f%% CIs\n",
            BOOTSTRAP_B, boot_results$B_used, CONF_LEVEL * 100))
cat("\n  Output files:\n")
cat("    Figures : Fig7_1_cohort_bar.png\n")
cat("              Fig7_2_km_curves.png\n")
cat("              Fig7_3_forest_CI.png        <- NEW: includes CI error bars\n")
cat("              Fig7_4_deviation_tile.png\n")
cat("              Fig7_5_mcp_vs_scad.png\n")
cat("              Fig7_6_selection_sizes.png\n")
cat("    Tables  : Table7_1_cohort_chars.csv\n")
cat("              Table7_2_biomarkers_CI.csv  <- NEW: SE and 95% CI columns\n")
cat("              Table7_3_method_summary.csv\n")
cat("              Table7_4_deviations_CI.csv  <- NEW: epsilon CIs per study\n")
cat("================================================================\n")