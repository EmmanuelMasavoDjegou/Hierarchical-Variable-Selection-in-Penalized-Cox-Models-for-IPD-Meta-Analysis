################################################################################
##                                                                            ##
##   SECTION 7: REAL DATA APPLICATION                                         ##
##   "Additive Hierarchical Variable Selection for Nonconvex Penalized        ##
##    Cox Models in Individual-Participant-Data Meta-Analysis"                ##
##                                                                            ##
##   Datasets:                                                                ##
##     1. curatedOvarianData  – multi-cohort ovarian cancer microarray        ##
##        survival data (Ganzfried et al., 2013)                              ##                ##
##                                                                            ##
##   Analysis pipeline:                                                       ##
##     7.1  Data description & cohort characteristics                         ##
##     7.2  Preprocessing (QC, imputation, scaling)                           ##
##     7.3  Analysis setup (method fitting, CV tuning)                        ##
##     7.4  Results (selected biomarkers, coefficient profiles)               ##
##     7.5  Cross-study heterogeneity (deviation profiles, forest plots)      ##
##                                                                            ##
##   OUTPUT FILES:                                                            ##
##     Figures  : Fig7_1_cohort_summary.pdf  Fig7_2_survival_curves.pdf       ##
##                Fig7_3_selection_heatmap.pdf  Fig7_4_coef_dotplot.pdf       ##
##                Fig7_5_heterogeneity_forest.pdf  Fig7_6_deviation_tile.pdf  ##
##                Fig7_7_theta_comparison.pdf  Fig7_8_method_comparison.pdf   ##
##     Tables   : Table7_1_cohort_chars.csv  Table7_2_selected_genes.csv      ##
##                Table7_3_method_comparison.csv                              ##
##                                                                            ##
################################################################################

## ── 0. Packages ──────────────────────────────────────────────────────────────
# Install Bioconductor packages the first time:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("curatedOvarianData")
#install.packages(c("survival","ncvreg","glmnet","MASS","ggplot2","tidyr",
#                    "dplyr","patchwork","ggrepel","RColorBrewer","scales",
#                    "survminer","stringr","kableExtra","gridExtra"))

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
})

## ── source simulation helpers ─────────────────────────────────────────────────
## The functions below are copied from meta-Cox_model.R so this script is
## self-contained.  If you already have meta-Cox_model.R in your working
## directory you can replace this block with:  source("meta-Cox_model.R")

# ---------- paste the full body of meta-Cox_model.R here, OR just source it:
# source("meta-Cox_model.R")   # <-- adjust path if needed

################################################################################
##  SECTION 1: DATA GENERATING PROCESS                                        ##
################################################################################

#' Generate IPD survival data for K studies under the additive hierarchical model
#'
#' Model:  h_ki(t|x_ki) = h_k0(t) * exp(x_ki' * theta_k)
#'         theta_kj = alpha_j + epsilon_kj          [Eq. 1 in paper]
#'
#' Baseline hazard:  Weibull(shape, lambda_k),  lambda_k ~ Unif(range)
#' Event times:      Inverse-CDF method (Bender et al. 2005)
#'                   T = (-log(U) / (lambda_k * exp(eta)))^(1/shape)
#' Covariates:       AR(1)  Sigma_jl = rho^|j-l|   (Fan & Li 2001 benchmark)
#' Heterogeneity:    epsilon_kj ~ N(0, epsilon_sd^2) for j = 1..s_epsilon
#'
#' @param K                   Number of studies
#' @param n_per_study         Sample sizes (scalar = equal; vector = unequal)
#' @param p                   Total number of covariates
#' @param s_alpha             Number of truly active global covariates
#' @param alpha_true          True global effects (length s_alpha)
#' @param epsilon_sd          SD of study-specific deviations (heterogeneity)
#' @param s_epsilon           Number of active covariates with study heterogeneity
#' @param rho                 AR(1) correlation parameter
#' @param cens_rate           Target censoring proportion
#' @param weibull_shape       Shape parameter for Weibull baseline hazard
#' @param weibull_scale_range Uniform range for study-specific Weibull scales
#'
#' @return List: studies (K data frames), alpha_true, epsilon_true,
#'         theta_true, active_set, and DGP metadata
generate_ipd_data <- function(
    K                   = 5,
    n_per_study         = 200,
    p                   = 50,
    s_alpha             = 6,
    alpha_true          = c(1.0, -0.8, 0.6, -0.5, 0.4, -0.3),
    epsilon_sd          = 0.3,
    s_epsilon           = 4,
    rho                 = 0.5,
    cens_rate           = 0.30,
    weibull_shape       = 1.5,
    weibull_scale_range = c(0.05, 0.15)
) {
  if (length(n_per_study) == 1L) n_per_study <- rep(n_per_study, K)
  stopifnot(length(n_per_study) == K, length(alpha_true) == s_alpha,
            s_alpha <= p, s_epsilon <= s_alpha)
  
  # True sparse global effect vector (p-vector)
  alpha_full <- c(alpha_true, rep(0, p - s_alpha))
  
  # Study-specific Weibull scales: creates diversity in baseline hazards
  lambda_k <- runif(K, weibull_scale_range[1], weibull_scale_range[2])
  
  # AR(1) covariance: Sigma_jl = rho^|j-l|
  idx   <- seq_len(p)
  Sigma <- rho ^ abs(outer(idx, idx, "-"))
  
  # Study-specific deviations: N(0, epsilon_sd^2) for j <= s_epsilon
  # Mirrors classical random-effects meta-analysis (DerSimonian & Laird 1986)
  epsilon_list <- lapply(seq_len(K), function(k) {
    eps <- rep(0, p)
    if (s_epsilon > 0L)
      eps[seq_len(s_epsilon)] <- rnorm(s_epsilon, 0, epsilon_sd)
    eps
  })
  
  # Generate K study datasets
  study_data <- vector("list", K)
  for (k in seq_len(K)) {
    nk <- n_per_study[k]
    
    # Covariates: multivariate normal, column-standardised (Algorithm 1, Step 1)
    Xk <- mvrnorm(nk, mu = rep(0, p), Sigma = Sigma)
    Xk <- scale(Xk)
    colnames(Xk) <- paste0("X", seq_len(p))
    
    # Study-specific coefficient: theta_k = alpha + epsilon_k  (Eq. 1)
    theta_k <- alpha_full + epsilon_list[[k]]
    
    # Linear predictor
    eta <- as.numeric(Xk %*% theta_k)
    
    # Event times via Weibull inverse-CDF (Bender et al. 2005):
    #   T = (-log(U) / (lambda_k * exp(eta)))^(1/shape),  U ~ Unif(0,1)
    U       <- runif(nk)
    T_event <- (-log(U) / (lambda_k[k] * exp(eta))) ^ (1 / weibull_shape)
    
    # Exponential censoring calibrated to approximate target cens_rate
    lambda_c <- cens_rate / mean(T_event)
    C_time   <- rexp(nk, rate = lambda_c)
    
    Y_obs  <- pmin(T_event, C_time)
    delta  <- as.integer(T_event <= C_time)
    
    study_data[[k]] <- data.frame(
      study_id = k, time = Y_obs, status = delta, Xk
    )
  }
  
  list(
    studies      = study_data,
    alpha_true   = alpha_full,
    epsilon_true = epsilon_list,
    theta_true   = lapply(epsilon_list, function(eps) alpha_full + eps),
    active_set   = which(alpha_full != 0),
    K = K, p = p, s_alpha = s_alpha, s_epsilon = s_epsilon,
    n_per_study  = n_per_study
  )
}


################################################################################
##  SECTION 2: CORE ESTIMATION – HIERARCHICAL META-COX FITTER                 ##
################################################################################
##
##  Architecture (Algorithm 1 in the paper):
##
##  Stage 1 – Global selection:
##    Stack all K studies with per-subject offsets o_ki = x_ki' * eps_hat_k.
##    Fit penalised Cox:
##      SCAD / MCP  → ncvsurv() (Breheny & Huang 2011)
##      LASSO / ENet → glmnet()
##    Result: alpha_hat and active set A = {j : |alpha_j| > 0}
##
##  Stage 2 – Study-specific adjustment:
##    For each study k, fit penalised Cox with offset = X_k %*% alpha_hat.
##    Restrict to active columns A only (hierarchical constraint).
##    epsilon_kj = 0 for j not in A.
##
##  Stage 3 – Iterate until convergence.
##
################################################################################

#' Build stacked design matrix + survival outcome + offsets from K studies.
build_pooled_data <- function(data_list, offset_list = NULL) {
  K <- length(data_list)
  X_list <- y_list <- off_list <- sid_list <- vector("list", K)
  for (k in seq_len(K)) {
    dk            <- data_list[[k]]
    Xk            <- as.matrix(dk[, -(1:3), drop = FALSE])
    X_list[[k]]   <- Xk
    y_list[[k]]   <- cbind(time = dk$time, status = dk$status)
    off_list[[k]] <- if (is.null(offset_list)) rep(0, nrow(dk)) else offset_list[[k]]
    sid_list[[k]] <- rep(k, nrow(dk))
  }
  list(
    X        = do.call(rbind, X_list),
    y        = Surv(do.call(c, lapply(y_list, `[`, , "time")),
                    do.call(c, lapply(y_list, `[`, , "status"))),
    offset   = do.call(c, off_list),
    study_id = do.call(c, sid_list)
  )
}


#' Select lambda via CV for ncvsurv or glmnet (SCAD/MCP/LASSO/Ridge/ENet).
cv_select_lambda_ncv_glm <- function(X, y, offset, pen_upper,
                                     alpha_enet = 1, nfolds = 5L) {
  if (pen_upper %in% c("SCAD", "MCP")) {
    cv_fit <- tryCatch(
      cv.ncvsurv(X, y, penalty = pen_upper, offset = offset,
                 nfolds = nfolds, warn = FALSE),
      error = function(e) NULL)
    if (!is.null(cv_fit)) return(cv_fit$lambda.min)
    fit_path <- ncvsurv(X, y, penalty = pen_upper, offset = offset)
    nnz      <- colSums(fit_path$beta != 0)
    return(fit_path$lambda[which.min(abs(nnz - 6L))])
  } else {
    cv_fit <- tryCatch(
      cv.glmnet(X, y, family = "cox", alpha = alpha_enet, offset = offset,
                nfolds = nfolds, type.measure = "C"),
      error = function(e) NULL)
    if (!is.null(cv_fit)) return(cv_fit$lambda.min)
    fit_path <- glmnet(X, y, family = "cox", alpha = alpha_enet, offset = offset)
    return(fit_path$lambda[max(1L, length(fit_path$lambda) %/% 2L)])
  }
}


#' Fit the hierarchical penalised IPD meta-Cox model (Algorithm 1).
#'
#' @param data_list    List of K study data frames (study_id, time, status, X1..Xp)
#' @param penalty      "scad", "mcp", "lasso", or "enet"
#' @param lambda_alpha Tuning for global alpha  (NULL = CV)
#' @param lambda_eps   Tuning for study epsilons (NULL = CV)
#' @param max_iter     Maximum outer alternating iterations
#' @param conv_tol     Convergence: max|theta_new - theta_old| < conv_tol
#' @param active_tol   Hard-zero threshold after fitting
#' @param nfolds_cv    CV folds
#' @param verbose      Print iteration progress
#' @return List: alpha_hat, eps_hat (K x p), theta_hat (K x p),
#'         active_set, lambda_alpha_used, lambda_eps_used, n_iter
fit_hierarchical_meta_cox <- function(
    data_list,
    penalty      = "scad",
    lambda_alpha = NULL,
    lambda_eps   = NULL,
    max_iter     = 30L,
    conv_tol     = 1e-4,
    active_tol   = 1e-6,
    nfolds_cv    = 5L,
    verbose      = TRUE
) {
  K         <- length(data_list)
  p         <- ncol(data_list[[1]]) - 3L
  pen_upper <- toupper(penalty)
  
  use_ncvreg <- pen_upper %in% c("SCAD", "MCP")
  alpha_enet <- switch(pen_upper, "LASSO" = 1.0, "ENET" = 0.5, 1.0)
  
  # Initialise (Algorithm 1, Step 2)
  alpha_hat         <- rep(0, p)
  eps_hat           <- matrix(0, K, p)
  theta_hat         <- matrix(0, K, p)
  lambda_alpha_used <- NA_real_
  lambda_eps_used   <- NA_real_
  
  if (verbose)
    cat(sprintf("  Fitting: %-6s | lambda_alpha=%s | lambda_eps=%s\n",
                pen_upper,
                ifelse(is.null(lambda_alpha), "BIC/CV",
                       sprintf("%.4f", lambda_alpha)),
                ifelse(is.null(lambda_eps), "BIC/CV",
                       sprintf("%.4f", lambda_eps))))
  
  for (iter in seq_len(max_iter)) {
    theta_old <- theta_hat
    
    ## ── STAGE 1: Update alpha (global effects) ────────────────────────────
    ## Offset = X_k %*% eps_hat_k (current epsilon contribution per subject)
    ## Fit penalised Cox on stacked data from all K studies.
    
    offset_eps <- lapply(seq_len(K), function(k) {
      Xk <- as.matrix(data_list[[k]][, -(1:3), drop = FALSE])
      as.numeric(Xk %*% eps_hat[k, ])
    })
    
    pd <- build_pooled_data(data_list, offset_eps)
    
    if (use_ncvreg) {
      if (is.null(lambda_alpha))
        lambda_alpha_used <- cv_select_lambda_ncv_glm(
          pd$X, pd$y, pd$offset, pen_upper, nfolds = nfolds_cv)
      else
        lambda_alpha_used <- lambda_alpha
      
      fit_a     <- ncvsurv(pd$X, pd$y, penalty = pen_upper,
                           offset = pd$offset, lambda = lambda_alpha_used,
                           warn = FALSE)
      alpha_hat <- as.numeric(coef(fit_a))
      
    } else {
      if (is.null(lambda_alpha))
        lambda_alpha_used <- cv_select_lambda_ncv_glm(
          pd$X, pd$y, pd$offset, pen_upper, alpha_enet, nfolds = nfolds_cv)
      else
        lambda_alpha_used <- lambda_alpha
      
      fit_a     <- glmnet(pd$X, pd$y, family = "cox", alpha = alpha_enet,
                          offset = pd$offset, lambda = lambda_alpha_used)
      alpha_hat <- as.numeric(coef(fit_a))
    }
    
    # Hard-zero numerical dust
    alpha_hat[abs(alpha_hat) < active_tol] <- 0
    # Active set A = {j : |alpha_j| > 0}  (Algorithm 1, Step 3)
    active_set <- which(alpha_hat != 0)
    
    ## ── STAGE 2: Update epsilon_k (study-specific deviations) ────────────
    ## Hierarchical constraint: epsilon_kj = 0 for j NOT in active_set.
    ## Offset = X_k %*% alpha_hat (global effects absorbed as offset).
    
    eps_hat <- matrix(0, K, p)
    
    if (length(active_set) > 0L) {
      for (k in seq_len(K)) {
        dk          <- data_list[[k]]
        Xk          <- as.matrix(dk[, -(1:3), drop = FALSE])
        nk          <- nrow(dk)
        surv_k      <- Surv(dk$time, dk$status)
        off_alpha_k <- as.numeric(Xk %*% alpha_hat)
        Xk_active   <- Xk[, active_set, drop = FALSE]
        
        if (ncol(Xk_active) == 1L) {
          # Single active variable: plain coxph most stable
          df_tmp <- data.frame(time = dk$time, status = dk$status,
                               v1 = Xk_active[, 1])
          fm <- tryCatch(
            coxph(Surv(time, status) ~ v1 + offset(off_alpha_k),
                  data = df_tmp, control = coxph.control(iter.max = 50L)),
            error = function(e) NULL)
          if (!is.null(fm)) {
            est <- coef(fm)["v1"]
            eps_hat[k, active_set] <- ifelse(abs(est) < active_tol, 0, est)
          }
          
        } else if (use_ncvreg) {
          if (is.null(lambda_eps)) {
            lam_e <- tryCatch(
              cv_select_lambda_ncv_glm(
                Xk_active, surv_k, off_alpha_k, pen_upper,
                nfolds = min(nfolds_cv, max(3L, nk %/% 20L))),
              error = function(e) lambda_alpha_used * 0.5)
          } else {
            lam_e <- lambda_eps
          }
          lambda_eps_used <- lam_e
          fit_e <- tryCatch(
            ncvsurv(Xk_active, surv_k, penalty = pen_upper,
                    offset = off_alpha_k, lambda = lam_e, warn = FALSE),
            error = function(e) NULL)
          if (!is.null(fit_e)) {
            coefs <- as.numeric(coef(fit_e))
            coefs[abs(coefs) < active_tol] <- 0
            eps_hat[k, active_set] <- coefs
          }
          
        } else {
          if (is.null(lambda_eps)) {
            lam_e <- tryCatch(
              cv_select_lambda_ncv_glm(
                Xk_active, surv_k, off_alpha_k, pen_upper, alpha_enet,
                nfolds = min(nfolds_cv, max(3L, nk %/% 20L))),
              error = function(e) lambda_alpha_used * 0.5)
          } else {
            lam_e <- lambda_eps
          }
          lambda_eps_used <- lam_e
          fit_e <- tryCatch(
            glmnet(Xk_active, surv_k, family = "cox", alpha = alpha_enet,
                   offset = off_alpha_k, lambda = lam_e),
            error = function(e) NULL)
          if (!is.null(fit_e)) {
            coefs <- as.numeric(coef(fit_e))
            coefs[abs(coefs) < active_tol] <- 0
            eps_hat[k, active_set] <- coefs
          }
        }
      }
    }
    
    # Update theta_kj = alpha_j + epsilon_kj  (Eq. 1)
    theta_hat <- sweep(eps_hat, 2L, alpha_hat, "+")
    
    # Convergence check
    max_change <- max(abs(theta_hat - theta_old))
    if (verbose && (iter == 1L || iter %% 5L == 0L))
      cat(sprintf("    iter %2d | active=%2d | max_change=%.2e\n",
                  iter, length(active_set), max_change))
    
    if (max_change < conv_tol) {
      if (verbose)
        cat(sprintf("    Converged at iter %d (max_change=%.2e)\n",
                    iter, max_change))
      break
    }
    if (iter == max_iter && verbose)
      cat(sprintf("    Reached max_iter=%d (max_change=%.2e)\n",
                  max_iter, max_change))
  }
  
  list(
    alpha_hat         = alpha_hat,
    eps_hat           = eps_hat,
    theta_hat         = theta_hat,
    active_set        = which(alpha_hat != 0),
    lambda_alpha_used = lambda_alpha_used,
    lambda_eps_used   = lambda_eps_used,
    n_iter            = iter
  )
}


#' Fit POOLED penalised Cox (competitor) — naive stacking, no stratification.
fit_pooled_penalized <- function(data_list, method = "lasso",
                                 lambda = NULL, nfolds = 5L) {
  combined  <- do.call(rbind, data_list)
  X_pool    <- as.matrix(combined[, -(1:3), drop = FALSE])
  surv_pool <- Surv(combined$time, combined$status)
  K         <- length(data_list)
  alpha_g   <- switch(method,
                      "lasso" = 1.0, "enet" = 0.5,
                      stop("method must be 'lasso' or 'enet'"))
  
  if (is.null(lambda)) {
    cv_fit <- cv.glmnet(X_pool, surv_pool, family = "cox",
                        alpha = alpha_g, nfolds = nfolds, type.measure = "C")
    lambda <- cv_fit$lambda.min
  }
  fit      <- glmnet(X_pool, surv_pool, family = "cox",
                     alpha = alpha_g, lambda = lambda)
  beta_hat <- as.numeric(coef(fit))
  list(
    alpha_hat  = beta_hat,
    eps_hat    = matrix(0, K, length(beta_hat)),
    theta_hat  = matrix(rep(beta_hat, K), nrow = K, byrow = TRUE),
    active_set = which(abs(beta_hat) > 1e-8),
    lambda     = lambda
  )
}


#' Fit POOLED nonconvex penalised Cox via ncvsurv (competitor).
fit_ncvreg_pooled <- function(data_list, penalty = "SCAD",
                              lambda = NULL, nfolds = 5L) {
  combined  <- do.call(rbind, data_list)
  X_pool    <- as.matrix(combined[, -(1:3), drop = FALSE])
  surv_pool <- Surv(combined$time, combined$status)
  K         <- length(data_list)
  
  if (is.null(lambda)) {
    cv_fit <- tryCatch(
      cv.ncvsurv(X_pool, surv_pool, penalty = penalty, nfolds = nfolds),
      error = function(e) NULL)
    lambda <- if (!is.null(cv_fit)) cv_fit$lambda.min else 0.05
  }
  fit      <- ncvsurv(X_pool, surv_pool, penalty = penalty,
                      lambda = lambda, warn = FALSE)
  beta_hat <- as.numeric(coef(fit))
  list(
    alpha_hat  = beta_hat,
    eps_hat    = matrix(0, K, length(beta_hat)),
    theta_hat  = matrix(rep(beta_hat, K), nrow = K, byrow = TRUE),
    active_set = which(abs(beta_hat) > 1e-8),
    lambda     = lambda
  )
}


################################################################################
##  SECTION 4: PERFORMANCE METRICS                                             ##
################################################################################

#' Compute variable selection and estimation metrics.
#'
#' TPR        - True Positive Rate (sensitivity):  TP / s
#' TNR        - True Negative Rate (specificity):  TN / (p - s)
#' FDR        - False Discovery Rate:              FP / |selected|
#' F1         - 2*TP / (2*TP + FP + FN)
#' MCC        - Matthews Correlation Coefficient
#' alpha_bias - Mean |alpha_hat_j - alpha_j| over true signals
#' fp_bias    - Mean |alpha_hat_j| over true zeros
#' MSE_theta  - Mean squared error of study-specific theta_kj
#' exact_match- 1 if selected set == true active set
compute_metrics <- function(est, truth) {
  p        <- truth$p;  K <- truth$K
  active   <- truth$active_set
  inactive <- setdiff(seq_len(p), active)
  s        <- length(active)
  selected <- est$active_set
  
  TP <- length(intersect(selected, active))
  FP <- length(setdiff(selected, active))
  TN <- length(intersect(inactive, setdiff(seq_len(p), selected)))
  FN <- length(setdiff(active, selected))
  
  TPR <- TP / max(s, 1L)
  TNR <- TN / max(length(inactive), 1L)
  FDR <- FP / max(length(selected), 1L)
  F1  <- 2 * TP / max(2L * TP + FP + FN, 1L)
  
  denom_mcc <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  MCC       <- if (denom_mcc > 0) (TP * TN - FP * FN) / denom_mcc else 0
  
  alpha_bias <- if (s > 0L)
    mean(abs(est$alpha_hat[active] - truth$alpha_true[active])) else NA_real_
  fp_bias <- if (length(inactive) > 0L)
    mean(abs(est$alpha_hat[inactive])) else NA_real_
  
  mse_theta <- 0
  if (!is.null(est$theta_hat) && !is.null(truth$theta_true)) {
    for (k in seq_len(K))
      mse_theta <- mse_theta +
        mean((est$theta_hat[k, ] - truth$theta_true[[k]])^2)
    mse_theta <- mse_theta / K
  }
  
  data.frame(TPR = TPR, TNR = TNR, FDR = FDR, F1 = F1, MCC = MCC,
             alpha_bias = alpha_bias, fp_bias = fp_bias,
             MSE_theta = mse_theta,
             exact_match = as.integer(setequal(selected, active)))
}


################################################################################
##  SECTION 5: SIMULATION SCENARIOS                                            ##
################################################################################

#' Define all simulation scenarios.
#'
#' Factors varied (Breheny & Huang 2011; Fan & Li 2001; Zhang 2010 benchmarks):
#'   A: Correlation structure (rho = 0, 0.5, 0.8)
#'   B: Heterogeneity level  (epsilon_sd = 0, 0.3, 0.6)
#'   C: Censoring rate       (10%, 30%, 50%)
#'   D: Signal strength      (weak, moderate, strong)
#'   E: Dimensionality       (p = 50, 100, 200)
#'   F: Number of studies    (K = 5, 10)
#'   G: Unequal sample sizes
define_scenarios <- function() {
  base <- list(
    K = 5, n_per_study = 200, p = 50, s_alpha = 6, s_epsilon = 4,
    alpha_true = c(1.0, -0.8, 0.6, -0.5, 0.4, -0.3),
    epsilon_sd = 0.3, rho = 0.5, cens_rate = 0.30
  )
  list(
    S01_base     = base,
    S02_indep    = modifyList(base, list(rho = 0.0)),
    S03_highcor  = modifyList(base, list(rho = 0.8)),
    S04_homog    = modifyList(base, list(epsilon_sd = 0.0, s_epsilon = 0L)),
    S05_highhet  = modifyList(base, list(epsilon_sd = 0.6)),
    S06_lowcens  = modifyList(base, list(cens_rate = 0.10)),
    S07_highcens = modifyList(base, list(cens_rate = 0.50)),
    S08_weak     = modifyList(base, list(alpha_true = c(0.4, -0.3, 0.25, -0.2, 0.15, -0.10))),
    S09_strong   = modifyList(base, list(alpha_true = c(2.0, -1.5, 1.0, -0.8, 0.6, -0.4))),
    S10_highd    = modifyList(base, list(p = 100L)),
    S11_ultrad   = modifyList(base, list(p = 200L, n_per_study = 100L)),
    S12_moreK    = modifyList(base, list(K = 10L, n_per_study = 100L)),
    S13_unequal  = modifyList(base, list(n_per_study = c(100, 200, 300, 150, 250)))
  )
}


################################################################################
##  SECTION 6: SIMULATION RUNNER                                               ##
################################################################################

#' Run one replicate: generate data, fit all methods, compute metrics.
run_one_replicate <- function(
    scenario_params,
    methods      = c("scad", "mcp", "lasso", "enet", "scad_pool", "mcp_pool"),
    lambda_alpha = 0.08,
    lambda_eps   = 0.04,
    nfolds_cv    = 5L
) {
  dat     <- do.call(generate_ipd_data, scenario_params)
  results <- vector("list", length(methods))
  names(results) <- methods
  
  for (method in methods) {
    fit <- tryCatch({
      if (method %in% c("scad", "mcp")) {
        fit_hierarchical_meta_cox(
          dat$studies, penalty = method,
          lambda_alpha = lambda_alpha, lambda_eps = lambda_eps,
          nfolds_cv = nfolds_cv, max_iter = 30L, verbose = FALSE)
      } else if (method %in% c("lasso", "enet")) {
        fit_pooled_penalized(dat$studies, method = method)
      } else if (method %in% c("scad_pool", "mcp_pool")) {
        fit_ncvreg_pooled(dat$studies,
                          penalty = toupper(sub("_pool", "", method)))
      } else stop(paste("Unknown method:", method))
    }, error = function(e) {
      message(sprintf("  [WARN] '%s' failed: %s", method, conditionMessage(e)))
      NULL
    })
    
    if (!is.null(fit)) {
      m <- compute_metrics(fit, dat); m$method <- method
      results[[method]] <- m
    }
  }
  do.call(rbind, Filter(Negate(is.null), results))
}


#' Run full simulation: nsim replicates x scenarios x methods.
run_simulation_study <- function(
    nsim         = 100L,
    scenarios    = NULL,
    methods      = c("scad", "mcp", "lasso", "enet"),
    n_cores      = 1L,
    lambda_alpha = 0.08,
    lambda_eps   = 0.04,
    seed         = 2024L
) {
  if (is.null(scenarios)) scenarios <- define_scenarios()
  cat("================================================================\n")
  cat("  Hierarchical Nonconvex Penalised Meta-Cox Simulation\n")
  cat(sprintf("  nsim=%d | scenarios=%d | methods=%d | cores=%d\n",
              nsim, length(scenarios), length(methods), n_cores))
  cat("================================================================\n\n")
  
  worker_pkgs <- c("survival", "ncvreg", "glmnet", "MASS")
  worker_fns  <- c("run_one_replicate", "generate_ipd_data",
                   "fit_hierarchical_meta_cox", "fit_pooled_penalized",
                   "fit_ncvreg_pooled", "compute_metrics",
                   "build_pooled_data", "cv_select_lambda_ncv_glm")
  
  all_results <- vector("list", length(scenarios))
  names(all_results) <- names(scenarios)
  
  for (sc_name in names(scenarios)) {
    cat(sprintf(">>> Scenario: %s\n", sc_name))
    sc_params <- scenarios[[sc_name]]
    
    if (n_cores > 1L) {
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      rep_results <- foreach(
        rep = seq_len(nsim), .combine = rbind,
        .packages = worker_pkgs, .export = worker_fns
      ) %dopar% {
        set.seed(seed + rep)
        r <- run_one_replicate(sc_params, methods = methods,
                               lambda_alpha = lambda_alpha,
                               lambda_eps = lambda_eps)
        if (!is.null(r)) { r$rep <- rep; r$scenario <- sc_name; r }
      }
      stopCluster(cl)
    } else {
      rep_results <- NULL
      for (rep in seq_len(nsim)) {
        set.seed(seed + rep)
        if (rep %% 10L == 0L) cat(sprintf("  Replicate %d / %d\n", rep, nsim))
        r <- run_one_replicate(sc_params, methods = methods,
                               lambda_alpha = lambda_alpha,
                               lambda_eps = lambda_eps)
        if (!is.null(r)) {
          r$rep <- rep; r$scenario <- sc_name
          rep_results <- rbind(rep_results, r)
        }
      }
    }
    
    all_results[[sc_name]] <- rep_results
    cat(sprintf("  Done: %d replicates\n\n", nsim))
  }
  do.call(rbind, all_results)
}


################################################################################
##  SECTION 7: SUMMARY AND VISUALISATION                                       ##
################################################################################

summarize_results <- function(sim_results) {
  metrics_cols <- c("TPR", "TNR", "FDR", "F1", "MCC",
                    "alpha_bias", "fp_bias", "MSE_theta", "exact_match")
  sim_results %>%
    group_by(scenario, method) %>%
    summarise(
      across(all_of(metrics_cols),
             list(mean = ~mean(.x, na.rm = TRUE),
                  sd   = ~sd(.x,   na.rm = TRUE)),
             .names = "{.col}_{.fn}"),
      n_reps = n(), .groups = "drop")
}

print_comparison_table <- function(summary_df, scenario = "S01_base") {
  sub_df <- dplyr::filter(summary_df, scenario == !!scenario)
  cat(sprintf("\n=== Scenario: %s ===\n", scenario))
  cat(sprintf("%-12s %7s %7s %7s %10s %10s %8s\n",
              "Method", "TPR", "FDR", "MCC", "AlphaBias", "MSE_theta", "Exact%"))
  cat(strrep("-", 68), "\n")
  for (i in seq_len(nrow(sub_df))) {
    r <- sub_df[i, ]
    cat(sprintf("%-12s %7.3f %7.3f %7.3f %10.4f %10.4f %8.1f\n",
                r$method, r$TPR_mean, r$FDR_mean, r$MCC_mean,
                r$alpha_bias_mean, r$MSE_theta_mean, r$exact_match_mean * 100))
  }
  cat("\n")
}

plot_simulation_results <- function(sim_results, metric = "MCC",
                                    title = "Variable Selection Accuracy (MCC)") {
  method_colors <- c(
    scad = "#E41A1C", mcp = "#FF7F00",
    lasso = "#377EB8", enet = "#A65628",
    scad_pool = "#F781BF", mcp_pool = "#999999")
  method_labels <- c(
    scad = "SCAD (Hierarchical)", mcp = "MCP (Hierarchical)",
    lasso = "LASSO (Pooled)", enet = "Elastic Net (Pooled)",
    scad_pool = "SCAD (Pooled)", mcp_pool = "MCP (Pooled)")
  
  plot_df <- sim_results %>%
    group_by(scenario, method) %>%
    summarise(mean_val = mean(.data[[metric]], na.rm = TRUE),
              se_val   = sd(.data[[metric]], na.rm = TRUE) / sqrt(n()),
              .groups = "drop") %>%
    mutate(method_label = method_labels[method],
           method_type = ifelse(method %in% c("scad", "mcp"),
                                "Proposed", "Competitor"))
  
  ggplot(plot_df, aes(x = scenario, y = mean_val, color = method,
                      group = method, shape = method_type)) +
    geom_line(aes(linetype = method_type), linewidth = 0.9) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = mean_val - 1.96 * se_val,
                      ymax = mean_val + 1.96 * se_val),
                  width = 0.25, alpha = 0.5) +
    scale_color_manual(values = method_colors, labels = method_labels,
                       name = "Method") +
    scale_shape_manual(values = c(Proposed = 16L, Competitor = 17L),
                       name = "Type") +
    scale_linetype_manual(values = c(Proposed = "solid", Competitor = "dashed"),
                          name = "Type") +
    labs(title = title, x = "Scenario", y = metric,
         caption = "Error bars: 95% CI across replicates") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right",
          plot.title = element_text(face = "bold"))
}


################################################################################
##  SECTION 8: DIAGNOSTIC PLOTS                                                ##
################################################################################

plot_coefficient_comparison <- function(fit, truth, method_name = "SCAD") {
  K <- truth$K; p <- truth$p
  df <- do.call(rbind, lapply(seq_len(K), function(k)
    data.frame(covariate = seq_len(p), true_theta = truth$theta_true[[k]],
               est_theta = fit$theta_hat[k, ], study = paste0("Study ", k))))
  ggplot(df, aes(x = true_theta, y = est_theta, color = study)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
    geom_point(alpha = 0.6, size = 1.8) +
    facet_wrap(~study, nrow = 2L) +
    labs(title = sprintf("%s: Estimated vs True theta_kj", method_name),
         x = "True theta_kj", y = "Estimated theta_kj") +
    theme_bw(base_size = 11) + theme(legend.position = "none")
}

plot_heterogeneity <- function(fit, truth, top_j = 10L) {
  active <- sort(truth$active_set); K <- truth$K
  show_j <- head(active, top_j)
  df <- do.call(rbind, lapply(show_j, function(j)
    do.call(rbind, lapply(seq_len(K), function(k)
      data.frame(covariate = paste0("X", j), study = paste0("Study ", k),
                 True = truth$epsilon_true[[k]][j],
                 Estimated = fit$eps_hat[k, j])))))
  df_long <- tidyr::pivot_longer(df, c(True, Estimated),
                                 names_to = "type", values_to = "value")
  ggplot(df_long, aes(x = covariate, y = value, fill = type)) +
    geom_col(position = "dodge", alpha = 0.85) +
    geom_hline(yintercept = 0) +
    facet_wrap(~study, nrow = 2L) +
    scale_fill_manual(values = c(True = "#2C7BB6", Estimated = "#D7191C"),
                      name = "") +
    labs(title = "Study-Specific Deviations: True vs Estimated",
         x = "Covariate", y = "epsilon_kj") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


################################################################################
##                                                                            ##
##   SECTION 7: REAL DATA APPLICATION                                         ##
##   "Additive Hierarchical Variable Selection for Nonconvex Penalized        ##
##    Cox Models in Individual-Participant-Data Meta-Analysis"                ##
##                                                                            ##
##   Dataset:  curatedOvarianData                                             ##
##     Multi-cohort high-grade serous ovarian carcinoma microarray data.      ##
##     Outcome: overall survival (days → years).                              ##
##     Ganzfried et al. (2013), Bioinformatics.                               ##
##     Gene retention: top 2,000 most variable genes per cohort (by IQR),    ##
##     intersected across all 6 cohorts to yield a common high-dimensional   ##
##     feature set (p ≈ 1,500–3,000 genes).                                  ##
##                                                                            ##
##   Pipeline:                                                                ##
##     7.0  Packages, helpers, and shared estimation functions               ##
##     7.1  Data loading and cohort characteristics                           ##
##     7.2  Preprocessing (QC, standardisation, gene intersection)           ##
##     7.3  Method fitting and cross-validation tuning                        ##
##     7.4  Results: biomarker selection and reproducibility                  ##
##     7.5  Cross-study heterogeneity analysis                                ##
##                                                                            ##
##   OUTPUT FILES (written to Section7_outputs/):                             ##
##     Figures : Fig7_1_cohort_bar.png      Fig7_2_km_curves.png             ##
##               Fig7_3_forest.png          Fig7_4_deviation_tile.png        ##
##               Fig7_5_mcp_vs_scad.png     Fig7_6_selection_sizes.png       ##
##     Tables  : Table7_1_cohort_chars.csv  Table7_2_biomarkers.csv          ##
##               Table7_3_method_summary.csv                                  ##
##                                                                            ##
################################################################################


################################################################################
##  7.0  PACKAGES, HELPERS, AND SHARED ESTIMATION FUNCTIONS                  ##
################################################################################

## ── Install once (uncomment on first run) ────────────────────────────────────
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("curatedOvarianData", "Biobase"))
# install.packages(c("survival","ncvreg","glmnet","MASS","ggplot2","tidyr",
#                    "dplyr","patchwork","ggrepel","RColorBrewer","scales",
#                    "survminer","stringr","gridExtra"))

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
  library(Biobase)
  library(curatedOvarianData)
})

## ── Output directory ─────────────────────────────────────────────────────────
out_dir <- "Section7_outputs"
dir.create(out_dir, showWarnings = FALSE)
fig <- function(name) file.path(out_dir, name)

## ── Global plot theme ────────────────────────────────────────────────────────
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


## ── Build stacked data for pooled methods ────────────────────────────────────
build_pooled_data <- function(data_list, offset_list = NULL) {
  K <- length(data_list)
  X_list <- y_list <- off_list <- vector("list", K)
  for (k in seq_len(K)) {
    dk            <- data_list[[k]]
    Xk            <- as.matrix(dk[, -(1:3), drop = FALSE])
    X_list[[k]]   <- Xk
    y_list[[k]]   <- cbind(time = dk$time, status = dk$status)
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

## ── CV lambda selection ──────────────────────────────────────────────────────
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

## ── Hierarchical Meta-Cox fitter (Algorithm 1 in paper) ─────────────────────
fit_hierarchical_meta_cox <- function(
    data_list,
    penalty      = "mcp",
    lambda_alpha = NULL,
    lambda_eps   = NULL,
    max_iter     = 40L,
    conv_tol     = 1e-4,
    active_tol   = 1e-6,
    nfolds_cv    = 5L,
    verbose      = TRUE
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
    cat(sprintf("  [%s]  lambda_alpha = %s   lambda_eps = %s\n",
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
      cat(sprintf("    iter %2d | active = %3d | delta = %.2e\n",
                  iter, length(active_set), max_change))
    if (max_change < conv_tol) {
      if (verbose) cat(sprintf("    Converged at iter %d\n", iter))
      break
    }
    if (iter == max_iter && verbose)
      cat(sprintf("    Max iter reached (delta = %.2e)\n", max_change))
  }
  
  list(alpha_hat = alpha_hat, eps_hat = eps_hat, theta_hat = theta_hat,
       active_set        = which(alpha_hat != 0),
       lambda_alpha_used = lam_a_used,
       lambda_eps_used   = lam_e_used,
       n_iter            = iter)
}

## ── Pooled penalized Cox (LASSO / Elastic Net) ───────────────────────────────
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

## ── Pooled nonconvex Cox (SCAD / MCP) ────────────────────────────────────────
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
##  7.1  DATA LOADING AND COHORT CHARACTERISTICS                              ##
################################################################################

cat("\n================================================================\n")
cat("  7.1  Loading curatedOvarianData\n")
cat("================================================================\n\n")

## Candidate cohorts — six well-characterised HGSOC studies
candidate_eset_names <- c(
  "TCGA_eset",             # TCGA ovarian       (N ~ 489)
  "GSE32062.GPL6480_eset", # Japanese cohort    (N ~ 260)
  "GSE9891_eset",          # Australian cohort  (N ~ 276)
  "GSE26712_eset",         # Duke cohort        (N ~ 185)
  "GSE17260_eset",         # Tothill cohort     (N ~ 110)
  "GSE30161_eset"          # Phase III cohort   (N ~  58)
)

## ── Cohort extraction helper ─────────────────────────────────────────────────
#' Extract survival data and top variable genes from an ExpressionSet.
#'
#' @param eset_name     Name of the ExpressionSet object in curatedOvarianData
#' @param top_var_genes Number of most variable genes to retain per cohort
#'                      (ranked by IQR; set to 2000 for high-dimensional analysis)
extract_cohort <- function(eset_name, top_var_genes = 2000L) {
  tryCatch({
    data(list = eset_name, package = "curatedOvarianData", envir = environment())
    eset <- get(eset_name, envir = environment())
    pd   <- pData(eset)
    
    ## Identify survival columns (column names vary across datasets)
    time_col   <- intersect(c("days_to_death", "os_time", "survival_time",
                              "t.os", "T.OS"), colnames(pd))[1]
    status_col <- intersect(c("vital_status", "os_status", "censored",
                              "e.os", "E.OS"), colnames(pd))[1]
    if (is.na(time_col) || is.na(status_col)) return(NULL)
    
    time_raw   <- as.numeric(pd[[time_col]])
    status_raw <- pd[[status_col]]
    
    ## Harmonise event indicator to 0/1
    status_bin <- if (is.character(status_raw) || is.factor(status_raw))
      as.integer(tolower(as.character(status_raw)) %in%
                   c("deceased", "dead", "1", "died", "yes"))
    else as.integer(as.numeric(status_raw))
    
    ## Retain subjects with valid, positive follow-up
    keep <- !is.na(time_raw) & !is.na(status_bin) & time_raw > 0
    if (sum(keep) < 30 || sum(status_bin[keep]) < 10) return(NULL)
    
    ## Convert days to years; transpose to samples × genes
    time_use <- time_raw[keep] / 365.25
    expr_mat <- t(exprs(eset)[, keep])
    
    ## Select top genes by IQR (robust to outliers) then standardise
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

## ── Load all six cohorts retaining top 2,000 variable genes each ─────────────
cat("Loading cohorts (top_var_genes = 2000)...\n")
cohort_list <- Filter(Negate(is.null),
                      lapply(candidate_eset_names, extract_cohort,
                             top_var_genes = 2000L))

cat(sprintf("  Loaded %d / %d cohorts successfully\n",
            length(cohort_list), length(candidate_eset_names)))
for (ch in cohort_list)
  cat(sprintf("    %-35s  N = %3d  events = %3d  cens = %.0f%%\n",
              ch$eset_name, ch$n, ch$n_events, ch$cens_rate * 100))


################################################################################
##  7.2  PREPROCESSING: GENE INTERSECTION AND RE-STANDARDISATION             ##
################################################################################

cat("\n")
cat("================================================================\n")
cat("  7.2  Preprocessing\n")
cat("================================================================\n\n")

## Intersect gene sets across all cohorts so every cohort has the same features
gene_sets    <- lapply(cohort_list, `[[`, "gene_names")
common_genes <- Reduce(intersect, gene_sets)

cat(sprintf("  Genes per cohort (before intersection): %s\n",
            paste(sapply(gene_sets, length), collapse = ", ")))
cat(sprintf("  Common genes after intersection:  p = %d\n",
            length(common_genes)))

if (length(common_genes) < 200)
  warning("Fewer than 200 common genes — consider reducing the number of cohorts.")

## Subset each cohort to the common gene set and re-standardise within cohort
study_data_list <- lapply(cohort_list, function(ch) {
  df        <- ch$data
  keep_cols <- c("study_id", "time", "status",
                 intersect(common_genes, colnames(df)))
  df_sub    <- df[, keep_cols, drop = FALSE]
  ecols     <- keep_cols[-(1:3)]
  scaled    <- scale(df_sub[, ecols])
  scaled[is.nan(scaled)] <- 0          # replace NaN from zero-variance genes
  df_sub[, ecols] <- scaled
  df_sub
})

K        <- length(study_data_list)
p        <- length(common_genes)
N_total  <- sum(sapply(study_data_list, nrow))
n_events <- sum(sapply(study_data_list, function(d) sum(d$status)))

cat(sprintf("\n  Final dataset:  K = %d studies | N = %d patients | p = %d genes\n",
            K, N_total, p))
cat(sprintf("  Total events:   %d  (%.1f%%)\n",
            n_events, 100 * n_events / N_total))


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
print(table_cohort[, -7])


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
## FIGURE 7.2  Kaplan–Meier Curves by Cohort                                ##
################################################################################

km_data <- do.call(rbind, lapply(seq_len(K), function(k) {
  d <- study_data_list[[k]]
  d$study_label <- paste0("Study ", LETTERS[k])
  d
}))

km_fit <- survfit(Surv(time, status) ~ study_label, data = km_data)

p_km <- ggsurvplot(
  km_fit,
  data = km_data,
  palette = c(
    "#E41A1C",  # red
    "#377EB8",  # blue
    "#4DAF4A",  # green
    "#984EA3",  # purple
    "#FF7F00",  # orange
    "#A65628"   # brown
  ),
  conf.int          = FALSE,
  risk.table        = TRUE,
  risk.table.height = 0.28,
  xlab              = "Time (years)",
  ylab              = "Overall survival probability",
  title             = "Kaplan–Meier Curves by Cohort",
  legend.title      = "",
  ggtheme           = theme_paper(),
  fontsize          = 3.5
)

ggsave(
  fig("Fig7_2_km_curves.png"),
  print(p_km),
  width  = 9,
  height = 7
)

cat("Figure 7.2 saved.\n")


################################################################################
##  7.3  METHOD FITTING AND CROSS-VALIDATION TUNING                          ##
################################################################################

cat("\n================================================================\n")
cat("  7.3  Fitting all six methods\n")
cat("================================================================\n\n")
cat(sprintf("  Dataset: K=%d | N=%d | p=%d\n\n", K, N_total, p))

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
##  7.4  RESULTS: BIOMARKER SELECTION AND REPRODUCIBILITY                    ##
################################################################################

cat("\n================================================================\n")
cat("  7.4  Biomarker selection results\n")
cat("================================================================\n\n")

## ── TABLE 7.2  Biomarker comparison table ────────────────────────────────────
# For each gene selected by MCP (Hier.), report the estimated global
# log-hazard ratio from all six methods. Genes are ordered by the number
# of methods that selected them (reproducibility) and then by |alpha_hat|.

mcp_active <- fit_mcp_hier$active_set
if (length(mcp_active) > 0) {
  
  # Order by |alpha_hat| descending
  ord       <- order(abs(fit_mcp_hier$alpha_hat[mcp_active]), decreasing = TRUE)
  top_idx   <- mcp_active[ord]
  top_genes <- common_genes[top_idx]
  
  tbl2 <- data.frame(Gene = top_genes, stringsAsFactors = FALSE)
  for (nm in names(all_fits)) {
    vals <- round(all_fits[[nm]]$alpha_hat[top_idx], 4)
    vals[abs(vals) < 1e-8] <- 0
    tbl2[[nm]] <- vals
  }
  tbl2$N_methods <- rowSums(tbl2[, -1] != 0)
  tbl2 <- tbl2[order(-tbl2$N_methods, -abs(tbl2[["MCP (Hier.)"]])), ]
  
  write.csv(tbl2, fig("Table7_2_biomarkers.csv"), row.names = FALSE)
  cat("Table 7.2 saved.\n\n")
  cat("  Top reproducible biomarkers:\n")
  print(head(tbl2, 10))
  
} else {
  cat("  No genes selected by MCP (Hier.).\n")
  tbl2 <- NULL
}

## ── FIGURE 7.3  Forest plot: global + study-specific estimates ───────────────
if (length(mcp_active) > 0) {
  alpha_active <- fit_mcp_hier$alpha_hat[mcp_active]
  n_show   <- min(12, length(mcp_active))
  show_idx <- mcp_active[order(abs(alpha_active), decreasing = TRUE)][
    seq_len(n_show)]
  show_genes_f <- common_genes[show_idx]
  
  forest_df <- do.call(rbind, lapply(seq_along(show_idx), function(gi) {
    j    <- show_idx[gi]
    gene <- common_genes[j]
    rbind(
      do.call(rbind, lapply(seq_len(K), function(k)
        data.frame(gene     = gene,
                   label    = paste0("Study ", LETTERS[k]),
                   estimate = fit_mcp_hier$theta_hat[k, j],
                   type     = "Study-specific",
                   stringsAsFactors = FALSE))),
      data.frame(gene     = gene,
                 label    = "Global (meta)",
                 estimate = fit_mcp_hier$alpha_hat[j],
                 type     = "Global",
                 stringsAsFactors = FALSE))
  }))
  forest_df$gene  <- factor(forest_df$gene, levels = rev(show_genes_f))
  forest_df$label <- factor(forest_df$label,
                            levels = c(paste0("Study ", LETTERS[seq_len(K)]),
                                       "Global (meta)"))
  
  p_forest <- ggplot(forest_df,
                     aes(x = estimate, y = label,
                         colour = type, size = type, shape = type)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.85) +
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
      subtitle = "Red diamond = global meta-analytic effect; blue circles = study-specific \u03b8\u0302kj",
      x = "Log-hazard ratio", y = NULL) +
    theme_paper(base_size = 10) + theme(legend.position = "top")
  
  ggsave(fig("Fig7_3_forest.png"), p_forest,
         width = 11, height = max(6, ceiling(n_show / 3) * 3), dpi = 300)
  cat("\nFigure 7.3 saved.\n")
}


################################################################################
##  7.5  CROSS-STUDY HETEROGENEITY ANALYSIS                                  ##
################################################################################

cat("\n================================================================\n")
cat("  7.5  Cross-study heterogeneity\n")
cat("================================================================\n\n")

if (length(mcp_active) > 0) {
  
  alpha_active <- fit_mcp_hier$alpha_hat[mcp_active]
  active_genes <- common_genes[mcp_active]
  eps_active   <- fit_mcp_hier$eps_hat[, mcp_active, drop = FALSE]
  
  ## tau^2 analogue: cross-study variance of epsilon_kj
  tau2 <- apply(eps_active, 2, var)
  I2   <- tau2 / (tau2 + 1e-6)   # scale-free heterogeneity measure in [0,1]
  
  het_summary <- data.frame(
    gene      = active_genes,
    alpha_hat = alpha_active,
    tau2      = tau2,
    I2        = I2,
    het_level = cut(I2, c(-Inf, 0.25, 0.50, 0.75, Inf),
                    labels = c("Low", "Moderate", "Substantial", "High"))
  ) %>% arrange(desc(abs(alpha_hat)))
  
  cat("Heterogeneity summary (MCP Hierarchical):\n")
  print(het_summary)
  
  ## ── FIGURE 7.4  Deviation heatmap ────────────────────────────────────────
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
         width = max(6, K * 1.1 + 2),
         height = max(5, n_tile * 0.45 + 2), dpi = 300)
  cat("Figure 7.4 saved.\n")
  
  ## ── FIGURE 7.5  MCP vs SCAD study-specific estimate comparison ───────────
  n_genes   <- ncol(fit_mcp_hier$theta_hat)   # avoids name clash with ggplot2 p
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
                              "comparison across all studies and " * p ~ "genes"),
      x = expression("MCP (Hier.)  " * hat(theta)[kj]),
      y = expression("SCAD (Hier.) " * hat(theta)[kj])) +
    theme_paper(base_size = 10) + theme(legend.position = "top")
  
  ggsave(fig("Fig7_5_mcp_vs_scad.png"), p_cmp,
         width = max(7, ceiling(K / 3) * 3.5 + 1), height = 6, dpi = 300)
  cat("Figure 7.5 saved.\n")
}

## ── FIGURE 7.6  Selected set sizes ───────────────────────────────────────────
sel_sizes <- data.frame(
  Method   = names(all_fits),
  Selected = sapply(all_fits, function(f) length(f$active_set)),
  Type     = c("Proposed", "Proposed", "Competitor",
               "Competitor", "Competitor", "Competitor"),
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
##  TABLE 7.3  Method Summary Table                                           ##
################################################################################

## Reproducible core: genes selected by at least 4 of 6 methods
sel_bin <- sapply(all_fits, function(fit)
  as.integer(common_genes %in% common_genes[fit$active_set]))
rownames(sel_bin) <- common_genes
repro_genes <- rownames(sel_bin)[rowSums(sel_bin) >= 4]

tbl3 <- data.frame(
  Method = names(all_fits),
  p_total          = p,
  Selected_genes   = sapply(all_fits, function(f) length(f$active_set)),
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

cat("\n================================================================\n")
cat(sprintf("  Section 7 complete.  Outputs in: %s/\n", out_dir))
cat(sprintf("  Dataset: curatedOvarianData | K=%d | N=%d | p=%d\n",
            K, N_total, p))
cat("================================================================\n")