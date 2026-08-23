################################################################################
##                                                                            ##
##   SIMULATION CODE FOR:                                                     ##
##   "Additive Hierarchical Variable Selection for Nonconvex Penalized        ##
##    Cox Models in Individual-Participant-Data Meta-Analysis"                ##
##                                                                            ##
##   Authors: Emmanuel DJEGOU, Bertin DEHIGBE, Jarrad BOTCHWAY                ##
##                                                                            ##
##                                                                            ##
##   REFERENCES:                                                              ##
##     Fan & Li (2001, 2002)      - SCAD oracle properties, Cox models        ##
##     Zhang (2010)               - MCP nearly unbiased selection             ##
##     Breheny & Huang (2011)     - ncvreg / ncvsurv backends                 ##
##     Simon et al. (2011)        - Penalized Cox coordinate descent          ##
##     DerSimonian & Laird (1986) - Random-effects meta-analysis              ##
##     Crowther et al. (2012)     - IPD meta-analysis survival models         ##
##     Bender et al. (2005)       - Simulating survival times (inverse CDF)   ##
##                                                                            ##
################################################################################

# ── Packages ──────────────────────────────────────────────────────────────────
# install.packages(c("survival","ncvreg","glmnet","MASS",
#                    "foreach","doParallel","ggplot2","tidyr","dplyr"))
suppressPackageStartupMessages({
  library(survival)
  library(ncvreg)
  library(glmnet)
  library(MASS)
  library(foreach)
  library(doParallel)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
})

set.seed(2024L)


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
##  SECTION 9: DEMO RUN                                                        ##
################################################################################

cat("\n")
cat("================================================================\n")
cat("  DEMO RUN                                                      \n")
cat("  (single replicate + nsim=20 over 3 scenarios)                 \n")
cat("  Set run_full=TRUE in Section 10 for the full paper sim.       \n")
cat("================================================================\n\n")

# ── Single replicate demonstration ───────────────────────────────────────────
cat("--- Single replicate demonstration ---\n\n")
set.seed(42L)
demo_data <- generate_ipd_data(
  K = 10, n_per_study = 100, p = 50, s_alpha = 6, s_epsilon = 4,
  alpha_true = c(1.0, -0.8, 0.6, -0.5, 0.4, -0.3),
  epsilon_sd = 0.3, rho = 0.5, cens_rate = 0.30
)
cat(sprintf("K=%d | N=%d | p=%d | True active: {%s}\n",
            demo_data$K, sum(demo_data$n_per_study), demo_data$p,
            paste(demo_data$active_set, collapse = ",")))
cat(sprintf("Censoring rates: %s\n\n",
            paste(round(sapply(demo_data$studies,
                               function(d) 1 - mean(d$status)), 3),
                  collapse = ", ")))

# Proposed hierarchical methods
cat("Fitting SCAD (hierarchical)...\n")
fit_scad <- fit_hierarchical_meta_cox(demo_data$studies, penalty = "scad",
                                      lambda_alpha = 0.08, lambda_eps = 0.04,
                                      max_iter = 30L, verbose = TRUE)

cat("\nFitting MCP (hierarchical)...\n")
fit_mcp  <- fit_hierarchical_meta_cox(demo_data$studies, penalty = "mcp",
                                      lambda_alpha = 0.08, lambda_eps = 0.04,
                                      max_iter = 30L, verbose = TRUE)

# Competitor methods
cat("\nFitting competitor methods...\n")
fit_lasso     <- fit_pooled_penalized(demo_data$studies, method = "lasso")
fit_enet      <- fit_pooled_penalized(demo_data$studies, method = "enet")
fit_scad_pool <- fit_ncvreg_pooled(demo_data$studies, penalty = "SCAD")
fit_mcp_pool  <- fit_ncvreg_pooled(demo_data$studies, penalty = "MCP")

all_fits <- list(scad = fit_scad, mcp = fit_mcp,
                 lasso = fit_lasso, enet = fit_enet,
                 scad_pool = fit_scad_pool, mcp_pool = fit_mcp_pool)

metrics_df <- do.call(rbind, lapply(names(all_fits), function(m) {
  met <- compute_metrics(all_fits[[m]], demo_data); met$method <- m; met
}))

cat("\n\n=== Single Replicate Results ===\n")
cat(sprintf("True active set: {%s}\n\n", paste(demo_data$active_set, collapse = ",")))
cat(sprintf("%-12s | %5s | %5s | %5s | %9s | %9s | Selected\n",
            "Method", "TPR", "FDR", "MCC", "AlphaBias", "MSE_theta"))
cat(strrep("-", 85), "\n")
for (i in seq_len(nrow(metrics_df))) {
  r <- metrics_df[i, ]
  m <- all_fits[[r$method]]
  cat(sprintf("%-12s | %.3f | %.3f | %.3f | %.4f    | %.4f    | {%s}\n",
              r$method, r$TPR, r$FDR, r$MCC,
              ifelse(is.na(r$alpha_bias), 0, r$alpha_bias),
              r$MSE_theta,
              paste(sort(m$active_set), collapse = ",")))
}

cat("\nTrue global effects:\n")
cat(sprintf("  %s\n",
            paste(sprintf("j=%d: %+.2f", demo_data$active_set,
                          demo_data$alpha_true[demo_data$active_set]),
                  collapse = "  ")))

for (nm in c("scad", "mcp")) {
  cat(sprintf("\n%s estimated alpha (first 8 covariates):\n", toupper(nm)))
  cat(sprintf("  %s\n",
              paste(sprintf("j=%d:%+.3f", 1:8, all_fits[[nm]]$alpha_hat[1:8]),
                    collapse = "  ")))
}

cat("\nSCAD study-specific deviations (active covariates):\n")
A <- fit_scad$active_set
if (length(A) > 0L) {
  for (k in seq_len(demo_data$K))
    cat(sprintf("  Study %d: %s\n", k,
                paste(sprintf("j=%d:%+.3f", A, fit_scad$eps_hat[k, A]),
                      collapse = "  ")))
}

cat("\nMCP study-specific deviations (active covariates):\n")
B <- fit_mcp$active_set
if (length(B) > 0L) {
  for (k in seq_len(demo_data$K))
    cat(sprintf("  Study %d: %s\n", k,
                paste(sprintf("j=%d:%+.3f", B, fit_mcp$eps_hat[k, B]),
                      collapse = "  ")))
}

# ── Small simulation ──────────────────────────────────────────────────────────
cat("\n\n--- Small simulation (nsim=20, 3 scenarios) ---\n\n")
small_scenarios <- define_scenarios()[c("S01_base", "S04_homog", "S05_highhet")]

sim_results <- run_simulation_study(
  nsim         = 20L,
  scenarios    = small_scenarios,
  methods      = c("scad", "mcp", "lasso", "enet", "scad_pool", "mcp_pool"),
  n_cores      = 1L,
  lambda_alpha = 0.08,
  lambda_eps   = 0.04,
  seed         = 2024L
)

summary_df <- summarize_results(sim_results)
for (sc in unique(summary_df$scenario)) print_comparison_table(summary_df, sc)

# Plots
cat("Saving plots...\n")
ggsave("plot_MCC.png",
       plot_simulation_results(sim_results, "MCC", "MCC by Scenario"),
       width = 10, height = 6)
ggsave("plot_bias.png",
       plot_simulation_results(sim_results, "alpha_bias", "Alpha Bias by Scenario"),
       width = 10, height = 6)
ggsave("plot_FDR.png",
       plot_simulation_results(sim_results, "FDR", "FDR by Scenario"),
       width = 10, height = 6)
ggsave("plot_MSE.png",
       plot_simulation_results(sim_results, "MSE_theta", "MSE(theta) by Scenario"),
       width = 10, height = 6)
ggsave("plot_coef_SCAD.png",
       plot_coefficient_comparison(fit_scad, demo_data, "SCAD (Hierarchical)"),
       width = 10, height = 7)
ggsave("plot_coef_MCP.png",
       plot_coefficient_comparison(fit_mcp, demo_data, "MCP (Hierarchical)"),
       width = 10, height = 7)
ggsave("plot_het_SCAD.png",
       plot_heterogeneity(fit_scad, demo_data),
       width = 10, height = 7)
ggsave("plot_het_MCP.png",
       plot_heterogeneity(fit_mcp, demo_data),
       width = 10, height = 7)
cat("All plots saved.\n")


################################################################################
##  SECTION 10: FULL PAPER SIMULATION                                          ##
################################################################################

run_full <- FALSE   # <-- Set to FALSE to skip the full simulation

if (run_full) {
  cat("\n\n=== FULL PAPER SIMULATION (nsim=500, all scenarios) ===\n\n")
  n_cores_full <- 1L
  cat(sprintf("Using %d parallel cores\n\n", n_cores_full))
  
  full_results <- run_simulation_study(
    nsim         = 500L,
    scenarios    = define_scenarios(),
    methods      = c("scad", "mcp", "lasso", "enet", "scad_pool", "mcp_pool"),
    n_cores      = n_cores_full,
    lambda_alpha = NULL,   # NULL = BIC/CV per replicate (principled)
    lambda_eps   = NULL,
    seed         = 2024L
  )
  
  saveRDS(full_results, "sim_results_full.rds")
  write.csv(full_results, "sim_results_full.csv", row.names = FALSE)
  
  full_summary <- summarize_results(full_results)
  write.csv(full_summary, "sim_summary_full.csv", row.names = FALSE)
  for (sc in names(define_scenarios())) print_comparison_table(full_summary, sc)
  
  for (metric in c("MCC", "TPR", "FDR", "F1", "alpha_bias",
                   "fp_bias", "MSE_theta", "exact_match")) {
    ggsave(sprintf("plot_full_%s.pdf", metric),
           plot_simulation_results(full_results, metric,
                                   paste(metric, "across Scenarios")),
           width = 14, height = 7)
  }
  cat("Full simulation complete.\n")
}

cat("\n=== Script complete ===\n")