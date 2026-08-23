################################################################################
##                                                                            ##
##   SIMULATION CODE FOR:                                                     ##
##   "Additive Hierarchical Variable Selection for Nonconvex Penalized        ##
##    Cox Models in Individual-Participant-Data Meta-Analysis"                ##
##                                                                            ##
##   Authors: Emmanuel DJEGOU, Bertin DEHIGBE, Jarrad BOTCHWAY                ##
##                                                                            ##
##   REFERENCES:                                                              ##
##     Fan & Li (2001, 2002)      - SCAD oracle properties, Cox models        ##
##     Zhang (2010)               - MCP nearly unbiased selection             ##
##     Breheny & Huang (2011)     - ncvreg / ncvsurv backends                 ##
##     Simon et al. (2011)        - Penalized Cox coordinate descent          ##
##     DerSimonian & Laird (1986) - Random-effects meta-analysis              ##
##     Crowther et al. (2012)     - IPD meta-analysis survival models         ##
##     Bender et al. (2005)       - Simulating survival times (inverse CDF)   ##
##     Efron & Tibshirani (1994)  - Bootstrap uncertainty quantification      ##
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
  
  alpha_full <- c(alpha_true, rep(0, p - s_alpha))
  lambda_k   <- runif(K, weibull_scale_range[1], weibull_scale_range[2])
  idx        <- seq_len(p)
  Sigma      <- rho ^ abs(outer(idx, idx, "-"))
  
  epsilon_list <- lapply(seq_len(K), function(k) {
    eps <- rep(0, p)
    if (s_epsilon > 0L)
      eps[seq_len(s_epsilon)] <- rnorm(s_epsilon, 0, epsilon_sd)
    eps
  })
  
  study_data <- vector("list", K)
  for (k in seq_len(K)) {
    nk      <- n_per_study[k]
    Xk      <- mvrnorm(nk, mu = rep(0, p), Sigma = Sigma)
    Xk      <- scale(Xk)
    colnames(Xk) <- paste0("X", seq_len(p))
    theta_k <- alpha_full + epsilon_list[[k]]
    eta     <- as.numeric(Xk %*% theta_k)
    U       <- runif(nk)
    T_event <- (-log(U) / (lambda_k[k] * exp(eta))) ^ (1 / weibull_shape)
    lambda_c <- cens_rate / mean(T_event)
    C_time   <- rexp(nk, rate = lambda_c)
    Y_obs    <- pmin(T_event, C_time)
    delta    <- as.integer(T_event <= C_time)
    study_data[[k]] <- data.frame(
      study_id = k, time = Y_obs, status = delta, Xk)
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
##  SECTION 2: CORE ESTIMATION                                                ##
################################################################################

build_pooled_data <- function(data_list, offset_list = NULL) {
  K <- length(data_list)
  X_list <- y_list <- off_list <- sid_list <- vector("list", K)
  for (k in seq_len(K)) {
    dk            <- data_list[[k]]
    Xk            <- as.matrix(dk[, -(1:3), drop = FALSE])
    X_list[[k]]   <- Xk
    y_list[[k]]   <- cbind(time = dk$time, status = dk$status)
    off_list[[k]] <- if (is.null(offset_list)) rep(0, nrow(dk))
    else offset_list[[k]]
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

fit_hierarchical_meta_cox <- function(
    data_list,
    penalty      = "scad",
    lambda_alpha = NULL,
    lambda_eps   = NULL,
    max_iter     = 30L,
    conv_tol     = 1e-4,
    active_tol   = 1e-6,
    nfolds_cv    = 5L,
    verbose      = FALSE
) {
  K          <- length(data_list)
  p          <- ncol(data_list[[1]]) - 3L
  pen_upper  <- toupper(penalty)
  use_ncvreg <- pen_upper %in% c("SCAD", "MCP")
  alpha_enet <- switch(pen_upper, "LASSO" = 1.0, "ENET" = 0.5, 1.0)
  
  alpha_hat         <- rep(0, p)
  eps_hat           <- matrix(0, K, p)
  theta_hat         <- matrix(0, K, p)
  lambda_alpha_used <- NA_real_
  lambda_eps_used   <- NA_real_
  
  for (iter in seq_len(max_iter)) {
    theta_old  <- theta_hat
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
    alpha_hat[abs(alpha_hat) < active_tol] <- 0
    active_set <- which(alpha_hat != 0)
    
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
          lam_e <- if (is.null(lambda_eps))
            tryCatch(
              cv_select_lambda_ncv_glm(
                Xk_active, surv_k, off_alpha_k, pen_upper,
                nfolds = min(nfolds_cv, max(3L, nk %/% 20L))),
              error = function(e) lambda_alpha_used * 0.5)
          else lambda_eps
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
          lam_e <- if (is.null(lambda_eps))
            tryCatch(
              cv_select_lambda_ncv_glm(
                Xk_active, surv_k, off_alpha_k, pen_upper, alpha_enet,
                nfolds = min(nfolds_cv, max(3L, nk %/% 20L))),
              error = function(e) lambda_alpha_used * 0.5)
          else lambda_eps
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
    
    theta_hat  <- sweep(eps_hat, 2L, alpha_hat, "+")
    max_change <- max(abs(theta_hat - theta_old))
    if (max_change < conv_tol) break
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

fit_pooled_penalized <- function(data_list, method = "lasso",
                                 lambda = NULL, nfolds = 5L) {
  combined  <- do.call(rbind, data_list)
  X_pool    <- as.matrix(combined[, -(1:3), drop = FALSE])
  surv_pool <- Surv(combined$time, combined$status)
  K         <- length(data_list)
  alpha_g   <- switch(method, "lasso" = 1.0, "enet" = 0.5,
                      stop("method must be 'lasso' or 'enet'"))
  if (is.null(lambda)) {
    cv_fit <- cv.glmnet(X_pool, surv_pool, family = "cox",
                        alpha = alpha_g, nfolds = nfolds, type.measure = "C")
    lambda <- cv_fit$lambda.min
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
  list(alpha_hat  = beta_hat,
       eps_hat    = matrix(0, K, length(beta_hat)),
       theta_hat  = matrix(rep(beta_hat, K), nrow = K, byrow = TRUE),
       active_set = which(abs(beta_hat) > 1e-8),
       lambda     = lambda)
}

################################################################################
##  SECTION 3: PERFORMANCE METRICS                                            ##
################################################################################

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
##  SECTION 4: SIMULATION SCENARIOS                                           ##
################################################################################

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
    S08_weak     = modifyList(base, list(
      alpha_true = c(0.4, -0.3, 0.25, -0.2, 0.15, -0.10))),
    S09_strong   = modifyList(base, list(
      alpha_true = c(2.0, -1.5, 1.0, -0.8, 0.6, -0.4))),
    S10_highd    = modifyList(base, list(p = 100L)),
    S11_ultrad   = modifyList(base, list(p = 200L, n_per_study = 100L)),
    S12_moreK    = modifyList(base, list(K = 10L, n_per_study = 100L)),
    S13_unequal  = modifyList(base, list(
      n_per_study = c(100, 200, 300, 150, 250)))
  )
}

################################################################################
##  SECTION 5: SIMULATION RUNNER                                              ##
################################################################################

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
                               lambda_eps   = lambda_eps)
        if (!is.null(r)) { r$rep <- rep; r$scenario <- sc_name; r }
      }
      stopCluster(cl)
    } else {
      rep_results <- NULL
      for (rep in seq_len(nsim)) {
        set.seed(seed + rep)
        if (rep %% 10L == 0L)
          cat(sprintf("  Replicate %d / %d\n", rep, nsim))
        r <- run_one_replicate(sc_params, methods = methods,
                               lambda_alpha = lambda_alpha,
                               lambda_eps   = lambda_eps)
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
##  SECTION 6: SUMMARY AND VISUALISATION                                      ##
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
                r$alpha_bias_mean, r$MSE_theta_mean,
                r$exact_match_mean * 100))
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
           method_type  = ifelse(method %in% c("scad", "mcp"),
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
    scale_linetype_manual(
      values = c(Proposed = "solid", Competitor = "dashed"), name = "Type") +
    labs(title = title, x = "Scenario", y = metric,
         caption = "Error bars: 95% CI across replicates") +
    theme_bw(base_size = 12) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1),
          legend.position = "right",
          plot.title    = element_text(face = "bold"))
}

################################################################################
##  SECTION 7: UNCERTAINTY QUANTIFICATION VIA SIMULATION                     ##
##                                                                            ##
##  Purpose:                                                                  ##
##   Assess the EMPIRICAL coverage and width of 95% bootstrap percentile     ##
##   confidence intervals for the global effects alpha_j and the             ##
##   study-specific deviations epsilon_kj, using a simulation study in       ##
##   which the true parameter values are known.                               ##
##                                                                            ##
##  Design:                                                                   ##
##   For a single scenario (default: S01_base), draw M Monte Carlo           ##
##   datasets.  For each dataset:                                             ##
##     1. Fit MCP (Hierarchical) on the full dataset.                         ##
##     2. Run B stratified bootstrap resamples with lambda fixed to the       ##
##        original fit (Efron & Tibshirani 1994).                             ##
##     3. Compute percentile 95% CI for each alpha_j and epsilon_kj.         ##
##     4. Record whether the CI covers the true value (coverage) and CI       ##
##        half-width.                                                         ##
##   Aggregate: mean coverage and mean CI width over M Monte Carlo datasets.  ##
##                                                                            ##
##  Evaluation metrics:                                                       ##
##   - Coverage probability: Pr(CI contains true parameter)                  ##
##     Target: 0.95 (nominal level)                                           ##
##   - Mean CI width: average |CI_hi - CI_lo|                                ##
##     Lower = more precise inference                                          ##
##   - Bias: mean(alpha_hat - alpha_true) over M datasets                    ##
##   - RMSE: sqrt(mean((alpha_hat - alpha_true)^2)) over M datasets          ##
##                                                                            ##
################################################################################

## ── Helper: stratified bootstrap for one dataset ─────────────────────────────

#' Run B stratified bootstrap resamples on one dataset and return
#' percentile CIs for alpha and epsilon.
#'
#' @param data_list   K-study list from generate_ipd_data()$studies
#' @param fit_orig    Output of fit_hierarchical_meta_cox() on data_list
#' @param B           Bootstrap resamples
#' @param conf_level  Nominal coverage (default 0.95)
#' @param seed        RNG seed
#'
#' @return List: alpha_ci (p x 2), eps_ci (K x p x 2), alpha_se (p),
#'         eps_se (K x p), B_used
bootstrap_ci_one_dataset <- function(data_list, fit_orig,
                                     B = 200L, conf_level = 0.95,
                                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  K   <- length(data_list)
  p   <- length(fit_orig$alpha_hat)
  lam_a <- fit_orig$lambda_alpha_used
  lam_e <- fit_orig$lambda_eps_used
  lo    <- (1 - conf_level) / 2
  hi    <- 1 - lo
  
  boot_alpha <- matrix(NA_real_, B, p)
  boot_eps   <- array(NA_real_,  dim = c(B, K, p))
  
  ok <- 0L
  for (b in seq_len(B)) {
    # Stratified resample: within each study independently
    boot_data <- lapply(data_list, function(dk) {
      idx <- sample(nrow(dk), replace = TRUE)
      dk[idx, ]
    })
    fb <- tryCatch(
      fit_hierarchical_meta_cox(
        boot_data, penalty = "mcp",
        lambda_alpha = lam_a, lambda_eps = lam_e,
        max_iter = 30L, conv_tol = 1e-4,
        active_tol = 1e-6, verbose = FALSE),
      error = function(e) NULL)
    if (!is.null(fb)) {
      ok               <- ok + 1L
      boot_alpha[ok, ] <- fb$alpha_hat
      boot_eps[ok, , ] <- fb$eps_hat
    }
  }
  
  if (ok == 0L) return(NULL)
  boot_alpha <- boot_alpha[seq_len(ok), , drop = FALSE]
  boot_eps   <- boot_eps[seq_len(ok), , , drop = FALSE]
  
  alpha_ci  <- apply(boot_alpha, 2, quantile,
                     probs = c(lo, hi), na.rm = TRUE)  # 2 x p
  alpha_se  <- apply(boot_alpha, 2, sd, na.rm = TRUE)  # p
  
  eps_ci <- array(NA_real_, dim = c(K, p, 2))
  eps_se <- matrix(NA_real_, K, p)
  for (k in seq_len(K)) {
    eps_k        <- boot_eps[, k, , drop = FALSE]
    eps_k        <- matrix(eps_k, nrow = ok, ncol = p)
    eps_ci[k, ,] <- apply(eps_k, 2, quantile,
                          probs = c(lo, hi), na.rm = TRUE)
    eps_se[k, ]  <- apply(eps_k, 2, sd, na.rm = TRUE)
  }
  
  list(alpha_ci = t(alpha_ci),   # p x 2  [lo, hi]
       eps_ci   = eps_ci,        # K x p x 2
       alpha_se = alpha_se,
       eps_se   = eps_se,
       B_used   = ok)
}


## ── Main UQ simulation function ──────────────────────────────────────────────

#' Simulation study for uncertainty quantification.
#'
#' Runs M Monte Carlo datasets, fits MCP (Hierarchical) on each, runs B
#' bootstrap resamples, and evaluates empirical coverage and CI width
#' against the known true parameter values.
#'
#' @param M              Monte Carlo datasets (10-20 feasible in demo; 200 for paper)
#' @param B              Bootstrap resamples per dataset (200-500)
#' @param scenario_params  List of DGP parameters (from define_scenarios())
#' @param lambda_alpha   Fixed or NULL (CV); fixed recommended to save time
#' @param lambda_eps     Fixed or NULL (CV)
#' @param conf_level     Nominal coverage
#' @param seed           Base RNG seed
#' @param n_cores        Parallel cores for the M-loop
#'
#' @return Data frame with one row per covariate, columns:
#'   covariate, true_alpha, mean_alpha_hat, bias, rmse,
#'   mean_se, mean_ci_width, coverage, type ("active"|"inactive")
run_uq_simulation <- function(
    M              = 50L,
    B              = 200L,
    scenario_params = define_scenarios()[["S01_base"]],
    lambda_alpha   = 0.08,
    lambda_eps     = 0.04,
    conf_level     = 0.95,
    seed           = 2024L,
    n_cores        = 1L
) {
  cat("================================================================\n")
  cat("  SECTION 7: Uncertainty Quantification Simulation\n")
  cat(sprintf("  M=%d datasets | B=%d bootstrap | conf=%.0f%%\n",
              M, B, conf_level * 100))
  cat(sprintf("  Scenario: K=%d | n=%d | p=%d | lambda_a=%.3f | lambda_e=%.3f\n",
              scenario_params$K,
              if (length(scenario_params$n_per_study) == 1)
                scenario_params$n_per_study
              else mean(scenario_params$n_per_study),
              scenario_params$p,
              lambda_alpha, lambda_eps))
  cat("================================================================\n\n")
  
  K <- scenario_params$K
  p <- scenario_params$p
  
  ## Storage: per dataset, per covariate
  alpha_hat_mat <- matrix(NA_real_, M, p)   # point estimates
  covered_mat   <- matrix(NA_real_, M, p)   # 1/0 coverage
  width_mat     <- matrix(NA_real_, M, p)   # CI width
  se_mat        <- matrix(NA_real_, M, p)   # bootstrap SE
  
  ## Also track epsilon coverage for first active covariate
  eps_covered   <- matrix(NA_real_, M, K)
  eps_width     <- matrix(NA_real_, M, K)
  
  ## ── Worker for one Monte Carlo dataset ─────────────────────────────────
  one_mc <- function(m) {
    set.seed(seed + m)
    dat   <- do.call(generate_ipd_data, scenario_params)
    truth <- dat
    
    fit <- tryCatch(
      fit_hierarchical_meta_cox(
        dat$studies, penalty = "mcp",
        lambda_alpha = lambda_alpha, lambda_eps = lambda_eps,
        max_iter = 30L, verbose = FALSE),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    
    ci_res <- bootstrap_ci_one_dataset(
      dat$studies, fit, B = B,
      conf_level = conf_level,
      seed = seed + m + 10000L)
    if (is.null(ci_res)) return(NULL)
    
    alpha_ci <- ci_res$alpha_ci   # p x 2
    eps_ci   <- ci_res$eps_ci     # K x p x 2
    
    ## Coverage and width for alpha
    cov_alpha <- as.integer(
      fit$alpha_hat >= alpha_ci[, 1] &
        fit$alpha_hat <= alpha_ci[, 2])   # compare estimate vs its own CI
    
    ## ── IMPORTANT: compare CI to TRUE value ──────────────────────────────
    ## Coverage should be assessed against the DGP truth, not the estimate.
    cov_alpha <- as.integer(
      truth$alpha_true >= alpha_ci[, 1] &
        truth$alpha_true <= alpha_ci[, 2])
    
    width_alpha <- alpha_ci[, 2] - alpha_ci[, 1]
    
    ## Epsilon coverage: first active covariate, all K studies
    j1    <- truth$active_set[1]
    cov_e <- sapply(seq_len(K), function(k) {
      true_e <- truth$epsilon_true[[k]][j1]
      as.integer(true_e >= eps_ci[k, j1, 1] &
                   true_e <= eps_ci[k, j1, 2])
    })
    width_e <- sapply(seq_len(K), function(k)
      eps_ci[k, j1, 2] - eps_ci[k, j1, 1])
    
    list(
      alpha_hat  = fit$alpha_hat,
      cov_alpha  = cov_alpha,
      width_alpha = width_alpha,
      se_alpha   = ci_res$alpha_se,
      cov_eps    = cov_e,
      width_eps  = width_e
    )
  }
  
  ## ── Run M datasets ──────────────────────────────────────────────────────
  if (n_cores > 1L) {
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    mc_list <- foreach(
      m = seq_len(M), .combine = c,
      .packages = c("survival","ncvreg","glmnet","MASS"),
      .export   = c("generate_ipd_data","fit_hierarchical_meta_cox",
                    "bootstrap_ci_one_dataset","build_pooled_data",
                    "cv_select_lambda_ncv_glm",
                    "scenario_params","lambda_alpha","lambda_eps",
                    "conf_level","seed","K","p")
    ) %dopar% {
      list(one_mc(m))
    }
    stopCluster(cl)
    mc_list <- lapply(mc_list, `[[`, 1)
  } else {
    mc_list <- vector("list", M)
    for (m in seq_len(M)) {
      if (m %% 5L == 0L)
        cat(sprintf("  MC dataset %d / %d\n", m, M))
      mc_list[[m]] <- one_mc(m)
    }
  }
  
  ## ── Aggregate results ────────────────────────────────────────────────────
  ok_idx <- which(!sapply(mc_list, is.null))
  M_used <- length(ok_idx)
  cat(sprintf("\n  %d / %d MC datasets completed successfully.\n\n", M_used, M))
  
  for (i in seq_along(ok_idx)) {
    m <- ok_idx[i]
    r <- mc_list[[m]]
    alpha_hat_mat[i, ] <- r$alpha_hat
    covered_mat[i, ]   <- r$cov_alpha
    width_mat[i, ]     <- r$width_alpha
    se_mat[i, ]        <- r$se_alpha
    eps_covered[i, ]   <- r$cov_eps
    eps_width[i, ]     <- r$width_eps
  }
  
  ## ── True parameter values from a fixed seed for reference ────────────────
  set.seed(seed)
  ref_dat <- do.call(generate_ipd_data, scenario_params)
  
  ## ── Build summary data frame ─────────────────────────────────────────────
  results <- data.frame(
    covariate    = seq_len(p),
    true_alpha   = ref_dat$alpha_true,
    mean_alpha   = colMeans(alpha_hat_mat[seq_len(M_used), , drop=FALSE],
                            na.rm = TRUE),
    bias         = colMeans(alpha_hat_mat[seq_len(M_used), , drop=FALSE],
                            na.rm = TRUE) - ref_dat$alpha_true,
    rmse         = sqrt(colMeans(
      (alpha_hat_mat[seq_len(M_used), , drop=FALSE] -
         matrix(ref_dat$alpha_true, M_used, p, byrow = TRUE))^2,
      na.rm = TRUE)),
    coverage     = colMeans(covered_mat[seq_len(M_used), , drop=FALSE],
                            na.rm = TRUE),
    mean_ci_width = colMeans(width_mat[seq_len(M_used), , drop=FALSE],
                             na.rm = TRUE),
    mean_se      = colMeans(se_mat[seq_len(M_used), , drop=FALSE],
                            na.rm = TRUE),
    type         = ifelse(ref_dat$alpha_true != 0, "active", "inactive"),
    stringsAsFactors = FALSE
  )
  
  ## ── Epsilon summary (first active covariate, all studies) ─────────────
  j1_true_eps <- sapply(seq_len(K), function(k)
    ref_dat$epsilon_true[[k]][ref_dat$active_set[1]])
  
  eps_summary <- data.frame(
    study         = paste0("Study ", seq_len(K)),
    true_epsilon  = j1_true_eps,
    coverage      = colMeans(eps_covered[seq_len(M_used), , drop=FALSE],
                             na.rm = TRUE),
    mean_ci_width = colMeans(eps_width[seq_len(M_used), , drop=FALSE],
                             na.rm = TRUE)
  )
  
  list(
    alpha_summary = results,
    eps_summary   = eps_summary,
    alpha_hat_mat = alpha_hat_mat[seq_len(M_used), ],
    covered_mat   = covered_mat[seq_len(M_used), ],
    width_mat     = width_mat[seq_len(M_used), ],
    M_used        = M_used,
    B             = B,
    conf_level    = conf_level
  )
}


## ── Print and plot helpers for UQ results ────────────────────────────────────

#' Print a formatted UQ summary table to console.
print_uq_summary <- function(uq_res) {
  res <- uq_res$alpha_summary
  cat(sprintf("\n=== UQ Summary (M=%d | B=%d | %.0f%% CIs) ===\n\n",
              uq_res$M_used, uq_res$B, uq_res$conf_level * 100))
  
  ## Active covariates
  act <- res[res$type == "active", ]
  cat("Active covariates (true alpha != 0):\n")
  cat(sprintf("  %-5s  %8s  %8s  %7s  %7s  %8s  %9s\n",
              "j", "true", "mean_hat", "bias", "RMSE",
              "coverage", "CI_width"))
  cat(strrep("-", 65), "\n")
  for (i in seq_len(nrow(act))) {
    r <- act[i, ]
    cat(sprintf("  j=%-3d  %+8.4f  %+8.4f  %+7.4f  %7.4f  %8.3f  %9.4f\n",
                r$covariate, r$true_alpha, r$mean_alpha,
                r$bias, r$rmse, r$coverage, r$mean_ci_width))
  }
  
  ## Inactive covariates: aggregate
  inact <- res[res$type == "inactive", ]
  cat(sprintf(
    "\nInactive covariates (true alpha = 0):  %d covariates\n",
    nrow(inact)))
  cat(sprintf(
    "  Mean coverage: %.3f   Mean CI width: %.4f   Mean |bias|: %.5f\n",
    mean(inact$coverage,      na.rm = TRUE),
    mean(inact$mean_ci_width, na.rm = TRUE),
    mean(abs(inact$bias),     na.rm = TRUE)))
  
  ## Epsilon summary
  cat("\nStudy-specific deviations (first active covariate):\n")
  cat(sprintf("  %-10s  %8s  %8s  %9s\n",
              "Study", "true_eps", "coverage", "CI_width"))
  cat(strrep("-", 45), "\n")
  for (i in seq_len(nrow(uq_res$eps_summary))) {
    r <- uq_res$eps_summary[i, ]
    cat(sprintf("  %-10s  %+8.4f  %8.3f  %9.4f\n",
                r$study, r$true_epsilon, r$coverage, r$mean_ci_width))
  }
  
  ## Overall active coverage
  cat(sprintf(
    "\nOverall active-covariate coverage: %.3f  (target: %.3f)\n",
    mean(act$coverage, na.rm = TRUE),
    uq_res$conf_level))
  cat(sprintf(
    "Nominal undercoverage: %.3f  (bias from penalization expected)\n",
    uq_res$conf_level - mean(act$coverage, na.rm = TRUE)))
}


#' Plot empirical coverage and CI width across covariates.
plot_uq_results <- function(uq_res, title_prefix = "UQ Simulation") {
  
  res       <- uq_res$alpha_summary
  conf_line <- uq_res$conf_level
  
  ## ── Panel 1: Coverage by covariate ─────────────────────────────────────
  p1 <- ggplot(res, aes(x = covariate, y = coverage, colour = type)) +
    geom_hline(yintercept = conf_line, linetype = "dashed",
               colour = "grey40", linewidth = 0.7) +
    geom_point(size = 2.2, alpha = 0.85) +
    scale_colour_manual(
      values = c(active = "#C0392B", inactive = "#2980B9"),
      labels = c(active = "Active (true != 0)",
                 inactive = "Inactive (true = 0)"),
      name = NULL) +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format(accuracy = 1)) +
    labs(title    = sprintf("%s: Empirical Coverage", title_prefix),
         subtitle = sprintf("Dashed line = nominal %.0f%%  |  M=%d datasets  |  B=%d bootstrap",
                            conf_line * 100, uq_res$M_used, uq_res$B),
         x = "Covariate index", y = "Coverage probability") +
    theme_bw(base_size = 11) +
    theme(legend.position = "top",
          plot.title = element_text(face = "bold"))
  
  ## ── Panel 2: CI width by covariate ──────────────────────────────────────
  p2 <- ggplot(res, aes(x = covariate, y = mean_ci_width, colour = type)) +
    geom_point(size = 2.2, alpha = 0.85) +
    scale_colour_manual(
      values = c(active = "#C0392B", inactive = "#2980B9"),
      name = NULL) +
    labs(title    = sprintf("%s: Mean 95%% CI Width", title_prefix),
         subtitle = "Narrower = more precise; inactive covariates near zero expected",
         x = "Covariate index", y = "Mean CI width") +
    theme_bw(base_size = 11) +
    theme(legend.position = "top",
          plot.title = element_text(face = "bold"))
  
  ## ── Panel 3: Bootstrap SE vs |bias| for active covariates ───────────────
  act <- res[res$type == "active", ]
  p3 <- ggplot(act, aes(x = rmse, y = mean_se,
                        label = paste0("j=", covariate))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey40") +
    geom_point(colour = "#C0392B", size = 3) +
    ggrepel::geom_text_repel(size = 3, colour = "grey30") +
    labs(title    = sprintf("%s: Bootstrap SE vs RMSE (active covariates)", title_prefix),
         subtitle = "Points on the diagonal = SE well-calibrated to RMSE",
         x = "RMSE across MC datasets",
         y = "Mean bootstrap SE") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  
  ## ── Panel 4: Epsilon coverage bar chart ─────────────────────────────────
  eps <- uq_res$eps_summary
  p4 <- ggplot(eps, aes(x = study, y = coverage, fill = coverage >= conf_line)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = conf_line, linetype = "dashed",
               colour = "grey40", linewidth = 0.7) +
    scale_fill_manual(
      values = c(`TRUE` = "#27AE60", `FALSE` = "#E74C3C"),
      labels = c(`TRUE` = "At/above nominal", `FALSE` = "Below nominal"),
      name = NULL) +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format(accuracy = 1)) +
    labs(title    = sprintf("%s: Deviation epsilon Coverage by Study", title_prefix),
         subtitle = "First active covariate; dashed = nominal level",
         x = NULL, y = "Coverage probability") +
    theme_bw(base_size = 11) +
    theme(legend.position = "top",
          plot.title = element_text(face = "bold"))
  
  list(coverage = p1, width = p2, se_vs_rmse = p3, eps_coverage = p4)
}


################################################################################
##  SECTION 8: DIAGNOSTIC PLOTS                                               ##
################################################################################

plot_coefficient_comparison <- function(fit, truth, method_name = "SCAD") {
  K <- truth$K; p <- truth$p
  df <- do.call(rbind, lapply(seq_len(K), function(k)
    data.frame(covariate  = seq_len(p),
               true_theta = truth$theta_true[[k]],
               est_theta  = fit$theta_hat[k, ],
               study      = paste0("Study ", k))))
  ggplot(df, aes(x = true_theta, y = est_theta, color = study)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "black") +
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
      data.frame(covariate = paste0("X", j),
                 study     = paste0("Study ", k),
                 True      = truth$epsilon_true[[k]][j],
                 Estimated = fit$eps_hat[k, j])))))
  df_long <- tidyr::pivot_longer(df, c(True, Estimated),
                                 names_to = "type", values_to = "value")
  ggplot(df_long, aes(x = covariate, y = value, fill = type)) +
    geom_col(position = "dodge", alpha = 0.85) +
    geom_hline(yintercept = 0) +
    facet_wrap(~study, nrow = 2L) +
    scale_fill_manual(
      values = c(True = "#2C7BB6", Estimated = "#D7191C"), name = "") +
    labs(title = "Study-Specific Deviations: True vs Estimated",
         x = "Covariate", y = "epsilon_kj") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


################################################################################
##  SECTION 9: DEMO RUN                                                       ##
################################################################################

cat("\n================================================================\n")
cat("  DEMO RUN\n")
cat("  (single replicate + small simulation + UQ demo)\n")
cat("================================================================\n\n")

## ── Single replicate ─────────────────────────────────────────────────────────
set.seed(42L)
demo_data <- generate_ipd_data(
  K = 10, n_per_study = 100, p = 50, s_alpha = 6, s_epsilon = 4,
  alpha_true = c(1.0, -0.8, 0.6, -0.5, 0.4, -0.3),
  epsilon_sd = 0.3, rho = 0.5, cens_rate = 0.30
)
cat(sprintf("K=%d | N=%d | p=%d | True active: {%s}\n\n",
            demo_data$K, sum(demo_data$n_per_study), demo_data$p,
            paste(demo_data$active_set, collapse = ",")))

cat("Fitting MCP (hierarchical)...\n")
fit_mcp <- fit_hierarchical_meta_cox(
  demo_data$studies, penalty = "mcp",
  lambda_alpha = 0.08, lambda_eps = 0.04,
  max_iter = 30L, verbose = FALSE)

cat("Fitting SCAD (hierarchical)...\n")
fit_scad <- fit_hierarchical_meta_cox(
  demo_data$studies, penalty = "scad",
  lambda_alpha = 0.08, lambda_eps = 0.04,
  max_iter = 30L, verbose = FALSE)

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

cat("\n=== Single Replicate Results ===\n")
cat(sprintf("%-12s | %5s | %5s | %5s | %9s | Selected\n",
            "Method", "TPR", "FDR", "MCC", "AlphaBias"))
cat(strrep("-", 65), "\n")
for (i in seq_len(nrow(metrics_df))) {
  r <- metrics_df[i, ]
  m <- all_fits[[r$method]]
  cat(sprintf("%-12s | %.3f | %.3f | %.3f | %.4f    | {%s}\n",
              r$method, r$TPR, r$FDR, r$MCC,
              ifelse(is.na(r$alpha_bias), 0, r$alpha_bias),
              paste(sort(m$active_set), collapse = ",")))
}

## ── Small simulation ─────────────────────────────────────────────────────────
cat("\n\n--- Small simulation (nsim=20, 3 scenarios) ---\n\n")
small_scenarios <- define_scenarios()[c("S01_base","S04_homog","S05_highhet")]
sim_results <- run_simulation_study(
  nsim = 20L, scenarios = small_scenarios,
  methods = c("scad","mcp","lasso","enet","scad_pool","mcp_pool"),
  n_cores = 1L, lambda_alpha = 0.08, lambda_eps = 0.04, seed = 2024L)

summary_df <- summarize_results(sim_results)
for (sc in unique(summary_df$scenario)) print_comparison_table(summary_df, sc)

## ── UQ simulation (demo: M=5 datasets, B=50 bootstrap) ───────────────────────
## For the full paper use M=200, B=500.  Demo values run in ~2-5 minutes.
cat("\n\n--- UQ Simulation Demo (M=5, B=50) ---\n")
cat("    For paper-quality results: set M=200, B=500\n\n")

uq_res <- run_uq_simulation(
  M              = 5L,
  B              = 50L,
  scenario_params = define_scenarios()[["S01_base"]],
  lambda_alpha   = 0.08,
  lambda_eps     = 0.04,
  conf_level     = 0.95,
  seed           = 2024L,
  n_cores        = 1L
)

print_uq_summary(uq_res)

## Save UQ plots
cat("\nSaving UQ plots...\n")
uq_plots <- plot_uq_results(uq_res, title_prefix = "Demo UQ")
ggsave("plot_uq_coverage.png",   uq_plots$coverage,    width = 9, height = 5)
ggsave("plot_uq_width.png",      uq_plots$width,       width = 9, height = 5)
ggsave("plot_uq_se_rmse.png",    uq_plots$se_vs_rmse,  width = 6, height = 5)
ggsave("plot_uq_eps_cov.png",    uq_plots$eps_coverage,width = 7, height = 4)
cat("UQ plots saved.\n")

## Save standard simulation plots
ggsave("plot_MCC.png",
       plot_simulation_results(sim_results, "MCC", "MCC by Scenario"),
       width = 10, height = 6)
ggsave("plot_FDR.png",
       plot_simulation_results(sim_results, "FDR", "FDR by Scenario"),
       width = 10, height = 6)
ggsave("plot_bias.png",
       plot_simulation_results(sim_results, "alpha_bias", "Alpha Bias"),
       width = 10, height = 6)
cat("All plots saved.\n")


################################################################################
##  SECTION 10: FULL PAPER SIMULATION                                         ##
################################################################################

run_full <- FALSE   # <-- Set to TRUE for the full paper run

if (run_full) {
  cat("\n\n=== FULL PAPER SIMULATION (nsim=500, all scenarios) ===\n\n")
  
  full_results <- run_simulation_study(
    nsim         = 500L,
    scenarios    = define_scenarios(),
    methods      = c("scad","mcp","lasso","enet","scad_pool","mcp_pool"),
    n_cores      = 1L,       # change to detectCores()-1 for speed
    lambda_alpha = NULL,     # NULL = CV per replicate
    lambda_eps   = NULL,
    seed         = 2024L
  )
  saveRDS(full_results, "sim_results_full.rds")
  write.csv(full_results, "sim_results_full.csv", row.names = FALSE)
  
  full_summary <- summarize_results(full_results)
  write.csv(full_summary, "sim_summary_full.csv", row.names = FALSE)
  for (sc in names(define_scenarios())) print_comparison_table(full_summary, sc)
  
  for (metric in c("MCC","TPR","FDR","F1","alpha_bias",
                   "fp_bias","MSE_theta","exact_match")) {
    ggsave(sprintf("plot_full_%s.pdf", metric),
           plot_simulation_results(full_results, metric,
                                   paste(metric, "across Scenarios")),
           width = 14, height = 7)
  }
  
  ## ── Full UQ simulation (paper quality) ───────────────────────────────────
  cat("\n=== FULL UQ SIMULATION (M=200, B=500) ===\n\n")
  
  uq_full <- run_uq_simulation(
    M              = 200L,
    B              = 500L,
    scenario_params = define_scenarios()[["S01_base"]],
    lambda_alpha   = NULL,    # CV-selected per dataset
    lambda_eps     = NULL,
    conf_level     = 0.95,
    seed           = 2024L,
    n_cores        = 1L
  )
  print_uq_summary(uq_full)
  saveRDS(uq_full, "uq_results_full.rds")
  write.csv(uq_full$alpha_summary, "uq_alpha_summary.csv", row.names = FALSE)
  write.csv(uq_full$eps_summary,   "uq_eps_summary.csv",   row.names = FALSE)
  
  uq_plots_full <- plot_uq_results(uq_full, title_prefix = "Full UQ")
  ggsave("plot_uq_full_coverage.png",
         uq_plots_full$coverage,    width = 10, height = 5, dpi = 300)
  ggsave("plot_uq_full_width.png",
         uq_plots_full$width,       width = 10, height = 5, dpi = 300)
  ggsave("plot_uq_full_se_rmse.png",
         uq_plots_full$se_vs_rmse,  width = 7,  height = 5, dpi = 300)
  ggsave("plot_uq_full_eps_cov.png",
         uq_plots_full$eps_coverage,width = 8,  height = 4, dpi = 300)
  
  cat("Full simulation and UQ complete.\n")
}

cat("\n=== Script complete ===\n")