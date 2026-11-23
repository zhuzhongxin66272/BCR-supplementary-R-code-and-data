# ================================================================
# Part1
# BCR_global_simulation_pilot_v2_full
# BCR / WQS / QGCOMP simulation framework for GLOBAL TESTING only
# Full version with incremental saving, resume support, and summaries
# ================================================================

required_packages <- c("MASS", "fda", "glmnet", "gWQS", "qgcomp", "dplyr")

load_required_packages <- function() {
  missing_pkgs <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing packages: ", paste(missing_pkgs, collapse = ", "),
      "\nPlease install them before running."
    )
  }
  invisible(lapply(required_packages, library, character.only = TRUE))
}

# ================================================================
# Scenario generation
# ================================================================

make_pilot_scenarios <- function() {
  expand.grid(
    scenario_type = c("null", "nonnull"),
    n = c(150, 300),
    M = c(50, 100),
    corr_type = c("cs", "block"),
    rho = c(0.5, 0.8),
    signal = c(0.01, 0.03),
    structure = c("sparse_same_sign", "sparse_mixed"),
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(
      signal = ifelse(scenario_type == "null", 0, signal),
      structure = ifelse(scenario_type == "null", "null", structure),
      scenario_id = dplyr::row_number()
    ) |>
    dplyr::distinct()
}

make_formal_scenarios_main <- function() {
  n_vals <- c(300, 800, 1500)
  M_vals <- c(50, 100, 200)
  
  dep_settings <- data.frame(
    corr_type = c("cs", "ar1", "block"),
    rho       = c(0.1,   0.3,   0.6),
    stringsAsFactors = FALSE
  )
  
  base_grid <- merge(
    expand.grid(
      n = n_vals,
      M = M_vals,
      stringsAsFactors = FALSE
    ),
    dep_settings,
    by = NULL
  )
  
  null_df <- base_grid
  null_df$scenario_type <- "null"
  null_df$signal <- 0.00
  null_df$structure <- "null"
  
  same_df <- merge(
    base_grid,
    data.frame(signal = c(0.05, 0.10), stringsAsFactors = FALSE),
    by = NULL
  )
  same_df$scenario_type <- "nonnull"
  same_df$structure <- "sparse_same_sign"
  
  mixed10_df <- base_grid
  mixed10_df$scenario_type <- "nonnull"
  mixed10_df$signal <- 0.10
  mixed10_df$structure <- "sparse_mixed"
  
  mixed05_df <- subset(base_grid, n >= 800 | M <= 100)
  mixed05_df$scenario_type <- "nonnull"
  mixed05_df$signal <- 0.05
  mixed05_df$structure <- "sparse_mixed"
  
  scenarios <- rbind(
    null_df,
    same_df,
    mixed05_df,
    mixed10_df
  )
  
  scenarios <- scenarios[order(
    scenarios$scenario_type,
    scenarios$structure,
    scenarios$signal,
    scenarios$n,
    scenarios$M,
    scenarios$corr_type,
    scenarios$rho
  ), ]
  
  rownames(scenarios) <- NULL
  scenarios$scenario_id <- seq_len(nrow(scenarios))
  
  scenarios[, c(
    "scenario_id", "scenario_type", "n", "M",
    "corr_type", "rho", "signal", "structure"
  )]
}

make_formal_scenarios_refine <- function(main_scenarios = make_formal_scenarios_main()) {
  ## Second-round refinement set:
  ## 1) all null scenarios (for more precise type I error)
  ## 2) all sparse_mixed, signal = 0.10 (key advantage scenarios)
  ## 3) sparse_mixed, signal = 0.05 with n >= 800 (moderately informative)
  ## 4) sparse_same_sign, signal = 0.10 with n >= 800 (stable comparison against WQS)
  
  refine <- subset(
    main_scenarios,
    structure == "null" |
      (structure == "sparse_mixed" & signal == 0.10) |
      (structure == "sparse_mixed" & signal == 0.05 & n >= 800) |
      (structure == "sparse_same_sign" & signal == 0.10 & n >= 800)
  )
  
  refine[order(refine$scenario_id), ]
}

# ================================================================
# Data generation
# ================================================================

generate_sigma <- function(M, corr_type = "cs", rho = 0.5, n_blocks = 5) {
  if (corr_type == "cs") {
    S <- matrix(rho, M, M)
    diag(S) <- 1
    return(S)
  }
  if (corr_type == "ar1") {
    idx <- seq_len(M)
    S <- outer(idx, idx, function(i, j) rho^(abs(i - j)))
    return(S)
  }
  if (corr_type == "block") {
    S <- matrix(0.1, M, M)
    diag(S) <- 1
    block_id <- rep(seq_len(n_blocks), length.out = M)
    for (b in unique(block_id)) {
      ii <- which(block_id == b)
      S[ii, ii] <- rho
      diag(S)[ii] <- 1
    }
    return(S)
  }
  stop("Unknown corr_type")
}

generate_theta <- function(M, structure = "sparse_same_sign") {
  theta <- rep(0, M)
  if (structure == "null") return(theta)
  
  n_active <- max(3, ceiling(0.1 * M))
  active <- sort(sample(seq_len(M), n_active, replace = FALSE))
  
  if (structure == "sparse_same_sign") {
    vals <- seq(0.2, 1, length.out = n_active)
  } else if (structure == "sparse_mixed") {
    vals <- seq(0.2, 1, length.out = n_active)
    signs <- rep(c(-1, 1), length.out = n_active)
    vals <- vals * signs
  } else if (structure == "clustered_mixed") {
    active <- unique(c(1:ceiling(n_active / 2), (M - ceiling(n_active / 2) + 1):M))
    vals <- c(
      seq(0.3, 1, length.out = floor(length(active) / 2)),
      -seq(0.3, 1, length.out = ceiling(length(active) / 2))
    )
  } else {
    stop("Unknown structure")
  }
  
  theta[active] <- vals
  theta
}

scale_theta_to_target <- function(theta, Sigma, gamma, Zmat, target_R2) {
  if (all(theta == 0) || target_R2 <= 0) {
    return(list(theta = theta, sigma2 = 1))
  }
  
  var_Z <- as.numeric(var(Zmat %*% gamma))
  qform <- as.numeric(t(theta) %*% Sigma %*% theta)
  if (qform <= 0) {
    return(list(theta = theta * 0, sigma2 = 1))
  }
  
  c2 <- (target_R2 * (var_Z + 1)) / ((1 - target_R2) * qform)
  c_scale <- sqrt(max(c2, 0))
  theta_scaled <- c_scale * theta
  var_theta <- as.numeric(t(theta_scaled) %*% Sigma %*% theta_scaled)
  sigma2 <- max(0.05, (var_theta + var_Z) * (1 - target_R2) / target_R2 - var_Z - var_theta)
  
  list(theta = theta_scaled, sigma2 = sigma2)
}

generate_dataset <- function(n, M, corr_type, rho, structure, signal,
                             gamma = c(0.5, 1.0, 0.3)) {
  Z1 <- rbinom(n, 1, 0.5)
  Z2 <- runif(n)
  Z3 <- rnorm(n)
  Z <- cbind(Z1, Z2, Z3)
  colnames(Z) <- paste0("Z", 1:3)
  
  Sigma <- generate_sigma(M, corr_type, rho)
  theta0 <- generate_theta(M, structure)
  scaled <- scale_theta_to_target(theta0, Sigma, gamma, Z, signal)
  theta <- scaled$theta
  sigma2 <- scaled$sigma2
  
  X <- MASS::mvrnorm(n, mu = rep(0, M), Sigma = Sigma)
  X <- scale(X)
  colnames(X) <- paste0("X", 1:M)
  
  Y <- as.numeric(X %*% theta + Z %*% gamma + rnorm(n, 0, sqrt(sigma2)))
  data <- data.frame(Y = Y, X, Z, check.names = FALSE)
  
  list(
    data = data,
    theta_true = theta,
    sigma2 = sigma2,
    Sigma = Sigma,
    theta_total_true = sum(theta * rep(1, M))
  )
}

# ================================================================
# BCR
# ================================================================

bcr_fit <- function(data, M, K = NULL) {
  if (is.null(K)) K <- max(4, ceiling(sqrt(M)))
  X <- as.matrix(data[, paste0("X", seq_len(M))])
  Z <- as.matrix(data[, paste0("Z", 1:3)])
  Y <- data$Y
  
  cvfit <- glmnet::cv.glmnet(cbind(X, Z), Y, alpha = 0, nfolds = 5)
  cf <- as.matrix(glmnet::coef.glmnet(cvfit, s = "lambda.min"))[-1, 1]
  ridge_beta <- cf[seq_len(M)]
  ordering <- order(ridge_beta)
  
  pos <- seq_len(M) / M
  basis_obj <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = K, norder = 4)
  B <- fda::eval.basis(pos, basis_obj)
  X_ord <- X[, ordering, drop = FALSE]
  U <- X_ord %*% B
  colnames(U) <- paste0("U", seq_len(K))
  
  dat_full <- data.frame(Y = Y, U, Z)
  dat_red <- data.frame(Y = Y, Z)
  fit_full <- lm(Y ~ ., data = dat_full)
  fit_red <- lm(Y ~ ., data = dat_red)
  
  rss_full <- sum(resid(fit_full)^2)
  rss_red <- sum(resid(fit_red)^2)
  df1 <- K
  df2 <- nrow(data) - K - ncol(Z) - 1
  F_stat <- ((rss_red - rss_full) / df1) / (rss_full / df2)
  p_naive <- pf(F_stat, df1 = df1, df2 = df2, lower.tail = FALSE)
  
  beta_hat <- coef(fit_full)[2:(K + 1)]
  theta_ord <- as.numeric(B %*% beta_hat)
  theta_hat <- numeric(M)
  theta_hat[ordering] <- theta_ord
  
  list(
    F_stat = F_stat,
    p_naive = p_naive,
    fit_full = fit_full,
    fit_red = fit_red,
    theta_hat = theta_hat,
    ordering = ordering,
    B = B,
    K = K
  )
}

bcr_global_test <- function(data, M, K = NULL, B = 199,
                            resample = c("permutation", "bootstrap")) {
  resample <- match.arg(resample)
  obs <- bcr_fit(data, M, K)
  F_obs <- obs$F_stat
  
  Y <- data$Y
  Z <- as.matrix(data[, paste0("Z", 1:3)])
  fit_red <- lm(Y ~ Z)
  yhat <- fitted(fit_red)
  ehat <- resid(fit_red)
  
  F_null <- numeric(B)
  for (b in seq_len(B)) {
    if (resample == "permutation") {
      e_star <- sample(ehat, replace = FALSE)
    } else {
      e_star <- sample(ehat, replace = TRUE)
    }
    data_b <- data
    data_b$Y <- yhat + e_star
    F_null[b] <- bcr_fit(data_b, M, K)$F_stat
  }
  
  p_resample <- (sum(F_null >= F_obs) + 1) / (B + 1)
  
  list(
    p_value = p_resample,
    reject = p_resample < 0.05,
    F_obs = F_obs,
    p_naive = obs$p_naive,
    fit = obs,
    B = B
  )
}

bcr_total_effect_ci <- function(bcr_obj, contrast = NULL, alpha = 0.05) {
  theta_hat <- bcr_obj$theta_hat
  M <- length(theta_hat)
  if (is.null(contrast)) contrast <- rep(1, M)
  
  sigma2 <- sum(resid(bcr_obj$fit_full)^2) / df.residual(bcr_obj$fit_full)
  Xd <- model.matrix(bcr_obj$fit_full)
  U_cols <- grep("^U", colnames(Xd))
  U <- Xd[, U_cols, drop = FALSE]
  Var_beta <- sigma2 * solve(t(U) %*% U)
  Var_theta_ord <- bcr_obj$B %*% Var_beta %*% t(bcr_obj$B)
  
  ord <- bcr_obj$ordering
  inv_ord <- order(ord)
  Var_theta <- Var_theta_ord[inv_ord, inv_ord, drop = FALSE]
  
  est <- sum(theta_hat * contrast)
  se <- sqrt(as.numeric(t(contrast) %*% Var_theta %*% contrast))
  z <- qnorm(1 - alpha / 2)
  
  c(est = est, se = se, lower = est - z * se, upper = est + z * se)
}

# ================================================================
# WQS and QGCOMP: global test only
# ================================================================

wqs_global_test <- function(data, M, wqs_b = 30) {
  X_vars <- paste0("X", seq_len(M))
  Z_vars <- paste0("Z", 1:3)
  fml <- as.formula(paste("Y ~ wqs +", paste(Z_vars, collapse = " + ")))
  
  fit <- gWQS::gwqs(
    formula = fml,
    mix_name = X_vars,
    data = data,
    q = 4,
    validation = 0.6,
    b = wqs_b,
    b_constr = FALSE,
    family = "gaussian"
  )
  
  cf <- coef(fit)
  idx <- grep("^wqs", names(cf), ignore.case = TRUE)
  if (length(idx) == 0) stop("WQS coefficient not found")
  vc <- vcov(fit)
  beta <- as.numeric(cf[idx[1]])
  se <- sqrt(vc[idx[1], idx[1]])
  z <- beta / se
  p <- 2 * (1 - pnorm(abs(z)))
  
  list(p_value = p, reject = p < 0.05)
}

qgc_global_test <- function(data, M, qgc_B = 50) {
  X_vars <- paste0("X", seq_len(M))
  Z_vars <- paste0("Z", 1:3)
  fml <- as.formula(paste("Y ~", paste(c(X_vars, Z_vars), collapse = " + ")))
  
  fit <- qgcomp::qgcomp.boot(
    fml,
    expnms = X_vars,
    data = data,
    q = 4,
    B = qgc_B,
    family = gaussian(),
    bayes = FALSE
  )
  
  p <- if (!is.null(fit$pval)) fit$pval else {
    psi <- fit$psi
    se <- sqrt(fit$var.psi)
    2 * (1 - pnorm(abs(psi / se)))
  }
  
  list(p_value = as.numeric(p), reject = as.numeric(p) < 0.05)
}

# ================================================================
# One replicate with explicit error logging
# ================================================================

run_one_replicate <- function(scenario, rep_id,
                              bcr_B = 199,
                              bcr_resample = "permutation",
                              wqs_b = 30,
                              qgc_B = 50) {
  n <- scenario$n
  M <- scenario$M
  corr_type <- scenario$corr_type
  rho <- scenario$rho
  signal <- scenario$signal
  structure <- scenario$structure
  scenario_type <- scenario$scenario_type
  
  dat <- generate_dataset(n, M, corr_type, rho, structure, signal)
  data <- dat$data
  
  out <- list(
    scenario_id = scenario$scenario_id,
    scenario_type = scenario_type,
    n = n,
    M = M,
    corr_type = corr_type,
    rho = rho,
    signal = signal,
    structure = structure,
    replicate = rep_id,
    theta_total_true = dat$theta_total_true,
    bcr_p = NA_real_,
    bcr_reject = NA,
    bcr_time = NA_real_,
    bcr_error = NA_character_,
    bcr_theta_total_hat = NA_real_,
    bcr_theta_total_lower = NA_real_,
    bcr_theta_total_upper = NA_real_,
    bcr_theta_rmse = NA_real_,
    wqs_p = NA_real_,
    wqs_reject = NA,
    wqs_time = NA_real_,
    wqs_error = NA_character_,
    qgc_p = NA_real_,
    qgc_reject = NA,
    qgc_time = NA_real_,
    qgc_error = NA_character_,
    fatal_error = NA_character_
  )
  
  # BCR
  t0 <- proc.time()[3]
  tryCatch({
    bcr_res <- bcr_global_test(data, M, B = bcr_B, resample = bcr_resample)
    ci <- bcr_total_effect_ci(bcr_res$fit, contrast = rep(1, M))
    out$bcr_p <- bcr_res$p_value
    out$bcr_reject <- bcr_res$reject
    out$bcr_theta_total_hat <- ci["est"]
    out$bcr_theta_total_lower <- ci["lower"]
    out$bcr_theta_total_upper <- ci["upper"]
    out$bcr_theta_rmse <- sqrt(mean((bcr_res$fit$theta_hat - dat$theta_true)^2))
  }, error = function(e) {
    out$bcr_error <<- conditionMessage(e)
  })
  out$bcr_time <- proc.time()[3] - t0
  
  # WQS
  t0 <- proc.time()[3]
  tryCatch({
    wqs_res <- wqs_global_test(data, M, wqs_b = wqs_b)
    out$wqs_p <- wqs_res$p_value
    out$wqs_reject <- wqs_res$reject
  }, error = function(e) {
    out$wqs_error <<- conditionMessage(e)
  })
  out$wqs_time <- proc.time()[3] - t0
  
  # QGCOMP
  t0 <- proc.time()[3]
  tryCatch({
    qgc_res <- qgc_global_test(data, M, qgc_B = qgc_B)
    out$qgc_p <- qgc_res$p_value
    out$qgc_reject <- qgc_res$reject
  }, error = function(e) {
    out$qgc_error <<- conditionMessage(e)
  })
  out$qgc_time <- proc.time()[3] - t0
  
  as.data.frame(out, stringsAsFactors = FALSE)
}

# ================================================================
# Summaries
# ================================================================

summarise_type1_error <- function(results) {
  results |>
    dplyr::filter(scenario_type == "null") |>
    dplyr::group_by(scenario_id, n, M, corr_type, rho, signal, structure) |>
    dplyr::summarise(
      n_reps_completed = dplyr::n(),
      bcr_type1 = mean(bcr_reject, na.rm = TRUE),
      wqs_type1 = mean(wqs_reject, na.rm = TRUE),
      qgc_type1 = mean(qgc_reject, na.rm = TRUE),
      bcr_time_mean = mean(bcr_time, na.rm = TRUE),
      wqs_time_mean = mean(wqs_time, na.rm = TRUE),
      qgc_time_mean = mean(qgc_time, na.rm = TRUE),
      bcr_nonmissing = sum(!is.na(bcr_reject)),
      wqs_nonmissing = sum(!is.na(wqs_reject)),
      qgc_nonmissing = sum(!is.na(qgc_reject)),
      .groups = "drop"
    )
}

summarise_power <- function(results) {
  results |>
    dplyr::filter(scenario_type == "nonnull") |>
    dplyr::group_by(scenario_id, n, M, corr_type, rho, signal, structure) |>
    dplyr::summarise(
      n_reps_completed = dplyr::n(),
      bcr_power = mean(bcr_reject, na.rm = TRUE),
      wqs_power = mean(wqs_reject, na.rm = TRUE),
      qgc_power = mean(qgc_reject, na.rm = TRUE),
      bcr_time_mean = mean(bcr_time, na.rm = TRUE),
      wqs_time_mean = mean(wqs_time, na.rm = TRUE),
      qgc_time_mean = mean(qgc_time, na.rm = TRUE),
      bcr_theta_total_bias_mean = mean(bcr_theta_total_hat - theta_total_true, na.rm = TRUE),
      bcr_theta_total_cover = mean(theta_total_true >= bcr_theta_total_lower & theta_total_true <= bcr_theta_total_upper, na.rm = TRUE),
      bcr_theta_rmse_mean = mean(bcr_theta_rmse, na.rm = TRUE),
      bcr_nonmissing = sum(!is.na(bcr_reject)),
      wqs_nonmissing = sum(!is.na(wqs_reject)),
      qgc_nonmissing = sum(!is.na(qgc_reject)),
      .groups = "drop"
    )
}

summarise_one_error <- function(x, method_name) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) {
    return(data.frame(
      method = method_name,
      error = NA_character_,
      n = 0,
      stringsAsFactors = FALSE
    ))
  }
  tab <- sort(table(x), decreasing = TRUE)
  data.frame(
    method = method_name,
    error = names(tab),
    n = as.integer(tab),
    stringsAsFactors = FALSE
  )
}

summarise_errors <- function(results) {
  out1 <- summarise_one_error(results$bcr_error, "BCR")
  out2 <- summarise_one_error(results$wqs_error, "WQS")
  out3 <- summarise_one_error(results$qgc_error, "QGC")
  dplyr::bind_rows(out1, out2, out3)
}
summarise_runtime <- function(results) {
  results |>
    dplyr::group_by(scenario_id, scenario_type, n, M, corr_type, rho, signal, structure) |>
    dplyr::summarise(
      n_reps_completed = dplyr::n(),
      bcr_time_mean = mean(bcr_time, na.rm = TRUE),
      wqs_time_mean = mean(wqs_time, na.rm = TRUE),
      qgc_time_mean = mean(qgc_time, na.rm = TRUE),
      .groups = "drop"
    )
}

combine_results_from_outdir <- function(outdir) {
  rds_files <- list.files(outdir, pattern = "^scenario_[0-9]{3}\\.rds$", full.names = TRUE)
  if (length(rds_files) == 0) return(data.frame())
  dplyr::bind_rows(lapply(rds_files, readRDS))
}

update_running_summaries <- function(outdir) {
  all_res <- combine_results_from_outdir(outdir)
  if (nrow(all_res) == 0) return(invisible(NULL))
  
  saveRDS(all_res, file.path(outdir, "all_results_current.rds"))
  write.csv(all_res, file.path(outdir, "all_results_current.csv"), row.names = FALSE)
  
  err_tab <- summarise_errors(all_res)
  type1_tab <- summarise_type1_error(all_res)
  power_tab <- summarise_power(all_res)
  runtime_tab <- summarise_runtime(all_res)
  
  write.csv(err_tab, file.path(outdir, "summary_errors_current.csv"), row.names = FALSE)
  write.csv(type1_tab, file.path(outdir, "summary_type1_current.csv"), row.names = FALSE)
  write.csv(power_tab, file.path(outdir, "summary_power_current.csv"), row.names = FALSE)
  write.csv(runtime_tab, file.path(outdir, "summary_runtime_current.csv"), row.names = FALSE)
  
  invisible(NULL)
}

# ================================================================
# Incremental runner with resume support
# ================================================================

run_simulation_incremental <- function(
    scenarios,
    n_reps = 300,
    outdir = "simulation_incremental",
    bcr_B = 199,
    bcr_resample = "permutation",
    wqs_b = 30,
    qgc_B = 50,
    save_every = 5,
    overwrite = FALSE,
    verbose = TRUE
) {
  load_required_packages()
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  write.csv(scenarios, file.path(outdir, "scenario_grid.csv"), row.names = FALSE)
  saveRDS(scenarios, file.path(outdir, "scenario_grid.rds"))
  
  n_scen <- nrow(scenarios)
  
  cat("====================================================\n")
  cat("Incremental simulation started\n")
  cat("Output directory:", outdir, "\n")
  cat("Number of scenarios:", n_scen, "\n")
  cat("Target replicates per scenario:", n_reps, "\n")
  cat("save_every:", save_every, "\n")
  cat("====================================================\n\n")
  
  for (i in seq_len(n_scen)) {
    sc <- as.list(scenarios[i, , drop = FALSE])
    sid <- sprintf("%03d", sc$scenario_id)
    
    scenario_rds <- file.path(outdir, paste0("scenario_", sid, ".rds"))
    scenario_csv <- file.path(outdir, paste0("scenario_", sid, ".csv"))
    
    if (verbose) {
      cat("----------------------------------------------------\n")
      cat(sprintf(
        "[%d/%d] scenario_id=%s | type=%s | n=%d | M=%d | %s rho=%.1f | signal=%.2f | %s\n",
        i, n_scen, sid, sc$scenario_type, sc$n, sc$M, sc$corr_type, sc$rho, sc$signal, sc$structure
      ))
    }
    
    if (file.exists(scenario_rds) && !overwrite) {
      res_s <- readRDS(scenario_rds)
      done_reps <- if ("replicate" %in% names(res_s)) sort(unique(res_s$replicate)) else integer(0)
      start_rep <- if (length(done_reps) == 0) 1 else max(done_reps) + 1
      if (verbose) cat("Existing file found. Completed reps:", length(done_reps), "\n")
    } else {
      res_s <- data.frame()
      start_rep <- 1
      if (verbose) cat("Starting new scenario.\n")
    }
    
    if (start_rep > n_reps) {
      if (verbose) cat("Scenario already completed for target n_reps. Skipping.\n")
      next
    }
    
    for (r in start_rep:n_reps) {
      
      bcr_B_use <- if (sc$scenario_type == "null") 499 else bcr_B
      
      one_res <- tryCatch(
        {
          run_one_replicate(
            scenario = sc,
            rep_id = r,
            bcr_B = bcr_B_use,
            bcr_resample = bcr_resample,
            wqs_b = wqs_b,
            qgc_B = qgc_B
          )
        },
        error = function(e) {
          data.frame(
            scenario_id = sc$scenario_id,
            scenario_type = sc$scenario_type,
            n = sc$n,
            M = sc$M,
            corr_type = sc$corr_type,
            rho = sc$rho,
            signal = sc$signal,
            structure = sc$structure,
            replicate = r,
            theta_total_true = NA_real_,
            bcr_p = NA_real_,
            bcr_reject = NA,
            bcr_time = NA_real_,
            bcr_error = NA_character_,
            bcr_theta_total_hat = NA_real_,
            bcr_theta_total_lower = NA_real_,
            bcr_theta_total_upper = NA_real_,
            bcr_theta_rmse = NA_real_,
            wqs_p = NA_real_,
            wqs_reject = NA,
            wqs_time = NA_real_,
            wqs_error = NA_character_,
            qgc_p = NA_real_,
            qgc_reject = NA,
            qgc_time = NA_real_,
            qgc_error = NA_character_,
            fatal_error = as.character(e$message),
            stringsAsFactors = FALSE
          )
        }
      )
      
      res_s <- dplyr::bind_rows(res_s, one_res)
      
      if ((r %% save_every == 0) || (r == n_reps)) {
        saveRDS(res_s, scenario_rds)
        write.csv(res_s, scenario_csv, row.names = FALSE)
        if (verbose) cat(sprintf("  saved after rep %d / %d\n", r, n_reps))
      }
    }
    
    update_running_summaries(outdir)
    if (verbose) cat("Scenario finished and master summaries updated.\n")
  }
  
  cat("\n====================================================\n")
  cat("Incremental simulation completed.\n")
  cat("Results saved in:", outdir, "\n")
  cat("====================================================\n")
  
  invisible(outdir)
}


##Part2#########################################################################
##--generate_main_scenarios--##
make_main_scenarios <- function() {
  
  n_vals <- c(300, 800, 1500)
  M_vals <- c(50, 100)
  corr_types <- c("cs", "ar1", "block")
  rho_vals <- c(0.1, 0.3, 0.6)
  signal_vals <- c(0.05, 0.10, 0.20)
  
  base_grid <- expand.grid(
    n = n_vals,
    M = M_vals,
    corr_type = corr_types,
    rho = rho_vals,
    stringsAsFactors = FALSE
  )
  
  ## -------------------------
  ## Null scenarios
  ## -------------------------
  null_df <- base_grid
  null_df$scenario_type <- "null"
  null_df$signal <- 0.00
  null_df$structure <- "null"
  
  ## -------------------------
  ## Non-null: sparse same-sign
  ## -------------------------
  same_df <- merge(
    base_grid,
    data.frame(signal = signal_vals, stringsAsFactors = FALSE),
    by = NULL
  )
  same_df$scenario_type <- "nonnull"
  same_df$structure <- "sparse_same_sign"
  
  ## -------------------------
  ## Non-null: sparse mixed-sign
  ## -------------------------
  mixed_df <- merge(
    base_grid,
    data.frame(signal = signal_vals, stringsAsFactors = FALSE),
    by = NULL
  )
  mixed_df$scenario_type <- "nonnull"
  mixed_df$structure <- "sparse_mixed"
  
  ## -------------------------
  ## Combine all
  ## -------------------------
  scenarios <- rbind(null_df, same_df, mixed_df)
  
  scenarios <- scenarios[order(
    scenarios$scenario_type,
    scenarios$structure,
    scenarios$signal,
    scenarios$n,
    scenarios$M,
    scenarios$corr_type,
    scenarios$rho
  ), ]
  
  rownames(scenarios) <- NULL
  scenarios$scenario_id <- seq_len(nrow(scenarios))
  
  scenarios <- scenarios[, c(
    "scenario_id", "scenario_type", "n", "M",
    "corr_type", "rho", "signal", "structure"
  )]
  
  cat("=====================================\n")
  cat("Main analysis scenarios generated.\n")
  cat("Total scenarios:", nrow(scenarios), "\n")
  cat("Expected total = 54 null + 162 same-sign + 162 mixed-sign = 378\n")
  cat("=====================================\n")
  print(table(scenarios$structure, scenarios$signal))
  
  return(scenarios)
}


split_main_scenarios <- function(scenarios) {
  list(
    null  = subset(scenarios, scenario_type == "null"),
    same  = subset(scenarios, structure == "sparse_same_sign"),
    mixed = subset(scenarios, structure == "sparse_mixed")
  )
}


#Part3#########################################################################
# ============================================================
# run_simulation_parallel_universal.R
# Parallel simulation runner for null / same / mixed scenarios
# ============================================================

setwd("C:/Users/zhang/Desktop/BCR_parallel")

# ---- source main functions ----
source("BCR_global_simulation_pilot_v2_full3.R")
source("generate_main_scenarios3.R")

# ---- packages (install if missing) ----
# install.packages(c("cowplot", "gWQS", "foreach", "doParallel", "doRNG"))
library(cowplot)
library(gWQS)
library(foreach)
library(doParallel)
library(doRNG)

# ---- force single-thread BLAS to avoid oversubscription ----
Sys.setenv(OMP_NUM_THREADS   = 1)
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS   = 1)

# ---- main parallel function ----
run_simulation_parallel <- function(
    scenarios,
    n_reps,
    outdir,
    bcr_B = 199,
    wqs_b = 30,
    qgc_B = 50,
    n_cores = 8,
    batch_size = 100
) {
  cl <- makePSOCKcluster(n_cores, useXDR = FALSE)
  registerDoParallel(cl)
  
  clusterEvalQ(cl, {
    library(MASS)
    library(glmnet)
    library(gWQS)
    library(qgcomp)
    library(dplyr)
    source("C:/Users/zhang/Desktop/BCR_parallel/BCR_global_simulation_pilot_v2_full3.R")
  })
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # save scenario table for reference
  write.csv(scenarios, file.path(outdir, "scenario_grid.csv"), row.names = FALSE)
  saveRDS(scenarios, file.path(outdir, "scenario_grid.rds"))
  
  n_scen <- nrow(scenarios)
  
  cat("\n====================================================\n")
  cat("Parallel simulation on", n_cores, "cores\n")
  cat("Output:", outdir, "\n")
  cat("Scenarios:", n_scen, " | Reps per scenario:", n_reps, "\n")
  cat("Batch size:", batch_size, "\n")
  cat("====================================================\n\n")
  
  for (i in seq_len(n_scen)) {
    sc <- as.list(scenarios[i, , drop = FALSE])
    sid <- sprintf("%03d", sc$scenario_id)
    scenario_rds <- file.path(outdir, paste0("scenario_", sid, ".rds"))
    scenario_csv <- file.path(outdir, paste0("scenario_", sid, ".csv"))
    
    cat(sprintf("[%d/%d] scenario_id=%s | n=%d M=%d %s rho=%.1f signal=%.2f %s\n",
                i, n_scen, sid, sc$n, sc$M, sc$corr_type, sc$rho, sc$signal, sc$structure))
    
    # ---- resume or start fresh ----
    if (file.exists(scenario_rds)) {
      res_s <- readRDS(scenario_rds)
      # remove any duplicate rows that may have been written
      if (any(duplicated(res_s$replicate))) {
        cat("  (Fixing duplicated replicates)\n")
        res_s <- res_s[!duplicated(res_s$replicate), ]
        saveRDS(res_s, scenario_rds)
        write.csv(res_s, scenario_csv, row.names = FALSE)
      }
      done_reps <- sort(unique(res_s$replicate))
      start_rep <- if (length(done_reps) == 0) 1 else max(done_reps) + 1
      cat("  Existing file, completed:", length(done_reps), "reps\n")
    } else {
      res_s <- data.frame()
      start_rep <- 1
      cat("  Starting new scenario.\n")
    }
    
    if (start_rep > n_reps) {
      cat("  Already completed. Skipping.\n\n")
      next
    }
    
    # ---- run remaining replicates in parallel batches ----
    rep_seq <- start_rep:n_reps
    batch_starts <- seq(1, length(rep_seq), by = batch_size)
    
    for (bs in batch_starts) {
      be <- min(bs + batch_size - 1, length(rep_seq))
      batch_reps <- rep_seq[bs:be]
      
      batch_res <- foreach(
        r = batch_reps,
        .packages = c("glmnet", "gWQS", "qgcomp", "dplyr"),
        .export = c("run_one_replicate", "generate_dataset",
                    "generate_sigma", "generate_theta",
                    "scale_theta_to_target", "bcr_fit",
                    "bcr_global_test", "bcr_total_effect_ci",
                    "wqs_global_test", "qgc_global_test")
      ) %dorng% {
        bcr_B_use <- if (sc$scenario_type == "null") 499 else bcr_B
        
        tryCatch(
          run_one_replicate(scenario = sc, rep_id = r,
                            bcr_B = bcr_B_use, wqs_b = wqs_b, qgc_B = qgc_B),
          error = function(e) {
            data.frame(
              scenario_id = sc$scenario_id,
              scenario_type = sc$scenario_type,
              n = sc$n, M = sc$M,
              corr_type = sc$corr_type,
              rho = sc$rho,
              signal = sc$signal,
              structure = sc$structure,
              replicate = r,
              fatal_error = as.character(e$message),
              stringsAsFactors = FALSE
            )
          }
        )
      }
      
      new_batch <- dplyr::bind_rows(batch_res)
      res_s <- dplyr::bind_rows(res_s, new_batch)
      
      saveRDS(res_s, scenario_rds)
      write.csv(res_s, scenario_csv, row.names = FALSE)
      cat(sprintf("  saved up to rep %d\n", be))
    }
    
    cat(sprintf("  Scenario %s finished (%d replicates)\n\n", sid, nrow(res_s)))
    
    # update global summaries after each scenario
    update_running_summaries(outdir)
  }
  
  stopCluster(cl)
  cat("\nAll scenarios completed in:", outdir, "\n")
  invisible(outdir)
}


# ============================================================
# Generate all scenarios and split by type
# ============================================================
scenarios_main <- make_main_scenarios()
batches <- split_main_scenarios(scenarios_main)

# ============================================================
# Run simulations (resume automatically from saved progress)
# ============================================================

# 1. null scenarios (54)
run_simulation_parallel(
  scenarios = batches$null,
  n_reps = 1000,
  outdir = "main_null",
  bcr_B = 199,
  wqs_b = 30,
  qgc_B = 50,
  n_cores = 8,
  batch_size = 100
)

# 2. sparse same-sign scenarios (162)
run_simulation_parallel(
  scenarios = batches$same,
  n_reps = 500,
  outdir = "main_same_sign",
  bcr_B = 199,
  wqs_b = 30,
  qgc_B = 50,
  n_cores = 8,
  batch_size = 100
)

# 3. sparse mixed-sign scenarios (162)
run_simulation_parallel(
  scenarios = batches$mixed,
  n_reps = 500,
  outdir = "main_mixed_sign",
  bcr_B = 199,
  wqs_b = 30,
  qgc_B = 50,
  n_cores = 8,
  batch_size = 100
)

cat("All tasks submitted. Check output directories for results.\n")
