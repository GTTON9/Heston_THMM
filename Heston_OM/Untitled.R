initial_heuristic <- function(data, states, positive_mu) {
  
  # ============================
  # 1. K-means Clustering (with NA handling)
  # ============================
  valid_idx <- which(!is.na(data))
  data_clean <- data[valid_idx]
  
  if (length(data_clean) < states * 5) {
    warning("Insufficient data for clustering. Using simple quantile-based initialization.")
    # Fallback: assign states by quantiles
    quantiles <- quantile(data_clean, probs = seq(0, 1, length.out = states + 1))
    cluster_clean <- cut(data_clean, breaks = quantiles, labels = FALSE, include.lowest = TRUE)
  } else {
    cluster_clean <- stats::kmeans(
      data_clean, 
      centers = states, 
      iter.max = 100, 
      nstart = 100
    )$cluster
  }
  
  # Create full cluster vector
  cluster <- rep(NA, length(data))
  cluster[valid_idx] <- cluster_clean
  

  
  # ============================
  # 2. FIXED: Transition Matrix Estimation
  # ============================
  Gamma <- matrix(NA, states, states)
  smooth_factor <- 0.01  # Laplace smoothing parameter
  
  for (i in 1:states) {
    indices_i <- which(cluster[1:(length(cluster) - 1)] == i)
    
    if (length(indices_i) == 0) {
      # No observations in this state, use uniform
      Gamma[i, ] <- 1 / states
      next
    }
    
    # Count transitions
    numerator_row <- numeric(states)
    for (j in 1:states) {
      N_ij <- sum(cluster[indices_i + 1] == j, na.rm = TRUE)
      numerator_row[j] <- N_ij  # ✅ FIXED: No more max()!
    }
    
    # Apply Laplace smoothing
    numerator_row <- numerator_row + smooth_factor
    
    # Normalize
    S_i <- sum(numerator_row)
    Gamma[i, ] <- numerator_row / S_i
  }
  
  # ============================
  # 3. Parameter Estimation
  # ============================
  kappa_est <- numeric(states)
  theta_est <- numeric(states)
  sigma_est <- numeric(states)
  
  dt <- 1/252
  Delta_Xt <- diff(data)
  Xt_lag <- data[-length(data)]
  cluster_lag <- cluster[-length(cluster)]
  
  for (s in seq_len(states)) {
    
    Xt_s <- data[cluster == s]
    
    # Get indices for this state
    indices_s <- which(cluster_lag == s & !is.na(cluster_lag))
    
    if (length(indices_s) < 3) {
      warning(paste("State", s, "has too few observations (<3). Using default parameters."))
      kappa_est[s] <- 2.0
      theta_est[s] <- max(mean(Xt_s, na.rm = TRUE), 1e-6)
      sigma_est[s] <- 0.1
      next
    }
    
    Delta_Xt_s <- Delta_Xt[indices_s]
    Xt_lag_s <- Xt_lag[indices_s]
    
    # ============================
    # OLS Regression: ΔX = α + β·X_{t-1} + ε
    # ============================
    ols_model_s <- tryCatch({
      lm(Delta_Xt_s ~ Xt_lag_s)
    }, error = function(e) {
      warning(paste("OLS failed for state", s, ":", e$message))
      return(NULL)
    })
    
    if (is.null(ols_model_s)) {
      # Fallback to simple estimates
      kappa_est[s] <- 2.0
      theta_est[s] <- max(mean(Xt_s, na.rm = TRUE), 1e-6)
      sigma_est[s] <- max(sd(Xt_s, na.rm = TRUE), 0.01)
      next
    }
    
    beta_hat <- coef(ols_model_s)["Xt_lag_s"]
    alpha_hat <- coef(ols_model_s)["(Intercept)"]
    
    # ============================
    # FIXED: Extract parameters from OLS
    # ============================
    
    # From: ΔX ≈ κ(θ-X)Δt + noise
    #     = κθΔt - κXΔt + noise
    #     = α + βX + noise
    # We have: α = κθΔt, β = -κΔt
    
    # Extract κ
    kappa_s_raw <- -beta_hat / dt
    kappa_est[s] <- max(kappa_s_raw, 0.1)  # Ensure positive
    
    # ✅ FIXED: Extract θ from OLS results
    # θ = α / (κΔt) = -α / β
    if (abs(beta_hat) > 1e-10) {
      theta_s_ols <- -alpha_hat / beta_hat
      theta_est[s] <- max(theta_s_ols, 1e-6)
    } else {
      # Fallback if β ≈ 0
      theta_est[s] <- max(mean(Xt_s, na.rm = TRUE), 1e-6)
    }
    
    # ============================
    # FIXED: Sigma estimation from residuals
    # ============================
    residuals_s <- residuals(ols_model_s)
    
    # For CIR: Var(ΔX | X_{t-1}) = σ²·X_{t-1}·Δt
    # So: σ² ≈ residual² / (X_{t-1}·Δt)
    
    # Use median for robustness
    sigma_sq_candidates <- residuals_s^2 / (Xt_lag_s * dt)
    
    # Remove outliers and non-positive values
    sigma_sq_candidates <- sigma_sq_candidates[
      is.finite(sigma_sq_candidates) & sigma_sq_candidates > 0
    ]
    
    if (length(sigma_sq_candidates) > 0) {
      sigma_sq_est <- median(sigma_sq_candidates, na.rm = TRUE)
      sigma_est[s] <- sqrt(max(sigma_sq_est, 1e-8))
    } else {
      # Fallback
      sigma_est[s] <- max(sd(Delta_Xt_s, na.rm = TRUE) / sqrt(dt * mean(Xt_lag_s)), 0.01)
    }
    
    # Sanity check: σ should be reasonable (0.01 to 1.0)
    sigma_est[s] <- pmax(pmin(sigma_est[s], 1.0), 0.01)
  }

  if (verbose) {
    cat("\nInitial parameter estimates:\n")
    for (s in 1:states) {
      cat(sprintf("  State %d: κ=%.4f, θ=%.4f, σ=%.4f\n", 
                  s, kappa_est[s], theta_est[s], sigma_est[s]))
    }
    cat("\nTransition Matrix:\n")
    print(round(Gamma, 4))
  }
  
  list(
    "cluster" = cluster,
    "pars" = list(
      "kappa" = kappa_est,
      "theta" = theta_est,
      "sigma" = sigma_est,
      "Gamma" = Gamma
    )
  )
}


batch_results <- calculate_expected_moves_batch(S_daily, V_daily)
plot_expected_moves_with_shift(batch_results, H_days = 5)
result_with_paths <- construct_linear_paths(batch_results, H = 5/252, data_with_IV)

controls_em <- Heston_set_controls(
  states = 2,
  sdds = "Heston",
  horizon = nrow(result_with_paths)
)

model_em <- Heston_fit_model_EM(
  result_with_paths = result_with_paths,
  controls = controls_em,
  max_iter = 50,         
  tol = 1e-4,             
  lambda_om = 1.0,        
  H = 5/252,
  N_integration = 100,
  seed = 123,
  verbose = TRUE
)



Gamma_est <- model_em$estimate$Gamma

regime_params_BW <- list(
  list(kappa = model_em$estimate$kappa[1], 
       theta = model_em$estimate$theta[1], 
       sigma = model_em$estimate$sigma[1]),
  list(kappa = model_em$estimate$kappa[2], 
       theta = model_em$estimate$theta[2], 
       sigma = model_em$estimate$sigma[2])
)


result_complete <- calculate_OM_batch(result_with_paths, regime_params_BW)

OM_viterbi_result <- viterbi_om_pure(OM_upper= result_complete[, c("OM_upper_R1", "OM_upper_R2")],
                                     OM_lower= result_complete[, c("OM_lower_R1", "OM_lower_R2")],
                                     nstates,
                                     Gamma_est,
                                     lambda_om = 1.0,
                                     use_upper = TRUE,
                                     use_lower = TRUE)


plot2 <- plot_viterbi_om(
  batch_results_complete = result_complete,
  nstates = 2,
  Gamma = Gamma_est,
  kappa = c(regime_params_BW[[1]]$kappa, regime_params_BW[[2]]$kappa),
  theta = c(regime_params_BW[[1]]$theta, regime_params_BW[[2]]$theta),
  sigma = c(regime_params_BW[[1]]$sigma, regime_params_BW[[2]]$sigma),
  Reg_chain = Reg_chain_year ,
  lambda_om = 1,
  normalize_method = "log",
  show_om_contribution = TRUE
)


result <- plot_cir_confidence_simple(
  true_vol = V_daily,
  param = model_em$estimate,
  states_estimate = OM_viterbi_result$states_estimate,
  dt = 1/252,
  interval_step = 10
)

print(result$plot)


