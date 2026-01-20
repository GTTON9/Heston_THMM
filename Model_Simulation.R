simulate_heston <- function(S0, v0, Reg_series, Reg_param, T, N, M,
                            method = "E", interp = TRUE, substeps = 100,
                            seed = 999, min_var = 1e-6) {
  if (!is.null(seed)) set.seed(seed)
  
  # Check Feller (just warn)
  if (2 * Reg_param[1,2] * Reg_param[1,3] < Reg_param[1,4]^2 ||
      2 * Reg_param[2,2] * Reg_param[2,3] < Reg_param[2,4]^2) {
    message("Warning: Feller condition (2*kappa*theta > sigma^2) NOT satisfied in at least one regime.")
  }
  
  if (interp) {
    N <- N * substeps
    Reg_series <- rep(Reg_series, each = substeps)
  }
  
  dt <- (T / N)
  sqrt_dt <- sqrt(dt)
  method <- toupper(method)
  
  S_paths <- matrix(0, nrow = M, ncol = N + 1)
  V_paths <- matrix(0, nrow = M, ncol = N + 1)
  S_paths[, 1] <- S0
  V_paths[, 1] <- v0
  
  for (i in 1:N) {
    reg_index <- Reg_series[i] + 1
    mu_curr <- Reg_param[reg_index, 1]
    kappa_curr <- Reg_param[reg_index, 2]
    theta_curr <- Reg_param[reg_index, 3]
    sigma_curr <- Reg_param[reg_index, 4]
    rho_curr <- Reg_param[reg_index, 5]
    
    rho_prime <- sqrt(max(0, 1 - rho_curr^2))
    
    Z1 <- rnorm(M)
    Z2 <- rnorm(M)
    dW1 <- Z1 * sqrt_dt
    dW2 <- (rho_curr * Z1 + rho_prime * Z2) * sqrt_dt
    
    V_prev_raw <- V_paths[, i]
    V_prev_pos <- pmax(V_prev_raw, 0)
    sqrt_V_prev <- sqrt(V_prev_pos)
    
    if (method == "E_C") {
      kappa <- kappa_curr; theta <- theta_curr; sigma <- sigma_curr
      exp_k_dt <- exp(-kappa * dt)
      C_val <- (2 * kappa) / ((1 - exp_k_dt) * (sigma^2))
      df_val <- 4 * kappa * theta / sigma^2
      lambda_val <- 2 * C_val * V_prev_pos * exp_k_dt
      Y_sample <- rchisq(M, df = df_val, ncp = lambda_val)
      V_next <- Y_sample / (2 * C_val)
    } else if (method == "E") {
      V_next <- V_prev_raw + kappa_curr * (theta_curr - V_prev_raw) * dt +
        sigma_curr * sqrt_V_prev * dW2
    } else if (method == "M") {
      V_euler_term <- V_prev_raw + kappa_curr * (theta_curr - V_prev_raw) * dt +
        sigma_curr * sqrt_V_prev * dW2
      V_milstein_term <- 0.25 * sigma_curr^2 * (dW2^2 - dt)
      V_next <- V_euler_term + V_milstein_term
    } else {
      stop("Unknown method. Use 'E', 'M', or 'E_C'.")
    }
    
    V_next <- pmax(V_next, min_var)
    V_paths[, i + 1] <- V_next
    
    S_prev <- S_paths[, i]
    S_paths[, i + 1] <- S_prev * exp((mu_curr - 0.5 * V_prev_pos) * dt + sqrt_V_prev * dW1)
  }
  
  if (interp) {
    # return intraday paths (do NOT downsample here; let caller choose)
    S_paths <- as.vector(S_paths)
    V_paths <- as.vector(V_paths)
  } else {
    S_paths <- as.vector(S_paths)
    V_paths <- as.vector(V_paths)
  }
  
  return(list(S_paths = S_paths, V_paths = V_paths, method_used = method))
}




simulate_heston_with_IV <- function(S0 = 100,
                                    v0 = 0.04,
                                    regime_params,
                                    Gamma,
                                    T = 1,
                                    N = 252,
                                    substeps = 1000,  
                                    H = 5/252,
                                    r = 0.02,
                                    strikes_relative = seq(0.85, 1.15, by = 0.001),
                                    method = "E",
                                    seed = 999) {
  
  library(NMOF)
  library(derivmkts)
  
  set.seed(seed)
  Reg_chain <- numeric(N)
  Reg_chain[1] <- 0
  for (t in 2:N) {
    Reg_chain[t] <- sample(c(0, 1), 1, prob = Gamma[Reg_chain[t-1] + 1, ])
  }
  
  Reg_param <- rbind(
    c(regime_params[[1]]$mu, regime_params[[1]]$kappa, regime_params[[1]]$theta,
      regime_params[[1]]$sigma, regime_params[[1]]$rho),
    c(regime_params[[2]]$mu, regime_params[[2]]$kappa, regime_params[[2]]$theta,
      regime_params[[2]]$sigma, regime_params[[2]]$rho)
  )
  
  sim <- simulate_heston(
    S0 = S0, v0 = v0,
    Reg_series = Reg_chain,
    Reg_param = Reg_param,
    T = T, N = N, M = 1,
    method = method,
    interp = TRUE,
    substeps = substeps,
    seed = seed
  )
  
  S <- sim$S_paths
  V <- sim$V_paths
  
  # ---- 4) Daily sampling: take end-of-day points ----
  # intraday length is N*substeps + 1
  sample_indices <- seq(1, N * substeps + 1, by = substeps)
  S_daily_all <- S[sample_indices]
  V_daily_all <- V[sample_indices]
  
  # drop t=0 to keep N daily observations
  S_daily <- S_daily_all[-1]
  V_daily <- V_daily_all[-1]
  
  # ---- 5) Compute IV surface each day ----
  n_strikes <- length(strikes_relative)
  call_IV <- matrix(NA, N, n_strikes)
  put_IV  <- matrix(NA, N, n_strikes)
  
  for (t in 1:N) {
    params <- regime_params[[Reg_chain[t] + 1]]
    K_grid <- S_daily[t] * strikes_relative
    
    for (k in 1:n_strikes) {
      call_price <- tryCatch(
        callHestoncf(S = S_daily[t], X = K_grid[k], tau = H, r = r, q = 0,
                     v0 = V_daily[t], vT = params$theta, rho = params$rho,
                     k = params$kappa, sigma = params$sigma),
        error = function(e) NA
      )
      
      if (!is.na(call_price) && call_price > 1e-6) {
        call_IV[t, k] <- tryCatch(
          bscallimpvol(s = S_daily[t], k = K_grid[k], price = call_price,
                       tt = H, r = r, d = 0),
          error = function(e) NA
        )
        
        put_price <- call_price - S_daily[t] + K_grid[k] * exp(-r * H)
        if (put_price > 1e-6) {
          put_IV[t, k] <- tryCatch(
            bsputimpvol(s = S_daily[t], k = K_grid[k], price = put_price,
                        tt = H, r = r, d = 0),
            error = function(e) NA
          )
        }
      }
    }
  }
  
  colnames(call_IV) <- paste0("CallIV_m", strikes_relative)
  colnames(put_IV)  <- paste0("PutIV_m", strikes_relative)
  
  # ---- 6) Realised variance from intraday S ----
  RV_V <- numeric(N)
  for (t in 1:N) {
    idx <- ((t - 1) * substeps + 1):(t * substeps + 1)  # include endpoints
    S_day <- S[idx]
    r_day <- diff(log(S_day))
    RV_V[t] <- sum(r_day^2) * 252
  }
  
  result <- cbind(
    data.frame(time = 1:N, S = S_daily, v = V_daily, regime = Reg_chain),
    call_IV,
    put_IV,
    RV_V
  )
  
  return(result)
}










simulate_heston_with_IV_Heston <- function(
    S0 = 100,
    v0 = 0.1,
    regime_params,
    Gamma,
    T = 1,
    N = 252,
    substeps = 1000,
    H = 5/252,
    r = 0.02,
    q = 0.00,
    strikes_relative = seq(0.85, 1.15, by = 0.01),
    method = "E",
    seed = 999
) {
  library(NMOF)
  library(derivmkts)
  
  set.seed(seed)
  
  # ---- 1) Regime chain ----
  Reg_chain <- numeric(N)
  Reg_chain[1] <- 0
  for (t in 2:N) {
    Reg_chain[t] <- sample(c(0, 1), 1, prob = Gamma[Reg_chain[t-1] + 1, ])
  }
  
  # ---- 2) Build Reg_param for simulate_heston ----
  # Use risk-neutral drift mu = r - q for consistency with option pricing under Q
  Reg_param <- rbind(
    c(r - q, regime_params[[1]]$kappa, regime_params[[1]]$theta,
      regime_params[[1]]$sigma, regime_params[[1]]$rho),
    c(r - q, regime_params[[2]]$kappa, regime_params[[2]]$theta,
      regime_params[[2]]$sigma, regime_params[[2]]$rho)
  )
  
  # ---- 3) Intraday simulation ----
  sim <- simulate_heston(
    S0 = S0, v0 = v0,
    Reg_series = Reg_chain,
    Reg_param = Reg_param,
    T = T, N = N, M = 1,
    method = method,
    interp = TRUE,
    substeps = substeps,
    seed = seed
  )
  
  S <- sim$S_paths
  V <- sim$V_paths
  
  # ---- 4) Daily sampling (end of each day) ----
  sample_indices <- seq(1, N * substeps + 1, by = substeps)
  S_daily_all <- S[sample_indices]
  V_daily_all <- V[sample_indices]
  
  S_daily <- S_daily_all[-1]
  V_daily <- V_daily_all[-1]
  
  # ---- 5) IV surface each day via Heston CF pricing ----
  n_strikes <- length(strikes_relative)
  call_IV <- matrix(NA_real_, N, n_strikes)
  put_IV  <- matrix(NA_real_, N, n_strikes)
  
  # Helper: safe inversion
  safe_call_iv <- function(S, K, price, tt, r, d) {
    if (is.na(price) || !is.finite(price) || price <= 1e-12) return(NA_real_)
    tryCatch(
      derivmkts::bscallimpvol(s = S, k = K, price = price, tt = tt, r = r, d = d),
      error = function(e) NA_real_
    )
  }
  safe_put_iv <- function(S, K, price, tt, r, d) {
    if (is.na(price) || !is.finite(price) || price <= 1e-12) return(NA_real_)
    tryCatch(
      derivmkts::bsputimpvol(s = S, k = K, price = price, tt = tt, r = r, d = d),
      error = function(e) NA_real_
    )
  }
  
  for (t in 1:N) {
    params <- regime_params[[Reg_chain[t] + 1]]
    
    S_t <- S_daily[t]
    v_t <- V_daily[t]
    K_grid <- S_t * strikes_relative
    
    # --- Heston CF call pricing ---
    # Try vectorized pricing first
    call_price_vec <- tryCatch(
      NMOF::callHestoncf(
        S = S_t, X = K_grid, tau = H, r = r, q = q,
        v0 = v_t, vT = params$theta,
        rho = params$rho, k = params$kappa, sigma = params$sigma
      ),
      error = function(e) NA
    )
    
    # If vectorized call failed, fall back to per-strike
    if (length(call_price_vec) == 1 && is.na(call_price_vec)) {
      call_price_vec <- sapply(K_grid, function(K) {
        tryCatch(
          NMOF::callHestoncf(
            S = S_t, X = K, tau = H, r = r, q = q,
            v0 = v_t, vT = params$theta,
            rho = params$rho, k = params$kappa, sigma = params$sigma
          ),
          error = function(e) NA_real_
        )
      })
    }
    
    # Put by parity (with dividend yield q)
    # P = C - S*exp(-qH) + K*exp(-rH)
    put_price_vec <- call_price_vec - S_t * exp(-q * H) + K_grid * exp(-r * H)
    
    # Invert to IV (BS)
    for (k in 1:n_strikes) {
      call_IV[t, k] <- safe_call_iv(S_t, K_grid[k], call_price_vec[k], H, r, d = q)
      put_IV[t, k]  <- safe_put_iv( S_t, K_grid[k], put_price_vec[k],  H, r, d = q)
    }
  }
  
  colnames(call_IV) <- paste0("CallIV_m", strikes_relative)
  colnames(put_IV)  <- paste0("PutIV_m", strikes_relative)
  
  # ---- 6) Realized variance from intraday S ----
  RV_V <- numeric(N)
  for (t in 1:N) {
    idx <- ((t - 1) * substeps + 1):(t * substeps + 1)
    S_day <- S[idx]
    r_day <- diff(log(S_day))
    RV_V[t] <- sum(r_day^2) * 252
  }
  
  result <- cbind(
    data.frame(time = 1:N, S = S_daily, v = V_daily, regime = Reg_chain),
    call_IV, put_IV,
    RV_V = RV_V
  )
  
  return(result)
}




