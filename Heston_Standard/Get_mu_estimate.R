estimate_mu_improved <- function(S_daily, S_intraday, n_intraday,
                                 states_estimate, nstates = 2, 
                                 dt = 1/252,
                                 use_weights = TRUE) {
  
  # Calculate RV
  V_daily <- calculate_daily_rv(S_intraday, n_intraday, 
                                annualize = TRUE, method = "log")
  V_daily_smoothed <- lowess(V_daily, f = 0.1)$y
  
  # Log returns
  log_returns <- diff(log(S_daily))
  n <- length(log_returns)
  
  # ===== 改进1: 时间对齐 =====
  # log_returns[t]是从t到t+1，应该用t时刻的variance
  V_aligned <- V_daily_smoothed[1:n]
  states_aligned <- states_estimate[1:n]
  
  mu_estimates <- numeric(nstates)
  mu_std_errors <- numeric(nstates)
  
  for (state in 1:nstates) {
    state_idx <- which(states_aligned == state)
    
    if (length(state_idx) == 0) {
      mu_estimates[state] <- NA
      mu_std_errors[state] <- NA
      next
    }
    
    state_returns <- log_returns[state_idx]
    state_vols <- V_aligned[state_idx]
    
    # ===== 改进2: 异常值处理 =====
    # 移除极端RV值（可能是估计误差）
    valid_idx <- state_vols > 0 & state_vols < quantile(state_vols, 0.99)
    state_returns <- state_returns[valid_idx]
    state_vols <- state_vols[valid_idx]
    
    if (use_weights) {
      # ===== 改进3: 加权估计 =====
      weights <- 1 / (state_vols * dt + 1e-8)
      weights <- weights / sum(weights)
      
      weighted_mean_return <- sum(weights * state_returns)
      weighted_mean_vol <- sum(weights * state_vols)
      
      mu_estimates[state] <- (weighted_mean_return / dt) + 0.5 * weighted_mean_vol
      
      # Weighted standard error
      residuals <- state_returns - (mu_estimates[state] - 0.5*state_vols)*dt
      weighted_var <- sum(weights * residuals^2)
      mu_std_errors[state] <- sqrt(weighted_var / length(state_returns)) / dt
      
    } else {
      # 原始方法（简单平均）
      mean_return <- mean(state_returns)
      mean_vol <- mean(state_vols)
      
      mu_estimates[state] <- (mean_return / dt) + (mean_vol / 2)
      
      sd_return <- sd(state_returns)
      mu_std_errors[state] <- (sd_return / sqrt(length(state_returns))) / dt
    }
    
    cat(sprintf("State %d:\n", state))
    cat(sprintf("  μ̂ = %.4f\n", mu_estimates[state]))
    cat(sprintf("  SE = %.4f\n", mu_std_errors[state]))
    cat(sprintf("  95%% CI: [%.4f, %.4f]\n", 
                mu_estimates[state] - 1.96*mu_std_errors[state],
                mu_estimates[state] + 1.96*mu_std_errors[state]))
    cat(sprintf("  N = %d observations (after filtering)\n\n", 
                length(state_returns)))
  }
  
  return(list(
    mu = mu_estimates,
    se = mu_std_errors,
    V_daily = V_daily_smoothed,
    method = ifelse(use_weights, "Improved (Weighted)", "Improved (Unweighted)")
  ))
}



a <- estimate_mu_improved(S_daily, S_simulated, n_intraday,
                                 states_estimate, nstates = 2, 
                                 dt = 1/252,
                                 use_weights = TRUE)
