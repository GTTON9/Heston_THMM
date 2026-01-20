# ============================================================
# 修复版：估计ρ（添加NA检查）
# ============================================================

#' 估计ρ - 价格和波动率过程的相关系数（修复版）
#'
#' @description
#' 修复：添加NA/NaN检查，确保数值稳定性
estimate_rho_from_residuals <- function(S_daily, V_daily, states_estimate,
                                        mu_estimates, param_cir,
                                        nstates = 2, dt = 1/252,
                                        method = "spearman") {
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("估计ρ（基于价格和波动率残差）\n")
  cat(rep("=", 70), "\n\n")
  
  cat("方法：", method, "相关系数\n")
  cat("理论：ρ = Corr(dW^S, dW^V)\n\n")
  
  # 计算收益率和波动率变化
  log_returns <- diff(log(S_daily))
  V_changes <- diff(V_daily)
  n <- length(log_returns)
  
  # 对齐数据
  V_daily_lag <- V_daily[1:n]
  V_daily_lead <- V_daily[2:(n+1)]
  states_estimate <- states_estimate[1:n]
  
  # 为每个regime估计ρ
  rho_estimates <- numeric(nstates)
  rho_std_errors <- numeric(nstates)
  
  for (state in 1:nstates) {
    
    cat(sprintf("=== Regime %d ===\n", state))
    
    # 找到该regime的所有时间点
    state_indices <- which(states_estimate == state)
    
    if (length(state_indices) < 20) {
      warning(sprintf("Regime %d: 样本量太少 (%d < 20)，ρ设为0", 
                      state, length(state_indices)))
      rho_estimates[state] <- 0
      rho_std_errors[state] <- NA
      cat(sprintf("⚠️ 样本量不足，跳过\n\n"))
      next
    }
    
    cat(sprintf("样本量: %d\n", length(state_indices)))
    
    # 检查μ是否为NA
    if (is.na(mu_estimates[state])) {
      warning(sprintf("Regime %d: μ为NA，ρ设为0", state))
      rho_estimates[state] <- 0
      rho_std_errors[state] <- NA
      cat(sprintf("⚠️ μ为NA，跳过\n\n"))
      next
    }
    
    # 提取该regime的参数
    mu_state <- mu_estimates[state]
    kappa_state <- param_cir$kappa[state]
    theta_state <- param_cir$theta[state]
    sigma_state <- param_cir$sigma[state]
    
    cat(sprintf("参数: μ=%.4f, κ=%.4f, θ=%.4f, σ=%.4f\n",
                mu_state, kappa_state, theta_state, sigma_state))
    
    # 初始化残差向量
    price_residuals <- numeric(length(state_indices))
    vol_residuals <- numeric(length(state_indices))
    
    # 计算残差
    for (idx in seq_along(state_indices)) {
      i <- state_indices[idx]
      
      # ========== 价格过程残差 ==========
      expected_return <- (mu_state - V_daily_lag[i]/2) * dt
      actual_return <- log_returns[i]
      price_res_raw <- actual_return - expected_return
      
      # 标准化（添加保护）
      price_sd <- sqrt(max(V_daily_lag[i] * dt, 1e-10))
      price_residuals[idx] <- price_res_raw / price_sd
      
      # ========== 波动率过程残差 ==========
      expected_vol_change <- kappa_state * (theta_state - V_daily_lag[i]) * dt
      actual_vol_change <- V_changes[i]
      vol_res_raw <- actual_vol_change - expected_vol_change
      
      # 标准化（添加保护）
      vol_sd <- sigma_state * sqrt(max(V_daily_lag[i] * dt, 1e-10))
      vol_residuals[idx] <- vol_res_raw / vol_sd
    }
    
    # ========== 检查残差有效性 ==========
    
    # 移除NA/NaN/Inf
    valid_idx <- is.finite(price_residuals) & is.finite(vol_residuals)
    
    if (sum(valid_idx) < 20) {
      warning(sprintf("Regime %d: 有效残差太少 (%d < 20)，ρ设为0", 
                      state, sum(valid_idx)))
      rho_estimates[state] <- 0
      rho_std_errors[state] <- NA
      cat(sprintf("⚠️ 有效残差不足，跳过\n\n"))
      next
    }
    
    # 使用有效残差
    price_residuals <- price_residuals[valid_idx]
    vol_residuals <- vol_residuals[valid_idx]
    
    # 检查标准差
    sd_price <- sd(price_residuals, na.rm = TRUE)
    sd_vol <- sd(vol_residuals, na.rm = TRUE)
    
    # 添加NA检查
    if (is.na(sd_price) || is.na(sd_vol)) {
      warning(sprintf("Regime %d: 标准差为NA，ρ设为0", state))
      rho_estimates[state] <- 0
      rho_std_errors[state] <- NA
      cat(sprintf("⚠️ 标准差为NA，跳过\n\n"))
      next
    }
    
    if (sd_price < 1e-10 || sd_vol < 1e-10) {
      warning(sprintf("Regime %d: 残差标准差过小 (sd_price=%.2e, sd_vol=%.2e)，ρ设为0",
                      state, sd_price, sd_vol))
      rho_estimates[state] <- 0
      rho_std_errors[state] <- NA
      cat(sprintf("⚠️ 标准差过小，跳过\n\n"))
      next
    }
    
    # ========== 计算相关系数 ==========
    
    tryCatch({
      
      if (method == "spearman") {
        rho_estimates[state] <- cor(price_residuals, vol_residuals, 
                                    method = "spearman", use = "complete.obs")
      } else if (method == "pearson") {
        rho_estimates[state] <- cor(price_residuals, vol_residuals, 
                                    method = "pearson", use = "complete.obs")
      } else if (method == "kendall") {
        rho_estimates[state] <- cor(price_residuals, vol_residuals, 
                                    method = "kendall", use = "complete.obs")
      }
      
      # 检查结果是否为NA
      if (is.na(rho_estimates[state])) {
        warning(sprintf("Regime %d: 相关系数计算返回NA，设为0", state))
        rho_estimates[state] <- 0
        rho_std_errors[state] <- NA
        next
      }
      
      # 限制在 [-1, 1] 范围
      rho_estimates[state] <- max(min(rho_estimates[state], 0.99), -0.99)
      
      # 估计标准误
      n_obs <- length(price_residuals)
      if (method == "pearson") {
        rho_std_errors[state] <- 1 / sqrt(n_obs - 3)
      } else {
        rho_std_errors[state] <- 1 / sqrt(n_obs)
      }
      
      # 输出统计
      cat("\n残差统计:\n")
      cat(sprintf("  有效样本: %d\n", n_obs))
      cat(sprintf("  价格残差: mean=%.4f, sd=%.4f, range=[%.2f, %.2f]\n",
                  mean(price_residuals, na.rm = TRUE), 
                  sd(price_residuals, na.rm = TRUE),
                  min(price_residuals, na.rm = TRUE), 
                  max(price_residuals, na.rm = TRUE)))
      cat(sprintf("  波动率残差: mean=%.4f, sd=%.4f, range=[%.2f, %.2f]\n",
                  mean(vol_residuals, na.rm = TRUE), 
                  sd(vol_residuals, na.rm = TRUE),
                  min(vol_residuals, na.rm = TRUE), 
                  max(vol_residuals, na.rm = TRUE)))
      
      # ρ估计结果
      cat(sprintf("\nρ估计 (%s): %.4f ± %.4f\n", 
                  method, rho_estimates[state], rho_std_errors[state]))
      
      # 95%置信区间
      if (!is.na(rho_std_errors[state]) && method == "pearson") {
        z <- 0.5 * log((1 + rho_estimates[state]) / (1 - rho_estimates[state]))
        se_z <- rho_std_errors[state]
        z_lower <- z - 1.96 * se_z
        z_upper <- z + 1.96 * se_z
        
        ci_lower <- (exp(2*z_lower) - 1) / (exp(2*z_lower) + 1)
        ci_upper <- (exp(2*z_upper) - 1) / (exp(2*z_upper) + 1)
        
        cat(sprintf("95%% CI: [%.4f, %.4f]\n", ci_lower, ci_upper))
      }
      
    }, error = function(e) {
      warning(sprintf("Regime %d: 相关系数计算出错: %s，ρ设为0", state, e$message))
      rho_estimates[state] <<- 0
      rho_std_errors[state] <<- NA
    })
    
    cat("\n")
  }
  
  return(list(
    rho = rho_estimates,
    se = rho_std_errors,
    method = method
  ))
}


# ============================================================
# 重新运行（使用修复版）
# ============================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("重新运行ρ估计（修复版）\n")
cat(rep("=", 80), "\n", sep = "")

# 首先检查result_mle的结果
cat("\n检查μ估计结果:\n")
cat("μ estimates:", result_mle$mu, "\n")
cat("μ有NA?:", any(is.na(result_mle$mu)), "\n")
cat("V_daily长度:", length(result_mle$V_daily), "\n")
cat("V_daily有NA?:", any(is.na(result_mle$V_daily)), "\n")
cat("states_estimate长度:", length(states_estimate), "\n\n")

# 运行修复版
rho_result <- estimate_rho_from_residuals(
  S_daily = S_daily,
  V_daily = result_mle$V_daily,
  states_estimate = states_estimate,
  mu_estimates = result_mle$mu,
  param_cir = param,
  nstates = 2,
  dt = 1/252,
  method = "spearman"
)

# 显示结果
cat("\n", rep("=", 80), "\n", sep = "")
cat("ρ估计结果\n")
cat(rep("=", 80), "\n\n")

for (i in 1:2) {
  cat(sprintf("Regime %d: ρ = %.4f", i, rho_result$rho[i]))
  if (!is.na(rho_result$se[i])) {
    cat(sprintf(" ± %.4f", rho_result$se[i]))
  }
  cat("\n")
}

cat("\n真实ρ:", c(-0.1, -0.1), "\n")
cat("估计误差:", rho_result$rho - c(-0.1, -0.1), "\n")


# ============================================================
# 如果仍然出错，使用更保守的方法
# ============================================================

if (all(rho_result$rho == 0)) {
  
  cat("\n", rep("=", 80), "\n", sep = "")
  cat("⚠️ 标准方法失败，尝试保守方法\n")
  cat(rep("=", 80), "\n\n")
  
  # 保守方法：不分regime，计算全局ρ
  log_returns <- diff(log(S_daily))
  V_changes <- diff(result_mle$V_daily)
  n <- length(log_returns)
  
  V_daily_lag <- result_mle$V_daily[1:n]
  
  # 使用全局μ估计
  global_mu <- mean(result_mle$mu, na.rm = TRUE)
  
  # 使用全局CIR参数
  global_kappa <- mean(param$kappa)
  global_theta <- mean(param$theta)
  global_sigma <- mean(param$sigma)
  
  # 计算全局残差
  price_residuals <- numeric(n)
  vol_residuals <- numeric(n)
  
  for (i in 1:n) {
    # 价格残差
    expected_return <- (global_mu - V_daily_lag[i]/2) * (1/252)
    price_res_raw <- log_returns[i] - expected_return
    price_sd <- sqrt(max(V_daily_lag[i] / 252, 1e-10))
    price_residuals[i] <- price_res_raw / price_sd
    
    # 波动率残差
    expected_vol_change <- global_kappa * (global_theta - V_daily_lag[i]) / 252
    vol_res_raw <- V_changes[i] - expected_vol_change
    vol_sd <- global_sigma * sqrt(max(V_daily_lag[i] / 252, 1e-10))
    vol_residuals[i] <- vol_res_raw / vol_sd
  }
  
  # 移除NA
  valid_idx <- is.finite(price_residuals) & is.finite(vol_residuals)
  price_residuals <- price_residuals[valid_idx]
  vol_residuals <- vol_residuals[valid_idx]
  
  cat(sprintf("有效残差数: %d\n", length(price_residuals)))
  
  # 计算全局ρ
  if (length(price_residuals) > 20) {
    global_rho <- cor(price_residuals, vol_residuals, 
                      method = "spearman", use = "complete.obs")
    
    cat(sprintf("\n全局ρ (Spearman): %.4f\n", global_rho))
    cat("将此值用于两个regime\n")
    
    # 更新结果
    rho_result$rho <- rep(global_rho, 2)
    rho_result$se <- rep(1/sqrt(length(price_residuals)), 2)
    rho_result$method <- "spearman (global)"
  }
}

# 保存结果
cat("\n最终ρ估计:\n")
print(data.frame(
  Regime = 1:2,
  rho = rho_result$rho,
  se = rho_result$se,
  True = c(-0.1, -0.1),
  Error = rho_result$rho - c(-0.1, -0.1)
))
cat("\n=== 详细诊断 ===\n\n")

# 1. 检查μ
cat("1. μ估计:\n")
cat("  值:", result_mle$mu, "\n")
cat("  类型:", class(result_mle$mu), "\n")
cat("  长度:", length(result_mle$mu), "\n")
cat("  有NA?:", any(is.na(result_mle$mu)), "\n")
cat("  有Inf?:", any(is.infinite(result_mle$mu)), "\n\n")

# 2. 检查V_daily
cat("2. V_daily:\n")
cat("  长度:", length(result_mle$V_daily), "\n")
cat("  范围:", range(result_mle$V_daily, na.rm = TRUE), "\n")
cat("  有NA?:", sum(is.na(result_mle$V_daily)), "/", length(result_mle$V_daily), "\n")
cat("  有负值?:", sum(result_mle$V_daily < 0, na.rm = TRUE), "\n")
cat("  有Inf?:", sum(is.infinite(result_mle$V_daily)), "\n\n")

# 3. 检查states
cat("3. states_estimate:\n")
cat("  长度:", length(states_estimate), "\n")
cat("  唯一值:", unique(states_estimate), "\n")
cat("  分布:\n")
print(table(states_estimate))
cat("\n")

# 4. 检查S_daily
cat("4. S_daily:\n")
cat("  长度:", length(S_daily), "\n")
cat("  范围:", range(S_daily, na.rm = TRUE), "\n")
cat("  有NA?:", sum(is.na(S_daily)), "\n\n")

# 5. 试算第一个regime的残差
cat("5. 试算Regime 1的残差:\n")
state1_idx <- which(states_estimate[1:(length(S_daily)-1)] == 1)
cat("  样本量:", length(state1_idx), "\n")

if (length(state1_idx) > 0 && !is.na(result_mle$mu[1])) {
  log_returns <- diff(log(S_daily))
  V_lag <- result_mle$V_daily[1:length(log_returns)]
  
  # 第一个残差
  i <- state1_idx[1]
  expected_return <- (result_mle$mu[1] - V_lag[i]/2) * (1/252)
  price_res <- (log_returns[i] - expected_return) / sqrt(V_lag[i] / 252)
  
  cat("  第一个价格残差:", price_res, "\n")
  cat("  是否有限?:", is.finite(price_res), "\n")
}

cat("\n=== 诊断完成 ===\n")