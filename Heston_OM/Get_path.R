

# ============================================================
# Onsager-Machlup Functional using IV-based Linear Paths
# ============================================================

# ------------------------------------------------------------
# 1. 从IV数据构建线性路径 φ(t) = a + bt
# ------------------------------------------------------------
construct_linear_paths <- function(result, H = 5/252, m_grid) {
 
  m_upper <- result$K_upper / result$S  # 应该是 > 1 的值
  m_lower <- result$K_lower / result$S  # 应该是 < 1 的值
  print(m_lower)
  get_nearest_iv <- function(target_m, iv_type, result_df) {
    pattern <- paste0("^", iv_type, "IV_m")
    cols <- grep(pattern, colnames(result_df), value = TRUE)
    grid_vals <- as.numeric(gsub(pattern, "", cols))
    matched_ivs <- sapply(1:length(target_m), function(i) {
      nearest_idx <- which.min(abs(grid_vals - target_m[i]))
      col_to_use <- cols[nearest_idx]
      cat(sprintf("[%d] Target: %.4f -> Matched: %s\n", i, target_m[i], col_to_use))
      print(col_to_use)
      return(result_df[i, col_to_use])
    })
     # print(matched_ivs)
    return(matched_ivs)
  }
  

  result$sigma_upper_IV <- as.numeric(get_nearest_iv(m_upper, "Call", m_grid))
  result$sigma_lower_IV <- as.numeric(get_nearest_iv(m_lower, "Put", m_grid))


  # 当前方差作为起点
  result$v_start <- result$v

  result$v_end_upper <- (result$sigma_upper_IV)^2

  result$v_end_lower <- (result$sigma_lower_IV)^2

  # 线性路径参数: φ(τ) = a + b*τ
  result$path_upper_a <- result$v_start
  result$path_upper_b <- (result$v_end_upper - result$v_start) / H

  
  result$path_lower_a <- result$v_start
  result$path_lower_b <- (result$v_end_lower - result$v_start) / H
  
  cat("✓ Upper path: v_start → v_end_upper (using CallIV_m1.03)\n")
  cat("✓ Lower path: v_start → v_end_lower (using PutIV_m0.97)\n\n")
  
  return(result)
}

# ------------------------------------------------------------
# 2. 计算单条路径的 O-M Functional
# ------------------------------------------------------------
calculate_OM_single <- function(a, b, kappa, theta, sigma, 
                                H = 5/252, N = 100) {
  # φ(τ) = a + b*τ, τ ∈ [0, H]
  # L[φ] = (1/2) ∫₀ᴴ [(φ̇ - κ(θ-φ))² / (σ²φ)] dτ
  
  # 检查参数有效性
  if (is.na(a) || is.na(b) || a <= 0) {
    return(Inf)
  }
  
  delta_tau <- H / N
  tau <- seq(0, H, by = delta_tau)
  
  f <- numeric(N + 1)
  
  for (i in 1:(N + 1)) {
    tau_i <- tau[i]
    phi_i <- a + b * tau_i
    
    # 避免非正值
    if (phi_i <= 1e-6) {
      return(Inf)
    }
    
    phi_dot <- b  # 导数是常数
    drift <- kappa * (theta - phi_i)
    
    numerator <- (phi_dot - drift)^2
    denominator <- sigma^2 * phi_i
    
    f[i] <- numerator / denominator
  }

  # 梯形积分
  integral <- (delta_tau / 2) * (f[1] + 2 * sum(f[2:N]) + f[N + 1])
  L <- 0.5 * integral

  return(L)
}

# ------------------------------------------------------------
# 3. 批量计算所有路径和regime的 O-M Functional
# ------------------------------------------------------------
calculate_OM_batch <- function(result_with_paths,
                               regime_params,
                               H = 5/252,
                               N = 100) {
  
  T_periods <- nrow(result_with_paths)
  
  cat("\n=== Calculating O-M Functional ===\n")
  cat("Time points:", T_periods, "\n")
  cat("Integration steps:", N, "\n")
  cat("Time window H:", H, "\n\n")
  
  # 预分配
  result_with_paths$OM_upper_R1 <- numeric(T_periods)
  result_with_paths$OM_upper_R2 <- numeric(T_periods)
  result_with_paths$OM_lower_R1 <- numeric(T_periods)
  result_with_paths$OM_lower_R2 <- numeric(T_periods)
  
  pb <- txtProgressBar(min = 0, max = T_periods, style = 3)
  
  for (t in 1:T_periods) {
    # 提取路径参数
    a_upper <- result_with_paths$path_upper_a[t]
    b_upper <- result_with_paths$path_upper_b[t]
    a_lower <- result_with_paths$path_lower_a[t]
    b_lower <- result_with_paths$path_lower_b[t]
    
    # Regime 1 参数
    kappa1 <- regime_params[[1]]$kappa
    theta1 <- regime_params[[1]]$theta
    sigma1 <- regime_params[[1]]$sigma
    print(kappa1)
    # Regime 2 参数
    kappa2 <- regime_params[[2]]$kappa
    theta2 <- regime_params[[2]]$theta
    sigma2 <- regime_params[[2]]$sigma
    
    # 计算4个O-M值 (2条路径 × 2个regime)
    result_with_paths$OM_upper_R1[t] <- calculate_OM_single(
      a_upper, b_upper, kappa1, theta1, sigma1, H, N
    )
    
    result_with_paths$OM_upper_R2[t] <- calculate_OM_single(
      a_upper, b_upper, kappa2, theta2, sigma2, H, N
    )
    
    result_with_paths$OM_lower_R1[t] <- calculate_OM_single(
      a_lower, b_lower, kappa1, theta1, sigma1, H, N
    )
    
    result_with_paths$OM_lower_R2[t] <- calculate_OM_single(
      a_lower, b_lower, kappa2, theta2, sigma2, H, N
    )
    
    setTxtProgressBar(pb, t)
  }
  
  close(pb)
  
  # 统计摘要
  cat("\n\n=== O-M Functional Summary ===\n\n")
  
  cat("Regime 1 (Calm):\n")
  cat("  Upper path - Mean:", sprintf("%.4f", mean(result_with_paths$OM_upper_R1, na.rm = TRUE)), "\n")
  cat("  Lower path - Mean:", sprintf("%.4f", mean(result_with_paths$OM_lower_R1, na.rm = TRUE)), "\n\n")
  
  cat("Regime 2 (Turbulent):\n")
  cat("  Upper path - Mean:", sprintf("%.4f", mean(result_with_paths$OM_upper_R2, na.rm = TRUE)), "\n")
  cat("  Lower path - Mean:", sprintf("%.4f", mean(result_with_paths$OM_lower_R2, na.rm = TRUE)), "\n\n")
  
  # 检查Inf值
  n_inf <- sum(is.infinite(c(result_with_paths$OM_upper_R1, 
                             result_with_paths$OM_upper_R2,
                             result_with_paths$OM_lower_R1,
                             result_with_paths$OM_lower_R2)))
  
  if (n_inf > 0) {
    cat("⚠️  Warning:", n_inf, "O-M values are Inf\n\n")
  }
  
  return(result_with_paths)
}

# ------------------------------------------------------------
# 4. 可视化O-M Functional结果
# ------------------------------------------------------------
plot_OM_results <- function(result_complete) {
  par(mfrow = c(2, 2))
  
  # Plot 1: Upper path O-M
  with(result_complete, {
    plot(time, OM_upper_R1, type = "l", col = "blue", lwd = 2,
         main = "O-M Functional: Upper Path (Call IV 103%)",
         ylab = "L[φ]", xlab = "Time")
    lines(time, OM_upper_R2, col = "red", lwd = 2)
    legend("topright", legend = c("Regime 1 (Calm)", "Regime 2 (Turbulent)"),
           col = c("blue", "red"), lwd = 2, cex = 0.3)
  })
  
  # Plot 2: Lower path O-M
  with(result_complete, {
    plot(time, OM_lower_R1, type = "l", col = "blue", lwd = 2,
         main = "O-M Functional: Lower Path (Put IV 97%)",
         ylab = "L[φ]", xlab = "Time")
    lines(time, OM_lower_R2, col = "red", lwd = 2)
    legend("topright", legend = c("Regime 1 (Calm)", "Regime 2 (Turbulent)"),
           col = c("blue", "red"), lwd = 2, cex = 0.3)
  })
  
  # Plot 3: O-M差异 (upper)
  with(result_complete, {
    diff_upper <- OM_upper_R1 - OM_upper_R2
    plot(time, diff_upper, type = "l", col = "darkgreen", lwd = 2,
         main = "O-M Difference: Upper Path (R1 - R2)",
         ylab = "ΔL[φ]", xlab = "Time")
    abline(h = 0, lty = 2, col = "gray")
  })
  
  # Plot 4: True regime
  with(result_complete, {
    plot(regime, type = 's', col = 'black', lwd = 2,
         main = "True Regime",
         ylim = c(-0.1, 1.1), ylab = "Regime", xlab = "Time")
  })
  
  par(mfrow = c(1, 1))
}
