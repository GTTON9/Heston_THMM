

#' Forward Filter for Pure O-M Functional (THMM)
#' 
#' Following PDF Section 11.3.1 - Forward Probabilities
#' Formula: α_k(j) = [Σ_i α_{k-1}(i) P_ij] exp(-L_j[φ_k])
#'
#' @param OM_upper Matrix [T x nstates] of O-M functionals for upper paths
#' @param OM_lower Matrix [T x nstates] of O-M functionals for lower paths
#' @param nstates Number of hidden states
#' @param Gamma Transition probability matrix
#' @param lambda_om Weight for O-M contribution (default 1.0)
#' @param use_upper Whether to use upper path (default TRUE)
#' @param use_lower Whether to use lower path (default TRUE)
#'
#' @return List containing:
#'   - posterior: Matrix [nstates x T] of P(s_t = i | φ_{1:t})
#'   - alpha: Matrix [nstates x T] of forward probabilities
#'   - states_estimate: Vector [T] of most likely states at each time
#'   - L_total: Matrix [nstates x T] of total O-M costs

forward_filter_om_pure <- function(OM_upper,
                                   OM_lower,
                                   nstates,
                                   Gamma,
                                   lambda_om = 1.0,
                                   use_upper = TRUE,
                                   use_lower = TRUE) {
  
  T <- nrow(OM_upper)
  
  cat("=================================================================\n")
  cat("Pure O-M Forward Filter (THMM)\n")
  cat("=================================================================\n\n")

  
  cat("Computing O-M functionals...\n")
  
  L_total <- matrix(0, T, nstates)
  
  for (t in 1:T) {
    for (i in 1:nstates) {
      
      om_value <- 0
      
      # Upper path
      if (use_upper && !is.na(OM_upper[t, i]) && is.finite(OM_upper[t, i])) {
        om_value <- om_value + lambda_om * OM_upper[t, i]
      }
      
      # Lower path
      if (use_lower && !is.na(OM_lower[t, i]) && is.finite(OM_lower[t, i])) {
        om_value <- om_value + lambda_om * OM_lower[t, i]
      }
      
      L_total[t, i] <- om_value
    }
  }
  
  cat("  L_total statistics:\n")
  cat("    Mean:", mean(L_total), "\n")
  cat("    Range: [", min(L_total), ",", max(L_total), "]\n\n")
  
  delta <- oeli::stationary_distribution(Gamma) # initial distribution
  
  cat("Initial distribution (π):\n")
  for (i in 1:nstates) {
    cat(sprintf("  State %d: %.4f\n", i, delta[i]))
  }
  cat("\n")
  
  log_alpha <- matrix(-Inf, nstates, T)

  
  for (i in 1:nstates) {
    log_alpha[i, 1] <- log(delta[i]) - L_total[1, i]
  }
  
  pb <- txtProgressBar(min = 2, max = T, style = 3)
  
  log_Gamma <- log(Gamma)
  
  for (t in 2:T) {
    for (j in 1:nstates) {
      
      log_terms <- numeric(nstates)
      for (i in 1:nstates) {
        log_terms[i] <- log_alpha[i, t-1] + log_Gamma[i, j]
      }
      log_sum <- logsumexp(log_terms)
      
      log_alpha[j, t] <- log_sum - L_total[t, j]
    }
    
    setTxtProgressBar(pb, t)
  }
  
  close(pb)
  
  cat("\n\nComputing posterior probabilities...\n")
  
  posterior <- matrix(0, nstates, T)
  
  for (t in 1:T) {
    log_norm <- logsumexp(log_alpha[, t])
    
    for (i in 1:nstates) {
      posterior[i, t] <- exp(log_alpha[i, t] - log_norm)
    }
  }
  
  states_estimate <- apply(posterior, 2, which.max)
  
  
  transitions <- diff(states_estimate)
  n_switches <- sum(transitions != 0)
  
  cat("\n=== Forward Filter Results (Pure O-M) ===\n")
  cat("Total time windows:", T, "\n")
  cat("Number of regime switches:", n_switches, "\n")
  
  for (i in 1:nstates) {
    n_in_regime <- sum(states_estimate == i)
    pct <- 100 * n_in_regime / T
    avg_posterior <- mean(posterior[i, states_estimate == i], na.rm = TRUE)
    
    cat(sprintf("  Regime %d:\n", i))
    cat(sprintf("    Frequency: %d windows (%.1f%%)\n", n_in_regime, pct))
    cat(sprintf("    Avg posterior when assigned: %.4f\n", avg_posterior))
  }
  
  cat("\nDetected switches:\n")
  switch_idx <- which(transitions != 0)
  if (length(switch_idx) > 0) {
    for (idx in switch_idx) {
      cat(sprintf("  t=%d: R%d → R%d (prob: %.3f → %.3f)\n", 
                  idx, 
                  states_estimate[idx], 
                  states_estimate[idx+1],
                  posterior[states_estimate[idx], idx],
                  posterior[states_estimate[idx+1], idx+1]))
    }
  } else {
    cat("  (No regime switches detected)\n")
  }
  
  # Compute average posterior confidence
  avg_confidence <- mean(apply(posterior, 2, max))
  cat(sprintf("\nAverage posterior confidence: %.4f\n", avg_confidence))
  
  cat("=================================================================\n")
  
  return(list(
    posterior = posterior,
    alpha = exp(log_alpha),
    log_alpha = log_alpha,
    states_estimate = states_estimate,
    L_total = L_total,
    n_switches = n_switches,
    switch_points = switch_idx,
    avg_confidence = avg_confidence
  ))
}



logsumexp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}





viterbi_om_pure <- function(OM_upper,
                            OM_lower,
                            nstates,
                            Gamma,
                            lambda_om = 1.0,
                            use_upper = TRUE,
                            use_lower = TRUE) {
  
  T_len <- nrow(OM_upper)
  log_Gamma <- log(Gamma)
  
  L_total <- matrix(0, T_len, nstates)
  for (t in 1:T_len) {
    for (i in 1:nstates) {
      val <- 0
      if (use_upper) val <- val + lambda_om * OM_upper[t, i]
      if (use_lower) val <- val + lambda_om * OM_lower[t, i]
      L_total[t, i] <- val
    }
  }
  
  v <- matrix(-Inf, nstates, T_len)
  psi <- matrix(0, nstates, T_len)
  delta <- oeli::stationary_distribution(Gamma)
  
  for (i in 1:nstates) {
    v[i, 1] <- log(delta[i]) - L_total[1, i]
  }
  
  for (t in 2:T_len) {
    for (j in 1:nstates) {
      
      trans_probs <- v[, t-1] + log_Gamma[, j]
      best_prev_state <- which.max(trans_probs)
      
      v[j, t] <- trans_probs[best_prev_state] - L_total[t, j]
      psi[j, t] <- best_prev_state
    }
  }
  
  states_viterbi <- numeric(T_len)
  states_viterbi[T_len] <- which.max(v[, T_len])
  
 
  for (t in (T_len - 1):1) {
    states_viterbi[t] <- psi[states_viterbi[t + 1], t + 1]
  }
  
  cat("Viterbi decoding complete.\n")
  return(list(
    states_estimate = states_viterbi,
    viterbi_score = v,
    path_matrix = psi
  ))
}











# ================================================================
# Plot Viterbi Results with O-M Functional (Probability Scale)
# ================================================================




plot_viterbi_om_path <- function(batch_results_complete,
                                 nstates, 
                                 Gamma, 
                                 Reg_chain,
                                 lambda_om = 1.0,
                                 use_upper = TRUE,
                                 use_lower = TRUE,
                                 normalize_method = "log") {
  
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  
  
  cat("Preparing O-M matrices...\n")
  
  if (normalize_method == "log") {
    OM_upper <- log(batch_results_complete[, c("OM_upper_R1", "OM_upper_R2")] + 1)
    OM_lower <- log(batch_results_complete[, c("OM_lower_R1", "OM_lower_R2")] + 1)
  } else if (normalize_method == "zscore") {
    OM_all <- c(batch_results_complete$OM_upper_R1, batch_results_complete$OM_upper_R2,
                batch_results_complete$OM_lower_R1, batch_results_complete$OM_lower_R2)
    OM_finite <- OM_all[is.finite(OM_all)]
    mean_om <- mean(OM_finite); sd_om <- sd(OM_finite)
    OM_upper <- (batch_results_complete[, c("OM_upper_R1", "OM_upper_R2")] - mean_om) / sd_om
    OM_lower <- (batch_results_complete[, c("OM_lower_R1", "OM_lower_R2")] - mean_om) / sd_om
  } else {
    OM_upper <- batch_results_complete[, c("OM_upper_R1", "OM_upper_R2")]
    OM_lower <- batch_results_complete[, c("OM_lower_R1", "OM_lower_R2")]
  }
  
  OM_upper <- as.matrix(OM_upper)
  OM_lower <- as.matrix(OM_lower)
  
  # ===== Step 2: 运行 Viterbi 算法 (替换原有的 Forward Filter) =====
  cat("Running Viterbi Algorithm with O-M...\n")
  
  # 调用之前定义的 viterbi_om_pure 函数
  viterbi_results <- viterbi_om_pure(
    OM_upper = OM_upper,
    OM_lower = OM_lower,
    nstates = nstates,
    Gamma = Gamma,
    lambda_om = lambda_om,
    use_upper = use_upper,
    use_lower = use_lower
  )
  
  # Viterbi 得到的是 1, 2 序列，转为 0, 1 以匹配 Reg_chain
  viterbi_path <- viterbi_results$states_estimate - 1
  
  # ===== Step 3: 准备绘图数据 =====
  time_index <- 1:length(viterbi_path)
  
  state_data <- data.frame(
    Time = time_index,
    Viterbi_State = viterbi_path,
    True_State = as.factor(Reg_chain[1:length(viterbi_path)])
  )
  
  # 背景矩形数据 (保持原有逻辑)
  regime_background <- state_data %>%
    mutate(change = True_State != lag(True_State, default = first(True_State))) %>%
    filter(change | row_number() == 1) %>%
    dplyr::select(Time, True_State) %>%
    mutate(Time_End = lead(Time, default = max(time_index) + 1)) %>%
    rename(Time_Start = Time)
  
  # ===== Step 4: 创建 Viterbi 路径图 =====
  
  p <- ggplot() +
    # 背景着色：代表真实的 Regime
    geom_rect(
      data = regime_background,
      aes(xmin = Time_Start, xmax = Time_End, ymin = -0.05, ymax = 1.05, fill = True_State),
      alpha = 0.15
    ) +
    
    # 真实的 Regime 阶梯线 (虚线)
    geom_step(
      data = state_data,
      aes(x = Time, y = as.numeric(as.character(True_State)), linetype = "True Path"),
      color = "gray40", linewidth = 0.8
    ) +
    
    # Viterbi 识别出的最可能路径 (粗实线)
    geom_step(
      data = state_data,
      aes(x = Time, y = Viterbi_State, color = "Viterbi Path"),
      linewidth = 1.2
    ) +
    
    # 格式设置
    scale_fill_manual(
      values = c("0" = "skyblue", "1" = "salmon"),
      labels = c("0" = "True Regime 1", "1" = "True Regime 2"),
      name = "Background (Ground Truth)"
    ) +
    scale_color_manual(
      values = c("Viterbi Path" = "darkblue"),
      name = "Inference"
    ) +
    scale_linetype_manual(
      values = c("True Path" = "dotted"),
      name = "Reference"
    ) +
    scale_y_continuous(
      breaks = c(0, 1),
      labels = c("Regime 1", "Regime 2"),
      limits = c(-0.1, 1.1)
    ) +
    labs(
      title = "Regime Detection: Viterbi Global Decoding with O-M Functional",
      subtitle = sprintf("Lambda_OM = %.2f | Normalization = %s", lambda_om, normalize_method),
      x = "Time Step",
      y = "Decoded State"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(list(
    plot = p,
    viterbi_path = viterbi_path,
    accuracy = mean(viterbi_path == as.numeric(as.character(state_data$True_State)))
  ))
}
