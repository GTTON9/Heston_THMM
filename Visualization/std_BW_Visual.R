plot_cir_confidence_improved <- function(true_vol, 
                                              param, 
                                              states_estimate, 
                                              Reg_chain,
                                              dt = 1/252, 
                                              interval_step = 10,
                                              subtitle_text = "") {
  
  library(ggplot2)
  library(dplyr)
  
  # -------------------------
  # 0) Align lengths
  # -------------------------
  min_length <- min(length(true_vol), length(states_estimate), length(Reg_chain))
  true_vol <- true_vol[1:min_length]
  states_estimate <- states_estimate[1:min_length]
  Reg_chain <- Reg_chain[1:min_length]
  time_points <- 1:min_length
  
  # -------------------------
  # 1) V grid for density calculation
  # -------------------------
  V_min <- max(0.0001, min(true_vol, na.rm = TRUE) * 0.1)
  V_max <- max(true_vol, na.rm = TRUE) * 3
  V_grid <- seq(V_min, V_max, length.out = 500)
  dV <- diff(V_grid)[1]
  
  # -------------------------
  # 2) Find regime change points (based on states_estimate)
  # -------------------------
  changepoints <- c(1, which(diff(states_estimate) != 0) + 1, length(states_estimate) + 1)
  
  all_ci_data <- data.frame()
  
  # 固定预测步长：interval_step
  k_fixed <- interval_step * dt
  
  # -------------------------
  # 3) Compute rolling CIs for each regime segment
  # -------------------------
  for (seg in 1:(length(changepoints) - 1)) {
    start_idx <- changepoints[seg]
    end_idx <- changepoints[seg + 1] - 1
    
    if (end_idx - start_idx < interval_step) next
    
    current_regime <- states_estimate[start_idx]
    kappa_reg <- param$kappa[current_regime]
    theta_reg <- param$theta[current_regime]
    sigma_reg <- param$sigma[current_regime]
    
    # 在该段内：每隔 interval_step 计算一次 CI
    t_seq <- seq(start_idx + interval_step, end_idx, by = interval_step)
    
    for (t_idx in t_seq) {
      # ✅ rolling reference: 用前 interval_step 的真实值作为条件
      v_ref <- true_vol[t_idx - interval_step]
      
      log_densities <- sapply(V_grid, function(v) {
        ln_d_Heston(
          V_t = v_ref, 
          V_t_plus_k = v, 
          k = k_fixed, 
          kappa = kappa_reg, 
          theta = theta_reg, 
          sigma = sigma_reg
        )
      })
      
      densities <- exp(log_densities)
      if (any(!is.finite(densities)) || all(densities == 0)) next
      
      cdf <- cumsum(densities * dV)
      if (!is.finite(max(cdf)) || max(cdf) <= 0) next
      cdf <- cdf / max(cdf)
      
      idx_lower <- which.min(abs(cdf - 0.025))
      idx_mean  <- which.min(abs(cdf - 0.5))
      idx_upper <- which.min(abs(cdf - 0.975))
      
      lower_val <- V_grid[idx_lower]
      mean_val  <- V_grid[idx_mean]
      upper_val <- V_grid[idx_upper]
      
      if (is.finite(lower_val) && is.finite(upper_val) &&
          lower_val > 0 && upper_val > lower_val) {
        
        all_ci_data <- rbind(all_ci_data, data.frame(
          time = t_idx,
          lower = lower_val,
          mean  = mean_val,
          upper = upper_val,
          regime = current_regime
        ))
      }
    }
  }
  
  if (nrow(all_ci_data) == 0) {
    warning("No confidence intervals computed")
    return(NULL)
  }
  
  # -------------------------
  # 4) Prepare data
  # -------------------------
  vol_data <- data.frame(
    time = time_points, 
    variance = true_vol,
    regime_estimate = states_estimate
  )
  
  state_data <- data.frame(
    Time = time_points,
    Estimated_State = states_estimate - 1,
    True_State = as.factor(Reg_chain)
  )
  
  # Background rectangles for true regime
  max_time <- max(state_data$Time)
  regime_background <- state_data %>%
    mutate(change = True_State != lag(True_State, default = first(True_State))) %>%
    filter(change | row_number() == 1) %>%
    dplyr::select(Time, True_State) %>%
    mutate(Time_End = lead(Time, default = max_time + 1)) %>%
    rename(Time_Start = Time)
  
  # y-axis limits
  y_min <- min(c(vol_data$variance, all_ci_data$lower), na.rm = TRUE) * 0.95
  y_max <- max(c(vol_data$variance, all_ci_data$upper), na.rm = TRUE) * 1.05
  
  # accuracy
  accuracy <- mean(state_data$Estimated_State == as.numeric(as.character(state_data$True_State))) * 100
  
  # Create segments for variance line colored by regime estimate
  vol_segments <- data.frame()
  for (i in 1:(nrow(vol_data) - 1)) {
    vol_segments <- rbind(vol_segments, data.frame(
      x = vol_data$time[i],
      xend = vol_data$time[i + 1],
      y = vol_data$variance[i],
      yend = vol_data$variance[i + 1],
      regime = vol_data$regime_estimate[i]
    ))
  }
  
  # -------------------------
  # 5) Plot
  # -------------------------
  p <- ggplot() +
    geom_rect(
      data = regime_background,
      aes(xmin = Time_Start, xmax = Time_End, ymin = y_min, ymax = y_max, fill = True_State),
      alpha = 0.15,
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = vol_segments,
      aes(x = x, xend = xend, y = y, yend = yend, color = factor(regime)),
      linewidth = 0.8, alpha = 0.9
    ) +
    geom_segment(
      data = all_ci_data,
      aes(x = time, xend = time, y = lower, yend = upper, linetype = "95% CI"),
      color = "gray40",
      alpha = 0.5, linewidth = 0.8
    ) +
    geom_point(
      data = all_ci_data,
      aes(x = time, y = lower),
      color = "gray30",
      size = 0.8, alpha = 0.6
    ) +
    geom_point(
      data = all_ci_data,
      aes(x = time, y = upper),
      color = "gray30",
      size = 0.8, alpha = 0.6
    ) +
    scale_fill_manual(
      values = c("0" = "skyblue", "1" = "salmon"),
      labels = c("0" = "Regime 1\n(Calm)", "1" = "Regime 2\n(Turbulent)"),
      name = "True Regime"
    ) +
    scale_color_manual(
      values = c("1" = "blue", "2" = "red"),
      labels = c("1" = "Regime 1", "2" = "Regime 2"),
      name = "Variance Process\n(Estimated Regime)"
    ) +
    scale_linetype_manual(values = c("95% CI" = "solid"), name = "") +
    labs(
      title = sprintf("%s | Accuracy: %.2f%%", subtitle_text, accuracy),
      x = "Time Step",
      y = "Variance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), # 加大标题
      plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold"),
      
      # 核心修改部分：
      legend.position = "bottom",          # 图例置于底部
      legend.box = "vertical",             # 强制图例组垂直堆叠，确保分成三排
      legend.margin = margin(t = 10),      # 增加图例与主图的间距
      legend.spacing.y = unit(0.2, "cm"),  # 增加排与排之间的垂直间距
      
      legend.text = element_text(size = 12),    # 加大图例文字
      legend.title = element_text(size = 11, face = "bold"), # 加大图例标题
      legend.key.size = unit(1.2, "lines")      # 稍微加大图例图标大小，防止文字挤压
    ) +
    guides(
      # 通过 ncol = 1 确保每个图例内部是纵向的，
      # 配合 theme 里的 legend.box = "vertical" 即可实现三排布局
      fill = guide_legend(order = 1, ncol = 2),   # True Regime 如果选项多可以用 ncol=2
      color = guide_legend(order = 2, ncol = 1),  
      shape = guide_legend(order = 3, ncol = 1)
    )
  return(list(plot = p, ci_data = all_ci_data))
}




# ===== 更新的 Baum-Welch Triple Analysis =====
baum_welch_triple_analysis <- function(V_daily, RV_V, RV_smooth, 
                                       Reg_chain_year, 
                                       date_sequence,
                                       nstates = 2,
                                       runs_V = 10,
                                       runs_RV = 10,
                                       runs_RV_smooth = 10,
                                       dt = 1/252,
                                       interval_step = 10,
                                       save_viterbi_path = "baum_welch_viterbi.png",
                                       save_ci_path = "baum_welch_confidence.png",
                                       width = 18, height = 6, dpi = 300) {
  
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(patchwork)
  
  cat("=================================================================\n")
  cat("Baum-Welch Triple Analysis\n")
  cat("=================================================================\n\n")
  
  BW_all_param <- list()
  all_states_estimate <- list()
  
  # ===== 1. Real Variance (V_daily) =====
  cat("Processing Real Variance...\n")
  
  my_data_df <- data.frame(Date = date_sequence, Var = V_daily)
  series_control <- Heston_set_controls( 
    states = nstates, sdds = "Heston", date_column = "Date",
    file = my_data_df, data_column = "Var", logreturns = FALSE,
    from = date_sequence[1], to = date_sequence[length(date_sequence)],
    runs = runs_V
  )
  
  data_hmm <- prepare_data(series_control)
  model_hmm <- Heston_fit_model(data_hmm) 
  final_model <- decode_states_heston(model_hmm) 
  
  param_true_V <- parUncon2par_heston(final_model$estimate, series_control, 
                                      FALSE, numerical_safeguard = TRUE)
  BW_all_param[[1]] <- param_true_V
  all_states_estimate[[1]] <- final_model$decoding
  
  cat("  Real Variance - Done\n\n")
  
  # ===== 2. Realized Variance (RV_V) =====
  cat("Processing Realized Variance...\n")
  
  my_data_df <- data.frame(Date = date_sequence, Var = RV_V)
  series_control <- Heston_set_controls( 
    states = nstates, sdds = "Heston", date_column = "Date",
    file = my_data_df, data_column = "Var", logreturns = FALSE,
    from = date_sequence[1], to = date_sequence[length(date_sequence)],
    runs = runs_RV
  )
  
  data_hmm <- prepare_data(series_control)
  model_hmm <- Heston_fit_model(data_hmm) 
  final_model <- decode_states_heston(model_hmm) 
  
  param_RV <- parUncon2par_heston(final_model$estimate, series_control, 
                                  FALSE, numerical_safeguard = TRUE)
  BW_all_param[[2]] <- param_RV
  all_states_estimate[[2]] <- final_model$decoding
  
  cat("  Realized Variance - Done\n\n")
  
  # ===== 3. Smoothed RV =====
  cat("Processing Smoothed Realized Variance...\n")
  
  my_data_df <- data.frame(Date = date_sequence, Var = RV_smooth)
  series_control <- Heston_set_controls( 
    states = nstates, sdds = "Heston", date_column = "Date",
    file = my_data_df, data_column = "Var", logreturns = FALSE,
    from = date_sequence[1], to = date_sequence[length(date_sequence)],
    runs = runs_RV_smooth
  )
  
  data_hmm <- prepare_data(series_control)
  model_hmm <- Heston_fit_model(data_hmm) 
  final_model <- decode_states_heston(model_hmm) 
  
  param_RV_smooth <- parUncon2par_heston(final_model$estimate, series_control, 
                                         FALSE, numerical_safeguard = TRUE)
  BW_all_param[[3]] <- param_RV_smooth
  all_states_estimate[[3]] <- final_model$decoding
  
  cat("  Smoothed RV - Done\n\n")
  
  # ===== Internal function for single Viterbi plot =====
  plot_single_viterbi <- function(V_simulated, param, Reg_chain, subtitle_text) {
    
    prob <- forward_filter(V_simulated, nstates, param$Gamma, 
                          param$kappa, param$theta, param$sigma)
    states_estimate <- prob$states_estimate
    a <- prob$posterior
    
    time_index <- 1:ncol(a)
    
    prob_data <- data.frame(
      Time = time_index,
      Prob_State1 = a[1,],
      Prob_State2 = 1 - a[1,]
    )
    
    state_data <- data.frame(
      Time = time_index,
      Estimated_State = states_estimate - 1, 
      True_State = as.factor(Reg_chain[1:ncol(a)])
    )
    
    prob_data_long <- pivot_longer(
      prob_data,
      cols = starts_with("Prob_"),
      names_to = "Probability_Type",
      values_to = "Probability"
    )
    
    max_time <- max(state_data$Time)
    regime_start_times <- state_data %>%
      mutate(change = True_State != lag(True_State, default = first(True_State))) %>%
      filter(change | row_number() == 1) %>%
      dplyr::select(Time, True_State) %>%
      mutate(Time_End = lead(Time, default = max_time + 1)) %>%
      rename(Time_Start = Time)
    
    accuracy <- mean(state_data$Estimated_State == as.numeric(as.character(state_data$True_State))) * 100
    
    plot <- ggplot() +
      geom_rect(
        data = regime_start_times,
        aes(xmin = Time_Start, xmax = Time_End, ymin = -0.05, ymax = 1.05, fill = True_State),
        alpha = 0.2, 
        inherit.aes = FALSE
      ) +
      geom_line(
        data = prob_data_long, 
        aes(x = Time, y = Probability, color = Probability_Type),
        linewidth = 0.9
      ) +
      geom_point(
        data = state_data,
        aes(x = Time, y = Estimated_State, shape = "Estimated State"), 
        size = 1,
        color = "black",
        alpha = 0.8
      ) +
      scale_fill_manual(
        values = c("0" = "skyblue", "1" = "salmon"), 
        labels = c("0" = "Regime 1\n(Calm)", "1" = "Regime 2\n(Turbulent)"),
        name = "True Regime"
      ) +
      scale_color_manual(
        values = c("Prob_State1" = "blue", "Prob_State2" = "red"),
        labels = c("Prob_State1" = expression(paste(italic(P), "(", bold(V)[paste("1:", italic(t))], ", ", italic(R)[italic(t)], " = 1 | ", bold(theta), ")")),
                   "Prob_State2" = expression(paste(italic(P), "(", bold(V)[paste("1:", italic(t))], ", ", italic(R)[italic(t)], " = 2 | ", bold(theta), ")"))),
        name = "Filtered\nProbability"
      ) +
      scale_shape_manual(
        values = c("Estimated State" = 1),
        labels = c("Estimated State" = ""),
        name = "Regime Estimate"
      ) +
      scale_y_continuous(
        breaks = c(0, 0.5, 1),
        labels = c("0", "0.5", "1"),
        name = "Probability / State"
      ) +
      labs(
        title = sprintf("Accuracy: %.2f%%", accuracy),
        subtitle = subtitle_text,
        x = "Time Step"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12, face = "bold"),
      
        legend.position = "bottom",          
        legend.box = "vertical",             
        legend.margin = margin(t = 10),     
        legend.spacing.y = unit(0.2, "cm"),  
        
        legend.text = element_text(size = 12),   
        legend.title = element_text(size = 11, face = "bold"),
        legend.key.size = unit(1.2, "lines")     
      ) +
      guides(
        fill = guide_legend(order = 1, ncol = 2),  
        color = guide_legend(order = 2, ncol = 1),  
        shape = guide_legend(order = 3, ncol = 1)
      )
    
    return(plot)
  }
  
  # ===== 4. Create Viterbi Plots =====
  cat("Creating Viterbi plots...\n")
  
  p1_vit <- plot_single_viterbi(V_daily, param_true_V, Reg_chain_year, "Real Variance")
  p2_vit <- plot_single_viterbi(RV_V, param_RV, Reg_chain_year, "Realized Variance")
  p3_vit <- plot_single_viterbi(RV_smooth, param_RV_smooth, Reg_chain_year, "Smoothed RV")
  
  combined_viterbi_plot <- p1_vit | p2_vit | p3_vit
  combined_viterbi_plot <- combined_viterbi_plot + 
    plot_annotation(
      title = "Baum-Welch HMM: Estimated Regime Sequence",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    )
  
  if (!is.null(save_viterbi_path)) {
    ggsave(save_viterbi_path, combined_viterbi_plot, width = width, height = height, dpi = dpi)
    cat(sprintf("Viterbi plot saved to: %s\n", save_viterbi_path))
  }
  
  # ===== 5. Create Confidence Interval Plots =====
  cat("Creating confidence interval plots...\n")
  
  result1 <- plot_cir_confidence_improved(
    true_vol = V_daily,
    param = param_true_V,
    states_estimate = all_states_estimate[[1]],
    Reg_chain = Reg_chain_year,
    dt = dt,
    interval_step = interval_step,
    subtitle_text = "Real Variance"
  )
  
  result2 <- plot_cir_confidence_improved(
    true_vol = V_daily,
    param = param_RV,
    states_estimate = all_states_estimate[[2]],
    Reg_chain = Reg_chain_year,
    dt = dt,
    interval_step = interval_step,
    subtitle_text = "Realized Variance"
  )
  
  result3 <- plot_cir_confidence_improved(
    true_vol = V_daily,
    param = param_RV_smooth,
    states_estimate = all_states_estimate[[3]],
    Reg_chain = Reg_chain_year,
    dt = dt,
    interval_step = interval_step,
    subtitle_text = "Smoothed RV"
  )
  
  # Combine CI plots
  combined_ci_plot <- result1$plot | result2$plot | result3$plot
  combined_ci_plot <- combined_ci_plot + 
    plot_annotation(
      title = "CIR Model Confidence Intervals: Baum-Welch Estimates",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    )
  
  if (!is.null(save_ci_path)) {
    ggsave(save_ci_path, combined_ci_plot, width = width, height = height, dpi = dpi)
    cat(sprintf("CI plot saved to: %s\n", save_ci_path))
  }
  
  cat("\n=================================================================\n")
  cat("Analysis Complete!\n")
  cat("=================================================================\n")
  
  return(list(
    params = BW_all_param,
    states_estimate = all_states_estimate,
    viterbi_plot = combined_viterbi_plot,
    ci_plot = combined_ci_plot,
    ci_data = list(
      real_var = result1$ci_data,
      rv = result2$ci_data,
      rv_smooth = result3$ci_data
    )
  ))
}

