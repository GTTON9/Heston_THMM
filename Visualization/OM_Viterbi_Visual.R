plot_viterbi_om_triple <- function(V_daily, RV_V, RV_smooth,
                                   data_with_IV,
                                   regime_params,
                                   nstates, 
                                   Gamma, 
                                   Reg_chain_year,
                                   H = 5/252,
                                   lambda_om = 1.0,
                                   normalize_method = "log",
                                   save_path = "om_viterbi_comparison.png",
                                   width = 18, height = 6, dpi = 300) {
  
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(patchwork)
  
  cat("=================================================================\n")
  cat("OM Viterbi Triple Analysis\n")
  cat("=================================================================\n\n")
  
  # 内部函数：单个 OM Viterbi 绘图
  plot_single_om_viterbi <- function(variance_data, subtitle_text) {
    
    cat(sprintf("Processing %s...\n", subtitle_text))
    
    # Step 1: Calculate expected moves
    batch_results <- calculate_expected_moves_batch(data_with_IV$S, variance_data)
    
    # Step 2: Construct linear paths
    result_with_paths <- construct_linear_paths(batch_results, H, data_with_IV)
    
    # Step 3: Calculate OM functional
    result_complete <- calculate_OM_batch(result_with_paths, regime_params)
    
    # Step 4: Prepare OM matrices
    if (normalize_method == "log") {
      OM_upper <- log(result_complete[, c("OM_upper_R1", "OM_upper_R2")] + 1)
      OM_lower <- log(result_complete[, c("OM_lower_R1", "OM_lower_R2")] + 1)
    } else if (normalize_method == "zscore") {
      OM_all <- c(result_complete$OM_upper_R1, result_complete$OM_upper_R2,
                  result_complete$OM_lower_R1, result_complete$OM_lower_R2)
      OM_finite <- OM_all[is.finite(OM_all)]
      mean_om <- mean(OM_finite)
      sd_om <- sd(OM_finite)
      OM_upper <- (result_complete[, c("OM_upper_R1", "OM_upper_R2")] - mean_om) / sd_om
      OM_lower <- (result_complete[, c("OM_lower_R1", "OM_lower_R2")] - mean_om) / sd_om
    } else {
      OM_upper <- result_complete[, c("OM_upper_R1", "OM_upper_R2")]
      OM_lower <- result_complete[, c("OM_lower_R1", "OM_lower_R2")]
    }
    
    OM_upper <- as.matrix(OM_upper)
    OM_lower <- as.matrix(OM_lower)
    

    
    # Step 6: Run forward filter for posterior probabilities
    forward_results <- forward_filter_om_pure(
      OM_upper = OM_upper,
      OM_lower = OM_lower,
      nstates = nstates,
      Gamma = Gamma,
      lambda_om = lambda_om,
      use_upper = TRUE,
      use_lower = TRUE
    )
    
    # Prepare data
    viterbi_path <- forward_results$states_estimate - 1
    posterior <- forward_results$posterior
    time_index <- 1:length(viterbi_path)
    
    # Probability data
    prob_data <- data.frame(
      Time = time_index,
      Prob_State1 = posterior[1,],
      Prob_State2 = posterior[2,]
    )
    
    # State data
    state_data <- data.frame(
      Time = time_index,
      Estimated_State = viterbi_path,
      True_State = as.factor(Reg_chain_year[1:length(viterbi_path)])
    )
    
    # Long format probability data
    prob_data_long <- pivot_longer(
      prob_data,
      cols = starts_with("Prob_"),
      names_to = "Probability_Type",
      values_to = "Probability"
    )
    
    # Regime background rectangles
    max_time <- max(state_data$Time)
    regime_start_times <- state_data %>%
      mutate(change = True_State != lag(True_State, default = first(True_State))) %>%
      filter(change | row_number() == 1) %>%
      dplyr::select(Time, True_State) %>%
      mutate(Time_End = lead(Time, default = max_time + 1)) %>%
      rename(Time_Start = Time)
    
    # Calculate accuracy
    accuracy <- mean(state_data$Estimated_State == as.numeric(as.character(state_data$True_State))) * 100
    
    # Plot
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
        labels = c("Prob_State1" = expression(paste(italic(P), "(", phi[paste("1:", italic(t))], ", ", italic(R)[italic(t)], " = 1 | ", bold(theta), ")")),
                   "Prob_State2" = expression(paste(italic(P), "(", phi[paste("1:", italic(t))], ", ", italic(R)[italic(t)], " = 2 | ", bold(theta), ")"))),
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
    
    cat(sprintf("  %s - Done (Accuracy: %.2f%%)\n\n", subtitle_text, accuracy))
    
    return(list(
      plot = plot,
      viterbi_path = viterbi_path,
      accuracy = accuracy
    ))
  }
  
  # Generate three plots
  result1 <- plot_single_om_viterbi(V_daily, "Real Variance")
  result2 <- plot_single_om_viterbi(RV_V, "Realized Volatility")
  result3 <- plot_single_om_viterbi(RV_smooth, "Smoothed RV")
  
  # Combine plots
  combined_plot <- result1$plot | result2$plot | result3$plot
  combined_plot <- combined_plot +
    plot_annotation(
      title = "OM-HMM: Estimated Regime Sequence with O-M Functional",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    )
  
  # Save
  if (!is.null(save_path)) {
    ggsave(save_path, combined_plot, width = width, height = height, dpi = dpi)
    cat(sprintf("Plot saved to: %s\n", save_path))
  }
  
  cat("=================================================================\n")
  cat("OM Viterbi Triple Analysis Complete!\n")
  cat(sprintf("Average Accuracy: %.2f%%\n", 
              mean(c(result1$accuracy, result2$accuracy, result3$accuracy))))
  cat("=================================================================\n\n")
  
  return(list(
    combined_plot = combined_plot,
    results = list(
      real_var = result1,
      rv = result2,
      rv_smooth = result3
    )
  ))
}
