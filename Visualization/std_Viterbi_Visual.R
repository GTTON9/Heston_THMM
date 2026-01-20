plot_viterbi_triple <- function(V_daily, RV_V, RV_smooth, 
                                nstates, Gamma, kappa, theta, sigma, 
                                Reg_chain_year,
                                save_path = "viterbi_comparison.png",
                                width = 18, height = 6, dpi = 300) {
  
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(patchwork)
  
  plot_single_viterbi <- function(V_simulated, nstates, Gamma, kappa, theta, sigma, Reg_chain, subtitle_text) {
    

    prob <- forward_filter(V_simulated, nstates, Gamma, kappa, theta, sigma)
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
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), # 加大标题
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
  
  p1 <- plot_single_viterbi(V_daily, nstates, Gamma, kappa, theta, sigma, 
                            Reg_chain_year, "Real Variance")
  
  p2 <- plot_single_viterbi(RV_V, nstates, Gamma, kappa, theta, sigma, 
                            Reg_chain_year, "Realized Volatility")
  
  p3 <- plot_single_viterbi(RV_smooth, nstates, Gamma, kappa, theta, sigma, 
                            Reg_chain_year, "Smoothed RV")
  

  combined_plot <- p1 | p2 | p3
  combined_plot <- combined_plot + 
    plot_annotation(
      title = "Standard HMM: Estimated Regime Sequence",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    )
  

  if (!is.null(save_path)) {
    ggsave(save_path, combined_plot, width = width, height = height, dpi = dpi)
    cat(sprintf("Plot saved to: %s\n", save_path))
  }
  
  return(combined_plot)
}

