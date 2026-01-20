plot_viterbi_om <- function(batch_results_complete,
                            nstates, 
                            Gamma, 
                            kappa, 
                            theta, 
                            sigma,
                            Reg_chain,
                            lambda_om = 1.0,
                            use_upper = TRUE,
                            use_lower = TRUE,
                            normalize_method = "log",
                            show_om_contribution = FALSE) {
  
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  
  # ===== Step 1: Prepare O-M matrices =====
  cat("Preparing O-M matrices...\n")
  
  # Normalize O-M values if needed
  if (normalize_method == "log") {
    OM_upper <- log(batch_results_complete[, c("OM_upper_R1", "OM_upper_R2")] + 1)
    OM_lower <- log(batch_results_complete[, c("OM_lower_R1", "OM_lower_R2")] + 1)
  } else if (normalize_method == "zscore") {
    OM_all <- c(
      batch_results_complete$OM_upper_R1,
      batch_results_complete$OM_upper_R2,
      batch_results_complete$OM_lower_R1,
      batch_results_complete$OM_lower_R2
    )
    OM_finite <- OM_all[is.finite(OM_all)]
    mean_om <- mean(OM_finite)
    sd_om <- sd(OM_finite)
    
    OM_upper <- (batch_results_complete[, c("OM_upper_R1", "OM_upper_R2")] - mean_om) / sd_om
    OM_lower <- (batch_results_complete[, c("OM_lower_R1", "OM_lower_R2")] - mean_om) / sd_om
  } else {
    # Raw values
    OM_upper <- batch_results_complete[, c("OM_upper_R1", "OM_upper_R2")]
    OM_lower <- batch_results_complete[, c("OM_lower_R1", "OM_lower_R2")]
  }
  
  OM_upper <- as.matrix(OM_upper)
  OM_lower <- as.matrix(OM_lower)
  
  # ===== Step 2: Run Forward Filter with O-M =====
  cat("Running Forward Filter with O-M...\n")
  
  observations <- batch_results_complete$v
  prob <-forward_filter_om_pure(OM_upper,OM_lower,
                                     nstates,
                                     Gamma,
                                     lambda_om = 1.0,
                                     use_upper = TRUE,
                                     use_lower = TRUE)

  
  # prob <- forward_filter_om(
  #   observations = observations,
  #   OM_upper = OM_upper,
  #   OM_lower = OM_lower,
  #   nstates = nstates,
  #   Gamma = Gamma,
  #   kappa = kappa,
  #   theta = theta,
  #   sigma = sigma,
  #   lambda_om = lambda_om,
  #   use_upper = use_upper,
  #   use_lower = use_lower
  # )
  
  states_estimate <- prob$states_estimate
  posterior <- prob$posterior
  print(posterior[1,])
  # ===== Step 3: Prepare data for plotting =====
  
  time_index <- 1:ncol(posterior)
  
  # Probability data
  prob_data <- data.frame(
    Time = time_index,
    Prob_State1 = posterior[1, ],
    Prob_State2 = posterior[2, ]
  )
  
  prob_data_long <- pivot_longer(
    prob_data,
    cols = starts_with("Prob_"),
    names_to = "Probability_Type",
    values_to = "Probability"
  )
  
  # State data
  state_data <- data.frame(
    Time = time_index,
    Estimated_State = states_estimate - 1,
    True_State = as.factor(Reg_chain[1:ncol(posterior)])
  )
  
  # Regime background rectangles
  max_time <- max(state_data$Time)
  regime_start_times <- state_data %>%
    mutate(change = True_State != lag(True_State, default = first(True_State))) %>%
    filter(change | row_number() == 1) %>%
    dplyr::select(Time, True_State) %>%
    mutate(Time_End = lead(Time, default = max_time + 1)) %>%
    rename(Time_Start = Time)
  
  # ===== Step 4: Create main plot =====
  
  main_plot <- ggplot() +
    # Background shading for true regimes
    geom_rect(
      data = regime_start_times,
      aes(
        xmin = Time_Start,
        xmax = Time_End,
        ymin = -0.05,
        ymax = 1.05,
        fill = True_State
      ),
      alpha = 0.2,
      inherit.aes = FALSE
    ) +
    
    # Posterior probability lines
    geom_line(
      data = prob_data_long,
      aes(x = Time, y = Probability, color = Probability_Type),
      linewidth = 0.9
    ) +
    
    # Estimated states
    geom_point(
      data = state_data,
      aes(x = Time, y = Estimated_State, shape = "Estimated State"),
      size = 1,
      color = "black",
      alpha = 0.8
    ) +
    
    # Scales
    scale_fill_manual(
      values = c("0" = "skyblue", "1" = "salmon"),
      labels = c("0" = "Regime 1 (Calm)", "1" = "Regime 2 (Turbulent)"),
      name = "True Regime Background"
    ) +
    scale_color_manual(
      values = c("Prob_State1" = "blue", "Prob_State2" = "red"),
      labels = c(
        "Prob_State1" = expression(P(s[t] == 1 * " | " * bold(y)[1:t], bold(phi)[1:t])),
        "Prob_State2" = expression(P(s[t] == 2 * " | " * bold(y)[1:t], bold(phi)[1:t]))
      ),
      name = "Forward Probability (O-M Enhanced)"
    ) +
    scale_shape_manual(
      values = c("Estimated State" = 1),
      name = "Regime Estimate"
    ) +
    scale_y_continuous(
      breaks = c(0, 0.5, 1),
      labels = c("0", "0.5", "1"),
      name = "Posterior Probability / State"
    ) +
    labs(
      title = sprintf("HMM with O-M Functional (norm=%s)", 
                      normalize_method),
      subtitle = sprintf("Global Optimal Path | Accuracy: %.2f%%", 
                         mean(state_data$Estimated_State == as.numeric(as.character(state_data$True_State))) * 100),
      x = "Time Step"
    ) +
    theme_minimal() +
    guides(
      fill = guide_legend(order = 1),
      color = guide_legend(order = 2),
      shape = guide_legend(order = 3)
    )
  
  # ===== Step 5: Optional O-M contribution plot (PROBABILITY SCALE) =====
  
  if (show_om_contribution) {
    
    # Calculate O-M contribution as PROBABILITY (not log)
    # Emission probability: exp(-λ × L[φ])
    
    # Total O-M functional per regime
    L_total_R1 <- OM_upper[, 1] + OM_lower[, 1]
    L_total_R2 <- OM_upper[, 2] + OM_lower[, 2]
    
    # Convert to probability: P ∝ exp(-λ × L)
    prob_R1 <- exp(-lambda_om * L_total_R1)
    prob_R2 <- exp(-lambda_om * L_total_R2)
    
    # Normalize to make probabilities comparable
    # (Each regime's probability at each time point)
    total_prob <- prob_R1 + prob_R2
    prob_R1_norm <- prob_R1 / total_prob
    prob_R2_norm <- prob_R2 / total_prob
    
    om_contrib_data <- data.frame(
      Time = 1:nrow(batch_results_complete),
      Prob_R1 = prob_R1_norm,
      Prob_R2 = prob_R2_norm
    )
    
    om_contrib_long <- pivot_longer(
      om_contrib_data,
      cols = starts_with("Prob_"),
      names_to = "Regime",
      values_to = "OM_Probability"
    )
    
    om_plot <- ggplot(om_contrib_long, aes(x = Time, y = OM_Probability, color = Regime)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50", alpha = 0.5) +
      scale_color_manual(
        values = c("Prob_R1" = "blue", "Prob_R2" = "red"),
        labels = c("Prob_R1" = "Regime 1 (Calm)", "Prob_R2" = "Regime 2 (Turbulent)"),
        name = "Regime"
      ) +
      scale_y_continuous(
        limits = c(0, 1),
        breaks = c(0, 0.25, 0.5, 0.75, 1),
        labels = c("0", "0.25", "0.5", "0.75", "1")
      ) +
      labs(
        title = "O-M Based Path Probability (Normalized)",
        subtitle = expression(
          "P(path | regime)"
        ),
        x = "Time Step",
        y = "Probability"
      ) +
      theme_minimal() +
      theme(
        plot.subtitle = element_text(size = 9, color = "gray30")
      )
    
    # Combine plots
    library(gridExtra)
    combined_plot <- grid.arrange(main_plot, om_plot, ncol = 1, heights = c(2, 1))
    
    return(list(
      plot = combined_plot,
      main_plot = main_plot,
      om_plot = om_plot,
      states_estimate = states_estimate,
      posterior = posterior,
      om_probabilities = om_contrib_data  # 返回概率数据
    ))
    
  } else {
    return(list(
      plot = main_plot,
      states_estimate = states_estimate,
      posterior = posterior
    ))
  }
}