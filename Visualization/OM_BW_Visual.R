plot_cir_confidence_om <- function(true_vol,
                                        param_om,
                                        states_estimate_om,
                                        Reg_chain,
                                        dt = 1/252,
                                        interval_step = 10,
                                        subtitle_text = "",
                                        V_grid_length = 500,
                                        V_min_floor = 1e-4,
                                        V_min_scale = 0.1,
                                        V_max_scale = 3,
                                        true_fill = c("0"="skyblue","1"="salmon"),
                                        est_color = c("1"="blue","2"="red")) {
  library(ggplot2)
  library(dplyr)
  
  # ---- align length
  min_length <- min(length(true_vol), length(states_estimate_om), length(Reg_chain))
  true_vol <- true_vol[1:min_length]
  states_estimate_om <- states_estimate_om[1:min_length]
  Reg_chain <- Reg_chain[1:min_length]
  time_points <- 1:min_length
  
  # ---- states to 1/2 indexing
  uniq_states <- sort(unique(states_estimate_om))
  if (all(uniq_states %in% c(0,1))) states_idx <- states_estimate_om + 1
  else if (all(uniq_states %in% c(1,2))) states_idx <- states_estimate_om
  else stop("states_estimate_om must be {0,1} or {1,2}")
  
  true_state_01 <- as.numeric(Reg_chain)
  est_state_01  <- states_idx - 1
  accuracy <- mean(est_state_01 == true_state_01) * 100
  
  # ---- grid
  V_min <- max(V_min_floor, min(true_vol, na.rm=TRUE) * V_min_scale)
  V_max <- max(true_vol, na.rm=TRUE) * V_max_scale
  V_grid <- seq(V_min, V_max, length.out = V_grid_length)
  dV <- diff(V_grid)[1]
  
  # ---- segment by OM regimes
  changepoints <- c(1, which(diff(states_idx) != 0) + 1, length(states_idx) + 1)
  
  all_ci_data <- data.frame()
  
  # ---- rolling CI: fixed horizon = interval_step*dt, reference updates each time
  k_fixed <- interval_step * dt
  
  for (seg in 1:(length(changepoints)-1)) {
    start_idx <- changepoints[seg]
    end_idx <- changepoints[seg+1] - 1
    
    if (end_idx - start_idx < interval_step) next
    
    current_regime <- states_idx[start_idx]  # 1/2
    kappa_reg <- param_om$kappa[current_regime]
    theta_reg <- param_om$theta[current_regime]
    sigma_reg <- param_om$sigma[current_regime]
    
    # t_idx runs inside this segment
    t_seq <- seq(start_idx + interval_step, end_idx, by = interval_step)
    
    for (t_idx in t_seq) {
      v_ref <- true_vol[t_idx - interval_step]   # ✅ rolling reference
      
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
      
      dens <- exp(log_densities)
      if (any(!is.finite(dens)) || all(dens == 0)) next
      
      cdf <- cumsum(dens * dV)
      if (!is.finite(max(cdf)) || max(cdf) <= 0) next
      cdf <- cdf / max(cdf)
      
      idx_lower <- which.min(abs(cdf - 0.025))
      idx_mean  <- which.min(abs(cdf - 0.5))
      idx_upper <- which.min(abs(cdf - 0.975))
      
      all_ci_data <- rbind(all_ci_data, data.frame(
        time = t_idx,
        lower = V_grid[idx_lower],
        mean  = V_grid[idx_mean],
        upper = V_grid[idx_upper],
        regime = current_regime
      ))
    }
  }
  
  if (nrow(all_ci_data) == 0) {
    warning("No confidence intervals computed (OM roll version).")
    return(NULL)
  }
  
  # ---- plot prep
  vol_data <- data.frame(time=time_points, variance=true_vol, regime_estimate=states_idx)
  
  state_data <- data.frame(
    Time=time_points,
    Estimated_State=states_idx - 1,
    True_State=as.factor(true_state_01)
  )
  
  max_time <- max(state_data$Time)
  regime_background <- state_data %>%
    mutate(change = True_State != lag(True_State, default = first(True_State))) %>%
    filter(change | row_number()==1) %>%
    select(Time, True_State) %>%
    mutate(Time_End = lead(Time, default = max_time+1)) %>%
    rename(Time_Start=Time)
  
  y_min <- min(c(vol_data$variance, all_ci_data$lower), na.rm=TRUE) * 0.95
  y_max <- max(c(vol_data$variance, all_ci_data$upper), na.rm=TRUE) * 1.05
  
  vol_segments <- do.call(rbind, lapply(1:(nrow(vol_data)-1), function(i){
    data.frame(x=vol_data$time[i], xend=vol_data$time[i+1],
               y=vol_data$variance[i], yend=vol_data$variance[i+1],
               regime=vol_data$regime_estimate[i])
  }))
  
  p <- ggplot() +
    geom_rect(data=regime_background,
              aes(xmin=Time_Start, xmax=Time_End, ymin=y_min, ymax=y_max, fill=True_State),
              alpha=0.15, inherit.aes=FALSE) +
    geom_segment(data=vol_segments,
                 aes(x=x, xend=xend, y=y, yend=yend, color=factor(regime)),
                 linewidth=0.8, alpha=0.9) +
    geom_segment(data=all_ci_data,
                 aes(x=time, xend=time, y=lower, yend=upper, linetype="95% CI"),
                 color="gray40", alpha=0.5, linewidth=0.8) +
    geom_point(data=all_ci_data, aes(x=time, y=lower),
               color="gray30", size=0.8, alpha=0.6) +
    geom_point(data=all_ci_data, aes(x=time, y=upper),
               color="gray30", size=0.8, alpha=0.6) +
    scale_fill_manual(values=true_fill,
                      labels=c("0"="Regime 1\n(Calm)","1"="Regime 2\n(Turbulent)"),
                      name="True Regime") +
    scale_color_manual(values=est_color,
                       labels=c("1"="Regime 1","2"="Regime 2"),
                       name="Variance Process\n(OM Estimated Regime)") +
    scale_linetype_manual(values=c("95% CI"="solid"), name="") +
    labs(title=sprintf("%s | OM Accuracy: %.2f%%", subtitle_text, accuracy),
         x="Time Step", y="Variance") +
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
  
  list(plot=p, ci_data=all_ci_data, accuracy=accuracy)
}






plot_viterbi_om <- function(batch_results_complete,
                            nstates,
                            Gamma,
                            Reg_chain,
                            lambda_om = 1.0,
                            normalize_method = c("log", "zscore", "none"),
                            use_upper = TRUE,
                            use_lower = TRUE,
                            show_om_contribution = FALSE,
                            subtitle_text = NULL) {
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  normalize_method <- match.arg(normalize_method)
  
  # ------------------------------------------------------------
  # 1) Extract OM columns
  # ------------------------------------------------------------
  stopifnot(all(c("OM_upper_R1", "OM_upper_R2", "OM_lower_R1", "OM_lower_R2") %in%
                  colnames(batch_results_complete)))
  
  OM_upper_raw <- as.matrix(batch_results_complete[, c("OM_upper_R1", "OM_upper_R2")])
  OM_lower_raw <- as.matrix(batch_results_complete[, c("OM_lower_R1", "OM_lower_R2")])
  
  # ------------------------------------------------------------
  # 2) Normalize OM (optional)
  # ------------------------------------------------------------
  if (normalize_method == "log") {
    OM_upper <- log(OM_upper_raw + 1)
    OM_lower <- log(OM_lower_raw + 1)
  } else if (normalize_method == "zscore") {
    OM_all <- c(OM_upper_raw, OM_lower_raw)
    OM_finite <- OM_all[is.finite(OM_all)]
    m <- mean(OM_finite)
    s <- sd(OM_finite)
    if (is.finite(s) && s > 0) {
      OM_upper <- (OM_upper_raw - m) / s
      OM_lower <- (OM_lower_raw - m) / s
    } else {
      OM_upper <- OM_upper_raw
      OM_lower <- OM_lower_raw
    }
  } else {
    OM_upper <- OM_upper_raw
    OM_lower <- OM_lower_raw
  }
  
  OM_upper <- as.matrix(OM_upper)
  OM_lower <- as.matrix(OM_lower)
  
  # ------------------------------------------------------------
  # 3) Viterbi + Forward filter
  # ------------------------------------------------------------
  vit <- viterbi_om_pure(
    OM_upper = OM_upper,
    OM_lower = OM_lower,
    nstates = nstates,
    Gamma = Gamma,
    lambda_om = lambda_om,
    use_upper = use_upper,
    use_lower = use_lower
  )
  
  ff <- forward_filter_om_pure(
    OM_upper = OM_upper,
    OM_lower = OM_lower,
    nstates = nstates,
    Gamma = Gamma,
    lambda_om = lambda_om,
    use_upper = use_upper,
    use_lower = use_lower
  )
  
  states_estimate <- vit$states_estimate - 1   # assume {1,2} -> {0,1}
  posterior <- ff$posterior                    # K x T
  
  Tlen <- length(states_estimate)
  Reg_chain <- Reg_chain[1:Tlen]
  
  # ------------------------------------------------------------
  # 4) Prepare plot data
  # ------------------------------------------------------------
  prob_data <- data.frame(
    Time = 1:Tlen,
    Prob_State1 = posterior[1,],
    Prob_State2 = posterior[2,]
  )
  
  state_data <- data.frame(
    Time = 1:Tlen,
    Estimated_State = states_estimate,
    True_State = as.factor(Reg_chain)
  )
  
  prob_data_long <- pivot_longer(
    prob_data,
    cols = starts_with("Prob_"),
    names_to = "Probability_Type",
    values_to = "Probability"
  )
  
  max_time <- max(state_data$Time)
  regime_rects <- state_data %>%
    mutate(change = True_State != lag(True_State, default = first(True_State))) %>%
    filter(change | row_number() == 1) %>%
    dplyr::select(Time, True_State) %>%
    mutate(Time_End = lead(Time, default = max_time + 1)) %>%
    rename(Time_Start = Time)
  
  accuracy <- mean(
    state_data$Estimated_State == as.numeric(as.character(state_data$True_State)),
    na.rm = TRUE
  ) * 100
  
  # Optional OM contribution (difference between regimes)
  # This is just a diagnostic curve to show separability; scaled to [0,1] for plotting.
  if (show_om_contribution) {
    om_sep <- rowSums(cbind(
      OM_upper[,2] - OM_upper[,1],
      OM_lower[,2] - OM_lower[,1]
    ), na.rm = TRUE)
    om_sep_scaled <- (om_sep - min(om_sep, na.rm = TRUE)) /
      (max(om_sep, na.rm = TRUE) - min(om_sep, na.rm = TRUE) + 1e-12)
    om_df <- data.frame(Time = 1:Tlen, OM_sep = om_sep_scaled)
  }
  
  # ------------------------------------------------------------
  # 5) Plot
  # ------------------------------------------------------------
  if (is.null(subtitle_text)) subtitle_text <- "OM-Viterbi Decoding"
  
  p <- ggplot() +
    geom_rect(
      data = regime_rects,
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
      labels = c("Prob_State1" = "Filtered P(Regime=1)",
                 "Prob_State2" = "Filtered P(Regime=2)"),
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
  
  if (show_om_contribution) {
    p <- p + geom_line(
      data = om_df,
      aes(x = Time, y = OM_sep, linetype = "OM separability"),
      linewidth = 0.7,
      color = "black",
      alpha = 0.7,
      inherit.aes = FALSE
    ) +
      scale_linetype_manual(values = c("OM separability" = "solid"),
                            name = "",
                            labels = c("OM separability (scaled)"))
  }
  
  return(list(
    plot = p,
    states_estimate = states_estimate,
    posterior = posterior,
    accuracy = accuracy
  ))
}














# ===== Topological HMM Triple Analysis (OM Version) =====
om_triple_analysis <- function(S_daily,
                               data_with_IV,
                               V_daily, 
                               RV_V, 
                               RV_smooth,
                               Reg_chain_year,
                               date_sequence,
                               nstates = 2,
                               runs_trueV = 10,
                               runs_RV = 10,
                               runs_smooth = 10,
                               H = 5/252,
                               lambda_fit = 1.0,
                               lambda_viterbi = 1.0,
                               sigma_penalty_trueV = 0.0,
                               sigma_penalty_RV = 10.0,
                               sigma_penalty_smooth = 0.0,
                               use_upper = TRUE,
                               use_lower = TRUE,
                               use_normalization = TRUE,
                               N_integration = 100,
                               optimizer = "nlm",
                               normalize_method = "log",
                               dt = 1/252,
                               interval_step = 10,
                               save_viterbi_path = "om_viterbi.png",
                               save_ci_path = "om_confidence.png",
                               width = 18, 
                               height = 6, 
                               dpi = 300,
                               seed = 999,
                               verbose = TRUE) {
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  
  cat("=================================================================\n")
  cat("Topological HMM (OM Functional) Triple Analysis\n")
  cat("=================================================================\n\n")
  
  # ---------------------------------------------------------------
  # Internal: Process single volatility series
  # ---------------------------------------------------------------
  process_single <- function(vol_series, 
                             runs, 
                             sigma_penalty, 
                             label) {
    
    cat(sprintf("Processing: %s\n", label))
    
    # 1) Calculate expected moves
    batch_results <- calculate_expected_moves_batch(S_daily, vol_series)
    
    # 2) Construct linear paths
    result_with_paths <- construct_linear_paths(batch_results, H, data_with_IV)
    
    # 3) Set controls
    controls_om <- Heston_set_controls(
      states = nstates,
      sdds = "Heston",
      horizon = nrow(result_with_paths),
      runs = runs
    )
    
    # 4) Fit model using OM functional
    model_om <- Heston_fit_model_om(
      result_with_paths = result_with_paths,
      controls = controls_om,
      runs = runs,
      lambda_om = lambda_fit,
      sigma_penalty = sigma_penalty,
      use_normalization = use_normalization,
      H = H,
      N_integration = N_integration,
      optimizer = optimizer,
      seed = seed,
      verbose = verbose
    )
    
    # 5) Extract parameters
    param_om <- parUncon2par_heston(
      model_om$estimate,
      controls_om,
      FALSE,
      numerical_safeguard = TRUE
    )
    
    Gamma_est <- param_om$Gamma
    
    regime_params <- list(
      list(kappa = param_om$kappa[1], 
           theta = param_om$theta[1], 
           sigma = param_om$sigma[1]),
      list(kappa = param_om$kappa[2], 
           theta = param_om$theta[2], 
           sigma = param_om$sigma[2])
    )
    
    # 6) Calculate OM contributions
    result_complete <- calculate_OM_batch(result_with_paths, regime_params)
    
    # 7) Normalize OM
    OM_upper_raw <- as.matrix(result_complete[, c("OM_upper_R1", "OM_upper_R2")])
    OM_lower_raw <- as.matrix(result_complete[, c("OM_lower_R1", "OM_lower_R2")])
    
    if (normalize_method == "log") {
      OM_upper <- log(OM_upper_raw + 1)
      OM_lower <- log(OM_lower_raw + 1)
    } else if (normalize_method == "zscore") {
      OM_all <- c(OM_upper_raw, OM_lower_raw)
      OM_finite <- OM_all[is.finite(OM_all)]
      m <- mean(OM_finite)
      s <- sd(OM_finite)
      if (is.finite(s) && s > 0) {
        OM_upper <- (OM_upper_raw - m) / s
        OM_lower <- (OM_lower_raw - m) / s
      } else {
        OM_upper <- OM_upper_raw
        OM_lower <- OM_lower_raw
      }
    } else {
      OM_upper <- OM_upper_raw
      OM_lower <- OM_lower_raw
    }
    
    OM_upper <- as.matrix(OM_upper)
    OM_lower <- as.matrix(OM_lower)
    
    # 8) Viterbi decoding
    vit <- viterbi_om_pure(
      OM_upper = OM_upper,
      OM_lower = OM_lower,
      nstates = nstates,
      Gamma = Gamma_est,
      lambda_om = lambda_viterbi,
      use_upper = use_upper,
      use_lower = use_lower
    )
    
    # 9) Forward filter
    ff <- forward_filter_om_pure(
      OM_upper = OM_upper,
      OM_lower = OM_lower,
      nstates = nstates,
      Gamma = Gamma_est,
      lambda_om = lambda_viterbi,
      use_upper = use_upper,
      use_lower = use_lower
    )
    
    states_estimate <- vit$states_estimate
    posterior <- ff$posterior
    
    # 10) Plot Viterbi
    viterbi_plot <- plot_viterbi_single(
      states_estimate = states_estimate,
      posterior = posterior,
      Reg_chain = Reg_chain_year,
      subtitle_text = label
    )
    
    # 11) Plot CI
  
    ci_result <- plot_cir_confidence_om(
      true_vol = V_daily,
      param_om = param_om,
      states_estimate_om = states_estimate, # or OM_viterbi_result$states_estimate - 1
      Reg_chain = Reg_chain_year,
      dt = dt,
      interval_step = interval_step,
      subtitle_text = label
    )
    
    cat(sprintf("  %s - Done\n\n", label))
    
    return(list(
      label = label,
      param = param_om,
      Gamma = Gamma_est,
      states_estimate = states_estimate,
      posterior = posterior,
      viterbi_plot = viterbi_plot,
      ci_plot = ci_result$plot,
      ci_data = ci_result$ci_data
    ))
  }
  
  # ---------------------------------------------------------------
  # Internal: Plot single Viterbi (similar to standard version)
  # ---------------------------------------------------------------
  plot_viterbi_single <- function(states_estimate, 
                                  posterior, 
                                  Reg_chain, 
                                  subtitle_text) {
    
    states_estimate <- states_estimate - 1  # {1,2} -> {0,1}
    Tlen <- length(states_estimate)
    Reg_chain <- Reg_chain[1:Tlen]
    
    prob_data <- data.frame(
      Time = 1:Tlen,
      Prob_State1 = posterior[1,],
      Prob_State2 = posterior[2,]
    )
    
    state_data <- data.frame(
      Time = 1:Tlen,
      Estimated_State = states_estimate,
      True_State = as.factor(Reg_chain)
    )
    
    prob_data_long <- pivot_longer(
      prob_data,
      cols = starts_with("Prob_"),
      names_to = "Probability_Type",
      values_to = "Probability"
    )
    
    max_time <- max(state_data$Time)
    regime_rects <- state_data %>%
      mutate(change = True_State != lag(True_State, default = first(True_State))) %>%
      filter(change | row_number() == 1) %>%
      dplyr::select(Time, True_State) %>%
      mutate(Time_End = lead(Time, default = max_time + 1)) %>%
      rename(Time_Start = Time)
    
    accuracy <- mean(
      state_data$Estimated_State == as.numeric(as.character(state_data$True_State)),
      na.rm = TRUE
    ) * 100
    
    p <- ggplot() +
      geom_rect(
        data = regime_rects,
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
        labels = c("Prob_State1" = "Filtered P(Regime=1)",
                   "Prob_State2" = "Filtered P(Regime=2)"),
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
    
    return(p)
  }
  
  # ---------------------------------------------------------------
  # Process three scenarios
  # ---------------------------------------------------------------
  res_trueV <- process_single(
    vol_series = V_daily,
    runs = runs_trueV,
    sigma_penalty = sigma_penalty_trueV,
    label = "True Variance"
  )
  
  res_RV <- process_single(
    vol_series = RV_V,
    runs = runs_RV,
    sigma_penalty = sigma_penalty_RV,
    label = "Realized Variance"
  )
  
  res_smooth <- process_single(
    vol_series = RV_smooth,
    runs = runs_smooth,
    sigma_penalty = sigma_penalty_smooth,
    label = "Smoothed RV"
  )
  
  # ---------------------------------------------------------------
  # Combine plots
  # ---------------------------------------------------------------
  cat("Combining plots...\n")
  
  viterbi_combined <- res_trueV$viterbi_plot | res_RV$viterbi_plot | res_smooth$viterbi_plot
  viterbi_combined <- viterbi_combined +
    plot_annotation(
      title = "Topological HMM (OM Functional): Estimated Regime Sequence",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    )
  
  ci_combined <- res_trueV$ci_plot | res_RV$ci_plot | res_smooth$ci_plot
  ci_combined <- ci_combined +
    plot_annotation(
      title = "CIR Confidence Intervals: Topological HMM Estimates",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    )
  
  # ---------------------------------------------------------------
  # Save plots
  # ---------------------------------------------------------------
  if (!is.null(save_viterbi_path)) {
    ggsave(save_viterbi_path, viterbi_combined, width = width, height = height, dpi = dpi)
    cat(sprintf("Viterbi plot saved to: %s\n", save_viterbi_path))
  }
  
  if (!is.null(save_ci_path)) {
    ggsave(save_ci_path, ci_combined, width = width, height = height, dpi = dpi)
    cat(sprintf("CI plot saved to: %s\n", save_ci_path))
  }
  
  cat("\n=================================================================\n")
  cat("Topological HMM Analysis Complete!\n")
  cat("=================================================================\n")
  
  # ---------------------------------------------------------------
  # Return results
  # ---------------------------------------------------------------
  return(list(
    results = list(
      trueV = res_trueV,
      RV = res_RV,
      smooth = res_smooth
    ),
    params = list(
      trueV = res_trueV$param,
      RV = res_RV$param,
      smooth = res_smooth$param
    ),
    states_estimate = list(
      trueV = res_trueV$states_estimate,
      RV = res_RV$states_estimate,
      smooth = res_smooth$states_estimate
    ),
    viterbi_plot = viterbi_combined,
    ci_plot = ci_combined,
    ci_data = list(
      trueV = res_trueV$ci_data,
      RV = res_RV$ci_data,
      smooth = res_smooth$ci_data
    )
  ))
}

