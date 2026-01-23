library(ggplot2)
library(tidyr)
library(dplyr)


forward_filter <- function(observations, nstates, Gamma, kappa, theta, sigma) {
  
  T <- length(observations) - 1
  logGamma <- log(Gamma)
  delta <- log(oeli::stationary_distribution(Gamma))
  
  # emission log-prob
  allprobs <- matrix(0, nstates, T)
  for (i in 1:nstates) {
    allprobs[i,] <- get_transition_density_heston_ln(observations, kappa[i], theta[i], sigma[i])
  }
  
  logalpha <- matrix(-Inf, nstates, T)
  
  # initialization
  logalpha[,1] <- delta + allprobs[,1]
  
  # recursion
  for (t in 2:T) {
    for (j in 1:nstates) {
      logalpha[j,t] <- logsumexp(logalpha[,t-1] + logGamma[,j]) + allprobs[j,t]
    }
  }
  
  # normalized posterior: p(s_t=j | y_1:t)
  posterior <- matrix(0, nstates, T)
  reg <- matrix(0, 1, T)
  for (t in 1:T) {
    c <- logsumexp(logalpha[,t])
    posterior[,t] <- exp(logalpha[,t] - c)
    reg[,t] <- which.max(posterior[,t])
  }
  
  return(list(posterior = posterior,
              states_estimate = as.vector(reg)))
}

logsumexp <- function(x) {
  
  m <- max(x)
  
  m + log(sum(exp(x - m)))
  
}

plot_viterbi <- function(V_simulated, nstates, Gamma, kappa, theta, sigma, Reg_chain){
  
  
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
    Estimated_State = states_estimate -1, 
    True_State = as.factor(Reg_chain[1:ncol(a)])        
  )
  
  
print("np")
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

  
plot <- ggplot() +
    
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
    
  
    geom_line(
      data = prob_data_long, 
      aes(x = Time, y = Probability, color = Probability_Type), # Probability_Type 映射到颜色图例
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
      labels = c("0" = "Regime 1 (Calm)", "1" = "Regime 2 (Turbulent)"),
      name = "True Regime Background"
    ) +
  scale_color_manual(
    values = c("Prob_State1" = "blue", "Prob_State2" = "red"),
    labels = c("Prob_State1" = expression(P(X[1:t], S[t] == 1 * " | " * bold(theta))),
               "Prob_State2" = expression(P(X[1:t], S[t] == 2 * " | " * bold(theta)))),
               name = "Forward Probability" 
    )+
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
      title = "HMM Regime Switching Heston: Estimation",
      subtitle = sprintf("Global Optimal Path | Accuracy: %.2f%%", 
                         mean(state_data$Estimated_State == as.numeric(as.character(state_data$True_State))) * 100),
      
      x = "Time Step",
    ) +
    theme_minimal() +

    guides(
      fill = guide_legend(order = 1),       
      color = guide_legend(order = 2),
      shape = guide_legend(order = 3) 
    )
  return(plot)
}










plot_cir_confidence_simple <- function(true_vol, 
                                       param, 
                                       states_estimate, 
                                       dt = 1/252, 
                                       interval_step = 10,
                                       V_grid = NULL) {
  
  min_length <- min(length(true_vol), length(states_estimate))
  true_vol <- true_vol[1:min_length]
  states_estimate <- states_estimate[1:min_length]
  time_points <- 1:min_length
  
  
  if (is.null(V_grid)) {
    V_min <- max(0.0001, min(true_vol) * 0.1)
    V_max <- max(true_vol) * 3
    V_grid <- seq(V_min, V_max, length.out = 500)
  }
  
  
  changepoints <- c(1, which(diff(states_estimate) != 0) + 1, length(states_estimate) + 1)
  
  
  all_ci_data <- data.frame()
  
  
  for (i in 1:(length(changepoints) - 1)) {
    start_idx <- changepoints[i]
    end_idx <- changepoints[i + 1] - 1
    
    if (end_idx - start_idx < interval_step) next
    
    current_regime <- states_estimate[start_idx]
    kappa_reg <- param$kappa[current_regime]
    theta_reg <- param$theta[current_regime]
    sigma_reg <- param$sigma[current_regime]
    
    v_ref <- true_vol[start_idx]
    
    
    for (j in seq(interval_step, (end_idx - start_idx + 1), by = interval_step)) {
      t_idx <- start_idx + j - 1
      if (t_idx > end_idx) break
      
      k <- j * dt 
      
      log_densities <- sapply(V_grid, function(v) {
        ln_d_Heston(V_t = v_ref, 
                    V_t_plus_k = v, 
                    k = k, 
                    kappa = kappa_reg, 
                    theta = theta_reg, 
                    sigma = sigma_reg)
      })
      
      densities <- exp(log_densities)
      
      if (any(is.na(densities)) || all(densities == 0)) {
        warning(sprintf("Invalid densities at time %d, regime %d", t_idx, current_regime))
        next
      }
      
      dV <- diff(V_grid)[1]
      cdf <- cumsum(densities * dV)
      cdf <- cdf / max(cdf) 
      
      
      idx_lower <- which.min(abs(cdf - 0.025))
      idx_mean <- which.min(abs(cdf - 0.5))
      idx_upper <- which.min(abs(cdf - 0.975))
      
      lower_val <- V_grid[idx_lower]
      mean_val <- V_grid[idx_mean]
      upper_val <- V_grid[idx_upper]
      
      
      if (!is.na(lower_val) && !is.na(upper_val) && 
          lower_val > 0 && upper_val > lower_val) {
        
        new_row <- data.frame(
          time = t_idx,
          lower = lower_val,
          mean = mean_val,
          upper = upper_val,
          regime = current_regime
        )
        all_ci_data <- rbind(all_ci_data, new_row)
      }
    }
  }
  
  if (nrow(all_ci_data) == 0) {
    warning("No confidence intervals computed")
    return(NULL)
  }
  

  vol_data <- data.frame(
    time = time_points, 
    true = true_vol
  )
  
  print(head(all_ci_data, 10))
  

  p <- ggplot() +

    geom_line(data = vol_data, 
              aes(x = time, y = true), 
              color = "black", size = 0.5, alpha = 0.8) +
    

    geom_point(data = all_ci_data, 
               aes(x = time, y = lower, fill = factor(regime)), 
               shape = 21, size = 0.7, alpha = 0.7, stroke = 1) +
    geom_point(data = all_ci_data, 
               aes(x = time, y = mean, fill = factor(regime)), 
               shape = 21, size = 0.7, alpha = 0.7, stroke = 1)+
    geom_point(data = all_ci_data, 
               aes(x = time, y = upper, fill = factor(regime)), 
               shape = 21, size = 0.7, alpha = 0.7, stroke = 1) +
    

    geom_segment(data = all_ci_data,
                 aes(x = time, xend = time, y = lower, yend = upper, 
                     color = factor(regime)),
                 alpha = 0.3, size = 0.8) +
    
    scale_color_manual(
      values = c("1" = "blue", "2" = "red"),
      name = "Regime"
    ) +
    scale_fill_manual(
      values = c("1" = "blue", "2" = "red"),
      name = "Regime (95% CI)"
    ) +
    labs(
      title = "CIR Model Confidence Intervals",
      subtitle = sprintf("95%% CI | Every %d days | black line: True volatility", interval_step),
      x = "Time (Days)", 
      y = "Volatility"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )
  
  
  return(list(plot = p, ci_data = all_ci_data))
}







