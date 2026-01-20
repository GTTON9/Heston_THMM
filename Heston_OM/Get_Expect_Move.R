# ================================================================
# calculate_expected_move
# ================================================================
#' Calculate Expected Move Strikes using Black-Scholes Delta
#' This function finds the strike prices K_upper and K_lower such that
#' the Black-Scholes delta equals 0.1 (for calls) and -0.1 (for puts)

#' @param S Current asset price
#' @param v Current variance (volatility^2)
#' @param r Risk-free rate (annualized)
#' @param H Time to maturity in years (default = 5/252 for 5 trading days)
#' @param strike_step Grid search step size (default = 1.0)
#' @param search_range_pct Percentage range to search around S (default = 0.3 = 30%)

#' @return A list containing:
#'   - K_upper: Strike price with delta = 0.1 (call)
#'   - K_lower: Strike price with delta = -0.1 (put)
#'   - delta_upper: Actual delta at K_upper
#'   - delta_lower: Actual delta at K_lower
calculate_expected_move <- function(S, v, r = 0.02, H = 5/252, 
                                    strike_step = 1.0, 
                                    search_range_pct = 0.3) {
  
  # Convert variance to volatility
  sigma <- sqrt(v)
  
  # Helper function: Black-Scholes d1
  calculate_d1 <- function(S, K, r, sigma, T) {
    (log(S/K) + (r + 0.5 * sigma^2) * T) / (sigma * sqrt(T))
  }
  
  # Helper function: Black-Scholes Delta
  bs_delta <- function(S, K, r, sigma, T, option_type = "call") {
    d1 <- calculate_d1(S, K, r, sigma, T)
    
    if (option_type == "call") {
      delta <- pnorm(d1)
    } else if (option_type == "put") {
      delta <- pnorm(d1) - 1
    } else {
      stop("option_type must be 'call' or 'put'")
    }
    
    return(delta)
  }
  
  # ============================================================
  # STEP 1: Find K_upper (Call option with delta = 0.1)
  # ============================================================
  
  # Search range: from current price to 30% above
  K_upper_range <- seq(from = S, 
                       to = S * (1 + search_range_pct), 
                       by = strike_step)
  
  # Calculate delta for each strike
  deltas_upper <- sapply(K_upper_range, function(K) {
    bs_delta(S, K, r, sigma, H, "call")
  })
  
  # Find the strike closest to delta = 0.1
  target_delta_upper <- 0.1
  errors_upper <- abs(deltas_upper - target_delta_upper)
  best_idx_upper <- which.min(errors_upper)
  
  K_upper <- K_upper_range[best_idx_upper]
  delta_upper <- deltas_upper[best_idx_upper]
  
  # ============================================================
  # STEP 2: Find K_lower (Put option with delta = -0.1)
  # ============================================================
  
  # Search range: from 30% below to current price
  K_lower_range <- seq(from = S * (1 - search_range_pct), 
                       to = S, 
                       by = strike_step)
  
  # Calculate delta for each strike
  deltas_lower <- sapply(K_lower_range, function(K) {
    bs_delta(S, K, r, sigma, H, "put")
  })
  
  # Find the strike closest to delta = -0.1
  target_delta_lower <- -0.1
  errors_lower <- abs(deltas_lower - target_delta_lower)
  best_idx_lower <- which.min(errors_lower)
  
  K_lower <- K_lower_range[best_idx_lower]
  delta_lower <- deltas_lower[best_idx_lower]
  
  # ============================================================
  # RETURN RESULTS
  # ============================================================
  
  return(list(
    K_upper = K_upper,
    K_lower = K_lower,
    delta_upper = delta_upper,
    delta_lower = delta_lower,
    error_upper = errors_upper[best_idx_upper],
    error_lower = errors_lower[best_idx_lower]
  ))
}


# ================================================================
# calculate_expected_moves_batch
# ================================================================
#' Calculate Expected Moves for Multiple Time Points

#' @param S_vec Vector of asset prices at each time point
#' @param v_vec Vector of variances at each time point
#' @param r Risk-free rate
#' @param H Time to maturity (default = 5/252)

#' @return A data frame with columns:
#'   - time: Time index
#'   - S: Asset price
#'   - v: Variance
#'   - K_upper: Upper strike
#'   - K_lower: Lower strike
#'   - delta_upper: Actual delta at K_upper
#'   - delta_lower: Actual delta at K_lower

calculate_expected_moves_batch <- function(S_vec, v_vec, r = 0.02, H = 5/252) {
  
  T <- length(S_vec)
  
  # Pre-allocate results
  results <- data.frame(
    time = 1:T,
    S = S_vec,
    v = v_vec,
    K_upper = numeric(T),
    K_lower = numeric(T),
    delta_upper = numeric(T),
    delta_lower = numeric(T),
    error_upper = numeric(T),
    error_lower = numeric(T)
  )
  
  # Calculate for each time point
  for (t in 1:T) {
    em <- calculate_expected_move(S_vec[t], v_vec[t], r, H)
    
    results$K_upper[t] <- em$K_upper
    results$K_lower[t] <- em$K_lower
    results$delta_upper[t] <- em$delta_upper
    results$delta_lower[t] <- em$delta_lower
    results$error_upper[t] <- em$error_upper
    results$error_lower[t] <- em$error_lower
  }
  
  return(results)
}



plot_expected_moves_with_shift <- function(batch_results,
                                           H_days = 5,
                                           title = "Expected Moves: Current vs 1-Week Ahead") {
  
  T <- nrow(batch_results)
  
  batch_results_shifted <- batch_results %>%
    mutate(
      time_future = time + H_days,
      K_upper_shifted = K_upper,
      K_lower_shifted = K_lower
    )
  
  # Create plot
  p <- ggplot() +
    # Current price
    geom_line(
      data = batch_results,
      aes(x = time, y = S, color = "Current Price"),
      linewidth = 1.2
    ) +
    # Expected upper move (shifted to future)
    geom_line(
      data = batch_results_shifted,
      aes(x = time_future, y = K_upper_shifted, color = "Expected Upper (in 1 week)"),
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    # Expected lower move (shifted to future)
    geom_line(
      data = batch_results_shifted,
      aes(x = time_future, y = K_lower_shifted, color = "Expected Lower (in 1 week)"),
      linewidth = 0.8,
      linetype = "dashed"
    ) +
    # Shaded confidence region (shifted)
    geom_ribbon(
      data = batch_results_shifted,
      aes(x = time_future, ymin = K_lower_shifted, ymax = K_upper_shifted),
      fill = "lightblue", alpha = 0.2
    ) +
    # Add arrows to show "1 week ahead"
    geom_segment(
      data = batch_results[1:3, ],  # Just show a few examples
      aes(x = time, xend = time + H_days, y = S, yend = S),
      arrow = arrow(length = unit(0.2, "cm")),
      color = "gray50", linetype = "dotted", alpha = 0.5
    ) +
    scale_color_manual(
      name = "Series",
      values = c(
        "Current Price" = "darkblue",
        "Expected Upper (in 1 week)" = "darkgreen",
        "Expected Lower (in 1 week)" = "darkred"
      )
    ) +
    labs(
      title = title,
      subtitle = sprintf("Expected moves are shifted %d days forward to show future alignment", H_days),
      x = "Time Index",
      y = "Price ($)"
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

