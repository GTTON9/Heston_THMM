plot_price_and_Variance <- function(S_daily, V_daily, True_state) {
  
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  
  N <- length(S_daily)
  time_index <- 1:N
  
  plot_data <- data.frame(
    Time = time_index,
    Price = S_daily,
    Variance = V_daily,
    True_State = as.factor(True_state)
  )
  

  max_time <- max(plot_data$Time)
  regime_rects <- plot_data %>%
    mutate(change = True_State != lag(True_State, default = first(True_State))) %>%
    filter(change | row_number() == 1) %>%
    dplyr::select(Time, True_State) %>%
    mutate(Time_End = lead(Time, default = max_time + 1)) %>%
    rename(Time_Start = Time)
  

  y_min_price <- min(S_daily) * 0.95
  y_max_price <- max(S_daily) * 1.05
  
  p1 <- ggplot() +
    geom_rect(
      data = regime_rects,
      aes(xmin = Time_Start, xmax = Time_End, 
          ymin = y_min_price, ymax = y_max_price, fill = True_State),
      alpha = 0.2
    ) +
    geom_line(
      data = plot_data,
      aes(x = Time, y = Price),
      color = "darkblue",
      linewidth = 0.8
    ) +
    scale_fill_manual(
      values = c("0" = "skyblue", "1" = "salmon"),
      labels = c("0" = "Regime 1 (Calm)", "1" = "Regime 2 (Turbulent)")
    ) +
    labs(title = "Asset Price Path", x = "Time Step", y = "Price (S)") +
    theme_minimal() +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))
  

  y_min_vol <- min(V_daily) * 0.95
  y_max_vol <- max(V_daily) * 1.05
  
  p2 <- ggplot() +
    geom_rect(
      data = regime_rects,
      aes(xmin = Time_Start, xmax = Time_End, 
          ymin = y_min_vol, ymax = y_max_vol, fill = True_State),
      alpha = 0.2
    ) +
    geom_line(
      data = plot_data,
      aes(x = Time, y = Variance),
      color = "darkred",
      linewidth = 0.8
    ) +
    scale_fill_manual(
      values = c("0" = "skyblue", "1" = "salmon"),
      labels = c("0" = "Regime 1 (Calm)", "1" = "Regime 2 (Turbulent)")
    ) +
    labs(title = "Variance Path", x = "Time Step", y = "Variance (v)") +
    theme_minimal() +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))
  
  grid.arrange(p1, p2, ncol = 2)
}




















plot_price_and_variance_with_RV <- function(S_daily, V_daily, RV_V, RV_smooth, True_state) {
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(tidyr)
  
  # ---- Align lengths safely ----
  N <- min(length(S_daily), length(V_daily), length(True_state), length(RV_V), length(RV_smooth))
  S_daily   <- S_daily[1:N]
  V_daily   <- V_daily[1:N]
  True_state <- True_state[1:N]
  RV_V      <- RV_V[1:N]
  RV_smooth <- RV_smooth[1:N]
  
  time_index <- 1:N
  
  plot_data <- data.frame(
    Time = time_index,
    Price = S_daily,
    Variance = V_daily,
    RV_V = RV_V,
    RV_smooth = RV_smooth,
    True_State = as.factor(True_state)
  )
  
  # ---- Regime rectangles (background shading) ----
  max_time <- max(plot_data$Time)
  regime_rects <- plot_data %>%
    mutate(change = True_State != lag(True_State, default = first(True_State))) %>%
    filter(change | row_number() == 1) %>%
    dplyr::select(Time, True_State) %>%
    mutate(Time_End = lead(Time, default = max_time + 1)) %>%
    rename(Time_Start = Time)
  
  # ---- Determine y ranges ----
  y_min_price <- min(plot_data$Price, na.rm = TRUE) * 0.95
  y_max_price <- max(plot_data$Price, na.rm = TRUE) * 1.05
  
  y_min_vol <- min(c(plot_data$Variance, plot_data$RV_V, plot_data$RV_smooth), na.rm = TRUE) * 0.95
  y_max_vol <- max(c(plot_data$Variance, plot_data$RV_V, plot_data$RV_smooth), na.rm = TRUE) * 1.2
  
  # ---- Fill legend mapping (assume states are 0/1; still works if 1/2 but labels may need tweak) ----
  fill_vals <- c("0" = "skyblue", "1" = "salmon")
  fill_labs <- c("0" = "Regime 1 (Calm)", "1" = "Regime 2 (Turbulent)")
  
  # If your states are 1/2 instead of 0/1, uncomment:
  # fill_vals <- c("1"="skyblue","2"="salmon")
  # fill_labs <- c("1"="Regime 1 (Calm)","2"="Regime 2 (Turbulent)")
  
  # ---- Left: Price path ----
  p1 <- ggplot() +
    geom_rect(
      data = regime_rects,
      aes(xmin = Time_Start, xmax = Time_End, ymin = y_min_price, ymax = y_max_price, fill = True_State),
      alpha = 0.2
    ) +
    geom_line(
      data = plot_data,
      aes(x = Time, y = Price),
      color = "black",
      linewidth = 0.8
    ) +
    scale_fill_manual(values = fill_vals, labels = fill_labs) +
    labs(title = "Asset Price Path", x = "Time Step", y = "Price (S)") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # ---- Right: Variance + RV + Smoothed RV (same panel) ----
  vol_long <- plot_data %>%
    dplyr::select(Time, Variance, RV_V, RV_smooth) %>%
    pivot_longer(cols = c(Variance, RV_V, RV_smooth),
                 names_to = "Series", values_to = "Value")
  
  series_labs <- c(
    "Variance"  = "True variance v(t)",
    "RV_V"      = "Realised variance (from S)",
    "RV_smooth" = "Smoothed realised variance"
  )
  
  p2 <- ggplot() +
    geom_rect(
      data = regime_rects,
      aes(xmin = Time_Start, xmax = Time_End, ymin = y_min_vol, ymax = y_max_vol, fill = True_State),
      alpha = 0.2
    ) +
    # 1) Realised variance: noisy, put at bottom, semi-transparent
    geom_line(
      data = plot_data,
      aes(x = Time, y = RV_V, color = "Realised variance (from S)"),
      linewidth = 0.9, alpha = 0.55
    ) +
    # 2) Smoothed realised variance: middle
    geom_line(
      data = plot_data,
      aes(x = Time, y = RV_smooth, color = "Smoothed realised variance"),
      linewidth = 0.9, alpha = 0.75
    ) +
    # 3) True variance: draw last so it is always visible
    geom_line(
      data = plot_data,
      aes(x = Time, y = Variance, color = "True variance v(t)"),
      linewidth = 0.5, alpha = 1
    ) +
    scale_fill_manual(values = fill_vals) +
    scale_color_manual(
      values = c("True variance v(t)" = "black",
        "Realised variance (from S)" = "red",
        "Smoothed realised variance" = "green"
      ),
      name = NULL
    ) +
    labs(title = "Variance Path", x = "Time Step", y = "Variance / RV") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = c(0.6, 0.9),
      legend.justification = c(1, 1),
      legend.text = element_text(size = 9),
      legend.key.size = grid::unit(0.22, "cm")
    ) +
    guides(fill = "none")
  
  return(grid.arrange(p1, p2, ncol = 2))
}

