# 加载必要的包
library(NMOF)
library(derivmkts)
library(ggplot2)

# 原始模拟函数
simulate_heston_with_IV <- function(S0 = 100,
                                    v0 = 0.04,
                                    regime_params,
                                    Gamma,
                                    T = 1,
                                    N = 252,
                                    H = 5/252,
                                    r = 0.02,
                                    strikes_relative = seq(0.85, 1.15, by = 0.001),
                                    method = "E",
                                    seed = 999) {
  
  set.seed(seed)
  Reg_chain <- numeric(N)
  Reg_chain[1] <- 0
  for (t in 2:N) {
    Reg_chain[t] <- sample(c(0, 1), 1, prob = Gamma[Reg_chain[t-1] + 1, ])
  }
  
  Reg_param <- rbind(
    c(regime_params[[1]]$mu, regime_params[[1]]$kappa, regime_params[[1]]$theta,
      regime_params[[1]]$sigma, regime_params[[1]]$rho),
    c(regime_params[[2]]$mu, regime_params[[2]]$kappa, regime_params[[2]]$theta,
      regime_params[[2]]$sigma, regime_params[[2]]$rho)
  )
  
  N_H <- N*1000
  Reg_chain_H <- rep(Reg_chain, each = 1000)
  dt <- T / N_H
  sqrt_dt <- sqrt(dt)
  S <- numeric(N_H + 1)
  V <- numeric(N_H + 1)
  S[1] <- S0
  V[1] <- v0
  
  for (i in 1:N_H) {
    reg <- Reg_chain_H[i] + 1
    mu <- Reg_param[reg, 1]
    kappa <- Reg_param[reg, 2]
    theta <- Reg_param[reg, 3]
    sigma <- Reg_param[reg, 4]
    rho <- Reg_param[reg, 5]
    
    Z1 <- rnorm(1)
    Z2 <- rnorm(1)
    dW1 <- Z1 * sqrt_dt
    dW2 <- (rho * Z1 + sqrt(1 - rho^2) * Z2) * sqrt_dt
    
    V_pos <- max(V[i], 0)
    V[i + 1] <- max(V[i] + kappa * (theta - V[i]) * dt + sigma * sqrt(V_pos) * dW2, 1e-6)
    S[i + 1] <- S[i] * exp((mu - 0.5 * V_pos) * dt + sqrt(V_pos) * dW1)
  }
  
  sample_indices <- seq(1, N_H + 1, by = 1000)
  S_daily <- S[sample_indices]
  V_daily <- V[sample_indices]
  S_daily <- S_daily[-1]
  V_daily <- V_daily[-1]
  
  n_strikes <- length(strikes_relative)
  call_IV <- matrix(NA, N, n_strikes)
  put_IV <- matrix(NA, N, n_strikes)
  
  for (t in 1:N) {
    if (t %% 50 == 0) print(paste("Processing time:", t))
    params <- regime_params[[Reg_chain[t] + 1]]
    K_grid <- S_daily[t] * strikes_relative
    
    for (k in 1:n_strikes) {
      call_price <- tryCatch(
        callHestoncf(S = S_daily[t], X = K_grid[k], tau = H, r = r, q = 0,
                     v0 = V_daily[t], vT = params$theta, rho = params$rho,
                     k = params$kappa, sigma = params$sigma),
        error = function(e) NA
      )
      
      if (!is.na(call_price) && call_price > 1e-6) {
        call_IV[t, k] <- tryCatch(
          bscallimpvol(s = S_daily[t], k = K_grid[k], price = call_price, tt = H, r = r, d = 0),
          error = function(e) NA
        )
        
        put_price <- call_price - S_daily[t] + K_grid[k] * exp(-r * H)
        if (put_price > 1e-6) {
          put_IV[t, k] <- tryCatch(
            bsputimpvol(s = S_daily[t], k = K_grid[k], price = put_price, tt = H, r = r, d = 0),
            error = function(e) NA
          )
        }
      }
    }
  }
  
  colnames(call_IV) <- paste0("CallIV_m", strikes_relative)
  colnames(put_IV) <- paste0("PutIV_m", strikes_relative)
  
  RV_V <- numeric(N)
  for (t in 1:N) {
    idx <- ((t - 1) * 1000 + 1):(t * 1000)
    S_day <- S[idx]
    r_day <- diff(log(S_day))      
    RV_V[t] <- sum(r_day^2) * 252 
  }
  
  result <- cbind(
    data.frame(time = 1:N, S = S_daily, v = V_daily, regime = Reg_chain),
    call_IV,
    put_IV,
    RV_V
  )
  
  return(result)
}

# 设置参数
regime_params <- list(
  list(mu = 0.05, kappa = 3, theta = 0.04, sigma = 0.3, rho = -0.7),
  list(mu = -0.02, kappa = 5, theta = 0.09, sigma = 0.5, rho = -0.8)
)

Gamma <- matrix(c(0.95, 0.05, 0.10, 0.90), nrow = 2, byrow = TRUE)

# 运行模拟（使用较少的strikes加快速度）
strikes_relative <- seq(0.90, 1.10, by = 0.01)
result <- simulate_heston_with_IV(
  S0 = 100,
  v0 = 0.04,
  regime_params = regime_params,
  Gamma = Gamma,
  strikes_relative = strikes_relative,
  seed = 999
)

# 提取time=1的数据
t1_data <- result[1, ]
call_cols <- grep("^CallIV_m", names(result))
moneyness <- strikes_relative
call_iv_t1 <- as.numeric(t1_data[call_cols])

# 创建数据框
plot_data <- data.frame(
  moneyness = moneyness,
  implied_vol = call_iv_t1
)
plot_data <- plot_data[!is.na(plot_data$implied_vol), ]

# 绘制volatility smile
ggplot(plot_data, aes(x = moneyness, y = implied_vol)) +
  geom_line(color = "blue", size = 1.2) +
  geom_point(color = "red", size = 2) +
  labs(
    title = "Volatility Smile at Time = 1",
    subtitle = paste("Spot Price:", round(t1_data$S, 2), 
                     "| Regime:", t1_data$regime,
                     "| Variance:", round(t1_data$v, 4)),
    x = "Moneyness (K/S)",
    y = "Implied Volatility"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.title = element_text(size = 12)
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
  annotate("text", x = 1, y = max(plot_data$implied_vol, na.rm = TRUE) * 0.95, 
           label = "ATM", size = 3, color = "gray30")

print(paste("Time 1 - S:", round(t1_data$S, 2), 
            "| v:", round(t1_data$v, 4), 
            "| Regime:", t1_data$regime))

