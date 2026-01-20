
calculate_parkinson_from_series <- function(S_simulated, obs_per_day) {
  
  total_obs <- length(S_simulated)
  rows_per_day <- obs_per_day + 1 
  
  # 1. 确定可以完整覆盖的天数
  # 总价格观测点数量 - 1 = 总收益率数量
  # 总收益率数量 / 每日收益率数量 = 完整天数
  # S_simulated[1] 是第一天的开始价格，S_simulated[obs_per_day + 1] 是第一天的结束价格。
  
  # 完整天数的计算：
  num_full_days <- floor((total_obs - 1) / obs_per_day)
  
  if (num_full_days == 0) {
    stop("观测数据太少，不足以形成一整天的完整观测。")
  }
  
  # 2. 裁剪价格序列，使其长度恰好匹配矩阵维度
  # 我们需要 (num_full_days * obs_per_day) 个收益率，即 (num_full_days * obs_per_day + 1) 个价格点。
  
  required_length <- num_full_days * obs_per_day + 1 
  S_trimmed <- S_simulated[1:required_length]
  print(length(S_trimmed))
  # 3. 将连续价格转换为每日观测矩阵
  # 矩阵的行数是 obs_per_day + 1，列数是 num_full_days
  # 每次分块包含 (obs_per_day + 1) 个价格点，代表一个交易日的价格路径。
  daily_price_segments <- matrix(
    S_trimmed, 
    nrow = rows_per_day,        # 11 行
    ncol = num_full_days,       # 227 列
    byrow = FALSE
  )
  # 4. 提取每日极值并计算 Parkinson 估计值 (略...)
  daily_high <- apply(daily_price_segments, 2, max)
  daily_low <- apply(daily_price_segments, 2, min)
  
  # [Parkinson 估计值的计算公式]
  scaling_factor <- 1 / (4 * log(2))
  ln_ratio_squared <- log(daily_high / daily_low)^2
  parkinson_volatility_daily <- sqrt(scaling_factor * ln_ratio_squared)
  
  return(parkinson_volatility_daily)
}



sim_series  <- simulate_heston(S0, v0, Reg_chain, Reg_param, T = 0.5, N, M=1, method = "E")
S_simulated <- sim_series$S_paths
plot(S_simulated,  type = "l")

N = 250
Parkinson_estimates <- calculate_parkinson_from_series(
  S_simulated = S_simulated, 
  obs_per_day = 100
)


# 
# 
# plot(V_simulated, type = 'l', ylim = c(0,1))
# lines(volatility_proxy, col = "blue")
# 
# 
# print(paste("总价格观测点数:", length(S_simulated_example)))
# print(paste("每日观测次数:", obs_per_day_input))
# print(paste("计算出的天数:", length(Parkinson_estimates)))
# cat("\n每日 Parkinson 波动率估计值:\n")
# print(Parkinson_estimates)




