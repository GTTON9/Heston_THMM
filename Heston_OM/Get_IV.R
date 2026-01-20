# ================================================================
# Volatility Skew via Dividend Yield Differential
# ================================================================

library(derivmkts)

#' Calculate Implied Volatility with Skew
#' 
#' @param batch_results Data frame with S, K_upper, K_lower, v, time
#' @param r Risk-free rate
#' @param H Time to maturity (years)
#' @param d_call Dividend yield for Call pricing
#' @param d_put Dividend yield for Put pricing

calculate_IV_with_dividend_skew <- function(batch_results,
                                            r = 0.02,
                                            H = 5/252,
                                            dividend_yield = 0.03,
                                            days_in_week = 5) {

  library(derivmkts)

  weekly_yield <- dividend_yield / 52
  compound_factor <- (1 + weekly_yield)^days_in_week

  # Upper (Call)
  batch_results$C_theoretical_upper <- sapply(1:nrow(batch_results), function(i) {
    derivmkts::bscall(
      s = batch_results$S[i],
      k = batch_results$K_upper[i],
      v = sqrt(batch_results$v[i]),
      tt = H, r = r, d = 0
    )
  })

  batch_results$sigma_upper_IV <- sapply(1:nrow(batch_results), function(i) {
    tryCatch({
      derivmkts::bscallimpvol(
        s = batch_results$S[i],
        k = batch_results$K_upper[i],
        price = batch_results$C_theoretical_upper[i],
        tt = H, r = r, d = 0
      )
    }, error = function(e) NA)
  })


  P_base <- sapply(1:nrow(batch_results), function(i) {
    derivmkts::bsput(  # ← 改为 bsput
      s = batch_results$S[i],
      k = batch_results$K_lower[i],
      v = sqrt(batch_results$v[i]),
      tt = H, r = r, d = 0
    )
  })

  batch_results$P_theoretical_lower <- P_base * compound_factor

  batch_results$sigma_lower_IV <- sapply(1:nrow(batch_results), function(i) {
    tryCatch({
      derivmkts::bsputimpvol(
        s = batch_results$S[i],
        k = batch_results$K_lower[i],
        price = batch_results$P_theoretical_lower[i],
        tt = H, r = r, d = 0
      )
    }, error = function(e) NA)
  })

  batch_results$current_sigma <- sqrt(batch_results$v)
  batch_results$skew <- batch_results$sigma_lower_IV - batch_results$sigma_upper_IV

  return(batch_results)
}

# 
# calculate_IV_with_dividend_skew <- function(batch_results,
#                                             r = 0.02,
#                                             H = 5/252,
#                                             d_call = 0,      
#                                             d_put = 0.03,    
#                                             d_iv = 0) {
# 
#   library(derivmkts)
#   
#   # Calculate Call prices (no dividend)
#   batch_results$C_theoretical_upper <- sapply(1:nrow(batch_results), function(i) {
#     derivmkts::bscall(
#       s = batch_results$S[i],
#       k = batch_results$K_upper[i],
#       v = sqrt(batch_results$v[i]),
#       tt = H,
#       r = r,
#       d = d_call
#     )
#   })
#   
#   # Calculate Put prices (with dividend)
#   batch_results$P_theoretical_lower <- sapply(1:nrow(batch_results), function(i) {
#     derivmkts::bsput(
#       s = batch_results$S[i],
#       k = batch_results$K_lower[i],
#       v = sqrt(batch_results$v[i]),
#       tt = H,
#       r = r,
#       d = d_put
#     )
#   })
#   
#   # Invert for Upper IV
#   batch_results$sigma_upper_IV <- sapply(1:nrow(batch_results), function(i) {
#     tryCatch({
#       derivmkts::bscallimpvol(
#         s = batch_results$S[i],
#         k = batch_results$K_upper[i],
#         price = batch_results$C_theoretical_upper[i],
#         tt = H,
#         r = r,
#         d = d_iv
#       )
#     }, error = function(e) NA)
#   })
#   
#   # Invert for Lower IV
#   batch_results$sigma_lower_IV <- sapply(1:nrow(batch_results), function(i) {
#     tryCatch({
#       derivmkts::bsputimpvol(
#         s = batch_results$S[i],
#         k = batch_results$K_lower[i],
#         price = batch_results$P_theoretical_lower[i],
#         tt = H,
#         r = r,
#         d = d_iv
#       )
#     }, error = function(e) NA)
#   })
#   
#   # Calculate statistics
#   batch_results$current_sigma <- sqrt(batch_results$v)
#   batch_results$skew <- batch_results$sigma_lower_IV - batch_results$sigma_upper_IV
#   
#   return(batch_results)
# }



#' Generate Volatility Smile using dividend skew logic
#' Based on calculate_IV_with_dividend_skew
batch_results_with_IV <- calculate_IV_with_dividend_skew(
  batch_results,
  r = 0.02,
  H = 5/252,
  d_call = 0,      
  d_put = 0.03,    
  d_iv = 0         
)


batch_results_with_IV <- calculate_IV_with_dividend_skew(
batch_results,
r = 0.02,
H = 5/252,
dividend_yield = 0.03,
days_in_week = 5
)




# ================================================================
# Plot
# ================================================================
with(batch_results_with_IV[1:50,], {
  
  n <- length(time)
  shift <- 5
  
  plot(time, current_sigma, 
       type = "l", 
       col = "blue", 
       lwd = 2,
       main = "Volatility: True vs Implied (5-Day Forward)",
       ylab = "Volatility (σ)", 
       xlab = "Time",
       xlim = range(time),
       ylim = range(c(current_sigma, sigma_upper_IV, sigma_lower_IV), na.rm = TRUE))
  
  # Shift IV by 5 days forward
  if (n > shift) {
    lines(time[(shift+1):n], sigma_upper_IV[1:(n-shift)], 
          col = "darkgreen", lwd = 2, lty = 2)
    lines(time[(shift+1):n], sigma_lower_IV[1:(n-shift)], 
          col = "darkred", lwd = 2, lty = 2)
  }
  
  # Add vertical line to show the shift
  abline(v = time[shift], lty = 3, col = "gray50")
  text(time[shift], max(current_sigma, na.rm = TRUE), 
       "5-day horizon", pos = 4, cex = 0.7, col = "gray50")
  
  legend("topright",
         legend = c("True σ (current)", 
                    "Upper IV (5d fwd)", 
                    "Lower IV (5d fwd)"),
         col = c("blue", "darkgreen", "darkred"),
         lwd = 2,
         lty = c(1, 2, 2),
         bg = "white",
         cex = 0.5,        
         y.intersp = 0.8,  
         x.intersp = 0.8,  
         seg.len = 1.5     
  )
  
  grid()
})




# View results
head(batch_results_with_IV[, c("time", "current_sigma", 
                               "sigma_upper_IV", "sigma_lower_IV", "skew")])














with(batch_results_with_IV[1:20,], {
  
  n <- length(time)
  shift <- 5
  
  plot(time, current_sigma, 
       type = "n",  
       main = "IV Forecast: Should Cover Future Volatility Path",
       ylab = "Volatility (σ)", 
       xlab = "Time",
       ylim = range(c(current_sigma, sigma_upper_IV, sigma_lower_IV), na.rm = TRUE))
  
  # 为每个时间点画出IV预测的路径
  for (t in 1:min(40, n - shift)) {
    
    # ===== 真实的未来波动率路径 =====
    future_path <- current_sigma[t:(t + shift)]
    time_path <- time[t:(t + shift)]
    lines(time_path, future_path, col = rgb(0, 0, 1, 0.3), lwd = 1)
    
    # ===== Upper IV 预测路径 =====
    # 从 current_sigma[t] 线性到 sigma_upper_IV[t]
    segments(time[t], current_sigma[t], 
             time[t + shift], sigma_upper_IV[t],
             col = "darkgreen", lwd = 1.5, lty = 2)
    
    # ===== Lower IV 预测路径 =====
    # 从 current_sigma[t] 线性到 sigma_lower_IV[t]
    segments(time[t], current_sigma[t], 
             time[t + shift], sigma_lower_IV[t],
             col = "darkred", lwd = 1.5, lty = 2)
  }
  
  # 当前波动率（粗线）
  lines(time, current_sigma, col = "blue", lwd = 2)
  
  legend("topright",
         legend = c("Current Vol", 
                    "Future Vol Paths", 
                    "Upper IV Path", 
                    "Lower IV Path"),
         col = c("blue", rgb(0,0,1,0.3), "darkgreen", "darkred"),
         lwd = c(2, 1, 1.5, 1.5),
         lty = c(1, 1, 2, 2),
         cex = 0.7)
  
  grid()
})



