# ================================================================
# calculate_IV_with_dividend_skew
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


