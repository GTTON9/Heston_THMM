source("Model_Simulation.R")
source("Heston_controls.R")
source("Heston_data.R")
source("Heston_fit_model.R")
source("Heston_parameters.R")
source("Heston_get_init.R")
source("Heston_Conversion.R")
source("LL_HMM_R.R")
source("Heston_reorder_states.R")
source("fHMM_model.R")
source("Heston_decode.R")
source("Heston_likelihood.R")
source("Viterbi_Visual.R")


Gen_fit <- function(Gamma, mu, kappa, theta, sigma, rho, input_ = "V_", gen_length = 252, v0 = 0.03, S0 = 100,  V_ = TRUE , plot_path = TRUE, seed = 999){
  # V_: whether we know the volatility process
  
  set.seed(999)
  
  
  if(input_ == "V_"){
    N <- gen_length
    Reg_chain <- simulate_Reg(series_length = N, Reg_tran = Gamma)
    Reg_param <- cbind(mu, kappa, theta, sigma, rho)
    sim_series  <- simulate_heston(S0, v0, Reg_chain, Reg_param, T, N, M=1, method = "E", seed = seed)
    S_simulated <- sim_series$S_paths
    print(S_simulated)
    V_simulated <- sim_series$V_paths
    
    
    par(mfrow = c(1, 2))
    if(plot_path){
      plot(
        S_simulated,
        type = "l",
        main = "S path",
        col = "blue")
      
      plot(
        V_simulated,
        type = "l",
        main = "V path",
        col = "red"
      )
    }
    par(mfrow = c(1, 1))
    
    
    series_length <- length(S_simulated * 100)
    start_date <- as.Date("2024-01-01")
    date_sequence <- seq(from = start_date, by = "day", length.out = series_length)
    
    
    
    my_data_df <- data.frame(
      Date = date_sequence,
      Var = V_simulated
    )
    
    
    series_control <- Heston_set_controls( 
      states      = 2,     # 2 state
      sdds        = "Heston",         
      date_column = "Date",
      file        = my_data_df, 
      data_column = "Var",      
      logreturns  = FALSE,         
      from        = date_sequence[1],             
      to          = date_sequence[length(date_sequence)],
      runs = 10
    )
    
    
    data_hmm <- prepare_data(series_control)
    model_hmm <- Heston_fit_model(data_hmm) 
    
    final_model <- decode_states_heston(model_hmm) 
    states_estimate <- final_model$decoding
    states_estimate
    
    plot(states_estimate, col = "blue")
    lines(Reg_chain+1)
    param <- parUncon2par_heston(final_model$estimate, series_control, FALSE, numerical_safeguard = TRUE)
    
    
    p <- plot_viterbi(V_simulated, nstates, param$Gamma, 
                      param$kappa, param$theta, 
                      param$sigma, Reg_chain)
    plot(p)
    
    result <- plot_cir_confidence_simple(
      true_vol = V_simulated,
      param = param,
      states_estimate = states_estimate,
      dt = 1/252,
      interval_step = 10
    )
    
    print(result$plot)
    
  }else if(input_ == "S_"){
    
    N <- gen_length * 100
    Reg_chain_year <- simulate_Reg(series_length = N/100, Reg_tran = Gamma)
    Reg_chain <- rep(Reg_chain_year, each = 100)
    
    Reg_param <- cbind(mu, kappa, theta, sigma, rho)
    sim_series  <- simulate_heston(S0, v0, Reg_chain, Reg_param, T, N, M=1, method = "E",interp = T, seed = seed)
    S_simulated <- sim_series$S_paths

    V_simulated <- sim_series$V_paths
    
    
    par(mfrow = c(1, 2))
    if(plot_path){
      plot(
        S_simulated,
        type = "l",
        main = "S path",
        col = "blue")
      
      plot(
        V_simulated,
        type = "l",
        main = "V path",
        col = "red"
      )
    }
    par(mfrow = c(1, 1))
    
    
    
    
    S <- S_simulated
    n_days <- gen_length
    n_intraday <- 100
    
    RV_V <- numeric(n_days)
    
    for (t in 1:n_days) {
      idx <- ((t - 1) * n_intraday + 1):(t * n_intraday)
      S_day <- S[idx]
      r_day <- diff(log(S_day))      # length 3
      RV_V[t] <- sum(r_day^2) * 250 
    }
    
    
    plot(V_simulated[seq(1, length(V_simulated), by = 100)], type = 'l', col = "black")
    lines(RV_V, col = "blue")
    lines(lowess(RV_V, f = 0.1), col = "red")

    
    
    # RV_V <- lowess(RV_V, f = 0.1)$y

    
    
    
    series_length <- length(RV_V)
    start_date <- as.Date("2024-01-01")
    date_sequence <- seq(from = start_date, by = "day", length.out = series_length)
    
    
    my_data_df <- data.frame(
      Date = date_sequence,
      Var = RV_V
    )
    colnames(my_data_df) <- c("Date", "Var")


    
    # followed the example
    series_control <- Heston_set_controls( 
      states      = 2,     # 2 state
      sdds        = "Heston",         
      date_column = "Date",
      file        = my_data_df, 
      data_column = "Var",      
      logreturns  = FALSE,         
      from        = date_sequence[1],             
      to          = date_sequence[length(date_sequence)],
      runs = 10
    )
    
    
    data_hmm <- prepare_data(series_control)
    model_hmm <- Heston_fit_model(data_hmm) 
    
    final_model <- decode_states_heston(model_hmm) 
    states_estimate <- final_model$decoding
    states_estimate
    
    plot(states_estimate, col = "blue")
    lines(Reg_chain_year+1)
    param <- parUncon2par_heston(final_model$estimate, series_control, FALSE, numerical_safeguard = TRUE)
    
    
    
    
    
    p <- plot_viterbi(RV_V, nstates, param$Gamma, 
                      param$kappa, param$theta, 
                      param$sigma, Reg_chain_year)
    plot(p)
    
    
    result <- plot_cir_confidence_simple(
      true_vol = V_simulated,
      param = param,
      states_estimate = states_estimate,
      dt = 1/252,
      interval_step = 10
    )

    print(result$plot)
    
    
  }
  
  
  
  
  
  return(list(states_estimate = states_estimate, 
              param = param,
              fisher_inverse = final_model$inverse_fisher,
              S_path = S_simulated))
}








MC_Gen_fit <- function(
    n_rep = 20,
    seeds = 1:20,
    Gamma, mu, kappa, theta, sigma, rho,
    input_ = "S_",
    gen_length = 252,
    v0 = 0.03,
    S0 = 100,
    V_ = TRUE,
    plot_path = FALSE   # Monte Carlo 时一定关掉
) {
  
  stopifnot(length(seeds) == n_rep)
  
  results <- vector("list", n_rep)
  
  for (i in seq_len(n_rep)) {
    cat("Replication:", i, " Seed:", seeds[i], "\n")
    
    fit_i <- Gen_fit(
      Gamma  = Gamma,
      mu     = mu,
      kappa  = kappa,
      theta  = theta,
      sigma  = sigma,
      rho    = rho,
      input_ = input_,
      gen_length = gen_length,
      v0 = v0,
      S0 = S0,
      V_ = V_,
      plot_path = FALSE,
      seed = seeds[i]
    )
    
    results[[i]] <- list(
      seed   = seeds[i],
      Gamma  = fit_i$param$Gamma,
      kappa  = fit_i$param$kappa,
      theta  = fit_i$param$theta,
      sigma  = fit_i$param$sigma,
      fisher_inverse = fit_i$fisher_inverse,
      S_path = fit_i$S_path
    )
  }
  
  return(results)
}





MC_results <- MC_Gen_fit(
  n_rep = 20,
  seeds = 1:20,
  Gamma = Gamma,
  mu    = mu,
  kappa = kappa,
  theta = theta,
  sigma = sigma,
  rho   = rho,
  input_ = "V_"
)









# Gamma
Gamma11 <- sapply(MC_results, function(x) x$Gamma[1,1])
Gamma22 <- sapply(MC_results, function(x) x$Gamma[2,2])

# kappa, theta, sigma (each has 2 states)
kappa_mat <- do.call(rbind, lapply(MC_results, function(x) x$kappa))
theta_mat <- do.call(rbind, lapply(MC_results, function(x) x$theta))
sigma_mat <- do.call(rbind, lapply(MC_results, function(x) x$sigma))



par(mfrow = c(1,5), mar = c(5,4,3,1))  # 下边距略大

# 1. Gamma[11]
boxplot(Gamma11,
        main = expression(Gamma[11]),
        ylab = "",
        xaxt = "n")
abline(h = Gamma[1,1], col = "red", lwd = 2)
axis(1, at = 1, labels = "Reg 1", cex.axis = 0.7)
mtext(paste0("true = ", round(Gamma[1,1],3)),
      side = 1, line = 2.5, at = 1, cex = 0.7)

# 2. Gamma[22]
boxplot(Gamma22,
        main = expression(Gamma[22]),
        ylab = "",
        xaxt = "n")
abline(h = Gamma[2,2], col = "red", lwd = 2)
axis(1, at = 1, labels = "Reg 2", cex.axis = 0.7)
mtext(paste0("true = ", round(Gamma[2,2],3)),
      side = 1, line = 2.5, at = 1, cex = 0.7)

# 3. kappa (两状态)
bp <- boxplot(kappa_mat,
              main = expression(kappa),
              ylab = "",
              xaxt = "n")
segments(0.7, kappa[1], 1.3, kappa[1], col = "red", lwd = 2)
segments(1.7, kappa[2], 2.3, kappa[2], col = "red", lwd = 2)
axis(1, at = 1:2, labels = c("Reg 1","Reg 2"), cex.axis = 0.58)
mtext(paste0("true = (", round(kappa[1],3), ", ", round(kappa[2],3), ")"),
      side = 1, line = 2.5, at = 1.5, cex = 0.7)

# 4. theta (两状态)
bp <- boxplot(theta_mat,
              main = expression(theta),
              ylab = "",
              xaxt = "n")
segments(0.7, theta[1], 1.3, theta[1], col = "red", lwd = 2)
segments(1.7, theta[2], 2.3, theta[2], col = "red", lwd = 2)
axis(1, at = 1:2, labels = c("Reg 1","Reg 2"), cex.axis = 0.58)
mtext(paste0("true = (", round(theta[1],3), ", ", round(theta[2],3), ")"),
      side = 1, line = 2.5, at = 1.5, cex = 0.7)

# 5. sigma (两状态)
bp <- boxplot(sigma_mat,
              main = expression(sigma),
              ylab = "",
              xaxt = "n")
segments(0.7, sigma[1], 1.3, sigma[1], col = "red", lwd = 2)
segments(1.7, sigma[2], 2.3, sigma[2], col = "red", lwd = 2)
axis(1, at = 1:2, labels = c("Reg 1","Reg 2"), cex.axis = 0.58)
mtext(paste0("true = (", round(sigma[1],3), ", ", round(sigma[2],3), ")"),
      side = 1, line = 2.5, at = 1.5, cex = 0.7)

# 恢复默认画布
par(mfrow = c(1, 1))



# 构建矩阵，每列是一条路径
S_paths_mat <- sapply(MC_results[1:20], function(x) x$S_path)

# matplot 每列作为一条线
matplot(S_paths_mat, type="l", lty=1, col=rainbow(20),
        main="20 Simulated Paths", xlab="Time", ylab="S")
