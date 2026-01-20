

regime_params <- list(
  list(mu = 0.05, kappa = 10, theta = 0.1, sigma = 0.1, rho = -0.7),
  list(mu = 0.10, kappa = 5, theta = 0.5, sigma = 0.1, rho = -0.7)
)

Gamma <- matrix(c(0.99, 0.01, 0.01, 0.99), 2, 2)
N <- 252
v0 <- 0.1
S0 <- 100
nstates = 2
T = 1
n_days <- 252
n_intraday <- 2400




result <- simulate_heston_with_IV(
  S0 = S0,
  v0 = v0,
  regime_params = regime_params,
  Gamma = Gamma,
  N = N,
  seed = 999
)


S_daily <- result$S
V_daily <- result$v
colnames(result)
batch_results <- calculate_expected_moves_batch(S_daily, V_daily)

plot_expected_moves_with_shift(batch_results, H_days = 5)


result_with_paths <- construct_linear_paths(batch_results, H = 5/252, result)
result_complete <- calculate_OM_batch(result_with_paths, regime_params_for_OM)












OM_viterbi_result <- viterbi_om_pure(OM_upper= result_complete[, c("OM_upper_R1", "OM_upper_R2")],
                            OM_lower= result_complete[, c("OM_lower_R1", "OM_lower_R2")],
                            nstates,
                            Gamma,
                            lambda_om = 1.0,
                            use_upper = TRUE,
                            use_lower = TRUE)

plot(OM_viterbi_result$states_estimate, col="blue")
lines(Reg_chain_year+1)
#' Viterbi Algorithm for Pure O-M Functional (THMM)
#' 
#' Finds the single most likely sequence of hidden states by maximizing 
#' the joint probability P(S_1:T, φ_1:T).
#'
#' @param OM_upper Matrix [T x nstates] of O-M functionals
#' @param OM_lower Matrix [T x nstates] of O-M functionals
#' @param nstates Number of hidden states
#' @param Gamma Transition probability matrix
#' @param lambda_om Weight for O-M contribution
#' 
#' @return List with states_viterbi, log_prob, and path_matrix



plot2 <- plot_viterbi_om(
  batch_results_complete = batch_results_complete,
  nstates = 2,
  Gamma = Gamma,
  kappa = c(10,5),
  theta = c(0.1, 0.5),
  sigma = c(0.1, 0.1),
  Reg_chain = Reg_chain_year ,
  lambda_om = 1,
  normalize_method = "log",
  show_om_contribution = TRUE
)






