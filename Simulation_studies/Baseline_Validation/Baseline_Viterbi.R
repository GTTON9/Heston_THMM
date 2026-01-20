# source("Model_Simulation.R")
# source("Heston_likelihood.R")
# source("Heston_decode.R")
# 

# easy case

Gamma <- matrix(c(0.99, 0.01, 0.01, 0.99), 2, 2)
N <- 252

v0 <- 0.1
S0 <- 100
nstates = 2
T = 1
n_days <- 252
n_intraday <- 2400
mu <- c(0.3, 0.3)
kappa <- c(10, 5)
theta <- c(0.1, 0.5)
sigma <- c(0.1, 0.1)
rho <- c(-0.1, -0.1)





set.seed(999)
Reg_chain <- simulate_Reg(series_length = N)
Reg_param <- cbind(mu, kappa, theta, sigma, rho)
sim_series  <- simulate_heston(S0, v0, Reg_chain, Reg_param, T, N, M=1, method = "E_C")
S_simulated <- sim_series$S_paths
plot(S_simulated, type = "l")
V_simulated <- sim_series$V_paths
plot(V_simulated, type = "l")


states_estimate <- viterbi(V_simulated, nstates, Gamma, kappa, theta, sigma)
states_estimate

plot(states_estimate, col = "blue")
lines(Reg_chain+1)

p <- plot_viterbi(RV_V, nstates, Gamma, 
                  kappa, theta, 
                  sigma, Reg_chain_year)
plot(p)





# hard case

Gamma <- matrix(c(0.99, 0.01, 0.01, 0.99), 2, 2)
N <- 250
v0 <- 1
S0 <- 100
nstates = 2
mu <- c(0.5, 0.5)
kappa <- c(2, 1)
theta <- c(0.1, 0.2)
sigma <- c(0.1, 0.1)
rho <- c(-0.1, -0.1)




set.seed(999)
Reg_chain <- simulate_Reg(series_length = N)
Reg_param <- cbind(mu, kappa, theta, sigma, rho)
sim_series  <- simulate_heston(S0, v0, Reg_chain, Reg_param, T, N, M=1, method = "E_C")
S_simulated <- sim_series$S_paths
plot(S_simulated, type = "l")
V_simulated <- sim_series$V_paths
plot(V_simulated, type = "l")




states_estimate <- viterbi(V_simulated, nstates, Gamma, kappa, theta, sigma)
states_estimate

plot(states_estimate, col = "blue")
lines(Reg_chain+1)






