source("Gen_fit.R")

# easy case
# jitters 10

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


BW_result <- Gen_fit(Gamma, mu, kappa, theta, sigma, rho, 
                     gen_length = N, v0 = v0, S0 = S0)
BW_result$param
BW_result$states_estimate





# hard case
# jitters 5
source("Gen_fit.R")
N <- 252
v0 <- 1
S0 <- 100
nstates = 2
Gamma <- matrix(c(0.99, 0.01, 0.01, 0.99), 2, 2)
mu <- c(0.5, 0.5)
kappa <- c(2, 1)
theta <- c(0.1, 0.2)
sigma <- c(0.05, 0.05)
rho <- c(-0.1, -0.1)





BW_result <- Gen_fit(Gamma, mu, kappa, theta, sigma, rho, 
                     gen_length = N, v0 = v0, S0 = S0)
BW_result$param
BW_result$states_estimate





