
# Viterbi Gaussian Approximation
get_gaussian_approx_density <- function(x_next, x_prev, kappa, theta, sigma, dt) {
  # Euler-Maruyama Expected Mean
  mu <- x_prev + kappa * (theta - x_prev) * dt
  
  # Euler-Maruyama Variance (Sigma^2 * Xt-1 * dt)
  # We use max(..., epsilon) to prevent division by zero or negative variance
  variance <- (sigma^2) * x_prev * dt
  sd <- sqrt(max(variance, 1e-8)) 
  
  # Calculate Log Probability Density
  return(dnorm(x_next, mean = mu, sd = sd, log = TRUE))
}


viterbi_mdp_heston <- function(obs, nstates, trans_matrix, params, dt, G_approx = F) {
  
  T_len <- length(obs)
  
  # --- 1. Initialization ---
  
  # Value Function Matrix (xi / delta): Stores the max accumulated log-prob up to time t
  # Rows = States, Cols = Time
  value_func <- matrix(-Inf, nrow = nstates, ncol = T_len)
  
  # Backpointer Matrix (psi): Stores the optimal previous state index
  backpointer <- matrix(0, nrow = nstates, ncol = T_len)
  
  # Initial State Probabilities (Assuming stationary or uniform start)
  # For MDP, we initialize t=1 based on the first observation likelihood
  pi_init <- rep(1/nstates, nstates) # Or use stationary distribution of trans_matrix
  
  for (j in 1:nstates) {
    # Note: We cannot calculate transition density for t=1 without t=0.
    # We assume a steady state density or simple initialization. 
    # Here we initialize based on log(pi) assuming likelihood is handled in steps t > 1
    value_func[j, 1] <- log(pi_init[j]) 
  }
  
  # --- 2. Recursion (The Forward Pass) ---
  
  for (t in 2:T_len) {
    for (j in 1:nstates) { # Current State j
      
      # A. Calculate Transition Rewards from all possible previous states i
      # transition_rewards[i] = V_{t-1}(i) + log(Transition i->j)
      prev_values <- value_func[, t-1]
      trans_probs <- log(trans_matrix[, j]) # Column j: prob of moving TO j
      
      transition_rewards <- prev_values + trans_probs
      
      # B. Maximization Step (The Bellman Optimality)
      # Find the best previous state 'i' that maximizes the path to 'j'
      best_prev_state <- which.max(transition_rewards)
      max_prev_val <- transition_rewards[best_prev_state]
      
      # C. Calculate Emission Reward (Approximated Conditional Expectation)
      # P(y_t | y_{t-1}, State_j) using Gaussian Approximation
      if (G_approx)
        emission_reward <- get_gaussian_approx_density(
          x_next = obs[t], 
          x_prev = obs[t-1], 
          kappa = params$kappa[j], 
          theta = params$theta[j], 
          sigma = params$sigma[j], 
          dt = dt
        )else{
          emission_reward <- ln_d_Heston(
            obs[t], 
            obs[t-1], 
            k = 1, 
            params$kappa[j], 
            params$theta[j], 
            params$sigma[j])
        }
      # D. Update Value Function
      value_func[j, t] <- max_prev_val + emission_reward
      
      # E. Store Backpointer
      backpointer[j, t] <- best_prev_state
    }
  }
  
  # --- 3. Termination & Backtracking (The Backward Pass) ---
  
  # Find the best final state
  best_path <- numeric(T_len)
  best_path[T_len] <- which.max(value_func[, T_len])
  
  # Backtrack
  for (t in T_len:2) {
    best_path[t-1] <- backpointer[best_path[t], t]
  }
  
  return(best_path)
}


params <- list(
  kappa = c(2, 1),   # Fast vs Slow mean reversion
  theta = c(0.03, 0.6), # Means are closer together (Challenging)
  sigma = c(0.05, 0.05)    # Different volatilities
)

# params <- list(
#   kappa = c(6.120589e-05, 6.634052e-02 ),   # Fast vs Slow mean reversion
#   theta = c(0.03421763, 0.19811660 ), # Means are closer together (Challenging)
#   sigma = c(0.005357511, 0.018877040 )    # Different volatilities
# )

# 
# mu <- c(0.5, 0.5)
# kappa <- c(2, 1)
# theta <- c(0.03, 0.6)
# sigma <- c(0.05, 0.05)
# rho <- c(-0.1, -0.1)


# Run Viterbi
estimated_states <- viterbi_mdp_heston(V_simulated, nstates=2, Gamma, params, dt=1/252, T)
estimated_states

plot(estimated_states, col = "blue")
lines(Reg_chain+1)



