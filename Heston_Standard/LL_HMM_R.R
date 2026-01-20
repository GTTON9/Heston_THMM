#' @title Log-Sum-Exp Trick
#' @description Computes log(sum(exp(x))) in a numerically stable way.
log_sum_exp <- function(x) {
  c <- max(x)
  return(c + log(sum(exp(x - c))))
}


#' @title Standard HMM Log-Likelihood (Forward Algorithm)
#' @description Calculates the log-likelihood of the observations using the Forward Algorithm.
#'
#' @param allprobs A matrix of emission probabilities (state-dependent observation probabilities),
#'        where allprobs[i, t] = P(observation_t | state=i).
#' @param Gamma The state transition probability matrix, where Gamma[j, i] = P(state_t=i | state_t-1=j).
#' @param delta The vector of initial state probabilities, delta[i] = P(state_1=i).
#' @return The total log-likelihood of the observation sequence.
LL_HMM_R <- function(allprobs, Gamma, delta) {
  
  # N: Number of states (from delta vector length)
  N <- length(delta)
  # T: Length of sequence (from number of columns in allprobs)
  T <- ncol(allprobs) 
  
  # Initialize phi matrix (log-forward variables: log P(o_1..t, s_t=i))
  # R indexing is 1-based, so phi[, 1] corresponds to t=0 in C++
  phi <- matrix(0, nrow = N, ncol = T)
  
  # --- Initialization (t = 1 in R, t=0 in C++) ---
  # phi(i, 0) = log(delta(i) * allprobs(i, 0));
  phi[, 1] <- log(delta) + allprobs[, 1]
  
  # --- Recursion (t = 2 to T in R, t=1 to T-1 in C++) ---
  # Outer loop iterates through time steps (columns)
  for (t in 2:T) {
    
    # phi[, t-1] is the log-forward vector from the previous time step
    prev_phi <- phi[, t - 1]
    
    # Inner loop iterates through the current state (i)
    for (i in 1:N) {
      
      # C++ equivalent of: phi.col(t - 1) + log(Gamma.col(i))
      # This calculates log P(o_1..t-1, s_t-1=j) + log P(s_t=i | s_t-1=j) for all j
      log_sum_terms <- prev_phi + log(Gamma[, i])
      
      # Apply Log-Sum-Exp Trick: log(sum(exp(...)))
      log_sum_exp_result <- log_sum_exp(log_sum_terms)
      
      # The final forward variable:
      # phi(i, t) = log(sum) + log(allprobs(i, t))
      phi[i, t] <- log_sum_exp_result + allprobs[i, t]
    }
  }
  
  # --- Termination (Calculate final log-likelihood) ---
  # log(sum(exp(phi.col(T - 1) - c))) + c
  final_log_likelihood <- log_sum_exp(phi[, T])
  
  return(final_log_likelihood)
}



