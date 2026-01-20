source("Heston_likelihood.R")
#' Decode the underlying hidden state sequence
#'
#' @description
#' This function decodes the (most likely) underlying hidden state sequence by
#' applying the Viterbi algorithm for global decoding.
#'
#' @references
#' <https://en.wikipedia.org/wiki/Viterbi_algorithm>
#'
#' @param x
#' An object of class \code{\link{fHMM_model}}.
#' @param verbose
#' Set to \code{TRUE} to print progress messages.
#'
#' @return
#' An object of class \code{\link{fHMM_model}} with decoded state sequence 
#' included.
#'
#' @export
#'
#' @examples
#' decode_states(dax_model_3t)
#' plot(dax_model_3t, type = "ts")

decode_states_heston <- function(x, verbose = TRUE) {
  
  ### check input
  if (!inherits(x,"fHMM_model")) {
    stop("'x' must be of class 'fHMM_model'.", call. = FALSE)
  }
  if (!isTRUE(verbose) && !isFALSE(verbose)) {
    stop("'verbose' must be either TRUE or FALSE.", call. = FALSE)
  }
  
  ### apply Viterbi algorithm
  par <- parUncon2par_heston(x$estimate, x$data$controls)

  # decoding <- viterbi(
  #   observations = x$data$data,
  #   nstates = x$data$controls$states[1],
  #   Gamma = par$Gamma, kappa = par$kappa,
  #   theta = par$theta, sigma = par$sigma)

  decoding <- forward_filter(x$data$data, x$data$controls$states[1], Gamma = par$Gamma, kappa = par$kappa,
                 theta = par$theta, sigma = par$sigma)$states_estimate
    
  ### save decoding in 'x' and return 'x'
  if (verbose) message("Decoded states")
  x$decoding <- decoding
  return(x)
}




#' @rdname decode_states
#' @param observations
#' A \code{numeric} \code{vector} of state-dependent observations.
#' @param nstates
#' The number of states.
#' @param sdd
#' A \code{character}, specifying the state-dependent distribution. One of 
#' \itemize{
#'   \item \code{"normal"} (the normal distribution),
#'   \item \code{"lognormal"} (the log-normal distribution),
#'   \item \code{"t"} (the t-distribution),
#'   \item \code{"gamma"} (the gamma distribution),
#'   \item \code{"poisson"} (the Poisson distribution).
#' }
#' @param Gamma
#' A transition probability \code{matrix} of dimension \code{nstates}.
#' @param mu
#' A \code{numeric} vector of expected values for the state-dependent 
#' distribution in the different states of length \code{nstates}.
#' 
#' For the gamma- or Poisson-distribution, \code{mu} must be positive.
#' 
#' @param sigma
#' A positive \code{numeric} vector of standard deviations for the 
#' state-dependent distribution in the different states of length \code{nstates}. 
#' 
#' Not relevant in case of a state-dependent Poisson distribution.
#' 
#' @param df
#' A positive \code{numeric} vector of degrees of freedom for the 
#' state-dependent distribution in the different states of length \code{nstates}. 
#' 
#' Only relevant in case of a state-dependent t-distribution.
#' 
#' @examples
#' viterbi(
#'   observations = c(1, 1, 1, 10, 10, 10),
#'   nstates      = 2,
#'   sdd          = "poisson",
#'   Gamma        = matrix(0.5, 2, 2),
#'   mu           = c(1, 10)
#' )
#' 
#' @export

viterbi <- function(
    observations, nstates, Gamma, kappa, theta, sigma, mu = NULL, rho = NULL, G_approx = TRUE
) {

  T <- length(observations)-1
  delta <- oeli::stationary_distribution(Gamma, soft_fail = TRUE)
  # allprobs <- matrix(0, nstates, T)
  allprobs <- matrix(NA_real_, nstates, T)
  for (i in seq_len(nstates)) {

    allprobs[i, ] <- get_transition_density_heston_ln(observations, kappa[i], theta[i], sigma[i])
  }

  xi <- matrix(0, nstates, T)
  
  for (n in seq_len(nstates)) { 
    xi[n, 1] <- delta[n] + allprobs[n, 1]
  }
  for (t in seq_len(T)[-1]) {
    for (n in seq_len(nstates)) {
      xi[n, t] <- max(xi[, t - 1] + Gamma[, n]) + allprobs[n, t]
    }
  }

  iv <- numeric(T)
  iv[T] <- which.max(xi[, T])
  for (t in rev(seq_len(T - 1))) {
    iv[t] <- which.max(xi[, t] + Gamma[, iv[t + 1]])
  }
  
  
  return(iv)
}



forward_filter <- function(observations, nstates, Gamma, kappa, theta, sigma) {
  
  T <- length(observations) - 1
  logGamma <- log(Gamma)
  delta <- log(oeli::stationary_distribution(Gamma))
  
  # emission log-prob
  allprobs <- matrix(0, nstates, T)
  for (i in 1:nstates) {
    allprobs[i,] <- get_transition_density_heston_ln(observations, kappa[i], theta[i], sigma[i])
  }
  
  logalpha <- matrix(-Inf, nstates, T)
  
  # initialization
  logalpha[,1] <- delta + allprobs[,1]
  
  # recursion
  for (t in 2:T) {
    for (j in 1:nstates) {
      logalpha[j,t] <- logsumexp(logalpha[,t-1] + logGamma[,j]) + allprobs[j,t]
    }
  }
  
  # normalized posterior: p(s_t=j | y_1:t)
  posterior <- matrix(0, nstates, T)
  reg <- matrix(0, 1, T)
  for (t in 1:T) {
    c <- logsumexp(logalpha[,t])
    posterior[,t] <- exp(logalpha[,t] - c)
    reg[,t] <- which.max(posterior[,t])
  }
  
  return(list(posterior = posterior,
              states_estimate = as.vector(reg)))
}

logsumexp <- function(x) {
  
  m <- max(x)
  
  m + log(sum(exp(x - m)))
  
}






