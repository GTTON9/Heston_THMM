

#' Negative log-likelihood function of an HMM
#'
#' @description
#' This function computes the negative log-likelihood of an HMM.
#'
#' @param parUncon
#' An object of class \code{parUncon}.
#' @param observations
#' The vector of the simulated or empirical data used for estimation.
#' @param controls
#' An object of class \code{Heston_controls}.
#'
#' @return
#' The negative log-likelihood value.
#'
#' @keywords internal

# install.packages("Bessel")
# library(Bessel)
Heston_nLL_hmm <- function(parUncon, observations, controls) {

  class(parUncon) <- "parUncon"
  T <- length(observations)
  
  nstates <- controls[["states"]][1]

  par <- parUncon2par_heston(parUncon, controls, FALSE, numerical_safeguard = TRUE)
  
  
  # print("----- likelihood par-----")
  # print(par)
  sdd <- controls[["sdds"]][[1]]$name

  Gamma <- par[["Gamma"]]

  delta <- try(
    solve(t(diag(nstates) - Gamma + 1), rep(1, nstates)),
    silent = TRUE
  )

  if (inherits(delta, "try-error")) {
    delta <- rep(1, nstates) / nstates
  }
  
  kappa <- par[["kappa"]]
  theta <- par[["theta"]]
  sigma <- par[["sigma"]]
  # kappa <- c(10,5)

  allprobs <- matrix(NA_real_, nstates, T-1)
  for (i in 1:nstates) {
    if (sdd == "Heston") {
        
        allprobs[i, ] <- get_transition_density_heston_ln(observations, kappa[i], theta[i], sigma[i])
        } 
    
    else {
      stop("Unknown state-dependent distribution", call. = FALSE)
    }
  }
  # print(allprobs)
  
  ll <- - LL_HMM_R(allprobs, Gamma, delta)

  # print(paste("negative_log_likelihood:",ll))

  return(ll)
  # -LL_HMM_Rcpp(allprobs, Gamma, delta, nstates, T)

}
get_transition_density_Gaussian <- function(observations, kappa, theta, sigma){
  T_obs <- length(observations)
  all_ln_densities <- numeric(T_obs - 1) 
  
  for (t in 1:(T_obs - 1)) {
    V_t <- observations[t]
    V_t_plus_k <- observations[t + 1]
    transition_density <- get_gaussian_approx_density(
      V_t, V_t_plus_k, 
      kappa = kappa, theta = theta, sigma = sigma, 
      dt = 1/252)
    
    all_ln_densities[t] <- transition_density
  }
  return(all_ln_densities)
}


get_transition_density_heston_ln <- function(observations, kappa, theta, sigma) {
  T_obs <- length(observations)

  all_ln_densities <- numeric(T_obs - 1) 
  
  for (t in 1:(T_obs - 1)) {
 
    V_t <- observations[t]
    V_t_plus_k <- observations[t + 1]
    
    ln_transition_density <- ln_d_Heston(
      V_t = V_t, 
      V_t_plus_k = V_t_plus_k, 
      k = 1/252, kappa = kappa, theta = theta, sigma = sigma
    )
    
    all_ln_densities[t] <- (ln_transition_density)
  }

  return(all_ln_densities)
}


# 
# ln_d_Heston <- function(V_t, V_t_plus_k, k = 1/252, kappa, theta, sigma) {
#   
#   # --- 1. 参数计算 ---
#   sigma <- max(sigma, 1e-8)
#   exp_k <- exp(-kappa * k)
#   den_factor <- (1 - exp_k) * (sigma^2)
#   
#   if(abs(den_factor) < 1e-10) return(-Inf) 
#   
#   C <- 2 * kappa / den_factor
#   q <- (2 * kappa * theta) / (sigma^2) - 1 # 
#   u <- C * V_t * exp_k
#   v <- C * V_t_plus_k
#   z_arg <- 2 * sqrt(u * v) 
#   
#   bessel_val <- tryCatch({
#     BesselI(z_arg, nu = q, expon.scaled = TRUE) 
#   }, error = function(e) {
#     besselI(z_arg, nu = q, expon.scaled = TRUE)
#   })
#   
#   
#   if (is.finite(bessel_val) && bessel_val > 0) {
#     ln_I_q <- z_arg + log(bessel_val)
#     
#   } else {
#     
#     
#     q_abs <- abs(q)
#     
#     if (z_arg > 100 && q_abs < z_arg / 5) {
#       # besselIasym for large z regular q
#       ln_I_q <- tryCatch({
#         besselIasym(z_arg, nu = q, log = TRUE)
#       }, error = function(e) {
#         # DLMF 10.40.1: I_q(z) ≈ exp(z) / sqrt(2*pi*z) (z -> ∞)
#         warning(paste("besselIasym failed for z =", z_arg, 
#                       ". Switching to large-z asymptotic approximation (DLMF 10.40.1). Original error:", 
#                       conditionMessage(e)))
#         z_arg - 0.5 * log(2 * pi * z_arg)
#       })
#       
#       
#     } else if (q_abs > 10 && q_abs > z_arg * 1.2) {
#       # for large q regular z
#       # besselI.nuAsym (nu >= 0)
#       ln_I_q <- tryCatch({
#         besselI.nuAsym(z_arg, nu = q_abs, log = TRUE, k.max = 5)
#       }, error = function(e) {
#         # q * (1 + log(z / (2 * q))) - 0.5 * log(2 * pi * q)
#         # DLMF 10.41.1
#         warning(paste("besselI.nuAsym failed for q =", q_abs, 
#                       ". Switching to asymptotic approximation (DLMF 10.41.1). Original error:", 
#                       conditionMessage(e)))
#         
#         q_abs * (1 + log(z_arg / (2 * q_abs))) - 0.5 * log(2 * pi * q_abs)
#       })
#       
#     } else {
#       # ratio_abs <- z_arg / q_abs
#       # r <- sqrt(1 + ratio_abs^2)
#       # # eta = r + log(ratio_abs / (1 + r))
#       # eta <- r + log(ratio_abs) - log(1 + r)
#       # 
#       # # ln_I_q ≈ |q| * eta - 0.5 * log(2 * pi * |q| * r)
#       # ln_I_q <- q_abs * eta - 0.5 * log(2 * pi * q_abs * r)
#       # 
#       # 
#       # 
#       mu_approx <- V_t + kappa * (theta - V_t) * k
# 
#       # Variance: Var[V_{t+k} | V_t] ≈ diffusion^2 * k = (sigma^2 * V_t) * k
#       variance_approx <- (sigma^2) * V_t * k
# 
#       # Ensure standard deviation is positive
#       sd_approx <- sqrt(pmax(variance_approx, 1e-8))
# 
#       # Calculate the log density using the Normal distribution
#       ln_density <- dnorm(V_t_plus_k, mean = mu_approx, sd = sd_approx, log = TRUE)
#       return(ln_density)
#     }
#   }
# 
#   log_vu_ratio <- 0
#   if(u > 1e-10 && v > 1e-10) {
#     log_vu_ratio <- log(v / u)
#   }
#   
#   ln_density <- log(C) - u - v + (q / 2) * log_vu_ratio + ln_I_q
#   
#   return(ln_density)
# }



ln_d_Heston <- function(V_t, V_t_plus_k, 
                                   k = 1/252, 
                                   kappa, theta, sigma,
                                   min_val = 1e-300) {
  # ---- 1. Compute Constants ----
  ek <- exp(-kappa * k)
  
  c  <- (sigma^2 * (1 - ek)) / (4 * kappa)
  d  <- 4 * kappa * theta / sigma^2
  lambda <- (4 * kappa * ek / (sigma^2 * (1 - ek))) * V_t
  
  # Avoid numerical blowups
  c      <- max(c,  min_val)
  d      <- max(d,  min_val)
  lambda <- max(lambda, min_val)
  V_ratio <- V_t_plus_k / c
  
  # ---- 2. Compute noncentral chi-square density ----
  dens <- dchisq(V_ratio, df = d, ncp = lambda, log = FALSE)
  
  # Adjust for variable transformation (divide by c)
  return(log(dens / c))
}
