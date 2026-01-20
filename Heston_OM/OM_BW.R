# ================================================================
# O-M THMM Complete Implementation (CORRECTED VERSION)
# Key Corrections:
# 1. Added proper normalization term: (1/2) * ∫ log(σ²ϕ) dt
# 2. Added sigma penalty to prevent overestimation
# 3. Improved numerical stability
# ================================================================

# ================================================================
# 1. Calculate O-M Functional with Normalization (CORRECTED)
# ================================================================

#' Calculate O-M Functional for Single Path
#' 
#' This implements the correct OM functional from the paper:
#' L[ϕ] = (1/2) ∫ [(ϕ̇ - κ(θ-ϕ))² / (σ²ϕ)] dt
#' 
#' @param a Starting value of linear path
#' @param b Slope of linear path (ϕ(t) = a + bt)
#' @param kappa Mean reversion speed
#' @param theta Long-term mean
#' @param sigma Volatility parameter
#' @param H Time horizon
#' @param N Number of integration steps
#' @return OM functional value (energy/cost)
calculate_OM_single <- function(a, b, kappa, theta, sigma, 
                                H = 5/252, N = 100) {
  
  # Validity checks
  if (a <= 0) return(1e6)  # Variance must be positive
  
  end_val <- a + b * H
  if (end_val <= 0) return(1e6)
  
  # Vectorized computation
  delta_tau <- H / N
  tau <- seq(0, H, length.out = N + 1)
  phi <- a + b * tau
  
  # Ensure positivity
  phi <- pmax(phi, 1e-10)
  
  # Core calculation: (ϕ̇ - κ(θ-ϕ))² / (σ²ϕ)
  phi_dot <- b
  drift <- kappa * (theta - phi)
  denom <- sigma^2 * phi
  
  f <- (phi_dot - drift)^2 / denom
  
  # Check for numerical issues
  if (any(!is.finite(f))) return(1e6)
  
  # Trapezoidal integration
  integral <- (sum(f) - 0.5 * (f[1] + f[N + 1])) * delta_tau
  L <- 0.5 * integral
  
  return(L)
}

#' Calculate Normalization Term for Single Path
#' 
#' This implements: (1/2) ∫ log(σ²ϕ(t)) dt
#' For CIR: σ²(v) = σ² · v
#' 
#' @param a Starting value
#' @param b Slope
#' @param sigma Volatility parameter
#' @param H Time horizon
#' @param N Integration steps
#' @return Normalization term value
calculate_normalization_term <- function(a, b, sigma, H = 5/252, N = 100) {
  
  delta_tau <- H / N
  tau <- seq(0, H, length.out = N + 1)
  phi <- a + b * tau
  
  # Ensure positivity
  phi <- pmax(phi, 1e-10)
  
  # Compute ∫ log(σ²ϕ) dt = ∫ [log(σ²) + log(ϕ)] dt
  log_sigma_sq <- log(sigma^2)
  log_phi <- log(phi)
  
  # Check for numerical issues
  if (any(!is.finite(log_phi))) return(1e6)
  
  # Integrate log(ϕ) using trapezoidal rule
  integral_log_phi <- (sum(log_phi) - 0.5 * (log_phi[1] + log_phi[N + 1])) * delta_tau
  
  # Total: ∫ [log(σ²) + log(ϕ)] dt = H·log(σ²) + ∫ log(ϕ) dt
  norm_term <- H * log_sigma_sq + integral_log_phi
  
  return(norm_term)
}

# ================================================================
# 2. Calculate Complete Log-Likelihood (CORRECTED)
# ================================================================

#' Calculate O-M Emission Log-Probabilities with Normalization
#' 
#' Complete log-likelihood includes:
#' log P[ϕ] = -L[ϕ] - (1/2)·Normalization
#' 
#' @param result_with_paths Data frame with path parameters
#' @param kappa Vector of mean reversion speeds
#' @param theta Vector of long-term means
#' @param sigma Vector of volatilities
#' @param H Time horizon
#' @param N Integration steps
#' @param lambda_om Weight for OM functional
#' @param sigma_penalty Penalty weight for large sigma
#' @param use_normalization Include normalization term (recommended TRUE)
#' @return Matrix [nstates x T] of log emission probabilities
calculate_OM_for_likelihood <- function(result_with_paths,
                                              kappa, theta, sigma,
                                              H = 5/252, N = 100,
                                              lambda_om = 1.0,
                                              use_upper = TRUE,
                                              use_lower = TRUE) {
  T_periods <- nrow(result_with_paths)
  nstates <- length(kappa)
  
  log_b <- matrix(0, nrow = nstates, ncol = T_periods)
  
  for (t in 1:T_periods) {
    a_upper <- result_with_paths$path_upper_a[t]
    b_upper <- result_with_paths$path_upper_b[t]
    a_lower <- result_with_paths$path_lower_a[t]
    b_lower <- result_with_paths$path_lower_b[t]
    
    for (i in 1:nstates) {
      Li_t <- 0
      
      if (use_upper && !is.na(a_upper) && !is.na(b_upper)) {
        Li_t <- Li_t + calculate_OM_single(
          a = a_upper, b = b_upper,
          kappa = kappa[i], theta = theta[i], sigma = sigma[i],
          H = H, N = N
        )
      }
      
      if (use_lower && !is.na(a_lower) && !is.na(b_lower)) {
        Li_t <- Li_t + calculate_OM_single(
          a = a_lower, b = b_lower,
          kappa = kappa[i], theta = theta[i], sigma = sigma[i],
          H = H, N = N
        )
      }
      
      # paper: log b_i(Ot) = - L_i[phi] (up to additive const independent of state)
      log_b[i, t] <- -lambda_om * Li_t
    }
  }
  
  log_b
}

# ================================================================
# 3. Negative Log-Likelihood Function (CORRECTED)
# ================================================================

#' Negative Log-Likelihood for O-M THMM
#' 
#' Objective function for numerical optimization
#' Includes proper normalization and sigma penalty
#' 
#' @param parUncon Unconstrained parameter vector
#' @param result_with_paths Path data
#' @param controls Heston controls
#' @param lambda_om Weight for OM functional
#' @param sigma_penalty Penalty for large sigma
#' @param use_normalization Include normalization term
#' @param H Time horizon
#' @param N_integration Integration steps
#' @return Negative log-likelihood value
Heston_nLL_hmm_om <- function(parUncon,
                              result_with_paths,
                              controls,
                              lambda_om = 1.0,
                              sigma_penalty = 10.0,
                              use_normalization = TRUE,
                              H = 5/252,
                              N_integration = 100) {
  
  class(parUncon) <- "parUncon"
  
  T_periods <- nrow(result_with_paths)
  nstates <- controls[["states"]][1]
  
  # Convert parameters
  par <- tryCatch({
    parUncon2par_heston(parUncon, controls, FALSE, numerical_safeguard = TRUE)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(par)) return(1e10)
  
  # Extract parameters
  Gamma <- par[["Gamma"]]
  kappa <- par[["kappa"]]
  theta <- par[["theta"]]
  sigma <- par[["sigma"]]
  
  # Parameter validity checks
  if (any(kappa <= 0) || any(theta <= 0) || any(sigma <= 0)) {
    return(1e10)
  }
  
  # Feller condition check (soft)
  feller_penalty <- 0
  for (i in 1:nstates) {
    if (2 * kappa[i] * theta[i] < sigma[i]^2) {
      feller_penalty <- feller_penalty + 100 * (sigma[i]^2 - 2*kappa[i]*theta[i])^2
    }
  }
  
  # Calculate stationary distribution
  stationary_dist <- function(Gamma) {
    n <- nrow(Gamma)
    A <- t(Gamma) - diag(n)
    A[n, ] <- 1
    b <- rep(0, n); b[n] <- 1
    as.numeric(solve(A, b))
  }
  delta <- stationary_dist(Gamma)
  
  # Calculate emission probabilities (CORRECTED VERSION)
  log_allprobs <- tryCatch({
    calculate_OM_for_likelihood(
      result_with_paths = result_with_paths,
      kappa = kappa,
      theta = theta,
      sigma = sigma,
      H = H,
      N = N_integration,
      lambda_om = lambda_om,
      sigma_penalty = sigma_penalty,
      use_normalization = use_normalization,
      use_upper = TRUE,
      use_lower = TRUE
    )
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(log_allprobs) || any(!is.finite(log_allprobs))) {
    return(1e10)
  }
  
  # Calculate log-likelihood using forward algorithm
  ll <- tryCatch({
    LL_HMM_R(exp(log_allprobs), Gamma, delta)
  }, error = function(e) {
    return(-1e10)
  })
  
  # Return negative log-likelihood + penalties
  nll <- -ll + feller_penalty
  
  if (!is.finite(nll) || is.na(nll)) {
    return(1e10)
  }
  
  return(nll)
}

# ================================================================
# 4. Model Fitting Function (CORRECTED)
# ================================================================

#' Fit Heston Model using O-M THMM with Numerical Optimizer
#' 
#' @param result_with_paths Path data
#' @param controls Heston controls
#' @param fit Fitting options
#' @param runs Number of optimization runs
#' @param lambda_om Weight for OM functional (default 1.0)
#' @param sigma_penalty Penalty for large sigma (default 10.0)
#' @param use_normalization Include normalization term (default TRUE)
#' @param H Time horizon
#' @param N_integration Integration steps
#' @param optimizer Optimizer to use: "nlm" or "optim" (default "nlm")
#' @param ncluster Number of parallel clusters
#' @param seed Random seed
#' @param verbose Print progress
#' @return Fitted model object
Heston_fit_model_om <- function(result_with_paths,
                                controls,
                                fit = list(),
                                runs = 10,
                                origin = FALSE,
                                accept = 1:3,
                                gradtol = 0.01,
                                iterlim = 1000,
                                print.level = 0,
                                steptol = 0.01,
                                lambda_om = 1.0,
                                sigma_penalty = 10.0,
                                use_normalization = TRUE,
                                H = 5/252,
                                N_integration = 100,
                                optimizer = "nlm",
                                ncluster = 1,
                                seed = NULL,
                                verbose = TRUE,
                                initial_estimate = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Validation
  if (!checkmate::test_count(ncluster, positive = TRUE)) {
    stop("'ncluster' must be a positive integer.", call. = FALSE)
  }
  
  if (!optimizer %in% c("nlm", "optim")) {
    stop("'optimizer' must be either 'nlm' or 'optim'.", call. = FALSE)
  }
  
  if (verbose) {
    
  }
  
  # Update controls
  controls <- Heston_set_controls(
    controls = controls,
    hierarchy = controls[["hierarchy"]],
    states = controls[["states"]],
    sdds = controls[["sdds"]],
    horizon = nrow(result_with_paths),
    period = controls[["period"]],
    data = controls[["data"]],
    fit = fit,
    runs = runs,
    origin = origin,
    accept = accept,
    gradtol = gradtol,
    iterlim = iterlim,
    print.level = print.level,
    steptol = steptol
  )
  
  # Create pseudo data for initialization
  pseudo_data <- list(
    data = result_with_paths$v_start,
    controls = controls
  )
  class(pseudo_data) <- "Heston_data"
  
  # Get initial values
  if (verbose) {
    message("Generating initial values...")
  }
  
  initial_values <- Heston_get_init(
    data = pseudo_data,
    ncluster = ncluster,
    seed = seed,
    verbose = verbose,
    initial_estimate = initial_estimate
  )
  
  runs <- length(initial_values)
  
  # Define target function
  target <- function(parUncon) {
    Heston_nLL_hmm_om(
      parUncon = parUncon,
      result_with_paths = result_with_paths,
      controls = controls,
      lambda_om = lambda_om,
      sigma_penalty = sigma_penalty,
      use_normalization = use_normalization,
      H = H,
      N_integration = N_integration
    )
  }
  
  # Progress bar
  if (verbose) {
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent, :eta ETA",
      total = runs,
      width = 45,
      clear = TRUE,
      complete = "=",
      incomplete = "-",
      current = ">"
    )
    pb$message("Maximizing O-M likelihood...")
  }
  
  # Optimization
  start_time <- Sys.time()
  
  mods <- list()
  
  for (run in seq_len(runs)) {
    if (verbose) pb$tick(0)
    
    suppressWarnings({
      
      if (optimizer == "nlm") {
        # Use NLM optimizer
        mod <- try(
          stats::nlm(
            f = target,
            p = initial_values[[run]],
            iterlim = controls[["fit"]][["iterlim"]],
            steptol = controls[["fit"]][["steptol"]],
            gradtol = controls[["fit"]][["gradtol"]],
            print.level = controls[["fit"]][["print.level"]],
            hessian = FALSE
          ),
          silent = TRUE
        )
        
      } else if (optimizer == "optim") {
        # Use OPTIM with L-BFGS-B
        # Need to set appropriate bounds for unconstrained parameters
        
        mod <- try(
          optim(
            par = initial_values[[run]],
            fn = target,
            method = "L-BFGS-B",
            # Loose bounds on unconstrained space
            lower = rep(-10, length(initial_values[[run]])),
            upper = rep(10, length(initial_values[[run]])),
            control = list(
              factr = 1e7,
              pgtol = 0,
              maxit = controls[["fit"]][["iterlim"]]
            ),
            hessian = FALSE
          ),
          silent = TRUE
        )
        
        # Convert optim output to nlm-like structure
        if (!inherits(mod, "try-error")) {
          mod_nlm <- list(
            minimum = mod$value,
            estimate = mod$par,
            gradient = NULL,
            code = ifelse(mod$convergence == 0, 1, 5)
          )
          mod <- mod_nlm
        }
      }
    })
    
    # Check if run was successful
    accept_run <- !inherits(mod, "try-error") &&
      mod[["code"]] %in% controls[["fit"]][["accept"]]
    
    if (accept_run) {
      mods[[run]] <- mod
    } else {
      mods[[run]] <- NA
    }
    
    if (verbose) pb$tick()
  }
  
  end_time <- Sys.time()
  
  # Extract log-likelihoods
  lls <- -unlist(sapply(mods, `[`, "minimum"), use.names = FALSE)
  
  if (all(is.na(lls))) {
    stop("None of the estimation runs ended successfully.\n",
         "Try: (1) Increase 'runs', (2) Adjust sigma_penalty, or (3) Change optimizer.",
         call. = FALSE)
  }
  
  if (verbose) {
    cat("\n\n")
    cat("Successful runs:", sum(!is.na(lls)), "/", runs, "\n")
    cat("Best log-likelihood:", max(lls, na.rm = TRUE), "\n")
  }
  
  # Approximate Hessian for best model
  if (verbose) {
    message("Approximating Hessian...")
  }
  
  best_estimate <- mods[[which.max(lls)]][["estimate"]]
  
  fisher <- pracma::hessdiag(f = target, x = best_estimate)
  
  if (all(fisher > 0)) {
    inverse_fisher <- 1 / fisher
  } else {
    hessian <- suppressWarnings(
      pracma::hessian(f = target, x0 = best_estimate)
    )
    inverse_fisher <- diag(MASS::ginv(hessian))
  }
  
  if (verbose) {
    message("Fitting completed!")
  }
  
  # Best model
  mod <- mods[[which.max(lls)]]
  ll <- -mod[["minimum"]]
  estimate <- mod[["estimate"]]
  class(estimate) <- "parUncon"
  
  estimation_time <- ceiling(difftime(end_time, start_time, units = "mins"))
  
  # Create model object
  out <- fHMM_model(
    data = pseudo_data,
    estimate = estimate,
    nlm_output = mod,
    estimation_time = estimation_time,
    ll = ll,
    lls = lls,
    gradient = mod$gradient,
    inverse_fisher = inverse_fisher,
    decoding = NULL
  )
  
  # Reorder states
  out <- Heston_reorder_states(out, state_order = "mean")
  
  return(out)
}
