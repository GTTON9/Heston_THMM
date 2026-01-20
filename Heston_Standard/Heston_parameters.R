

Heston_parameters <- function(
    controls = list(), 
    hierarchy = FALSE, 
    states = 2, 
    sdds = "Heston", 
    Gamma = NULL, 
    # Heston/CIR Parameters
    kappa = NULL, theta = NULL, sigma = NULL, 
 
    scale_par = c(1), 
    seed = NULL,
    check_controls = TRUE
) {
  
  ### check 'controls' and 'scale_par'
  if (isTRUE(check_controls)) {

    controls <- Heston_set_controls(
      controls = controls, hierarchy = FALSE, states = states, sdds = sdds
    )
  } else {
    controls <- structure(
      oeli::merge_lists(
        controls,

        list("hierarchy" = FALSE, "states" = states[1], "sdds" = sdds[1]) 
      ),
      class = "Heston_controls"
    )
  }
  

  if (!checkmate::test_numeric(scale_par[1], len = 1, lower = 0)) {
    stop(
      "'scale_par' must be a positive numeric vector of length at least 1.",
      call. = FALSE
    )
  }
  scale_par_c <- scale_par[1] 
  
  ### extract specifications
  M <- controls[["states"]][1] 
  sdds <- controls[["sdds"]]
  
  ### set seed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  ### specify missing parameters
  if (is.null(Gamma)) {
    print("random param")
    Gamma <- oeli::sample_transition_probability_matrix(
      dim = M, state_persistent = TRUE
    )
    colnames(Gamma) <- rownames(Gamma) <- paste0("state_", seq_len(M))
  }
  

  if (is.null(kappa)) {
    kappa <- stats::qunif((0:(M - 1) / M + stats::runif(1, 0, 1 / M)), 0.5, 3) * scale_par_c 
    kappa <- sort(kappa, decreasing = TRUE) 
  }

  if (is.null(theta)) {
    theta <- stats::qunif((0:(M - 1) / M + stats::runif(1, 0, 1 / M)), 0.01, 1) * scale_par_c
  }
  
  if (is.null(sigma)) {
    sigma <- stats::qunif((0:(M - 1) / M + stats::runif(1, 0, 1 / M)), 0.1, 1.5) * scale_par_c
  }
  # print("------------------Parameter Update------------------")
  # print(Gamma)
  # print(kappa)
  # print(theta)
  # print(sigma)
  
  # --- Standard HMM parameters (mu, df) are explicitly set to NULL ---
  mu <- NULL
  df <- NULL
  
  # --- Hierarchy parameters (*_star) are explicitly set to NULL ---
  Gamma_star <- mu_star <- sigma_star <- df_star <- NULL
  
  ### set fixed parameters (if any)
  sdd_pars <- sdds[[1]]$pars
  
  if ("kappa" %in% names(sdd_pars)) {
    kappa <- sdd_pars$kappa
    if (length(kappa) == 1) {
      kappa <- rep(kappa, M)
    }
  }
  
  if ("theta" %in% names(sdd_pars)) {
    theta <- sdd_pars$theta
    if (length(theta) == 1) {
      theta <- rep(theta, M)
    }
  }
  
  if ("sigma" %in% names(sdd_pars)) {
    sigma <- sdd_pars$sigma
    if (length(sigma) == 1) {
      sigma <- rep(sigma, M)
    }
  }
  
  ### check parameters
  oeli::assert_transition_probability_matrix(Gamma, dim = M)
  

  if (!checkmate::test_numeric(kappa, len = M, lower = 0, finite = TRUE)) {
    stop(
      paste("'kappa' must be a positive numeric vector of length", M),
      call. = FALSE
    )
  }
  

  if (!checkmate::test_numeric(theta, len = M, lower = 0, finite = TRUE)) {
    stop(
      paste("'theta' must be a positive numeric vector of length", M),
      call. = FALSE
    )
  }
  

  if (!checkmate::test_numeric(sigma, len = M, lower = 0, finite = TRUE)) {
    stop(
      paste("'sigma' must be a positive numeric vector of length", M),
      call. = FALSE
    )
  }
  
  
  out <- list("sdds" = sdds)
  
  if (!is.null(Gamma)) {
    colnames(Gamma) <- rownames(Gamma) <- paste0("state_", seq_len(M))
    out <- c(out, list("Gamma" = Gamma))
  }
  

  if (!is.null(kappa)) {
    names(kappa) <- paste0("kappaCon_", seq_along(kappa))
    out <- c(out, list("kappa" = kappa))
  }
  if (!is.null(theta)) {
    names(theta) <- paste0("thetaCon_", seq_along(theta))
    out <- c(out, list("theta" = theta))
  }
  if (!is.null(sigma)) {
    names(sigma) <- paste0("sigmaCon_", seq_along(sigma))
    out <- c(out, list("sigma" = sigma))
  }
  
  structure(out, class = c("Heston_parameters", "list"))
}




