

# par: the original parameters of Heston Model
# Con: parameters that have certain condition
#   Gamma: transition matrix, which 0 <= P_ij <= 1, logit to 
#   kappa: long term moving average rate, kappa > 0 
#   theta: long term mean, theta >0
#.  sigma: volatility of V, sigma > 0


# Convert raw parameter to unconditional parameters
par2parUncon_heston <- function(par, controls, use_parameter_labels = TRUE) {
  
  
  stopifnot(inherits(par, "Heston_parameters"))
  stopifnot(inherits(controls, "Heston_controls"))
  
  sdds <- controls[["sdds"]]
  
  
  parUncon <- Gamma2gammasUncon(
    par[["Gamma"]], 
    prefix = "gammasUncon_",
    use_parameter_labels = use_parameter_labels
  )
  
  
  if (!"kappa" %in% names(sdds[[1]]$pars)) {
    parUncon <- c(
      parUncon,
      kappaCon2kappaUncon(
        kappaCon = par[["kappa"]],
        prefix = "kappaUncon_",
        use_parameter_labels = use_parameter_labels
      )
    )
  }
  
  
  if (!"theta" %in% names(sdds[[1]]$pars)) {
    parUncon <- c(
      parUncon,
      thetaCon2thetaUncon(
        thetaCon = par[["theta"]],
        prefix = "thetaUncon_",
        use_parameter_labels = use_parameter_labels
      )
    )
  }
  
  
  if (!"sigma" %in% names(sdds[[1]]$pars)) {
    parUncon <- c(
      parUncon,
      sigmaCon2sigmaUncon(
        sigmaCon = par[["sigma"]], 
        prefix = "sigmaUncon_",
        use_parameter_labels = use_parameter_labels
      )
    )
  }
  
  structure(parUncon, class = c("parUncon", "numeric"))
}

kappaCon2kappaUncon <- function(kappaCon, prefix = "kappaUncon_", use_parameter_labels = TRUE) {
  stopifnot(is.numeric(kappaCon))
  
  kappaUncon <- log(kappaCon)
  
  if (isTRUE(use_parameter_labels)) {
    
    names(kappaUncon) <- paste0(prefix, seq_along(kappaUncon))
  }
  return(kappaUncon)
}

#' @rdname kappaCon2kappaUncon
#' @export
thetaCon2thetaUncon <- function(thetaCon, prefix = "thetaUncon_", use_parameter_labels = TRUE) {
  stopifnot(is.numeric(thetaCon))
  
  thetaUncon <- log(thetaCon)
  
  if (isTRUE(use_parameter_labels)) {
    names(thetaUncon) <- paste0(prefix, seq_along(thetaUncon))
  }
  return(thetaUncon)
}

#' @rdname kappaCon2kappaUncon
#' @export
sigmaCon2sigmaUncon <- function(sigmaCon, prefix = "sigmaUncon_", use_parameter_labels = TRUE) {
  stopifnot(is.numeric(sigmaCon))
  
  sigmaUncon <- log(sigmaCon)
  
  if (isTRUE(use_parameter_labels)) {
    names(sigmaUncon) <- paste0(prefix, seq_along(sigmaUncon))
  }
  return(sigmaUncon)
}







# --- Heston Parameter Forward Conversions (Constrained -> Unconstrained) ---

#' @rdname heston_parameter_transformations
#' @export
kappaCon2kappaUncon <- function(kappaCon, prefix = "kappaUncon_", use_parameter_labels = TRUE) {
  kappaUncon <- log(kappaCon)
  if (isTRUE(use_parameter_labels)) {
    names(kappaUncon) <- paste0(prefix, seq_along(kappaUncon))
  }
  return(kappaUncon)
}

#' @rdname heston_parameter_transformations
#' @export
thetaCon2thetaUncon <- function(thetaCon, prefix = "thetaUncon_", use_parameter_labels = TRUE) {
  thetaUncon <- log(thetaCon)
  if (isTRUE(use_parameter_labels)) {
    names(thetaUncon) <- paste0(prefix, seq_along(thetaUncon))
  }
  return(thetaUncon)
}

#' @rdname heston_parameter_transformations
#' @export
sigmaCon2sigmaUncon_heston <- function(sigmaCon, prefix = "sigmaUncon_", use_parameter_labels = TRUE) {
  sigmaUncon <- log(sigmaCon)
  if (isTRUE(use_parameter_labels)) {
    names(sigmaUncon) <- paste0(prefix, seq_along(sigmaUncon))
  }
  return(sigmaUncon)
}

# --- Heston Parameter Inverse Conversions (Unconstrained -> Constrained) ---

#' @rdname heston_parameter_transformations
#' @export
kappaUncon2kappaCon <- function(kappaUncon, prefix = "kappaCon_", use_parameter_labels = TRUE) {
  kappaCon <- exp(kappaUncon)
  if (isTRUE(use_parameter_labels)) {
    names(kappaCon) <- paste0(prefix, seq_along(kappaCon))
  }
  return(kappaCon)
}

#' @rdname heston_parameter_transformations
#' @export
thetaUncon2thetaCon <- function(thetaUncon, prefix = "thetaCon_", use_parameter_labels = TRUE) {
  thetaCon <- exp(thetaUncon)
  if (isTRUE(use_parameter_labels)) {
    names(thetaCon) <- paste0(prefix, seq_along(thetaCon))
  }
  return(thetaCon)
}

#' @rdname heston_parameter_transformations
#' @export
sigmaUncon2sigmaCon_heston <- function(
    sigmaUncon, prefix = "sigmaCon_", use_parameter_labels = TRUE,
    numerical_safeguard = FALSE
) {
  sigmaCon <- exp(sigmaUncon)
  if (isTRUE(numerical_safeguard)) {
    
    sigmaCon[sigmaCon > 100] <- 100 
  }
  if (isTRUE(use_parameter_labels)) {
    names(sigmaCon) <- paste0(prefix, seq_along(sigmaCon))
  }
  return(sigmaCon)
}


# --- Heston Main Conversion Wrappers ---

#' @rdname heston_parameter_transformations
#' @export
parUncon2parCon_heston <- function(
    parUncon, controls, use_parameter_labels = TRUE, numerical_safeguard = FALSE
) {
  stopifnot(inherits(parUncon, "parUncon"))
  stopifnot(inherits(controls, "Heston_controls"))
  # print("p2")
  M <- controls[["states"]][1]
  
  
  Gamma_size <- (M - 1) * M

  parCon <- gammasUncon2gammasCon(
    parUncon[1:Gamma_size], 
    dim = M,
    prefix = "gammasCon_",
    use_parameter_labels = use_parameter_labels,
    numerical_safeguard = numerical_safeguard
  )
  parUncon <- parUncon[-(1:Gamma_size)]

  
  kappa_size <- M
  parCon <- c(
    parCon,
    kappaUncon2kappaCon(
      parUncon[1:kappa_size],
      prefix = "kappaCon_",
      use_parameter_labels = use_parameter_labels
    )
  )

  parUncon <- parUncon[-(1:kappa_size)]
  
  theta_size <- M
  parCon <- c(
    parCon,
    thetaUncon2thetaCon(
      parUncon[1:theta_size],
      prefix = "thetaCon_",
      use_parameter_labels = use_parameter_labels
    )
  )
  parUncon <- parUncon[-(1:theta_size)]
  
  sigma_size <- M
  parCon <- c(
    parCon,
    sigmaUncon2sigmaCon_heston(
      parUncon[1:sigma_size],
      prefix = "sigmaCon_",
      use_parameter_labels = use_parameter_labels,
      numerical_safeguard = numerical_safeguard
    )
  )
  parUncon <- parUncon[-(1:sigma_size)]

  structure(parCon, class = c("parCon", "numeric"))
}

#' @rdname heston_parameter_transformations
#' @export
parCon2par_heston <- function(parCon, controls, use_parameter_labels = TRUE) {

  stopifnot(inherits(parCon, "parCon"))
  stopifnot(inherits(controls, "Heston_controls"))
  # print("p3")
  M <- controls[["states"]][1]
  
  Gamma_size <- (M - 1) * M
  Gamma <- gammasCon2Gamma(
    parCon[1:Gamma_size], 
    M, 
    use_parameter_labels = use_parameter_labels
  )
  parCon <- parCon[-(1:Gamma_size)]
  
  kappa_size <- M
  kappa <- parCon[1:kappa_size]
  parCon <- parCon[-(1:kappa_size)]
  
  theta_size <- M
  theta <- parCon[1:theta_size]
  parCon <- parCon[-(1:theta_size)]
  
  sigma_size <- M
  sigma <- parCon[1:sigma_size]
  parCon <- parCon[-(1:sigma_size)]
  
  
  #print("return param")
  Heston_parameters(
    controls = controls,
    Gamma = Gamma, kappa = kappa, theta = theta, sigma = sigma,
    check_controls = FALSE
  )
}


#' @rdname heston_parameter_transformations
#' @export
par2parCon_heston <- function(par, controls, use_parameter_labels = TRUE) {
  stopifnot(inherits(par, "Heston_parameters"))
  stopifnot(inherits(controls, "Heston_controls"))
  parUncon2parCon_heston(
    par2parUncon_heston(par, controls, use_parameter_labels = use_parameter_labels), 
    controls,
    use_parameter_labels = use_parameter_labels
  )
}

#' @rdname heston_parameter_transformations
#' @export
parCon2parUncon_heston <- function(parCon, controls, use_parameter_labels = TRUE) { 
  stopifnot(inherits(parCon, "parCon"))
  stopifnot(inherits(controls, "Heston_controls"))
  par2parUncon_heston(
    parCon2par_heston(parCon, controls, use_parameter_labels = use_parameter_labels), 
    controls,
    use_parameter_labels = use_parameter_labels
  )
}

#' @rdname heston_parameter_transformations
#' @export
parUncon2par_heston <- function(
    parUncon, controls, use_parameter_labels = TRUE, numerical_safeguard = FALSE
) {
  
  stopifnot(inherits(parUncon, "parUncon"))
  stopifnot(inherits(controls, "Heston_controls"))
  
  parCon2par_heston(
    parUncon2parCon_heston(
      parUncon, controls, use_parameter_labels = use_parameter_labels,
      numerical_safeguard = numerical_safeguard
    ), 
    controls,
    use_parameter_labels = use_parameter_labels
  )
  
}








#' @rdname parameter_transformations
#' @return
#' For \code{Gamma2gammasCon}: a vector of constrained non-diagonal matrix 
#' elements (column-wise).

Gamma2gammasCon <- function(
    Gamma, prefix = "gammasCon_", use_parameter_labels = TRUE,
    numerical_safeguard = FALSE
) {
  gammasCon <- Gamma[row(Gamma) != col(Gamma)]
  if (isTRUE(numerical_safeguard)) {
    gammasCon <- replace(gammasCon, gammasCon == 0, 1e-3)
    gammasCon <- replace(gammasCon, gammasCon == 1, 1 - 1e-3)
  }
  if (isTRUE(use_parameter_labels)) {
    names(gammasCon) <- oeli::matrix_indices(
      Gamma, prefix = prefix, exclude_diagonal = TRUE
    )
  }
  return(gammasCon)
}

#' @rdname parameter_transformations
#' @return
#' For \code{Gamma2gammasUncon}: a vector of unconstrained non-diagonal matrix 
#' elements (column-wise).

Gamma2gammasUncon <- function(
    Gamma, prefix = "gammasUncon_", use_parameter_labels = TRUE
) {
  diag(Gamma) <- 0
  Gamma <- log(Gamma / (1 - rowSums(Gamma)))
  diag(Gamma) <- NA_real_
  gammasUncon <- Gamma[!is.na(Gamma)]
  if (isTRUE(use_parameter_labels)) {
    names(gammasUncon) <- oeli::matrix_indices(
      Gamma, prefix = prefix, exclude_diagonal = TRUE
    )
  }
  return(gammasUncon)
}

#' @rdname parameter_transformations
#' @return
#' For \code{gammasCon2Gamma}: a transition probability matrix.

gammasCon2Gamma <- function(
    gammasCon, dim, prefix = "state_", use_parameter_labels = TRUE
) {
  Gamma <- diag(dim)
  Gamma[!Gamma] <- gammasCon
  for (i in 1:dim) {
    Gamma[i, i] <- 1 - (rowSums(Gamma)[i] - 1)
  }
  if (isTRUE(use_parameter_labels)) {
    colnames(Gamma) <- rownames(Gamma) <- paste0(prefix, seq_len(dim))
  }
  return(Gamma)
}

#' @rdname parameter_transformations
#' @return
#' For \code{gammasCon2gammasUncon}: a vector of unconstrained non-diagonal 
#' elements of the transition probability matrix.

gammasCon2gammasUncon <- function(
    gammasCon, dim, prefix = "gammasUncon_", use_parameter_labels = TRUE
) {
  Gamma2gammasUncon(
    gammasCon2Gamma(gammasCon, dim, use_parameter_labels = use_parameter_labels), 
    prefix = prefix,
    use_parameter_labels = use_parameter_labels
  )
}

#' @rdname parameter_transformations
#' @return
#' For \code{gammasUncon2Gamma}: a transition probability matrix.

gammasUncon2Gamma <- function(
    gammasUncon, dim, prefix = "state_", use_parameter_labels = TRUE,
    numerical_safeguard = FALSE
) {
  Gamma <- diag(dim)
  Gamma[!Gamma] <- exp(gammasUncon)
  if (isTRUE(numerical_safeguard)) {
    Gamma[!is.finite(Gamma)] <- 100
  }
  Gamma <- Gamma / rowSums(Gamma)
  if (isTRUE(use_parameter_labels)) {
    colnames(Gamma) <- rownames(Gamma) <- paste0(prefix, seq_len(dim))
  }
  return(Gamma)
}

#' @rdname parameter_transformations
#' @return
#' For \code{gammasUncon2gammasCon}: a vector of constrained non-diagonal 
#' elements of a transition probability matrix.

gammasUncon2gammasCon <- function(
    gammasUncon, dim, prefix = "gammasCon_", use_parameter_labels = TRUE,
    numerical_safeguard = FALSE
) {
  Gamma2gammasCon(
    gammasUncon2Gamma(
      gammasUncon, dim, use_parameter_labels = use_parameter_labels,
      numerical_safeguard = numerical_safeguard
    ), 
    prefix = prefix,
    use_parameter_labels = use_parameter_labels,
    numerical_safeguard = numerical_safeguard
  )
}