
source("Heston_Standard/Heston_get_init.R")
source("Heston_Standard/fHMM_model.R")
library(iterators)
library(parallel)
#' Model fitting
#'
#' @description
#' This function fits a hidden Markov model via numerical likelihood 
#' maximization.
#'
#' @details
#' Multiple optimization runs starting from different initial values are 
#' computed in parallel if \code{ncluster > 1}.
#'
#' @param data
#' An object of class \code{\link{fHMM_data}}.
#' 
#' @param ncluster
#' Set the number of clusters for parallel optimization runs to reduce 
#' optimization time.
#' By default, \code{ncluster = 1} (no clustering).
#' 
#' @param verbose
#' Set to \code{TRUE} to print progress messages.
#' 
#' @inheritParams Heston_set_controls
#' @inheritParams get_initial_values
#' @inheritParams Heston_set_controls
#'
#' @return
#' An object of class \code{\link{fHMM_model}}.
#' 
#' @examples
#' ### 2-state HMM with normal distributions
#' 
#' # set specifications
#' controls <- Heston_set_controls(
#'   states = 2, sdds = "normal", horizon = 100, runs = 10
#' )
#' 
#' # define parameters
#' parameters <- fHMM_parameters(controls, mu = c(-1, 1), seed = 1)
#' 
#' # sample data
#' data <- prepare_data(controls, true_parameter = parameters, seed = 1)
#' 
#' # fit model
#' model <- fit_model(data, seed = 1)
#' 
#' # inspect fit
#' summary(model)
#' plot(model, "sdds")
#' 
#' # decode states
#' model <- decode_states(model)
#' plot(model, "ts")
#' 
#' # predict
#' predict(model, ahead = 5)
#'
#' @export

Heston_fit_model <- function (data, controls = data[["controls"]], fit = list(), 
                       runs = 10, origin = FALSE, accept = 1:3, gradtol = 0.01, 
                       iterlim = 1000, print.level = 0, steptol = 0.01, ncluster = 1, 
                       seed = NULL, verbose = TRUE, initial_estimate = NULL) 
{ 
  set.seed(999)
  if (!inherits(data, "Heston_data")) {
    stop("'data' is not of class 'Heston_data'.", call. = FALSE)
  }
  if (!checkmate::test_count(ncluster, positive = TRUE)) {
    stop("'ncluster' must be a positive integer.", call. = FALSE)
  }
  if (!checkmate::test_flag(verbose)) {
    stop("'verbose' must be either TRUE or FALSE.", call. = FALSE)
  }
  
  data[["controls"]] <- Heston_set_controls(controls = controls, hierarchy = data[["controls"]][["hierarchy"]], 
                                     states = data[["controls"]][["states"]], sdds = data[["controls"]][["sdds"]], 
                                     horizon = data[["controls"]][["horizon"]], period = data[["controls"]][["period"]], 
                                     data = data[["controls"]][["data"]], fit = fit, runs = runs, 
                                     origin = origin, accept = accept, gradtol = gradtol, 
                                     iterlim = iterlim, print.level = print.level, steptol = steptol)
  
  initial_values <- Heston_get_init(data = data, ncluster = ncluster, # problem
                                       seed = seed, verbose = verbose, initial_estimate = initial_estimate)
  
  # check nstate match the true regime
  runs <- length(initial_values)

  target <- Heston_nLL_hmm

  if (verbose) {
    pb <- progress::progress_bar$new(format = "[:bar] :percent, :eta ETA", 
                                     total = runs, width = 45, clear = TRUE, complete = "=", 
                                     incomplete = "-", current = ">")
    pb$message("Maximizing likelihood...")
  }

  
  

  start_time <- Sys.time()
  if (ncluster == 1) {
    mods <- list()
    for (run in seq_len(runs)) {
      if (verbose) 
        pb$tick(0)
      suppressWarnings({
        
        parameter_history <- list(
          par_uncon = list(),
          nll_value = numeric()
        )
        
        mod <- try(
          stats::nlm(f = target, p = initial_values[[run]], 
                              observations = data[["data"]], controls = data[["controls"]], 
                              iterlim = data[["controls"]][["fit"]][["iterlim"]],  # maximum number of iterations
                              steptol = data[["controls"]][["fit"]][["steptol"]], 
                              gradtol = data[["controls"]][["fit"]][["gradtol"]], 
                              print.level = data[["controls"]][["fit"]][["print.level"]], 
                              hessian = FALSE), silent = FALSE) # silent off
        
      })
      
      accept_run <- !inherits(mod, "try-error") && mod[["code"]] %in% 
        data[["controls"]][["fit"]][["accept"]]
     
      if (accept_run) {
        mods[[run]] <- mod
      }
      else {
        mods[[run]] <- NA
      }
      if (verbose) 
        pb$tick()
    }
  }

  end_time <- Sys.time()
  lls <- -unlist(sapply(mods, `[`, "minimum"), use.names = FALSE)
  if (all(is.na(lls))) {
    stop("None of the estimation runs ended successfully.\n", 
         "Adapt 'accept' or increase 'runs' in 'controls'.", 
         call. = FALSE)
  }
  if (verbose) 
    message("Approximating Hessian...")
  fisher <- pracma::hessdiag(f = target, x = mods[[which.max(lls)]][["estimate"]], 
                             observations = data[["data"]], controls = data[["controls"]])
  if (all(fisher > 0)) {
    inverse_fisher <- 1/fisher
  }
  else {
    hessian <- suppressWarnings(pracma::hessian(f = target, 
                                                x0 = mods[[which.max(lls)]][["estimate"]], observations = data[["data"]], 
                                                controls = data[["controls"]]))
    inverse_fisher <- diag(MASS::ginv(hessian))
  }
  if (verbose) 
    message("Fitting completed!")
  mod <- mods[[which.max(lls)]]
  ll <- -mod[["minimum"]]
  estimate <- mod[["estimate"]]
  class(estimate) <- "parUncon"
  estimation_time <- ceiling(difftime(end_time, start_time, 
                                      units = "mins"))
  out <- fHMM_model(data = data, estimate = estimate, nlm_output = mod, 
                    estimation_time = estimation_time, ll = ll, lls = lls, 
                    gradient = mod$gradient, inverse_fisher = inverse_fisher, 
                    decoding = NULL)

  out <- Heston_reorder_states(out, state_order = "mean")
  return(out)
}
