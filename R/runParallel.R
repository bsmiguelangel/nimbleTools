#' Run one NIMBLE chain
#'
#' Build, compile, configure and run one NIMBLE MCMC chain.
#'
#' @param seed Numeric value used as seed for the chain.
#' @param inits Optional function returning initial values for the chain.
#' @param control.model Optional list of arguments passed to
#'   `nimble::nimbleModel()`.
#' @param control.compile Optional list of arguments passed to
#'   `nimble::compileNimble()`.
#' @param control.configure Optional list of arguments passed to
#'   `nimble::configureMCMC()` or `nimbleHMC::configureHMC()`.
#' @param control.build Optional list of arguments passed to
#'   `nimble::buildMCMC()`.
#' @param HMC Logical value indicating whether HMC sampling should be used.
#' @param WAIC Logical value indicating whether the stochastic parent nodes of
#'   the data nodes should be added to the internal monitored variables for
#'   WAIC calculation.
#' @param replaceSamplers Optional list used to replace default NIMBLE samplers.
#' @param ... Additional arguments passed to `nimble::runMCMC()`.
#'
#' @returns Posterior samples returned by `nimble::runMCMC()`.
#'
#' @keywords internal
runParallel <- function(seed, inits = NULL, control.model, control.compile,
                        control.configure, control.build,
                        HMC, WAIC, replaceSamplers, ...) {

  # Set the seed for reproducibility of this chain
  set.seed(seed)

  # Check that nimbleHMC is available when HMC sampling is requested
  if (HMC && !requireNamespace("nimbleHMC", quietly = TRUE)) {
    stop("Package 'nimbleHMC' is required when HMC = TRUE.")
  }

  # Register the Leroux CAR distribution before building the model
  load_leroux()

  # Create empty control lists when they are not provided
  if (is.null(control.configure)) control.configure <- list()
  if (is.null(control.compile)) control.compile <- list()
  if (is.null(control.build)) control.build <- list()

  # Generate initial values for this chain, if provided
  if (!is.null(inits)) control.model$inits <- inits()

  # Build the NIMBLE model
  model.nimble <- do.call(nimble::nimbleModel,
                          c(control.model, buildDerivs = HMC),
                          quote = TRUE)

  # Compile the model before configuring the sampler
  model.precompiled <- do.call(nimble::compileNimble,
                               c(model.nimble, control.compile),
                               quote = TRUE)

  # Add the stochastic parent nodes of the data nodes to the internal monitors
  # when WAIC is requested
  if (WAIC) {
    control.configure$monitors <- unique(c(control.configure$monitors,
                                           model.nimble$getParents(
                                             model.nimble$getNodeNames(dataOnly = TRUE),
                                             stochOnly = TRUE)))
  }

  # Configure the sampler, using HMC or standard MCMC
  if (HMC) {
    model.configure <- do.call(nimbleHMC::configureHMC, c(model.nimble, control.configure), quote = TRUE)
  } else {
    model.configure <- do.call(nimble::configureMCMC, c(model.nimble, control.configure), quote = TRUE)
  }

  # Replace default samplers when requested
  if (!is.null(replaceSamplers)) {
    for (i in seq_along(replaceSamplers[[1]])) {
      theith <- lapply(replaceSamplers, function(y, j = i) {y[[j]]})
      theith <- Filter(Negate(anyNA), theith)
      do.call(model.configure$replaceSamplers, theith)
    }
    model.configure$printSamplers()
  }

  # Build and compile the configured MCMC
  model.build <- do.call(nimble::buildMCMC, c(model.configure, control.build), quote = TRUE)

  model.compiled <- do.call(nimble::compileNimble, c(model.build, control.compile), quote = TRUE)

  # Run a single chain and return posterior samples
  model.output <- nimble::runMCMC(mcmc = model.compiled, nchains = 1, setSeed = seed,
                                  progressBar = FALSE, samplesAsCodaMCMC = TRUE,
                                  summary = FALSE, WAIC = FALSE, ...)

  return(model.output)
}
