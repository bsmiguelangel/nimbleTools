#' Run NIMBLE models in parallel
#'
#' Run NIMBLE MCMC chains in parallel and return posterior samples, with
#' summaries and WAIC when requested.
#'
#' @param code NIMBLE model code, usually created with `nimble::nimbleCode()`.
#' @param data Named list with the observed data used by the model.
#' @param constants Named list with model constants.
#' @param inits Initial values for the MCMC chains.
#' @param nchains Number of MCMC chains to run. Defaults to 3.
#' @param seeds Numeric or integer vector with one seed per chain. If `NULL`,
#'   seeds are set to `1:nchains`.
#' @param summary If `TRUE`, posterior summaries are calculated when
#'   multiple chains are available.
#' @param monitors Character vector with the names of the model variables to
#'   return in the posterior samples. If `NULL`, the default monitored variables
#'   from the NIMBLE MCMC configuration are returned.
#' @param ntfyAccount Optional account or topic used by `notify()` to send a
#'   notification when the model run finishes.
#' @param email If `TRUE`, an email notification is sent through `notify()`.
#' @param WAIC If `TRUE`, the WAIC is calculated from the posterior samples.
#'   The stochastic parent nodes of the data nodes are automatically added to
#'   the internal monitored variables required for WAIC calculation.
#' @param HMC If `TRUE`, the model is configured to use HMC sampling
#'   inside `runParallel()`.
#' @param replaceSamplers Optional object passed to `runParallel()` to replace
#'   default NIMBLE samplers.
#' @param parallel If `TRUE`, chains are run in parallel using a cluster with
#'   one worker per chain. If `FALSE`, chains are run sequentially.
#' @param control.model Optional list of arguments passed to `nimble::nimbleModel()`.
#' @param control.compile Optional list of arguments passed to `nimble::compileNimble()`.
#' @param control.configure Optional list of arguments passed to
#'   `nimble::configureMCMC()` or `nimbleHMC::configureHMC()`.
#' @param control.build Optional list of arguments passed to
#'   `nimble::buildMCMC()`.
#' @param ... Additional arguments passed to `nimble::runMCMC()` through
#'   `runParallel()`.
#'
#' @returns A list with at least the following elements:
#' \describe{
#'   \item{samples}{Posterior samples as a `coda::mcmc.list` object.}
#'   \item{summary}{Posterior summary returned by `MCMCvis::MCMCsummary()`, if
#'   `summary = TRUE` and multiple chains are available.}
#'   \item{WAIC}{WAIC value, if `WAIC = TRUE`.}
#'   \item{total.time}{Total computation time in seconds.}
#' }
#'
#' @examples
#' \dontrun{
#' code <- nimble::nimbleCode({
#'   for (i in 1:N) {
#'     y[i] ~ dnorm(mu, sd = sigma)
#'   }
#'   mu ~ dnorm(0, sd = 10)
#'   sigma ~ dunif(0, 10)
#' })
#'
#' data <- list(y = rnorm(50))
#' constants <- list(N = 50)
#' inits <- list(mu = 0, sigma = 1)
#'
#' fit <- pNimble(code = code,
#'                data = data,
#'                constants = constants,
#'                inits = inits,
#'                nchains = 3,
#'                monitors = c("mu", "sigma"),
#'                parallel = TRUE)
#'
#' fit$samples
#' fit$summary
#' }
#'
#' @export
pNimble <- function(code = NULL, data = NULL, constants = NULL, inits = NULL,
                    nchains = 3, seeds = NULL, summary = TRUE, monitors = NULL,
                    ntfyAccount = NULL, email = FALSE, WAIC = FALSE, HMC = FALSE,
                    replaceSamplers = NULL, parallel = TRUE,
                    control.model = NULL, control.compile = NULL,
                    control.configure = NULL, control.build = NULL, ...) {

  # Initial arrangements and checkings
  time.start <- Sys.time()

  # Send notification when the function exits, even if an error occurs
  on.exit({
    time.end <- Sys.time()
    total.time <- as.numeric(difftime(time.end, time.start,
                                      units = "secs"))
    notify(time = total.time,
           email = email,
           ntfyAccount = ntfyAccount,
           model = "NIMBLE model")
  }, add = TRUE)

  # Create the control lists when they are not provided by the user
  if (is.null(control.model)) control.model <- list()

  # Add model components to the control.model list when provided directly
  if (!is.null(code)) control.model$code <- code
  if (!is.null(data)) control.model$data <- data
  if (!is.null(constants)) control.model$constants <- constants

  # Add monitored variables to the MCMC configuration
  if (!is.null(monitors)) {
    if (is.null(control.configure)) {
      control.configure <- list()
    }
    control.configure$monitors <- monitors
  }

  # Create the remaining control lists when they are not provided by the user
  if (is.null(control.configure)) control.configure <- list()
  if (is.null(control.compile)) control.compile <- list()
  if (is.null(control.build)) control.build <- list()

  # Use one seed per chain if seeds are not specified
  if (is.null(seeds)) seeds <- 1:nchains

  # Check that all input arguments are valid before running NIMBLE
  checkArguments(inits = inits, nchains = nchains, seeds = seeds, summary = summary,
                 monitors = monitors, ntfyAccount = ntfyAccount, email = email,
                 WAIC = WAIC, HMC = HMC, replaceSamplers = replaceSamplers, parallel = parallel,
                 control.model = control.model, control.compile = control.compile,
                 control.configure = control.configure, control.build = control.build)

  # Check that nimbleHMC is available when HMC sampling is requested
  if (HMC && !requireNamespace("nimbleHMC", quietly = TRUE)) {
    stop("Package 'nimbleHMC' is required when HMC = TRUE.")
  }

  # Parallelized call to Nimble
  if (parallel) {

    # A cluster with one worker per chain is created
    my.cluster <- parallel::makeCluster(nchains)

    # Ensure that the cluster is stopped when the function exits
    on.exit(parallel::stopCluster(my.cluster), add = TRUE)

    # Load NIMBLE in each worker
    parallel::clusterEvalQ(my.cluster, {
      library(nimble)
    })

    # Export the Leroux distribution functions to each worker so that NIMBLE
    # can find them when processing model code in parallel
    parallel::clusterExport(my.cluster,
                            varlist = c("dcar_leroux", "rcar_leroux"),
                            envir = environment(dcar_leroux))

    # Load nimbleHMC in each worker when HMC sampling is requested
    if (HMC) {
      parallel::clusterEvalQ(my.cluster, {
        library(nimbleHMC)
      })
    }

    # Run one independent NIMBLE chain for each seed in parallel
    resul <- parallel::parLapply(cl = my.cluster, X = seeds, fun = runParallel, inits = inits,
                                 control.model = control.model, control.compile = control.compile,
                                 control.configure = control.configure, control.build = control.build,
                                 HMC = HMC, replaceSamplers = replaceSamplers, WAIC = WAIC, ...)

  } else {

    # Run one independent NIMBLE chain for each seed sequentially
    resul <- lapply(X = seeds, FUN = runParallel, inits = inits,
                    control.model = control.model, control.compile = control.compile,
                    control.configure = control.configure, control.build = control.build,
                    HMC = HMC, replaceSamplers = replaceSamplers, WAIC = WAIC, ...)
  }

  # Output arrangement

  # Name the outputs according to the chain number
  names(resul) <- paste0("chain", 1:length(resul))
  resul2 <- list()

  # Convert the raw samples into a coda mcmc.list object
  samples.all <- coda::as.mcmc.list(lapply(resul, function(x) {coda::as.mcmc(x)}))

  # Keep only the variables requested by the user in the final posterior samples
  resul2$samples <- samples.all

  if (!is.null(monitors)) {

    sample.names <- colnames(samples.all[[1]])

    selected <- unlist(lapply(monitors, function(monitor) {
      sample.names[sample.names == monitor |
                     startsWith(sample.names, paste0(monitor, "["))]
    }), use.names = FALSE)

    selected <- unique(selected)

    if (length(selected) == 0) {
      stop("None of the variables in monitors were found in the posterior samples.")
    }

    resul2$samples <- coda::as.mcmc.list(lapply(samples.all, function(x) {
      coda::as.mcmc(x[, selected, drop = FALSE])
    }))
  }

  if (summary) {

    # Try to calculate posterior summaries, but return posterior samples if summary fails
    summary.result <- tryCatch({

      # Summary is only calculated when several chains are available
      if (length(resul2$samples) > 1) {

        # Combine all chains to detect problematic variables
        AllSamples <- do.call(rbind, resul2$samples)

        # Identify variables with zero variance or undefined variance
        Var0s <- which(apply(AllSamples, 2, stats::var) == 0)
        NAs <- which(is.na(apply(AllSamples, 2, stats::var)))

        if (length(Var0s) > 0 | length(NAs) > 0) {

          # Remove problematic variables before calculating the standard summary
          aux <- lapply(resul2$samples, function(x) {x[, -c(Var0s, NAs)]})
          summary.out <- MCMCvis::MCMCsummary(aux)

          if (length(Var0s) > 0) {

            # Add summary rows manually for variables with zero variance
            newSummaries <- cbind(apply(AllSamples[, Var0s], 2, mean),
                                  rep(0, length(Var0s)), rep(0, length(Var0s)),
                                  rep(0, length(Var0s)), rep(0, length(Var0s)),
                                  rep(NA, length(Var0s)), rep(NA, length(Var0s)))
            rownames(newSummaries) <- colnames(AllSamples)[Var0s]
            colnames(newSummaries) <- colnames(summary.out)
            summary.out <- rbind(summary.out, newSummaries)
          }

          if (length(NAs) > 0) {

            # Add empty summary rows for variables with undefined variance
            newSummaries <- matrix(rep(NA, 7 * length(NAs)), ncol = 7)
            rownames(newSummaries) <- colnames(AllSamples)[NAs]
            colnames(newSummaries) <- colnames(summary.out)
            summary.out <- rbind(summary.out, newSummaries)
          }

          # Reorder the summary to match the original variable order
          summary.out <- summary.out[match(colnames(AllSamples),
                                           rownames(summary.out)), ]
        } else {

          # Standard summary when no problematic variables are found
          summary.out <- MCMCvis::MCMCsummary(resul2$samples)
        }

        summary.out

      } else {
        warning("Summary cannot be calculated with a single chain. Posterior samples are returned without summary.")
        NULL
      }

    }, error = function(e) {

      warning("Summary could not be calculated. Posterior samples are returned without summary. ",
              "Original error: ", conditionMessage(e))

      NULL
    })

    if (!is.null(summary.result)) {
      resul2$summary <- summary.result
    }
  }

  if (WAIC) {

    # Try to calculate WAIC, but return posterior samples if WAIC fails
    WAIC.result <- tryCatch({

      # Register the Leroux CAR distribution before rebuilding the model for WAIC
      load_leroux()

      # Rebuild and compile the model to calculate WAIC from the posterior samples
      model.nimble <- nimble::nimbleModel(code = code, data = data, constants = constants,
                                          check = FALSE, calculate = FALSE, buildDerivs = FALSE)
      cmodel <- nimble::compileNimble(model.nimble)

      nimble::calculateWAIC(do.call(rbind, samples.all), cmodel)

    }, error = function(e) {

      warning("WAIC could not be calculated. Posterior samples are returned without WAIC. ",
              "Original error: ", conditionMessage(e))

      NA
    })

    resul2$WAIC <- WAIC.result
  }

  # Add total computation time to the output object
  time.end <- Sys.time()
  resul2$total.time <- as.numeric(difftime(time.end, time.start,
                                           units = "secs"))

  return(resul2)
}
