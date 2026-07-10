#' Identify problematic MCMC parameters
#'
#' Identify parameters with problematic MCMC behaviour and optionally produce
#' traceplots only for those parameters.
#'
#' The function can be used with posterior samples from `pNimble()` or from
#' other MCMC models. It calculates posterior summaries using
#' `MCMCvis::MCMCsummary()` and identifies parameters with `Rhat`
#' greater than `Rhat.max` or effective sample size `n.eff` lower than
#' `n.eff.min`. If requested, traceplots are produced only for the problematic
#' parameters using `MCMCvis::MCMCtrace()`.
#'
#' @param object Object returned by `pNimble()` containing a `samples` element,
#'   or posterior samples as a `coda::mcmc.list` object.
#' @param params Optional character vector with the names of the parameters to
#'   check. If `NULL`, all monitored parameters are checked.
#' @param Rhat.max Maximum acceptable `Rhat` value. Parameters with `Rhat`
#'   greater than this value are identified as problematic. The default is
#'   `1.1`.
#' @param n.eff.min Minimum acceptable effective sample size. Parameters with
#'   `n.eff` lower than this value are identified as problematic. The default
#'   is `100`.
#' @param plot Logical value. If `TRUE`, traceplots are produced for the
#'   problematic parameters. The default is `TRUE`.
#' @param pdf Logical value passed to `MCMCvis::MCMCtrace()`. If `TRUE`, plots
#'   are saved to a PDF file. The default is `FALSE`.
#' @param ind Logical value passed to `MCMCvis::MCMCtrace()`. If `TRUE`,
#'   density lines are shown separately for each chain. The default is `TRUE`.
#' @param round Number of decimal places used by `MCMCvis::MCMCsummary()`.
#'   The default is `4`.
#' @param ... Additional arguments passed to `MCMCvis::MCMCsummary()` and
#'   `MCMCvis::MCMCtrace()`.
#'
#' @returns A posterior summary table restricted to the parameters with
#'   problematic MCMC behaviour. If no problematic parameters are found, an
#'   empty summary table is returned.
#'
#' @examples
#' \dontrun{
#' MCMCproblems(object = fit,
#'              params = c("rho", "theta", "beta_age"),
#'              Rhat.max = 1.02, n.eff.min = 400)
#' }
#'
#' @export
MCMCproblems <- function(object, params = NULL,
                         Rhat.max = 1.1, n.eff.min = 100,
                         plot = TRUE, pdf = FALSE, ind = TRUE,
                         round = 4, ...) {

  # Extract posterior samples if a pNimble output object is provided
  if (is.list(object) && !is.null(object$samples)) {
    samples <- object$samples
  } else {
    samples <- object
  }

  # Calculate posterior summaries
  summary.out <- MCMCvis::MCMCsummary(object = samples, params = params,
                                      round = round, exact = TRUE,
                                      ISB = FALSE, ...)

  # Check that required diagnostic columns are available
  if (!all(c("Rhat", "n.eff") %in% colnames(summary.out))) {
    stop("MCMCsummary output must contain columns named 'Rhat' and 'n.eff'.")
  }

  # Identify parameters with problematic MCMC behaviour
  problematic <- summary.out[
    summary.out[, "Rhat"] > Rhat.max | summary.out[, "n.eff"] < n.eff.min,
    ,
    drop = FALSE
  ]

  # Return an empty table if no problematic parameters are found
  if (nrow(problematic) == 0) {
    message("No parameters with Rhat > ", Rhat.max,
            " or n.eff < ", n.eff.min, " were found.")
    return(problematic)
  }

  # Plot only problematic parameters when requested
  if (plot) {
    MCMCvis::MCMCtrace(object = samples,
                       params = rownames(problematic),
                       pdf = pdf, ind = ind,
                       exact = TRUE, ISB = FALSE,
                       Rhat = TRUE, n.eff = TRUE, ...)
  }

  return(problematic)
}
