#' Identify problematic MCMC parameters
#'
#' Identify parameters with problematic MCMC behaviour and optionally produce
#' traceplots only for those parameters.
#'
#' The function can be used with posterior samples from `pNimble()` or from
#' other MCMC models. It calculates posterior summaries using
#' `MCMCvis::MCMCsummary()` and identifies parameters with `Rhat`
#' greater than `Rhat.max` or effective sample size `n.eff` lower than
#' `n.eff.min`. Parameters with no variability, such as reference categories
#' of factors, are excluded when they have `NaN` in `Rhat` and `0` in `n.eff`.
#' If requested, traceplots are produced only for the problematic parameters
#' using `MCMCvis::MCMCtrace()`.
#'
#' @param object Object returned by `pNimble()` containing a `samples` element,
#'   or posterior samples as a `coda::mcmc.list` object.
#' @param params Optional character vector with the names of the parameters to
#'   check. If `NULL`, all monitored parameters are checked.
#' @param Rhat.max Maximum acceptable `Rhat` value. Parameters with `Rhat`
#'   greater than this value are identified as problematic. The default is
#'   `1.10`.
#' @param n.eff.min Minimum acceptable effective sample size. Parameters with
#'   `n.eff` lower than this value are identified as problematic. The default
#'   is `100`.
#' @param plot Logical value. If `TRUE`, traceplots are produced for the
#'   problematic parameters. The default is `TRUE`.
#' @param round Number of decimal places used by `MCMCvis::MCMCsummary()`.
#'   The default is `4`.
#' @param type Type of plot produced by `MCMCvis::MCMCtrace()`. The default is
#'   `"trace"`, which produces traceplots only.
#' @param pdf Logical value passed to `MCMCvis::MCMCtrace()`. If `TRUE`, plots
#'   are saved to a PDF file. The default is `FALSE`.
#' @param ind Logical value passed to `MCMCvis::MCMCtrace()`. The default is
#'   `TRUE`.
#' @param exact Logical value passed to `MCMCvis::MCMCtrace()`. The default is
#'   `TRUE`.
#' @param ISB Logical value passed to `MCMCvis::MCMCtrace()`. The default is
#'   `FALSE`.
#' @param Rhat Logical value passed to `MCMCvis::MCMCtrace()`. If `TRUE`, `Rhat`
#'   values are added to the plots. The default is `TRUE`.
#' @param n.eff Logical value passed to `MCMCvis::MCMCtrace()`. If `TRUE`,
#'   effective sample sizes are added to the plots. The default is `TRUE`.
#' @param ... Additional arguments passed to `MCMCvis::MCMCtrace()`.
#'
#' @returns A posterior summary table restricted to the parameters with
#'   problematic MCMC behaviour. If the result of `MCMCproblems()` is assigned
#'   to an object and problematic parameters are found, this summary table is
#'   stored in that object. If no problematic parameters are found, a message is
#'   printed and `NULL` is returned invisibly.
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
                         Rhat.max = 1.10, n.eff.min = 100,
                         plot = TRUE, round = 4,
                         type = "trace", pdf = FALSE, ind = TRUE,
                         exact = TRUE, ISB = FALSE,
                         Rhat = TRUE, n.eff = TRUE, ...) {

  # Extract posterior samples if a pNimble output object is provided
  if (is.list(object) && !is.null(object$samples)) {
    samples <- object$samples
  } else {
    samples <- object
  }

  # Calculate posterior summaries
  if (is.null(params)) {
    summary.out <- MCMCvis::MCMCsummary(object = samples,
                                        round = round)
  } else {
    summary.out <- MCMCvis::MCMCsummary(object = samples,
                                        params = params,
                                        round = round)
  }

  # Check that required diagnostic columns are available
  if (!all(c("Rhat", "n.eff") %in% colnames(summary.out))) {
    stop("MCMCsummary output must contain columns named 'Rhat' and 'n.eff'.")
  }

  # Remove parameters with no variability, such as reference categories of
  # factors, before identifying problematic parameters
  no.variability <- is.nan(summary.out[, "Rhat"]) &
    !is.na(summary.out[, "n.eff"]) &
    summary.out[, "n.eff"] == 0

  summary.out <- summary.out[!no.variability, , drop = FALSE]

  # Identify parameters with problematic MCMC behaviour
  problematic <- summary.out[
    summary.out[, "Rhat"] > Rhat.max | summary.out[, "n.eff"] < n.eff.min,
    ,
    drop = FALSE
  ]

  # Print a message if no problematic parameters are found
  if (nrow(problematic) == 0) {
    message("No parameters with Rhat > ", Rhat.max,
            " or n.eff < ", n.eff.min, " were found.")
    return(invisible(NULL))
  }

  # Plot only problematic parameters when requested
  if (plot) {
    MCMCvis::MCMCtrace(object = samples,
                       params = rownames(problematic),
                       type = type, pdf = pdf, ind = ind,
                       exact = exact, ISB = ISB, Rhat = Rhat,
                       n.eff = n.eff, ...)
  }

  return(problematic)
}
