#' Monitor MCMC diagnostics
#'
#' Monitor MCMC diagnostics and optionally produce traceplots for selected
#' parameters.
#'
#' The function can be used with posterior samples from `pNimble()` or from
#' other MCMC models. When available, the posterior summary stored in the
#' `pNimble()` output is used directly. Otherwise, posterior summaries are
#' calculated using `MCMCvis::MCMCsummary()`.
#'
#' Parameters with `Rhat` greater than `Rhat.max` or effective sample size
#' `n.eff` lower than `n.eff.min` are identified as problematic. Parameters with
#' no variability or missing diagnostic values, such as reference categories of
#' factors, are not treated as problematic.
#'
#' If `only.problematic = TRUE`, only parameters with problematic MCMC behaviour
#' are returned and plotted. If `only.problematic = FALSE`, all selected
#' parameters are returned and plotted, regardless of their `Rhat` or `n.eff`
#' values.
#'
#' @param object Object returned by `pNimble()` containing a `samples` element,
#'   or posterior samples as a `coda::mcmc.list` object.
#' @param params Optional character vector with the names of the parameters to
#'   check. If `NULL`, all monitored parameters are checked.
#' @param only.problematic Logical value. If `TRUE`, only parameters with
#'   problematic MCMC behaviour are returned and plotted. If `FALSE`, all
#'   selected parameters are returned and plotted. The default is `TRUE`.
#' @param Rhat.max Maximum acceptable `Rhat` value. Parameters with `Rhat`
#'   greater than this value are identified as problematic. The default is
#'   `1.10`.
#' @param n.eff.min Minimum acceptable effective sample size. Parameters with
#'   `n.eff` lower than this value are identified as problematic. The default
#'   is `100`.
#' @param plot Logical value. If `TRUE`, traceplots are produced. The default is
#'   `TRUE`.
#' @param round Number of decimal places used for the posterior summary returned
#'   by the function. When the posterior summary has to be calculated, this
#'   value is also passed to `MCMCvis::MCMCsummary()`. The default is `4`.
#' @param type Type of plot produced by `MCMCvis::MCMCtrace()`. The default is
#'   `"trace"`, which produces traceplots only.
#' @param pdf Logical value passed to `MCMCvis::MCMCtrace()`. If `TRUE`, plots
#'   are saved to a PDF file. The default is `FALSE`.
#' @param ind Logical value passed to `MCMCvis::MCMCtrace()`. If `TRUE`,
#'   density plots are shown separately for each chain. The default is `TRUE`.
#' @param exact Logical value used when selecting parameters through `params`.
#'   If `TRUE`, parameter names are matched exactly after applying the `ISB`
#'   rule. If `FALSE`, parameter names are matched using regular expressions.
#'   The default is `TRUE`.
#' @param ISB Logical value used when selecting parameters through `params`.
#'   If `TRUE`, square brackets and indices are ignored when matching parameter
#'   names. If `FALSE`, square brackets are kept. The default is `FALSE`.
#' @param Rhat Logical value passed to `MCMCvis::MCMCtrace()`. If `TRUE`, `Rhat`
#'   values are added to the plots. The default is `TRUE`.
#' @param n.eff Logical value passed to `MCMCvis::MCMCtrace()`. If `TRUE`,
#'   effective sample sizes are added to the plots. The default is `TRUE`.
#' @param ... Additional arguments passed to `MCMCvis::MCMCtrace()`.
#'
#' @returns A posterior summary table restricted either to the problematic
#'   parameters or to all selected parameters, depending on `only.problematic`.
#'   If no problematic parameters are found when `only.problematic = TRUE`, a
#'   message is printed and `NULL` is returned invisibly.
#'
#' @examples
#' \dontrun{
#' MCMCmonitor(object = fit,
#'             params = c("rho", "theta", "beta_age"),
#'             ISB = TRUE,
#'             Rhat.max = 1.02,
#'             n.eff.min = 400)
#'
#' MCMCmonitor(object = fit,
#'             params = "beta_age[2]",
#'             only.problematic = FALSE,
#'             ISB = FALSE)
#' }
#'
#' @export
MCMCmonitor <- function(object, params = NULL, only.problematic = TRUE,
                        Rhat.max = 1.10, n.eff.min = 100,
                        plot = TRUE, round = 4,
                        type = "trace", pdf = FALSE, ind = TRUE,
                        exact = TRUE, ISB = FALSE,
                        Rhat = TRUE, n.eff = TRUE, ...) {

  # Extract posterior samples and posterior summary if a pNimble output object
  # is provided
  if (is.list(object) && !is.null(object$samples)) {
    samples <- object$samples
    summary.out <- object$summary
  } else {
    samples <- object
    summary.out <- NULL
  }

  # Use the existing posterior summary when available. Otherwise, calculate it
  # from the posterior samples.
  if (is.null(summary.out)) {

    if (is.null(params)) {
      summary.out <- MCMCvis::MCMCsummary(object = samples,
                                          round = round)
    } else {
      summary.out <- MCMCvis::MCMCsummary(object = samples,
                                          params = params,
                                          exact = exact,
                                          ISB = ISB,
                                          round = round)
    }

  } else {

    # If parameters are selected and an existing summary is available, use
    # MCMCvis only to identify the selected MCMC columns. The posterior summary
    # itself is not recalculated.
    if (!is.null(params)) {

      selected.samples <- MCMCvis::MCMCchains(object = samples,
                                              params = params,
                                              exact = exact,
                                              ISB = ISB,
                                              mcmc.list = TRUE)

      selected <- colnames(selected.samples[[1]])
      selected <- selected[selected %in% rownames(summary.out)]

      if (length(selected) == 0) {
        stop("None of the parameters in params were found in the posterior summary.")
      }

      summary.out <- summary.out[selected, , drop = FALSE]
    }
  }

  # Check that required diagnostic columns are available
  if (!all(c("Rhat", "n.eff") %in% colnames(summary.out))) {
    stop("MCMCsummary output must contain columns named 'Rhat' and 'n.eff'.")
  }

  # Remove parameters with no variability, such as reference categories of
  # factors, before identifying problematic parameters. The original summary is
  # kept so that all selected parameters can still be returned when
  # only.problematic = FALSE.
  if ("sd" %in% colnames(summary.out)) {
    no.variability <- !is.na(summary.out[, "sd"]) &
      summary.out[, "sd"] == 0 &
      is.na(summary.out[, "Rhat"])
  } else {
    no.variability <- rep(FALSE, nrow(summary.out))
  }

  summary.check <- summary.out[!no.variability, , drop = FALSE]

  # Identify parameters with problematic MCMC behaviour
  is.problematic <- summary.check[, "Rhat"] > Rhat.max |
    summary.check[, "n.eff"] < n.eff.min

  # Parameters with missing diagnostic values are not treated as problematic
  is.problematic[is.na(is.problematic)] <- FALSE

  problematic <- summary.check[is.problematic, , drop = FALSE]

  # Round the summaries returned by the function
  summary.return <- summary.out
  problematic.return <- problematic

  if (!is.null(round)) {
    summary.return <- base::round(summary.return, digits = round)
    problematic.return <- base::round(problematic.return, digits = round)
  }

  # If only.problematic = FALSE, return and plot all selected parameters
  if (!only.problematic) {

    if (plot) {
      MCMCvis::MCMCtrace(object = samples,
                         params = rownames(summary.out),
                         type = type, pdf = pdf, ind = ind,
                         exact = TRUE, ISB = FALSE,
                         Rhat = Rhat, n.eff = n.eff, ...)
    }

    return(summary.return)
  }

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
                       exact = TRUE, ISB = FALSE,
                       Rhat = Rhat, n.eff = n.eff, ...)
  }

  return(problematic.return)
}
