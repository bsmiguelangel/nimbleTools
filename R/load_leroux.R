#' Load Leroux distribution
#'
#' Register the Leroux CAR distribution so that it can be used in NIMBLE model
#' code.
#'
#' @returns No return value.
#'
#' @keywords internal
load_leroux <- function() {

  # Deregister previous Leroux definition, if any
  suppressMessages(
    suppressWarnings(
      try(nimble::deregisterDistributions("dcar_leroux"), silent = TRUE)
    )
  )

  # Register the Leroux distribution in NIMBLE
  suppressMessages(
    nimble::registerDistributions(list(
      dcar_leroux = list(
        BUGSdist = "dcar_leroux(rho, sd, Lambda, from.to, zero_mean)",
        types = c(
          "value = double(1)",
          "rho = double(0)",
          "sd = double(0)",
          "Lambda = double(1)",
          "from.to = double(2)",
          "zero_mean = double(0)"
        ),
        pqAvail = FALSE
      )
    ))
  )
}
