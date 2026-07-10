#' Random generation function for the Leroux CAR distribution
#'
#' Random generation function associated with [dcar_leroux()].
#'
#' This function is included for compatibility with NIMBLE user-defined
#' distributions. Random generation from the Leroux CAR distribution is not
#' implemented, and the function stops if called.
#'
#' @param n Number of observations to simulate. This argument is included for
#'   compatibility with NIMBLE user-defined distributions.
#' @param rho Spatial dependence parameter. Values close to 0 correspond to
#'   weak spatial dependence, while values close to 1 correspond to strong
#'   spatial dependence.
#' @param sd Marginal standard deviation parameter of the spatial random
#'   effects.
#' @param Lambda Numeric vector containing the eigenvalues of
#'   \eqn{\boldsymbol{D} - \boldsymbol{W}}, where \eqn{\boldsymbol{D}} is the
#'   diagonal matrix of the numbers of neighbours and \eqn{\boldsymbol{W}} is
#'   the neighbourhood matrix. This object can be constructed using
#'   [lerouxObjects()].
#' @param from.to Matrix with two columns defining the distinct neighbouring
#'   pairs. Each row contains the indices of two neighbouring spatial units.
#'   The implementation assumes that each neighbouring pair is included once.
#'   This object can be constructed using [lerouxObjects()].
#' @param zero_mean Numeric indicator. If `zero_mean = 1`, a zero-mean
#'   constraint is added to the spatial random effects. If `zero_mean = 0`,
#'   no zero-mean constraint is added. The default is `0`.
#'
#' @returns This function does not return simulated values. It stops with an
#'   error because random generation is not implemented.
#'
#' @seealso [dcar_leroux()], [lerouxObjects()]
#'
#' @export
rcar_leroux <- nimble::nimbleFunction(
  name = "rcar_leroux",
  run = function(n = integer(0),
                 rho = double(0),
                 sd = double(0),
                 Lambda = double(1),
                 from.to = double(2),
                 zero_mean = double(0, default = 0)) {

    returnType(double(1))
    nimStop("user-defined distribution dcar_leroux provided without random generation function.")
    x <- nimNumeric(length(Lambda))
    return(x)
  }
)
