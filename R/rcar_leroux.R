utils::globalVariables(c("rnorm", "nimNumeric", "returnType"))

#' Random generation function associated with the Leroux CAR distribution
#'
#' Simple random generation function associated with [dcar_leroux()].
#'
#' This function is included for compatibility with NIMBLE user-defined
#' distributions. It generates independent normal values with mean zero and
#' standard deviation `0.1`. Therefore, it should be understood as a simple
#' generator for initial or simulated values, not as an exact simulation from
#' the Leroux CAR distribution.
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
#' @param zero_mean Numeric indicator. If `zero_mean = 1`, generated values are
#'   centred to have mean zero. If `zero_mean = 0`, generated values are not
#'   centred. The default is `0`.
#'
#' @returns A numeric vector of simulated values.
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

    # returnType(double(1))
    # nimStop("user-defined distribution dcar_leroux provided without random generation function.")
    # x <- nimNumeric(length(Lambda))
    # return(x)

    # Number of small areas
    NAreas <- length(Lambda)

    # Simulated values
    x <- nimNumeric(NAreas)

    for (i in 1:NAreas) {
      x[i] <- rnorm(1, mean = 0, sd = 0.1)
    }

    # Centre generated values when requested
    if (zero_mean == 1) {
      x.mean <- mean(x[1:NAreas])
      for (i in 1:NAreas) {
        x[i] <- x[i] - x.mean
      }
    }

    returnType(double(1))
    return(x)
  }
)
