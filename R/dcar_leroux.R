utils::globalVariables(c("ADbreak", "dnorm", "nimInteger", "nimDim",
                         "nimNumeric", "pow", "returnType"))

#' Density function for the Leroux CAR distribution
#'
#' Density function for the Leroux conditional autoregressive distribution for
#' use in NIMBLE models.
#'
#' The Leroux CAR prior is a spatial prior for a vector of random effects
#' \eqn{\boldsymbol{\theta} = (\theta_1, \dots, \theta_K)}. Given a spatial
#' neighbourhood matrix \eqn{\boldsymbol{W}}, the model can be written through
#' the full conditional distributions
#'
#' \deqn{
#' \theta_k \mid \boldsymbol{\theta}_{-k}, \sigma^2, \rho
#' \sim
#' N\left(
#' \frac{\rho}{1 - \rho + \rho w_{k+}}
#' \sum_{j \sim k} \theta_j,
#' \frac{\sigma^2}{1 - \rho + \rho w_{k+}}
#' \right),
#' }
#'
#' where \eqn{w_{k+}} is the number of neighbours of spatial unit \eqn{k}.
#' The parameter \eqn{\rho} controls the amount of spatial dependence.
#'
#' Equivalently, the joint distribution can be written as
#'
#' \deqn{
#' \boldsymbol{\theta} \mid \sigma^2, \rho
#' \sim
#' N\left(
#' \boldsymbol{0},
#' \sigma^2
#' \left[
#' \rho(\boldsymbol{D} - \boldsymbol{W}) + (1 - \rho)\boldsymbol{I}
#' \right]^{-1}
#' \right),
#' }
#'
#' where \eqn{\boldsymbol{D}} is the diagonal matrix of row sums of
#' \eqn{\boldsymbol{W}}. The implementation uses the eigenvalues of
#' \eqn{\boldsymbol{D} - \boldsymbol{W}} and the list of neighbouring pairs to
#' compute the log-density efficiently.
#'
#' @param x Numeric vector of spatial random effects.
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
#' @param log Logical or integer indicator. If `TRUE`, the log-density is
#'   returned. If `FALSE`, the density is returned.
#'
#' @returns The density or log-density of the Leroux CAR distribution.
#'
#' @seealso [rcar_leroux()], [lerouxObjects()]
#'
#' @examples
#' \dontrun{
#' W <- matrix(c(0, 1, 0,
#'               1, 0, 1,
#'               0, 1, 0),
#'             nrow = 3, byrow = TRUE)
#'
#' leroux.obj <- lerouxObjects(W)
#' Lambda <- leroux.obj$Lambda
#' from.to <- leroux.obj$from.to
#' NDist <- leroux.obj$NDist
#' NAreas <- nrow(W)
#'
#' code <- nimble::nimbleCode({
#'   theta[1:NAreas] ~ dcar_leroux(rho = rho,
#'                            sd = 1,
#'                            Lambda = Lambda[1:NAreas],
#'                            from.to = from.to[1:NDist, 1:2],
#'                            zero_mean = 0)
#' })
#' }
#'
#' @export
dcar_leroux <- nimble::nimbleFunction(
  name = "dcar_leroux",
  run = function(x = double(1),                      # Spatial random effect (vector)
                 rho = double(0),                    # Amount of spatial dependence (scalar)
                 sd = double(0),                     # Standard deviation (scalar)
                 Lambda = double(1),                 # Eigenvalues of matrix D - W
                 from.to = double(2),                # Matrix of distinct pairs of neighbours from.to[, 1] < from.to[, 2]
                 zero_mean = double(0, default = 0), # Apply zero-mean constraint when set to 1
                 log = integer(0, default = 0)) {

    # Number of small areas
    NAreas <- dim(x)[1]

    # Number of distinct pairs of neighbours
    NDist <- dim(from.to)[1]

    # Required vectors
    from <- nimInteger(NDist)
    to <- nimInteger(NDist)
    x.from <- nimNumeric(NDist)
    x.to <- nimNumeric(NDist)

    for (Dist in 1:NDist) {
      from[Dist] <- ADbreak(from.to[Dist, 1])
      to[Dist] <- ADbreak(from.to[Dist, 2])
      x.from[Dist] <- x[from[Dist]]
      x.to[Dist] <- x[to[Dist]]
    }

    logDens <- sum(dnorm(x[1:NAreas], mean = 0, sd = sd * pow(1 - rho, -1/2), log = TRUE)) -
      NAreas/2 * log(1 - rho) +  sum(log(rho * (Lambda[1:NAreas] - 1) + 1))/2 -
      pow(sd, -2) * rho * sum(pow(x.from[1:NDist] - x.to[1:NDist], 2))/2

    # Add a zero-mean constraint when requested
    if (zero_mean == 1) {
      logDens <- logDens + dnorm(mean(x[1:NAreas]),
                                 mean = 0,
                                 sd = sd / (10 * pow(NAreas, 1/2)),
                                 log = TRUE)
    }

    returnType(double(0))
    if (log) {
      return(logDens)
    } else {
      return(exp(logDens))
    }

  },
  buildDerivs = list(run = list(ignore = c("Dist", "NAreas", "NDist")))
)
