#' Construct Leroux CAR objects from a neighbourhood matrix
#'
#' Construct the objects required by `dcar_leroux()` from a binary
#' neighbourhood matrix.
#'
#' The function takes a binary adjacency matrix `W` and returns the matrix of
#' distinct neighbouring pairs `from.to`, the number of distinct neighbouring
#' pairs `NDist`, and the eigenvalues `Lambda` of
#' \eqn{\boldsymbol{D} - \boldsymbol{W}}, where \eqn{\boldsymbol{D}} is the
#' diagonal matrix of row sums of `W`.
#'
#' @param W Binary neighbourhood matrix. Rows and columns must correspond to
#'   the same spatial units. The matrix must be square, symmetric, contain only
#'   `0` and `1` values, and have a zero diagonal.
#'
#' @returns A list with the following elements:
#' \describe{
#'   \item{from.to}{Matrix with two columns defining the distinct neighbouring
#'   pairs. Each row contains the indices of two neighbouring spatial units.}
#'   \item{Lambda}{Numeric vector containing the eigenvalues of
#'   \eqn{\boldsymbol{D} - \boldsymbol{W}}.}
#'   \item{NDist}{Number of distinct neighbouring pairs.}
#' }
#'
#' @seealso [dcar_leroux()]
#'
#' @examples
#' W <- matrix(c(0, 1, 0,
#'               1, 0, 1,
#'               0, 1, 0),
#'             nrow = 3, byrow = TRUE)
#'
#' leroux.obj <- lerouxObjects(W)
#' leroux.obj$from.to
#' leroux.obj$Lambda
#' leroux.obj$NDist
#'
#' @export
lerouxObjects <- function(W) {

  # Check that W is a matrix
  if (!is.matrix(W)) {
    stop("W must be a matrix.")
  }

  # Check that W is square
  if (nrow(W) != ncol(W)) {
    stop("W must be a square matrix.")
  }

  # Check that W does not contain missing values
  if (anyNA(W)) {
    stop("W must not contain missing values.")
  }

  # Check that W is binary
  if (!all(W %in% c(0, 1))) {
    stop("W must contain only 0 and 1 values.")
  }

  # Check that W has a zero diagonal
  if (any(diag(W) != 0)) {
    stop("W must have a zero diagonal.")
  }

  # Check that W is symmetric
  if (!isTRUE(all.equal(W, t(W)))) {
    stop("W must be symmetric.")
  }

  # Number of spatial units
  N <- nrow(W)

  # Matrix D - W
  Q <- diag(rowSums(W), nrow = N) - W

  # Number of neighbours of each area
  nadj <- apply(W, 1, sum)

  # Neighbours of each area
  map <- unlist(apply(W, 1, function(x) which(x != 0)))

  # Sum of all the neighbour numbers of all areas
  nadj.tot <- length(map)

  # Cumulative sums of the number of neighbours of each area
  index <- c(0, cumsum(nadj))

  # All the neighbourhoods j ~ i where i < j
  from.to <- cbind(rep(1:N, times = nadj), map)
  colnames(from.to) <- c("from", "to")
  from.to <- from.to[which(from.to[, 1] < from.to[, 2]), ]

  # Number of distinct pairs of neighbours
  NDist <- nrow(from.to)

  # Eigenvalues of D - W
  Lambda <- eigen(Q)$values

  return(list(from.to = from.to, Lambda = Lambda, NDist = NDist))
}
