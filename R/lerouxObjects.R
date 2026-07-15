#' Construct Leroux CAR objects from neighbourhood structures
#'
#' Construct the objects required by `dcar_leroux()` from different
#' neighbourhood representations.
#'
#' The function can take a binary adjacency matrix `W`, the usual `WinBUGS`
#' neighbourhood objects `adj` and `num`, or a graph object.
#'
#' The argument `W` should be a binary adjacency matrix. It can be created
#' manually or obtained by converting a neighbourhood object using, for example,
#' `spdep::nb2mat(Neigh, style = "B", zero.policy = TRUE)`, where `Neigh` is an
#' object produced by `spdep::poly2nb()`. If the `poly2nb()` neighbourhood
#' object itself is used directly, it should be passed through the `graph`
#' argument, not through `W`.
#'
#' The function returns the matrix of distinct neighbouring pairs `from.to`,
#' the number of distinct neighbouring pairs `NDist`, and the eigenvalues
#' `Lambda` of \eqn{\boldsymbol{D} - \boldsymbol{W}}, where
#' \eqn{\boldsymbol{D}} is the diagonal matrix of row sums of `W`.
#'
#' Exactly one neighbourhood representation must be provided: either `W`,
#' both `adj` and `num`, or `graph`.
#'
#' @param W Optional binary neighbourhood matrix. Rows and columns must
#'   correspond to the same spatial units. The matrix must be square,
#'   symmetric, contain only `0` and `1` values, and have a zero diagonal.
#' @param adj Optional vector containing the ID numbers of the adjacent areas
#'   for each area, using the usual `WinBUGS` neighbourhood representation.
#'   This argument must be provided together with `num`.
#' @param num Optional vector giving the number of neighbours of each area,
#'   using the usual `WinBUGS` neighbourhood representation. Its length must be
#'   equal to the total number of areas. This argument must be provided together
#'   with `adj`.
#' @param graph Optional graph representation. This can be an object of class
#'   `nb`, such as those produced by `spdep::poly2nb()`, a dense or sparse
#'   adjacency matrix used in INLA models, or an INLA-like graph object with a
#'   list of neighbours.
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
#' leroux.obj <- lerouxObjects(W = W)
#' leroux.obj$from.to
#' leroux.obj$Lambda
#' leroux.obj$NDist
#'
#' adj <- c(2, 1, 3, 2)
#' num <- c(1, 2, 1)
#'
#' leroux.obj <- lerouxObjects(adj = adj, num = num)
#' leroux.obj$from.to
#' leroux.obj$Lambda
#' leroux.obj$NDist
#'
#' graph <- list(nbs = list(2, c(1, 3), 2))
#' class(graph) <- "inla.graph"
#'
#' leroux.obj <- lerouxObjects(graph = graph)
#' leroux.obj$from.to
#' leroux.obj$Lambda
#' leroux.obj$NDist
#'
#' \dontrun{
#' Neigh <- spdep::poly2nb(cartography)
#'
#' W <- spdep::nb2mat(Neigh, style = "B", zero.policy = TRUE)
#' leroux.obj <- lerouxObjects(W = W)
#' leroux.obj$from.to
#' leroux.obj$Lambda
#' leroux.obj$NDist
#'
#' leroux.obj <- lerouxObjects(graph = Neigh)
#' leroux.obj$from.to
#' leroux.obj$Lambda
#' leroux.obj$NDist
#'
#' graph <- INLA::inla.read.graph("graph.file")
#' leroux.obj <- lerouxObjects(graph = graph)
#' leroux.obj$from.to
#' leroux.obj$Lambda
#' leroux.obj$NDist
#' }
#'
#' @export
lerouxObjects <- function(W = NULL, adj = NULL, num = NULL, graph = NULL) {

  # Check that exactly one neighbourhood representation is provided
  input.sources <- c(!is.null(W),
                     !is.null(adj) || !is.null(num),
                     !is.null(graph))

  if (sum(input.sources) != 1) {
    stop("Provide exactly one neighbourhood representation: W, adj and num, or graph.")
  }

  # Check that adj and num are provided together
  if (xor(is.null(adj), is.null(num))) {
    stop("adj and num must be provided together.")
  }

  # Internal function to check W
  checkW <- function(W) {

    if (!is.matrix(W)) {
      stop("W must be a matrix.")
    }

    if (nrow(W) != ncol(W)) {
      stop("W must be a square matrix.")
    }

    if (anyNA(W)) {
      stop("W must not contain missing values.")
    }

    if (!all(W %in% c(0, 1))) {
      stop("W must contain only 0 and 1 values.")
    }

    if (any(diag(W) != 0)) {
      stop("W must have a zero diagonal.")
    }

    if (!isTRUE(all.equal(W, t(W)))) {
      stop("W must be symmetric.")
    }

    return(W)
  }

  # Internal function to construct W from adj and num
  adjNumToW <- function(adj, num) {

    if (!is.numeric(num) || anyNA(num)) {
      stop("num must be a numeric vector without missing values.")
    }

    if (!is.numeric(adj) || anyNA(adj)) {
      stop("adj must be a numeric vector without missing values.")
    }

    if (any(num < 0) || any(num %% 1 != 0)) {
      stop("num must contain non-negative integers.")
    }

    if (length(adj) != sum(num)) {
      stop("length of adj must be equal to sum(num).")
    }

    NAreas <- length(num)

    if (length(adj) > 0 && any(adj < 1 | adj > NAreas | adj %% 1 != 0)) {
      stop("adj must contain integer area indices between 1 and length(num).")
    }

    W <- matrix(0, nrow = NAreas, ncol = NAreas)
    index <- c(0, cumsum(num))

    for (i in 1:NAreas) {
      if (num[i] > 0) {
        neighbours <- adj[(index[i] + 1):index[i + 1]]
        W[i, neighbours] <- 1
      }
    }

    return(W)
  }

  # Internal function to construct W from a graph object
  graphToW <- function(graph) {

    # Graph supplied as a matrix, including sparse adjacency matrices used in
    # INLA models
    if (is.matrix(graph) || inherits(graph, "Matrix")) {

      W <- as.matrix(graph)
      diag(W) <- 0
      W <- 1 * (W != 0)

      return(W)
    }

    # Graph supplied as an spdep nb object, such as from spdep::poly2nb()
    if (inherits(graph, "nb")) {

      NAreas <- length(graph)
      W <- matrix(0, nrow = NAreas, ncol = NAreas)

      for (i in 1:NAreas) {
        neighbours <- graph[[i]]
        neighbours <- neighbours[neighbours > 0]

        if (length(neighbours) > 0) {
          W[i, neighbours] <- 1
        }
      }

      return(W)
    }

    # Graph supplied as an INLA-like graph object with a list of neighbours
    if (is.list(graph) && !is.null(graph$nbs)) {

      NAreas <- length(graph$nbs)
      W <- matrix(0, nrow = NAreas, ncol = NAreas)

      for (i in 1:NAreas) {
        neighbours <- graph$nbs[[i]]
        neighbours <- neighbours[neighbours > 0]

        if (length(neighbours) > 0) {
          W[i, neighbours] <- 1
        }
      }

      return(W)
    }

    stop("graph must be an nb object, an adjacency matrix, or an INLA-like graph object with a list of neighbours.")
  }

  # Construct W from the selected input
  if (!is.null(W)) {

    W <- checkW(W)

  } else if (!is.null(adj) && !is.null(num)) {

    W <- adjNumToW(adj = adj, num = num)
    W <- checkW(W)

  } else if (!is.null(graph)) {

    W <- graphToW(graph)
    W <- checkW(W)
  }

  # Number of spatial units
  NAreas <- nrow(W)

  # Matrix D - W
  Q <- diag(rowSums(W), nrow = NAreas) - W

  # Number of neighbours of each area
  num <- rowSums(W)

  # Neighbours of each area
  adj <- unlist(lapply(1:NAreas, function(i) which(W[i, ] != 0)),
                use.names = FALSE)

  # All the neighbourhoods j ~ i where i < j
  from.to <- cbind(rep(1:NAreas, times = num), adj)
  colnames(from.to) <- c("from", "to")
  from.to <- from.to[from.to[, 1] < from.to[, 2], , drop = FALSE]

  # Number of distinct pairs of neighbours
  NDist <- nrow(from.to)

  # Eigenvalues of D - W
  Lambda <- eigen(Q)$values

  return(list(from.to = from.to,
              Lambda = Lambda,
              NDist = NDist))
}
