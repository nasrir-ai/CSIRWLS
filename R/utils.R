## Small internal helpers shared by the estimators.

#' @noRd
center_cols <- function(X) sweep(X, 2L, colMeans(X), "-")

## Silverman-type bandwidth for an r-dimensional product Gaussian kernel
#' @noRd
silverman_bw <- function(Z) {
  n <- nrow(Z); r <- ncol(Z)
  s <- apply(Z, 2L, stats::sd)
  s[!(s > 0)] <- 1
  mean(s) * n^(-1 / (r + 4))
}

## H - 1 interior cut points by quantiles of y, padded with -Inf and Inf
#' @noRd
slice_cuts <- function(y, H) {
  q <- stats::quantile(y, seq(0, 1, length.out = H + 1L)[2:H], names = FALSE, type = 7)
  c(-Inf, q, Inf)
}

## symmetric inverse square root of a covariance matrix
#' @noRd
inv_sqrt <- function(S, ridge = 1e-6) {
  e <- eigen(S, symmetric = TRUE)
  w <- pmax(e$values, ridge)
  e$vectors %*% (t(e$vectors) * w^(-0.5))
}

## basic argument checks shared by the user-facing functions
#' @noRd
check_xyd <- function(X, y, delta = NULL) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.numeric(X)) stop("'X' must be a numeric matrix", call. = FALSE)
  if (anyNA(X)) stop("'X' contains missing values", call. = FALSE)
  if (length(y) != nrow(X)) stop("'y' must have nrow(X) elements", call. = FALSE)
  if (anyNA(y)) stop("'y' contains missing values", call. = FALSE)
  if (!is.null(delta)) {
    if (length(delta) != nrow(X)) stop("'delta' must have nrow(X) elements", call. = FALSE)
    if (!all(delta %in% c(0, 1))) stop("'delta' must be coded 1 = event, 0 = censored", call. = FALSE)
  }
  dimnames(X) <- NULL
  X
}

## default floor multiplier: option, then environment variable, then 0
#' @noRd
default_floor_mult <- function() {
  o <- getOption("csirwls.floor_mult")
  if (!is.null(o)) return(as.numeric(o))
  as.numeric(Sys.getenv("CSIRWLS_FLOOR_MULT", "0"))
}
