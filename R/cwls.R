#' Censored weighted leverage scores (CWLS)
#'
#' Computes the censored weighted leverage score \eqn{\hat\omega^c_j} of every
#' predictor. The design matrix is centred and decomposed as
#' \eqn{X = U D V^\top}; the score of predictor \eqn{j} is
#' \eqn{V_{(j)} \hat\Pi V_{(j)}^\top}, where \eqn{\hat\Pi} is built from the
#' censoring-corrected slice means of the rows of \eqn{U}
#' (\code{\link{censored_slice_moments}}). Unlike marginal screening statistics,
#' the score depends on all predictors jointly.
#'
#' @param X numeric matrix (n x p) of predictors.
#' @param y observed times (censored or not).
#' @param delta event indicator, 1 = event, 0 = censored.
#' @param H number of slices.
#' @param r dimension of the projection used by the kernel weight function.
#' @param k number of singular components retained (default \code{min(n, p)}).
#' @param h bandwidth of the kernel weight function (\code{NULL}: automatic).
#' @param floor_mult floor multiplier, see \code{\link{lambda_weights}}.
#' @param return_parts if \code{TRUE}, return the intermediate quantities too.
#' @return The vector of the \code{p} scores, or, when \code{return_parts = TRUE},
#'   a list with \code{omega}, \code{U}, \code{s}, \code{V}, \code{Pi},
#'   \code{phat} and \code{Ubar}.
#' @seealso \code{\link{bic_cut}}, \code{\link{csir_wls}}
#' @references Zhong, W., Liu, Y. and Zeng, P. (2023). A model-free variable
#'   screening method based on leverage score. \emph{Journal of the American
#'   Statistical Association}, 118, 135-146.
#' @examples
#' set.seed(1)
#' dat <- simulate_csir_data(n = 200, p = 100, dgp = 1, censoring = 0.25)
#' om <- cwls(dat$X, dat$y, dat$delta)
#' head(order(om, decreasing = TRUE))   # the active predictors are 1:5
#' @export
cwls <- function(X, y, delta, H = 5L, r = 2L, k = NULL, h = NULL,
                 floor_mult = NULL, return_parts = FALSE) {
  X <- check_xyd(X, y, delta)
  n <- nrow(X); p <- ncol(X)
  sv <- svd(center_cols(X))
  if (is.null(k)) k <- min(n, p)
  k <- min(k, length(sv$d))
  U <- sv$u[, 1:k, drop = FALSE]
  s <- sv$d[1:k]
  V <- sv$v[, 1:k, drop = FALSE]

  Z <- U[, 1:r, drop = FALSE] * rep(s[1:r], each = n)   # r-dim projection of x
  cm <- censored_slice_moments(U, y, delta, Z, H = H, h = h, floor_mult = floor_mult)

  Pi <- matrix(0, k, k)
  for (l in seq_along(cm$phat)) Pi <- Pi + cm$phat[l] * tcrossprod(cm$Ubar[l, ])

  omega <- rowSums((V %*% Pi) * V)
  if (return_parts) return(list(omega = omega, U = U, s = s, V = V, Pi = Pi,
                                phat = cm$phat, Ubar = cm$Ubar))
  omega
}

#' Stopping rule G_c for the CWLS ranking
#'
#' Modified BIC of Zhong et al. (2023), censored version \eqn{G_c(r)}: the number
#' of top-ranked predictors \eqn{r \le r_{\max}} minimising
#' \eqn{-\log\sum_{j \le r}\omega_{(j)} + r\{\log n + c_{n2}\log p\}/\max(n, p)}.
#'
#' @param omega vector of scores (for example from \code{\link{cwls}}).
#' @param n sample size.
#' @param p number of predictors.
#' @param c_n2 constant multiplying \eqn{\log p} in the penalty.
#' @param rmax largest admissible size (default \eqn{\min(p, \lfloor n/\log n\rfloor)}).
#' @return The selected number of predictors (an integer).
#' @examples
#' set.seed(1)
#' dat <- simulate_csir_data(n = 200, p = 100, dgp = 1, censoring = 0.25)
#' om <- cwls(dat$X, dat$y, dat$delta)
#' bic_cut(om, n = 200, p = 100)
#' @export
bic_cut <- function(omega, n, p, c_n2 = 1, rmax = NULL) {
  o <- sort(omega, decreasing = TRUE)
  if (is.null(rmax)) rmax <- min(p, as.integer(n / log(n)))
  pen <- (log(n) + c_n2 * log(p)) / max(n, p)
  cs <- cumsum(o[1:rmax])
  vals <- -log(pmax(cs, 1e-300)) + seq_len(rmax) * pen
  which.min(vals)
}

#' CSIR-WLS: censored SIR after weighted leverage screening
#'
#' The two-stage procedure of the paper. The predictors are ranked by the
#' censored weighted leverage score (\code{\link{cwls}}), the first
#' \code{n_keep} are retained (by default chosen by \code{\link{bic_cut}}), and
#' censored sliced inverse regression (\code{\link{csir}}) estimates the EDR
#' directions from the retained predictors.
#'
#' @inheritParams cwls
#' @param d number of EDR directions to estimate.
#' @param n_keep number of predictors to retain; \code{NULL} uses the stopping
#'   rule \eqn{G_c} of \code{\link{bic_cut}}. At least \code{d + 1} predictors
#'   are always retained.
#' @param c_n2 constant of the stopping rule, see \code{\link{bic_cut}}.
#' @return An object of class \code{"csirwls"}: a list with
#'   \item{beta}{a \code{p x d} matrix of unit-norm directions, zero outside the
#'     retained set;}
#'   \item{A}{the sorted indices of the retained predictors;}
#'   \item{omega}{the CWLS scores;}
#'   \item{n_keep}{the number of retained predictors.}
#' @references Li, K.-C., Wang, J.-L. and Chen, C.-H. (1999). Dimension
#'   reduction for censored regression data. \emph{The Annals of Statistics},
#'   27, 1-23.
#' @seealso \code{\link{cwls}}, \code{\link{csir}}, \code{\link{cv_survival}}
#' @examples
#' set.seed(1)
#' dat <- simulate_csir_data(n = 200, p = 250, dgp = 1, censoring = 0.25)
#' fit <- csir_wls(dat$X, dat$y, dat$delta)
#' fit
#' proximity(fit$beta, dat$beta)
#' @export
csir_wls <- function(X, y, delta, H = 5L, r = 2L, d = 1L, n_keep = NULL, c_n2 = 1,
                     floor_mult = NULL) {
  X <- check_xyd(X, y, delta)
  n <- nrow(X); p <- ncol(X)
  omega <- cwls(X, y, delta, H = H, r = r, floor_mult = floor_mult)
  if (is.null(n_keep)) n_keep <- bic_cut(omega, n, p, c_n2 = c_n2)
  A <- order(omega, decreasing = TRUE)[seq_len(min(p, max(n_keep, d + 1L)))]
  B_A <- csir(X[, A, drop = FALSE], y, delta, H = H, d = d)
  B <- matrix(0, p, d)
  B[A, ] <- B_A
  structure(list(beta = B, A = sort(A), omega = omega, n_keep = length(A)),
            class = "csirwls")
}

#' @export
print.csirwls <- function(x, digits = 3, ...) {
  cat("CSIR-WLS fit:", nrow(x$beta), "predictors,", length(x$A), "retained,",
      ncol(x$beta), if (ncol(x$beta) > 1L) "directions\n" else "direction\n")
  top <- x$A[order(-rowSums(abs(x$beta[x$A, , drop = FALSE])))]
  show <- utils::head(top, 10L)
  B <- x$beta[show, , drop = FALSE]
  dimnames(B) <- list(paste0("x", show), paste0("dir", seq_len(ncol(B))))
  cat("Largest loadings:\n")
  print(round(B, digits))
  invisible(x)
}

#' Estimated indices of a CSIR-WLS fit
#'
#' @param object a \code{"csirwls"} object.
#' @param newdata numeric matrix with the same columns as the fitted data.
#' @param ... unused.
#' @return The matrix of estimated indices \code{newdata \%*\% beta}.
#' @export
predict.csirwls <- function(object, newdata, ...) {
  newdata <- as.matrix(newdata)
  if (ncol(newdata) != nrow(object$beta)) stop("'newdata' has the wrong number of columns")
  newdata %*% object$beta
}
