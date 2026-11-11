## Breslow negative partial log-likelihood (divided by n) and its gradient
#' @noRd
cox_nll_grad <- function(beta, X, y, delta) {
  n <- length(y)
  o <- order(y, decreasing = TRUE)
  Xo <- X[o, , drop = FALSE]; do_ <- delta[o]
  eta <- pmin(pmax(as.vector(Xo %*% beta), -30), 30)
  e <- exp(eta)
  risk <- pmax(cumsum(e), 1e-300)                       # sum over y_j >= y_i
  nll <- -sum(do_ * (eta - log(risk)))
  wx <- apply(Xo * e, 2L, cumsum) / risk
  ie <- which(do_ == 1L)
  grad <- -colSums(Xo[ie, , drop = FALSE] - wx[ie, , drop = FALSE])
  list(nll = nll / n, grad = grad / n)
}

#' Cox proportional hazards fit, with optional L1 and L2 penalties
#'
#' Minimises the Breslow negative partial log-likelihood (divided by \eqn{n})
#' plus \eqn{l_1\lVert\beta\rVert_1 + (l_2/2)\lVert\beta\rVert_2^2} by FISTA
#' with backtracking. With \code{l1 = 0} and a small \code{l2} this is the
#' classical Cox fit used as a benchmark in the paper. No dependency on the
#' \pkg{survival} package.
#'
#' @param X numeric matrix (n x p).
#' @param y observed times.
#' @param delta event indicator.
#' @param l1 lasso penalty.
#' @param l2 ridge penalty.
#' @param n_iter maximal number of iterations.
#' @param step0 initial step size.
#' @param warm optional starting value.
#' @return The vector of the \code{p} coefficients.
#' @seealso \code{\link{cox_lasso_bic}}
#' @examples
#' set.seed(1)
#' dat <- simulate_csir_data(n = 300, p = 10, dgp = 3, censoring = 0.25)
#' round(cox_fit(dat$X, dat$y, dat$delta, l2 = 1e-3), 2)
#' @export
cox_fit <- function(X, y, delta, l1 = 0, l2 = 1e-4, n_iter = 400L, step0 = 1, warm = NULL) {
  X <- check_xyd(X, y, delta)
  p <- ncol(X)
  b <- if (is.null(warm)) numeric(p) else warm
  z <- b; tk <- 1; step <- step0; f_prev <- Inf
  cand <- b; fc <- Inf
  for (it in seq_len(n_iter)) {
    fg <- cox_nll_grad(z, X, y, delta)
    f <- fg$nll; g <- fg$grad + l2 * z
    ok <- FALSE
    for (bt in 1:10) {
      cand <- z - step * g
      if (l1 > 0) cand <- sign(cand) * pmax(abs(cand) - step * l1, 0)
      fc <- cox_nll_grad(cand, X, y, delta)$nll
      if (fc <= f + sum(g * (cand - z)) + sum((cand - z)^2) / (2 * step) + 1e-12) { ok <- TRUE; break }
      step <- step * 0.5
    }
    if (ok) step <- step * 1.3                          # let the step grow back
    b_new <- cand
    tk_new <- (1 + sqrt(1 + 4 * tk^2)) / 2
    z <- b_new + ((tk - 1) / tk_new) * (b_new - b)
    if (max(abs(b_new - b)) < 1e-7) { b <- b_new; break }
    b <- b_new; tk <- tk_new
    if (abs(f_prev - fc) < 1e-11) break
    f_prev <- fc
  }
  b
}

#' Cox-lasso with a BIC choice of the penalty
#'
#' Computes a Cox-lasso path on a logarithmic grid of \code{n_lambda} penalties,
#' from the smallest penalty giving the null model down to one hundredth of it,
#' with warm starts, and returns the fit minimising
#' \eqn{2n\,\mathrm{nll} + k\log(\#\mathrm{events})}, \eqn{k} being the number of
#' non-zero coefficients (fits with more than \eqn{\min(n/2, p)} non-zero
#' coefficients are discarded).
#'
#' @inheritParams cox_fit
#' @param n_lambda length of the penalty grid.
#' @return The vector of the \code{p} coefficients (all zero if no fit on the
#'   path is admissible).
#' @export
cox_lasso_bic <- function(X, y, delta, n_lambda = 12L, l2 = 1e-4, n_iter = 120L) {
  X <- check_xyd(X, y, delta)
  n <- nrow(X); p <- ncol(X)
  lmax <- max(abs(cox_nll_grad(numeric(p), X, y, delta)$grad))
  grid <- lmax * 10^seq(0, -2, length.out = n_lambda)
  best <- numeric(p); best_bic <- Inf; warm <- numeric(p)
  for (lam in grid) {
    b <- cox_fit(X, y, delta, l1 = lam, l2 = l2, n_iter = n_iter, warm = warm)
    warm <- b
    k <- sum(abs(b) > 1e-8)
    if (k == 0 || k > min(n %/% 2, p)) next
    nll <- cox_nll_grad(b, X, y, delta)$nll
    bic <- 2 * n * nll + k * log(max(sum(delta), 2))
    if (bic < best_bic) { best_bic <- bic; best <- b }
  }
  best
}
