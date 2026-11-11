## Kernel matrix of double-slicing censored SIR (internal): returns the
## standardisation, the kernel G and the number of non-empty slices.
#' @noRd
csir_kernel <- function(X, y, delta, H = 5L, ridge = 1e-6) {
  n <- nrow(X); p <- ncol(X)
  Xc <- center_cols(X)
  S <- stats::cov(Xc) + ridge * diag(p)
  Sih <- inv_sqrt(S, ridge)
  Xt <- Xc %*% Sih

  G <- matrix(0, p, p); nslices <- 0L
  for (dd in c(1L, 0L)) {                          # double slicing on delta
    m <- (delta == dd)
    if (sum(m) < H) next
    yy <- y[m]
    cuts <- stats::quantile(yy, seq(0, 1, length.out = H + 1L), names = FALSE, type = 7)
    cuts[1] <- -Inf; cuts[length(cuts)] <- Inf
    for (l in seq_len(H)) {
      sl <- m & (y > cuts[l]) & (y <= cuts[l + 1L])
      if (!any(sl)) next
      nslices <- nslices + 1L
      mu <- colMeans(Xt[sl, , drop = FALSE])
      G <- G + (sum(sl) / n) * tcrossprod(mu)
    }
  }
  list(G = G, Sih = Sih, nslices = nslices)
}

#' Censored sliced inverse regression (CSIR)
#'
#' Double-slicing censored SIR of Li, Wang and Chen (1999): the observations are
#' split by the event indicator, each group is sliced by quantiles of the
#' observed time, and the directions are the leading eigenvectors of the
#' weighted covariance of the standardised slice means, mapped back to the
#' original scale. Requires \code{p < n}; in high dimension use
#' \code{\link{csir_wls}}.
#'
#' @param X numeric matrix (n x p).
#' @param y observed times.
#' @param delta event indicator, 1 = event, 0 = censored.
#' @param H number of slices within each group.
#' @param d number of directions.
#' @param ridge ridge added to the covariance matrix before inversion.
#' @return A \code{p x d} matrix of unit-norm directions.
#' @references Li, K.-C., Wang, J.-L. and Chen, C.-H. (1999). Dimension
#'   reduction for censored regression data. \emph{The Annals of Statistics},
#'   27, 1-23.
#' @seealso \code{\link{csir_dim_test}}, \code{\link{csir_wls}}
#' @examples
#' set.seed(1)
#' dat <- simulate_csir_data(n = 300, p = 10, dgp = 1, censoring = 0.25)
#' proximity(csir(dat$X, dat$y, dat$delta), dat$beta)
#' @export
csir <- function(X, y, delta, H = 5L, d = 1L, ridge = 1e-6) {
  X <- check_xyd(X, y, delta)
  ker <- csir_kernel(X, y, delta, H = H, ridge = ridge)
  ev <- eigen(ker$G, symmetric = TRUE)
  B <- ker$Sih %*% ev$vectors[, 1:d, drop = FALSE]
  sweep(B, 2L, sqrt(colSums(B^2)), "/")
}

#' Sequential chi-squared test for the number of CSIR directions
#'
#' The sequential test of Li (1991), applied to the eigenvalues
#' \eqn{\hat\lambda_1 \ge \dots \ge \hat\lambda_p} of the double-slicing CSIR
#' kernel. For \eqn{k = 0, 1, \dots}, the statistic
#' \eqn{n\sum_{j>k}\hat\lambda_j} is compared with a chi-squared distribution on
#' \eqn{(p-k)(h-k-1)} degrees of freedom, \eqn{h} being the number of non-empty
#' slices; \eqn{\hat K} is the first \eqn{k} for which \eqn{H_0: K = k} is not
#' rejected. The asymptotics require \eqn{p} fixed and \eqn{n \to \infty}
#' (normal predictors in Li, 1991); the test should not be used when \eqn{p} is
#' comparable to \eqn{n}.
#'
#' @inheritParams csir
#' @param alpha level of each test.
#' @param kmax largest dimension tested.
#' @return A list with \code{K} (the estimated dimension) and \code{table}, a data
#'   frame with columns \code{k}, \code{statistic}, \code{df} and \code{p.value}.
#' @references Li, K.-C. (1991). Sliced inverse regression for dimension
#'   reduction. \emph{Journal of the American Statistical Association}, 86,
#'   316-342.
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   pbc <- survival::pbc
#'   pbc <- pbc[stats::complete.cases(pbc), ]
#'   X <- scale(as.matrix(pbc[, c("age", "edema", "bili", "albumin",
#'                                "platelet", "protime")]))
#'   csir_dim_test(X, pbc$time, as.integer(pbc$status == 2), alpha = 0.01)
#' }
#' @export
csir_dim_test <- function(X, y, delta, H = 5L, alpha = 0.05, kmax = NULL, ridge = 1e-6) {
  X <- check_xyd(X, y, delta)
  n <- nrow(X); p <- ncol(X)
  ker <- csir_kernel(X, y, delta, H = H, ridge = ridge)
  lam <- pmax(eigen(ker$G, symmetric = TRUE, only.values = TRUE)$values, 0)
  h <- ker$nslices
  if (is.null(kmax)) kmax <- min(p, h - 1L) - 1L
  tab <- data.frame(k = 0:kmax, statistic = NA_real_, df = NA_real_, p.value = NA_real_)
  for (i in seq_len(nrow(tab))) {
    k <- tab$k[i]
    tab$statistic[i] <- n * sum(lam[(k + 1L):p])
    tab$df[i] <- (p - k) * (h - k - 1L)
    tab$p.value[i] <- stats::pchisq(tab$statistic[i], tab$df[i], lower.tail = FALSE)
  }
  acc <- which(tab$p.value >= alpha)
  K <- if (length(acc)) tab$k[acc[1L]] else kmax + 1L
  list(K = K, table = tab, eigenvalues = lam, slices = h)
}

#' Sliced inverse regression on the observed time
#'
#' Classical SIR of Li (1991) applied to the observed (possibly censored) time,
#' ignoring the event indicator. Used as a naive benchmark.
#'
#' @param X numeric matrix (n x p).
#' @param y observed times.
#' @param H number of slices.
#' @param d number of directions.
#' @param ridge ridge added to the covariance matrix before inversion.
#' @return A \code{p x d} matrix of unit-norm directions.
#' @references Li, K.-C. (1991). Sliced inverse regression for dimension
#'   reduction. \emph{Journal of the American Statistical Association}, 86,
#'   316-342.
#' @export
sir <- function(X, y, H = 5L, d = 1L, ridge = 1e-6) {
  X <- check_xyd(X, y)
  n <- nrow(X); p <- ncol(X)
  Xc <- center_cols(X)
  Sih <- inv_sqrt(stats::cov(Xc) + ridge * diag(p), ridge)
  Xt <- Xc %*% Sih
  cuts <- stats::quantile(y, seq(0, 1, length.out = H + 1L), names = FALSE, type = 7)
  cuts[1] <- -Inf; cuts[length(cuts)] <- Inf
  G <- matrix(0, p, p)
  for (l in seq_len(H)) {
    sl <- (y > cuts[l]) & (y <= cuts[l + 1L])
    if (!any(sl)) next
    mu <- colMeans(Xt[sl, , drop = FALSE])
    G <- G + (sum(sl) / n) * tcrossprod(mu)
  }
  ev <- eigen(G, symmetric = TRUE)
  B <- Sih %*% ev$vectors[, 1:d, drop = FALSE]
  sweep(B, 2L, sqrt(colSums(B^2)), "/")
}

#' Counting-process SIR (CP-SIR)
#'
#' Algorithm 2 of Sun, Zhu, Wang and Zeng (2019). The covariates are
#' standardised; for every failure time \eqn{Y_i},
#' \eqn{\hat\phi(Y_i)} is the mean of \eqn{z_j} over the failures in
#' \eqn{[Y_i, Y_i + h)} minus the mean over the risk set, and the directions are
#' the leading left singular vectors of
#' \eqn{n^{-1}\sum_{\delta_i = 1}(z_i - \bar z_{R(Y_i)})\hat\phi(Y_i)^\top},
#' mapped back through the standardisation. Requires \code{p < n}.
#'
#' @param X numeric matrix (n x p).
#' @param y observed times.
#' @param delta event indicator.
#' @param d number of directions.
#' @param h window width on the log-time scale (\code{NULL}: rule of thumb).
#' @param ridge ridge added to the covariance matrix.
#' @param time_is_log set to \code{TRUE} when \code{y} is already on the log scale.
#' @return A \code{p x d} matrix of unit-norm directions.
#' @references Sun, Q., Zhu, R., Wang, T. and Zeng, D. (2019). Counting process
#'   based dimension reduction methods for censored outcomes.
#'   \emph{Biometrika}, 106, 181-196.
#' @export
cp_sir <- function(X, y, delta, d = 1L, h = NULL, ridge = 1e-6, time_is_log = FALSE) {
  X <- check_xyd(X, y, delta)
  n <- nrow(X); p <- ncol(X)
  Xc <- center_cols(X)
  Sih <- inv_sqrt(stats::cov(Xc) + ridge * diag(p), ridge)
  Z <- Xc %*% Sih
  ly <- if (time_is_log) as.numeric(y) else log(pmax(y, 1e-300))
  ev <- delta == 1
  if (is.null(h)) {
    s_ev <- if (sum(ev) > 1) stats::sd(ly[ev]) else 1
    h <- (4 / 3)^0.2 * sum(ev)^(-0.2) * s_ev
  }
  o <- order(ly)
  Zs <- Z[o, , drop = FALSE]; lys <- ly[o]; ds <- ev[o]

  ## suffix means for the at-risk sets {j : Y_j >= Y_i}
  csum <- apply(Zs[n:1, , drop = FALSE], 2L, cumsum)[n:1, , drop = FALSE]
  risk_mean <- csum / (n:1)

  ev_idx <- which(ds)
  le <- lys[ev_idx]
  ce <- apply(Zs[ev_idx, , drop = FALSE], 2L, cumsum)
  M <- matrix(0, p, p)
  for (i in ev_idx) {
    a <- findInterval(lys[i], le, left.open = TRUE)              # count of le < lys[i]
    b <- findInterval(lys[i] + h, le, left.open = TRUE)          # count of le < lys[i] + h
    if (b <= a) next
    win_mean <- (ce[b, ] - if (a > 0) ce[a, ] else 0) / (b - a)
    phi <- win_mean - risk_mean[i, ]
    M <- M + tcrossprod(Zs[i, ] - risk_mean[i, ], phi)
  }
  M <- M / n
  U <- svd(M)$u
  B <- Sih %*% U[, 1:d, drop = FALSE]
  sweep(B, 2L, sqrt(colSums(B^2)), "/")
}
