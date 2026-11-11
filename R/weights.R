#' Kernel weight function of Li, Wang and Chen (1999)
#'
#' Computes \eqn{\exp\{-\hat\Lambda(y_i, t \mid z_i)\}} for every observation
#' \eqn{i} and every threshold \eqn{t} of \code{t_grid}. This is the weight that
#' corrects the slice means for right censoring in CSIR and in the censored
#' weighted leverage score.
#'
#' The kernel smoothing is carried out on a low-dimensional projection \code{Z}
#' of the predictors, with a product Gaussian kernel and a Silverman-type
#' bandwidth. The denominator \eqn{\hat f} can be floored at \eqn{b_n}, equal to
#' \code{floor_mult} times \eqn{n^{-1/(r+4)}} times the modal density of a
#' Gaussian with the same scales, which is the floor required by condition
#' (C.9)(c) of the paper. The default \code{floor_mult = 0} is the estimator
#' used for every table of the paper.
#'
#' @param Z numeric matrix (n x r), the projection on which the kernel acts.
#' @param y numeric vector of observed times (any monotone transformation, for
#'   example the log time, gives the same result).
#' @param delta event indicator, 1 = event, 0 = censored.
#' @param t_grid numeric vector of thresholds.
#' @param h bandwidth; \code{NULL} uses a Silverman-type rule.
#' @param floor_mult multiplier of the floor \eqn{b_n}. \code{NULL} takes the
#'   option \code{csirwls.floor_mult}, then the environment variable
#'   \code{CSIRWLS_FLOOR_MULT}, then 0.
#' @return An \code{n x length(t_grid)} matrix \code{W} with
#'   \code{W[i, m]} \eqn{= \exp\{-\hat\Lambda(y_i, t_m \mid z_i)\}}.
#' @references Li, K.-C., Wang, J.-L. and Chen, C.-H. (1999). Dimension
#'   reduction for censored regression data. \emph{The Annals of Statistics},
#'   27, 1-23.
#' @examples
#' set.seed(1)
#' Z <- matrix(rnorm(100), 50, 2)
#' y <- rexp(50); d <- rbinom(50, 1, 0.7)
#' W <- lambda_weights(Z, y, d, t_grid = quantile(y, c(0.25, 0.5, 0.75)))
#' dim(W)
#' @export
lambda_weights <- function(Z, y, delta, t_grid, h = NULL, floor_mult = NULL) {
  if (is.null(floor_mult)) floor_mult <- default_floor_mult()
  Z <- as.matrix(Z)
  n <- length(y)
  if (is.null(h)) h <- silverman_bw(Z)

  ## kernel matrix K[i, j] = K_h(z_j - z_i), product Gaussian kernel
  D2 <- as.matrix(stats::dist(Z / h))^2
  K <- exp(-0.5 * D2)
  fhat <- rowSums(K)                              # n h^r fhat(z_i), up to a constant
  bn <- floor_mult * n^(3 / (ncol(Z) + 4))        # floor of condition (C.9)(c)
  fhat <- pmax(fhat, max(bn, 1e-12))

  ## S_hat(y_i | z_i) = sum_{j : y_j >= y_i} K[i, j] / sum_j K[i, j]
  ge <- outer(y, y, "<=")                         # ge[i, j] = y_j >= y_i
  S <- rowSums(K * ge) / fhat
  S <- pmin(pmax(S, 1e-6), 1)

  ## Lambda_hat(y_i, t | z_i)
  ##   = sum_{j : y_i < y_j <= t, delta_j = 1} K[i, j] / S(y_j | z_j) / fhat_i
  invS <- ifelse(delta == 1, 1 / S, 0)
  gt <- outer(y, y, "<")                          # gt[i, j] = y_j > y_i
  KG <- K * gt
  W <- matrix(0, n, length(t_grid))
  for (m in seq_along(t_grid)) {
    v <- invS * (y <= t_grid[m])                  # only uncensored, only y_j <= t
    W[, m] <- exp(-as.vector(KG %*% v) / fhat)
  }
  W
}

#' Censoring-corrected slice moments
#'
#' Slice probabilities \eqn{\hat p_\ell} and censoring-corrected slice means
#' \eqn{\bar U^c_\ell} of the columns of \code{U}, slicing the (unobserved)
#' survival time into \code{H} slices by quantiles of the observed time. The
#' contribution of a censored observation to the slices after its censoring time
#' is weighted by \code{\link{lambda_weights}}.
#'
#' @param U numeric matrix (n x k) whose slice means are wanted.
#' @param y observed times.
#' @param delta event indicator, 1 = event.
#' @param Z projection used for the kernel weights (n x r).
#' @param H number of slices.
#' @param h bandwidth passed to \code{\link{lambda_weights}}.
#' @param floor_mult floor multiplier passed to \code{\link{lambda_weights}}.
#' @return A list with \code{phat} (length \code{H}) and \code{Ubar}
#'   (\code{H x k}).
#' @export
censored_slice_moments <- function(U, y, delta, Z, H = 5L, h = NULL, floor_mult = NULL) {
  n <- length(y); k <- ncol(U)
  cuts <- slice_cuts(y, H)
  t_int <- cuts[2:(length(cuts) - 1L)]             # finite thresholds

  W <- if (length(t_int)) lambda_weights(Z, y, delta, t_int, h = h, floor_mult = floor_mult)
       else matrix(0, n, 0)

  ## E_hat[u 1{y0 >= t}] and P_hat[y0 >= t]; m indexes t within t_int
  tail_moments <- function(t, m) {
    if (is.infinite(t) && t < 0) return(list(mom = colSums(U) / n, prob = 1))
    if (is.infinite(t))          return(list(mom = rep(0, k),      prob = 0))
    obs <- y >= t
    cen <- (y < t) & (delta == 0)
    w <- W[, m]
    num <- colSums(U[obs, , drop = FALSE]) +
           colSums(U[cen, , drop = FALSE] * w[cen])
    den <- sum(obs) + sum(w[cen])
    list(mom = num / n, prob = den / n)
  }

  mom <- vector("list", length(cuts)); prob <- numeric(length(cuts))
  for (a in seq_along(cuts)) {
    m <- match(cuts[a], t_int, nomatch = NA_integer_)
    tm <- tail_moments(cuts[a], m)
    mom[[a]] <- tm$mom; prob[a] <- tm$prob
  }

  Ubar <- matrix(0, H, k); phat <- numeric(H)
  for (l in seq_len(H)) {
    p_l <- prob[l] - prob[l + 1L]
    num <- mom[[l]] - mom[[l + 1L]]
    if (p_l <= 1e-8) { Ubar[l, ] <- 0; phat[l] <- 0 } else { Ubar[l, ] <- num / p_l; phat[l] <- p_l }
  }
  list(phat = phat, Ubar = Ubar)
}
