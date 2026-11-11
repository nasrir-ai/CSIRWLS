#' Concordance screening statistic (CSIS)
#'
#' Marginal concordance utility \eqn{|c_j - 1/2|} of Cheng et al. (2024), where
#' \eqn{c_j} is the proportion of comparable pairs (an event followed by a later
#' observed time) in which the predictor \eqn{x_j} is larger for the later time.
#'
#' @param X numeric matrix (n x p).
#' @param y observed times.
#' @param delta event indicator, 1 = event.
#' @return The vector of the \code{p} screening statistics (larger = stronger).
#' @references Cheng, X., Li, G. and Wang, H. (2024). The concordance filter:
#'   an adaptive model-free feature screening procedure. \emph{Computational
#'   Statistics}, 39, 2413-2436.
#' @seealso \code{\link{cris}}, \code{\link{cmv_sis}}, \code{\link{cwls}}
#' @export
csis <- function(X, y, delta) {
  X <- check_xyd(X, y, delta)
  n <- nrow(X); p <- ncol(X)
  o <- order(y)
  Xs <- X[o, , drop = FALSE]; ys <- y[o]; ds <- delta[o]
  conc <- numeric(p); tot <- 0
  for (i in seq_len(n - 1L)) {
    if (ds[i] != 1L) next
    start <- findInterval(ys[i], ys) + 1L            # strictly later times
    if (start > n) next
    later <- Xs[start:n, , drop = FALSE]
    conc <- conc + colSums(sweep(later, 2L, Xs[i, ], "<"))
    tot <- tot + (n - start + 1L)
  }
  if (tot == 0) return(numeric(p))
  abs(conc / tot - 0.5)
}

## Kaplan-Meier estimate G(t-) = P(C >= t) of the censoring survival at every
## observed time (left limit); failures precede censorings at tied times.
#' @noRd
km_censoring <- function(y, delta) {
  n <- length(y)
  o <- order(y)
  ys <- y[o]; ds <- delta[o]
  Gs <- numeric(n); s <- 1; i <- 1L
  while (i <= n) {
    j <- i
    while (j <= n && ys[j] == ys[i]) j <- j + 1L
    Gs[i:(j - 1L)] <- s
    n_risk <- n - i + 1L
    n_cens <- sum(ds[i:(j - 1L)] == 0L)
    if (n_risk > 0) s <- s * (1 - n_cens / n_risk)
    i <- j
  }
  G <- numeric(n); G[o] <- Gs
  G
}

#' Censored rank independence screening (CRIS)
#'
#' Statistic of Song, Lu, Ma and Jeng (2014):
#' \eqn{\tau_j = \{n(n-1)\}^{-1}\sum_{i \ne k}\delta_i 1\{y_i < y_k\}
#' 1\{x_{ij} < x_{kj}\}/\hat G(y_i-)^2}, with \eqn{\hat G} the Kaplan-Meier
#' estimate of the censoring survival; predictors are ranked by
#' \eqn{|\tau_j - 1/4|}.
#'
#' @inheritParams csis
#' @return The vector of the \code{p} screening statistics.
#' @references Song, R., Lu, W., Ma, S. and Jeng, X. J. (2014). Censored rank
#'   independence screening of high-dimensional survival data.
#'   \emph{Biometrika}, 101, 799-814.
#' @export
cris <- function(X, y, delta) {
  X <- check_xyd(X, y, delta)
  n <- nrow(X); p <- ncol(X)
  G <- km_censoring(y, delta)
  o <- order(y)
  Xs <- X[o, , drop = FALSE]; ys <- y[o]; ds <- delta[o]; Gs <- G[o]
  acc <- numeric(p)
  for (i in seq_len(n - 1L)) {
    if (ds[i] != 1L) next
    start <- findInterval(ys[i], ys) + 1L
    if (start > n) next
    w <- 1 / max(Gs[i], 1e-8)^2
    acc <- acc + w * colSums(sweep(Xs[start:n, , drop = FALSE], 2L, Xs[i, ], ">"))
  }
  abs(acc / (n * (n - 1)) - 0.25)
}

#' @noRd
km_curve_sorted <- function(ds, member) {
  at_risk <- rev(cumsum(rev(as.numeric(member))))
  ev <- as.numeric(ds == 1L & member)
  haz <- ifelse(at_risk > 0, ev / pmax(at_risk, 1), 0)
  cumprod(1 - haz)
}

#' Censored mean-variance screening (cMV-SIS)
#'
#' Statistic of Zhong, Wang and Chen (2021):
#' \eqn{\sum_R\sum_r p_r\int\{S_r(t) - S(t)\}^2 dF(t)}, where \eqn{S_r} and
#' \eqn{S} are Kaplan-Meier curves within the \eqn{r}th of \eqn{R} slices of
#' \eqn{x_j} and overall, summed over the numbers of slices in \code{slices}.
#'
#' @inheritParams csis
#' @param slices numbers of slices of each predictor.
#' @return The vector of the \code{p} screening statistics.
#' @references Zhong, W., Wang, J. and Chen, X. (2021). Censored mean variance
#'   sure independence screening for ultrahigh dimensional survival data.
#'   \emph{Computational Statistics & Data Analysis}, 159, 107206.
#' @export
cmv_sis <- function(X, y, delta, slices = c(3L, 4L, 5L)) {
  X <- check_xyd(X, y, delta)
  n <- nrow(X); p <- ncol(X)
  o <- order(y)
  Xs <- X[o, , drop = FALSE]; ds <- delta[o]
  S_all <- km_curve_sorted(ds, rep(TRUE, n))
  S_prev <- c(1, S_all[-n])
  dF <- ifelse(ds == 1L, S_prev - S_all, 0)
  out <- numeric(p)
  for (j in seq_len(p)) {
    ranks <- rank(Xs[, j], ties.method = "first") - 1L
    for (R in slices) {
      grp <- (ranks * R) %/% n
      v <- 0
      for (r in 0:(R - 1L)) {
        member <- grp == r
        pr <- mean(member)
        if (pr == 0) next
        Sr <- km_curve_sorted(ds, member)
        v <- v + pr * sum(dF * (Sr - S_all)^2)
      }
      out[j] <- out[j] + v
    }
  }
  out
}
