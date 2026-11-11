#' Accuracy criteria for estimated directions
#'
#' \code{proximity} is the proximity criterion PC, the absolute cosine between an
#' estimated and a true direction. \code{subspace_metrics} returns the recovery
#' criterion \eqn{R = \mathrm{tr}(P P_{\hat B})/d} and the validity criterion
#' \eqn{V = \mathrm{tr}(P P_{\hat B})/\hat d} for subspaces of possibly different
#' dimensions (both equal \eqn{PC^2} when \eqn{d = \hat d = 1}).
#' \code{index_corr} is the index correlation IC, the absolute cosine in the
#' \eqn{\Sigma}-metric, i.e. the correlation between the estimated and the true
#' index, for a diagonal covariance with diagonal \code{Sigma_diag}.
#'
#' @param bhat,Bhat estimated direction(s).
#' @param btrue,Btrue true direction(s).
#' @param Sigma_diag diagonal of the covariance matrix of the predictors.
#' @return A number (\code{proximity}, \code{index_corr}) or a named vector
#'   (\code{subspace_metrics}).
#' @examples
#' proximity(c(1, 1, 0), c(1, 0, 0))
#' subspace_metrics(cbind(c(1, 0, 0), c(0, 1, 0)), c(1, 0, 0))
#' @export
proximity <- function(bhat, btrue) {
  a <- as.vector(bhat); b <- as.vector(btrue)
  a <- a / max(sqrt(sum(a^2)), 1e-300)
  b <- b / max(sqrt(sum(b^2)), 1e-300)
  abs(sum(a * b))
}

#' @rdname proximity
#' @export
subspace_metrics <- function(Bhat, Btrue) {
  Bh <- as.matrix(Bhat); Bt <- as.matrix(Btrue)
  Qh <- qr.Q(qr(Bh)); Qt <- qr.Q(qr(Bt))
  tr <- sum((crossprod(Qt, Qh))^2)
  c(recovery = tr / ncol(Bt), validity = tr / ncol(Bh))
}

#' @rdname proximity
#' @export
index_corr <- function(bhat, btrue, Sigma_diag) {
  a <- as.vector(bhat); b <- as.vector(btrue)
  num <- abs(sum(a * b * Sigma_diag))
  den <- sqrt(sum(a * a * Sigma_diag) * sum(b * b * Sigma_diag))
  num / max(den, 1e-300)
}

#' Screening criteria
#'
#' \code{mms} is the minimum model size: the number of top-ranked predictors
#' needed to include every active one. It does not depend on a stopping rule, so
#' it compares rankings at a matched model size. \code{screening_rates} gives the
#' true positive rate, false positive rate, false discovery rate and size of a
#' selected set.
#'
#' @param scores screening statistics (larger = more important).
#' @param active indices of the active predictors.
#' @param sel indices of the selected predictors.
#' @param p number of predictors.
#' @return \code{mms}: an integer; \code{screening_rates}: a named vector
#'   (\code{tpr}, \code{fpr}, \code{fdr}, \code{size}).
#' @examples
#' mms(c(5, 4, 1, 3, 2), active = c(1, 4))
#' screening_rates(c(1, 2, 4), active = c(1, 4), p = 5)
#' @export
mms <- function(scores, active) {
  rank_ <- order(scores, decreasing = TRUE)
  max(match(active, rank_))
}

#' @rdname mms
#' @export
screening_rates <- function(sel, active, p) {
  tp <- length(intersect(active, sel)); fp <- length(setdiff(sel, active))
  c(tpr = tp / length(active), fpr = fp / (p - length(active)),
    fdr = if (length(sel)) fp / length(sel) else 0, size = length(sel))
}

#' Predictive criteria for survival data
#'
#' \code{c_index} is Harrell's concordance index (larger risk = earlier event).
#' \code{km} is the Kaplan-Meier estimator. \code{breslow_surv} gives the
#' survival curves predicted by a linear risk score through the Breslow baseline
#' hazard. \code{integrated_brier} is Graf's integrated Brier score with inverse
#' probability of censoring weights. \code{logrank_test} is the k-group log-rank
#' test.
#'
#' @param risk risk score (larger = higher risk).
#' @param y,delta observed times and event indicators.
#' @param beta coefficients of the linear risk score.
#' @param Xtr,ytr,dtr training data.
#' @param Xte covariates of the observations whose curves are wanted.
#' @param t_grid time grid.
#' @param S matrix of predicted survival probabilities (observations x
#'   \code{t_grid}).
#' @param yte,dte test times and event indicators.
#' @param group factor or integer vector of groups.
#' @return \code{c_index}, \code{integrated_brier}: a number; \code{km}: a list
#'   with \code{times} and \code{surv}; \code{breslow_surv}: a matrix;
#'   \code{logrank_test}: a list with \code{statistic}, \code{df} and
#'   \code{p.value}.
#' @references Graf, E., Schmoor, C., Sauerbrei, W. and Schumacher, M. (1999).
#'   Assessment and comparison of prognostic classification schemes for survival
#'   data. \emph{Statistics in Medicine}, 18, 2529-2545.
#' @examples
#' set.seed(1)
#' dat <- simulate_csir_data(n = 200, p = 10, dgp = 3, censoring = 0.25)
#' b <- cox_fit(dat$X, dat$y, dat$delta, l2 = 1e-3)
#' c_index(dat$X %*% b, dat$y, dat$delta)
#' @export
c_index <- function(risk, y, delta) {
  risk <- as.vector(risk)
  comp <- outer(y, y, "<") & (delta == 1L)
  ii <- which(comp, arr.ind = TRUE)
  if (nrow(ii) == 0L) return(NA_real_)
  ri <- risk[ii[, 1]]; rj <- risk[ii[, 2]]
  (sum(ri > rj) + 0.5 * sum(ri == rj)) / nrow(ii)
}

#' @rdname c_index
#' @export
km <- function(y, delta) {
  t <- sort(unique(y[delta == 1L]))
  s <- 1; out <- numeric(length(t))
  for (a in seq_along(t)) {
    n_risk <- sum(y >= t[a])
    dd <- sum((y == t[a]) & (delta == 1L))
    s <- s * (1 - dd / max(n_risk, 1))
    out[a] <- s
  }
  list(times = t, surv = out)
}

#' @noRd
step_fun <- function(times, vals, t) {
  if (length(times) == 0L) return(rep(1, length(t)))
  idx <- findInterval(t, times)
  ifelse(idx >= 1L, vals[pmax(pmin(idx, length(vals)), 1L)], 1)
}

#' @rdname c_index
#' @export
breslow_surv <- function(beta, Xtr, ytr, dtr, Xte, t_grid) {
  eta_tr <- pmin(pmax(as.vector(Xtr %*% beta), -30), 30)
  e <- exp(eta_tr)
  evt <- sort(unique(ytr[dtr == 1L]))
  H0 <- numeric(length(evt)); acc <- 0
  for (a in seq_along(evt)) {
    acc <- acc + sum((ytr == evt[a]) & (dtr == 1L)) / max(sum(e[ytr >= evt[a]]), 1e-300)
    H0[a] <- acc
  }
  base <- step_fun(evt, H0, t_grid)
  base[!is.finite(base)] <- 0
  eta_te <- pmin(pmax(as.vector(Xte %*% beta), -30), 30)
  exp(-outer(exp(eta_te), base))
}

#' @rdname c_index
#' @export
integrated_brier <- function(S, t_grid, yte, dte, ytr, dtr) {
  g <- km(ytr, 1L - dtr)                                  # censoring KM
  G_t <- pmax(step_fun(g$times, g$surv, t_grid), 1e-3)
  G_y <- pmax(step_fun(g$times, g$surv, yte), 1e-3)
  bs <- numeric(length(t_grid))
  for (m in seq_along(t_grid)) {
    t <- t_grid[m]
    a <- (yte <= t) & (dte == 1L)
    b <- (yte > t)
    term <- ifelse(a, (0 - S[, m])^2 / G_y, ifelse(b, (1 - S[, m])^2 / G_t[m], 0))
    bs[m] <- mean(term)
  }
  ## trapezoidal average over the grid
  d <- diff(t_grid)
  sum(d * (bs[-1] + bs[-length(bs)]) / 2) / (t_grid[length(t_grid)] - t_grid[1])
}

#' @rdname c_index
#' @export
logrank_test <- function(y, delta, group) {
  g <- as.integer(factor(group)); G <- max(g)
  O <- numeric(G); E <- numeric(G); V <- matrix(0, G, G)
  for (t in sort(unique(y[delta == 1L]))) {
    n_j <- tabulate(g[y >= t], G)
    d_j <- tabulate(g[y == t & delta == 1L], G)
    N <- sum(n_j); D <- sum(d_j)
    if (N <= 1 || D == 0) next
    O <- O + d_j; E <- E + D * n_j / N
    V <- V + D * (N - D) / (N - 1) * (diag(n_j / N, G) - tcrossprod(n_j) / N^2)
  }
  U <- (O - E)[-G]
  Vr <- V[-G, -G, drop = FALSE]
  stat <- as.numeric(t(U) %*% solve(Vr, U))
  list(statistic = stat, df = G - 1L,
       p.value = stats::pchisq(stat, G - 1L, lower.tail = FALSE),
       observed = O, expected = E)
}
