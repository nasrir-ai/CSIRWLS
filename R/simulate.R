## Data-generating processes of the simulation study of the paper.

.chol_cache <- new.env(parent = emptyenv())

#' @noRd
ar1_chol <- function(p, rho) {
  key <- paste0("ar1_", p, "_", rho)
  if (is.null(.chol_cache[[key]])) {
    i <- seq_len(p)
    S <- rho^abs(outer(i, i, "-"))
    .chol_cache[[key]] <- chol(S + 1e-10 * diag(p))
  }
  .chol_cache[[key]]
}

#' @noRd
block_chol <- function(p, alpha, s = 5L) {
  key <- paste0("blk_", p, "_", paste(alpha, collapse = "_"))
  if (is.null(.chol_cache[[key]])) {
    .chol_cache[[key]] <- chol(block_sigma(p, alpha, s) + 1e-9 * diag(p))
  }
  .chol_cache[[key]]
}

#' Simulation designs of the paper
#'
#' Covariate designs used in the simulation study.
#' \describe{
#'   \item{\code{design_ar1}}{Gaussian, AR(1) correlation \eqn{\rho^{|i-j|}}:
#'     Settings (1) (\eqn{\rho = 0.5}) and (2) (\eqn{\rho = 0.8});}
#'   \item{\code{design_block}}{Gaussian, block-exchangeable correlation with
#'     parameters \eqn{\alpha = (\alpha_1, \alpha_2, \alpha_3)} (within the first
#'     five predictors, between them and the others, among the others):
#'     Setting (5);}
#'   \item{\code{design_spiked}}{independent Gaussian columns with the spiked
#'     standard deviations of \code{spiked_lambda}: Settings (3) and (4);}
#'   \item{\code{design_skew}}{skew-normal (\code{kind = "skewnormal"}) or skew-t
#'     (\code{kind = "skewt"}) designs of Azzalini type built on the AR(1)
#'     design, used in the robustness study.}
#' }
#'
#' @param n sample size.
#' @param p number of predictors.
#' @param rho AR(1) correlation.
#' @param alpha block-exchangeable correlations.
#' @param setting 3 or 4.
#' @param kind \code{"skewnormal"} or \code{"skewt"}.
#' @param delta skewness parameter of the Azzalini designs.
#' @param nu degrees of freedom of the skew-t design.
#' @param s size of the first block.
#' @return An \code{n x p} matrix (\code{block_sigma}: the \code{p x p}
#'   covariance matrix; \code{spiked_lambda}: the \code{p} standard deviations).
#' @examples
#' X <- design_ar1(100, 20, rho = 0.5)
#' round(cor(X)[1:3, 1:3], 2)
#' @export
design_ar1 <- function(n, p, rho) matrix(stats::rnorm(n * p), n, p) %*% ar1_chol(p, rho)

#' @rdname design_ar1
#' @export
design_block <- function(n, p, alpha) matrix(stats::rnorm(n * p), n, p) %*% block_chol(p, alpha)

#' @rdname design_ar1
#' @export
block_sigma <- function(p, alpha, s = 5L) {
  S <- matrix(alpha[3], p, p)
  S[1:s, 1:s] <- alpha[1]
  S[1:s, (s + 1):p] <- alpha[2]
  S[(s + 1):p, 1:s] <- alpha[2]
  diag(S) <- 1
  S
}

#' @rdname design_ar1
#' @export
spiked_lambda <- function(n, p, setting) {
  cc <- ceiling(p / sqrt(n))
  top <- if (setting == 3L) 70 else 50
  c(seq(top + cc, cc, by = -1), rep(1, p - (top + 1)))
}

#' @rdname design_ar1
#' @export
design_spiked <- function(n, p, setting) {
  lam <- spiked_lambda(n, p, setting)
  matrix(stats::rnorm(n * p), n, p) * rep(lam, each = n)
}

#' @rdname design_ar1
#' @export
design_skew <- function(n, p, kind, rho = 0.5, delta = 0.8, nu = 5) {
  u <- design_ar1(n, p, rho)
  z0 <- abs(stats::rnorm(n))
  x <- delta * z0 + sqrt(1 - delta^2) * u
  if (kind == "skewt") x <- x / sqrt(stats::rchisq(n, nu) / nu)
  x
}

#' Survival times and censoring of the simulation study
#'
#' \code{log_time} draws the log survival time given the linear index
#' \code{eta}: DGP 1 is a log-normal accelerated failure time model, DGP 2 a
#' Weibull accelerated failure time model, DGP 3 a proportional hazards model.
#' Everything is generated and analysed on the log scale; every estimator of the
#' package is invariant to monotone transformations of the time.
#' \code{censor_log} draws \eqn{\log C = \mathrm{const} + N(0, 1)}.
#' \code{calibrate_const} finds by bisection the constant giving a target
#' censoring rate for a generator \code{gen()} of log survival times.
#' \code{draw_beta} draws the active coefficients uniformly on \eqn{[0.5, 1.5]}.
#'
#' @param dgp 1, 2 or 3.
#' @param eta linear index \eqn{\beta^\top x}.
#' @param n sample size.
#' @param const constant of the log censoring time.
#' @param gen function without arguments returning a vector of log survival
#'   times.
#' @param target target censoring rate.
#' @param tries,lo,hi,reps bisection settings.
#' @param p number of predictors.
#' @param active indices of the active predictors.
#' @return A numeric vector (\code{calibrate_const}: a number).
#' @export
log_time <- function(dgp, eta) {
  n <- length(eta)
  if (dgp == 1L) return(eta + stats::rnorm(n))
  if (dgp == 2L) return(-2.5 + eta + log(-log(1 - stats::runif(n))))
  if (dgp == 3L) return(log(-log(stats::runif(n))) - eta)
  stop("dgp must be 1, 2 or 3")
}

#' @rdname log_time
#' @export
draw_beta <- function(p, active = 1:5) {
  b <- numeric(p); b[active] <- stats::runif(length(active), 0.5, 1.5); b
}

#' @rdname log_time
#' @export
censor_log <- function(n, const) const + stats::rnorm(n)

#' @rdname log_time
#' @export
calibrate_const <- function(gen, target, tries = 28L, lo = -1500, hi = 1500, reps = 6L) {
  for (k in seq_len(tries)) {
    mid <- 0.5 * (lo + hi)
    r <- replicate(reps, {
      ly0 <- gen()
      mean(ly0 > censor_log(length(ly0), mid))
    })
    if (mean(r) > target) lo <- mid else hi <- mid
  }
  0.5 * (lo + hi)
}

#' Simulate censored data from the designs of the paper
#'
#' Convenience wrapper: draws a design, coefficients on the active set, a log
#' survival time from DGP 1, 2 or 3 and an independent log censoring time
#' calibrated to the requested censoring rate.
#'
#' @param n sample size.
#' @param p number of predictors.
#' @param dgp 1, 2 or 3 (see \code{\link{log_time}}).
#' @param design \code{"ar1"}, \code{"block"}, \code{"spiked3"} or
#'   \code{"spiked4"}.
#' @param rho AR(1) correlation.
#' @param alpha block-exchangeable correlations.
#' @param censoring target censoring rate.
#' @return A list with \code{X}, \code{y} (observed log time), \code{delta},
#'   \code{beta} and \code{active}.
#' @examples
#' set.seed(1)
#' dat <- simulate_csir_data(n = 200, p = 50)
#' mean(1 - dat$delta)
#' @export
simulate_csir_data <- function(n = 500L, p = 80L, dgp = 1L, design = c("ar1", "block", "spiked3", "spiked4"),
                               rho = 0.5, alpha = c(0.5, 0.2, 0.4), censoring = 0.25) {
  design <- match.arg(design)
  active <- switch(design, spiked3 = c(1, 10, 15, 20, 25, 30),
                   spiked4 = c(1, 10, 15, 20, 25, 30, 40, 45, 50), 1:5)
  mk <- function() switch(design,
    ar1 = design_ar1(n, p, rho), block = design_block(n, p, alpha),
    spiked3 = design_spiked(n, p, 3L), spiked4 = design_spiked(n, p, 4L))
  gen <- function() log_time(dgp, as.vector(mk() %*% draw_beta(p, active)))
  const <- calibrate_const(gen, censoring)
  X <- mk(); b <- draw_beta(p, active)
  ly0 <- log_time(dgp, as.vector(X %*% b)); lc <- censor_log(n, const)
  list(X = X, y = pmin(ly0, lc), delta = as.integer(ly0 <= lc), beta = b, active = active)
}
