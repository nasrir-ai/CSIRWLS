## =============================================================================
## designs.R -- data-generating processes of Section 5 and of the Supplement
##
## Everything is generated on the log-time scale (as in R/simulate.R): every
## estimator used here is invariant to a monotone transformation of the time.
##
##   design (covariates)            spec_*()            Section 5 / Supplement
##   AR(1), rho = 0.5 / 0.8         spec_ar1(rho)       Settings (1)-(2)
##   spiked, Upsilon diagonal       spec_spiked(3|4)    Settings (3)-(4)
##   block-exchangeable             spec_block(alpha)   Setting (5), Table S12
##   skew-normal / skew-t (Azzalini) spec_skew(kind)    Sect. 5.8.2, S18-S21
##   Gaussian with an Azzalini covariance  spec_ellctrl(base)   Table 7 (controls)
##   independent skew components, AR(1)    spec_nonellctrl(base) Table 7 (controls)
##   block design, two negative signs      spec_block(alpha, neg = 2)  Table 8
##
##   model (log survival time)                       model =
##   DGP 5.1  eta + N(0,1)              (lognormal AFT)   "dgp1"
##   DGP 5.2  -2.5 + eta + log(-log(1-U))  (Weibull AFT) "dgp2"
##   DGP 5.3  log(-log U) - eta         (PH, exponential) "dgp3"
##   DGP 5.4  3 z1/(0.5+(z2+1.5)^2)+0.5 e (d = 2)       "dgp4"   (*)
##   Table S13 links on the standardised index z:
##            3 tanh(z) + 0.5 e                           "saturating"
##            2.5 1(z > 0) + 0.5 e                        "threshold"
##            z + 0.5 exp(z/2) e                          "heteroscedastic"
##
##   censoring  log C = c + shift + N(0,1)           mech =
##            shift = 0                                   "indep"
##            shift = 0.6 * standardised beta'x           "edr"
##            shift = 0.8 * standardised sum_{41..45} x_j "misspec"
##
## (*) DGP 5.4 (two-dimensional EDR space): log y0 = 3 z1 / (0.5 + (z2 + 1.5)^2)
##     + sigma e, sigma = 0.5, with z_k = beta_k'x / sd(beta_k'x) the
##     standardised indices (Li's (1991) example on the log survival time).
##     The noise level can be changed with CSIRWLS_DGP4_SIGMA.
##
## The constant c is calibrated ONCE PER CELL (calibrate_const_cell) as the
## (1 - target) quantile of log T - shift - e' pooled over 400 independent draws
## of (beta, x); this is exact up to Monte-Carlo error, instead of the noisy
## bisection of R/simulate.R. The realised censoring rate is stored for every
## replication.
## =============================================================================

DGP4_SIGMA <- as.numeric(Sys.getenv("CSIRWLS_DGP4_SIGMA", "0.5"))

spec_ar1    <- function(rho) list(type = "ar1", rho = rho, label = sprintf("AR1(%.1f)", rho))
## neg = number of active coefficients whose sign is flipped, at positions
## redrawn in every replication (Table 8); neg = 0 is the design of Setting (5)
spec_block  <- function(alpha, neg = 0L) list(type = "block", alpha = alpha, neg = as.integer(neg),
                                    label = sprintf("block(%s)%s", paste(alpha, collapse = ","),
                                                    if (neg > 0) sprintf("-neg%d", neg) else ""))
spec_spiked <- function(setting) list(type = "spiked", setting = as.integer(setting),
                                      label = sprintf("spiked%d", setting))
spec_skew   <- function(kind = c("skewnormal", "skewt"), rho = 0.5, delta = 0.8, nu = 5) {
  kind <- match.arg(kind)
  list(type = kind, rho = rho, delta = delta, nu = nu, label = kind)
}

## Table 7 controls built on an Azzalini design 'base' (spec_skew):
##  * elliptical control: Gaussian with exactly the covariance of 'base'
##    (condition C.1 holds, same unbounded leading eigenvalue);
##  * non-elliptical control: x = Sigma_rho^{1/2} w, w with INDEPENDENT
##    standardised components of the same skew law (skew-normal or skew-t_nu),
##    so cov(x) = Sigma_rho (the reference AR(1) design) while C.1 fails.
##    Sigma_rho^{1/2} is taken as the Cholesky factor of the AR(1) matrix,
##    applied by the recursion x_j = rho x_{j-1} + sqrt(1 - rho^2) w_j.
spec_ellctrl    <- function(base) list(type = "ellctrl", base = base, rho = base$rho,
                                       label = paste0("ell-", base$label))
spec_nonellctrl <- function(base) list(type = "nonellctrl", base = base, rho = base$rho,
                                       label = paste0("nonell-", base$label))

## active set of a cell
active_set <- function(spec, model = "dgp1") {
  if (model == "dgp4") return(1:10)
  if (spec$type == "spiked") {
    if (spec$setting == 3L) return(c(1, 10, 15, 20, 25, 30))
    return(c(1, 10, 15, 20, 25, 30, 40, 45, 50))
  }
  1:5
}

## ------------------------------------------------------------------ covariates
## AR(1) Gaussian design: identical to design_ar1() of R/simulate.R (same random
## numbers, X = Z %*% chol(Sigma)), computed by the recursion
## x_j = rho x_{j-1} + sqrt(1 - rho^2) z_j in O(np) instead of O(np^2).
fast_ar1 <- function(n, p, rho) {
  X <- matrix(stats::rnorm(n * p), n, p)
  if (p > 1L) {
    s <- sqrt(1 - rho^2)
    for (j in 2:p) X[, j] <- rho * X[, j - 1L] + s * X[, j]
  }
  X
}

## n x p design. For the spiked designs the standard deviations depend on the
## sample size and the dimension of the study (n_design, p_design); this matters
## only for calibrate_const_cell(), which generates the leading columns only.
gen_X <- function(spec, n, p, n_design = n, p_design = p) {
  switch(spec$type,
    ar1    = fast_ar1(n, p, spec$rho),
    block  = matrix(stats::rnorm(n * p), n, p) %*% block_chol(p, spec$alpha),
    spiked = {
      lam <- spiked_lambda(n_design, p_design, spec$setting)[seq_len(p)]
      matrix(stats::rnorm(n * p), n, p) * rep(lam, each = n)
    },
    skewnormal = , skewt = {
      u  <- fast_ar1(n, p, spec$rho)
      z0 <- abs(stats::rnorm(n))
      x  <- spec$delta * z0 + sqrt(1 - spec$delta^2) * u
      if (spec$type == "skewt") x <- x / sqrt(stats::rchisq(n, spec$nu) / spec$nu)
      x
    },
    ellctrl = matrix(stats::rnorm(n * p), n, p) %*% ctrl_chol(spec, p),
    nonellctrl = {
      w <- skew_components(n, p, spec$base)
      s <- sqrt(1 - spec$rho^2)
      x <- w
      if (p > 1L) for (j in 2:p) x[, j] <- spec$rho * x[, j - 1L] + s * w[, j]
      x
    },
    stop("unknown design ", spec$type))
}

## n x p matrix of independent, standardised (mean 0, variance 1) skew components
## delta |z0| + sqrt(1 - delta^2) u  (skew-normal), divided by sqrt(W / nu) for
## the skew-t; the exact mean and variance are used for the standardisation.
skew_components <- function(n, p, base) {
  d <- base$delta
  x <- d * abs(stats::rnorm(n * p)) + sqrt(1 - d^2) * stats::rnorm(n * p)
  m <- d * sqrt(2 / pi); v <- 1 - m^2
  if (base$type == "skewt") {
    nu <- base$nu; ER <- sqrt(nu / 2) * gamma((nu - 1) / 2) / gamma(nu / 2)
    x <- x / sqrt(stats::rchisq(n * p, nu) / nu)
    m <- m * ER; v <- nu / (nu - 2) - m^2
  }
  matrix((x - m) / sqrt(v), n, p)
}

## Cholesky factor of the covariance of a Gaussian control (cached)
.ctrl_cache <- new.env()
ctrl_chol <- function(spec, p) {
  key <- paste(spec$label, p)
  if (is.null(.ctrl_cache[[key]]))
    .ctrl_cache[[key]] <- chol(sigma_sub(spec, NA, p, seq_len(p)) + 1e-10 * diag(p))
  .ctrl_cache[[key]]
}

## population covariance of the design restricted to the indices idx
sigma_sub <- function(spec, n, p, idx) {
  idx <- as.integer(idx)
  ar1 <- function(rho) rho^abs(outer(idx, idx, "-"))
  switch(spec$type,
    ar1   = ar1(spec$rho),
    block = {
      a <- spec$alpha; inA <- idx <= 5L
      S <- matrix(a[3], length(idx), length(idx))
      S[inA, inA] <- a[1]; S[inA, !inA] <- a[2]; S[!inA, inA] <- a[2]
      diag(S) <- 1; S
    },
    spiked = diag(spiked_lambda(n, p, spec$setting)[idx]^2, length(idx)),
    ellctrl    = sigma_sub(spec$base, n, p, idx),
    nonellctrl = ar1(spec$rho),
    skewnormal = spec$delta^2 * (1 - 2 / pi) + (1 - spec$delta^2) * ar1(spec$rho),
    skewt = {
      nu <- spec$nu; ER <- sqrt(nu / 2) * gamma((nu - 1) / 2) / gamma(nu / 2)
      nu / (nu - 2) * (spec$delta^2 + (1 - spec$delta^2) * ar1(spec$rho)) -
        spec$delta^2 * (2 / pi) * ER^2
    })
}

## population standard deviation of b'x
index_sd <- function(spec, n, p, b) {
  s <- which(b != 0)
  sqrt(as.numeric(crossprod(b[s], sigma_sub(spec, n, p, s) %*% b[s])))
}

## --------------------------------------------------------------- coefficients
## p x K matrix of EDR directions (K = 2 for DGP 5.4)
draw_coef <- function(p, spec, model) {
  if (model == "dgp4") return(cbind(draw_beta(p, 1:5), draw_beta(p, 6:10)))
  b <- draw_beta(p, active_set(spec, model))
  if (!is.null(spec$neg) && spec$neg > 0L) {          # Table 8: sign cancellation
    act <- active_set(spec, model)
    flip <- act[sample.int(length(act), spec$neg)]
    b[flip] <- -b[flip]
  }
  matrix(b, p, 1)
}

## Square root of the block covariance and the two population ranking margins of
## Table 8:  m_half = min_S |(Sigma^{1/2} beta)_j| / max_{S^c} |(Sigma^{1/2} beta)_j|,
## m_one  = same with Sigma beta (the target of the marginal statistics).
## Sigma^{1/2} is the symmetric root, computed once per (p, alpha).
.half_cache <- new.env()
sigma_half_block <- function(p, alpha) {
  key <- paste(p, paste(alpha, collapse = "_"))
  if (is.null(.half_cache[[key]])) {
    S <- block_sigma(p, alpha); e <- eigen(S, symmetric = TRUE)
    .half_cache[[key]] <- list(S = S, H = e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors)))
  }
  .half_cache[[key]]
}
block_margins <- function(spec, p, b, act = 1:5) {
  sh <- sigma_half_block(p, spec$alpha)
  mg <- function(v) { v <- abs(as.vector(v)); min(v[act]) / max(v[-act]) }
  c(m_half = mg(sh$H %*% b), m_one = mg(sh$S %*% b))
}

## ------------------------------------------------------------ survival times
log_time_model <- function(model, X, B, spec, n_design, p_design) {
  n <- nrow(X)
  zstd <- function(k) as.vector(X %*% B[, k]) / index_sd(spec, n_design, p_design, B[, k])
  switch(model,
    dgp1 = , dgp2 = , dgp3 = log_time(as.integer(substr(model, 4, 4)), as.vector(X %*% B[, 1])),
    dgp4 = { z1 <- zstd(1); z2 <- zstd(2); 3 * z1 / (0.5 + (z2 + 1.5)^2) + DGP4_SIGMA * stats::rnorm(n) },
    saturating      = { z <- zstd(1); 3 * tanh(z) + 0.5 * stats::rnorm(n) },
    threshold       = { z <- zstd(1); 2.5 * (z > 0) + 0.5 * stats::rnorm(n) },
    heteroscedastic = { z <- zstd(1); z + 0.5 * exp(z / 2) * stats::rnorm(n) },
    symmetric       = as.vector(X %*% B[, 1])^2 + 0.5 * stats::rnorm(n),
    stop("unknown model ", model))
}

## covariate-dependent part of log C
cens_shift <- function(mech, X, B, spec, n_design, p_design) {
  switch(mech,
    indep = rep(0, nrow(X)),
    edr   = 0.6 * as.vector(X %*% B[, 1]) / index_sd(spec, n_design, p_design, B[, 1]),
    misspec = {
      w <- numeric(p_design); w[41:45] <- 1
      0.8 * rowSums(X[, 41:45, drop = FALSE]) / index_sd(spec, n_design, p_design, w)
    },
    stop("unknown censoring mechanism ", mech))
}

## columns that determine (log T, shift): the calibration only needs those
needed_cols <- function(spec, model, mech) {
  m <- max(active_set(spec, model), if (mech == "misspec") 45 else 0)
  max(m, 6L)          # block_sigma() needs more than the 5 active predictors
}

## ------------------------------------------------------------- calibration
calibrate_const_cell <- function(spec, model, mech, n, p, target, key,
                                 n_beta = 400L, rows = 500L) {
  old <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv))
  set.seed(hash_seed(paste0(key, "|calibration")))
  pc <- min(p, needed_cols(spec, model, mech))
  v <- unlist(lapply(seq_len(n_beta), function(i) {
    Xc <- gen_X(spec, rows, pc, n_design = n, p_design = p)
    Bc <- draw_coef(p, spec, model)[seq_len(pc), , drop = FALSE]
    lt <- log_time_model(model, Xc, Bc, spec, n, p)
    sh <- if (mech == "indep") 0 else {
      Xf <- matrix(0, rows, p); Xf[, seq_len(pc)] <- Xc
      Bf <- matrix(0, p, ncol(Bc)); Bf[seq_len(pc), ] <- Bc
      cens_shift(mech, Xf, Bf, spec, n, p)
    }
    lt - sh - stats::rnorm(rows)
  }))
  as.numeric(stats::quantile(v, 1 - target, names = FALSE))
}

## ----------------------------------------------------------- one data set
## const: calibrated constant of log C. Order of the random draws: design,
## coefficients, survival noise, censoring noise (as in simulate_csir_data()).
gen_data <- function(spec, model, mech, n, p, const) {
  X  <- gen_X(spec, n, p)
  B  <- draw_coef(p, spec, model)
  lt <- log_time_model(model, X, B, spec, n, p)
  lc <- const + cens_shift(mech, X, B, spec, n, p) + stats::rnorm(n)
  list(X = X, y = pmin(lt, lc), delta = as.integer(lt <= lc), B = B,
       active = active_set(spec, model), cens_rate = mean(lt > lc))
}

## warm the Cholesky cache in the parent process before forking
warm_design <- function(spec, p) {
  if (spec$type == "block") { block_chol(p, spec$alpha); sigma_half_block(p, spec$alpha) }
  if (spec$type == "ellctrl") ctrl_chol(spec, p)
  invisible(NULL)
}

## ----------------------------------------------------------- Table 7 layout
## rows of Table 7: label, role, condition C.1, driver and key prefix of the
## stored cells (the reference and the Azzalini designs are the cells of Tables
## 1-2 and S18-S19, shared with run_main.R and run_nonelliptical.R)
table7_rows <- function() {
  sn <- spec_skew("skewnormal"); st <- spec_skew("skewt")
  list(
    list(name = "Gaussian reference design", spec = spec_ar1(0.5), c1 = "holds", driver = "main"),
    list(name = "skew-normal",     spec = sn,                 c1 = "fails", driver = "nonelliptical"),
    list(name = "elliptical control",     spec = spec_ellctrl(sn),    c1 = "holds", driver = "table7"),
    list(name = "non-elliptical control", spec = spec_nonellctrl(sn), c1 = "fails", driver = "table7"),
    list(name = "skew-t5",         spec = st,                 c1 = "fails", driver = "nonelliptical"),
    list(name = "elliptical control",     spec = spec_ellctrl(st),    c1 = "holds", driver = "table7"),
    list(name = "non-elliptical control", spec = spec_nonellctrl(st), c1 = "fails", driver = "table7"))
}
lambda_max <- function(spec, n, p)
  max(eigen(sigma_sub(spec, n, p, seq_len(p)), symmetric = TRUE, only.values = TRUE)$values)
