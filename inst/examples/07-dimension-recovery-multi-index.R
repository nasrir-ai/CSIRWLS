## =============================================================================
## Compact simulation study: recovery of the TRUE structural dimension
## d = 2 in a MULTI-INDEX censored model, comparing:
##   (a) the classical sequential chi-squared test (Li 1991), via
##       estimate_dim_chisq() (wraps sir.test())
##   (b) the BIC-type criterion of Section 5.3.2, via select_d_bic()
##
## Meant to substantiate the "multiple index model" claim in Section 5.3.2
## with simulation evidence (not just the PBC real-data example), and to
## give empirical support for a regime-dependent recommendation:
##   - chi-squared test  : better suited for p << n
##   - BIC-type criterion: better suited for p >> n
##
## Kept intentionally small (2 values of p, 2 censoring rates, N = 100
## replications each) so it adds roughly half a page, not a new section.
##
## History / bug fixes (worked out interactively before landing here):
##  1. The screening step originally called cen.wls() with arguments it
##     does not have (n.slice1, n.slice0, ndim) and read a "$w" field
##     that cen.wls() never returns -- every p >> n replication silently
##     errored out (100% failure). Fixed: call cen.wls() with its real
##     signature and screen on cwls_scr$select directly.
##  2. The BIC criterion was re-implemented locally ("estimate_K()")
##     with eigvals[eigvals > 0] (no numerical tolerance). With only
##     n.slice1 slices, sigma.eta has rank <= n.slice1 - 1; beyond that,
##     "eigenvalues" are floating-point noise (~1e-15), not exact zeros,
##     and passed that filter -- inflating the candidate space and
##     making d_hat collapse to (almost) the same value regardless of
##     the data. Fixed: use select_d_bic() (R/04-variable-selection-wls.R),
##     which filters with a numerical tolerance (> 1e-10).
##  3. The kernel bandwidth was fixed at h = 0.2; switched to a
##     data-driven Silverman bandwidth (as in cen.wls.sir2()) for better
##     adaptivity to the projected data's scale.
##  4. Added diagnostics: table(d_chisq) / table(d_bic) printed for each
##     configuration, since the aggregate recovery rate alone hid the
##     fact that estimate_K() was stuck on a single value.
##
## Once you've confirmed with the diagnostics that BOTH criteria behave
## as expected under an "easy" DGP (see Section 6 at the bottom), this
## script's main configs are the ones worth reporting in the paper.
## =============================================================================

library(CSIRWLS)
library(MASS)

set.seed(2024)

## ---------------------------------------------------------------------
## 1) Multi-index DGP with a KNOWN true structural dimension d = 2
##    Two distinct active sets S1, S2 (disjoint) define two genuinely
##    different EDR directions beta1, beta2. The ratio-form test
##    function ensures BOTH directions matter for y0 -- but note this
##    also makes beta2's signal weak/non-monotone from a SIR viewpoint
##    (see Section 6 below for a sanity check against an "easy" additive
##    DGP where both directions are trivially detectable).
## ---------------------------------------------------------------------
make_betas <- function(p, S1 = 1:5, S2 = 6:10) {
  beta1 <- rep(0, p)
  beta1[S1] <- runif(length(S1), 0.5, 1.5)
  beta2 <- rep(0, p)
  beta2[S2] <- runif(length(S2), 0.5, 1.5)
  list(beta1 = beta1, beta2 = beta2)
}

generate_data <- function(n, p, rho, cens_c, beta1, beta2, sigma = 1) {
  Sigma <- rho^abs(outer(1:p, 1:p, "-"))
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  lin1 <- as.vector(X %*% beta1)
  lin2 <- as.vector(X %*% beta2)
  eps <- rnorm(n)
  ## classical multi-index test function (ratio form, Li 1991-style)
  log_y0 <- lin1 / (0.5 + (lin2 + 1.5)^2) + sigma * eps
  y0 <- exp(log_y0)
  Cc <- exp(cens_c + rnorm(n))
  time <- pmin(y0, Cc)
  status <- as.numeric(y0 <= Cc)
  list(X = X, time = time, status = status, cens_rate = mean(1 - status))
}

## ---------------------------------------------------------------------
## 2) Censoring-rate calibration (small pilot search), same principle
##    used to hit the 15% / 40% targets elsewhere in the paper.
## ---------------------------------------------------------------------
calibrate_c <- function(target_rate, n, p, rho, beta1, beta2,
                         M = 20, c_grid = seq(-2, 4, by = 0.25)) {
  rates <- sapply(c_grid, function(c_val) {
    mean(replicate(M, generate_data(n, p, rho, c_val, beta1, beta2)$cens_rate))
  })
  c_grid[which.min(abs(rates - target_rate))]
}

## ---------------------------------------------------------------------
## 3) One replication: estimate d via both criteria
##    For p >> n, CWLS pre-screening is applied first, exactly as in the
##    real TRANSBIG pipeline, before computing the eigenvalues that both
##    estimate_dim_chisq() and select_d_bic() rely on.
## ---------------------------------------------------------------------
run_one_rep <- function(n, p, rho, cens_c, beta1, beta2,
                         n.slice1 = 5, n.slice0 = 5, c_floor = 0.05,
                         screen_first = FALSE, cn1 = 0.1, cn2 = 1, n_sel = 20) {
  dat <- generate_data(n, p, rho, cens_c, beta1, beta2)
  X <- dat$X
  time <- dat$time
  status <- dat$status

  tryCatch(
    {
      if (screen_first) {
        cwls_scr <- cen.wls(X, time, status,
          c = c_floor, n.slice = n.slice1 + n.slice0,
          cn1 = cn1, cn2 = cn2, choose.dir = FALSE
        )
        sel <- cwls_scr$select
        if (length(sel) > n_sel) sel <- sel[1:n_sel]
        Xs <- X[, sel, drop = FALSE]
      } else {
        Xs <- X
      }

      beta_init <- prcomp(Xs, rank. = 1)$rotation[, 1]
      h <- bw.nrd(as.vector(Xs %*% beta_init))

      ds <- double.slice(time, status, Xs, n.slice1, n.slice0)
      je <- edr.n(ds, 2)
      fit <- cen.sir(time, status, Xs, n.slice = n.slice1, joint.edrs = je, h = h, c = c_floor)

      c(
        d_chisq = estimate_dim_chisq(fit, alpha = 0.01),
        d_bic   = select_d_bic(fit$eval, n = n, cn1 = cn1)
      )
    },
    error = function(e) {
      cat("    [erreur capturee] ", conditionMessage(e), "\n")
      c(d_chisq = NA_real_, d_bic = NA_real_)
    }
  )
}

## ---------------------------------------------------------------------
## 4) Simulation driver : N replications -> recovery rate P(d_hat = 2),
##    plus the full distribution of d_hat for diagnostics.
## ---------------------------------------------------------------------
run_simulation <- function(n = 500, p, rho = 0.5, cens_c, N = 100, screen_first = FALSE) {
  betas <- make_betas(p)
  res <- t(sapply(1:N, function(i) {
    run_one_rep(n, p, rho, cens_c, betas$beta1, betas$beta2, screen_first = screen_first)
  }))
  res <- as.data.frame(res)

  cat("  -> distribution d_chisq :\n")
  print(table(res$d_chisq, useNA = "ifany"))
  cat("  -> distribution d_bic   :\n")
  print(table(res$d_bic, useNA = "ifany"))

  data.frame(
    p = p,
    censoring = NA, # filled by caller
    recovery_chisq = mean(res$d_chisq == 2, na.rm = TRUE),
    recovery_bic = mean(res$d_bic == 2, na.rm = TRUE),
    n_failed = sum(is.na(res$d_chisq))
  )
}

## ---------------------------------------------------------------------
## 5) Run : p = 80 (p << n, low-dim regime) and p = 1500 (p >> n),
##    each at two censoring rates (15%, 40%), N = 100 replications.
## ---------------------------------------------------------------------
n <- 500
rho <- 0.5
N <- 100
configs <- expand.grid(p = c(80, 1500), target_cens = c(0.15, 0.40))

results <- do.call(rbind, lapply(seq_len(nrow(configs)), function(i) {
  p_i <- configs$p[i]
  target_i <- configs$target_cens[i]
  betas_i <- make_betas(p_i)
  c_val <- calibrate_c(target_i, n, p_i, rho, betas_i$beta1, betas_i$beta2)
  cat(sprintf(
    "Running p = %d, target censoring = %.0f%% (c = %.2f) ...\n",
    p_i, target_i * 100, c_val
  ))
  out <- run_simulation(
    n = n, p = p_i, rho = rho, cens_c = c_val, N = N,
    screen_first = (p_i > n)
  )
  out$censoring <- target_i
  out
}))

results$censoring <- paste0(results$censoring * 100, "%")
print(results)

## ---------------------------------------------------------------------
## 6) Sanity check : an "easy" additive DGP (no ratio form)
##    If recovery is still poor here, the issue is in the pipeline, not
##    in the intrinsic difficulty of the ratio-form DGP above. Run this
##    FIRST when debugging a new change to run_one_rep()/select_d_bic().
## ---------------------------------------------------------------------
generate_data_easy <- function(n, p, rho, cens_c, beta1, beta2, sigma = 1) {
  Sigma <- rho^abs(outer(1:p, 1:p, "-"))
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  lin1 <- as.vector(X %*% beta1)
  lin2 <- as.vector(X %*% beta2)
  eps <- rnorm(n)
  log_y0 <- lin1 + lin2 + sigma * eps # additive, both directions "easy" for SIR
  y0 <- exp(log_y0)
  Cc <- exp(cens_c + rnorm(n))
  time <- pmin(y0, Cc)
  status <- as.numeric(y0 <= Cc)
  list(X = X, time = time, status = status, cens_rate = mean(1 - status))
}

## Decommenter pour lancer le sanity check (p = 80, censure 15%, N = 50) :
# betas_easy <- make_betas(80)
# c_easy <- calibrate_c(0.15, n = 500, p = 80, rho = 0.5, betas_easy$beta1, betas_easy$beta2)
# res_easy <- t(sapply(1:50, function(i) {
#   dat <- generate_data_easy(500, 80, 0.5, c_easy, betas_easy$beta1, betas_easy$beta2)
#   beta_init <- prcomp(dat$X, rank. = 1)$rotation[, 1]
#   h <- bw.nrd(as.vector(dat$X %*% beta_init))
#   ds <- double.slice(dat$time, dat$status, dat$X, 5, 5)
#   je <- edr.n(ds, 2)
#   fit <- cen.sir(dat$time, dat$status, dat$X, n.slice = 5, joint.edrs = je, h = h, c = 0.05)
#   c(d_chisq = estimate_dim_chisq(fit, alpha = 0.01),
#     d_bic   = select_d_bic(fit$eval, n = 500, cn1 = 0.1))
# }))
# print(table(res_easy[, "d_chisq"]))
# print(table(res_easy[, "d_bic"]))
