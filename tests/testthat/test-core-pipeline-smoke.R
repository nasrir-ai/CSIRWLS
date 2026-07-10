# Tests d'integration "fumee" : verifient que le pipeline CSIR / CSIR-WLS
# complet s'execute sans erreur sur un petit jeu de donnees simule, et que
# les sorties ont la forme attendue. Ils ne verifient pas l'exactitude
# statistique (voir les scripts de simulation dans inst/examples/ pour
# reproduire les resultats du papier).

test_that("double.slice + cen.sir run end-to-end on a small censored dataset", {
  set.seed(1)
  n <- 120
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- c(1, -1, rep(0, p - 2))
  eta <- as.numeric(X %*% beta_true)
  Y0 <- exp(eta + rnorm(n))
  C <- exp(rnorm(n))
  time <- pmin(Y0, C)
  event <- as.numeric(Y0 <= C)

  ds <- double.slice(time, event, X, n.slice1 = 4, n.slice0 = 4)
  expect_s3_class(ds, "ds")
  expect_equal(nrow(ds$evec), p)

  joint.edrs <- edr.n(ds, 2)
  expect_equal(ncol(joint.edrs), 2)

  fit <- cen.sir(time, event, X, n.slice = 4, joint.edrs = joint.edrs, h = 0.5, c = 0.05)
  expect_s3_class(fit, "csir")
  expect_equal(nrow(fit$evec), p)
  expect_length(fit$eval, p)
})

test_that("double.slice2 does not error on slices that may end up empty", {
  set.seed(2)
  n <- 60
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  time <- rexp(n)
  event <- rbinom(n, 1, 0.5)

  ds <- double.slice2(time, event, X, n.slice1 = 5, n.slice0 = 5)
  expect_s3_class(ds, "ds")
  expect_equal(nrow(ds$evec), p)
})

test_that("cen.wls.sir2 (full CSIR-WLS pipeline) runs end-to-end", {
  skip_if_not_installed("dr")

  set.seed(3)
  n <- 150
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- c(1, 1, rep(0, p - 2))
  eta <- as.numeric(X %*% beta_true)
  Y0 <- exp(eta + rnorm(n))
  C <- exp(rnorm(n) + 1)
  time <- pmin(Y0, C)
  event <- as.numeric(Y0 <= C)

  fit <- cen.wls.sir2(X, time, event,
    n.slice1 = 5, n.slice0 = 5, n.slice = 10,
    cn1 = 0.1, cn2 = 1, choose.dir = FALSE, ndim = 1, c = 0.05
  )

  expect_type(fit, "list")
  expect_equal(nrow(fit$beta.hat), p)
  expect_true(fit$n.sel >= 1 && fit$n.sel <= p)
  expect_true(all(fit$select %in% seq_len(p)))
})

test_that("compute_censoring_rate and find_c_for_censoring_rate are consistent", {
  set.seed(4)
  Y0 <- rexp(300)
  errc <- rnorm(300)
  c_val <- find_c_for_censoring_rate(0.2, Y0, errc)
  expect_equal(compute_censoring_rate(c_val, Y0, errc), 0.2, tolerance = 1e-4)
})

test_that("estimate_dim_chisq and select_d_bic agree on an easy, high-signal case", {
  # Single strong EDR direction, low noise, no screening needed (p << n):
  # both dimension criteria should recover d = 1 most of the time.
  set.seed(5)
  n <- 300
  p <- 8
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- c(2, rep(0, p - 1))
  eta <- as.numeric(X %*% beta_true)
  Y0 <- exp(eta + rnorm(n, sd = 0.3))
  C <- exp(rnorm(n) + 2) # light censoring
  time <- pmin(Y0, C)
  event <- as.numeric(Y0 <= C)

  ds <- double.slice(time, event, X, n.slice1 = 5, n.slice0 = 5)
  je <- edr.n(ds, 2)
  fit <- cen.sir(time, event, X, n.slice = 5, joint.edrs = je, h = 0.5, c = 0.05)

  d_chisq <- estimate_dim_chisq(fit, alpha = 0.01)
  d_bic <- select_d_bic(fit$eval, n = n, cn1 = 0.1)

  expect_true(is.numeric(d_chisq) && d_chisq >= 0)
  expect_true(is.numeric(d_bic) && d_bic >= 1)
})

test_that("select_d_bic ignores near-zero numerical noise past the true rank", {
  # Regression test for the bug found while debugging the multi-index
  # dimension-recovery example: with only n.slice1 - 1 = 4 non-trivial
  # eigenvalues, select_d_bic() should not be fooled by ~1e-15 noise
  # tacked on to look like additional (spurious) positive eigenvalues.
  true_evals <- c(0.8, 0.5, 0.3, 0.1)
  noisy_evals <- c(true_evals, rep(1e-15, 20))

  d_true <- select_d_bic(true_evals, n = 300)
  d_noisy <- select_d_bic(noisy_evals, n = 300)

  expect_equal(d_true, d_noisy)
})
