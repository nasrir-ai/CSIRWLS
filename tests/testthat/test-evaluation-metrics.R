test_that("trace_corr is 1 for identical subspaces and near 0 for orthogonal ones", {
  B1 <- matrix(c(1, 0, 0, 1, 0, 0), ncol = 2)
  expect_equal(trace_corr(B1, B1), 1, tolerance = 1e-10)

  B2 <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2)
  expect_equal(trace_corr(B1, B2), 0, tolerance = 1e-10)
})

test_that("subspace_dist is 0 for identical subspaces", {
  B1 <- matrix(c(1, 0, 0, 1, 0, 0), ncol = 2)
  expect_equal(subspace_dist(B1, B1), 0, tolerance = 1e-10)
})

test_that("safe_cor returns 0 instead of NA for degenerate input", {
  a <- c(1, 1, 1, 1)
  b <- c(1, 2, 3, 4)
  expect_equal(safe_cor(a, b), 0) # cor() with a constant vector is NA
})

test_that("compute_TPR_FPR recovers perfect detection", {
  out <- compute_TPR_FPR(estimated_active = 1:5, true_active = 1:5, total_features = 20)
  expect_equal(out$TPR, 1)
  expect_equal(out$FPR, 0)
  expect_equal(out$FDR, 0)
})

test_that("compute_metrics matches compute_TPR_FPR on the same inputs", {
  est <- c(1, 2, 3, 10)
  true <- 1:5
  p <- 20

  full <- compute_TPR_FPR(est, true, p)
  compact <- compute_metrics(est, true, p)

  expect_equal(compact$TPR, full$TPR)
  expect_equal(compact$FPR, full$FPR)
  expect_equal(compact$FDR, full$FDR)
})
