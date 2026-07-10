test_that("costeta gives 1 for identical vectors and 0 for orthogonal vectors", {
  v <- c(1, 2, 3, 4)
  expect_equal(as.numeric(costeta(v, v)), 1, tolerance = 1e-10)

  a <- c(1, 0)
  b <- c(0, 1)
  expect_equal(as.numeric(costeta(a, b)), 0, tolerance = 1e-10)
})

test_that("normalize_vector returns a unit-norm vector", {
  v <- c(3, 4)
  nv <- normalize_vector(v)
  expect_equal(sqrt(sum(nv^2)), 1, tolerance = 1e-10)
  expect_equal(nv, c(0.6, 0.8), tolerance = 1e-10)
})

test_that("normalize_vector leaves the zero vector unchanged", {
  v <- c(0, 0, 0)
  expect_equal(normalize_vector(v), v)
})

test_that("Dproj is zero for two collinear vectors and positive for orthogonal ones", {
  a <- c(1, 0, 0)
  b <- c(2, 0, 0) # colinear with a
  cc <- c(0, 1, 0) # orthogonal to a

  expect_equal(Dproj(a, b), 0, tolerance = 1e-10)
  expect_gt(Dproj(a, cc), 0)
})

test_that("cov1 uses a denominator of n (biased) instead of n - 1", {
  set.seed(1)
  x <- matrix(rnorm(50), ncol = 2)
  biased <- cov1(x)
  unbiased <- cov(x)
  n <- nrow(x)
  expect_equal(biased, unbiased * (n - 1) / n, tolerance = 1e-10)
})

test_that("eigen.decomp solves a simple generalized eigenvalue problem", {
  # m2 = identity -> the generalized problem reduces to a standard
  # eigendecomposition of m1.
  m1 <- diag(c(3, 1))
  m2 <- diag(c(1, 1))
  out <- eigen.decomp(m1, m2)
  expect_equal(sort(out$values, decreasing = TRUE), c(3, 1), tolerance = 1e-8)
})
