test_that("X_rand returns a matrix of the right dimensions (autoR)", {
  set.seed(1)
  X <- X_rand(n = 50, p = 5, q = 2, type = "autoR", rho = 0.5)
  expect_equal(dim(X), c(50, 5))
})

test_that("X_rand returns a matrix of the right dimensions (block)", {
  set.seed(1)
  X <- X_rand(n = 50, p = 10, q = 3, type = "block", alpha = c(0.5, 0.2, 0.1))
  expect_equal(dim(X), c(50, 10))
})

test_that("generate_X supports the normal and autoR types", {
  set.seed(1)
  Xn <- generate_X(n = 40, p = 5, type = "normal")
  expect_equal(dim(Xn), c(40, 5))

  Xr <- generate_X(n = 40, p = 5, type = "autoR", rho = 0.3)
  expect_equal(dim(Xr), c(40, 5))
})

test_that("generate_X errors clearly on an unsupported type", {
  expect_error(generate_X(n = 10, p = 2, type = "not_a_real_type"), "Unsupported X_type")
})

test_that("generate_survival_times works with a vector Beta (single index)", {
  set.seed(1)
  X <- matrix(rnorm(100 * 3), 100, 3)
  Beta <- c(1, -1, 0)
  y <- generate_survival_times(X, Beta, model_type = "PHM")
  expect_length(y, 100)
  expect_true(all(y > 0))
})

test_that("generate_survival_times works with a matrix Beta (multi index)", {
  set.seed(1)
  X <- matrix(rnorm(100 * 4), 100, 4)
  Beta <- matrix(0, 4, 2)
  Beta[1, 1] <- 1
  Beta[2, 2] <- 1
  y <- generate_survival_times(X, Beta, model_type = "PHM2")
  expect_length(y, 100)
  expect_true(all(y > 0))
})

test_that("generate_survival_times_multi is an alias of generate_survival_times", {
  set.seed(42)
  X <- matrix(rnorm(60 * 3), 60, 3)
  Beta <- c(1, 0, -1)

  set.seed(1)
  y1 <- generate_survival_times(X, Beta, "PHM")
  set.seed(1)
  y2 <- generate_survival_times_multi(X, Beta, "PHM")
  expect_equal(y1, y2)
})

test_that("find_c_for_censoring_rate calibrates the target censoring rate", {
  set.seed(1)
  Y0 <- rexp(500, rate = 1)
  errc <- rnorm(500)
  target <- 0.3
  c_val <- find_c_for_censoring_rate(target, Y0, errc)
  achieved_rate <- compute_censoring_rate(c_val, Y0, errc)
  expect_equal(achieved_rate, target, tolerance = 1e-4)
})
