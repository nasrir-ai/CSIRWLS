## =============================================================================
## Section 6.3 -- optimism of the in-sample C-index after screening, on data of
## the TRANSBIG shape without any signal: n = 198, p = 2000 independent N(0,1)
## covariates, survival time independent of x, 69 % censoring, 40 replications.
## Both pipelines (CSIS -> Cox keeping floor(n / log n) = 37 covariates, and
## CWLS -> CSIR with G_c) are applied
##   (a) to the full sample and evaluated on it (in-sample C-index);
##   (b) to a random half and evaluated on the other half (held-out C-index).
## The correct value is 0.5. Paper: in sample 0.902 (0.002) and 0.718 (0.013),
## held out 0.508 and 0.505.
## =============================================================================
source(file.path("reproduction", "applications", "app_common.R"))
n <- 198L; p <- 2000L; R <- n_rep(40)
const <- sqrt(2) * stats::qnorm(1 - 0.69)        # P(log T > log C) = 0.69
one <- function(r) {
  set.seed(rep_seed("null-optimism", r))
  X <- matrix(stats::rnorm(n * p), n, p)
  lt <- stats::rnorm(n); lc <- const + stats::rnorm(n)
  y <- pmin(lt, lc); d <- as.integer(lt <= lc)
  fit <- function(m, i) fit_risk_index(m, X[i, ], y[i], d[i])
  te <- sample(n, n %/% 2); tr <- setdiff(seq_len(n), te)
  out <- c(cens = mean(d == 0))
  for (m in c("CSIS+Cox", "CSIR-WLS")) {
    b <- fit(m, seq_len(n))
    out[paste0("in_", m)] <- c_index(as.vector(X %*% b), y, d)
    b <- fit(m, tr)
    out[paste0("out_", m)] <- c_index(as.vector(X[te, ] %*% b), y[te], d[te])
  }
  out
}
res <- do.call(rbind, if (CFG$cores > 1) parallel::mclapply(seq_len(R), one, mc.cores = CFG$cores)
                      else lapply(seq_len(R), one))
summ <- data.frame(pipeline = c("CSIS -> Cox", "CWLS -> CSIR"),
  in_sample = sprintf("%.3f (%.3f)", colMeans(res[, c("in_CSIS+Cox", "in_CSIR-WLS")]),
                      apply(res[, c("in_CSIS+Cox", "in_CSIR-WLS")], 2, stats::sd) / sqrt(R)),
  held_out = sprintf("%.3f (%.3f)", colMeans(res[, c("out_CSIS+Cox", "out_CSIR-WLS")]),
                     apply(res[, c("out_CSIS+Cox", "out_CSIR-WLS")], 2, stats::sd) / sqrt(R)))
cat(sprintf("replications: %d, mean censoring: %.3f\n", R, mean(res[, "cens"])))
print(summ); write_csv(summ, "transbig_null_optimism")
