## =============================================================================
## Section 6.2 -- PBC: dimension test, Table 10, Table 11, Figure 3
##
## Data: survival::pbc, complete cases (n = 276, 111 deaths); death = event,
## transplant and alive = censored. Six covariates (age, edema, bili, albumin,
## platelet, protime), standardised to zero mean and unit variance.
##  * sequential chi-squared test of Li (1991) on the double-slicing CSIR
##    eigenvalues, alpha = 0.01 (csir_dim_test)            -> K_hat
##  * Table 10: CSIR and CSIR-WLS with d = K_hat (= 2), B = 500 bootstrap
##    resamples. Each bootstrap direction is sign-aligned with the full-sample
##    direction. Reported: full-sample estimate, bootstrap SD, percentile 95 %
##    interval (the bootstrap mean is in the CSV too).
##    Sign convention of the full-sample directions: a larger projection means
##    a longer survival (which gives a negative bilirubin loading on Dir 1,
##    as in the paper).
##  * Table 11: in-sample and 10-fold CV C-index and IBS, 25 random partitions,
##    screening repeated inside every training fold (cv_survival).
##  * Figure 3: the two estimated directions versus time.
## =============================================================================
source(file.path("reproduction", "applications", "app_common.R"))
B_BOOT <- n_rep(500); R_CV <- n_rep(25)
pbc <- survival::pbc
pbc <- pbc[stats::complete.cases(pbc), ]
vars <- c("age", "edema", "bili", "albumin", "platelet", "protime")
X <- scale(as.matrix(pbc[, vars])); dimnames(X) <- list(NULL, vars)
y <- pbc$time; delta <- as.integer(pbc$status == 2)
cat(sprintf("PBC: n = %d, deaths = %d (%.1f %%)\n", nrow(X), sum(delta), 100 * mean(delta)))

## dimension of the EDR space
dt <- csir_dim_test(X, y, delta, H = 5L, alpha = 0.01)
print(dt$table); cat("K_hat =", dt$K, "\n")
write_csv(dt$table, "pbc_dimension_test")
K <- max(1L, dt$K)

orient <- function(B) {          # larger projection = longer survival
  for (k in seq_len(ncol(B))) if (c_index(-as.vector(X %*% B[, k]), y, delta) < 0.5) B[, k] <- -B[, k]
  B
}
fit_both <- function(Xb, yb, db) list(
  CSIR = csir(Xb, yb, db, H = 5L, d = K),
  `CSIR-WLS` = csir_wls(Xb, yb, db, H = 5L, d = K)$beta)
full <- lapply(fit_both(X, y, delta), orient)

set.seed(2026)
boot <- lapply(seq_len(B_BOOT), function(b) {
  i <- sample(nrow(X), replace = TRUE)
  fb <- fit_both(X[i, ], y[i], delta[i])
  lapply(names(fb), function(m) {           # sign alignment, direction by direction
    Bb <- fb[[m]]
    for (k in seq_len(K)) if (sum(Bb[, k] * full[[m]][, k]) < 0) Bb[, k] <- -Bb[, k]
    Bb
  })
})
tab10 <- do.call(rbind, lapply(seq_along(full), function(mi) {
  m <- names(full)[mi]
  do.call(rbind, lapply(seq_len(K), function(k) {
    draws <- sapply(boot, function(z) z[[mi]][, k])      # p x B
    data.frame(method = m, direction = paste("Dir", k), variable = vars,
               value = round(full[[m]][, k], 2), boot_mean = round(rowMeans(draws), 2),
               sd = round(apply(draws, 1, stats::sd), 2),
               lower = round(apply(draws, 1, stats::quantile, 0.025), 2),
               upper = round(apply(draws, 1, stats::quantile, 0.975), 2))
  }))
}))
tab10$cell <- sprintf("%.2f (%.2f) [%.2f, %.2f]", tab10$value, tab10$sd, tab10$lower, tab10$upper)
for (m in names(full)) {
  w <- tab10[tab10$method == m, ]
  cat("\n", m, "\n")
  print(data.frame(variable = vars, `Dir 1` = w$cell[w$direction == "Dir 1"],
                   `Dir 2` = if (K > 1) w$cell[w$direction == "Dir 2"] else "", check.names = FALSE))
}
## With p = 6, G_c keeps the six covariates on the full sample, so the
## full-sample CSIR-WLS directions coincide with the CSIR ones; the two
## bootstrap distributions differ only through the resamples in which G_c
## drops a covariate.
cat("covariates retained by G_c on the full sample:", length(csir_wls(X, y, delta, d = K)$A), "\n")
write_csv(tab10, "table10_pbc_bootstrap")

## Table 11
cv <- cv_survival(X, y, delta, methods = c("Cox", "Cox-lasso", "CSIR", "CSIR-WLS"),
                  K = 10L, reps = R_CV, insample = TRUE)
tab11 <- data.frame(Method = c(Cox = "COX", `Cox-lasso` = "Cox-lasso", CSIR = "CSIR", `CSIR-WLS` = "CSIR-WLS")[cv$method],
                    `C in sample` = round(cv$c_insample, 3), `IBS in sample` = round(cv$ibs_insample, 3),
                    `C 10-fold CV` = sprintf("%.3f (%.3f)", cv$c_index, cv$c_sd),
                    `IBS 10-fold CV` = sprintf("%.3f (%.3f)", cv$ibs, cv$ibs_sd), check.names = FALSE)
print(tab11); write_csv(tab11, "table11_pbc")

## Figure 3: standardised projections on the two directions (CSIR in red, CSIR-WLS in blue)
plot_fig3 <- function() {
  Z <- lapply(full, function(B) apply(X %*% B, 2, function(v) as.vector(scale(v))))
  yl <- range(unlist(Z))
  op <- par(mfrow = c(2, K), mar = c(4, 4, 3, 1)); on.exit(par(op))
  for (m in names(full)) for (k in seq_len(K)) {
    col <- METHOD_COLS[[m]]; z <- Z[[m]][, k]
    plot(y, z, type = "n", ylim = yl, xlab = "Time (days)", ylab = "Standardised projection",
         main = sprintf("%s, direction %d", m, k), font.main = 1)
    grid(col = "grey90", lty = 1); abline(h = 0, col = "grey70")
    points(y[delta == 1L], z[delta == 1L], pch = 16, cex = 0.9, col = col)
    points(y[delta == 0L], z[delta == 0L], pch = 3, cex = 0.9, col = col, lwd = 1.2)
    legend("bottomright", c("death", "censored"), pch = c(16, 3), col = col, bty = "n", cex = 0.85)
  }
}
grDevices::pdf(file.path(FIG, "figure3_pbc_projections.pdf"), width = 9, height = 7); plot_fig3(); grDevices::dev.off()
grDevices::png(file.path(FIG, "figure3_pbc_projections.png"), width = 9, height = 7, units = "in", res = 200); plot_fig3(); grDevices::dev.off()
