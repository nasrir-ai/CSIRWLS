## =============================================================================
## app_common.R -- helpers shared by the real-data scripts of Section 6
## =============================================================================
source(file.path("reproduction", "lib", "common.R"))
FIG <- Sys.getenv("CSIRWLS_FIGURES", file.path(REPO, "figures"))
TAB <- Sys.getenv("CSIRWLS_TABLES", file.path(REPO, "tables"))
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)

## time grid used by cv_survival(): 40 points between the 10 % and 90 %
## quantiles of the event times
ibs_grid <- function(y, delta) {
  qq <- stats::quantile(y[delta == 1L], c(0.1, 0.9), names = FALSE)
  seq(qq[1], qq[2], length.out = 40)
}

## in-sample Harrell C-index and Graf IBS of a linear risk index
## (b oriented so that a larger index means a higher risk, as fit_risk_index())
insample_eval <- function(b, X, y, delta) {
  tg <- ibs_grid(y, delta)
  S <- breslow_surv(b, X, y, delta, X, tg)
  c(c_index = c_index(as.vector(X %*% b), y, delta),
    ibs = integrated_brier(S, tg, y, delta, y, delta))
}

## three risk groups cut at the 33 % and 66 % quantiles of the index
risk_groups <- function(eta, probs = c(1/3, 2/3)) {
  q <- stats::quantile(eta, probs, names = FALSE)
  factor(findInterval(eta, q, rightmost.closed = TRUE) + 1L, levels = 1:3,
         labels = c("low", "medium", "high"))
}

## projections versus time: failures as filled dots, censored as "+"
plot_projection <- function(eta, y, delta, main, ylab = "projection") {
  plot(y, eta, type = "n", xlab = "Time", ylab = ylab, main = main)
  points(y[delta == 1L], eta[delta == 1L], pch = 16, cex = 0.8)
  points(y[delta == 0L], eta[delta == 0L], pch = 3, cex = 0.8)
}

## Kaplan-Meier curves of the three risk groups with the log-rank test
plot_km_groups <- function(eta, y, delta, main) {
  g <- risk_groups(eta)
  lr <- logrank_test(y, delta, g)
  cols <- c(low = "#2a78d6", medium = "#2e9e44", high = "#d6402a")
  plot(NA, xlim = c(0, max(y)), ylim = c(0, 1), xlab = "Time", ylab = "Survival probability",
       main = main)
  for (lv in levels(g)) {
    s <- g == lv
    k <- km(y[s], delta[s])
    if (!length(k$times)) { segments(0, 1, max(y[s]), 1, col = cols[[lv]], lwd = 2); next }
    lines(stats::stepfun(k$times, c(1, k$surv)), col = cols[[lv]], do.points = FALSE, lwd = 2,
          xlim = c(0, max(y[s])))
  }
  legend("bottomleft", c("< 33 %", "33-66 %", "> 66 %"), col = cols, lwd = 2, bty = "n")
  mtext(sprintf("log-rank chi2 = %.1f, p = %.1e", lr$statistic, lr$p.value), cex = 0.8)
  invisible(lr)
}

## -----------------------------------------------------------------------------
## Figures 1-2 of Section 6.1. Base R only.
##  * Figure 1: the risk indices are STANDARDISED (mean 0, variance 1) and
##    oriented so that a larger value means a higher risk (fit_risk_index), so
##    the three panels share one scale and can be compared; one colour per
##    method; the number of covariates is in the panel title.
##  * Figure 2: Kaplan-Meier curves of the tertile groups of the index, with the
##    log-rank statistic and the number of covariates in the panel title.
## The time unit of ipred::DLBCL is the month.
## -----------------------------------------------------------------------------
METHOD_COLS <- c(Cox = "#2e9e44", COX = "#2e9e44", CSIR = "#d6402a", `CSIR-WLS` = "#2a78d6")

plot_index_panels <- function(B, X, y, delta, time_lab = "Time (months)") {
  Z <- lapply(B, function(b) as.vector(scale(as.vector(X %*% b))))
  yl <- range(unlist(Z))
  op <- par(mfrow = c(1, length(B)), mar = c(4, 4, 3, 1)); on.exit(par(op))
  for (m in names(B)) {
    col <- METHOD_COLS[[m]]; nv <- sum(B[[m]] != 0)
    plot(y, Z[[m]], type = "n", ylim = yl, xlab = time_lab, ylab = "Standardised risk index",
         main = sprintf("%s (%d covariates)", m, nv), font.main = 1)
    grid(col = "grey90", lty = 1); abline(h = 0, col = "grey70")
    points(y[delta == 1L], Z[[m]][delta == 1L], pch = 16, cex = 1.1, col = col)
    points(y[delta == 0L], Z[[m]][delta == 0L], pch = 3, cex = 1.1, col = col, lwd = 1.5)
    legend("bottomleft", c("death", "censored"), pch = c(16, 3), col = col, bty = "n", cex = 0.9)
  }
}

plot_km_panels <- function(B, X, y, delta, time_lab = "Time (months)") {
  cols <- c(low = "#2a78d6", medium = "#2e9e44", high = "#d6402a")
  op <- par(mfrow = c(1, length(B)), mar = c(4, 4, 3.5, 1)); on.exit(par(op))
  out <- list()
  for (m in names(B)) {
    g <- risk_groups(as.vector(X %*% B[[m]]))
    lr <- logrank_test(y, delta, g); out[[m]] <- lr
    plot(NA, xlim = c(0, max(y)), ylim = c(0, 1), xlab = time_lab, ylab = "Survival probability",
         main = sprintf("%s (%d covariates)\nlog-rank %.1f, p-value %s", m, sum(B[[m]] != 0),
                        lr$statistic, format(signif(lr$p.value, 1), scientific = TRUE)),
         font.main = 1, cex.main = 0.95)
    grid(col = "grey90", lty = 1)
    for (lv in levels(g)) {
      s <- g == lv
      k <- km(y[s], delta[s])
      xs <- c(0, k$times, max(y[s])); ys <- c(1, k$surv, if (length(k$surv)) tail(k$surv, 1) else 1)
      lines(xs, ys, type = "s", col = cols[[lv]], lwd = 2)
      cz <- y[s][delta[s] == 0L]                       # censoring marks
      if (length(cz)) points(cz, stats::approx(xs, ys, cz, method = "constant", f = 0, rule = 2)$y,
                             pch = 3, col = cols[[lv]], cex = 0.8)
    }
    if (m == names(B)[1])                              # one legend, first panel
      legend("bottomright", c("low risk (< 33%)", "medium risk (33% to 66%)", "high risk (> 66%)"),
             col = cols, lwd = 2, bty = "n", cex = 0.85)
  }
  invisible(out)
}

write_csv <- function(d, name) {
  utils::write.csv(d, file.path(TAB, paste0(name, ".csv")), row.names = FALSE)
  message("wrote ", file.path(TAB, paste0(name, ".csv")))
}
