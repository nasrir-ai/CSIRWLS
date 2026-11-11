#' Fit a linear risk index by one of the methods of the paper
#'
#' Fits a single risk index \eqn{\hat\beta^\top x} on training data and orients
#' it so that a larger index means a higher risk, using the training data only.
#' This is the building block of \code{\link{cv_survival}}.
#'
#' Methods:
#' \describe{
#'   \item{\code{"CSIR-WLS"}}{CWLS screening with \eqn{G_c}, then CSIR
#'     (\code{\link{csir_wls}});}
#'   \item{\code{"CSIR"}}{censored SIR on all predictors (requires \eqn{p < n});}
#'   \item{\code{"Cox"}}{Cox fit with a small ridge (requires \eqn{p < n});}
#'   \item{\code{"Cox-lasso"}}{\code{\link{cox_lasso_bic}}; if the path selects
#'     nothing, a ridge Cox fit when \eqn{p < n} and \code{"CSIS+Cox"}
#'     otherwise;}
#'   \item{\code{"CWLS+Cox"}}{CWLS screening with \eqn{G_c}, then a Cox fit on
#'     the retained predictors;}
#'   \item{\code{"CSIS+Cox"}}{the \code{csis_keep} predictors with the largest
#'     CSIS statistics, then a Cox fit.}
#' }
#'
#' @param method one of the names above.
#' @param Xtr,ytr,dtr training data.
#' @param csis_keep number of predictors kept by CSIS (default
#'   \eqn{\lfloor n/\log n\rfloor} on the training data).
#' @param lasso_args list of further arguments for \code{\link{cox_lasso_bic}}.
#' @param n_keep number of predictors retained by CWLS in \code{"CSIR-WLS"}
#'   and \code{"CWLS+Cox"} (\code{NULL}: the stopping rule \eqn{G_c}).
#' @return The vector of the \code{p} coefficients of the index.
#' @export
fit_risk_index <- function(method, Xtr, ytr, dtr, csis_keep = NULL, lasso_args = list(),
                           n_keep = NULL) {
  n <- nrow(Xtr); p <- ncol(Xtr)
  if (is.null(csis_keep)) csis_keep <- as.integer(n / log(n))
  screen_cox <- function(A) {
    b <- numeric(p)
    b[A] <- cox_fit(Xtr[, A, drop = FALSE], ytr, dtr, l2 = 1e-2, n_iter = 150)
    b
  }
  b <- switch(method,
    `CSIR-WLS`  = as.vector(csir_wls(Xtr, ytr, dtr, d = 1L, n_keep = n_keep)$beta),
    CSIR        = as.vector(csir(Xtr, ytr, dtr, d = 1L)),
    Cox         = cox_fit(Xtr, ytr, dtr, l2 = 1e-3, n_iter = 300),
    `Cox-lasso` = {
      bb <- do.call(cox_lasso_bic, c(list(Xtr, ytr, dtr), lasso_args))
      if (any(bb != 0)) bb
      else if (p < n) cox_fit(Xtr, ytr, dtr, l2 = 1e-2, n_iter = 200)
      else screen_cox(order(csis(Xtr, ytr, dtr), decreasing = TRUE)[seq_len(csis_keep)])
    },
    `CWLS+Cox`  = {
      om <- cwls(Xtr, ytr, dtr)
      screen_cox(order(om, decreasing = TRUE)[seq_len(if (is.null(n_keep)) bic_cut(om, n, p) else n_keep)])
    },
    `CSIS+Cox`  = screen_cox(order(csis(Xtr, ytr, dtr), decreasing = TRUE)[seq_len(csis_keep)]),
    stop("unknown method '", method, "'", call. = FALSE))
  if (c_index(as.vector(Xtr %*% b), ytr, dtr) < 0.5) b <- -b
  b
}

#' Repeated K-fold cross-validation of survival risk indices
#'
#' Cross-validated Harrell C-index and integrated Brier score of the risk
#' indices of \code{\link{fit_risk_index}}. Every step that uses the response,
#' including the screening step of CSIR-WLS, is repeated inside every training
#' fold, so the cross-validated values are not affected by the optimism of an
#' in-sample evaluation after screening. Optionally the predictors are
#' standardised with the mean and standard deviation of the training fold.
#'
#' For each repetition the fold predictions are pooled and the C-index and the
#' integrated Brier score (on a grid of 40 times between the 10\% and 90\%
#' quantiles of the event times) are computed on the whole sample.
#'
#' @param X numeric matrix (n x p).
#' @param y observed times.
#' @param delta event indicator.
#' @param methods methods, see \code{\link{fit_risk_index}}.
#' @param K number of folds.
#' @param reps number of random fold assignments (a vector of length
#'   \code{length(methods)} is allowed, for example to repeat a slow method
#'   fewer times).
#' @param standardise standardise the predictors within each training fold.
#' @param seeds seeds of the fold assignments (default \code{seq_len(reps)}).
#' @param insample also return the in-sample values.
#' @param ... further arguments for \code{\link{fit_risk_index}}.
#' @return A data frame with one row per method: the mean of the
#'   cross-validated C-index and IBS over the repetitions, their standard
#'   deviations (\code{c_sd}, \code{ibs_sd}) and standard errors (\code{c_se},
#'   \code{ibs_se}), the average number of non-zero coefficients
#'   (\code{n_selected}) and, if requested, the in-sample values. The values of
#'   every repetition are in the attribute \code{"reps"}.
#' @examples
#' \donttest{
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   pbc <- survival::pbc
#'   pbc <- pbc[stats::complete.cases(pbc), ]
#'   X <- scale(as.matrix(pbc[, c("age", "edema", "bili", "albumin",
#'                                "platelet", "protime")]))
#'   cv_survival(X, pbc$time, as.integer(pbc$status == 2),
#'               methods = c("Cox", "CSIR-WLS"), reps = 2)
#' }
#' }
#' @export
cv_survival <- function(X, y, delta, methods = c("Cox", "Cox-lasso", "CSIR", "CSIR-WLS"),
                        K = 10L, reps = 20L, standardise = FALSE, seeds = NULL,
                        insample = TRUE, ...) {
  X <- check_xyd(X, y, delta)
  n <- nrow(X)
  reps <- rep_len(reps, length(methods))
  if (is.null(seeds)) seeds <- seq_len(max(reps))
  qq <- stats::quantile(y[delta == 1L], c(0.1, 0.9), names = FALSE)
  tg <- seq(qq[1], qq[2], length.out = 40)
  stdz <- function(Xtr, Xte) {
    if (!standardise) return(list(tr = Xtr, te = Xte))
    mu <- colMeans(Xtr); sd <- apply(Xtr, 2L, stats::sd); sd[!(sd > 0)] <- 1
    list(tr = sweep(sweep(Xtr, 2L, mu), 2L, sd, "/"), te = sweep(sweep(Xte, 2L, mu), 2L, sd, "/"))
  }
  per_rep <- list(); rows <- list()
  for (mi in seq_along(methods)) {
    m <- methods[mi]
    cvals <- ibs <- nsel <- numeric(0)
    for (r in seq_len(reps[mi])) {
      set.seed(seeds[r])
      idx <- sample(n)
      folds <- split(idx, cut(seq_along(idx), K, labels = FALSE))
      risk <- numeric(n); S <- matrix(0.5, n, length(tg)); ns <- numeric(0)
      for (f in folds) {
        tr <- setdiff(idx, f)
        ok <- try({
          sx <- stdz(X[tr, , drop = FALSE], X[f, , drop = FALSE])
          b <- fit_risk_index(m, sx$tr, y[tr], delta[tr], ...)
          risk[f] <- as.vector(sx$te %*% b)
          S[f, ] <- breslow_surv(b, sx$tr, y[tr], delta[tr], sx$te, tg)
          ns <- c(ns, sum(abs(b) > 1e-10))
        }, silent = TRUE)
        if (inherits(ok, "try-error")) { risk[f] <- 0; S[f, ] <- 0.5 }
      }
      cvals <- c(cvals, c_index(risk, y, delta))
      ibs <- c(ibs, integrated_brier(S, tg, y, delta, y, delta))
      nsel <- c(nsel, if (length(ns)) mean(ns) else NA_real_)
    }
    per_rep[[m]] <- data.frame(rep = seq_along(cvals), c_index = cvals, ibs = ibs, n_selected = nsel)
    sdv <- function(v) if (length(v) > 1L) stats::sd(v) else NA_real_
    row <- data.frame(method = m, reps = length(cvals),
                      c_index = mean(cvals), c_sd = sdv(cvals), c_se = sdv(cvals) / sqrt(length(cvals)),
                      ibs = mean(ibs), ibs_sd = sdv(ibs), ibs_se = sdv(ibs) / sqrt(length(ibs)),
                      n_selected = mean(nsel, na.rm = TRUE), stringsAsFactors = FALSE)
    if (insample) {
      sx <- stdz(X, X)
      b <- fit_risk_index(m, sx$tr, y, delta, ...)
      Sin <- breslow_surv(b, sx$tr, y, delta, sx$tr, tg)
      row$c_insample <- c_index(as.vector(sx$tr %*% b), y, delta)
      row$ibs_insample <- integrated_brier(Sin, tg, y, delta, y, delta)
    }
    rows[[m]] <- row
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL
  attr(out, "reps") <- per_rep
  out
}
