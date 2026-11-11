## =============================================================================
## methods.R -- the estimators of Section 5 wrapped for the simulation drivers
##
## Conventions (Sections 5.2-5.3 of the paper):
##   * H = 5 slices for every SIR-type method; r = 2 for the kernel weights;
##   * CWLS keeps the G_c-minimiser (bic_cut) with at least K + 1 predictors,
##     exactly as csir_wls(); the marginal statistics keep floor(n / log n);
##   * SIR, CSIR, CP-SIR and Cox are run only when p < n;
##   * CP-SIR receives the log time (time_is_log = TRUE): its window is
##     Silverman's rule on the log failure times;
##   * Cox is cox_fit() with its defaults (l2 = 1e-4), Cox-lasso is
##     cox_lasso_bic() (BIC over a path of 12 penalties).
## =============================================================================

H_SLICES <- 5L

## ------------------------------------------------------------------ screening
## CWLS + stopping rule G_c; returns the scores, the G_c size and the set A
cwls_select <- function(X, y, delta, K = 1L, floor_mult = NULL) {
  n <- nrow(X); p <- ncol(X)
  om <- cwls(X, y, delta, H = H_SLICES, r = 2L, floor_mult = floor_mult)
  k  <- bic_cut(om, n, p)
  A  <- order(om, decreasing = TRUE)[seq_len(min(p, max(k, K + 1L)))]
  list(omega = om, k = k, A = A)
}

## TPR / FPR / FDR / size of a selected set
sel_rates <- function(sel, active, p) screening_rates(sel, active, p)

## top-'keep' predictors of a score vector
top_k <- function(scores, keep) order(scores, decreasing = TRUE)[seq_len(keep)]

## minimum model size of a ranking given by |b|, ties (the zero coefficients of
## a lasso fit) broken at random -- Section 5.4
mms_random_ties <- function(b, active) {
  rk <- order(-abs(b), stats::runif(length(b)))
  max(match(active, rk))
}

## ---------------------------------------------------------- second stage (EDR)
## p x K matrix of CSIR directions computed on the predictors in A
csir_on <- function(X, y, delta, A, K = 1L) {
  B <- matrix(0, ncol(X), K)
  B[A, ] <- csir(X[, A, drop = FALSE], y, delta, H = H_SLICES, d = K)
  B
}
cpsir_on <- function(X, y, delta, A, K = 1L) {
  B <- matrix(0, ncol(X), K)
  B[A, ] <- cp_sir(X[, A, drop = FALSE], y, delta, d = K, time_is_log = TRUE)
  B
}

## one EDR estimator; returns a p x K matrix (p x 1 for Cox and Cox-lasso)
edr_fit <- function(method, X, y, delta, K = 1L) {
  switch(method,
    SIR        = sir(X, y, H = H_SLICES, d = K),
    CSIR       = csir(X, y, delta, H = H_SLICES, d = K),
    `CP-SIR`   = cp_sir(X, y, delta, d = K, time_is_log = TRUE),
    Cox        = matrix(cox_fit(X, y, delta), ncol = 1),
    `Cox-lasso`= matrix(cox_lasso_bic(X, y, delta), ncol = 1),
    stop("unknown method ", method))
}

## run a set of EDR estimators and return their accuracy
##   single index (K = 1): proximity criterion PC (and index correlation IC if
##   Sigma_diag is given); K = 2: recovery R and validity V.
edr_block <- function(dat, methods, A, K = 1L, Sigma_diag = NULL, lasso = TRUE,
                      prefix = "") {
  X <- dat$X; y <- dat$y; d <- dat$delta
  n <- nrow(X); p <- ncol(X)
  out <- list()
  score <- function(name, Bh) {
    if (K == 1L) {
      out[[paste0(prefix, "pc_", name)]] <<- proximity(Bh[, 1], dat$B[, 1])
      if (!is.null(Sigma_diag))
        out[[paste0(prefix, "ic_", name)]] <<- index_corr(Bh[, 1], dat$B[, 1], Sigma_diag)
    } else {
      sm <- if (all(Bh == 0)) c(recovery = 0, validity = 0) else subspace_metrics(Bh, dat$B)
      out[[paste0(prefix, "R_", name)]] <<- sm[["recovery"]]
      out[[paste0(prefix, "V_", name)]] <<- sm[["validity"]]
    }
  }
  for (m in methods) {
    if (m %in% c("SIR", "CSIR", "CP-SIR", "Cox") && p >= n) next
    if (m == "Cox-lasso" && !lasso) next
    t0 <- proc.time()[["elapsed"]]
    Bh <- switch(m,
      `CSIR-WLS`  = csir_on(X, y, d, A, K),
      `CP-SIRWLS` = cpsir_on(X, y, d, A, K),
      edr_fit(m, X, y, d, K))
    out[[paste0(prefix, "time_", m)]] <- proc.time()[["elapsed"]] - t0
    if (m == "Cox-lasso") {
      out[[paste0(prefix, "mms_lasso")]] <- mms_random_ties(Bh[, 1], dat$active)
      out[[paste0(prefix, "nz_lasso")]]  <- sum(Bh[, 1] != 0)
    }
    score(m, Bh)
  }
  out
}

## -------------------------------------------------------- ablation (Table S15)
## weighted leverage score of Zhong et al. (2023) on the observed time: the
## "naive" score that slices y and omits the factor exp(-Lambda)
wls_observed <- function(X, y, H = H_SLICES) {
  sv <- svd(center_cols(X))
  k <- length(sv$d)
  U <- sv$u[, seq_len(k), drop = FALSE]; V <- sv$v[, seq_len(k), drop = FALSE]
  cuts <- slice_cuts(y, H)
  Pi <- matrix(0, k, k)
  for (l in seq_len(H)) {
    sl <- y > cuts[l] & y <= cuts[l + 1L]
    if (!any(sl)) next
    Pi <- Pi + mean(sl) * tcrossprod(colMeans(U[sl, , drop = FALSE]))
  }
  rowSums((V %*% Pi) * V)
}
## complete-case score: the same score on the uncensored observations only
wls_complete_case <- function(X, y, delta, H = H_SLICES)
  wls_observed(X[delta == 1L, , drop = FALSE], y[delta == 1L], H)

## ------------------------------------------------- split / cross-fit (Table 6)
## D1 = first floor(n/2) observations, D2 = the others (Section 4.5)
split_crossfit <- function(dat) {
  X <- dat$X; y <- dat$y; d <- dat$delta
  n <- nrow(X); p <- ncol(X)
  D1 <- seq_len(n %/% 2L); D2 <- setdiff(seq_len(n), D1)
  half <- function(Dr, De) {                  # rank on Dr, CSIR on De
    s <- cwls_select(X[Dr, , drop = FALSE], y[Dr], d[Dr], K = 1L)
    list(A = s$A, b = csir_on(X[De, , drop = FALSE], y[De], d[De], s$A, 1L)[, 1])
  }
  h1 <- half(D1, D2); h2 <- half(D2, D1)
  b2 <- if (sum(h1$b * h2$b) < 0) -h2$b else h2$b    # sign alignment
  bcf <- (h1$b / sqrt(sum(h1$b^2)) + b2 / sqrt(sum(b2^2))) / 2
  list(pc_split = proximity(h1$b, dat$B[, 1]),
       pc_cross = proximity(bcf, dat$B[, 1]),
       size_split = length(h1$A),
       cover_split = as.numeric(all(dat$active %in% h1$A)))
}
