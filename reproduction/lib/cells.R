## =============================================================================
## cells.R -- the generic simulation cell used by most drivers
##
## One replication = one data set, then
##   (i)  screening: CWLS + G_c (TPR, FPR, FDR, size, MMS) and the marginal
##        statistics CSIS, CRIS, cMV-SIS at floor(n / log n) (TPR, MMS);
##  (ii)  EDR: the requested estimators (PC, or IC, or R/V when K = 2);
##        Cox-lasso only on the first NL replications of the cell.
## =============================================================================

N_DESIGN <- 500L

cell_key <- function(prefix, model, spec, mech, p, target)
  sprintf("%s|%s|%s|%s|p%d|c%.2f", prefix, model, spec$label, mech, p, target)

marginal_scores <- function(m, X, y, d) switch(m,
  CSIS = csis(X, y, d), CRIS = cris(X, y, d), `cMV-SIS` = cmv_sis(X, y, d),
  stop("unknown screening statistic ", m))

screen_edr_cell <- function(driver, prefix, spec, model, mech = "indep", p, target,
                            N, NL = 0L, n = N_DESIGN,
                            marginal = c("CSIS", "CRIS", "cMV-SIS"),
                            edr = c("SIR", "CSIR", "CP-SIR", "Cox", "Cox-lasso", "CSIR-WLS"),
                            K = 1L, ic = FALSE, extra = NULL) {
  key <- cell_key(prefix, model, spec, mech, p, target)
  N <- n_rep(N); NL <- if (NL > 0) n_lasso(NL) else 0L
  warm_design(spec, p)
  const <- calibrate_const_cell(spec, model, mech, n, p, target, key)
  Sdiag <- if (ic) diag(sigma_sub(spec, n, p, seq_len(p))) else NULL
  keep <- as.integer(n / log(n))
  fun <- function(r) {
    dat <- gen_data(spec, model, mech, n, p, const)
    X <- dat$X; y <- dat$y; d <- dat$delta; act <- dat$active
    out <- list(cens_rate = dat$cens_rate)
    t0 <- proc.time()[["elapsed"]]
    s <- cwls_select(X, y, d, K = K)
    out$time_CWLS <- proc.time()[["elapsed"]] - t0
    rt <- sel_rates(s$A, act, p)
    out$cwls_tpr <- rt[["tpr"]]; out$cwls_fpr <- rt[["fpr"]]; out$cwls_fdr <- rt[["fdr"]]
    out$cwls_size <- rt[["size"]]; out$cwls_k <- s$k
    out$mms_CWLS <- mms(s$omega, act)
    for (m in marginal) {
      t0 <- proc.time()[["elapsed"]]
      sc <- marginal_scores(m, X, y, d)
      out[[paste0("time_", m)]] <- proc.time()[["elapsed"]] - t0
      out[[paste0("tpr_", m)]] <- sel_rates(top_k(sc, keep), act, p)[["tpr"]]
      out[[paste0("mms_", m)]] <- mms(sc, act)
    }
    if (length(edr))
      out <- c(out, edr_block(dat, edr, s$A, K = K, Sigma_diag = Sdiag, lasso = r <= NL))
    if (!is.null(extra)) out <- c(out, extra(dat, s, r))
    out
  }
  run_cell(driver, key, N, fun,
           meta = list(key = key, model = model, design = spec$label, mech = mech, p = p,
                       target = target, const = const))
}

## restrict the dimensions run, e.g. CSIRWLS_P="80,300"
p_grid <- function(ps) {
  only <- Sys.getenv("CSIRWLS_P", "")
  if (!nzchar(only)) return(ps)
  intersect(ps, as.integer(strsplit(only, ",")[[1]]))
}
