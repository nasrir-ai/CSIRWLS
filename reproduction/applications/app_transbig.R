## =============================================================================
## Section 6.3 -- TRANSBIG breast cancer data: Table 12, the in-sample C-index
## as a function of the number of retained probes, and the bootstrap stability
## (Jaccard) of the CWLS and CSIS selections.
##
## Data: Bioconductor experiment package breastCancerTRANSBIG (Desmedt et al.,
## 2007): 198 patients, 22,283 probes (hgu133a) + age, tumour size, grade
## (p = 22,286); endpoint distant metastasis (t.dmfs, e.dmfs; 62 events).
##   BiocManager::install(c("Biobase", "breastCancerTRANSBIG"))
## Alternatively CSIRWLS_TRANSBIG_RDS=<file.rds> with list(X, y, delta).
## Missing clinical values (if any) are replaced by the median; the script
## prints how many.
##
##  * Table 12: 10-fold CV (20 fold assignments, 2 for the Cox-lasso),
##    screening and standardisation inside every training fold (cv_survival,
##    standardise = TRUE); in-sample C-index and log-rank statistic of three
##    risk groups (33rd / 67th percentiles, 2 df) on the full standardised
##    sample; |A| = average number of covariates entering the second stage.
##    Pipelines: CWLS -> CSIR (CSIR-WLS), CWLS -> Cox, CSIS -> Cox
##    (floor(n / log n) = 37 covariates), Cox-lasso on all covariates.
##  * CWLS -> CSIR with G_c, 100 and 300 retained probes: in sample and CV.
##  * Stability: B = 500 bootstrap resamples, both rules keeping 37
##    covariates; mean pairwise Jaccard similarity, with the standard error of
##    a U-statistic of order 2 (2 sd(J_b) / sqrt(B), J_b the mean similarity of
##    resample b with the others).
## Cost: one Cox-lasso path is about 80 s at this dimension; the whole script
## takes a few hours on one core (CSIRWLS_CORES is used by the bootstrap).
## =============================================================================
source(file.path("reproduction", "applications", "app_common.R"))
R_CV <- n_rep(20); R_CV_LASSO <- n_rep(2); B_BOOT <- n_rep(500)

load_transbig <- function() {
  f <- Sys.getenv("CSIRWLS_TRANSBIG_RDS", "")
  if (nzchar(f)) return(readRDS(f))
  for (pk in c("Biobase", "breastCancerTRANSBIG"))
    if (!requireNamespace(pk, quietly = TRUE))
      stop("package '", pk, "' is required: BiocManager::install(c('Biobase', 'breastCancerTRANSBIG'))")
  e <- new.env()
  utils::data("transbig", package = "breastCancerTRANSBIG", envir = e)
  ex <- t(Biobase::exprs(e$transbig)); cl <- Biobase::pData(e$transbig)
  ok <- !is.na(cl$t.dmfs) & !is.na(cl$e.dmfs)
  clin <- as.matrix(cl[, c("age", "size", "grade")])
  for (j in seq_len(ncol(clin))) {
    nas <- is.na(clin[, j])
    if (any(nas & ok)) message(sum(nas & ok), " missing value(s) of ", colnames(clin)[j], " replaced by the median")
    clin[nas, j] <- stats::median(clin[ok, j], na.rm = TRUE)
  }
  X <- cbind(ex, clin)[ok, ]
  list(X = X, y = as.numeric(cl$t.dmfs[ok]), delta = as.integer(cl$e.dmfs[ok]))
}
dat <- load_transbig()
X <- as.matrix(dat$X); y <- as.numeric(dat$y); delta <- as.integer(dat$delta)
keep_out <- apply(X, 2, stats::sd) > 0                 # constant columns cannot be standardised
X <- X[, keep_out, drop = FALSE]
n <- nrow(X); p <- ncol(X)
cat(sprintf("TRANSBIG: n = %d, p = %d, events = %d (%.1f %%)\n", n, p, sum(delta), 100 * mean(delta)))
Xs <- scale(X)
csis_keep <- as.integer(n / log(n))

## ------------------------------------------------------------ Table 12
pipes <- c(`CWLS -> CSIR (CSIR-WLS)` = "CSIR-WLS", `CWLS -> Cox` = "CWLS+Cox",
           `CSIS -> Cox` = "CSIS+Cox", `Cox-lasso` = "Cox-lasso")
if (Sys.getenv("CSIRWLS_TRANSBIG_LASSO", "1") == "0") pipes <- pipes[pipes != "Cox-lasso"]   # quick tests
cv <- cv_survival(X, y, delta, methods = unname(pipes), K = 10L,
                  reps = c(R_CV, R_CV, R_CV, R_CV_LASSO)[seq_along(pipes)], standardise = TRUE, insample = FALSE)
ins <- do.call(rbind, lapply(unname(pipes), function(m) {
  if (m == "Cox-lasso") return(data.frame(c_in = NA, chi2 = NA, p_lr = NA))   # "--" in the paper
  b <- fit_risk_index(m, Xs, y, delta)
  lr <- logrank_test(y, delta, risk_groups(as.vector(Xs %*% b)))
  data.frame(c_in = c_index(as.vector(Xs %*% b), y, delta), chi2 = lr$statistic, p_lr = lr$p.value)
}))
tab12 <- data.frame(Pipeline = names(pipes),
                    `C-index CV` = sprintf("%.3f (%.3f)", cv$c_index, cv$c_se),
                    `IBS CV` = sprintf("%.3f (%.3f)", cv$ibs, cv$ibs_se),
                    `C-index in sample` = ifelse(is.na(ins$c_in), "--", sprintf("%.3f", ins$c_in)),
                    `log-rank chi2` = ifelse(is.na(ins$chi2), "--", sprintf("%.1f (p = %.2g)", ins$chi2, ins$p_lr)),
                    `|A|` = round(cv$n_selected), check.names = FALSE)
print(tab12); write_csv(tab12, "table12_transbig")

## ---------------------------------- CWLS -> CSIR with G_c, 100 and 300 probes
nk <- list(`G_c` = NULL, `100` = 100L, `300` = 300L)
size_tab <- do.call(rbind, lapply(names(nk), function(k) {
  b <- fit_risk_index("CSIR-WLS", Xs, y, delta, n_keep = nk[[k]])
  cvk <- cv_survival(X, y, delta, methods = "CSIR-WLS", K = 10L, reps = R_CV,
                     standardise = TRUE, insample = FALSE, n_keep = nk[[k]])
  data.frame(retained = k, size = sum(b != 0), c_insample = round(c_index(as.vector(Xs %*% b), y, delta), 3),
             c_cv = sprintf("%.3f (%.3f)", cvk$c_index, cvk$c_se))
}))
print(size_tab); write_csv(size_tab, "transbig_retained_size")

## ------------------------------------------------------------- stability
jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
one_boot <- function(b) {
  set.seed(rep_seed("transbig-jaccard", b))
  i <- sample(n, replace = TRUE)
  Xb <- scale(X[i, ]); Xb[!is.finite(Xb)] <- 0
  list(cwls = top_k(cwls(Xb, y[i], delta[i]), csis_keep),
       csis = top_k(csis(Xb, y[i], delta[i]), csis_keep))
}
sel <- if (CFG$cores > 1) {
  parallel::mclapply(seq_len(B_BOOT), one_boot, mc.cores = CFG$cores)
} else lapply(seq_len(B_BOOT), one_boot)
stab <- do.call(rbind, lapply(c("cwls", "csis"), function(rule) {
  S <- lapply(sel, `[[`, rule); B <- length(S)
  J <- matrix(NA_real_, B, B)
  for (a in seq_len(B - 1)) for (b in (a + 1):B) J[a, b] <- J[b, a] <- jaccard(S[[a]], S[[b]])
  Jb <- rowMeans(J, na.rm = TRUE)
  data.frame(rule = toupper(rule), B = B, retained = csis_keep,
             mean_jaccard = round(mean(J[upper.tri(J)]), 3), se = round(2 * stats::sd(Jb) / sqrt(B), 3))
}))
print(stab); write_csv(stab, "transbig_jaccard")
