## =============================================================================
## Section 6.1 -- DLBCL: Table 9, Figures 1 and 2, the CWLS selected set
##
## Data: ipred::DLBCL (40 patients; 38 complete cases after removing the two
## missing IPI values; 18 censored). Covariates: IPI and MGEc.1-MGEc.10
## (p = 11). delta = cens (1 = death).
## Fits: Cox (all 11 covariates), CSIR (all 11), CSIR-WLS (CWLS + G_c, then
## CSIR), each through fit_risk_index(), so that a larger index means a higher
## risk. Table 9 gives the in-sample and the cross-validated C-index and IBS
## (10 folds, 20 repetitions, screening repeated inside each training fold):
## tables/table9_dlbcl.csv and .tex. Figure 1 plots the standardised
## indices on a common scale; Figure 2 the Kaplan-Meier curves of the tertile
## groups. The bandwidth of the kernel weights follows rule (27) of the paper.
##
## Covariates are used on their original scale by default; CSIRWLS_DLBCL_STD=1
## standardises them. Both choices give the same CWLS ranking.
##
## G_c retains 4 covariates (IPI, MGEc.2, MGEc.6, MGEc.8).
## CSIRWLS_DLBCL_KEEP=k forces the k top-ranked covariates.
## =============================================================================
source(file.path("reproduction", "applications", "app_common.R"))
if (!requireNamespace("ipred", quietly = TRUE))
  stop("package 'ipred' is required: install.packages('ipred')")
utils::data("DLBCL", package = "ipred", envir = environment())
vars <- c("IPI", paste0("MGEc.", 1:10))
d <- DLBCL[stats::complete.cases(DLBCL[, c("time", "cens", vars)]), ]
X <- as.matrix(d[, vars])
if (Sys.getenv("CSIRWLS_DLBCL_STD", "0") == "1") X <- scale(X)
y <- d$time; delta <- as.integer(d$cens)
cat(sprintf("DLBCL: n = %d, censored = %d (%.0f %%)\n", nrow(X), sum(delta == 0), 100 * mean(delta == 0)))

keep <- Sys.getenv("CSIRWLS_DLBCL_KEEP", "")
n_keep <- if (nzchar(keep)) as.integer(keep) else NULL
om <- cwls(X, y, delta)
A <- order(om, decreasing = TRUE)[seq_len(if (is.null(n_keep)) bic_cut(om, nrow(X), ncol(X)) else n_keep)]
cat("CWLS ranking:", vars[order(om, decreasing = TRUE)], "\n")
cat("retained by", if (is.null(n_keep)) "G_c" else "n_keep", ":", vars[sort(A)], "\n")

B <- list(COX = fit_risk_index("Cox", X, y, delta),
          CSIR = fit_risk_index("CSIR", X, y, delta),
          `CSIR-WLS` = fit_risk_index("CSIR-WLS", X, y, delta, n_keep = n_keep))
tab9 <- do.call(rbind, lapply(names(B), function(m) {
  ev <- insample_eval(B[[m]], X, y, delta)
  data.frame(Method = m, `C-index` = round(ev[["c_index"]], 3), IBS = round(ev[["ibs"]], 3),
             covariates = sum(B[[m]] != 0), check.names = FALSE)
}))
print(tab9)

## log-rank statistics of the three risk groups (Figure 2)
lr <- do.call(rbind, lapply(names(B), function(m) {
  t <- logrank_test(y, delta, risk_groups(as.vector(X %*% B[[m]])))
  data.frame(Method = m, chi2 = round(t$statistic, 1), df = t$df, p.value = signif(t$p.value, 2))
}))
print(lr); write_csv(lr, "dlbcl_logrank")
write_csv(data.frame(variable = vars, omega = om, rank = rank(-om),
                     retained = seq_along(vars) %in% A,
                     beta_cox = B$COX, beta_csir = B$CSIR, beta_csirwls = B$`CSIR-WLS`), "dlbcl_coefficients")

## Figure 1: standardised risk indices on a common scale, in colour
B_fig <- setNames(B, c("Cox", "CSIR", "CSIR-WLS"))
grDevices::pdf(file.path(FIG, "figure1_dlbcl_projections.pdf"), width = 10, height = 3.8)
plot_index_panels(B_fig, X, y, delta); grDevices::dev.off()
grDevices::png(file.path(FIG, "figure1_dlbcl_projections.png"), width = 2000, height = 760, res = 200)
plot_index_panels(B_fig, X, y, delta); grDevices::dev.off()
## Figure 2: Kaplan-Meier curves of the tertile groups, log-rank in the titles
grDevices::pdf(file.path(FIG, "figure2_dlbcl_km.pdf"), width = 11, height = 3.9)
plot_km_panels(B_fig, X, y, delta); grDevices::dev.off()
grDevices::png(file.path(FIG, "figure2_dlbcl_km.png"), width = 2200, height = 780, res = 200)
plot_km_panels(B_fig, X, y, delta); grDevices::dev.off()

## cross-validation (10 folds, 20 repetitions, screening repeated inside each
## training fold): reported in Table 9
cv <- cv_survival(X, y, delta, methods = c("Cox", "CSIR", "CSIR-WLS"), K = 10L,
                  reps = n_rep(20), n_keep = n_keep)
print(cv); write_csv(cv, "dlbcl_cv")

## Table 9: in sample and cross-validated C-index and IBS
f3 <- function(x) sprintf("%.3f", x)
t9 <- data.frame(Method = c("Cox", "CSIR", "CSIR-WLS"), Covariates = tab9$covariates,
                 C_in = f3(tab9$`C-index`), IBS_in = f3(tab9$IBS),
                 C_cv = sprintf("%.3f (%.3f)", cv$c_index, cv$c_se),
                 IBS_cv = sprintf("%.3f (%.3f)", cv$ibs, cv$ibs_se))
write_csv(t9, "table9_dlbcl")
tex <- c("\\begin{tabular}{lccccc}", "\\toprule",
         " & & \\multicolumn{2}{c}{In sample} & \\multicolumn{2}{c}{Cross-validated} \\\\",
         "Method & Covariates & C-index & IBS & C-index & IBS \\\\", "\\midrule",
         apply(t9, 1, function(r) paste0(paste(r, collapse = " & "), " \\\\")),
         "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(TAB, "table9_dlbcl.tex"))
print(t9)
