## Tables S21 and S22: the floor b_n of condition (C.9)(c) in the kernel weight
## estimator (argument floor_mult of cwls(); 0 = no floor, the default).
##  S21: DGP 5.1, 25 % censoring, mult in {0, 0.1, 1}: Gaussian AR(1) p = 80, 300,
##       skew-normal p = 80, skew-t p = 80. N = 200 (p = 80), 120 (p = 300),
##       as in the other secondary designs.
##  S22: the main cells (DGPs 5.1-5.3, rho = 0.5/0.8, p = 80/300, 25/40 %), same
##       seeds and N as Tables 1, 2, S1-S10 (same cell keys as run_main.R), with
##       mult = 0 ("none") and mult = 1 ("b_n"). The "none" columns reproduce the
##       CWLS / CSIR-WLS entries of run_main.R exactly.
source(file.path("reproduction", "lib", "common.R"))
args <- commandArgs(trailingOnly = TRUE)
part <- if (length(args)) args else c("S21", "S22")
floor_extra <- function(mults) function(dat, s, r) {
  out <- list()
  for (m in mults) {
    sm <- cwls_select(dat$X, dat$y, dat$delta, K = 1L, floor_mult = m)
    tag <- gsub("\\.", "", sprintf("m%g", m))
    out[[paste0("size_", tag)]] <- length(sm$A)
    out[[paste0("tpr_", tag)]]  <- sel_rates(sm$A, dat$active, ncol(dat$X))[["tpr"]]
    out[[paste0("mms_", tag)]]  <- mms(sm$omega, dat$active)
    out[[paste0("pc_", tag)]]   <- proximity(csir_on(dat$X, dat$y, dat$delta, sm$A)[, 1], dat$B[, 1])
  }
  out
}
if ("S21" %in% part) {
  cells <- list(list(spec_ar1(0.5), 80), list(spec_ar1(0.5), 300),
                list(spec_skew("skewnormal"), 80), list(spec_skew("skewt"), 80))
  for (cl in cells) {
    p <- cl[[2]]
    if (!(p %in% p_grid(p))) next
    screen_edr_cell("floor_S21", "floor", cl[[1]], "dgp1", "indep", p, 0.25,
                    N = if (p == 80) 200 else 120, marginal = character(0), edr = character(0),
                    extra = floor_extra(c(0, 0.1, 1)))
  }
}
if ("S22" %in% part) {
  NN <- c(`80` = 500, `300` = 300)
  for (p in p_grid(c(80, 300)))
    for (model in c("dgp1", "dgp2", "dgp3"))
      for (rho in c(0.5, 0.8))
        for (target in c(0.25, 0.40))
          screen_edr_cell("floor_S22", "main", spec_ar1(rho), model, "indep", p, target,
                          N = NN[[as.character(p)]], marginal = character(0), edr = character(0),
                          extra = floor_extra(c(0, 1)))
}
