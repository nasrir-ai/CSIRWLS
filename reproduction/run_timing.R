## Tables S16 (screening) and S17 (EDR estimation): average CPU time in seconds,
## DGP 5.3, Setting (1), n = 500, N = 100, p in {80, 1500}, 15/25/40 %.
## Run SEQUENTIALLY (one core) so that the timings are not distorted. SIR is
## timed at p = 1500 too (with the ridge of sir()), CSIR and Cox only at p = 80.
## Absolute times depend on the machine and on the BLAS library.
source(file.path("reproduction", "lib", "common.R"))
tm <- function(expr) { t0 <- proc.time(); force(expr); (proc.time() - t0)[["elapsed"]] }
for (p in p_grid(c(80, 1500)))
  for (target in c(0.15, 0.25, 0.40)) {
    spec <- spec_ar1(0.5); key <- cell_key("timing", "dgp3", spec, "indep", p, target)
    const <- calibrate_const_cell(spec, "dgp3", "indep", N_DESIGN, p, target, key)
    run_cell("timing", key, n_rep(100), cores = 1L,
      meta = list(model = "dgp3", design = spec$label, mech = "indep", p = p, target = target),
      fun = function(r) {
        dat <- gen_data(spec, "dgp3", "indep", N_DESIGN, p, const)
        X <- dat$X; y <- dat$y; d <- dat$delta
        out <- list(cens_rate = dat$cens_rate)
        out$t_CWLS <- tm(cwls(X, y, d, H = H_SLICES, r = 2L))
        out$t_CSIS <- tm(csis(X, y, d))
        out$t_SIR  <- tm(sir(X, y, H = H_SLICES, d = 1L))
        if (p < N_DESIGN) {
          out$t_CSIR <- tm(csir(X, y, d, H = H_SLICES, d = 1L))
          out$t_Cox  <- tm(cox_fit(X, y, d))
        }
        out$t_CSIRWLS <- tm({ s <- cwls_select(X, y, d); csir_on(X, y, d, s$A) })
        out
      })
  }
