## Table S14: symmetric link log y0 = (beta'x)^2 + 0.5 e (Remark 5.1),
## Setting (1), n = 500, p = 80, 25 % censoring. N = 300; Cox-lasso 60.
## Reported: minimum model size of CWLS, CSIS and the Cox-lasso, and the
## proximity criterion of CSIR, Cox, Cox-lasso and CSIR-WLS.
source(file.path("reproduction", "lib", "common.R"))
for (p in p_grid(80))
  screen_edr_cell("symmetric", "symmetric", spec_ar1(0.5), "symmetric", "indep", p, 0.25,
                  N = 300, NL = 60, marginal = "CSIS",
                  edr = c("CSIR", "Cox", "Cox-lasso", "CSIR-WLS"))
