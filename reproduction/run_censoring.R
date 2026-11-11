## Table 5: censoring mechanism. DGP 5.1, Setting (1) (rho = 0.5), n = 500,
## 25 % censoring, N = 400; independent, EDR-dependent and misspecified
## censoring (see lib/designs.R).
source(file.path("reproduction", "lib", "common.R"))
for (p in p_grid(c(80, 300)))
  for (mech in c("indep", "edr", "misspec"))
    screen_edr_cell("censoring", "censoring", spec_ar1(0.5), "dgp1", mech, p, 0.25,
                    N = 400, marginal = "CSIS", edr = "CSIR-WLS")
