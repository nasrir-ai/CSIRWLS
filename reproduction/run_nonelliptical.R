## Tables S18 (screening, 25 % censoring) and S19 (EDR, 15/25/40 %) under the
## non-elliptical designs of Section 5.8.2: skew-normal x = delta |z0| 1 +
## sqrt(1 - delta^2) u and skew-t (the skew-normal vector divided by
## sqrt(W / nu), W ~ chi2_nu), u ~ N(0, AR(1) rho = 0.5), delta = 0.8, nu = 5.
## DGP 5.1, n = 500, N = 200 / 120 / 60, Cox-lasso on the first 30 / 20 / 6.
## The Gaussian rows of S18-S19 are those of Tables 1-2 (run_main.R).
source(file.path("reproduction", "lib", "common.R"))
NN <- c(`80` = 200, `300` = 120, `1500` = 60)
NL <- c(`80` = 30,  `300` = 20,  `1500` = 6)
for (p in p_grid(c(80, 300, 1500)))
  for (kind in c("skewnormal", "skewt"))
    for (target in c(0.15, 0.25, 0.40))
      screen_edr_cell("nonelliptical", "nonelliptical", spec_skew(kind), "dgp1", "indep", p, target,
                      N = NN[[as.character(p)]], NL = NL[[as.character(p)]], marginal = "CSIS")
