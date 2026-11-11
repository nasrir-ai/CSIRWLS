## Tables 3 and 4: two-dimensional EDR space, DGP 5.4 (d = 2), Settings (1)-(2),
## n = 500, 25 % censoring. Recovery R and validity V of CSIR, CP-SIR, Cox,
## Cox-lasso, CP-SIRWLS (CP-SIR on the CWLS set) and CSIR-WLS, all given d = 2;
## minimum model size of CWLS, CSIS and the Cox-lasso ranking.
## N = 300 / 200 / 100; Cox-lasso on the first 60 / 50 / 30.
## The noise level of DGP 5.4 (default 0.5, see lib/designs.R) is set by
## CSIRWLS_DGP4_SIGMA and is part of the cell key, so runs with different
## values are stored side by side; make_tables.R uses the current value.
source(file.path("reproduction", "lib", "common.R"))
NN <- c(`80` = 300, `300` = 200, `1500` = 100)
NL <- c(`80` = 60,  `300` = 50,  `1500` = 30)
for (p in p_grid(c(80, 300, 1500)))
  for (rho in c(0.5, 0.8))
    screen_edr_cell("dgp4", sprintf("dgp4-sigma%g", DGP4_SIGMA), spec_ar1(rho), "dgp4", "indep", p, 0.25,
                    N = NN[[as.character(p)]], NL = NL[[as.character(p)]],
                    marginal = "CSIS", K = 2L,
                    edr = c("CSIR", "CP-SIR", "Cox", "Cox-lasso", "CP-SIRWLS", "CSIR-WLS"))
