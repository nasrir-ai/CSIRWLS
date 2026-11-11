## Table 7 (Section 5.8.2): which departure from the theory is responsible for
## the degradation under the Azzalini designs?  DGP 5.1, n = 500, 25 % censoring,
## p in {80, 300}; N = 200 / 120, Cox-lasso on the first 30 / 20 replications.
##
## Seven designs per dimension (see table7_rows() in lib/designs.R):
##   Gaussian AR(1) reference ............ cell of run_main.R        (shared)
##   skew-normal, skew-t5 ................ cells of run_nonelliptical.R (shared)
##   for each Azzalini design, two controls (computed here):
##     elliptical control     Gaussian with exactly the Azzalini covariance (C.1 holds)
##     non-elliptical control x = Sigma_rho^{1/2} w, independent standardised skew
##                            components (C.1 fails, covariance = reference)
## The shared cells are run with the SAME arguments as their own drivers, so
## nothing is computed twice whatever the order in which the drivers are run.
source(file.path("reproduction", "lib", "common.R"))
NN  <- c(`80` = 200, `300` = 120); NL <- c(`80` = 30, `300` = 20)
NNr <- c(`80` = 500, `300` = 300); NLr <- c(`80` = 60, `300` = 50)   # run_main.R
for (p in p_grid(c(80, 300))) {
  P <- as.character(p)
  for (r in table7_rows()) {
    if (r$driver == "main")
      screen_edr_cell("main", "main", r$spec, "dgp1", "indep", p, 0.25,
                      N = NNr[[P]], NL = NLr[[P]])
    else if (r$driver == "nonelliptical")
      screen_edr_cell("nonelliptical", "nonelliptical", r$spec, "dgp1", "indep", p, 0.25,
                      N = NN[[P]], NL = NL[[P]], marginal = "CSIS")
    else
      screen_edr_cell("table7", "table7", r$spec, "dgp1", "indep", p, 0.25,
                      N = NN[[P]], NL = NL[[P]], marginal = character(0),
                      edr = c("CSIR-WLS", "Cox-lasso"))
  }
}
