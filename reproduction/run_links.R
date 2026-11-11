## Table S13: single-index links away from proportional hazards (saturating,
## threshold, heteroscedastic), Setting (1), n = 500, 25 % censoring.
## Reported: squared proximity criterion PC^2. N = 300 / 200; Cox-lasso 60 / 50.
source(file.path("reproduction", "lib", "common.R"))
NN <- c(`80` = 300, `300` = 200); NL <- c(`80` = 60, `300` = 50)
for (p in p_grid(c(80, 300)))
  for (link in c("saturating", "threshold", "heteroscedastic"))
    screen_edr_cell("links", "links", spec_ar1(0.5), link, "indep", p, 0.25,
                    N = NN[[as.character(p)]], NL = NL[[as.character(p)]],
                    marginal = character(0),
                    edr = c("CSIR", "CP-SIR", "Cox", "Cox-lasso", "CSIR-WLS"))
