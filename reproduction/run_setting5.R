## Tables S8, S9, S10: block-exchangeable Setting (5), configurations
## (A1) alpha = (0.5,0.2,0.4), (A2) (0.8,0.2,0.4), (A3) (0.8,0.2,0.6), DGPs 5.1-5.3.
## Table S12: four further configurations under DGP 5.1, p in {80, 300}.
## N = 200 / 120 / 60; Cox-lasso on the first 30 / 20 / 6.
source(file.path("reproduction", "lib", "common.R"))
args   <- commandArgs(trailingOnly = TRUE)
part   <- if (length(args)) args else c("S8-S10", "S12")
NN <- c(`80` = 200, `300` = 120, `1500` = 60)
NL <- c(`80` = 30,  `300` = 20,  `1500` = 6)
if ("S8-S10" %in% part)
  for (p in p_grid(c(80, 300, 1500)))
    for (model in c("dgp1", "dgp2", "dgp3"))
      for (alpha in list(c(0.5, 0.2, 0.4), c(0.8, 0.2, 0.4), c(0.8, 0.2, 0.6)))
        for (target in c(0.15, 0.25, 0.40))
          screen_edr_cell("setting5", "setting5", spec_block(alpha), model, "indep", p, target,
                          N = NN[[as.character(p)]], NL = NL[[as.character(p)]])
if ("S12" %in% part)
  for (p in p_grid(c(80, 300)))
    for (alpha in list(c(0.2, 0.2, 0.2), c(0.8, 0.2, 0.2), c(0.8, 0.6, 0.6), c(0.5, 0.2, 0.8)))
      for (target in c(0.15, 0.25, 0.40))
        screen_edr_cell("setting5_more", "setting5", spec_block(alpha), "dgp1", "indep", p, target,
                        N = NN[[as.character(p)]], NL = NL[[as.character(p)]])
