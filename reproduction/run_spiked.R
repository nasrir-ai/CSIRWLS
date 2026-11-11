## Tables S3, S4, S5 (screening, DGPs 5.1-5.3) and S11 (EDR, PC and index
## correlation IC) under the spiked Settings (3)-(4).
## N = 200 / 120 / 60; Cox-lasso on the first 30 / 20 / 6.
## Table S11 shows DGP 5.1 at 25 % censoring; the other cells are computed too
## ("in the code output" in the paper).
source(file.path("reproduction", "lib", "common.R"))
args   <- commandArgs(trailingOnly = TRUE)
models <- if (length(args)) args else c("dgp1", "dgp2", "dgp3")
NN <- c(`80` = 200, `300` = 120, `1500` = 60)
NL <- c(`80` = 30,  `300` = 20,  `1500` = 6)
for (p in p_grid(c(80, 300, 1500)))
  for (model in models)
    for (setting in 3:4)
      for (target in c(0.15, 0.25, 0.40)) {
        do_edr <- target == 0.25           # S11 is at 25 % censoring
        screen_edr_cell("spiked", "spiked", spec_spiked(setting), model, "indep", p, target,
                        N = NN[[as.character(p)]], NL = if (do_edr) NL[[as.character(p)]] else 0,
                        edr = if (do_edr) c("SIR", "CSIR", "CP-SIR", "Cox", "Cox-lasso", "CSIR-WLS") else character(0),
                        ic = TRUE)
      }
