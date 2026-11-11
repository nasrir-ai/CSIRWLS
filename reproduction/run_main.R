## Tables 1, 2 (DGP 5.1), S1 + S6 (DGP 5.2), S2 + S7 (DGP 5.3), and the
## Cox-lasso minimum model size of Section 5.4.
## Settings (1)-(2): AR(1) Gaussian design, rho = 0.5 / 0.8, n = 500,
## p in {80, 300, 1500}, censoring 15/25/40 %.
## N = 500 / 300 / 150 replications; Cox-lasso on the first 60 / 50 / 20.
## Usage: Rscript reproduction/run_main.R [dgp1 dgp2 dgp3]
source(file.path("reproduction", "lib", "common.R"))
args   <- commandArgs(trailingOnly = TRUE)
models <- if (length(args)) args else c("dgp1", "dgp2", "dgp3")
NN <- c(`80` = 500, `300` = 300, `1500` = 150)
NL <- c(`80` = 60,  `300` = 50,  `1500` = 20)
for (p in p_grid(c(80, 300, 1500)))
  for (model in models)
    for (rho in c(0.5, 0.8))
      for (target in c(0.15, 0.25, 0.40))
        screen_edr_cell("main", "main", spec_ar1(rho), model, "indep", p, target,
                        N = NN[[as.character(p)]], NL = NL[[as.character(p)]])
