## Table 6 (Section 5.7): full-sample, split-sample and cross-fitted CSIR-WLS.
## DGP 5.1, Setting (1), n = 500, 25 % censoring, N = 300 / 200 / 100.
## Split: CWLS ranking and G_c on the first n/2 observations, CSIR on the other
## n/2. Cross-fitted (Remark 4.3): the two split estimators with the roles of the
## halves exchanged, sign-aligned and averaged.
source(file.path("reproduction", "lib", "common.R"))
NN <- c(`80` = 300, `300` = 200, `1500` = 100)
for (p in p_grid(c(80, 300, 1500)))
  screen_edr_cell("split", "split", spec_ar1(0.5), "dgp1", "indep", p, 0.25,
                  N = NN[[as.character(p)]], marginal = character(0), edr = "CSIR-WLS",
                  extra = function(dat, s, r) split_crossfit(dat))
