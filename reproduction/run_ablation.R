## Table S15: ablation of the censoring correction. Minimum model size of the
## proposed CWLS score, of the naive score (slices the observed time, omits
## exp(-Lambda)) and of the complete-case score (uncensored observations only).
## DGP 5.1, Setting (1), n = 500, N = 300, target censoring 25 / 40 / 60 %
## (achieved 24-26 / 39-40 / 60 % in the paper).
source(file.path("reproduction", "lib", "common.R"))
for (p in p_grid(c(80, 300)))
  for (target in c(0.25, 0.40, 0.60))
    screen_edr_cell("ablation", "ablation", spec_ar1(0.5), "dgp1", "indep", p, target,
                    N = 300, marginal = character(0), edr = character(0),
                    extra = function(dat, s, r) list(
                      mms_naive = mms(wls_observed(dat$X, dat$y), dat$active),
                      mms_cc    = mms(wls_complete_case(dat$X, dat$y, dat$delta), dat$active)))
