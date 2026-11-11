## Table 8 (Section 5.8.3): sign cancellation.  Setting (5), DGP 5.1, n = 500,
## 25 % censoring, |S| = 5, magnitudes U(0.5, 1.5); configurations
## (A1) alpha = (0.5, 0.2, 0.4), (A2) (0.8, 0.2, 0.4), (A3) (0.8, 0.2, 0.6),
## (A4) (0.2, 0.2, 0.2); signs: all positive, or two negative at positions
## redrawn in every replication.  p = 300 (N = 100) and p = 1500 (N = 50).
##
## For every replication the two POPULATION margins are stored
##   m_half = min_S |(Sigma^{1/2} beta)_j| / max_{S^c} |(Sigma^{1/2} beta)_j|
##   m_one  = min_S |(Sigma beta)_j|       / max_{S^c} |(Sigma beta)_j|
## together with the MMS and the TPR of CWLS (stopping rule G_c) and of CSIS,
## CRIS, cMV-SIS (floor(n / log n) = 80 predictors).  No EDR step is needed.
source(file.path("reproduction", "lib", "common.R"))
NN <- c(`300` = 100, `1500` = 50)
configs <- list(A1 = c(0.5, 0.2, 0.4), A2 = c(0.8, 0.2, 0.4), A3 = c(0.8, 0.2, 0.6),
                A4 = c(0.2, 0.2, 0.2))
for (p in p_grid(c(300, 1500)))
  for (cfg in names(configs))
    for (neg in c(0L, 2L)) {
      spec <- spec_block(configs[[cfg]], neg)
      screen_edr_cell("signs", "signs", spec, "dgp1", "indep", p, 0.25,
                      N = NN[[as.character(p)]], edr = character(0),
                      extra = function(dat, s, r) as.list(block_margins(spec, ncol(dat$X), dat$B[, 1])))
    }
