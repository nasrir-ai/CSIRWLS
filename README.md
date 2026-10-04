# CSIRWLS

Censored Sliced Inverse Regression with Weighted Least Squares variable
selection (**CSIR-WLS**) — the R package accompanying the CSIR-WLS
paper (Nadia Asrir).

CSIR-WLS is a dimension-reduction method for right-censored survival
data. It combines Censored SIR (CSIR) with a Weighted Least Squares
(WLS) leverage-score variable-selection step, making it usable in the
high-dimensional setting (p > n). This package also bundles the
competing methods used for comparison in the paper (SIR, Cox, CSIS,
  SIRS, SIS), the simulation-data generators, evaluation metrics,
parallel   simulation runners, and the LaTeX table helpers
used to reproduce the paper's tables and figures.

## Installation

```r
# install.packages("devtools")
devtools::install_github("nasrir-ai/CSIRWLS")
```

Some real-data example scripts (`inst/examples/`) use additional
packages not required for the core methods; install everything with:

```r
devtools::install_github("nasrir-ai/CSIRWLS", dependencies = TRUE)
```

## Quick example

```r
library(CSIRWLS)
library(survival)

set.seed(1)
n <- 200; p <- 10
X <- matrix(rnorm(n * p), n, p)
beta_true <- c(1, -1, rep(0, p - 2))
time_true <- exp(as.numeric(X %*% beta_true) + rnorm(n))
cens <- exp(rnorm(n) + 1)
time <- pmin(time_true, cens)
event <- as.numeric(time_true <= cens)

fit <- cen.wls.sir2(X, time, event,
  n.slice1 = 5, n.slice0 = 5, n.slice = 10,
  cn1 = 0.1, cn2 = 1, choose.dir = FALSE, ndim = 1, c = 0.05
)

fit$select     # variables retained by the WLS screening step
fit$beta.hat   # estimated EDR direction (0 outside selected variables)
```

## Package layout

```
R/
  01-linear-algebra-utils.R   costeta, Dproj, projec, cov1, eigen.decomp, ...
  02-core-censored-sir.R      cen.sir, double.slice, double.slice2, edr.n, ...
  03-plotting.R               plot.csir/plot.ds, Kaplan-Meier helpers
  04-variable-selection-wls.R cen.wls, cen.wls.sir2, select_d_bic,
                              estimate_dim_chisq, ...
  05-screening-methods.R      CSIS, DCSIS, SIRS, SIS, apply_cwls, apply_csis
  06-data-simulation.R        X_rand, generate_X, generate_survival_times
  07-evaluation-metrics.R     trace_corr, subspace_dist, TPR/FPR/FDR, ...
  08-simulation-runners.R     parallel Monte-Carlo loops (EDR + screening)
  09-reporting-tables.R       xtable LaTeX table builders
  10-experimental.R           incomplete/experimental code (see file header)
inst/examples/                 cleaned, runnable real-data & simulation scripts
inst/extdata/                  GSE7390 data for the breast-cancer example
tests/testthat/                unit + smoke tests
```

## Real-data & simulation examples

Runnable, commented scripts reproducing the paper's analyses live in
`inst/examples/` (find them locally after installing with
`system.file("examples", package = "CSIRWLS")`, or just open them in
this repo):

- `01-dlbcl-real-data.R` — Cox vs CSIR vs CSIR-WLS on the DLBCL dataset
  (`ipred::DLBCL`). Table 15 / Figures 2-3.
- `02-pbc-bootstrap.R` — CSIR vs CSIR-WLS on `survival::pbc`, with
  bootstrap standard errors and 95% CIs. The EDR dimension is no longer
  hardcoded to 2: it's estimated first via `select_d_bic()` (with
  `estimate_dim_chisq()` printed alongside for comparison), and the
  rest of the script (tables, plots) generalizes to whatever dimension
  gets estimated.
- `03-breast-cancer-gse7390.R` — CSIS vs CWLS gene screening on GSE7390
  (needs `inst/extdata/GSE7390_transbig2006affy.RData`, not pushed to
  GitHub — see note below).
- `04-simulation-single-index.R` — Monte-Carlo comparison SIR/CSIR/Cox/
  CSIR-WLS, single EDR direction (d = 1).
- `05-simulation-multi-index.R` — same comparison for d > 1, with fixed
  and automatically-selected dimension.
- `06-screening-cwls-vs-csis.R` — TPR/FPR/FDR comparison of CWLS vs
  CSIS in very high dimension (p up to 10000), including robustness
  checks under non-elliptical designs and block-correlation structures.
- `07-dimension-recovery-multi-index.R` — Section 5.3.2 evidence:
  recovery rate of the true structural dimension d=2 in a multi-index
  censored model, comparing the chi-squared test (`estimate_dim_chisq()`,
  Li 1991) against the BIC-type criterion (`select_d_bic()`), across
  p << n and p >> n regimes and two censoring rates. Includes an
  "easy" additive-DGP sanity check to isolate pipeline bugs from
  genuine detection-power limits.








  ```r
  devtools::document()  # regenerates NAMESPACE/man from the roxygen comments
  devtools::load_all()
  devtools::test()
  devtools::check()
  ```

 
 
 
