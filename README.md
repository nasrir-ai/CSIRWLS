# CSIRWLS

Censored Sliced Inverse Regression with Weighted Least Squares variable
selection (**CSIR-WLS**) — the R package accompanying the CSIR-WLS
paper (Nadia Asrir).

CSIR-WLS is a dimension-reduction method for right-censored survival
data. It combines Censored SIR (CSIR) with a Weighted Least Squares
(WLS) leverage-score variable-selection step, making it usable in the
high-dimensional setting (p > n). This package also bundles the
competing methods used for comparison in the paper (SIR, Cox, CSIS,
DCSIS, SIRS, SIS), the simulation-data generators, evaluation metrics,
parallel Monte-Carlo simulation runners, and the LaTeX table helpers
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

## About the reorganization

Your original 10 files (`myfunctions.R` + 9 analysis/simulation
scripts) had the same function names defined differently across files
(e.g. three near-identical `run_simulation_parallel()`, three
`latex_table()`), one dead nested copy-paste inside `cen.wls()`, and
hardcoded personal paths (`~/Documents/phd /R codes /...`). This
package:

- Deduplicates every function to a single canonical version.
- Renames colliding names explicitly (see the top of
  `R/CSIRWLS-package.R` and `R/08-simulation-runners.R` for the full
  old-name -> new-name table).
- Removes the dead nested `cen.wls()` redefinition (it lived inside an
  `if` branch that could never be reached — the fix does not change any
  numeric output).
- Replaces every hardcoded path with `library(CSIRWLS)` — no more
  `source("myfunctions.R")`.
- Documents every exported function with roxygen2 (in French, matching
  your original comments).

`wls.sir.fdr()` (FDR-controlled variable selection) was left as-is in
`R/10-experimental.R`: it references objects (`t`, `Ta`, `alf`, `W`,
`q`) that are never defined inside the function itself in the original
code, so it will error unless you supply them from the calling
environment. See that file's documentation for what would need
clarifying if you want to actually use it.

## Dimension selection: use `select_d_bic()`, not an ad hoc reimplementation

While debugging `07-dimension-recovery-multi-index.R`, an early version
of the BIC-type dimension criterion filtered eigenvalues with
`eigvals[eigvals > 0]` (no numerical tolerance). With only a handful of
slices, the true rank of the slice-mean covariance matrix is bounded by
`n.slice - 1`; beyond that, "eigenvalues" coming out of
`eigen.decomp()` are floating-point noise (~1e-15), not exact zeros,
and slipped through that filter — inflating the candidate dimension
space and making the estimated dimension collapse to (almost) the same
value regardless of the data. `select_d_bic()` (already in
`R/04-variable-selection-wls.R`) filters with a `> 1e-10` tolerance and
does not have this problem — always prefer it over a hand-rolled
filter. `estimate_dim_chisq()` (the sequential chi-squared alternative,
Li 1991) is also exported for direct comparison; it tends to
under-estimate d when one EDR direction only enters the model through
a non-monotone/weakly-asymmetric term (see the ratio-form DGP in
`07-dimension-recovery-multi-index.R`), and loses more power as
censoring increases.

## Known limitations / things to double-check

- **R was not available in the environment used to build this
  package**, so none of this has been run through `R CMD check` or
  `devtools::test()`. Before you rely on it, please run, from the
  package root:

  ```r
  devtools::document()  # regenerates NAMESPACE/man from the roxygen comments
  devtools::load_all()
  devtools::test()
  devtools::check()
  ```

  and fix anything that comes up — I did a careful manual review (brace
  matching, consistent exports, no leftover `source()`/hardcoded paths)
  but a real R session is the only way to catch every typo.
- `X_rand(type = "block")` needs `q` (number of active variables); the
  multi-index simulation runners (`run_simulation_multi*`) don't thread
  a `q` argument through to `generate_X()`, exactly like your original
  script — `X_type = "block"` will fail there (caught by `tryCatch`).
  It works fine in the single-index runners and in
  `run_alpha_study()`.
- `inst/extdata/GSE7390_transbig2006affy.RData` (~33 MB) is excluded
  from git via `.gitignore` to keep the repository light. It's on your
  disk after cloning your own copy of this folder, but **will not be
  on GitHub** unless you set up [Git LFS](https://git-lfs.com) or host
  it elsewhere and update `03-breast-cancer-gse7390.R` accordingly.
- Two datasets referenced in your original scripts were never uploaded
  and are not reproduced anywhere here: `LymphomaData.rda` and the
  `hnscc` / `nki70` / `metabric` exploratory sections of
  `CSIR-realdata.R`. Those sections of the original file also had
  unresolved bugs (undefined variables such as `cox_scores`,
  `common_patients`, `demo` used out of context) and were not carried
  over — ask if you want them cleaned up too.

## License

MIT (see `LICENSE`) — change this in `DESCRIPTION`/`LICENSE.md` if you'd
rather keep the code closed until the paper is published.
