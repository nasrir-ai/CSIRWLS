# csirwls — Censored Sliced Inverse Regression with Weighted Leverage Scores

R code and reproduction scripts for

> N. Asrir, T. Bouezmarni and A. Mkhadri, *Enhancing censored sliced inverse
> regression through weighted leverage scores*.

CSIR-WLS is a two-stage procedure for right-censored survival data with many
predictors. A censored weighted leverage score (CWLS) ranks the predictors,
the stopping rule G_c chooses how many to keep, and censored sliced inverse
regression (CSIR) estimates the effective dimension reduction directions from
the retained predictors.

The code depends on base R only (plus `parallel`). The real-data scripts also
use `survival` and `ipred` (CRAN) and, for TRANSBIG, the Bioconductor packages
`Biobase` and `breastCancerTRANSBIG`.

## Contents

```
R/                         package code: cwls(), bic_cut(), csir(), csir_wls(),
                           the competing methods (sir, cp_sir, csis, cris,
                           cmv_sis, cox_fit, cox_lasso_bic), evaluation criteria,
                           cross-validation (cv_survival)
reproduction/
  lib/                     shared infrastructure, designs and data-generating
                           processes, wrappers of the estimators
  run_main.R               Tables 1, 2, S1, S2, S6, S7 (and Section 5.4)
  run_dgp4.R               Tables 3, 4
  run_censoring.R          Table 5
  run_split.R              Table 6
  run_table7.R             Table 7
  run_table8.R             Table 8
  run_spiked.R             Tables S3, S4, S5, S11
  run_setting5.R           Tables S8, S9, S10, S12
  run_links.R              Table S13
  run_symmetric.R          Table S14
  run_ablation.R           Table S15
  run_timing.R             Tables S16, S17
  run_nonelliptical.R      Tables S18, S19
  run_defect.R             Table S20
  run_floor.R              Tables S21, S22
  make_tables.R            builds tables/<name>.tex and tables/<name>.csv
  run_all.R                runs every simulation, then make_tables.R and the applications
  applications/
    app_dlbcl.R            Section 6.1: Table 9, Figures 1-2
    app_pbc.R              Section 6.2: dimension test, Tables 10-11, Figure 3
    app_null_optimism.R    Section 6.3: optimism of the in-sample C-index without signal
    app_transbig.R         Section 6.3: Table 12, number of retained probes, stability
    run_applications.R     runs the four application scripts
tables/                    tables of the real-data applications (CSV and LaTeX)
figures/                   figures of the real-data applications (PDF and PNG)
```

## Running

Every command is run from the repository root.

```sh
# quick check of the whole pipeline (2 % of the replications, p = 80 only)
CSIRWLS_SCALE=0.02 CSIRWLS_P=80 Rscript reproduction/run_all.R

# full run (number of replications of the paper)
Rscript reproduction/run_all.R

# one table, then the LaTeX/CSV tables
Rscript reproduction/run_dgp4.R
Rscript reproduction/make_tables.R

# real-data applications
Rscript reproduction/applications/run_applications.R
```

Environment variables:

| variable | default | role |
|---|---|---|
| `CSIRWLS_SCALE` | 1 | fraction of the replications of the paper |
| `CSIRWLS_CORES` | all cores | cores used by `parallel::mclapply` (1 on Windows) |
| `CSIRWLS_P` | all | dimensions to run, e.g. `80,300` |
| `CSIRWLS_RESULTS`, `CSIRWLS_TABLES`, `CSIRWLS_FIGURES` | `results/`, `tables/`, `figures/` | output folders |
| `CSIRWLS_REDO` | 0 | 1 = recompute cells already stored |

Some drivers take an argument: `run_main.R dgp1` runs one data-generating
process, `run_setting5.R S12` runs Table S12 only and `run_floor.R S22` runs
Table S22 only.

**Resuming.** Each simulation cell is saved as soon as it finishes; an
interrupted run restarts where it stopped. The work can be split across
machines with `CSIRWLS_P` or the driver arguments, and the `results/` folders
copied together.

**Reproducibility.** Replication *r* of a cell always uses the same seed,
derived from the study, the data-generating process, the design, the
censoring mechanism, *p* and the censoring rate. Results therefore do not
depend on the number of cores or on the order in which cells are run.
Monte-Carlo results agree with the tables of the paper up to Monte-Carlo error.

## Settings

* n = 500; the active coefficients are drawn from U(0.5, 1.5) in every replication.
* H = 5 slices for every SIR-type method; projection dimension r = 2 and
  bandwidth rule (27) for the kernel weights.
* CWLS keeps the minimiser of G_c; CSIS, CRIS and cMV-SIS keep ⌊n / log n⌋
  predictors.
* SIR, CSIR, CP-SIR and the Cox model are fitted only when p < n; the
  Cox-lasso uses `cox_lasso_bic()` (BIC over 12 penalties).
* Censoring: log C = c + shift + N(0, 1), the constant c being calibrated once
  per cell to the target censoring rate; the achieved rates are written to
  `tables/censoring_rates.csv`.

## Real data

* DLBCL: data set `DLBCL` of the R package `ipred`.
* PBC: data set `pbc` of the R package `survival`.
* TRANSBIG: Bioconductor package `breastCancerTRANSBIG`
  (`BiocManager::install(c("Biobase", "breastCancerTRANSBIG"))`); alternatively
  set `CSIRWLS_TRANSBIG_RDS` to an `.rds` file holding `list(X, y, delta)`.

The application scripts write their tables to `tables/` and their figures to
`figures/`; the outputs included in this repository were produced with these
scripts.

## Computing cost

On one core with the reference BLAS, one replication takes about 1.5 s at
p = 80, 7 s at p = 300 and 10 s at p = 1500, and the whole set of simulations
of the order of 100 CPU-hours. A machine with 16–32 cores, or several machines
using `CSIRWLS_P`, brings this down to a few hours. The TRANSBIG application
takes a few hours on one core.
