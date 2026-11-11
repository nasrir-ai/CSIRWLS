## Runs every driver, then builds the tables. Long: see README for the cost.
## Quick check of the whole pipeline (2-5 % of the replications):
##   CSIRWLS_SCALE=0.02 Rscript reproduction/run_all.R
drivers <- c("run_main.R", "run_split.R", "run_dgp4.R", "run_censoring.R",
             "run_spiked.R", "run_setting5.R", "run_links.R", "run_symmetric.R", "run_ablation.R",
             "run_nonelliptical.R", "run_table7.R", "run_table8.R", "run_defect.R", "run_floor.R",
             "run_timing.R", "make_tables.R", "applications/run_applications.R")
for (d in drivers) {
  message("\n==================== ", d, " ====================")
  st <- system2(file.path(R.home("bin"), "Rscript"), file.path("reproduction", d))
  if (st != 0) warning(d, " exited with status ", st)
}
