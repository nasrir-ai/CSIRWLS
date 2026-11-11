## Section 6: runs the four application scripts. TRANSBIG is skipped (with a
## message) when breastCancerTRANSBIG / Biobase are not installed and
## CSIRWLS_TRANSBIG_RDS is not set.
scripts <- c("app_dlbcl.R", "app_pbc.R", "app_null_optimism.R", "app_transbig.R")
for (s in scripts) {
  if (s == "app_transbig.R" && !nzchar(Sys.getenv("CSIRWLS_TRANSBIG_RDS")) &&
      !(requireNamespace("breastCancerTRANSBIG", quietly = TRUE) && requireNamespace("Biobase", quietly = TRUE))) {
    message("skip app_transbig.R: BiocManager::install(c('Biobase', 'breastCancerTRANSBIG'))"); next
  }
  message("\n==================== ", s, " ====================")
  st <- system2(file.path(R.home("bin"), "Rscript"), file.path("reproduction", "applications", s))
  if (st != 0) warning(s, " exited with status ", st)
}
