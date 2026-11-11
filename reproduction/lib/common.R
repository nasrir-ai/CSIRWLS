## =============================================================================
## common.R -- infrastructure shared by every driver
##
##  * locates the repository and sources the package code in R/
##  * run-time configuration through environment variables
##      CSIRWLS_SCALE   fraction of the replications to run (default 1 = paper)
##      CSIRWLS_CORES   number of cores for parallel::mclapply (default: all)
##      CSIRWLS_RESULTS directory of the .rds results (default <repo>/results)
##      CSIRWLS_REDO    1 = recompute cells whose results already exist
##  * deterministic seeds: replication r of a cell always uses the same seed,
##    whatever the number of cores, the order in which cells are run or the
##    driver that runs them (this is what makes Table S22 a paired comparison
##    with Tables 1, 2, S1-S10)
##  * run_cell(): runs the replications of one cell in parallel and stores them
##    in results/<driver>/<key>.rds, so an interrupted run resumes where it
##    stopped.
## =============================================================================

find_repo_root <- function() {
  cand <- c(Sys.getenv("CSIRWLS_ROOT", ""), getwd(), file.path(getwd(), ".."),
            file.path(getwd(), "..", ".."))
  for (d in cand) {
    if (nzchar(d) && file.exists(file.path(d, "R", "cwls.R"))) return(normalizePath(d))
  }
  stop("cannot find the repository root (folder containing R/cwls.R); ",
       "run the scripts from the repository root or set CSIRWLS_ROOT")
}

REPO <- find_repo_root()
for (f in list.files(file.path(REPO, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
source(file.path(REPO, "reproduction", "lib", "designs.R"))
source(file.path(REPO, "reproduction", "lib", "methods.R"))
source(file.path(REPO, "reproduction", "lib", "cells.R"))

CFG <- list(
  scale   = as.numeric(Sys.getenv("CSIRWLS_SCALE", "1")),
  cores   = as.integer(Sys.getenv("CSIRWLS_CORES", as.character(max(1L, parallel::detectCores())))),
  results = Sys.getenv("CSIRWLS_RESULTS", file.path(REPO, "results")),
  redo    = Sys.getenv("CSIRWLS_REDO", "0") == "1"
)
if (.Platform$OS.type == "windows") CFG$cores <- 1L    # mclapply forks; not on Windows

## ---------------------------------------------------------------- replications
## number of replications actually run: N of the paper times CFG$scale
n_rep <- function(N) max(2L, as.integer(ceiling(N * CFG$scale)))
## Cox-lasso is run on the first NL replications of a cell (paper: 60/50/20 ...)
n_lasso <- function(NL) max(1L, as.integer(ceiling(NL * CFG$scale)))

## ----------------------------------------------------------------------- seeds
## A stable string hash: the seed of a cell depends only on its key.
hash_seed <- function(key) {
  h <- 0
  for (v in utf8ToInt(key)) h <- (h * 131 + v) %% 2147483629
  as.integer(h)
}
rep_seed <- function(key, r) as.integer((hash_seed(key) + 104729 * r) %% 2147483629)

## ----------------------------------------------------------------- cell runner
## fun(r) must start by generating the data (it is called right after
## set.seed(rep_seed(key, r))) and returns a named list or a data frame row.
run_cell <- function(driver, key, N, fun, meta = list(), cores = CFG$cores) {
  dir <- file.path(CFG$results, driver)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file <- file.path(dir, paste0(gsub("[^A-Za-z0-9._=-]", "_", key), ".rds"))
  if (file.exists(file) && !CFG$redo) {
    old <- readRDS(file)
    if (old$N >= N) { message(sprintf("[%s] %s: already done (N = %d)", driver, key, old$N)); return(invisible(old)) }
  }
  t0 <- proc.time()[["elapsed"]]
  one <- function(r) {
    set.seed(rep_seed(key, r))
    out <- tryCatch(fun(r), error = function(e) structure(list(msg = conditionMessage(e)), class = "rep_error"))
    out
  }
  res <- if (cores > 1L) parallel::mclapply(seq_len(N), one, mc.cores = cores, mc.preschedule = FALSE)
         else lapply(seq_len(N), one)
  bad <- vapply(res, function(z) inherits(z, "rep_error") || inherits(z, "try-error"), logical(1))
  if (any(bad)) {
    msgs <- unique(vapply(res[bad], function(z) if (is.list(z)) z$msg else as.character(z), character(1)))
    warning(sprintf("[%s] %s: %d replication(s) failed: %s", driver, key, sum(bad), paste(msgs, collapse = " | ")))
  }
  reps <- rbind_fill(lapply(which(!bad), function(i) {
    z <- as.data.frame(res[[i]], stringsAsFactors = FALSE, check.names = FALSE); z$rep <- i
    for (m in names(meta)) z[[m]] <- meta[[m]]
    z
  }))
  out <- list(driver = driver, key = key, N = N, n_failed = sum(bad), reps = reps,
              seconds = proc.time()[["elapsed"]] - t0, date = Sys.time(),
              scale = CFG$scale)
  saveRDS(out, file)
  message(sprintf("[%s] %s: %d replications in %.0f s", driver, key, N, out$seconds))
  invisible(out)
}

## rbind data frames whose columns differ (missing columns filled with NA)
rbind_fill <- function(lst) {
  lst <- Filter(function(d) !is.null(d) && nrow(d) > 0, lst)
  if (!length(lst)) return(NULL)
  cols <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(d) {
    for (cc in setdiff(cols, names(d))) d[[cc]] <- NA
    d[cols]
  }))
}

## read every stored cell of a driver into one data frame (cell descriptors are
## stored as columns of 'reps' by the drivers)
load_driver <- function(driver) {
  dir <- file.path(CFG$results, driver)
  fs <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  if (!length(fs)) return(NULL)
  rbind_fill(lapply(fs, function(f) readRDS(f)$reps))
}

## small helpers used by the drivers
mean_se <- function(v) { v <- v[is.finite(v)]; c(mean = mean(v), se = if (length(v) > 1) stats::sd(v) / sqrt(length(v)) else NA_real_) }
cat_cell <- function(...) message(sprintf(...))
