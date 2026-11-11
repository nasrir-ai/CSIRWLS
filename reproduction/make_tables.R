## =============================================================================
## make_tables.R -- builds every table from the stored replications
##
## For each table: tables/<name>.tex (LaTeX, booktabs) and tables/<name>.csv.
## Tables whose drivers have not been run are skipped with a message.
## Also writes tables/censoring_rates.csv (achieved censoring rate per cell).
## Usage: Rscript reproduction/make_tables.R
## =============================================================================
source(file.path("reproduction", "lib", "common.R"))
OUT <- Sys.getenv("CSIRWLS_TABLES", file.path(REPO, "tables"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------------ formatting
f_num <- function(x, k) if (is.na(x)) "--" else formatC(x, format = "f", digits = k)
f_ms  <- function(v, k = 3) {                     # "mean (se)"
  v <- v[is.finite(v)]
  if (!length(v)) return("$\\dagger$")
  se <- if (length(v) > 1) stats::sd(v) / sqrt(length(v)) else NA
  sprintf("%s (%s)", f_num(mean(v), k), f_num(se, k))
}
q1 <- function(v, pr) as.numeric(stats::quantile(v, pr, type = 1, names = FALSE, na.rm = TRUE))
f_mms <- function(v) { v <- v[is.finite(v)]; if (!length(v)) return("$\\dagger$"); sprintf("%d [%d]", q1(v, .5), q1(v, .95)) }
f_med <- function(v) { v <- v[is.finite(v)]; if (!length(v)) return("--"); sprintf("%d", q1(v, .5)) }
col_or_na <- function(d, nm) if (nm %in% names(d)) d[[nm]] else rep(NA_real_, nrow(d))
pct <- function(x) sprintf("%d%%", round(100 * x))

## write a character data frame as LaTeX + CSV; 'group' columns are blanked
## when equal to the previous row (as in the paper's layout)
write_table <- function(tab, name, caption, group = character(0)) {
  csv <- tab; names(csv)[names(csv) == "pp"] <- "p"
  utils::write.csv(csv, file.path(OUT, paste0(name, ".csv")), row.names = FALSE)
  shown <- tab
  if (nrow(tab) > 1) {
    same <- rep(TRUE, nrow(tab) - 1)
    for (g in group) {            # blank a group label repeated from the row above
      same <- same & (tab[[g]][-1] == tab[[g]][-nrow(tab)])
      shown[[g]][-1][same] <- ""
    }
  }
  names(tab)[names(tab) == "pp"] <- "p"; names(shown) <- names(tab)
  esc <- function(s) gsub("%", "\\\\%", gsub("_", "\\\\_", s))
  esc_mixed <- function(s) {          # escape outside $...$ only
    parts <- strsplit(s, "$", fixed = TRUE)[[1]]
    odd <- seq_along(parts) %% 2 == 1
    parts[odd] <- esc(parts[odd]); paste(parts, collapse = "$")
  }
  hdr <- paste(gsub("|", "$|$", esc(names(tab)), fixed = TRUE), collapse = " & ")
  body <- apply(shown, 1, function(r) paste(ifelse(grepl("\\$", r), r, esc(r)), collapse = " & "))
  tex <- c("\\begin{table}[ht]", "\\centering", "\\scriptsize",
           sprintf("\\caption{%s}", esc_mixed(caption)), sprintf("\\label{tab:%s}", name),
           sprintf("\\begin{tabular}{%s}", paste(rep("l", ncol(tab)), collapse = "")),
           "\\toprule", paste0(hdr, " \\\\"), "\\midrule",
           paste0(body, " \\\\"), "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(tex, file.path(OUT, paste0(name, ".tex")))
  message("wrote ", name, " (", nrow(tab), " rows)")
}

## loop over the cells of a data frame in a fixed order, building one row each
by_cells <- function(d, by, row_fun) {
  if (is.null(d) || !nrow(d)) return(NULL)
  keys <- unique(d[by]); keys <- keys[do.call(order, unname(as.list(keys))), , drop = FALSE]
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    sel <- Reduce(`&`, lapply(by, function(b) d[[b]] == keys[[b]][i]))
    sub <- d[sel, , drop = FALSE]
    c(lapply(keys[i, , drop = FALSE], as.character), row_fun(sub))
  })
  as.data.frame(do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE, check.names = FALSE))),
                stringsAsFactors = FALSE, check.names = FALSE)
}
rho_of <- function(lbl) sub("AR1\\((.*)\\)", "\\1", lbl)
skip <- function(what) message("skip ", what, ": no results yet")

## ---------------------------------------------------- screening tables (1, S1-S5)
screen_row <- function(s) list(
  TPR = f_ms(s$cwls_tpr, 2), FPR = f_ms(s$cwls_fpr, 3), FDR = f_ms(s$cwls_fdr, 2),
  `|S|` = f_ms(s$cwls_size, 1), MMS = f_mms(s$mms_CWLS),
  `CSIS TPR` = f_ms(col_or_na(s, "tpr_CSIS"), 2), `CSIS MMS` = f_mms(col_or_na(s, "mms_CSIS")),
  `CRIS TPR` = f_ms(col_or_na(s, "tpr_CRIS"), 2), `CRIS MMS` = f_mms(col_or_na(s, "mms_CRIS")),
  `cMV-SIS TPR` = f_ms(col_or_na(s, "tpr_cMV-SIS"), 2), `cMV-SIS MMS` = f_mms(col_or_na(s, "mms_cMV-SIS")))
edr_row <- function(s) list(
  SIR = f_ms(col_or_na(s, "pc_SIR")), CSIR = f_ms(col_or_na(s, "pc_CSIR")),
  `CP-SIR` = f_ms(col_or_na(s, "pc_CP-SIR")), Cox = f_ms(col_or_na(s, "pc_Cox")),
  `Cox-lasso` = f_ms(col_or_na(s, "pc_Cox-lasso")), `CSIR-WLS` = f_ms(s[["pc_CSIR-WLS"]]),
  `|A(c)|` = f_ms(s$cwls_size, 1))

main <- load_driver("main")
if (!is.null(main)) {
  main$Setting <- ifelse(rho_of(main$design) == "0.5", "(1): rho=0.5", "(2): rho=0.8")
  main$pp <- sprintf("%04d", main$p); main$CR <- pct(main$target)
  for (tb in list(c("dgp1", "table1", "table2", "DGP 5.1"), c("dgp2", "tableS1", "tableS6", "DGP 5.2"),
                  c("dgp3", "tableS2", "tableS7", "DGP 5.3"))) {
    d <- main[main$model == tb[1], ]
    if (!nrow(d)) next
    t1 <- by_cells(d, c("Setting", "pp", "CR"), screen_row); t1$pp <- as.integer(t1$pp)
    write_table(t1, tb[2], sprintf("Screening performance under %s, Settings (1)-(2), n = 500. CWLS with G_c; CSIS, CRIS, cMV-SIS keep floor(n/log n) = 80. MMS: median [95%% quantile].", tb[4]), c("Setting", "pp"))
    t2 <- by_cells(d, c("Setting", "pp", "CR"), edr_row); t2$pp <- as.integer(t2$pp)
    write_table(t2, tb[3], sprintf("Proximity criterion under %s, Settings (1)-(2), n = 500 (Monte-Carlo s.e.). $\\dagger$: not applicable (p > n).", tb[4]), c("Setting", "pp"))
  }
  ## Section 5.4: quantiles of the minimum model size, DGP 5.1, 25 %
  d <- main[main$model == "dgp1" & main$target == 0.25, ]
  if (nrow(d)) {
    qs <- function(v) { v <- v[is.finite(v)]; if (!length(v)) "--" else paste(q1(v, c(.05, .5, .95)), collapse = "/") }
    t54 <- by_cells(d, c("Setting", "pp"), function(s) list(
      CWLS = qs(s$mms_CWLS), CSIS = qs(s$mms_CSIS), CRIS = qs(s$mms_CRIS),
      `cMV-SIS` = qs(s[["mms_cMV-SIS"]]), `Cox-lasso` = qs(col_or_na(s, "mms_lasso"))))
    t54$pp <- as.integer(t54$pp)
    write_table(t54, "section5_4_mms", "Section 5.4: minimum model size, 5%/50%/95% quantiles, DGP 5.1, 25% censoring (Cox-lasso on its replications).", "Setting")
  }
} else skip("Tables 1, 2, S1, S2, S6, S7")

## ------------------------------------------------------------ spiked (S3-S5, S11)
sp <- load_driver("spiked")
if (!is.null(sp)) {
  sp$Setting <- sub("spiked", "(", sp$design); sp$Setting <- paste0(sp$Setting, ")")
  sp$pp <- sprintf("%04d", sp$p); sp$CR <- pct(sp$target)
  for (tb in list(c("dgp1", "tableS3", "DGP 5.1"), c("dgp2", "tableS4", "DGP 5.2"), c("dgp3", "tableS5", "DGP 5.3"))) {
    d <- sp[sp$model == tb[1], ]
    if (!nrow(d)) next
    t <- by_cells(d, c("Setting", "pp", "CR"), screen_row); t$pp <- as.integer(t$pp)
    write_table(t, tb[2], sprintf("Screening under %s, spiked Settings (3) (|S| = 6) and (4) (|S| = 9), n = 500.", tb[3]), c("Setting", "pp"))
  }
  for (mdl in unique(sp$model)) {
    d <- sp[sp$model == mdl & sp$target == 0.25, ]
    if (!nrow(d)) next
    rows <- list()
    for (st in sort(unique(d$Setting))) for (p in sort(unique(d$p))) {
      s <- d[d$Setting == st & d$p == p, ]
      for (m in c("SIR", "CSIR", "CP-SIR", "Cox", "Cox-lasso", "CSIR-WLS"))
        rows[[length(rows) + 1]] <- data.frame(Setting = st, p = as.character(p), Method = m,
          PC = f_ms(col_or_na(s, paste0("pc_", m))), IC = f_ms(col_or_na(s, paste0("ic_", m))),
          `|A(c)|` = if (m == "CSIR-WLS") f_ms(s$cwls_size, 1) else "",
          check.names = FALSE, stringsAsFactors = FALSE)
    }
    nm <- if (mdl == "dgp1") "tableS11" else paste0("tableS11_", mdl)
    write_table(do.call(rbind, rows), nm, sprintf("EDR estimation under the spiked Settings (3)-(4), %s, 25%% censoring: PC and index correlation IC.", toupper(mdl)), c("Setting", "p"))
  }
} else skip("Tables S3-S5, S11")

## ----------------------------------------------------------- Setting (5) (S8-S10, S12)
cfg_name <- c("block(0.5,0.2,0.4)" = "(A1)", "block(0.8,0.2,0.4)" = "(A2)", "block(0.8,0.2,0.6)" = "(A3)")
s5 <- load_driver("setting5")
if (!is.null(s5)) {
  s5$Config <- unname(cfg_name[s5$design]); s5$pp <- sprintf("%04d", s5$p); s5$CR <- pct(s5$target)
  for (tb in list(c("dgp1", "tableS8", "DGP 5.1"), c("dgp2", "tableS9", "DGP 5.2"), c("dgp3", "tableS10", "DGP 5.3"))) {
    d <- s5[s5$model == tb[1], ]
    if (!nrow(d)) next
    t <- by_cells(d, c("Config", "pp", "CR"), edr_row); t$pp <- as.integer(t$pp)
    write_table(t, tb[2], sprintf("Proximity criterion under %s, Setting (5), configurations (A1)-(A3).", tb[3]), c("Config", "pp"))
    t <- by_cells(d, c("Config", "pp", "CR"), screen_row); t$pp <- as.integer(t$pp)
    write_table(t, paste0(tb[2], "_screening"), sprintf("Screening under %s, Setting (5) (not in the paper; code output).", tb[3]), c("Config", "pp"))
  }
} else skip("Tables S8-S10")
s12 <- load_driver("setting5_more")
if (!is.null(s12)) {
  s12$alpha <- sub("block", "", s12$design); s12$pp <- sprintf("%04d", s12$p); s12$CR <- pct(s12$target)
  t <- by_cells(s12, c("alpha", "pp", "CR"), edr_row); t$pp <- as.integer(t$pp)
  write_table(t, "tableS12", "Proximity criterion under Setting (5), four further configurations, DGP 5.1.", c("alpha", "pp"))
  t <- by_cells(s12, c("alpha", "pp", "CR"), screen_row); t$pp <- as.integer(t$pp)
  write_table(t, "tableS12_screening", "Screening, further Setting (5) configurations (Remark S2.1).", c("alpha", "pp"))
} else skip("Table S12")

## ------------------------------------------------------------- DGP 5.4 (Tables 3, 4)
g4 <- load_driver("dgp4")
if (!is.null(g4)) g4 <- g4[grepl(sprintf("^dgp4-sigma%g\\|", DGP4_SIGMA), g4$key), , drop = FALSE]
if (!is.null(g4) && nrow(g4)) {
  g4$rho <- rho_of(g4$design)
  rows <- list()
  for (r in sort(unique(g4$rho))) for (p in sort(unique(g4$p))) {
    s <- g4[g4$rho == r & g4$p == p, ]
    for (m in c("CSIR", "CP-SIR", "Cox", "Cox-lasso", "CP-SIRWLS", "CSIR-WLS"))
      rows[[length(rows) + 1]] <- data.frame(rho = r, p = as.character(p), Method = m,
        R = f_ms(col_or_na(s, paste0("R_", m))), V = f_ms(col_or_na(s, paste0("V_", m))),
        `|A(c)|` = if (m == "CSIR-WLS") f_ms(s$cwls_size, 1) else "", check.names = FALSE, stringsAsFactors = FALSE)
  }
  write_table(do.call(rbind, rows), "table3", "Two-dimensional EDR space, DGP 5.4 (d = 2), n = 500, 25% censoring: recovery R and validity V.", c("rho", "p"))
  g4$pp <- sprintf("%04d", g4$p)
  t4 <- by_cells(g4, c("rho", "pp"), function(s) list(CWLS = f_mms(s$mms_CWLS), CSIS = f_mms(s$mms_CSIS),
                                                       `Cox-lasso` = f_mms(col_or_na(s, "mms_lasso"))))
  t4$pp <- as.integer(t4$pp)
  write_table(t4, "table4", "Minimum model size under DGP 5.4: median [95% quantile].", "rho")
} else skip("Tables 3-4")

## ------------------------------------------------------------- Table 5 (censoring)
ce <- load_driver("censoring")
if (!is.null(ce)) {
  ce$Mechanism <- c(indep = "1 Independent", edr = "2 EDR-dependent", misspec = "3 Misspecified")[ce$mech]
  ce$pp <- sprintf("%04d", ce$p)
  t <- by_cells(ce, c("pp", "Mechanism"), function(s) list(
    `Median MMS (CWLS)` = f_med(s$mms_CWLS), `Median MMS (CSIS)` = f_med(s$mms_CSIS),
    `PC (CSIR-WLS)` = f_ms(s[["pc_CSIR-WLS"]]), `|A(c)|` = f_ms(s$cwls_size, 1),
    `censoring` = f_num(mean(s$cens_rate), 3)))
  t$pp <- as.integer(t$pp); t$Mechanism <- sub("^[0-9] ", "", t$Mechanism)
  write_table(t, "table5", "Effect of the censoring mechanism, DGP 5.1, Setting (1), 25% censoring.", "pp")
} else skip("Table 5")

## ---------------------------------------------------------------- Table 6 (split)
spl <- load_driver("split")
if (!is.null(spl)) {
  spl$pp <- sprintf("%04d", spl$p)
  t <- by_cells(spl, "pp", function(s) list(
    `Full sample` = f_ms(s[["pc_CSIR-WLS"]]), `Split sample` = f_ms(s$pc_split),
    `Cross-fitted` = f_ms(s$pc_cross), `|A(c)| full` = f_num(mean(s$cwls_size), 1),
    `|A(c)| split` = f_num(mean(s$size_split), 1), `Coverage (split)` = f_num(mean(s$cover_split), 3)))
  t$pp <- as.integer(t$pp)
  write_table(t, "table6", "Full-sample, split-sample and cross-fitted CSIR-WLS, DGP 5.1, Setting (1), 25% censoring.")
} else skip("Table 6")

## ---------------------------------------------------------------- Table S13 (links)
lk <- load_driver("links")
if (!is.null(lk)) {
  lk$Link <- c(saturating = "1 Saturating", threshold = "2 Threshold", heteroscedastic = "3 Heteroscedastic")[lk$model]
  lk$pp <- sprintf("%04d", lk$p)
  sq <- function(s, m) f_ms(col_or_na(s, paste0("pc_", m))^2)
  t <- by_cells(lk, c("Link", "pp"), function(s) list(CSIR = sq(s, "CSIR"), `CP-SIR` = sq(s, "CP-SIR"),
    Cox = sq(s, "Cox"), `Cox-lasso` = sq(s, "Cox-lasso"), `CSIR-WLS` = sq(s, "CSIR-WLS")))
  t$pp <- as.integer(t$pp); t$Link <- sub("^[0-9] ", "", t$Link)
  write_table(t, "tableS13", "Single-index links away from proportional hazards, Setting (1), 25% censoring: squared proximity criterion.", "Link")
} else skip("Table S13")

## ----------------------------------------------------------- Table S14 (symmetric)
sy <- load_driver("symmetric")
if (!is.null(sy)) {
  s <- sy
  t <- data.frame(Method = c("CWLS / CSIR-WLS", "CSIS", "CSIR", "Cox", "Cox-lasso"),
                  MMS = c(f_mms(s$mms_CWLS), f_mms(col_or_na(s, "mms_CSIS")), "--", "--", f_mms(col_or_na(s, "mms_lasso"))),
                  PC = c(f_ms(s[["pc_CSIR-WLS"]]), "--", f_ms(col_or_na(s, "pc_CSIR")), f_ms(col_or_na(s, "pc_Cox")),
                         f_ms(col_or_na(s, "pc_Cox-lasso"))), check.names = FALSE)
  write_table(t, "tableS14", sprintf("Symmetric link, Setting (1), p = 80, 25%% censoring (achieved %s). MMS: median [95%% quantile]; PC: mean (s.e.). Retained by G_c: %s.",
                                     pct(mean(s$cens_rate)), f_ms(s$cwls_size, 1)))
} else skip("Table S14")

## ------------------------------------------------------------ Table S15 (ablation)
ab <- load_driver("ablation")
if (!is.null(ab)) {
  ab$pp <- sprintf("%04d", ab$p); ab$tg <- sprintf("%.2f", ab$target)
  mm <- function(v) sprintf("%d (%d)", q1(v, .5), q1(v, .95))
  t <- by_cells(ab, c("pp", "tg"), function(s) list(`C.R.` = pct(mean(s$cens_rate)),
    Proposed = mm(s$mms_CWLS), Naive = mm(s$mms_naive), `Complete-case` = mm(s$mms_cc)))
  t$pp <- as.integer(t$pp); t$tg <- NULL
  write_table(t, "tableS14", "Ablation of the censoring correction: minimum model size, median (95th percentile). DGP 5.1, Setting (1).", "pp")
} else skip("Table S15")

## -------------------------------------------------------- Tables S16, S17 (timing)
tmg <- load_driver("timing")
if (!is.null(tmg)) {
  tmg$pp <- sprintf("%04d", tmg$p); tmg$CR <- pct(tmg$target)
  av <- function(v) { v <- v[is.finite(v)]; if (!length(v)) "$\\dagger$" else f_num(mean(v), 2) }
  t <- by_cells(tmg, c("pp", "CR"), function(s) list(CWLS = av(s$t_CWLS), CSIS = av(s$t_CSIS)))
  t$pp <- as.integer(t$pp)
  write_table(t, "tableS15", "Average computational time (seconds) for variable screening, DGP 5.3, Setting (1).", "pp")
  t <- by_cells(tmg, c("pp", "CR"), function(s) list(SIR = av(s$t_SIR), CSIR = av(col_or_na(s, "t_CSIR")),
    COX = av(col_or_na(s, "t_Cox")), `CSIR-WLS` = av(s$t_CSIRWLS)))
  t$pp <- as.integer(t$pp)
  write_table(t, "tableS16", "Average computational time (seconds) for EDR estimation, DGP 5.3, Setting (1).", "pp")
} else skip("Tables S16-S17")

## ------------------------------------------------------------ Table 7 (Section 5.8.2)
## separation of the departures from the theory: shared cells of Tables 1-2 and
## S18-S19 (reference, Azzalini designs) + the two controls of run_table7.R
t7 <- list(); for (drv in c("main", "nonelliptical", "table7")) t7[[drv]] <- load_driver(drv)
if (!is.null(t7$table7)) {
  f_mms2 <- function(v) { v <- v[is.finite(v)]; sprintf("%d/%d", q1(v, .5), q1(v, .95)) }
  rows <- list()
  for (p in c(80, 300)) for (r in table7_rows()) {
    d <- t7[[r$driver]]
    if (is.null(d)) next
    d <- d[d$key == cell_key(r$driver, "dgp1", r$spec, "indep", p, 0.25), , drop = FALSE]
    if (!nrow(d)) next
    rows[[length(rows) + 1L]] <- data.frame(
      p = p, design = r$name, `(C.1)` = r$c1,
      `lambda_max` = sprintf("%.1f", lambda_max(r$spec, N_DESIGN, p)),
      `|A|` = f_ms(d$cwls_size, 1), MMS = f_mms2(d$mms_CWLS),
      `PC CSIR-WLS` = f_ms(d[["pc_CSIR-WLS"]]), `PC Cox-lasso` = f_ms(d[["pc_Cox-lasso"]]),
      `N` = nrow(d), check.names = FALSE, stringsAsFactors = FALSE)
  }
  write_table(do.call(rbind, rows), "table7",
    "Separating the departures from the theory in the non-elliptical designs. DGP 5.1, n = 500, 25% censoring. |A| = predictors retained by G_c, MMS = minimum model size (median/95% quantile), PC = proximity criterion (Monte-Carlo s.e.). The TPR of CWLS is not shown.", "p")
} else skip("Table 7")

## ------------------------------------------------------------ Table 8 (Section 5.8.3)
sg <- load_driver("signs")
if (!is.null(sg)) {
  cfg_of <- c("0.5,0.2,0.4" = "(A1)", "0.8,0.2,0.4" = "(A2)", "0.8,0.2,0.6" = "(A3)", "0.2,0.2,0.2" = "(A4)")
  alpha_of <- sub("^block\\((.*)\\).*$", "\\1", sg$design)
  sg$Config <- cfg_of[alpha_of]
  sg$Signs  <- ifelse(grepl("-neg", sg$design), "two negative", "all positive")
  sg$pp <- sprintf("%04d", sg$p)
  f_marg <- function(v) { v <- v[is.finite(v)]; sprintf("%s / %.2f", format(signif(stats::median(v), 3)), min(v)) }
  t <- by_cells(sg, c("pp", "Config", "Signs"), function(s) list(
    `m_1/2` = f_marg(s$m_half), `m_1` = f_marg(s$m_one),
    `CWLS MMS` = f_mms(s$mms_CWLS), `CWLS TPR` = f_ms(s$cwls_tpr, 2),
    `CSIS MMS` = f_mms(s$mms_CSIS), `CSIS TPR` = f_ms(s$tpr_CSIS, 2),
    `CRIS MMS` = f_mms(s$mms_CRIS), `CRIS TPR` = f_ms(s$tpr_CRIS, 2),
    `cMV-SIS MMS` = f_mms(s[["mms_cMV-SIS"]]), `cMV-SIS TPR` = f_ms(s[["tpr_cMV-SIS"]], 2),
    N = nrow(s)))
  t$pp <- as.integer(t$pp)
  write_table(t, "table8",
    "Sign cancellation. Setting (5), DGP 5.1, n = 500, 25% censoring, |S| = 5. Magnitudes U(0.5, 1.5); all positive or two negative signs at random positions. m_1/2 and m_1: population ranking margins for (Sigma^{1/2} beta) and (Sigma beta), median / minimum over the replications; MMS: median [95% quantile]; TPR (Monte-Carlo s.e.) at the stopping rules of Table 1.", c("pp", "Config"))
} else skip("Table 8")

## ------------------------------------------------- Tables S18, S19 (non-elliptical)
ne <- load_driver("nonelliptical")
if (!is.null(ne)) {
  g <- if (!is.null(main)) main[main$model == "dgp1" & main$design == "AR1(0.5)", ] else NULL
  if (!is.null(g)) g <- g[intersect(names(g), names(ne))]
  all <- rbind_fill(list(g, ne))
  all$Distribution <- c("AR1(0.5)" = "1 AR(1) Gaussian", skewnormal = "2 Skew-normal", skewt = "3 Skew-t (nu = 5)")[all$design]
  all$pp <- sprintf("%04d", all$p); all$CR <- pct(all$target)
  keep <- as.integer(N_DESIGN / log(N_DESIGN))
  d <- all[all$target == 0.25, ]
  t <- by_cells(d, c("Distribution", "pp"), function(s) list(
    `CWLS TPR` = f_ms(s$cwls_tpr, 2), `CWLS FPR` = f_ms(s$cwls_fpr, 3), `CWLS |S|` = f_ms(s$cwls_size, 1),
    `CWLS MMS` = f_mms(s$mms_CWLS), `CSIS TPR` = f_ms(s$tpr_CSIS, 2),
    `CSIS FPR` = f_ms((keep - s$tpr_CSIS * 5) / (s$p - 5), 3), `CSIS |S|` = f_ms(rep(keep, nrow(s)), 1),
    `CSIS MMS` = f_mms(s$mms_CSIS)))
  t$pp <- as.integer(t$pp); t$Distribution <- sub("^[0-9] ", "", t$Distribution)
  write_table(t, "tableS17", "Variable screening under non-elliptical designs, DGP 5.1, rho = 0.5, 25% censoring. Gaussian rows: Table 1.", "Distribution")
  t <- by_cells(all, c("Distribution", "pp", "CR"), edr_row)
  t$pp <- as.integer(t$pp); t$Distribution <- sub("^[0-9] ", "", t$Distribution)
  write_table(t, "tableS18", "Proximity criterion under non-elliptical designs, DGP 5.1, rho = 0.5. Gaussian rows: Table 2.", c("Distribution", "pp"))
} else skip("Tables S18-S19")

## ------------------------------------------------------------ Table S20 (defect)
dfile <- file.path(CFG$results, "defect", "defect.rds")
if (file.exists(dfile)) {
  d <- readRDS(dfile)$reps
  t <- data.frame(design = d$design, p = d$p, `xi_rel` = sprintf("%.3f", d$xi_rel),
                  `|c|` = sprintf("%.4f", d$c_abs), `sd(nu)` = sprintf("%.4f", d$sd_nu),
                  `Delta_lin` = sprintf("%.4f", d$Delta_lin), Delta = sprintf("%.4f", d$Delta),
                  `lambda_max` = sprintf("%.1f", d$lambda_max), check.names = FALSE, stringsAsFactors = FALSE)
  write_table(t, "tableS19", sprintf("Population linearity defect, DGP 5.1, AR(1) rho = 0.5, H = 5 (beta = %s, whitening %s).", d$beta[1], d$whiten[1]), "design")
} else skip("Table S20")

## -------------------------------------------------------- Tables S21, S22 (floor)
f20 <- load_driver("floor_S21")
if (!is.null(f20)) {
  f20$Design <- c("AR1(0.5)" = "1 Gaussian", skewnormal = "2 skew-normal", skewt = "3 skew-t5")[f20$design]
  f20$pp <- sprintf("%04d", f20$p)
  rows <- list()
  for (dz in sort(unique(f20$Design))) for (p in sort(unique(f20[f20$Design == dz, "p"]))) {
    s <- f20[f20$Design == dz & f20$p == p, ]
    for (m in c("m0", "m01", "m1"))
      rows[[length(rows) + 1]] <- data.frame(design = sub("^[0-9] ", "", dz), p = as.character(p),
        mult = c(m0 = "0.0", m01 = "0.1", m1 = "1.0")[[m]], `|A|` = f_ms(s[[paste0("size_", m)]], 1),
        TPR = f_num(mean(s[[paste0("tpr_", m)]]), 2),
        MMS = sprintf("%d/%d", q1(s[[paste0("mms_", m)]], .5), q1(s[[paste0("mms_", m)]], .95)),
        `PC CSIR-WLS` = f_ms(s[[paste0("pc_", m)]], 4), check.names = FALSE, stringsAsFactors = FALSE)
  }
  write_table(do.call(rbind, rows), "tableS20", "Sensitivity to the floor of (C.9)(c), DGP 5.1, 25% censoring.", c("design", "p"))
} else skip("Table S21")
f21 <- load_driver("floor_S22")
if (!is.null(f21)) {
  f21$DGP <- sub("dgp", "", f21$model); f21$rho <- rho_of(f21$design)
  f21$pp <- sprintf("%04d", f21$p); f21$CR <- pct(f21$target)
  t <- by_cells(f21, c("pp", "DGP", "rho", "CR"), function(s) list(
    `|A| none` = f_num(mean(s$size_m0), 1), `|A| b_n` = f_num(mean(s$size_m1), 1),
    `MMS none` = f_med(s$mms_m0), `MMS b_n` = f_med(s$mms_m1),
    `PC none` = f_ms(s$pc_m0), `PC b_n` = f_ms(s$pc_m1)))
  t$pp <- as.integer(t$pp)
  t <- t[c("DGP", "rho", "pp", "CR", setdiff(names(t), c("DGP", "rho", "pp", "CR")))]
  write_table(t, "tableS21", "Main cells re-run with the floor of (C.9)(c) (same seeds as Tables 1, 2, S1-S10).")
} else skip("Table S22")

## ----------------------------------------------------- achieved censoring rates
cr <- list()
for (drv in c("main", "spiked", "setting5", "setting5_more", "dgp4", "censoring", "split",
              "links", "symmetric", "ablation", "nonelliptical", "floor_S21", "timing")) {
  d <- load_driver(drv); if (is.null(d)) next
  d$cell <- paste(d$model, d$design, d$mech, d$p, d$target)
  a <- stats::aggregate(cens_rate ~ cell, data = d, FUN = mean)
  a$driver <- drv; a$n_rep <- as.vector(table(d$cell)[a$cell])
  cr[[drv]] <- a
}
if (length(cr)) {
  cr <- do.call(rbind, cr)
  utils::write.csv(cr[c("driver", "cell", "n_rep", "cens_rate")], file.path(OUT, "censoring_rates.csv"), row.names = FALSE)
  message("wrote censoring_rates.csv")
}
