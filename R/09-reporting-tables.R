#' Mise en forme des resultats (tables LaTeX)
#'
#' @description
#' Fonctions qui transforment la sortie des boucles de simulation
#' (\code{08-simulation-runners.R}) en tables \code{xtable} pretes a
#' coller dans le papier LaTeX. Comme pour les "runners", plusieurs
#' \code{latex_table()} portant le meme nom existaient dans vos scripts
#' avec des colonnes differentes ; elles sont renommees explicitement
#' ici : \code{\link{latex_table_edr}} (PC/Cor/Frob, pour
#' \code{\link{run_simulation_edr}}), \code{\link{latex_table_screening}}
#' (TPR/FPR/FDR, pour \code{\link{run_simulation_screening}}).
#'
#' @name reporting-tables
NULL

#' Table LaTeX : comparaison EDR a index simple (PC / Cor / Frob)
#'
#' @param results Sortie \code{$results} de \code{\link{run_simulation_edr}}.
#' @param censor_rates Vecteur des taux de censure a inclure.
#' @return Invisible ; imprime la table LaTeX (via \code{xtable}) sur la
#'   console.
#' @export
latex_table_edr <- function(results, censor_rates) {
  if (!requireNamespace("xtable", quietly = TRUE)) {
    stop("Le package 'xtable' est necessaire pour latex_table_edr().")
  }
  method_names <- c("SIR", "CSIR", "COX", "CSIR-WLS")
  combined_results <- data.frame()

  for (censor_rate in censor_rates) {
    result <- results[[paste0("censor_rate_", censor_rate)]]
    result_df <- data.frame(
      CensorRate = paste0(100 * censor_rate, "%"),
      Method = method_names,
      PC = result$PC,
      Cor = result$COR,
      Frob = result$Frob,
      Time_s = result$Time
    )
    combined_results <- rbind(combined_results, result_df)
  }

  print(
    xtable::xtable(combined_results,
      caption = "Results for different censoring rates",
      label = "tab:combined_censor_rates"
    ),
    type = "latex",
    include.rownames = FALSE,
    caption.placement = "top",
    table.placement = "H"
  )
}

#' Table LaTeX : comparaison EDR a index multiple (TC / SD / Cor)
#'
#' @param results Sortie \code{$results} de \code{\link{run_simulation_multi}}.
#' @param censor_rates Vecteur des taux de censure a inclure.
#' @param d_true Dimension EDR vraie (utilisee dans la legende).
#' @return Invisible ; imprime la table LaTeX sur la console.
#' @export
latex_table_multi <- function(results, censor_rates, d_true = NULL) {
  if (!requireNamespace("xtable", quietly = TRUE)) {
    stop("Le package 'xtable' est necessaire pour latex_table_multi().")
  }
  method_names <- c("SIR", "CSIR", "COX", "CSIR-WLS")
  combined_results <- data.frame()

  for (censor_rate in censor_rates) {
    result <- results[[paste0("censor_rate_", censor_rate)]]
    result_df <- data.frame(
      CensorRate = paste0(100 * censor_rate, "%"),
      Method = method_names,
      TC = result$TC,
      SD = result$SD,
      Cor = result$Cor,
      Time_s = result$Time
    )
    combined_results <- rbind(combined_results, result_df)
  }

  print(
    xtable::xtable(combined_results,
      caption = paste0(
        "Multi-index model ($d = ", d_true, "$): ",
        "trace correlation (TC, $\\uparrow$), ",
        "subspace distance (SD, $\\downarrow$), ",
        "mean correlation (Cor, $\\uparrow$)."
      ),
      label = "tab:multi_index"
    ),
    type = "latex",
    include.rownames = FALSE,
    caption.placement = "top",
    table.placement = "H"
  )
}

#' Table LaTeX : comparaison EDR index multiple, avec d_hat (BIC)
#'
#' @param results Sortie \code{$results} de
#'   \code{\link{run_simulation_multi_auto_d}}.
#' @param censor_rates Vecteur des taux de censure a inclure.
#' @param d_true Dimension EDR vraie (utilisee dans la legende).
#' @return Invisible ; imprime la table LaTeX sur la console.
#' @export
latex_table_multi_auto_d <- function(results, censor_rates, d_true) {
  if (!requireNamespace("xtable", quietly = TRUE)) {
    stop("Le package 'xtable' est necessaire pour latex_table_multi_auto_d().")
  }
  method_names <- c("SIR", "CSIR", "COX", "CSIR-WLS")
  combined_results <- data.frame()

  for (censor_rate in censor_rates) {
    result <- results[[paste0("censor_rate_", censor_rate)]]
    result_df <- data.frame(
      CensorRate = paste0(100 * censor_rate, "%"),
      Method = method_names,
      d_hat = result$d_hat,
      TC = result$TC,
      SD = result$SD,
      Cor = result$Cor,
      Time_s = result$Time
    )
    combined_results <- rbind(combined_results, result_df)
  }

  print(
    xtable::xtable(combined_results,
      caption = paste0(
        "Multi-index model ($d_{\\text{true}} = ", d_true, "$, ",
        "estimated $\\hat{d}$ via BIC): ",
        "trace correlation (TC, $\\uparrow$), ",
        "subspace distance (SD, $\\downarrow$), ",
        "mean correlation (Cor, $\\uparrow$). ",
        "COX is restricted to $d=1$ by construction."
      ),
      label = "tab:multi_index_auto_d"
    ),
    type = "latex",
    include.rownames = FALSE,
    caption.placement = "top",
    table.placement = "H"
  )
}

#' Table LaTeX : comparaison de screening (TPR / FPR / FDR)
#'
#' @param results Sortie \code{$results} de
#'   \code{\link{run_simulation_screening}}.
#' @param censor_rates Vecteur des taux de censure a inclure.
#' @return Invisible ; imprime la table LaTeX sur la console.
#' @export
latex_table_screening <- function(results, censor_rates) {
  if (!requireNamespace("xtable", quietly = TRUE)) {
    stop("Le package 'xtable' est necessaire pour latex_table_screening().")
  }
  method_names <- c("CWLS", "CSIS")
  combined_results <- data.frame()

  for (censor_rate in censor_rates) {
    result <- results[[paste0("censor_rate_", censor_rate)]]
    result_df <- data.frame(
      CensorRate = paste0(100 * censor_rate, "%"),
      Method = method_names,
      TPR = result$TPR,
      FPR = result$FPR,
      FDR = result$FDR,
      Time_s = result$Time
    )
    combined_results <- rbind(combined_results, result_df)
  }

  print(
    xtable::xtable(combined_results,
      caption = "Results for different censoring rates",
      label = "tab:combined_censor_rates"
    ),
    type = "latex",
    include.rownames = FALSE,
    caption.placement = "top",
    table.placement = "H"
  )
}

#' Table LaTeX : temps de calcul moyen (screening)
#'
#' @param results Sortie \code{$results} de
#'   \code{\link{run_simulation_screening}}.
#' @param censor_rates Vecteur des taux de censure a inclure.
#' @param n,p Taille d'echantillon et dimension (pour la legende).
#' @return Invisible ; imprime la table LaTeX sur la console.
#' @export
latex_time_table_screening <- function(results, censor_rates, n, p) {
  if (!requireNamespace("xtable", quietly = TRUE)) {
    stop("Le package 'xtable' est necessaire pour latex_time_table_screening().")
  }
  method_names <- c("CWLS", "CSIS")
  combined_results <- data.frame()

  for (censor_rate in censor_rates) {
    result <- results[[paste0("censor_rate_", censor_rate)]]
    result_df <- data.frame(
      CensorRate = paste0(100 * censor_rate, "%"),
      Method = method_names,
      Time_s = result$Time
    )
    combined_results <- rbind(combined_results, result_df)
  }

  print(
    xtable::xtable(combined_results,
      caption = paste0(
        "Average computational time in seconds per method ",
        "($n=", n, "$, $p=", p, "$)."
      ),
      label = "tab:time"
    ),
    type = "latex",
    include.rownames = FALSE,
    caption.placement = "top",
    table.placement = "H"
  )
}

#' Table LaTeX : robustesse du screening sous distributions non-elliptiques
#'
#' Table recapitulative (une ligne par distribution x taux de censure x
#' methode) pour repondre a un commentaire de reviewer sur la robustesse
#' de la Condition C.1 (voir \code{\link{generate_X}}, types
#' \code{"skewnormal"} / \code{"skewt"}).
#'
#' @param results_robustness Liste nommee par type de distribution de
#'   sorties \code{$results} de \code{\link{run_simulation_screening}}.
#' @param censor_rates Vecteur des taux de censure a inclure.
#' @return Invisible ; imprime la table LaTeX sur la console.
#' @export
latex_table_robustness_screening <- function(results_robustness, censor_rates) {
  if (!requireNamespace("xtable", quietly = TRUE)) {
    stop("Le package 'xtable' est necessaire pour latex_table_robustness_screening().")
  }
  method_names <- c("CWLS", "CSIS")
  combined_results <- data.frame()

  for (xt in names(results_robustness)) {
    res <- results_robustness[[xt]]
    for (censor_rate in censor_rates) {
      result <- res[[paste0("censor_rate_", censor_rate)]]
      result_df <- data.frame(
        Distribution = xt,
        CensorRate = paste0(100 * censor_rate, "%"),
        Method = method_names,
        TPR = result$TPR,
        FPR = result$FPR,
        FDR = result$FDR
      )
      combined_results <- rbind(combined_results, result_df)
    }
  }

  print(
    xtable::xtable(combined_results,
      caption = "Robustness of CWLS under non-elliptical distributions of x (response to Reviewer comment on Condition C.1)",
      label = "tab:robustness_C1"
    ),
    type = "latex",
    include.rownames = FALSE,
    caption.placement = "top",
    table.placement = "H"
  )
}

#' Table LaTeX : robustesse de l'estimation EDR (CSIR-WLS) sous
#' distributions non-elliptiques
#'
#' Meme principe que \code{\link{latex_table_robustness_screening}}, mais
#' pour la comparaison d'estimation EDR (\code{\link{run_simulation_edr}}).
#'
#' @param results_robustness Liste nommee par type de distribution de
#'   sorties \code{$results} de \code{\link{run_simulation_edr}}.
#' @param censor_rates Vecteur des taux de censure a inclure.
#' @return Invisible ; imprime la table LaTeX sur la console.
#' @export
latex_table_robustness_EDR <- function(results_robustness, censor_rates) {
  if (!requireNamespace("xtable", quietly = TRUE)) {
    stop("Le package 'xtable' est necessaire pour latex_table_robustness_EDR().")
  }
  method_names <- c("SIR", "CSIR", "COX", "CSIR-WLS")
  combined_results <- data.frame()

  for (xt in names(results_robustness)) {
    res <- results_robustness[[xt]]
    for (censor_rate in censor_rates) {
      result <- res[[paste0("censor_rate_", censor_rate)]]
      result_df <- data.frame(
        Distribution = xt,
        CensorRate = paste0(100 * censor_rate, "%"),
        Method = method_names,
        PC = result$PC,
        Cor = result$COR,
        Frob = result$Frob,
        Time_s = result$Time
      )
      combined_results <- rbind(combined_results, result_df)
    }
  }

  print(
    xtable::xtable(combined_results,
      caption = "Robustness of CSIR-WLS EDR estimation under non-elliptical distributions of x (response to Reviewer comment on Condition C.1)",
      label = "tab:robustness_C1_EDR"
    ),
    type = "latex",
    include.rownames = FALSE,
    caption.placement = "top",
    table.placement = "H"
  )
}
