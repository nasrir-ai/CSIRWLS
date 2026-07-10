#' Fonctions graphiques
#'
#' @description
#' Fonctions de visualisation utilisees dans les analyses du papier :
#' methodes \code{plot} pour les objets SIR/CSIR/double-slicing, et
#' fonctions plus specifiques (projection vs temps, Kaplan-Meier par
#' groupes de risque, temps de calcul).
#'
#' @name plotting
NULL

#' Nuage de points : chaque direction EDR vs temps de survie
#'
#' Trace, pour un objet SIR/csir/ds, un nuage de points de chaque
#' direction EDR selectionnee (\code{object$x \%*\% evec}) contre le
#' temps de survie (ou son log), en distinguant par couleur/symbole les
#' observations censurees et non censurees.
#'
#' @param object Objet contenant \code{x}, \code{evec}, \code{y}, \code{delta}.
#' @param which Indices des directions a tracer (maximum 6).
#' @param logY \code{TRUE} pour tracer \code{log(y)} au lieu de \code{y}.
#' @return Rien ; affiche des graphiques via \code{plot()}.
#' @export
plot.sir <- function(object, which = 1:ncol(object$evec), logY = FALSE) {
  if (length(which) > 6) {
    stop("no more than 6 directions can be plot at one time")
  }

  dirs <- object$x %*% object$evec[, which]
  if (logY) {
    surv <- log(object$y)
    labY <- "log(survival time)"
  } else {
    surv <- object$y
    labY <- "survival time"
  }
  labX <- colnames(object$evec)[which]

  if (length(which) == 2) {
    graphics::par(mfrow = c(1, 2))
  }
  if (length(which) == 3 || length(which) == 4) {
    graphics::par(mfrow = c(2, 2))
  }
  if (length(which) == 5 || length(which) == 6) {
    graphics::par(mfrow = c(2, 3))
  }

  for (i in 1:length(which)) {
    plot(dirs[, i][object$delta == 1], surv[object$delta == 1],
      col = "blue", type = "p", pch = 16,
      xlab = labX[i], ylab = labY, main = "",
      xlim = range(dirs[, i]), ylim = range(surv)
    )
    points(dirs[, i][object$delta == 0], surv[object$delta == 0],
      col = "red", type = "p", pch = 10
    )
    legend(min(dirs[, i]) + 0.05 * (max(dirs[, i]) - min(dirs[, i])),
      min(surv) + 0.95 * (max(surv) - min(surv)),
      legend = c("censored", "uncensored"),
      col = c("red", "blue"), pch = c(10, 16)
    )
  }
}

#' @rdname plot.sir
#' @export
plot.csir <- function(object, which = 1:ncol(object$evec), logY = FALSE) {
  plot.sir(object, which = which, logY = logY)
}

#' @rdname plot.sir
#' @export
plot.ds <- function(object, which = 1:ncol(object$evec), logY = FALSE) {
  plot.sir(object, which = which, logY = logY)
}

#' Nuage de points 3D pour deux directions EDR vs temps de survie
#'
#' Trace un nuage de points 3D (via \code{scatterplot3d}) de deux
#' directions EDR contre le temps de survie, sous differents angles de
#' vue, en distinguant censure/non-censure.
#'
#' @param object Objet contenant \code{x}, \code{evec}, \code{y}, \code{delta}.
#' @param which Exactement 2 indices de directions a tracer.
#' @param angles Vecteur d'angles de vue (maximum 6).
#' @param z.plane Hauteur d'un plan horizontal optionnel a tracer.
#' @param logY \code{TRUE} pour utiliser \code{log(y)} comme axe Z.
#' @return Rien ; affiche des graphiques 3D.
#' @export
plot.3d.sir <- function(object, which = 1:2, angles = c(60, 120), z.plane = NULL, logY = FALSE) {
  if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
    stop("Le package 'scatterplot3d' est necessaire pour plot.3d.sir().")
  }

  if (length(which) != 2) {
    stop("Please indicate 2 and only 2 directions for plot.3d\n")
  }
  if (length(angles) > 6) {
    stop("no more than 6 directions can be plot at one time")
  }

  dirs <- object$x %*% object$evec[, which]
  if (logY) {
    surv <- log(object$y)
    labZ <- "log(survival time)"
  } else {
    surv <- object$y
    labZ <- "survival time"
  }
  labX <- colnames(object$evec[, which])

  if (length(angles) == 2) {
    graphics::par(mfrow = c(1, 2))
  }
  if (length(angles) == 3 || length(which) == 4) {
    graphics::par(mfrow = c(2, 2))
  }
  if (length(angles) == 5 || length(which) == 6) {
    graphics::par(mfrow = c(2, 3))
  }

  d <- object$delta

  for (i in 1:length(angles)) {
    s3d <- scatterplot3d::scatterplot3d(dirs[, 1], dirs[, 2],
      xlab = labX[1], grid = TRUE,
      ylab = labX[2], zlab = labZ, surv, angle = angles[i], type = "n"
    )
    s3d$points3d(dirs[d == 1, 1], dirs[d == 1, 2], surv[d == 1],
      col = "blue", type = "p", pch = 16
    )
    s3d$points3d(dirs[d == 0, 1], dirs[d == 0, 2], surv[d == 0],
      col = "red", type = "p", pch = 10
    )
    if (!is.null(z.plane)) {
      s3d$plane3d(z.plane, 0, 0)
    }
  }
}

#' Nuage de points : projection EDR vs temps observe (censure/echec)
#'
#' Trace la projection \code{proj} contre le temps observe, avec un
#' symbole different pour les evenements observes ("." pch=16) et les
#' observations censurees ("+" pch=3). Utilise dans les analyses
#' DLBCL/PBC pour comparer Cox, CSIR et CSIR-WLS.
#'
#' @param proj Vecteur de projections (ex : \code{X \%*\% beta.hat}).
#' @param time Vecteur des temps observes.
#' @param delta Indicateur d'evenement (1 = evenement, 0 = censure).
#' @param title Titre du graphique.
#' @return Rien ; affiche un graphique.
#' @export
plot_projection <- function(proj, time, delta, title) {
  plot(proj[delta == 1], time[delta == 1],
    col = "blue", pch = 16,
    xlab = "Projected direction", ylab = "Time (days)",
    main = title,
    xlim = range(proj), ylim = range(time)
  )
  points(proj[delta == 0], time[delta == 0], col = "red", pch = 3)
  legend("bottomright",
    legend = c("Failure", "Censored"),
    col = c("blue", "red"), pch = c(16, 3), bty = "n"
  )
}

#' Decouper une projection en 3 groupes de risque (tertiles)
#'
#' Decoupe un vecteur de projection en trois groupes ("Low risk",
#' "Medium risk", "High risk") selon les tertiles empiriques. Utilise
#' pour construire les courbes de Kaplan-Meier par groupe de risque
#' (Figure 3 de l'analyse DLBCL).
#'
#' @param proj Vecteur de projections (ex : \code{X \%*\% beta.hat}).
#' @return Un facteur a 3 niveaux : "Low risk", "Medium risk", "High risk".
#' @export
make_risk_groups <- function(proj) {
  q <- quantile(proj, probs = c(1 / 3, 2 / 3))
  cut(proj,
    breaks = c(-Inf, q[1], q[2], Inf),
    labels = c("Low risk", "Medium risk", "High risk"),
    right = FALSE
  )
}

#' Courbe de Kaplan-Meier par groupe de risque (ggplot2)
#'
#' Construit un graphique ggplot2 des courbes de Kaplan-Meier pour une
#' methode donnee, a partir d'un data.frame long contenant \code{time},
#' \code{surv}, \code{strata} et \code{method}.
#'
#' @param data Data.frame long avec colonnes \code{time}, \code{surv},
#'   \code{strata}, \code{method}.
#' @param method_name Nom de la methode a filtrer (ex : "Cox", "CSIR",
#'   "CSIR-WLS").
#' @return Un objet ggplot.
#' @export
create_km_plot <- function(data, method_name) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Le package 'ggplot2' est necessaire pour create_km_plot().")
  }
  ggplot2::ggplot(
    data[data$method == method_name, ],
    ggplot2::aes(x = .data$time, y = .data$surv, color = .data$strata)
  ) +
    ggplot2::geom_step(linewidth = 1) +
    ggplot2::scale_color_manual(
      values = c("red", "green", "blue"),
      labels = c("High", "Medium", "Low"),
      name = "Risk Group"
    ) +
    ggplot2::labs(
      title = paste("Kaplan-Meier Curves by", method_name, "Projection"),
      x = "Time (days)",
      y = "Survival Probability"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}

#' Boxplot des temps de calcul par methode et taux de censure
#'
#' Construit un boxplot ggplot2 des temps de calcul (par repetition),
#' groupes par methode et taux de censure. Attend en entree la sortie de
#' \code{\link{run_simulation_screening}} (liste avec \code{$Time_by_rate}).
#'
#' @param sim_output Sortie de \code{\link{run_simulation_screening}}.
#' @param n,p Taille d'echantillon et dimension (pour le titre).
#' @return Un objet ggplot.
#' @export
plot_time_by_censor <- function(sim_output, n, p) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("reshape2", quietly = TRUE)) {
    stop("Les packages 'ggplot2' et 'reshape2' sont necessaires pour plot_time_by_censor().")
  }
  method_names <- c("CWLS", "CSIS")

  long_list <- list()
  for (cr_name in names(sim_output$Time_by_rate)) {
    cr_value <- as.numeric(gsub("censor_rate_", "", cr_name))
    time_df <- as.data.frame(t(sim_output$Time_by_rate[[cr_name]]))
    colnames(time_df) <- method_names
    time_df$Rep <- 1:nrow(time_df)
    time_df$CensorRate <- paste0(cr_value * 100, "%")
    long_list[[cr_name]] <- time_df
  }

  time_df_all <- do.call(rbind, long_list)

  time_long <- reshape2::melt(time_df_all,
    id.vars = c("Rep", "CensorRate"),
    variable.name = "Method",
    value.name = "Time_s"
  )

  ggplot2::ggplot(time_long, ggplot2::aes(x = .data$Method, y = .data$Time_s, fill = .data$CensorRate)) +
    ggplot2::geom_boxplot(alpha = 0.7, position = ggplot2::position_dodge(0.8)) +
    ggplot2::labs(
      title = paste0("Computational Time (n=", n, ", p=", p, ")"),
      x = "Method", y = "Time (seconds)"
    ) +
    ggplot2::theme_minimal()
}
