#' Criteres d'evaluation pour comparer des methodes
#'
#' @description
#' Fonctions de comparaison utilisees dans les etudes de simulation :
#' pour la selection de variables (TPR/FPR/FDR), et pour l'estimation de
#' directions EDR (cosinus / distance de projection pour d = 1, via
#' \code{\link{costeta}} / \code{\link{Dproj}} du fichier
#' \code{01-linear-algebra-utils.R} ; correlation de trace / distance de
#' sous-espace pour d > 1, ci-dessous).
#'
#' @name evaluation-metrics
NULL

#' TPR / FPR / FDR de la selection de variables (version detaillee)
#'
#' Compare l'ensemble des variables selectionnees par une methode
#' (\code{estimated_active}) a l'ensemble des vraies variables actives
#' (\code{true_active}), et calcule les indicateurs de performance
#' classiques de detection.
#'
#' @param estimated_active Indices des variables selectionnees.
#' @param true_active Indices des vraies variables actives.
#' @param total_features Nombre total de variables p (pour calculer TN).
#' @return Une liste : \code{TPR}, \code{FPR}, \code{FDR}, \code{TP},
#'   \code{FP}, \code{FN}, \code{TN}.
#' @export
compute_TPR_FPR <- function(estimated_active, true_active, total_features) {
  estimated_active <- unique(estimated_active)
  true_active <- unique(true_active)

  TP <- length(intersect(estimated_active, true_active))
  FP <- length(setdiff(estimated_active, true_active))
  FN <- length(setdiff(true_active, estimated_active))

  TN <- total_features - (TP + FP + FN)

  if (TP + FP + FN > total_features) {
    stop("The sum of TP, FP, and FN exceeds the total number of features (p). Check your data!")
  }

  TPR <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  FPR <- ifelse((FP + TN) > 0, FP / (FP + TN), 0)
  FDR <- ifelse((FP + TP) > 0, FP / (FP + TP), 0)

  return(list(TPR = TPR, FPR = FPR, FDR = FDR, TP = TP, FP = FP, FN = FN, TN = TN))
}

#' TPR / FPR / FDR de la selection de variables (version compacte)
#'
#' Version simplifiee de \code{\link{compute_TPR_FPR}}, utilisee dans les
#' boucles de simulation de screening (CWLS vs CSIS) : ne renvoie que
#' \code{TPR}, \code{FPR}, \code{FDR}.
#'
#' @inheritParams compute_TPR_FPR
#' @param p Nombre total de variables (equivalent a \code{total_features}).
#' @return Une liste : \code{TPR}, \code{FPR}, \code{FDR}.
#' @export
compute_metrics <- function(estimated_active, true_active, p) {
  estimated_active <- unique(estimated_active)

  TP <- sum(estimated_active %in% true_active)
  FP <- sum(!(estimated_active %in% true_active))
  FN <- length(true_active) - TP
  TN <- p - length(true_active) - FP

  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)
  FDR <- ifelse((FP + TP) == 0, 0, FP / (FP + TP))

  list(TPR = TPR, FPR = FPR, FDR = FDR)
}

#' Correlation de trace entre deux sous-espaces (generalisation de costeta)
#'
#' Generalisation de \code{\link{costeta}} au cas multi-directions (d > 1) :
#' orthonormalise \code{B1} et \code{B2} par QR, puis calcule
#' \code{sum((t(Q1) \%*\% Q2)^2) / d}. Vaut 1 si les deux sous-espaces
#' coincident, 0 s'ils sont orthogonaux.
#'
#' @param B1,B2 Matrices (ou vecteurs) p x d de directions EDR.
#' @return Un scalaire dans [0, 1].
#' @export
trace_corr <- function(B1, B2) {
  B1 <- as.matrix(B1)
  B2 <- as.matrix(B2)
  Q1 <- qr.Q(qr(B1))
  Q2 <- qr.Q(qr(B2))
  d <- ncol(Q1)
  M <- t(Q1) %*% Q2
  return(sum(M^2) / d)
}

#' Distance entre deux sous-espaces (generalisation de Dproj)
#'
#' Generalisation de \code{\link{Dproj}} au cas multi-directions (d > 1) :
#' orthonormalise \code{B1} et \code{B2} par QR, construit les matrices
#' de projection \code{P1}, \code{P2}, et renvoie
#' \code{norm(P1 - P2, "F") / sqrt(2 * d)}.
#'
#' @param B1,B2 Matrices (ou vecteurs) p x d de directions EDR.
#' @return Un scalaire >= 0 (0 = sous-espaces identiques).
#' @export
subspace_dist <- function(B1, B2) {
  B1 <- as.matrix(B1)
  B2 <- as.matrix(B2)
  Q1 <- qr.Q(qr(B1))
  Q2 <- qr.Q(qr(B2))
  P1 <- Q1 %*% t(Q1)
  P2 <- Q2 %*% t(Q2)
  d <- ncol(Q1)
  return(norm(P1 - P2, type = "F") / sqrt(2 * d))
}

#' Correlation, avec retour a 0 si NA (estimation degeneree)
#'
#' @param a,b Vecteurs numeriques.
#' @return \code{abs(cor(a, b))}, ou 0 si NA.
#' @export
safe_cor <- function(a, b) {
  out <- suppressWarnings(cor(a, b))
  ifelse(is.na(out), 0, abs(out))
}

#' Correlation moyenne entre plusieurs directions EDR vraies/estimees
#'
#' Pour chaque colonne k de \code{Beta_true} / \code{beta_hat}, calcule
#' \code{|cor(X \%*\% Beta_true[,k], X \%*\% beta_hat[,k])|} (0 si NA), et
#' renvoie la moyenne sur les d directions.
#'
#' @param X Matrice n x p de predicteurs.
#' @param Beta_true,beta_hat Matrices (ou vecteurs) p x d de directions
#'   EDR vraies / estimees.
#' @return Un scalaire dans [0, 1].
#' @export
safe_cor_multi <- function(X, Beta_true, beta_hat) {
  Beta_true <- as.matrix(Beta_true)
  beta_hat <- as.matrix(beta_hat)
  d <- ncol(Beta_true)
  cors <- numeric(d)
  for (k in 1:d) {
    c_val <- suppressWarnings(cor(
      X %*% Beta_true[, k],
      X %*% beta_hat[, k]
    ))
    cors[k] <- ifelse(is.na(c_val), 0, abs(c_val))
  }
  return(mean(cors))
}
