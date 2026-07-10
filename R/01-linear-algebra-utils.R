#' Fonctions d'algebre lineaire de base
#'
#' @description
#' Petites fonctions utilitaires utilisees dans tout le package : mesurer
#' la similarite entre deux directions (cosinus), normaliser un vecteur,
#' construire une matrice de projection, calculer une covariance biaisee,
#' et resoudre un probleme generalise aux valeurs propres. Elles viennent
#' toutes du fichier \code{myfunctions.R} d'origine.
#'
#' @name linear-algebra-utils
NULL

#' Cosinus de l'angle entre deux vecteurs (ou matrices)
#'
#' Calcule le produit scalaire de \code{x} et \code{y} divise par le
#' produit de leurs normes de Frobenius. Sert a mesurer la similarite
#' directionnelle entre deux estimations de direction EDR (ex : beta
#' estime vs beta vrai) dans le cas d = 1 (une seule direction). Pour
#' d > 1, utiliser \code{\link{trace_corr}} ou \code{\link{subspace_dist}}.
#'
#' @param x,y Deux vecteurs numeriques (ou matrices) de meme longueur.
#' @return Un scalaire = cos(theta) entre \code{x} et \code{y}.
#' @export
costeta <- function(x, y) {
  t <- crossprod(x, y) / (norm(as.matrix(x), "f") * norm(as.matrix(y), "f"))
  return(t)
}

#' Normaliser un vecteur (norme euclidienne 1)
#'
#' @param v Vecteur numerique.
#' @return Le vecteur \code{v} divise par sa norme euclidienne. Si
#'   \code{v} est nul, il est renvoye tel quel (pour eviter une division
#'   par 0).
#' @export
normalize_vector <- function(v) {
  magnitude <- sqrt(sum(v^2))
  if (magnitude == 0) {
    return(v)
  } else {
    return(v / magnitude)
  }
}

#' Normaliser chaque colonne d'une matrice
#'
#' Applique \code{\link{normalize_vector}} a chaque colonne de \code{mat}.
#'
#' @param mat Matrice numerique.
#' @return Matrice de meme dimension, colonnes normalisees.
#' @export
normalize_matrix <- function(mat) {
  apply(mat, MARGIN = 2, FUN = function(x) x / sqrt(sum(x^2)))
}

#' Matrice de projection orthogonale sur l'espace engendre par u
#'
#' P = u \%*\% t(u) / ||u||^2.
#'
#' @param u Vecteur ou matrice de directions.
#' @return Matrice de projection P (meme nombre de lignes que \code{u}).
#' @export
projec <- function(u) {
  P <- u %*% t(u) / norm(u, "2")^2
  return(P)
}

#' Distance de Frobenius entre deux espaces de projection (cas d = 1)
#'
#' Mesure la distance entre les espaces engendres par \code{u} et
#' \code{v}, via la norme de Frobenius de la difference des deux
#' matrices de projection (\code{\link{projec}}). Utile pour comparer
#' une direction EDR estimee a la direction vraie, dans le cas d = 1.
#' Pour d > 1, utiliser \code{\link{subspace_dist}}.
#'
#' @param u,v Deux vecteurs (ou matrices) de directions.
#' @return Un scalaire = distance de Frobenius entre les deux projections.
#' @export
Dproj <- function(u, v) {
  P1 <- projec(u)
  P2 <- projec(v)
  D <- norm(as.matrix(P1 - P2), "f")
  return(D)
}

#' Covariance empirique biaisee (denominateur n)
#'
#' Identique a \code{stats::cov(x)} mais avec un denominateur \code{n}
#' (biaise / population) au lieu de \code{n - 1}.
#'
#' @param x Matrice de donnees.
#' @param use Methode de gestion des valeurs manquantes, transmise a
#'   \code{stats::cov()} (defaut \code{"pair"} = pairwise.complete.obs).
#' @return Matrice de covariance p x p.
#' @export
cov1 <- function(x, use = "pair") {
  n <- nrow(x)
  c <- cov(x, use = use) * (n - 1) / n
  c
}

#' Probleme generalise aux valeurs propres via decomposition de Cholesky
#'
#' Resout m1 v = lambda * m2 v en utilisant une decomposition de Cholesky
#' de m2 (m2 = t(cm2) \%*\% cm2) suivie d'une SVD de la matrice
#' transformee. Utilise par \code{\link{cen.sir}}, \code{\link{double.slice}}
#' et \code{\link{double.slice2}} pour obtenir les directions EDR.
#'
#' @param m1,m2 Deux matrices carrees (\code{m2} doit etre definie positive).
#' @param tol Seuil numerique pour la stabilite de l'inversion (defaut 1e-5).
#' @return Une liste avec \code{values} (valeurs propres) et
#'   \code{vectors} (vecteurs propres).
#' @export
eigen.decomp <- function(m1, m2, tol = 1e-5) {
  # cholesky decomposition of m2 = t(cm2) %*% cm2
  # and inverse(t(cm2)) = t(inverse(cm2))
  cm2 <- chol(m2)
  inv.cm2 <- cm2
  idx.cm2 <- diag(cm2) > tol
  inv.cm2[idx.cm2, idx.cm2] <- solve(cm2[idx.cm2, idx.cm2])
  m3 <- t(inv.cm2) %*% m1 %*% (inv.cm2)
  svd.m3 <- svd(m3)
  eigen.value <- svd.m3$d
  eigen.vect <- inv.cm2 %*% svd.m3$v
  result <- list(values = eigen.value, vectors = eigen.vect)
  result
}
