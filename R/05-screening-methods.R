#' Methodes de screening de variables (concurrentes de CSIR-WLS)
#'
#' @description
#' Methodes de screening de variables pour la tres grande dimension,
#' utilisees dans le papier comme methodes concurrentes de CSIR-WLS :
#' CSIS (concordance index, Cheng, Li & Wang 2023), DCSIS (distance
#' correlation, Li, Zhong & Zhu 2012), SIRS (Zhu, Li, Li & Zhu 2011) et
#' SIS (Fan & Lv 2008). Deux fonctions "wrapper" (\code{apply_cwls},
#' \code{apply_csis}) simplifient leur appel dans les boucles de
#' simulation.
#'
#' @name screening-methods
NULL

#' CSIS : screening base sur l'indice de concordance
#'
#' Un screening de variables model-free et data-adaptive pour donnees en
#' tres grande dimension, y compris les donnees de survie. La methode
#' est basee sur l'indice de concordance qui mesure la concordance entre
#' vecteurs aleatoires, meme si l'un d'eux est un objet de survie
#' \code{Surv}. Cette methode basee sur une correlation de rang ne
#' necessite pas de specifier un modele de regression, et s'applique de
#' maniere robuste en presence de censure et de queues lourdes.
#'
#' Pour chaque variable, on calcule sa concordance avec \code{Y} ; si
#' elle est < 0.5 on prend (1 - concordance) pour se ramener a une
#' mesure d'association toujours dans [0.5, 1]. Les variables sont
#' triees par cette mesure et les \code{nsis} meilleures sont gardees.
#'
#' @param X Matrice n x p de predicteurs. Chaque ligne est un vecteur
#'   d'observation.
#' @param Y Vecteur reponse de dimension n. Pour les modeles de survie,
#'   \code{Y} doit etre un objet de classe \code{Surv} (package
#'   \pkg{survival}).
#' @param nsis Nombre de predicteurs a garder (defaut \code{n/log(n)}).
#' @return Les indices des \code{nsis} variables les plus concordantes.
#' @references Cheng X, Li G, Wang H. The concordance filter: an
#'   adaptive model-free feature screening procedure. Computational
#'   Statistics, 2023.
#' @export
CSIS <- function(X, Y, nsis = (dim(X)[1]) / log(dim(X)[1])) {
  if (dim(X)[1] != length(Y)) {
    stop("X and Y should have same number of rows!")
  }
  if (missing(X) | missing(Y)) {
    stop("The data is missing!")
  }
  if (TRUE %in% (is.na(X) | is.na(Y) | is.na(nsis))) {
    stop("The input vector or matrix cannot have NA!")
  }
  n <- dim(X)[1]
  p <- dim(X)[2]
  B <- vector(mode = "numeric", length = p)
  Cindex <- vector(mode = "numeric", length = p)
  if (n * p <= 2000000) {
    for (i in 1:p) {
      Cindex[i] <- survival::concordancefit(Y, X[, i])$concordance
    }
  } else {
    if (!requireNamespace("doParallel", quietly = TRUE) ||
      !requireNamespace("foreach", quietly = TRUE)) {
      stop("Les packages 'doParallel' et 'foreach' sont necessaires pour CSIS() en tres grande dimension.")
    }
    cores <- parallel::detectCores(logical = FALSE)
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl, cores = cores)
    j <- NULL
    Cindex <- foreach::foreach(
      j = 1:p, .combine = "c",
      .packages = c("survival")
    ) %dopar%
      survival::concordancefit(Y, X[, j])$concordance
    foreach::registerDoSEQ()
    parallel::stopCluster(cl)
  }
  num <- which(Cindex < 0.5)
  B[num] <- 1 - Cindex[num]
  B[-num] <- Cindex[-num]
  A <- order(B, decreasing = TRUE)
  return(A[1:nsis])
}

#' DCSIS : screening base sur la correlation de distance
#'
#' Screening de variables base sur la correlation de distance (distance
#' correlation) entre chaque colonne de \code{X} et \code{Y}. Ne
#' fonctionne pas si \code{Y} est un objet \code{Surv}.
#'
#' @param X Matrice n x p de predicteurs.
#' @param Y Reponse (vecteur, PAS un objet \code{Surv}).
#' @param nsis Nombre de variables a garder (defaut \code{n/log(n)}).
#' @return Les indices des \code{nsis} variables les plus correlees
#'   (distance correlation).
#' @references Li, R., Zhong, W. and Zhu, L. (2012). Feature screening
#'   via distance correlation learning. JASA 107(499), 1129-1139.
#' @export
DCSIS <- function(X, Y, nsis = (dim(X)[1]) / log(dim(X)[1])) {
  if (dim(X)[1] != length(Y)) {
    stop("X and Y should have same number of rows!")
  }
  if (missing(X) | missing(Y)) {
    stop("The data is missing!")
  }
  if (TRUE %in% (is.na(X) | is.na(Y) | is.na(nsis))) {
    stop("The input vector or matrix cannot have NA!")
  }
  if (inherits(Y, "Surv")) {
    stop("DCSIS can not implemented with object  of Surv")
  }
  n <- dim(X)[1]
  p <- dim(X)[2]
  B <- matrix(1, n, 1)
  C <- matrix(1, 1, p)
  sxy1 <- matrix(0, n, p)
  sxy2 <- matrix(0, n, p)
  sxy3 <- matrix(0, n, 1)
  sxx1 <- matrix(0, n, p)
  syy1 <- matrix(0, n, 1)
  for (i in 1:n) {
    XX1 <- abs(X - B %*% X[i, ])
    YY1 <- sqrt(apply((Y - B %*% Y[i])^2, 1, sum))
    sxy1[i, ] <- apply(XX1 * (YY1 %*% C), 2, mean)
    sxy2[i, ] <- apply(XX1, 2, mean)
    sxy3[i, ] <- mean(YY1)
    XX2 <- XX1^2
    sxx1[i, ] <- apply(XX2, 2, mean)
    YY2 <- YY1^2
    syy1[i, ] <- mean(YY2)
  }
  SXY1 <- apply(sxy1, 2, mean)
  SXY2 <- apply(sxy2, 2, mean) * apply(sxy3, 2, mean)
  SXY3 <- apply(sxy2 * (sxy3 %*% C), 2, mean)
  SXX1 <- apply(sxx1, 2, mean)
  SXX2 <- apply(sxy2, 2, mean)^2
  SXX3 <- apply(sxy2^2, 2, mean)
  SYY1 <- apply(syy1, 2, mean)
  SYY2 <- apply(sxy3, 2, mean)^2
  SYY3 <- apply(sxy3^2, 2, mean)
  dcovXY <- sqrt(SXY1 + SXY2 - 2 * SXY3)
  dvarXX <- sqrt(SXX1 + SXX2 - 2 * SXX3)
  dvarYY <- sqrt(SYY1 + SYY2 - 2 * SYY3)
  dcorrXY <- dcovXY / sqrt(dvarXX * dvarYY)
  A <- order(dcorrXY, decreasing = TRUE)
  return(A[1:nsis])
}

#' SIRS : screening base sur SIR (Sliced Inverse Regression Screening)
#'
#' Standardise \code{X}, trie les observations selon \code{Y}, puis
#' calcule pour chaque variable une statistique USIRS basee sur les
#' moyennes cumulees de x triees par y. Garde les \code{nsis} variables
#' avec la plus grande statistique. Ne fonctionne pas si \code{Y} est un
#' objet \code{Surv}.
#'
#' @param X Matrice n x p de predicteurs.
#' @param Y Reponse (vecteur, PAS un objet \code{Surv}).
#' @param nsis Nombre de variables a garder (defaut \code{n/log(n)}).
#' @return Les indices des \code{nsis} variables retenues.
#' @references Zhu, L.-P., Li, L., Li, R. and Zhu, L.-X. (2011).
#'   Model-free feature screening for ultrahigh-dimensional data. JASA
#'   106(496), 1464-1475.
#' @export
SIRS <- function(X, Y, nsis = (dim(X)[1]) / log(dim(X)[1])) {
  if (dim(X)[1] != length(Y)) {
    stop("X and Y should have same number of rows!")
  }
  if (missing(X) | missing(Y)) {
    stop("The data is missing!")
  }
  if (TRUE %in% (is.na(X) | is.na(Y) | is.na(nsis))) {
    stop("The input vector or matrix cannot have NA!")
  }
  if (inherits(Y, "Surv")) {
    stop("SIRS can not implemented with object  of Surv")
  }
  n <- dim(X)[1]
  posit <- order(Y, decreasing = FALSE)
  Y <- Y[posit]
  X <- X[posit, ]
  xx <- (X - as.matrix(rep(1, n)) %*% apply(X, 2, mean)) / (as.matrix(rep(1, n)) %*% apply(X, 2, sd))
  B <- matrix(1, n, n)
  B[!upper.tri(B, diag = TRUE)] <- 0
  USIRS <- apply((t(xx) %*% B / n)^2, 1, mean)
  A <- order(USIRS, decreasing = TRUE)
  B <- A[1:nsis]
  return(B)
}

#' SIS : Sure Independence Screening (Fan & Lv, 2008)
#'
#' Version la plus simple du screening : score de chaque variable =
#' \code{|t(X) \%*\% Y|} (produit scalaire entre la colonne de X et Y en
#' valeur absolue) ; garde les \code{nsis} variables avec le score le
#' plus eleve. Ne fonctionne pas si \code{Y} est un objet \code{Surv}.
#'
#' @param X Matrice n x p de predicteurs.
#' @param Y Reponse (vecteur, PAS un objet \code{Surv}).
#' @param nsis Nombre de variables a garder (defaut \code{n/log(n)}).
#' @return Les indices des \code{nsis} variables retenues.
#' @references Fan, J. and Lv, J. (2008). Sure independence screening
#'   for ultrahigh dimensional feature space. JRSS-B 70(5), 849-911.
#' @export
SIS <- function(X, Y, nsis = (dim(X)[1]) / log(dim(X)[1])) {
  if (dim(X)[1] != length(Y)) {
    stop("X and Y should have same number of rows!")
  }
  if (missing(X) | missing(Y)) {
    stop("The data is missing!")
  }
  if (TRUE %in% (is.na(X) | is.na(Y) | is.na(nsis))) {
    stop("The input vector or matrix cannot have NA!")
  }
  if (inherits(Y, "Surv")) {
    stop("SIS can not implemented with object  of Surv")
  }
  A <- order(abs(t(X) %*% Y), decreasing = TRUE)
  return(A[1:nsis])
}

#' Wrapper : selection de variables par CWLS (censored WLS)
#'
#' Applique \code{\link{cen.wls}} et renvoie uniquement les indices des
#' variables selectionnees. Utilise dans les boucles de simulation de
#' comparaison CWLS vs CSIS.
#'
#' @param X Matrice n x p de predicteurs.
#' @param time Temps de survie.
#' @param event Indicateur d'evenement (1 = evenement, 0 = censure).
#' @param c Plancher de la survie conditionnelle estimee (defaut 0.05).
#' @param n.slice Nombre de tranches (defaut 10).
#' @param cn1,cn2 Constantes de penalite BIC.
#' @return Vecteur des indices des variables selectionnees.
#' @export
apply_cwls <- function(X, time, event, c = 0.05, n.slice = 10,
                        cn1 = 0.1, cn2 = 1) {
  cwls_result <- cen.wls(X, time, event,
    c = c, n.slice = n.slice,
    cn1 = cn1, cn2 = cn2, choose.dir = FALSE
  )
  cwls_result$select
}

#' Wrapper : selection de variables par CSIS
#'
#' Applique \code{\link{CSIS}} avec une valeur par defaut de
#' \code{nsis = min(p, n/log(n))} si non fournie.
#'
#' @param X Matrice n x p de predicteurs.
#' @param time Temps de survie.
#' @param event Indicateur d'evenement (1 = evenement, 0 = censure).
#' @param nsis Nombre de predicteurs a garder (calcule automatiquement
#'   si \code{NULL}).
#' @return Vecteur des indices des variables selectionnees.
#' @export
apply_csis <- function(X, time, event, nsis = NULL) {
  if (is.null(nsis)) {
    nsis <- floor(min(ncol(X), nrow(X) / log(nrow(X))))
  }
  CSIS(X, survival::Surv(time, event), nsis = nsis)
}
