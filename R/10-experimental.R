#' Fonctions experimentales / incompletes
#'
#' @description
#' Ce fichier regroupe le code presents dans \code{myfunctions.R} qui
#' n'est pas termine ou qui depend de variables non definies dans la
#' fonction elle-meme. Il est garde tel quel (comme demande) mais
#' clairement isole ici pour eviter qu'il ne soit confondu avec le
#' pipeline principal, robuste, de \code{04-variable-selection-wls.R}.
#'
#' @name experimental
NULL

#' WLS-SIR avec controle du taux de fausses decouvertes (INCOMPLET)
#'
#' Variante de \code{\link{wls.sir}} censee ajouter un controle du taux
#' de fausses decouvertes (FDR) apres l'etape WLS+SIR habituelle.
#'
#' @section Attention :
#' La partie FDR en fin de fonction utilise des objets (\code{t},
#' \code{Ta}, \code{alf}, \code{W}, \code{q}) qui ne sont jamais
#' definis dans le corps de la fonction -- ils doivent exister deja
#' dans l'environnement appelant (typiquement un script de simulation
#' qui definit un seuil \code{alf}, une statistique de test \code{Ta},
#' etc.). Telle quelle, cette fonction \strong{plantera} si ces objets
#' n'existent pas deja en memoire au moment de l'appel. Le calcul de
#' FDP et AP est fait mais n'est pas renvoye dans la liste finale
#' (bug d'origine, non corrige ici pour rester fidele au code source).
#'
#' Si vous voulez utiliser reellement le controle FDR, il faudra
#' clarifier/re-ecrire cette fonction pour qu'elle prenne \code{t},
#' \code{Ta}, \code{alf}, \code{W}, \code{q} en parametres explicites.
#'
#' @inheritParams wls.sir
#' @return Une liste : \code{wls}, \code{select}, \code{betahat} (comme
#'   \code{\link{wls.sir}}).
#' @keywords internal
#' @export
wls.sir.fdr <- function(x, y, nslice = 10, cn1 = 0.1, cn2 = 1, choose.dir = FALSE,
                         categorical = FALSE, ndim = 1) {
  if (!requireNamespace("dr", quietly = TRUE)) {
    stop("Le package 'dr' est necessaire pour wls.sir.fdr().")
  }
  if (!requireNamespace("Rdimtools", quietly = TRUE)) {
    stop("Le package 'Rdimtools' est necessaire pour wls.sir.fdr().")
  }

  n <- dim(x)[1]
  p <- dim(x)[2]
  h <- nslice
  x <- scale(x, scale = FALSE)

  if (categorical == TRUE) {
    index <- as.numeric(factor(y))
    nh <- summary(factor(y))
    h <- length(nh)
  }
  if (categorical == FALSE) {
    slice <- dr::dr.slices.arc(y, h)
    index <- slice$slice.indicator
    nh <- slice$slice.sizes
  }
  ph <- nh / n

  wls <- c()

  svdx <- svd(x)
  u <- svdx$u
  d <- svdx$d
  v <- svdx$v

  if (choose.dir == FALSE) {
    dir <- min(n, p)
  } else {
    theta <- d^2 / (d[1])^2 + 1
    loglik <- penalty <- rep(0, length(d))
    for (i in 1:length(d)) {
      if (i < length(d)) {
        loglik[i] <- sum(log(theta[(i + 1):length(d)]) + 1 - theta[(i + 1):length(d)])
      } else {
        loglik[i] <- 0
      }
      penalty[i] <- i * cn1 / sqrt(n)
    }
    BIC <- -loglik + penalty
    dir <- which.min(BIC)
  }

  w <- matrix(ncol = dir, nrow = length(nh))
  for (j in 1:dir) {
    for (i in 1:length(nh)) w[i, j] <- sum(u[, j] * (index == i)) / nh[i]
  }

  uut <- array(dim = c(length(nh), dir, dir))
  for (i in 1:length(nh)) uut[i, , ] <- (nh[i]) * (w[i, ] %*% t(w[i, ]))

  for (j in 1:p) wls[j] <- t(v[j, 1:dir]) %*% colSums(uut) %*% v[j, 1:dir]

  wls.sort <- sort(wls, decreasing = TRUE)

  loglik <- penalty <- rep(0, min(n, p))
  for (k in 1:min(n, p)) {
    temp_loglik <- sum(wls.sort[1:k])
    loglik[k] <- -log(temp_loglik)
    penalty[k] <- (log(n) + cn2 * log(p)) * k / max(n, p)
  }
  BIC <- loglik + penalty
  sel.k <- which.min(BIC)

  select <- order(wls, decreasing = TRUE)[1:sel.k]

  z <- x[, select]
  outsir <- Rdimtools::do.rsir(z, y, ndim = ndim, regmethod = "Ridge")

  betahat <- matrix(0, nrow = ncol(x), ncol = ndim)
  for (j in 1:ndim) {
    theta <- Re(outsir$projection[, j])
    betahat[select, j] <- theta
  }

  # ── Controle FDR (INCOMPLET, voir @section Attention ci-dessus) ─────────
  # t, Ta, alf, W, q doivent deja exister dans l'environnement appelant.
  tmin <- min(t[which(Ta <= alf)])
  Aplus <- which(W >= tmin)
  betaselct <- rep(0, p)
  betaselct[Aplus] <- 1

  FDP <- (length(which(betaselct[(q + 1):p] != 0))) / max(1, length(Aplus)) # E(FDR)
  AP <- (length(which(betaselct[1:q] != 0))) / max(1, q) # E(TDP)

  return(list(wls = wls, select = select, betahat = betahat))
}
