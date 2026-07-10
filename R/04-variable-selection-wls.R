#' Selection de variables WLS et pipeline CSIR-WLS
#'
#' @description
#' Fonctions du "score de levier pondere" (Weighted Least Squares
#' leverage score, WLS) utilise pour selectionner les variables actives
#' avant d'appliquer SIR ou Censored SIR. \code{\link{cen.wls.sir2}} est
#' le pipeline complet CSIR-WLS presente dans le papier : selection WLS
#' censuree -> double slicing -> Censored SIR sur les variables
#' retenues.
#'
#' @section Nettoyage effectue :
#' La version d'origine de \code{cen.wls()} contenait une redefinition
#' imbriquee accidentelle (copier-coller) de la fonction a l'interieur
#' d'elle-meme, placee dans une branche \code{if (!is.matrix(sigma.eta))}
#' qui n'est en pratique jamais executee (une multiplication matricielle
#' renvoie toujours une matrice). Le code etait donc inoffensif mais mort
#' et illisible ; il a ete retire ici, le comportement numerique de la
#' fonction est inchange.
#'
#' \code{wls.sir()} etait definie trois fois de maniere identique dans
#' vos scripts (dernier "gagnant" silencieux) ; une seule version est
#' gardee ici.
#'
#' @name variable-selection-wls
NULL

#' Score de levier pondere (WLS) + SIR, reponse non censuree
#'
#' Methode de selection de variables "Weighted Leverage Score" (WLS)
#' combinee a SIR : (1) on centre \code{x} et on decoupe \code{y} en
#' tranches ; (2) on calcule un score de levier pondere pour chaque
#' variable, base sur la SVD de \code{x} et les moyennes par tranche ;
#' (3) on choisit le nombre de variables actives via un critere BIC
#' applique aux scores tries ; (4) sur les variables selectionnees, on
#' applique un SIR ridge (\code{Rdimtools::do.rsir}) pour estimer les
#' directions beta.
#'
#' @param x Matrice de donnees n x p.
#' @param y Reponse (vecteur de taille n).
#' @param nslice Nombre de tranches pour SIR (defaut 10).
#' @param cn1,cn2 Constantes de penalite pour les criteres BIC.
#' @param choose.dir Si \code{TRUE}, choisit automatiquement le nombre de
#'   directions de la SVD via un BIC ; sinon utilise \code{min(n, p)}.
#' @param categorical \code{TRUE} si \code{y} est categorielle (les
#'   tranches sont alors les classes).
#' @param ndim Dimension de l'espace EDR recherche (nombre de directions
#'   beta).
#' @return Une liste : \code{wls} (scores de levier pondere), \code{select}
#'   (indices des variables retenues), \code{betahat} (matrice p x ndim,
#'   0 hors des variables selectionnees).
#' @export
wls.sir <- function(x, y, nslice = 10, cn1 = 0.1, cn2 = 1, choose.dir = FALSE,
                     categorical = FALSE, ndim = 1) {
  if (!requireNamespace("dr", quietly = TRUE)) {
    stop("Le package 'dr' est necessaire pour wls.sir().")
  }
  if (!requireNamespace("Rdimtools", quietly = TRUE)) {
    stop("Le package 'Rdimtools' est necessaire pour wls.sir().")
  }

  n <- dim(x)[1]
  p <- dim(x)[2]
  h <- nslice
  x <- scale(x, scale = FALSE)

  ## slice
  if (categorical == TRUE) {
    index <- as.numeric(factor(y))
    nh <- summary(factor(y))
    h <- length(nh)
  }
  if (categorical == FALSE) {
    slice <- dr::dr.slices.arc(y, h)
    index <- slice$slice.indicator # Slice Index
    nh <- slice$slice.sizes # Observations per Slice
  }
  ph <- nh / n

  ## calculate wls
  wls <- c() # Weighted Leverage Score

  svdx <- svd(x)
  u <- svdx$u
  d <- svdx$d
  v <- svdx$v

  if (choose.dir == FALSE) {
    dir <- min(n, p)
  } else {
    ## selection of d
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

  ## calculate WLS
  w <- matrix(ncol = dir, nrow = length(nh))
  for (j in 1:dir) {
    for (i in 1:length(nh)) w[i, j] <- sum(u[, j] * (index == i)) / nh[i]
  }

  uut <- array(dim = c(length(nh), dir, dir)) # UUT Array
  for (i in 1:length(nh)) uut[i, , ] <- (nh[i]) * (w[i, ] %*% t(w[i, ]))

  ## LEVERAGE SCORES
  for (j in 1:p) wls[j] <- t(v[j, 1:dir]) %*% colSums(uut) %*% v[j, 1:dir]

  ## BIC
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
  return(list(wls = wls, select = select, betahat = betahat))
}

#' Score de levier pondere (WLS) seul, sans SIR
#'
#' Version "allegee" de \code{\link{wls.sir}} : calcule le meme score de
#' levier pondere et la meme selection BIC des variables, mais s'arrete
#' apres la selection (n'applique pas de SIR ensuite).
#'
#' @inheritParams wls.sir
#' @return Vecteur des indices de colonnes de \code{x} selectionnees.
#' @export
wls <- function(x, y, nslice = 10, cn1 = 0.1, cn2 = 1, choose.dir = FALSE, categorical) {
  if (!requireNamespace("dr", quietly = TRUE)) {
    stop("Le package 'dr' est necessaire pour wls().")
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

  return(select)
}

#' Score WLS pour donnees censurees (selection de variables CSIR-WLS)
#'
#' Version DONNEES CENSUREES (survie) du score WLS : au lieu de trancher
#' \code{y} directement, on estime une fonction de survie conditionnelle
#' par noyau (\code{\link{kernel.est}} / \code{\link{cond.sur.est}} /
#' \code{\link{weight.est}}) pour ponderer les observations censurees,
#' puis on calcule la matrice \code{sigma.eta} (variance inter-tranches
#' ponderee) et enfin le score wls de chaque variable = v_j' sigma.eta v_j.
#' Le nombre de variables actives est choisi par un critere BIC.
#' C'est le coeur de l'etape de selection de CSIR-WLS.
#'
#' @param x Matrice de donnees n x p.
#' @param y Temps de survie (ou reponse) de taille n.
#' @param delta Indicateur de censure (1 = evenement observe, 0 = censure).
#' @param c Borne inferieure pour la fonction de survie estimee (defaut 0.05).
#' @param n.slice Nombre de tranches sur le temps de survie.
#' @param h Largeur de bande du noyau (si \code{NULL}, estimee
#'   automatiquement via la regle de Silverman sur la 1ere composante
#'   principale de \code{x}).
#' @param cn1,cn2 Constantes de penalite BIC (choix de d_hat puis de p0_hat).
#' @param choose.dir \code{TRUE}/\code{FALSE} : choisir automatiquement le
#'   nombre de directions de la SVD (d_hat) ou prendre \code{min(n, p)}.
#' @return Une liste : \code{n.sel} (nombre de variables selectionnees),
#'   \code{select} (indices retenus), \code{dir_used} (nombre de
#'   directions d_hat effectivement utilisees).
#' @export
cen.wls <- function(x, y, delta, c = 0.05, n.slice = 10, h = NULL,
                     cn1 = 0.1, cn2 = 1, choose.dir = FALSE) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  x <- scale(x, scale = FALSE)

  svdx <- svd(x, nu = min(n, p), nv = min(n, p))
  u <- svdx$u
  d <- svdx$d
  v <- svdx$v # d = valeurs singulieres (papier: lambda_i)

  if (is.null(h)) {
    beta_init <- prcomp(x, rank. = 1)$rotation[, 1]
    h <- bw.nrd(as.vector(x %*% beta_init))
  }

  # ── ETAPE 1 : choix de d_hat (Section 4.1, Eq. 9 du papier) ─────────────
  if (choose.dir == FALSE) {
    d_hat <- min(n, p)
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
    BIC_d <- -loglik + penalty
    d_hat <- which.min(BIC_d)
  }

  # ── TRONCATURE de u et v a d_hat colonnes ────────────────────────────────
  u <- u[, 1:d_hat, drop = FALSE]
  v <- v[, 1:d_hat, drop = FALSE]
  u.bar <- colMeans(u, na.rm = TRUE)

  e.hat.u <- matrix(0, nrow = (n.slice + 1), ncol = d_hat)

  ordr <- order(y, -delta)
  y <- y[ordr]
  u <- u[ordr, , drop = FALSE]
  delta <- delta[ordr]

  kernel.mtr <- kernel.est(u, h)
  s.hat <- cond.sur.est(kernel.mtr, y, c)

  t.slice <- slice.time(y, n.slice)
  m <- matrix(0, nrow = n.slice, ncol = d_hat)
  p.hat <- rep(1, (n.slice + 1))
  P <- rep(1, n.slice)
  nobs <- length(y)

  for (j in 1:(n.slice - 1)) {
    weight.vect <- rep(1, nobs)
    idx <- (1:nobs)
    idx.case1 <- idx[y < t.slice[(j + 1)] & delta == 0 & s.hat > c]
    idx.case2 <- idx[y >= t.slice[(j + 1)]]
    for (i in idx.case1) {
      weight.vect[i] <- weight.est(y[i], t.slice[(j + 1)], i, y, delta, kernel.mtr, s.hat)
    }
    p.hat[(j + 1)] <- sum(weight.vect[idx.case1], na.rm = TRUE) / nobs +
      sum(weight.vect[idx.case2], na.rm = TRUE) / nobs

    if (length(idx.case1) > 1) {
      ave.case1 <- apply(u[idx.case1, , drop = FALSE] * weight.vect[idx.case1], 2, sum, na.rm = TRUE) / nobs
    } else if (length(idx.case1) == 1) {
      ave.case1 <- u[idx.case1, ] * weight.vect[idx.case1] / nobs
    } else {
      ave.case1 <- 0
    }

    if (length(idx.case2) > 1) {
      ave.case2 <- apply(u[idx.case2, , drop = FALSE], 2, sum, na.rm = TRUE) / nobs
    } else if (length(idx.case2) == 1) {
      ave.case2 <- u[idx.case2, ] / nobs
    } else {
      ave.case2 <- 0
    }

    e.hat.u[(j + 1), ] <- ave.case1 + ave.case2
    diff_p <- p.hat[j] - p.hat[(j + 1)]
    if (abs(diff_p) < 1e-10) {
      m[j, ] <- rep(0, d_hat)
    } else {
      m[j, ] <- (e.hat.u[j, ] - e.hat.u[(j + 1), ]) / diff_p - u.bar
    }
    P[j] <- max(0, p.hat[j] - p.hat[(j + 1)])
  }

  if (abs(p.hat[n.slice]) < 1e-10) {
    m[n.slice, ] <- rep(0, d_hat)
  } else {
    m[n.slice, ] <- e.hat.u[n.slice, ] / p.hat[n.slice] - u.bar
  }
  P[n.slice] <- p.hat[n.slice]

  mp <- m * sqrt(P)
  if (!is.matrix(mp)) mp <- matrix(mp, nrow = length(P))

  # sigma.eta est d x d
  sigma.eta <- t(mp) %*% mp
  if (!is.matrix(sigma.eta)) sigma.eta <- matrix(sigma.eta, nrow = d_hat, ncol = d_hat)

  # ── ETAPE 2 : calcul de omega_j sur les d_hat directions retenues ───────
  wls <- numeric(p)
  for (j in 1:p) {
    v_j <- matrix(v[j, ], ncol = 1) # v deja tronque a d_hat colonnes
    wls[j] <- as.numeric(t(v_j) %*% sigma.eta %*% v_j)
  }

  # ── ETAPE 3 : choix de p0_hat (Section 4.2, Eq. 10 du papier) ────────────
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

  return(list(n.sel = sel.k, select = select, dir_used = d_hat))
}

#' Pipeline CSIR-WLS (variante 1 : selection WLS non censuree + double slicing)
#'
#' Pipeline complet pour donnees censurees : (1) selection de variables
#' par le score WLS classique (\code{\link{wls.sir}}, applique a
#' \code{y} comme s'il n'etait pas censure) ; (2) sur les variables
#' retenues, calcul de directions initiales via
#' \code{\link{double.slice2}} (2 directions) ; (3) application de
#' \code{\link{cen.sir}} pour obtenir les directions finales beta,
#' replacees dans un vecteur de taille p (0 hors variables
#' selectionnees).
#'
#' @note Pour le pipeline utilise dans le papier (selection WLS
#'   *censuree*), voir plutot \code{\link{cen.wls.sir2}}.
#'
#' @param x,y,delta Donnees et censure.
#' @param n.slice1,n.slice0 Nombre de tranches (non-censure / censure)
#'   pour \code{\link{double.slice2}}.
#' @param n.slice Nombre de tranches pour le score WLS et pour
#'   \code{\link{cen.sir}}.
#' @param cn1,cn2 Constantes de penalite BIC.
#' @param choose.dir,categorical,ndim Comme dans \code{\link{wls.sir}}.
#' @return Une liste : \code{wls} (scores), \code{select} (variables
#'   retenues), \code{betahat} (matrice p x ndim), \code{n.sel} (nombre
#'   de variables retenues).
#' @export
cen.wls.sir <- function(x, y, delta, n.slice1, n.slice0, n.slice = 10, cn1 = 0.1,
                         cn2 = 1, choose.dir = FALSE, categorical = FALSE, ndim = 1) {
  if (!requireNamespace("dr", quietly = TRUE)) {
    stop("Le package 'dr' est necessaire pour cen.wls.sir().")
  }

  n <- dim(x)[1]
  p <- dim(x)[2]
  h <- n.slice
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
  ds <- double.slice2(y, delta, z, n.slice1, n.slice0)
  joint.edrs <- edr.n(ds, 2)
  sir <- cen.sir(y, delta, z, n.slice, joint.edrs, h, c = 0.05)

  betahat <- matrix(0, nrow = ncol(x), ncol = ndim)
  for (j in 1:ndim) {
    theta <- sir$evec[, j]
    betahat[select, j] <- theta
  }
  return(list(wls = wls, select = select, betahat = betahat, n.sel = sel.k))
}

#' Pipeline CSIR-WLS (variante utilisee dans le papier)
#'
#' Meme pipeline que \code{\link{cen.wls.sir}}, mais la selection de
#' variables utilise directement \code{\link{cen.wls}} (le score WLS
#' "version censuree" avec noyau/survie conditionnelle, plutot que le
#' score WLS "non censure"). La largeur de bande h pour
#' \code{\link{cen.sir}} est estimee automatiquement (regle de
#' Silverman) a partir de la 1ere composante principale de \code{x}.
#' C'est la fonction principale de la methode CSIR-WLS presentee dans
#' le papier.
#'
#' @param x,y,delta Donnees et censure.
#' @param n.slice1,n.slice0 Nombre de tranches pour \code{\link{double.slice2}}.
#' @param n.slice Nombre de tranches pour \code{\link{cen.wls}} /
#'   \code{\link{cen.sir}}.
#' @param cn1,cn2,choose.dir Comme dans \code{\link{cen.wls}}.
#' @param ndim Nombre de directions EDR souhaitees.
#' @param c Plancher de la survie conditionnelle estimee.
#' @return Une liste : \code{beta.hat} (matrice p x ndim), \code{eval.hat}
#'   (valeurs propres replacees dans un vecteur de taille p),
#'   \code{select} (variables retenues), \code{n.sel} (nombre de
#'   variables retenues).
#' @export
cen.wls.sir2 <- function(x, y, delta, n.slice1, n.slice0, n.slice = 10, cn1 = 0.1,
                          cn2 = 1, choose.dir = FALSE, ndim = 1, c = 0.05) {
  cwlsir <- cen.wls(x, as.vector(y), delta, c = 0.05, n.slice = n.slice, cn1 = 0.1, cn2 = 1, choose.dir = choose.dir)
  select <- cwlsir$select
  n.sel <- cwlsir$n.sel
  Z <- x[, select]
  ds <- double.slice2(y, delta, Z, n.slice1, n.slice0)
  joint.edrs <- edr.n(ds, 2)
  beta_init <- prcomp(x, rank. = 1)$rotation[, 1]
  # Rapide et justifie theoriquement -> Silverman
  h <- bw.nrd(as.vector(x %*% beta_init))
  # Plus precis si besoin -> Sheather-Jones : bw.SJ(as.vector(x %*% beta_init))
  censir <- cen.sir(y, delta, Z, n.slice, joint.edrs, h, c)
  beta.hat <- matrix(0, nrow = ncol(x), ncol = ndim)
  eval.hat <- rep(0, ncol(x))
  eval.hat[select] <- censir$eval
  for (j in 1:ndim) {
    theta <- censir$evec[, j]
    beta.hat[select, j] <- theta
  }
  rownames(beta.hat) <- colnames(x)
  return(list(beta.hat = beta.hat, eval.hat = eval.hat, select = select, n.sel = n.sel))
}

#' Choix de la dimension EDR d par un critere BIC
#'
#' Estime la dimension structurelle d en appliquant un critere BIC aux
#' valeurs propres renvoyees par \code{\link{cen.sir}} apres l'etape de
#' screening CSIR-WLS.
#'
#' @param evals Vecteur de valeurs propres (ordre quelconque, peut
#'   contenir des valeurs complexes/NA/negatives : elles sont nettoyees
#'   avant le calcul).
#' @param n Taille d'echantillon.
#' @param cn1 Constante de penalite (defaut 0.1, comme Zhong et al.).
#' @return \code{d_hat} : nombre estime de directions EDR (entier >= 1).
#' @export
select_d_bic <- function(evals, n, cn1 = 0.1) {
  # ── Nettoyage defensif ────────────────────────────────────────────────
  evals <- as.numeric(evals) # force numeric (au cas ou complexe)
  evals <- Re(evals) # partie reelle si complexe
  evals <- evals[!is.na(evals)] # retirer NA
  evals <- evals[is.finite(evals)] # retirer Inf/-Inf
  evals <- evals[evals > 1e-10] # garder seulement positifs non-nuls
  evals <- sort(evals, decreasing = TRUE) # ordre decroissant

  # si pas assez de valeurs propres exploitables -> d = 1 par defaut
  if (length(evals) <= 1) {
    return(1)
  }

  K <- length(evals)
  theta <- evals / evals[1] + 1 # theta_i = lambda_i/lambda_1 + 1 >= 1

  loglik <- numeric(K)
  penalty <- numeric(K)

  for (i in 1:K) {
    if (i < K) {
      terms <- log(theta[(i + 1):K]) + 1 - theta[(i + 1):K]
      loglik[i] <- sum(terms[is.finite(terms)])
    } else {
      loglik[i] <- 0
    }
    penalty[i] <- i * cn1 / sqrt(n)
  }

  BIC <- -loglik + penalty
  d_hat <- which.min(BIC)

  return(as.integer(d_hat))
}

#' Choix de la dimension EDR d par le test sequentiel du chi2 (Li, 1991)
#'
#' Alternative "classique" a \code{\link{select_d_bic}} : applique le test
#' sequentiel du chi2 renvoye par \code{\link{sir.test}} (\code{"D=0 vs
#' D>=1"}, \code{"D=1 vs D>=2"}, etc.) a un objet \code{"csir"} ou
#' \code{"ds"}, et s'arrete a la premiere hypothese non rejetee.
#'
#' D'apres les etudes de recuperation de dimension menees en simulation
#' (index multiple, d vrai = 2) : ce test tend a \strong{sous-estimer} d
#' des que l'une des directions EDR n'entre dans le modele que via une
#' forme non-lineaire/peu asymetrique (ex : un terme en ratio ou en carre
#' dans le predicteur), et perd encore plus de puissance quand le taux de
#' censure augmente. Comparer systematiquement a \code{\link{select_d_bic}}
#' sur les memes donnees avant de trancher.
#'
#' @param csir_obj Objet de classe \code{"csir"} (sortie de
#'   \code{\link{cen.sir}}) ou \code{"ds"} (sortie de
#'   \code{\link{double.slice}} / \code{\link{double.slice2}}).
#' @param alpha Seuil de signification pour "ne pas rejeter" (defaut 0.01,
#'   plus conservateur que le defaut usuel de 0.05, pour eviter de
#'   surestimer d sur de petits echantillons).
#' @return \code{d_hat} : nombre estime de directions EDR (entier >= 0).
#' @seealso \code{\link{select_d_bic}} pour le critere BIC-type utilise en
#'   parallele dans le papier (Section 5.3.2).
#' @export
estimate_dim_chisq <- function(csir_obj, alpha = 0.01) {
  test_tab <- sir.test(csir_obj)
  first_ns <- which(test_tab$p.value > alpha)[1]
  d_hat <- if (is.na(first_ns)) nrow(test_tab) - 1 else first_ns - 1
  d_hat
}

#' Pipeline CSIR-WLS avec choix automatique de la dimension d
#'
#' Version robuste de \code{\link{cen.wls.sir2}} qui estime elle-meme la
#' dimension EDR d via \code{\link{select_d_bic}}, au lieu de la fixer a
#' l'avance.
#'
#' @param x,y,delta Donnees et censure.
#' @param n.slice1,n.slice0,n.slice Comme dans \code{\link{cen.wls.sir2}}.
#' @param cn1,cn2 Constantes de penalite BIC.
#' @param c Plancher de la survie conditionnelle estimee.
#' @param d_max Dimension maximale autorisee.
#' @return Une liste : \code{beta.hat} (matrice p x d_hat), \code{eval.hat},
#'   \code{select}, \code{n.sel}, \code{d_hat} (dimension estimee).
#' @export
cen_wls_sir2_auto_d <- function(x, y, delta, n.slice1, n.slice0,
                                 n.slice = 10, cn1 = 0.1, cn2 = 1,
                                 c = 0.05, d_max = 10) {
  n <- nrow(x)
  p <- ncol(x)

  # Step 1: CWLS screening
  cwlsir <- cen.wls(x, as.vector(y), delta,
    c = c, n.slice = n.slice,
    cn1 = cn1, cn2 = cn2, choose.dir = FALSE
  )
  select <- cwlsir$select
  n.sel <- cwlsir$n.sel
  Z <- x[, select, drop = FALSE]

  # Step 2: double slicing
  ds <- double.slice2(y, delta, Z, n.slice1, n.slice0)
  joint.edrs <- edr.n(ds, 2)

  # Step 3: bandwidth
  beta_init <- prcomp(x, rank. = 1)$rotation[, 1]
  h <- bw.nrd(as.vector(x %*% beta_init))

  # Step 4: run cen.sir -> returns all eigenvalues
  censir <- cen.sir(y, delta, Z, n.slice, joint.edrs, h, c)

  # Step 5: estimate d via BIC -- defensive extraction
  evals_raw <- censir$eval
  evals_raw <- as.numeric(Re(evals_raw))
  evals_pos <- evals_raw[evals_raw > 1e-10]

  if (length(evals_pos) == 0) {
    d_hat <- 1 # fallback
  } else {
    d_hat <- select_d_bic(evals_pos, n = n, cn1 = cn1)
  }

  # cap d_hat
  d_hat <- as.integer(min(d_hat, d_max, length(select), ncol(censir$evec)))
  d_hat <- max(1L, d_hat)

  # Step 6: build beta.hat
  beta.hat <- matrix(0, nrow = p, ncol = d_hat)
  eval.hat <- rep(0, p)
  eval.hat[select] <- censir$eval

  for (j in 1:d_hat) {
    beta.hat[select, j] <- Re(censir$evec[, j])
  }
  rownames(beta.hat) <- colnames(x)

  return(list(
    beta.hat = beta.hat, eval.hat = eval.hat,
    select = select, n.sel = n.sel, d_hat = d_hat
  ))
}

#' Ajuster le nombre de colonnes d'une matrice de directions
#'
#' Complete par des zeros (si \code{d < d_target}) ou tronque (si
#' \code{d > d_target}) une matrice de directions beta pour qu'elle ait
#' exactement \code{d_target} colonnes. Utile pour comparer des beta
#' estimes de dimensions differentes (ex : d_hat variable) a un beta
#' vrai de dimension fixe.
#'
#' @param B Matrice (ou vecteur) de directions.
#' @param d_target Nombre de colonnes cible.
#' @return Matrice p x d_target.
#' @export
pad_beta <- function(B, d_target) {
  B <- as.matrix(Re(B))
  p <- nrow(B)
  d <- ncol(B)
  if (d == d_target) {
    return(B)
  }
  if (d < d_target) {
    return(cbind(B, matrix(0, p, d_target - d)))
  }
  return(B[, 1:d_target, drop = FALSE])
}
