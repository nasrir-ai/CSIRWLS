#' Coeur de la methode Censored SIR (CSIR)
#'
#' @description
#' Les fonctions de ce fichier implementent le "Censored SIR" (SIR pour
#' donnees censurees) tel qu'utilise dans le papier CSIR-WLS : estimation
#' a noyau de la survie conditionnelle, ponderation des observations
#' censurees (IPCW-like), decoupage en tranches ("slicing") tenant compte
#' de la censure, et resolution du probleme aux valeurs propres
#' generalise pour obtenir les directions EDR (Effective Dimension
#' Reduction). Toutes ces fonctions viennent de \code{myfunctions.R}.
#'
#' @name core-censored-sir
NULL

#' Matrice de poids a noyau gaussien entre observations projetees
#'
#' \code{kernel.x[i, j] = exp(-0.5 * distance_normalisee(Bx_i, Bx_j)^2)}.
#'
#' @param Bx Matrice n x r des directions projetees (ex : \code{x \%*\% beta}).
#' @param h Largeur de bande du noyau (meme largeur pour chaque colonne).
#' @return Matrice n x n des poids du noyau.
#' @export
kernel.est <- function(Bx, h) {
  n <- nrow(Bx)
  r <- ncol(Bx)
  b <- rep(h, r)

  kernel.x <- matrix(rep(0, n * n), nrow = n)
  for (i in 1:n) {
    Bx.temp <- t(matrix(rep(Bx[i, ], n), nrow = r, ncol = n))
    d2 <- ((Bx - Bx.temp) / b)^2
    ds <- apply(d2, 1, sum)
    fh <- exp(-0.5 * (ds))
    kernel.x[i, ] <- fh
  }
  kernel.x
}

#' Survie conditionnelle estimee a noyau
#'
#' Estime, pour chaque observation i, la probabilite de survie
#' conditionnelle S_hat(y_i | x_i) par une methode a noyau (moyenne
#' ponderee des indicatrices y_j > y_i), avec un plancher \code{c} (pour
#' eviter des divisions par des valeurs trop proches de 0 plus loin dans
#' les calculs).
#'
#' @param kernel.mtr Matrice de poids du noyau, sortie de
#'   \code{\link{kernel.est}}.
#' @param y Temps (survie ou reponse), deja tries.
#' @param c Valeur plancher pour \code{s.hat}.
#' @return Vecteur \code{s.hat} de longueur n (survie conditionnelle estimee).
#' @export
cond.sur.est <- function(kernel.mtr, y, c) {
  n <- length(y)
  s.hat <- rep(0, n)
  for (i in 1:n) {
    yi <- (1:n)[y > y[i]]
    fh <- kernel.mtr[i, ]
    if (length(yi) > 0) {
      s.hat[i] <- sum(fh[yi], na.rm = TRUE) / sum(fh, na.rm = TRUE)
      s.hat[i] <- max(s.hat[i], c, na.rm = TRUE)
    } else {
      s.hat[i] <- c
    }
  }
  s.hat
}

#' Poids de correction type IPCW pour une observation censuree
#'
#' Calcule le poids de correction (type Kaplan-Meier local / IPCW) pour
#' une observation censuree situee entre les temps \code{t1} et
#' \code{t2}, a partir du risque cumule estime a noyau. Ce poids sert a
#' "completer" l'information manquante des observations censurees dans
#' \code{\link{cen.wls}} / \code{\link{cen.sir}}.
#'
#' @param t1,t2 Bornes de temps (\code{t1 < t2}, sinon erreur).
#' @param index.x Indice de l'observation courante.
#' @param y,delta Temps et indicateur de censure (vecteurs).
#' @param kernel.mtr Matrice de poids du noyau.
#' @param s.hat Survie conditionnelle estimee, sortie de
#'   \code{\link{cond.sur.est}}.
#' @return Poids scalaire (entre 0 et 1).
#' @export
weight.est <- function(t1, t2, index.x, y, delta, kernel.mtr, s.hat) {
  if (t1 >= t2) {
    cat(paste(t1, "\t", t2, "\n"))
    stop("t1 >= t2!")
  }

  # equation (5.2)
  fhatx0 <- mean(kernel.mtr[index.x, ], na.rm = TRUE)
  lambdahat.ttx.temp <- kernel.mtr[index.x, ] * (t1 < y) * (y < t2) * delta
  lambdahat.ttx <- mean((lambdahat.ttx.temp / s.hat), na.rm = TRUE) / fhatx0

  # equation (5.1)
  weight.ttx <- exp(-lambdahat.ttx)
  weight.ttx
}

#' Decoupage d'un vecteur de temps en tranches de taille egale
#'
#' Decoupe un vecteur de temps tries \code{t} en \code{h} tranches
#' (slices) de taille a peu pres egale, et renvoie les h+1 bornes de
#' temps correspondantes.
#'
#' @param t Vecteur de temps (deja trie).
#' @param h Nombre de tranches souhaite.
#' @return Vecteur \code{t.slice} de longueur h+1 (bornes des tranches).
#' @export
slice.time <- function(t, h) {
  lt <- length(t)
  n1 <- ceiling(lt / h)
  n2 <- floor(lt / h)
  if (n1 == n2) {
    n <- lt / h
    index <- rep(n, h)
  } else if (n1 > n2) {
    h1 <- lt %% h
    h2 <- h - h1
    index <- rep(c(n1, n2), c(h1, h2))
  }
  t.slice <- numeric(h + 1) # (h+1) time points
  t.slice[1] <- min(t) - 0.000001
  end <- 0
  for (i in 1:h) {
    bgn <- end + 1
    end <- end + index[i]
    t.slice[i + 1] <- t[end]
  }
  t.slice
}

#' Censored SIR (SIR pour donnees censurees)
#'
#' Projette \code{x} sur des directions initiales (\code{joint.edrs}),
#' estime le noyau et la survie conditionnelle, calcule les moyennes
#' ponderees par tranche puis la matrice \code{sigma.eta}, et resout un
#' probleme aux valeurs propres generalise (\code{sigma.eta}, sigma.x)
#' pour obtenir les directions EDR triees par valeur propre decroissante.
#'
#' @param y Temps de survie (ou reponse).
#' @param delta Indicateur de censure (1 = evenement, 0 = censure).
#' @param x Matrice de donnees (deja reduite aux variables voulues).
#' @param n.slice Nombre de tranches.
#' @param joint.edrs Directions initiales utilisees pour projeter x
#'   (\code{Bx = x \%*\% joint.edrs}), typiquement issues de
#'   \code{\link{double.slice}} + \code{\link{edr.n}}.
#' @param h Largeur de bande du noyau.
#' @param c Plancher pour la survie conditionnelle estimee.
#' @return Un objet de classe \code{"csir"} (liste) contenant \code{y},
#'   \code{x}, \code{delta}, \code{eval} (valeurs propres), \code{evec}
#'   (vecteurs/directions EDR), \code{nslice}, \code{c}, \code{h}.
#' @export
cen.sir <- function(y, delta, x, n.slice, joint.edrs, h, c) {
  # order the data by survival time first
  ordr <- order(y, -delta)
  y <- y[ordr]
  x <- x[ordr, ]
  delta <- delta[ordr]

  # first, estimate the kernel matrix and conditional survival function
  Bx <- x %*% joint.edrs
  kernel.mtr <- kernel.est(Bx, h)
  s.hat <- cond.sur.est(kernel.mtr, y, c)

  # secondly, calculate the slice mean m_h
  t.slice <- slice.time(y, n.slice)
  m <- matrix(0, nrow = n.slice, ncol = ncol(x))
  p.hat <- rep(1, (n.slice + 1))
  p <- rep(1, n.slice)
  e.hat.x <- matrix(0, nrow = (n.slice + 1), ncol = ncol(x))
  x.bar <- apply(x, 2, mean, na.rm = TRUE)
  e.hat.x[1, ] <- x.bar
  nobs <- length(y)

  for (j in 1:(n.slice - 1)) {
    weight.vect <- rep(1, nobs)
    idx <- (1:nobs)
    idx.case1 <- idx[y < t.slice[(j + 1)] & delta == 0 & s.hat > c]
    idx.case2 <- idx[y >= t.slice[(j + 1)]]

    for (i in idx.case1) {
      weight.vect[i] <- weight.est(
        y[i], t.slice[(j + 1)], i, y, delta,
        kernel.mtr, s.hat
      )
    }

    p.hat[(j + 1)] <- sum(weight.vect[idx.case1], na.rm = TRUE) / nobs
    p.hat[(j + 1)] <- p.hat[(j + 1)] + sum(weight.vect[idx.case2], na.rm = TRUE) / nobs

    if (length(idx.case1) > 1) {
      x.case1 <- x[idx.case1, ] * weight.vect[idx.case1]
      ave.case1 <- apply(x.case1, 2, sum, na.rm = TRUE) / nobs
    } else if (length(idx.case1) == 1) {
      ave.case1 <- x[idx.case1, ] * weight.vect[idx.case1] / nobs
    } else if (length(idx.case1) == 0) {
      ave.case1 <- 0
    }

    if (length(idx.case2) > 1) {
      ave.case2 <- apply(x[idx.case2, ], 2, sum, na.rm = TRUE) / nobs
    } else if (length(idx.case2) == 1) {
      ave.case2 <- x[idx.case2, ] / nobs
    } else if (length(idx.case2) == 0) {
      ave.case2 <- 0
    }

    e.hat.x[(j + 1), ] <- ave.case1 + ave.case2

    m[j, ] <- (e.hat.x[j, ] - e.hat.x[(j + 1), ]) / (p.hat[j] - p.hat[(j + 1)]) - x.bar
    p[j] <- max(0, p.hat[j] - p.hat[(j + 1)])
  }

  m[n.slice, ] <- e.hat.x[n.slice, ] / p.hat[n.slice] - x.bar
  p[n.slice] <- p.hat[n.slice]

  # at last, compute the SIR direction
  mp <- m * sqrt(p)
  sigma.eta <- t(mp) %*% mp
  sigma.x <- cov1(x, use = "pair")
  eign <- eigen.decomp(sigma.eta, sigma.x)
  ordr <- order(eign$values, decreasing = TRUE)
  e.val <- eign$values[ordr]
  e.vec <- eign$vectors[, ordr]

  if (is.null(colnames(x))) {
    x.names <- paste("X", 1:ncol(x), sep = "")
  } else {
    x.names <- colnames(x)
  }
  dimnames(e.vec) <- list(x.names, paste("Dir", 1:ncol(e.vec), sep = ""))

  csir <- list(y = y, x = x, delta = delta, eval = e.val, evec = e.vec, nslice = n.slice)
  csir$c <- c
  csir$h <- h
  class(csir) <- "csir"
  csir
}

#' Extraire les n premieres directions EDR
#'
#' @param object Un objet contenant \code{$evec} (typiquement issu de
#'   \code{\link{cen.sir}}, \code{\link{double.slice}} ou
#'   \code{\link{double.slice2}}).
#' @param n Nombre de directions a extraire.
#' @return Matrice p x n des n premieres directions.
#' @export
edr.n <- function(object, n) {
  object$evec[, 1:n]
}

#' Selectionner les directions EDR expliquant un pourcentage de variance donne
#'
#' Selectionne le plus petit nombre j de directions EDR dont la somme
#' des valeurs propres explique au moins \code{percent} de la variance
#' totale.
#'
#' @param object Objet contenant \code{$eval} (valeurs propres) et
#'   \code{$evec}.
#' @param percent Proportion cible de variance expliquee (entre 0 et 1).
#' @return Matrice des j premiers vecteurs propres (directions EDR).
#' @export
edr.percent <- function(object, percent) {
  if (percent > 1 || percent < 0) {
    stop("percent should be between 0 and 1\n")
  }

  values <- object$eval
  vectors <- object$evec

  j <- 1
  while (sum(values[1:j]) / sum(values) < percent) {
    j <- j + 1 # take the first j eigenvalues
  }
  pct <- round(sum(values[1:j]) / sum(values), 2)
  cat(paste("eigen.value =", "\n"))
  cat(signif(values, 3))
  cat("\n")
  cat(paste("First", j, "eigen values explain", pct, "variance", "\n"))
  cat("\n")

  vectors.mtr <- vectors[, 1:j]
  vectors.mtr
}

#' Test asymptotique du chi2 pour le nombre de directions EDR
#'
#' Test sequentiel "D=0 vs D>=1", "D=1 vs D>=2", etc., base sur les
#' valeurs propres de l'objet SIR/csir/ds.
#'
#' @param object Objet de classe \code{"csir"} ou \code{"ds"}.
#' @param nd Nombre max de directions a tester (defaut = toutes).
#' @return Un \code{data.frame} avec les colonnes \code{Chisq}, \code{df},
#'   \code{p.value} (une ligne par test).
#' @export
sir.test <- function(object, nd = length(object$eval)) {
  if (class(object) != "csir" && class(object) != "ds") {
    stop("wrong class of object in sir.test\n")
  }
  e <- sort(object$eval)
  p <- length(object$eval)
  n <- length(object$y)
  l <- min(p, nd)
  chisq <- numeric(l)
  df <- numeric(l)
  p.val <- numeric(l)
  for (i in 0:(l - 1)) {
    j <- i + 1
    chisq[j] <- n * sum(e[1:(p - i)])
    df[j] <- (p - i) * (object$nslice - i - 1)
    p.val[j] <- 1 - pchisq(chisq[j], df[j])
  }
  test <- data.frame(Chisq = chisq, df = df, p.value = p.val)
  rnames <- character(l)
  for (i in 1:l) {
    rnames[i] <- paste("D=", i - 1, " vs. ", "D>=", i, sep = "")
  }
  rownames(test) <- rnames
  test
}

#' Double Slicing pour SIR censure (version originale)
#'
#' Les observations non censurees (\code{delta = 1}) et censurees
#' (\code{delta = 0}) sont tranchees separement (\code{n.slice1}
#' tranches pour les non-censurees, \code{n.slice0} pour les censurees).
#' On calcule la covariance moyenne intra-tranche, puis
#' \code{sigma.eta = sigma.x - ave.covx}, et on resout le probleme aux
#' valeurs propres generalise pour obtenir les directions EDR initiales
#' (utilisees ensuite comme point de depart de \code{\link{cen.sir}}).
#'
#' @note Cette version plante si une tranche est vide ou n'a qu'une
#'   observation. Pour un usage robuste, preferer
#'   \code{\link{double.slice2}} (meme algorithme, avec des gardes-fous).
#'
#' @param y Temps de survie (ou reponse).
#' @param delta Indicateur de censure.
#' @param x Matrice de donnees.
#' @param n.slice1 Nombre de tranches pour les non-censurees (delta = 1).
#' @param n.slice0 Nombre de tranches pour les censurees (delta = 0).
#' @return Un objet de classe \code{"ds"} (liste) : \code{y}, \code{x},
#'   \code{delta}, \code{eval}, \code{evec}, \code{nslice}, \code{nslice0},
#'   \code{nslice1}.
#' @export
double.slice <- function(y, delta, x, n.slice1, n.slice0) {
  slice.idx <- function(y, h, m) {
    ly <- length(y)
    if (m < max(y)) {
      stop("y is not a valid index vector")
    }
    n1 <- ceiling(ly / h)
    n2 <- floor(ly / h)
    split <- numeric(m)
    if (n1 == n2) {
      n <- ly / h
      index <- rep(n, h)
    } else if (n1 > n2) {
      h1 <- ly %% h
      h2 <- h - h1
      index <- rep(c(n1, n2), c(h1, h2))
    }
    end <- 0
    for (i in 1:h) {
      bgn <- end + 1
      end <- end + index[i]
      split[y[bgn:end]] <- i
    }
    split
  }

  # order the data by survival time first
  ordr <- order(y, -delta)
  y <- y[ordr]
  x <- x[ordr, ]
  delta <- delta[ordr]

  ds.idx <- numeric(nrow(x))

  y.idx0 <- (1:nrow(x))[delta == 0]
  y.idx0 <- y.idx0[order(y[delta == 0])]
  ds.idx0 <- slice.idx(y.idx0, n.slice0, nrow(x))

  y.idx1 <- (1:nrow(x))[delta == 1]
  y.idx1 <- y.idx1[order(y[delta == 1])]
  ds.idx1 <- slice.idx(y.idx1, n.slice1, nrow(x))

  ds.idx[ds.idx0 > 0] <- ds.idx0[ds.idx0 > 0]
  ds.idx[ds.idx1 > 0] <- ds.idx1[ds.idx1 > 0] + n.slice0

  n.per.slice <- as.vector(table(ds.idx))
  ave.covx <- matrix(0, ncol(x), ncol(x))
  probs <- n.per.slice / nrow(x)
  for (i in 1:length(n.per.slice)) {
    x.slice <- as.matrix(x[(1:nrow(x))[ds.idx == i], ])
    cov.slice <- cov1(x.slice, use = "pair")
    if (i == 1) {
      ave.covx <- probs[i] * cov.slice
    } else {
      ave.covx <- ave.covx + probs[i] * cov.slice
    }
  }

  sigma.x <- cov1(x, use = "pair")
  sigma.eta <- sigma.x - ave.covx

  eign <- eigen.decomp(sigma.eta, sigma.x)
  ordr <- order(eign$values, decreasing = TRUE)
  e.val <- eign$values[ordr]
  e.vec <- eign$vectors[, ordr]

  if (is.null(colnames(x))) {
    x.names <- paste("X", 1:ncol(x), sep = "")
  } else {
    x.names <- colnames(x)
  }
  dimnames(e.vec) <- list(x.names, paste("Dir", 1:ncol(e.vec), sep = ""))

  n.s <- n.slice1 + n.slice0
  ds <- list(y = y, x = x, delta = delta, eval = e.val, evec = e.vec, nslice = n.s)
  ds$nslice0 <- n.slice0
  ds$nslice1 <- n.slice1

  class(ds) <- "ds"
  ds
}

#' Double Slicing pour SIR censure (version robuste)
#'
#' Identique a \code{\link{double.slice}}, mais avec des gardes-fous pour
#' eviter les erreurs quand une tranche est vide ou ne contient qu'une
#' seule observation (dans ce cas, la covariance de la tranche est
#' remplacee par une matrice de zeros au lieu de faire planter
#' \code{cov()}). C'est la version utilisee en pratique par
#' \code{\link{cen.wls.sir2}}.
#'
#' @inheritParams double.slice
#' @return Identique a \code{\link{double.slice}}.
#' @export
double.slice2 <- function(y, delta, x, n.slice1, n.slice0) {
  slice.idx <- function(y, h, m) {
    ly <- length(y)
    if (m < max(y)) {
      stop("y is not a valid index vector")
    }
    n1 <- ceiling(ly / h)
    n2 <- floor(ly / h)
    split <- numeric(m)
    if (n1 == n2) {
      n <- ly / h
      index <- rep(n, h)
    } else if (n1 > n2) {
      h1 <- ly %% h
      h2 <- h - h1
      index <- rep(c(n1, n2), c(h1, h2))
    }
    end <- 0
    for (i in 1:h) {
      bgn <- end + 1
      end <- end + index[i]
      split[y[bgn:end]] <- i
    }
    split
  }

  cov1.safe <- function(x, use = "pair") {
    if (nrow(x) == 0) {
      return(matrix(0, ncol(x), ncol(x))) # Return a zero matrix if x is empty
    }
    n <- nrow(x)
    c <- cov(x, use = use) * (n - 1) / n
    c
  }

  # Order data by survival time
  ordr <- order(y, -delta)
  y <- y[ordr]
  x <- x[ordr, ]
  delta <- delta[ordr]

  ds.idx <- numeric(nrow(x))

  y.idx0 <- (1:nrow(x))[delta == 0]
  y.idx0 <- y.idx0[order(y[delta == 0])]
  ds.idx0 <- slice.idx(y.idx0, n.slice0, nrow(x))

  y.idx1 <- (1:nrow(x))[delta == 1]
  y.idx1 <- y.idx1[order(y[delta == 1])]
  ds.idx1 <- slice.idx(y.idx1, n.slice1, nrow(x))

  ds.idx[ds.idx0 > 0] <- ds.idx0[ds.idx0 > 0]
  ds.idx[ds.idx1 > 0] <- ds.idx1[ds.idx1 > 0] + n.slice0

  n.per.slice <- as.vector(table(ds.idx))
  probs <- n.per.slice / nrow(x)

  # Initialize ave.covx
  ave.covx <- matrix(0, ncol(x), ncol(x))

  for (i in 1:length(n.per.slice)) {
    idx <- (1:nrow(x))[ds.idx == i]
    if (length(idx) > 1) {
      x.slice <- as.matrix(x[idx, ])
      cov.slice <- cov1.safe(x.slice, use = "pair")
    } else {
      cov.slice <- matrix(0, ncol(x), ncol(x)) # zero matrix if slice empty/singleton
    }
    ave.covx <- ave.covx + probs[i] * cov.slice
  }

  sigma.x <- cov1.safe(x, use = "pair")
  sigma.eta <- sigma.x - ave.covx

  eign <- eigen.decomp(sigma.eta, sigma.x)
  ordr <- order(eign$values, decreasing = TRUE)
  e.val <- eign$values[ordr]
  e.vec <- eign$vectors[, ordr]

  if (is.null(colnames(x))) {
    x.names <- paste("X", 1:ncol(x), sep = "")
  } else {
    x.names <- colnames(x)
  }
  dimnames(e.vec) <- list(x.names, paste("Dir", 1:ncol(e.vec), sep = ""))

  n.s <- n.slice1 + n.slice0
  ds <- list(y = y, x = x, delta = delta, eval = e.val, evec = e.vec, nslice = n.s)
  ds$nslice0 <- n.slice0
  ds$nslice1 <- n.slice1

  class(ds) <- "ds"
  ds
}

#' Affichage d'un objet "ds" (double slicing)
#'
#' Methode \code{print} pour un objet issu de \code{\link{double.slice}}
#' / \code{\link{double.slice2}} : nombre d'observations, de censures,
#' de predicteurs, de tranches, valeurs propres (+ variance cumulee),
#' test du chi2 (\code{\link{sir.test}}) et vecteurs propres.
#'
#' @param object Objet de classe \code{"ds"}.
#' @param digits Nombre de decimales a afficher.
#' @param ... Non utilise.
#' @return \code{object}, invisiblement.
#' @export
print.ds <- function(object, digits = max(3, getOption("digits") - 3), ...) {
  cat("Double Slicing of Survival time\n\n")
  cat(paste("Number of observation:           ", length(object$y), "\n"))
  n.censor <- length(object$y[object$delta == 0])
  cat(paste("Number of censored observation:  ", n.censor, "\n"))
  cat(paste("Number of predictors:            ", ncol(object$x), "\n"))
  cat(paste("Number of slices (uncensored):   ", object$nslice1, "\n"))
  cat(paste("Number of slices (censored):     ", object$nslice0, "\n\n"))
  cat("Eigenvalues:\n")
  cum.e <- cumsum(object$eval / sum(object$eval))
  evals <- rbind(object$eval, cum.e)
  rownames(evals) <- c("Eigenvalues", "Cum.Sum.R^2")
  colnames(evals) <- colnames(object$evec)
  print(evals, digits = digits)
  cat("\nAsym Chi-square test of SIR directions:\n")
  print(sir.test(object), digits = digits)
  cat("\nEigenvectors:\n")
  print(object$evec, digits = digits)
  invisible(object)
}

#' Affichage d'un objet "csir" (sortie de cen.sir)
#'
#' @inheritParams print.ds
#' @param object Objet de classe \code{"csir"}.
#' @export
print.csir <- function(object, digits = max(3, getOption("digits") - 3), ...) {
  cat("Censored Sliced Inverse Regression Model\n\n")
  cat(paste("Number of observation:           ", length(object$y), "\n"))
  n.censor <- length(object$y[object$delta == 0])
  cat(paste("Number of censored observation:  ", n.censor, "\n"))
  cat(paste("Number of predictors:            ", ncol(object$x), "\n"))
  cat(paste("Number of slices:                ", (object$nslice), "\n"))
  cat(paste("Kernel width:                    ", (object$h), "\n\n"))
  cum.e <- cumsum(object$eval / sum(object$eval))
  evals <- rbind(object$eval, cum.e)
  rownames(evals) <- c("Eigenvalues", "Cum.Sum.R^2")
  colnames(evals) <- colnames(object$evec)
  cat("Eigenvalues:\n")
  print(evals, digits = digits)
  cat("\nAsym Chi-square test of SIR directions:\n")
  print(sir.test(object), digits = digits)
  cat("\nEigenvectors:\n")
  print(object$evec, digits = digits)
  invisible(object)
}
