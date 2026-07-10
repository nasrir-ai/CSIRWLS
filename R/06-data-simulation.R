#' Generateurs de donnees pour les etudes de simulation
#'
#' @description
#' Fonctions utilisees pour generer les jeux de donnees synthetiques des
#' etudes de simulation du papier : matrices de design \code{X} sous
#' differentes structures de covariance / distributions,
#' temps de survie simules (index simple ou multiple), et calibration du
#' taux de censure.
#'
#' @section Fusion de \code{generate_X} :
#' Vos scripts contenaient plusieurs versions de \code{generate_X()},
#' chacune supportant un sous-ensemble de types differents. La version
#' ci-dessous est la fusion de toutes les variantes rencontrees et
#' supporte : \code{"normal"}, \code{"autoR"}, \code{"block"},
#' \code{"spik"}, \code{"skewnormal"}, \code{"skewt"}, \code{"mixture"},
#' \code{"indep_nongauss"}.
#'
#' @section Fusion de \code{generate_survival_times} :
#' De meme, \code{generate_survival_times_multi()} (version "index
#' multiple", \code{Beta} = matrice p x d) et l'ancienne
#' \code{generate_survival_times()} (version "index simple", \code{Beta}
#' = vecteur) calculaient la meme chose : \code{eta = rowSums(X \%*\% Beta)}
#' fonctionne pour un \code{Beta} vecteur ou matrice. Une seule fonction
#' \code{\link{generate_survival_times}} est gardee ; l'ancien nom
#' \code{generate_survival_times_multi} reste disponible comme alias.
#'
#' @name data-simulation
NULL

#' Generer une matrice de design X avec une structure de covariance donnee
#'
#' Genere \code{X} (n x p) selon une structure de covariance par blocs
#' (\code{"block"} : correlation \code{alpha[1]} intra-actives,
#' \code{alpha[2]} actives/inactives, \code{alpha[3]} intra-inactives),
#' auto-regressive (\code{"autoR"} : \code{Cov(x_i, x_j) = rho^|i-j|}),
#' homogene (\code{"homog"} : correlation constante \code{rho} hors
#' diagonale), ou "homog2" (structure utilisee dans Wang & Yin, 2019,
#' "Sparse sliced inverse regression via Lasso").
#'
#' @param n Taille d'echantillon.
#' @param p Nombre de covariables.
#' @param q Nombre de coefficients non nuls (variables actives). Requis
#'   pour \code{type = "block"}.
#' @param type Structure de \code{Sigma} : \code{"block"} ou \code{"autoR"}
#'   (ou \code{"homog"} / \code{"homog2"}).
#' @param alpha Vecteur de 3 correlations (utilise avec \code{type = "block"}).
#' @param rho Correlation (mode auto-regressif ou homogene, utilise avec
#'   \code{type = "autoR"} / \code{"homog"} / \code{"homog2"}).
#' @param S Indices des variables actives (necessaire si les variables
#'   actives ne sont pas les q premieres ; utilise seulement pour
#'   \code{type = "homog2"}).
#' @return Matrice n x p tiree d'une loi normale multivariee N(0, Sigma).
#' @export
X_rand <- function(n, p, q, type = "block", alpha = c(0.5, 0.7, 0.9), rho = 0, S = S) {
  if (type == "block") {
    C11 <- matrix(alpha[1], nrow = q, ncol = q)
    diag(C11) <- 1

    C22 <- matrix(alpha[3], nrow = p - q, ncol = p - q)
    diag(C22) <- 1

    C12 <- matrix(alpha[2], nrow = q, ncol = p - q)

    sigma <- rbind(cbind(C11, C12), cbind(t(C12), C22))
  }

  ### Auto-regressive structure
  if (type == "autoR") {
    sigma <- rho^abs(outer(1:p, 1:p, "-"))
  }

  ### homogeneous structure
  if (type == "homog") {
    sigma <- diag(1, p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        if (i != j) sigma[i, j] <- rho
      }
    }
  }

  ### structure from "Sparse sliced inverse regression via lasso" (2019)
  if (type == "homog2") {
    sigma <- matrix(0, nrow = p, ncol = p)
    sigma[S, S] <- rho
    sigma[-S, -S] <- rho
    sigma[S, -S] <- 0.1
    sigma[-S, S] <- 0.1
    diag(sigma) <- 1
  }

  ### symetrisation
  if (isSymmetric(sigma) == FALSE) {
    sigma <- 0.5 * (sigma + t(sigma))
  }

  x <- mvtnorm::rmvnorm(n, rep(0, p), sigma)
  return(x)
}

#' Generer une matrice de design X (version canonique unifiee)
#'
#' Version fusionnee de \code{generate_X()} supportant tous les types de
#' distribution rencontres dans vos scripts de simulation.
#'
#' @param n Taille d'echantillon.
#' @param p Nombre de covariables.
#' @param q Nombre de variables actives (requis pour \code{type = "block"}).
#' @param type Un des \code{"normal"}, \code{"autoR"}, \code{"block"},
#'   \code{"spik"}, \code{"skewnormal"}, \code{"skewt"}, \code{"mixture"},
#'   \code{"indep_nongauss"}.
#' @param rho Correlation auto-regressive (types \code{"autoR"},
#'   \code{"skewnormal"}, \code{"skewt"}, \code{"mixture"}).
#' @param alpha Vecteur de 3 correlations (type \code{"block"}, voir
#'   \code{\link{X_rand}}).
#' @param spik Matrice de chargements (type \code{"spik"}).
#' @param skew_param Parametre d'asymetrie (types \code{"skewnormal"},
#'   \code{"skewt"}, defaut 5). Necessite le package \pkg{sn}.
#' @param df Degres de liberte (type \code{"skewt"}, defaut 5).
#' @return Matrice n x p.
#' @export
generate_X <- function(n, p, q = NULL, type, rho = NULL, alpha = NULL, spik = NULL,
                        skew_param = NULL, df = NULL) {
  if (type == "spik") {
    U <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = diag(1, p))
    X <- U %*% spik
  } else if (type == "block") {
    X <- X_rand(n = n, p = p, q = q, type = "block", alpha = alpha)
  } else if (type == "autoR") {
    X <- X_rand(n = n, p = p, q = q, type = "autoR", rho = rho)
  } else if (type == "normal") {
    X <- matrix(rnorm(n * p), n, p)
  } else if (type == "skewnormal") {
    # Skew-Normal Multivariee (Azzalini & Dalla Valle, 1996) : distribution
    # asymetrique (skewness > 0), mais a queues legeres. Viole la symetrie
    # de la Condition C.1 tout en conservant une structure de correlation
    # autoR(rho). Cas de violation moderee.
    if (!requireNamespace("sn", quietly = TRUE)) {
      stop("Le package 'sn' est necessaire pour generate_X(type = 'skewnormal').")
    }
    if (is.null(rho)) rho <- 0.5
    if (is.null(skew_param)) skew_param <- 5

    Sigma_corr <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
    xi <- rep(0, p)
    Omega <- Sigma_corr
    alpha_vec <- rep(skew_param, p)

    X <- sn::rmsn(n = n, xi = xi, Omega = Omega, alpha = alpha_vec)
    X <- scale(X, center = TRUE, scale = FALSE)
  } else if (type == "skewt") {
    # Skew-t Multivariee (Azzalini & Capitanio, 2003) : combine asymetrie
    # ET queues lourdes (df=5). Violation plus severe de la Condition C.1
    # que le skew-normal.
    if (!requireNamespace("sn", quietly = TRUE)) {
      stop("Le package 'sn' est necessaire pour generate_X(type = 'skewt').")
    }
    if (is.null(rho)) rho <- 0.5
    if (is.null(skew_param)) skew_param <- 5
    if (is.null(df)) df <- 5

    Sigma_corr <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
    xi <- rep(0, p)
    Omega <- Sigma_corr
    alpha_vec <- rep(skew_param, p)

    X <- sn::rmst(n = n, xi = xi, Omega = Omega, alpha = alpha_vec, nu = df)
    X <- scale(X, center = TRUE, scale = FALSE)
  } else if (type == "mixture") {
    # Melange bimodal de deux gaussiennes (50/50) avec shift=2 : viole
    # l'ellipticite de facon structurelle. NOTE : boucle for() sur n
    # observations -> lent pour n grand (voir vectorisation possible via
    # mvtnorm/MASS si besoin).
    if (is.null(rho)) rho <- 0.5
    Sigma_corr <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
    shift <- 2
    component <- rbinom(n, 1, 0.5)
    X <- matrix(0, n, p)
    for (i in 1:n) {
      mu_i <- if (component[i] == 1) rep(shift, p) else rep(-shift, p)
      X[i, ] <- MASS::mvrnorm(1, mu = mu_i, Sigma = Sigma_corr)
    }
    X <- scale(X, center = TRUE, scale = FALSE)
  } else if (type == "indep_nongauss") {
    # Composantes marginales independantes chi2(3) centrees : violation
    # radicale de la Condition C.1 (meme la forme des contours de
    # densite n'est plus elliptique).
    X <- matrix(0, n, p)
    for (j in 1:p) {
      X[, j] <- rchisq(n, df = 3) - 3
    }
  } else {
    stop("Unsupported X_type")
  }

  X
}

#' Generer des temps de survie simules (index simple ou multiple)
#'
#' Genere le predicteur lineaire (a index simple si \code{Beta} est un
#' vecteur, a index multiple si \code{Beta} est une matrice p x d, via
#' \code{eta = rowSums(X \%*\% Beta)}), puis simule un temps de survie
#' "vrai" selon le modele choisi.
#'
#' @param X Matrice n x p de predicteurs.
#' @param Beta Vecteur (p) ou matrice (p x d) de directions EDR vraies.
#' @param model_type Un des \code{"PHM"}, \code{"PHM2"}, \code{"cox"},
#'   \code{"sinus"}.
#' @return Vecteur de temps de survie simules (longueur n).
#' @export
generate_survival_times <- function(X, Beta, model_type) {
  n <- nrow(X)
  eta <- rowSums(as.matrix(X %*% Beta))

  if (model_type == "PHM") {
    Y0 <- exp(eta + rnorm(n))
  } else if (model_type == "PHM2") {
    u <- runif(n)
    err <- log(-log(1 - u))
    Y0 <- exp(-2.5 + eta + 0.25 * err)
  } else if (model_type == "cox") {
    if (any(is.na(eta))) stop("Linear predictor contains NA values.")
    if (any(is.infinite(eta))) stop("Linear predictor contains Inf values.")
    eta <- pmin(eta, 10) # clip to avoid Inf in rate
    Y0 <- rexp(n, rate = exp(eta))
  } else if (model_type == "sinus") {
    u <- runif(n)
    err <- log(-log(1 - u))
    Y0 <- exp(-2.5 + sin(0.1 * pi * eta) + 0.1 * (eta + 2)^2 + 0.25 * err)
  } else {
    stop("Unsupported model_type")
  }

  as.vector(Y0)
}

#' @rdname generate_survival_times
#' @export
generate_survival_times_multi <- function(X, Beta, model_type) {
  generate_survival_times(X, Beta, model_type)
}

#' Taux de censure resultant d'une constante c donnee
#'
#' Etant donne un temps de survie "vrai" \code{Y0}, un terme d'erreur
#' \code{errc} et une constante \code{c}, simule un temps de censure
#' \code{C = exp(errc + c)}, calcule le temps observe
#' \code{y = min(Y0, C)} et l'indicateur \code{delta = 1{Y0 <= C}}, puis
#' renvoie le taux de censure resultant.
#'
#' @param c Constante ajoutee dans \code{exp(errc + c)}.
#' @param Y0 Vrais temps de survie (non censures).
#' @param errc Terme d'erreur (meme longueur que \code{Y0}).
#' @return Proportion d'observations censurees (entre 0 et 1).
#' @export
compute_censoring_rate <- function(c, Y0, errc) {
  C <- exp(errc + c)
  y <- pmin(Y0, C)
  delta <- as.numeric(Y0 <= C)
  censoring_rate <- 1 - mean(delta)
  return(censoring_rate)
}

#' Calibrer c pour obtenir un taux de censure cible
#'
#' Recherche par resolution de racine (\code{uniroot}) la valeur de
#' \code{c} telle que \code{\link{compute_censoring_rate}(c, Y0, errc)}
#' soit egale au taux de censure cible \code{target_rate}. Sert a
#' calibrer les simulations de donnees censurees.
#'
#' @param target_rate Taux de censure cible (entre 0 et 1).
#' @param Y0,errc Comme dans \code{\link{compute_censoring_rate}}.
#' @return Valeur scalaire de c (racine trouvee par \code{uniroot}).
#' @export
find_c_for_censoring_rate <- function(target_rate, Y0, errc) {
  objective_function <- function(c) {
    current_rate <- compute_censoring_rate(c, Y0, errc)
    return(current_rate - target_rate)
  }

  result <- uniroot(objective_function, interval = c(-10, 10), extendInt = "yes")
  return(result$root)
}
