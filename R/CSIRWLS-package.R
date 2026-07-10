#' CSIRWLS : Censored Sliced Inverse Regression avec selection de variables
#' par Weighted Least Squares (CSIR-WLS)
#'
#' @description
#' Ce package regroupe l'ensemble des methodes developpees pour le papier
#' CSIR-WLS : une methode de reduction de dimension (Sliced Inverse
#' Regression) pour donnees de survie censurees, combinee a une selection
#' de variables par score de levier pondere (Weighted Least Squares, WLS).
#'
#' Le package est organise en 9 fichiers thematiques dans \code{R/} :
#' \itemize{
#'   \item \strong{01-linear-algebra-utils.R} : petites fonctions d'algebre
#'     lineaire utilisees partout (cosinus entre vecteurs, projections,
#'     covariance biaisee, decomposition generalisee en valeurs propres).
#'   \item \strong{02-core-censored-sir.R} : le coeur de la methode CSIR
#'     (noyau, survie conditionnelle ponderee, double slicing, cen.sir).
#'   \item \strong{03-plotting.R} : fonctions graphiques (SIR/CSIR, Kaplan-Meier
#'     par groupes de risque).
#'   \item \strong{04-variable-selection-wls.R} : la selection de variables
#'     WLS et le pipeline complet CSIR-WLS (cen.wls, cen.wls.sir2, choix
#'     automatique de la dimension d).
#'   \item \strong{05-screening-methods.R} : methodes de screening pour la
#'     tres grande dimension (CSIS, DCSIS, SIRS, SIS) utilisees comme
#'     methodes concurrentes dans le papier.
#'   \item \strong{06-data-simulation.R} : generateurs de donnees pour les
#'     etudes de simulation (X_rand, generate_X, temps de survie simules).
#'   \item \strong{07-evaluation-metrics.R} : criteres de comparaison des
#'     methodes (cosinus/PC, distance de Frobenius entre sous-espaces,
#'     TPR/FPR/FDR).
#'   \item \strong{08-simulation-runners.R} : boucles de simulation Monte-Carlo
#'     (paralleles), pour l'estimation EDR (SIR/CSIR/COX/CSIR-WLS) et pour
#'     le screening de variables (CWLS vs CSIS).
#'   \item \strong{09-reporting-tables.R} : mise en forme des resultats
#'     (tables LaTeX via xtable, graphiques de temps de calcul).
#' }
#'
#' @section Note de migration (important) :
#' Plusieurs de vos scripts d'origine definissaient des fonctions
#' \emph{portant le meme nom mais avec un contenu different} (ex :
#' \code{run_simulation_parallel}, \code{one_rep}, \code{latex_table},
#' \code{generate_X}). Dans le package, chaque variante a recu un nom
#' unique et explicite pour eviter tout ecrasement silencieux :
#' \itemize{
#'   \item \code{run_simulation_parallel()} (ancien, comparaison EDR
#'     SIR/CSIR/COX/CSIR-WLS, depuis \code{cwlsSIR parallel .R}) devient
#'     \code{\link{run_simulation_edr}()}.
#'   \item \code{run_simulation_parallel()} (ancien, comparaison de
#'     screening CWLS vs CSIS, depuis \code{CWLS-CSIS-bcorCISparalelle.R})
#'     devient \code{\link{run_simulation_screening}()}.
#'   \item \code{latex_table()} (PC/Cor/Frob) devient
#'     \code{\link{latex_table_edr}()} ; \code{latex_table()} (TPR/FPR/FDR)
#'     devient \code{\link{latex_table_screening}()}.
#'   \item \code{generate_X()} a ete fusionne en une seule version
#'     canonique (\code{\link{generate_X}}) qui supporte tous les types
#'     rencontres dans vos scripts : \code{"normal"}, \code{"autoR"},
#'     \code{"block"}, \code{"spik"}, \code{"skewnormal"}, \code{"skewt"},
#'     \code{"mixture"}, \code{"indep_nongauss"}.
#' }
#' Voir le README pour la table de correspondance complete ancien nom ->
#' nouveau nom, fichier par fichier.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom MASS mvrnorm
#' @importFrom stats coef cor cov cutree lm median na.omit predict prcomp
#' @importFrom stats quantile rbinom rchisq rexp rnorm runif sd uniroot bw.nrd
#' @importFrom stats chol pchisq
#' @importFrom survival Surv survfit survdiff coxph concordancefit
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
## usethis namespace: end
NULL
