#' csirwls: Censored Sliced Inverse Regression with Weighted Leverage Screening
#'
#' Implementation of the procedure of Asrir, Bouezmarni and Mkhadri,
#' "Enhancing censored sliced inverse regression through weighted leverage
#' scores", for right-censored survival data with many predictors. A censored
#' weighted leverage score (\code{\link{cwls}}) ranks the predictors, the
#' stopping rule \eqn{G_c} (\code{\link{bic_cut}}) chooses how many to keep, and
#' censored sliced inverse regression (\code{\link{csir}}) estimates the
#' effective dimension reduction directions from the retained ones
#' (\code{\link{csir_wls}}).
#'
#' The package also contains the competitors used in the paper (the marginal
#' screening statistics \code{\link{csis}}, \code{\link{cris}} and
#' \code{\link{cmv_sis}}; \code{\link{sir}}, \code{\link{cp_sir}},
#' \code{\link{cox_fit}} and \code{\link{cox_lasso_bic}}), the evaluation
#' criteria, the simulation designs and a cross-validation routine
#' (\code{\link{cv_survival}}). It depends on base R only. The scripts that
#' reproduce every table of the paper are in the \code{reproduction/} folder of
#' the source repository.
#'
#' @keywords internal
"_PACKAGE"
