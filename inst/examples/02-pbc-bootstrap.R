## =============================================================================
## Analyse de donnees reelles PBC (Primary Biliary Cirrhosis) avec bootstrap
## CSIR vs CSIR-WLS, avec ecarts-types et IC95% bootstrap
##
## Source d'origine : "CSIR real data app PBC with boostrap .R"
## Nettoyage effectue :
##  - suppression du chemin personnel
##    "~/Documents/phd /R codes /myfunctions.R"
##  - la fonction d'erreur du bootstrap CSIR-WLS referencait une variable
##    globale "num_vars" definie APRES l'appel a boot() dans le script
##    d'origine (donc inexistante si une erreur survenait sur la toute
##    premiere replication) ; elle utilise maintenant la variable locale
##    "nvar" deja calculee dans la fonction.
##  - la dimension EDR n'est plus fixee en dur a 2 : elle est estimee par
##    le critere BIC (Section 5.3.2 du papier) via select_d_bic(), qui
##    filtre correctement les valeurs propres proches de 0 (voir la note
##    ci-dessous). Le reste du script (tables, graphiques) se generalise
##    automatiquement a la dimension estimee, quelle qu'elle soit.
##
## Note sur estimate_K() vs select_d_bic() :
## une premiere version de ce script re-implementait ce critere BIC dans
## une fonction locale "estimate_K()" qui filtrait les valeurs propres
## avec eigvals[eigvals > 0] (sans tolerance numerique). Avec un petit
## nombre de tranches (n.slice), le rang theorique de la matrice
## sigma.eta est borne par (n.slice - 1) ; au-dela, les "valeurs propres"
## calculees ne sont que du bruit flottant (~1e-15), pas des zeros exacts,
## et passaient donc ce filtre -- ce qui gonflait artificiellement le
## nombre de candidats et rendait le d_hat quasi insensible aux donnees.
## select_d_bic() (dans R/04-variable-selection-wls.R) filtre avec un
## seuil de tolerance (> 1e-10), ce qui evite ce probleme. Utilisez
## toujours select_d_bic() plutot que de re-ecrire ce filtre a la main.
## =============================================================================

library(CSIRWLS)
library(survival)
library(boot)
library(xtable)
library(ggplot2)

## ---------------------------------------------------------------------------
## 1) Chargement et nettoyage des donnees
## ---------------------------------------------------------------------------
data("pbc", package = "survival")
str(pbc)
summary(pbc)

pbc_clean <- na.omit(pbc)
pbc_clean$sex <- ifelse(pbc_clean$sex == "m", 1, 0)

## ---------------------------------------------------------------------------
## 2) Variables : les 6 premieres comme dans CSIR 1999, puis toutes les
##    variables (17), avec les 6 premieres en tete pour comparabilite
## ---------------------------------------------------------------------------
base_vars <- c("age", "edema", "bili", "albumin", "platelet", "protime")
all_vars <- setdiff(colnames(pbc_clean)[4:20], base_vars)

## X17 : toutes les variables, 6 premieres en tete
X17 <- pbc_clean[, c(base_vars, all_vars)]

time <- pbc_clean$time
status <- ifelse(pbc_clean$status == 2, 1, 0) # 1 = deces, 0 = censure
X <- as.matrix(X17)

stopifnot(is.numeric(X))
dim(X)

## ---------------------------------------------------------------------------
## 3) Parametres CSIR / CSIR-WLS
## ---------------------------------------------------------------------------
n.slice1 <- 5
n.slice0 <- 5
n.slice <- 10
h <- 0.2
c <- 0.05
cn1 <- 0.1
cn2 <- 1

pbc_data <- data.frame(time = time, status = status, X)
num_vars <- ncol(X) # nombre de variables
variables <- colnames(X) # noms de variables

## ---------------------------------------------------------------------------
## 3 bis) Estimation de la dimension EDR d (Section 5.3.2 du papier)
##
## joint.edrs_prelim ci-dessous sert uniquement de base technique pour le
## lissage a noyau de cen.sir() (comme dans double.slice() partout
## ailleurs) -- elle n'a pas besoin de correspondre a la dimension
## structurelle finale, qui est justement ce qu'on estime ici.
## ---------------------------------------------------------------------------
ds_prelim <- double.slice(time, status, X, n.slice1, n.slice0)
joint.edrs_prelim <- edr.n(ds_prelim, 2)

csir_prelim <- cen.sir(time, status, X,
  n.slice = n.slice1,
  joint.edrs = joint.edrs_prelim, h = h, c = c
)

d_hat_bic <- select_d_bic(csir_prelim$eval, n = nrow(X), cn1 = cn1)
d_hat_chisq <- estimate_dim_chisq(csir_prelim, alpha = 0.01)
cat("Dimension estimee par le critere BIC (Section 5.3.2) :", d_hat_bic, "\n")
cat("Dimension estimee par le test sequentiel du chi2 (Li 1991) :", d_hat_chisq, "\n")

## Garde-fou contre un d_hat = 0 degenere (peut arriver si la courbe BIC
## est monotone) -- on retombe sur 1 direction, un modele EDR a d=0
## n'ayant pas de sens.
ndim <- max(1, d_hat_bic)
cat("ndim utilise pour CSIR / CSIR-WLS :", ndim, "\n")

## ---------------------------------------------------------------------------
## 4) Fonctions bootstrap
## ---------------------------------------------------------------------------
csir_boot <- function(data, indices) {
  d <- data[indices, ]
  Time <- d$time
  Status <- d$status
  Xmat <- as.matrix(d[, -(1:2)])
  ds1 <- double.slice(Time, Status, Xmat, n.slice1, n.slice0)
  joint.edrs1 <- edr.n(ds1, ndim)
  csir_result <- cen.sir(Time, Status, Xmat, n.slice = n.slice1, joint.edrs = joint.edrs1, h = h, c = c)
  as.vector(csir_result$evec[, 1:ndim])
}

csir_wls_boot <- function(data, indices) {
  nvar <- ncol(data) - 2 # data = time, status, puis les variables
  out <- tryCatch(
    {
      d <- data[indices, ]
      Time <- d$time
      Status <- d$status
      Xmat <- as.matrix(d[, -(1:2)])
      csir_wls_result <- cen.wls.sir2(Xmat, as.vector(Time), Status,
        n.slice1 = n.slice1,
        n.slice0 = n.slice0, n.slice = n.slice, cn1 = cn1,
        cn2 = cn2, choose.dir = FALSE, ndim = ndim
      )
      as.vector(Re(csir_wls_result$beta.hat[, 1:ndim]))
    },
    error = function(e) rep(NA_real_, nvar * ndim)
  )
  out
}

## ---------------------------------------------------------------------------
## 5) Lancer les bootstraps
## ---------------------------------------------------------------------------
set.seed(123)
boot_csir <- boot(data = pbc_data, statistic = csir_boot, R = 500)
boot_csir_wls <- boot(data = pbc_data, statistic = csir_wls_boot, R = 500)

## ---------------------------------------------------------------------------
## 6) Moyennes et ecarts-types bootstrap
## ---------------------------------------------------------------------------
csir_mean <- colMeans(boot_csir$t)
csir_se <- apply(boot_csir$t, 2, sd)
csirwls_mean <- colMeans(boot_csir_wls$t)
csirwls_se <- apply(boot_csir_wls$t, 2, sd)

format_coefs_se <- function(mean_vec, se_vec) {
  paste0(round(mean_vec, 2), " (", round(se_vec, 2), ")")
}

csir_formatted <- format_coefs_se(csir_mean, csir_se)
csirwls_formatted <- format_coefs_se(csirwls_mean, csirwls_se)

## Reshape en matrices (num_vars x ndim) -- fonctionne pour tout ndim
csir_mat <- matrix(csir_formatted, nrow = num_vars, ncol = ndim)
csirwls_mat <- matrix(csirwls_formatted, nrow = num_vars, ncol = ndim)

## Construire la table finale dynamiquement (une paire de colonnes par direction)
dir_cols <- lapply(seq_len(ndim), function(k) {
  setNames(
    data.frame(csir_mat[, k], csirwls_mat[, k]),
    c(paste0("CSIR_Dir", k), paste0("CSIRWLS_Dir", k))
  )
})
final_table <- data.frame(
  Variable = variables,
  do.call(cbind, dir_cols),
  stringsAsFactors = FALSE
)

print(xtable(final_table, caption = "Comparison of CSIR and CSIR-WLS coefficients with standard errors"))

## ---------------------------------------------------------------------------
## 7) Projections et graphique facette (ggplot2) -- generalise a tout ndim
## ---------------------------------------------------------------------------
csir_coefs <- lapply(seq_len(ndim), function(k) csir_mean[((k - 1) * num_vars + 1):(k * num_vars)])
csirwls_coefs <- lapply(seq_len(ndim), function(k) csirwls_mean[((k - 1) * num_vars + 1):(k * num_vars)])

plot_data <- do.call(rbind, lapply(seq_len(ndim), function(k) {
  rbind(
    data.frame(
      time = time, projection = as.vector(X %*% csir_coefs[[k]]),
      type = paste0("CSIR Projection ", k), status = status
    ),
    data.frame(
      time = time, projection = as.vector(X %*% csirwls_coefs[[k]]),
      type = paste0("CSIR-WLS Projection ", k), status = status
    )
  )
}))
plot_data$type <- factor(plot_data$type)

ggplot(plot_data, aes(x = time, y = projection)) +
  geom_point(aes(color = type, shape = factor(status)), size = 1) +
  scale_shape_manual(values = c(3, 16), labels = c("Censored", "Event")) +
  labs(
    title = "Projections Over Time",
    x = "Time",
    y = "Projection",
    shape = "Status",
    color = "Projection Type"
  ) +
  facet_wrap(~type, ncol = 2, scales = "free_y") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

## ---------------------------------------------------------------------------
## 8) Intervalles de confiance bootstrap (95%)
## ---------------------------------------------------------------------------
alpha <- 0.05

get_boot_ci <- function(boot_obj) {
  apply(boot_obj$t, 2, function(x) {
    quantile(x, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  })
}

csir_ci <- get_boot_ci(boot_csir)
csirwls_ci <- get_boot_ci(boot_csir_wls)

csir_ci_low <- csir_ci[1, ]
csir_ci_high <- csir_ci[2, ]
csirwls_ci_low <- csirwls_ci[1, ]
csirwls_ci_high <- csirwls_ci[2, ]

format_coefs_ci <- function(mean_vec, se_vec, low_vec, high_vec) {
  paste0(
    round(mean_vec, 2), " (", round(se_vec, 2), ") [",
    round(low_vec, 2), ", ", round(high_vec, 2), "]"
  )
}

csir_formatted_ci <- format_coefs_ci(csir_mean, csir_se, csir_ci_low, csir_ci_high)
csirwls_formatted_ci <- format_coefs_ci(csirwls_mean, csirwls_se, csirwls_ci_low, csirwls_ci_high)

csir_mat_ci <- matrix(csir_formatted_ci, nrow = num_vars, ncol = ndim)
csirwls_mat_ci <- matrix(csirwls_formatted_ci, nrow = num_vars, ncol = ndim)

## Construire la table CI dynamiquement (meme principe qu'a l'etape 6)
dir_cols_ci <- lapply(seq_len(ndim), function(k) {
  setNames(
    data.frame(csir_mat_ci[, k], csirwls_mat_ci[, k]),
    c(paste0("CSIR_Dir", k), paste0("CSIRWLS_Dir", k))
  )
})
final_table_ci <- data.frame(
  Variable = variables,
  do.call(cbind, dir_cols_ci),
  stringsAsFactors = FALSE
)

print(xtable(final_table_ci,
  caption = "Comparison of CSIR and CSIR-WLS coefficients with standard errors and 95% confidence intervals"
), include.rownames = FALSE)
