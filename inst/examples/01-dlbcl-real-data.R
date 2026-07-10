## =============================================================================
## Analyse de donnees reelles DLBCL : Cox vs CSIR vs CSIR-WLS
## Reproduit Table 15, Figure 2 et Figure 3 du papier CSIR-WLS
##
## Source d'origine : "CSIR-WLS real data DLBCL.R"
## Nettoyage effectue : suppression du chemin personnel
## "~/Documents/phd /R codes /myfunctions.R" -- desormais inutile, toutes
## les fonctions viennent du package CSIRWLS installe.
## =============================================================================

library(CSIRWLS)
library(ipred) # fournit le jeu de donnees DLBCL
library(survival)
library(xtable)

## ---------------------------------------------------------------------------
## 1) Chargement et nettoyage des donnees
## ---------------------------------------------------------------------------
data("DLBCL", package = "ipred")
str(DLBCL) # verification : noms/types de colonnes
DLBCL_clean <- na.omit(DLBCL)

time <- DLBCL_clean$time # temps de survie
delta <- DLBCL_clean$cens # 1 = evenement (deces), 0 = censure

## Selection des predicteurs PAR NOM (jamais par position) :
## IPI + MGEc.1 ... MGEc.10 -> 11 predicteurs, exactement comme dans le papier
pred_names <- c("IPI", paste0("MGEc.", 1:10))
X <- as.matrix(DLBCL_clean[, pred_names])
storage.mode(X) <- "double"

n <- nrow(X)
p <- ncol(X)
cat("n =", n, " p =", p, "\n") # attendu : n = 38 (apres na.omit), p = 11

## ---------------------------------------------------------------------------
## 2) Modele de Cox
## ---------------------------------------------------------------------------
cox_df <- as.data.frame(X)
cox_model <- coxph(Surv(time, delta) ~ ., data = cox_df)

cox_coef <- as.matrix(coef(cox_model))
rownames(cox_coef) <- pred_names
colnames(cox_coef) <- "COX"

## ---------------------------------------------------------------------------
## 3) Modele CSIR (index simple, d = 1, comme utilise pour DLBCL dans le papier)
## ---------------------------------------------------------------------------
n.slice <- 10
n.slice1 <- 5 # tranches pour les observations non censurees
n.slice0 <- 5 # tranches pour les observations censurees
h <- 0.2 # largeur de bande du noyau
c <- 0.05 # plancher de survie
ndim <- 1 # DLBCL utilise une seule direction EDR (contrairement a PBC, d = 2)

ds1 <- double.slice(time, delta, X, n.slice1, n.slice0)
joint.edrs1 <- edr.n(ds1, 2)

csir_result <- cen.sir(time, delta, X,
  n.slice = n.slice1,
  joint.edrs = joint.edrs1, h = h, c = c
)

beta_csir <- as.matrix(csir_result$evec[, 1:ndim])
rownames(beta_csir) <- pred_names
colnames(beta_csir) <- "CSIR"

## ---------------------------------------------------------------------------
## 4) Modele CSIR-WLS
## ---------------------------------------------------------------------------
csir.wls <- cen.wls.sir2(X, as.vector(time), delta,
  n.slice1 = n.slice1, n.slice0 = n.slice0,
  n.slice = n.slice,
  cn1 = 0.1, cn2 = 3, choose.dir = FALSE,
  ndim = ndim, c = c
)

beta_csirwls <- as.matrix(Re(csir.wls$beta.hat))
rownames(beta_csirwls) <- pred_names
colnames(beta_csirwls) <- "CSIR-WLS"

cat("Nombre de variables retenues par CWLS :", csir.wls$n.sel, "sur", p, "\n")

## ---------------------------------------------------------------------------
## 5) Table 15 : coefficients cote a cote
## ---------------------------------------------------------------------------
results_table <- round(cbind(cox_coef, beta_csir, beta_csirwls), 4)
print(results_table)

latex_table <- xtable(results_table,
  caption = "Estimated directional coefficients for the DLBCL dataset",
  label = "tab:dlbcl_coef"
)
print(latex_table, include.rownames = TRUE)

## ---------------------------------------------------------------------------
## 6) Projections
## ---------------------------------------------------------------------------
cox_projection <- X %*% cox_coef
csir_projection <- X %*% beta_csir
csirwls_projection <- X %*% beta_csirwls

## ---------------------------------------------------------------------------
## 7) Figure 2 : direction projetee vs temps observe
##    "." (pch = 16) = temps d'echec, "+" (pch = 3) = temps de censure
## ---------------------------------------------------------------------------
par(mfrow = c(1, 3))
plot_projection(cox_projection, time, delta, "COX")
plot_projection(csir_projection, time, delta, "CSIR")
plot_projection(csirwls_projection, time, delta, "CSIR-WLS")
par(mfrow = c(1, 1))

## ---------------------------------------------------------------------------
## 8) Figure 3 : courbes de Kaplan-Meier pour TROIS groupes de risque
##    Low risk  (< 33%)  = rouge
##    Medium risk (33-66%) = vert
##    High risk (> 66%)  = bleu
## ---------------------------------------------------------------------------
cox_group <- factor(make_risk_groups(cox_projection))
csir_group <- factor(make_risk_groups(csir_projection))
csirwls_group <- factor(make_risk_groups(csirwls_projection))

km_fit_cox <- survfit(Surv(time, delta) ~ cox_group)
km_fit_csir <- survfit(Surv(time, delta) ~ csir_group)
km_fit_csirwls <- survfit(Surv(time, delta) ~ csirwls_group)

plot_km3 <- function(km_fit, title) {
  plot(km_fit,
    col = c("red", "green", "blue"), lty = 1:3, lwd = 2,
    xlab = "Time (days)", ylab = "Survival Probability",
    main = title
  )
  legend("bottomleft",
    legend = c("Low risk", "Medium risk", "High risk"),
    col = c("red", "green", "blue"), lty = 1:3, bty = "n"
  )
}

par(mfrow = c(1, 3))
plot_km3(km_fit_cox, "Kaplan-Meier Curves by COX Projection")
plot_km3(km_fit_csir, "Kaplan-Meier Curves by CSIR Projection")
plot_km3(km_fit_csirwls, "Kaplan-Meier Curves by CSIR-WLS Projection")
par(mfrow = c(1, 1))

## ---------------------------------------------------------------------------
## 9) (Optionnel) tests du log-rank pour quantifier la separation par methode
## ---------------------------------------------------------------------------
cat("\nLog-rank test - COX groups:\n")
print(survdiff(Surv(time, delta) ~ cox_group))
cat("\nLog-rank test - CSIR groups:\n")
print(survdiff(Surv(time, delta) ~ csir_group))
cat("\nLog-rank test - CSIR-WLS groups:\n")
print(survdiff(Surv(time, delta) ~ csirwls_group))
