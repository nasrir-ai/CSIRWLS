## =============================================================================
## Analyse de donnees reelles GSE7390 (transbig2006, cancer du sein) :
## screening de genes CSIS vs CWLS (CSIR-WLS), puis courbes de Kaplan-Meier
##
## Source d'origine : "CSIS-CWLSIR-realdata.R" (section GSE7390)
## Nettoyage effectue : suppression des chemins personnels
## ("~/Documents/phd /R codes /...", "/Users/mac/Documents/phd /...").
##
## Donnees : ce script attend le fichier GSE7390_transbig2006affy.RData
## (dataset "transbig2006" du package Bioconductor 'breastCancerTRANSBIG',
## GEO accession GSE7390). Un exemplaire est fourni dans
## inst/extdata/GSE7390_transbig2006affy.RData -- retrouvez-le avec
## system.file("extdata", "GSE7390_transbig2006affy.RData",
## package = "CSIRWLS"). Ce fichier (~33 Mo) N'EST PAS versionne dans le
## depot GitHub (voir .gitignore) pour eviter d'alourdir le repo ; gardez
## une copie locale, ou utilisez Git LFS si vous voulez le versionner.
##
## IMPORTANT : les indices de colonnes ci-dessous (colonnes 29+ pour
## l'expression genique) dependent de la structure exacte de l'objet
## `demo` charge avec les donnees. Verifiez avec str(merged_data) avant
## de faire confiance aux indices sur un nouveau jeu de donnees.
## =============================================================================

library(CSIRWLS)
library(survival)
library(ggplot2)

data_path <- system.file("extdata", "GSE7390_transbig2006affy.RData", package = "CSIRWLS")
if (data_path == "") {
  stop(
    "Fichier GSE7390_transbig2006affy.RData introuvable. Placez-le dans ",
    "inst/extdata/ (voir l'en-tete de ce script) ou changez data_path ci-dessous."
  )
}
load(data_path) # charge les objets `data` (expression genique) et `demo` (clinique)

## ---------------------------------------------------------------------------
## 1) Preparation des donnees
## ---------------------------------------------------------------------------
data_df <- as.data.frame(data)
dim(data_df) # attendu : 198 echantillons x 22283 genes

data_standardized <- as.data.frame(scale(data_df))
merged_data <- cbind(demo, data_standardized)
str(merged_data) # verifiez ici la position reelle des colonnes d'expression

time <- merged_data$t.dmfs # temps jusqu'a metastase a distance ou censure
event <- merged_data$e.dmfs # 1 = metastase, 0 = censure

X <- as.matrix(merged_data[, 29:ncol(merged_data)]) # donnees d'expression genique
colnames(X) <- colnames(merged_data)[29:ncol(merged_data)]

## ---------------------------------------------------------------------------
## 2) Screening de genes : CSIS vs CWLS
## ---------------------------------------------------------------------------
csis_selected_genes <- apply_csis(X, time, event)
top_5_genes_csis <- colnames(X)[csis_selected_genes[1:5]]
top_10_genes_csis <- colnames(X)[csis_selected_genes[1:10]]

cwls_result <- cen.wls.sir2(X, time, event,
  n.slice1 = 10, n.slice0 = 5, n.slice = 5,
  cn1 = 0.1, cn2 = 1, choose.dir = FALSE
)
cwls_selected_genes <- cwls_result$select
top_5_genes_cwls <- cwls_selected_genes[1:5]
top_10_genes_cwls <- cwls_selected_genes[1:10]

beta_cwlsir <- as.matrix(round(Re(cwls_result$beta.hat), 4))

## ---------------------------------------------------------------------------
## 3) Comparaison avec Cox sur les genes selectionnes par CSIS
## ---------------------------------------------------------------------------
cox_model_csis <- coxph(Surv(time, event) ~ ., data = merged_data[, csis_selected_genes])
cox_coef_CSIS <- as.matrix(round(cox_model_csis$coefficients, 4))
X_CSIS <- as.matrix(merged_data[, csis_selected_genes])
cox_projection <- -X_CSIS %*% cox_coef_CSIS

CSIR_projection <- X %*% beta_cwlsir

## ---------------------------------------------------------------------------
## 4) Groupes de risque a 2 niveaux (median) et courbes de Kaplan-Meier
## ---------------------------------------------------------------------------
cox_group <- factor(ifelse(cox_projection > median(cox_projection), "High-Cox", "Low-Cox"))
km_fit_cox <- survfit(Surv(time, event) ~ cox_group)

CSIR_group <- factor(ifelse(CSIR_projection > median(CSIR_projection), "High-CWLS", "Low-CWLS"))
km_fit_CSIR <- survfit(Surv(time, event) ~ CSIR_group)

par(mfrow = c(1, 2))
plot(km_fit_cox,
  col = c("red", "blue"), lty = 1:2, lwd = 2,
  xlab = "Time (days)", ylab = "Survival Probability",
  main = "Kaplan-Meier Curves by Cox Projection (CSIS selected genes)"
)
legend("bottomleft", legend = c("High-Cox", "Low-Cox"), col = c("red", "blue"), lty = 1:2)

plot(km_fit_CSIR,
  col = c("red", "blue"), lty = 1:2, lwd = 2,
  xlab = "Time (days)", ylab = "Survival Probability",
  main = "Kaplan-Meier Curves by CSIR-WLS Projection (CWLS selected genes)"
)
legend("bottomleft", legend = c("High-CWLS", "Low-CWLS"), col = c("red", "blue"), lty = 1:2)
par(mfrow = c(1, 1))

## ---------------------------------------------------------------------------
## 5) Groupes de risque a 3 niveaux (tertiles)
## ---------------------------------------------------------------------------
risk_groups_cox <- make_risk_groups(cox_projection)
risk_groups_CSIR <- make_risk_groups(CSIR_projection)

km_fit_cox3 <- survfit(Surv(time, event) ~ risk_groups_cox)
km_fit_CSIR3 <- survfit(Surv(time, event) ~ risk_groups_CSIR)

par(mfrow = c(1, 2))
plot(km_fit_cox3,
  col = c("blue", "green", "red"), lty = 1:3, lwd = 2,
  xlab = "Time (days)", ylab = "Survival Probability",
  main = "Kaplan-Meier Curves by Cox Projection (CSIS selected genes)"
)
plot(km_fit_CSIR3,
  col = c("blue", "green", "red"), lty = 1:3, lwd = 2,
  xlab = "Time (days)", ylab = "Survival Probability",
  main = "Kaplan-Meier Curves by CSIR-WLS Projection (CWLS selected genes)"
)
par(mfrow = c(1, 1))

## ---------------------------------------------------------------------------
## 6) Version ggsurvplot (necessite le package 'survminer')
## ---------------------------------------------------------------------------
if (requireNamespace("survminer", quietly = TRUE)) {
  library(survminer)

  plot_cox <- ggsurvplot(km_fit_cox,
    data = merged_data,
    palette = c("red", "blue"), linetype = 1:2,
    xlab = "Time (days)", ylab = "Survival Probability",
    title = "Kaplan-Meier Curves by Cox Projection (CSIS selected genes)",
    legend.title = "Risk Group", legend.labs = c("High-Cox", "Low-Cox"),
    risk.table = TRUE, pval = TRUE, ggtheme = theme_minimal()
  )

  plot_CSIR <- ggsurvplot(km_fit_CSIR,
    data = merged_data,
    palette = c("red", "blue"), linetype = 1:2,
    xlab = "Time (days)", ylab = "Survival Probability",
    title = "Kaplan-Meier Curves by CSIR-WLS Projection (CWLS selected genes)",
    legend.title = "Risk Group", legend.labs = c("High-CWLS", "Low-CWLS"),
    risk.table = TRUE, pval = TRUE, ggtheme = theme_minimal()
  )

  print(plot_cox)
  print(plot_CSIR)

  if (requireNamespace("patchwork", quietly = TRUE) || requireNamespace("survminer", quietly = TRUE)) {
    combined_plots <- arrange_ggsurvplots(list(plot_cox, plot_CSIR), print = TRUE, ncol = 2, nrow = 1)
    print(combined_plots)
  }
} else {
  message("Package 'survminer' non installe : etape ggsurvplot ignoree.")
}

## ---------------------------------------------------------------------------
## 7) Nuage de points projection vs temps, colore par statut
## ---------------------------------------------------------------------------
plot_data <- data.frame(
  time = rep(time, 2),
  projection = c(CSIR_projection, cox_projection),
  status = rep(event, 2),
  type = factor(rep(c("CSIR-WLS", "Cox"), each = length(time)))
)

ggplot(plot_data, aes(x = time, y = projection)) +
  geom_point(size = 1, aes(color = type, shape = factor(status))) +
  scale_shape_manual(values = c(3, 16), labels = c("Censored", "Event")) +
  scale_color_manual(values = c("CSIR-WLS" = "brown", "Cox" = "purple")) +
  labs(
    title = "Projections Over Time",
    x = "Time", y = "Projection", color = "Method", shape = "Status"
  ) +
  facet_wrap(~type, ncol = 2, scales = "free_y") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )
