## =============================================================================
## Etude de simulation : estimation EDR a index simple (d = 1)
## SIR vs CSIR vs Cox vs CSIR-WLS, sous plusieurs taux de censure
##
## Source d'origine : "cwlsSIR parallel .R" / "CSIRWLS-withCPU time .R"
## (fonctions equivalentes fusionnees en une seule, voir run_simulation_edr()
## dans le package -- la version sequentielle de "CSIRWLS-withCPU time .R"
## faisait exactement le meme calcul, sans parallelisme).
## =============================================================================

library(CSIRWLS)
library(MASS)
library(survival)
library(Rdimtools) # methode SIR de reference

## ---------------------------------------------------------------------------
## Parametres de la simulation
## ---------------------------------------------------------------------------
n <- 500
p <- 80
s <- 1:5 # indices des variables actives
q <- length(s)

Beta <- rep(0, p)
Beta[s] <- runif(length(s), 0.5, 1.5)

rho <- 0.5
X_type <- "autoR" # structure de covariance auto-regressive
model_type <- "PHM" # modele a hasards proportionnels

censor_rates <- c(0.15, 0.25, 0.40)
num_repeats <- 100 # reduire pour un test rapide (ex : 10)

n.slice <- 10
n.slice1 <- 5
n.slice0 <- 5

n_cores <- max(1, parallel::detectCores() - 1)

## ---------------------------------------------------------------------------
## Lancer la simulation
## ---------------------------------------------------------------------------
set.seed(1)
total_time <- system.time({
  sim_output <- run_simulation_edr(
    n = n, p = p, q = q, Beta = Beta,
    censor_rates = censor_rates, num_repeats = num_repeats,
    X_type = X_type, model_type = model_type,
    rho = rho, n_cores = n_cores,
    n.slice = n.slice, n.slice1 = n.slice1, n.slice0 = n.slice0
  )
})
cat("Temps total :", total_time["elapsed"], "secondes\n")

## ---------------------------------------------------------------------------
## Table LaTeX des resultats (PC = |cosinus|, Cor, Frob = distance de
## projection, Time_s = temps de calcul moyen)
## ---------------------------------------------------------------------------
latex_table_edr(sim_output$results, censor_rates)

## ---------------------------------------------------------------------------
## (Optionnel) robustesse sous distributions non-elliptiques de X
## Reponse type a un commentaire de reviewer sur la Condition C.1
## ---------------------------------------------------------------------------
X_types_nonelliptic <- c("autoR", "skewnormal", "skewt")

results_robustness_EDR <- list()
for (xt in X_types_nonelliptic) {
  cat("=== X_type:", xt, "===\n")
  sim_out_xt <- tryCatch(
    run_simulation_edr(
      n = n, p = p, q = q, Beta = Beta,
      censor_rates = censor_rates, num_repeats = num_repeats,
      X_type = xt, model_type = model_type, rho = rho,
      n_cores = n_cores,
      n.slice = n.slice, n.slice1 = n.slice1, n.slice0 = n.slice0
    ),
    error = function(e) {
      cat("  Echec pour X_type =", xt, ":", conditionMessage(e), "\n")
      NULL
    }
  )
  if (!is.null(sim_out_xt)) results_robustness_EDR[[xt]] <- sim_out_xt$results
}

if (length(results_robustness_EDR) > 0) {
  latex_table_robustness_EDR(results_robustness_EDR, censor_rates)
}
