## =============================================================================
## Etude de simulation : screening de variables CWLS vs CSIS (TPR/FPR/FDR)
## en tres grande dimension, y compris robustesse a des structures de
## correlation par blocs et a des distributions non-elliptiques de X.
##
## Source d'origine : "CWLS-CSIS-bcorCISparalelle.R"
## =============================================================================

library(CSIRWLS)
library(survival)
library(MASS)

## ---------------------------------------------------------------------------
## Parametres
## ---------------------------------------------------------------------------
censor_rates <- c(0.15, 0.25, 0.40)
num_repeats <- 100 # reduire pour un test rapide (ex : 10)
n <- 500
p <- 10000 # essayer aussi : 80, 300, 1500

set.seed(1)

s <- 1:5
q <- length(s)

Beta <- rep(0, p)
Beta[s] <- runif(length(s), 0.5, 1.5)

rho <- 0.5
X_type <- "autoR"
model_type <- "cox"

n_cores <- max(1, parallel::detectCores() - 1)

stopifnot(length(Beta) == p, max(s) <= p)

## ---------------------------------------------------------------------------
## 1) Simulation principale (structure autoR, elliptique, baseline)
## ---------------------------------------------------------------------------
cat("Demarrage simulation parallele...\n")
total_time <- system.time({
  sim_output <- run_simulation_screening(
    n = n, p = p, q = q, Beta = Beta, s = s,
    censor_rates = censor_rates, num_repeats = num_repeats,
    X_type = X_type, model_type = model_type, rho = rho,
    n_cores = n_cores
  )
})
cat("Temps total de simulation :", total_time["elapsed"], "secondes\n")

latex_table_screening(sim_output$results, censor_rates)
latex_time_table_screening(sim_output$results, censor_rates, n, p)
plot_time_by_censor(sim_output, n, p)

## ---------------------------------------------------------------------------
## 2) Robustesse : distributions non-elliptiques de X (reponse type a un
##    commentaire de reviewer sur la Condition C.1)
## ---------------------------------------------------------------------------
model_type_robust <- "PHM"

X_types_nonelliptic <- c(
  "autoR", # reference : gaussien auto-regressif, satisfait C.1
  "skewnormal", # asymetrique, queues legeres : violation moderee de C.1
  "skewt" # asymetrique ET queues lourdes : violation plus severe
  # "mixture" et "indep_nongauss" sont aussi supportees par generate_X()
  # mais retirees ici pour le temps de calcul (voir generate_X() dans le
  # package pour les reactiver).
)

results_robustness <- list()
for (xt in X_types_nonelliptic) {
  cat("=== X_type:", xt, "===\n")
  sim_out_xt <- tryCatch(
    run_simulation_screening(
      n = n, p = p, q = q, Beta = Beta, s = s,
      censor_rates = censor_rates, num_repeats = num_repeats,
      X_type = xt, model_type = model_type_robust, rho = rho,
      n_cores = n_cores
    ),
    error = function(e) {
      cat("  Echec pour X_type =", xt, ":", conditionMessage(e), "\n")
      NULL
    }
  )
  if (!is.null(sim_out_xt)) results_robustness[[xt]] <- sim_out_xt$results
}

if (length(results_robustness) > 0) {
  latex_table_robustness_screening(results_robustness, censor_rates)
}

## ---------------------------------------------------------------------------
## 3) Impact de la structure de correlation par blocs (X_type = "block")
##    Reponse type a un commentaire de reviewer sur l'impact de la
##    correlation entre variables actives/inactives.
## ---------------------------------------------------------------------------
n_small <- 80
p_block <- 10000
s_block <- 1:20
q_block <- length(s_block)

Beta_block <- rep(0, p_block)
Beta_block[s_block] <- runif(length(s_block), 0.5, 1.5)

alpha_list <- list(
  "weak correlation" = c(0.2, 0.2, 0.2),
  "strong active-active" = c(0.8, 0.2, 0.2),
  "moderate active-inactive" = c(0.6, 0.5, 0.3),
  "strong inactive-inactive" = c(0.5, 0.2, 0.8),
  "mixed strong correlation" = c(0.8, 0.6, 0.6)
)

## Decommenter pour lancer l'etude complete (peut etre long) :
# alpha_results <- run_alpha_study(
#   alpha_list = alpha_list, n = n_small, p = p_block, q = q_block,
#   Beta = Beta_block, censor_rates = censor_rates,
#   num_repeats = num_repeats, model_type = "PHM", n_cores = n_cores
# )
# print(alpha_results)
