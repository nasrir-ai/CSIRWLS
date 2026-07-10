## =============================================================================
## Etude de simulation : estimation EDR a index multiple (d > 1)
## SIR vs CSIR vs Cox vs CSIR-WLS, dimension fixee puis dimension automatique
##
## Source d'origine : "CSIRWLS-Multiple indice.R" / "CSIRWLS multiple index1.R"
## =============================================================================

library(CSIRWLS)
library(MASS)
library(survival)
library(Rdimtools)

## ---------------------------------------------------------------------------
## Parametres
## ---------------------------------------------------------------------------
n <- 500
p <- 80 # essayer aussi : 300, 1500
d_true <- 2 # nombre vrai de directions EDR

set.seed(1)

## Beta vrai : p x d_true
## direction 1 : variables 1:5, direction 2 : variables 6:10
Beta <- matrix(0, p, d_true)
Beta[1:5, 1] <- 1
Beta[6:10, 2] <- 1
Beta <- qr.Q(qr(Beta)) # orthonormalisation

s <- 1:10 # union du support des deux directions
q <- length(s)

rho <- 0.5
X_type <- "autoR"
model_type <- "PHM"

censor_rates <- c(0.15, 0.25, 0.40)
num_repeats <- 20 # reduire pour un test rapide

n.slice <- 10
n.slice1 <- 5
n.slice0 <- 5

n_cores <- max(1, parallel::detectCores() - 1)

stopifnot(nrow(Beta) == p, ncol(Beta) == d_true)
cat("Verification d'orthogonalite de Beta (doit etre l'identite) :\n")
print(round(t(Beta) %*% Beta, 6))

## ---------------------------------------------------------------------------
## 1) Dimension fixee (d_true connue)
## ---------------------------------------------------------------------------
cat("Demarrage simulation multi-index (d =", d_true, ")...\n")
total_time <- system.time({
  sim_output <- run_simulation_multi(
    n = n, p = p, d_true = d_true, Beta = Beta,
    censor_rates = censor_rates, num_repeats = num_repeats,
    X_type = X_type, model_type = model_type, rho = rho,
    n_cores = n_cores,
    n.slice = n.slice, n.slice1 = n.slice1, n.slice0 = n.slice0
  )
})
cat("Temps total :", total_time["elapsed"], "secondes\n")

latex_table_multi(sim_output$results, censor_rates, d_true = d_true)

## ---------------------------------------------------------------------------
## 2) Dimension estimee automatiquement (BIC)
## ---------------------------------------------------------------------------
cat("Demarrage simulation multi-index avec d automatique...\n")
total_time_auto <- system.time({
  sim_output_auto <- run_simulation_multi_auto_d(
    n = n, p = p, d_true = d_true, Beta = Beta,
    censor_rates = censor_rates, num_repeats = num_repeats,
    X_type = X_type, model_type = model_type, rho = rho,
    n_cores = n_cores,
    n.slice = n.slice, n.slice1 = n.slice1, n.slice0 = n.slice0
  )
})
cat("Temps total :", total_time_auto["elapsed"], "secondes\n")

latex_table_multi_auto_d(sim_output_auto$results, censor_rates, d_true)

## ---------------------------------------------------------------------------
## 3) Boucle sur plusieurs structures de X et modeles de survie
##
## NOTE : run_simulation_multi()/one_rep_multi() n'acheminent pas de
## parametre "q" jusqu'a generate_X() (comme dans le script d'origine),
## donc X_type = "block" echouera ici (tryCatch l'attrape et l'ignore).
## Dans vos resultats publies, seul X_type = "autoR" a ete utilise pour
## les simulations a index multiple ; "block" fonctionne pour les
## simulations a index simple (voir 06-screening-cwls-vs-csis.R) ou q
## est bien transmis.
## ---------------------------------------------------------------------------
X_types_to_test <- c("autoR", "block")
model_types_to_test <- c("PHM", "PHM2", "cox")

results_all <- list()

for (xt in X_types_to_test) {
  for (mt in model_types_to_test) {
    key <- paste0(xt, "_", mt)
    cat("\n==================================\n")
    cat("X_type =", xt, "| model_type =", mt, "\n")
    cat("==================================\n")

    # alpha necessaire uniquement pour "block"
    alpha_val <- if (xt == "block") c(0.5, 0.2, 0.4) else NULL
    rho_val <- if (xt == "autoR") 0.5 else NULL

    sim_out <- tryCatch(
      run_simulation_multi(
        n = n, p = p, d_true = d_true, Beta = Beta,
        censor_rates = censor_rates, num_repeats = num_repeats,
        X_type = xt, model_type = mt,
        alpha = alpha_val, rho = rho_val,
        n_cores = n_cores,
        n.slice = n.slice, n.slice1 = n.slice1, n.slice0 = n.slice0
      ),
      error = function(e) {
        cat("  Echec:", conditionMessage(e), "\n")
        NULL
      }
    )

    if (!is.null(sim_out)) {
      results_all[[key]] <- sim_out$results
      latex_table_multi(sim_out$results, censor_rates, d_true = d_true)
    }
  }
}
