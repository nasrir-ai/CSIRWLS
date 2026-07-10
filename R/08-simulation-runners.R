#' Boucles de simulation Monte-Carlo (paralleles)
#'
#' @description
#' Ce fichier regroupe les boucles de simulation utilisees pour produire
#' les tableaux du papier. Il y a deux familles de comparaisons, chacune
#' declinee en version "index simple" (d=1) et parfois "index multiple"
#' (d>1) :
#' \enumerate{
#'   \item \strong{Estimation EDR} : SIR vs CSIR vs Cox vs CSIR-WLS,
#'     comparees via \code{\link{costeta}}/\code{\link{Dproj}} (d=1,
#'     \code{\link{run_simulation_edr}}) ou
#'     \code{\link{trace_corr}}/\code{\link{subspace_dist}} (d>1,
#'     \code{\link{run_simulation_multi}} /
#'     \code{\link{run_simulation_multi_auto_d}}).
#'   \item \strong{Screening de variables} : CWLS vs CSIS, comparees via
#'     TPR/FPR/FDR (\code{\link{run_simulation_screening}}).
#' }
#'
#' @section Renommage (voir aussi la note de migration dans
#'   \code{CSIRWLS-package.R}) :
#' Vos scripts definissaient plusieurs fonctions \code{run_simulation_parallel}
#' / \code{one_rep} differentes (memes noms, contenus differents selon le
#' script). Chaque variante a ete renommee explicitement :
#' \code{run_simulation_edr}/\code{one_rep_edr} (depuis
#' \code{cwlsSIR parallel .R}), \code{run_simulation_screening}/
#' \code{one_rep_screening} (depuis \code{CWLS-CSIS-bcorCISparalelle.R}).
#' La version sequentielle (non parallele) de la comparaison EDR, depuis
#' \code{CSIRWLS-withCPU time .R}, faisait exactement le meme calcul que
#' \code{one_rep_edr} : elle n'est pas dupliquee ici, il suffit d'appeler
#' \code{run_simulation_edr(..., n_cores = 1)}.
#'
#' @section Parallelisme :
#' Les workers utilisent desormais \code{.packages = "CSIRWLS"} dans
#' \code{foreach()} au lieu de \code{source("myfunctions.R")} +
#' \code{clusterExport()} de chaque fonction : une fois le package
#' installe, toutes les fonctions necessaires sont automatiquement
#' disponibles sur chaque worker.
#'
#' @name simulation-runners
NULL

# ============================================================================
# 1) ESTIMATION EDR, INDEX SIMPLE (d = 1) : SIR vs CSIR vs Cox vs CSIR-WLS
#    (source : cwlsSIR parallel .R)
# ============================================================================

#' Une repetition de la comparaison EDR a index simple
#'
#' Genere un jeu de donnees, puis estime la direction EDR (d=1) par SIR
#' (\code{Rdimtools::do.sir}), CSIR (\code{\link{cen.sir}}), Cox
#' (\code{survival::coxph}) et CSIR-WLS (\code{\link{cen.wls.sir2}}), et
#' calcule les criteres de comparaison (\code{\link{costeta}},
#' \code{\link{Dproj}}, \code{\link{safe_cor}}) ainsi que le temps de
#' calcul de chaque methode.
#'
#' @param rep_id Identifiant de repetition (utilise pour la graine
#'   aleatoire par l'appelant, pas dans cette fonction).
#' @param n,p,q Taille d'echantillon, dimension, nombre de variables actives.
#' @param Beta Vecteur p (direction EDR vraie).
#' @param censor_rate Taux de censure cible.
#' @param X_type,model_type,rho,alpha,spik Voir \code{\link{generate_X}} /
#'   \code{\link{generate_survival_times}}.
#' @param n.slice,n.slice1,n.slice0 Parametres de tranchage pour CSIR / CSIR-WLS.
#' @return Une liste : \code{abscos}, \code{Dist}, \code{Cor} (un scalaire
#'   par methode, ordre SIR/CSIR/COX/CSIR-WLS), \code{Time} (temps en
#'   secondes par methode).
#' @export
one_rep_edr <- function(rep_id, n, p, q, Beta, censor_rate,
                         X_type, model_type, rho, alpha, spik,
                         n.slice, n.slice1, n.slice0) {
  if (!requireNamespace("Rdimtools", quietly = TRUE)) {
    stop("Le package 'Rdimtools' est necessaire pour one_rep_edr() (methode SIR de reference).")
  }

  X <- generate_X(n, p, q, type = X_type, rho = rho, alpha = alpha, spik = spik)
  survival_times <- generate_survival_times(X, Beta, model_type)
  errc <- rnorm(n)
  c_val <- find_c_for_censoring_rate(censor_rate, survival_times, errc)
  censoring_times <- exp(c_val + errc)
  time <- pmin(survival_times, censoring_times)
  event <- as.numeric(survival_times <= censoring_times)

  beta_init <- prcomp(X, rank. = 1)$rotation[, 1]
  h <- bw.nrd(as.vector(X %*% beta_init))

  # --- SIR ---
  t1 <- system.time({
    sir_select <- Rdimtools::do.sir(as.matrix(X), time, ndim = 1)
    beta_sir <- as.vector(Re(sir_select$projection))
  })["elapsed"]

  # --- CSIR ---
  t2 <- system.time({
    if (p > n) {
      beta_censir <- rep(0, p)
    } else {
      ds <- double.slice(time, event, X, n.slice1, n.slice0)
      joint.edrs <- edr.n(ds, 2)
      censir_select <- cen.sir(time, event, X, n.slice, joint.edrs, h, c = 0.05)
      beta_censir <- censir_select$evec[, 1]
    }
  })["elapsed"]

  # --- COX ---
  t3 <- system.time({
    if (p < n) {
      cox_model <- survival::coxph(survival::Surv(time, event) ~ X)
      beta_cox <- cox_model$coef
    } else {
      beta_cox <- rep(0, p)
    }
  })["elapsed"]

  # --- CSIR-WLS ---
  t4 <- system.time({
    cwlsir_select <- cen.wls.sir2(
      X, time, event,
      n.slice1, n.slice0, n.slice,
      cn1 = 0.1, cn2 = 1,
      choose.dir = FALSE, c = 0.05
    )
    beta_cwlsir <- Re(cwlsir_select$beta.hat)
  })["elapsed"]

  abscos_rep <- c(
    abs(costeta(Beta, beta_sir)),
    abs(costeta(Beta, beta_censir)),
    abs(costeta(Beta, beta_cox)),
    abs(costeta(Beta, beta_cwlsir))
  )

  dist_rep <- c(
    Dproj(Beta, beta_sir),
    Dproj(Beta, beta_censir),
    Dproj(Beta, beta_cox),
    Dproj(Beta, beta_cwlsir)
  )

  cor_rep <- c(
    abs(safe_cor(X %*% beta_sir, survival_times)),
    abs(safe_cor(X %*% beta_censir, survival_times)),
    abs(safe_cor(X %*% beta_cox, survival_times)),
    abs(safe_cor(X %*% beta_cwlsir, survival_times))
  )

  time_rep <- c(t1, t2, t3, t4)

  list(abscos = abscos_rep, Dist = dist_rep, Cor = cor_rep, Time = time_rep)
}

#' Etude de simulation Monte-Carlo : estimation EDR a index simple
#'
#' Repete \code{\link{one_rep_edr}} \code{num_repeats} fois pour chaque
#' taux de censure, en parallele (\code{foreach}/\code{doParallel}), et
#' moyenne les resultats.
#'
#' @inheritParams one_rep_edr
#' @param censor_rates Vecteur des taux de censure a tester.
#' @param num_repeats Nombre de repetitions par taux de censure.
#' @param n_cores Nombre de coeurs (1 = sequentiel).
#' @return Une liste : \code{results} (moyennes par taux de censure :
#'   \code{PC}, \code{COR}, \code{Frob}, \code{Time}), \code{raw}
#'   (matrices brutes 4 x num_repeats par taux de censure).
#' @export
run_simulation_edr <- function(n, p, q, Beta, censor_rates, num_repeats,
                                X_type, model_type,
                                alpha = NULL, spik = NULL, rho = NULL,
                                n_cores = max(1, parallel::detectCores() - 1),
                                n.slice, n.slice1, n.slice0) {
  results <- list()
  raw_results <- list()

  for (censor_rate in censor_rates) {
    cat("-> Censoring rate:", censor_rate * 100, "%\n")

    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    rep_id <- NULL
    rep_results <- foreach::foreach(
      rep_id = 1:num_repeats,
      .combine = "list",
      .multicombine = TRUE,
      .packages = c("survival", "MASS", "CSIRWLS")
    ) %dopar% {
      set.seed(rep_id + round(10000 * censor_rate))
      one_rep_edr(
        rep_id = rep_id, n = n, p = p, q = q, Beta = Beta,
        censor_rate = censor_rate, X_type = X_type, model_type = model_type,
        rho = rho, alpha = alpha, spik = spik,
        n.slice = n.slice, n.slice1 = n.slice1, n.slice0 = n.slice0
      )
    }

    parallel::stopCluster(cl)

    if (num_repeats == 1) rep_results <- list(rep_results)

    abscos <- matrix(sapply(rep_results, function(r) r$abscos), nrow = 4, ncol = num_repeats)
    Dist <- matrix(sapply(rep_results, function(r) r$Dist), nrow = 4, ncol = num_repeats)
    Cor <- matrix(sapply(rep_results, function(r) r$Cor), nrow = 4, ncol = num_repeats)
    Time <- matrix(sapply(rep_results, function(r) r$Time), nrow = 4, ncol = num_repeats)

    results[[paste0("censor_rate_", censor_rate)]] <- list(
      PC = round(rowMeans(abscos, na.rm = TRUE), 3),
      COR = round(rowMeans(Cor, na.rm = TRUE), 3),
      Frob = round(rowMeans(Dist, na.rm = TRUE), 3),
      Time = round(rowMeans(Time, na.rm = TRUE), 4)
    )

    raw_results[[paste0("censor_rate_", censor_rate)]] <- list(
      abscos = abscos, Cor = Cor, Dist = Dist, Time = Time
    )
  }

  return(list(results = results, raw = raw_results))
}

# ============================================================================
# 2) SCREENING DE VARIABLES : CWLS vs CSIS (TPR/FPR/FDR)
#    (source : CWLS-CSIS-bcorCISparalelle.R)
# ============================================================================

#' Une repetition de la comparaison de screening CWLS vs CSIS
#'
#' @param rep_id Identifiant de repetition.
#' @param n,p,q Taille d'echantillon, dimension, nombre de variables actives.
#' @param Beta Vecteur p de coefficients vrais (0 hors variables actives).
#' @param s Indices des vraies variables actives.
#' @param censor_rate Taux de censure cible.
#' @param X_type,model_type,rho,alpha,spik Voir \code{\link{generate_X}}.
#' @return Une liste : \code{TPR}, \code{FPR}, \code{FDR} (un couple
#'   CWLS/CSIS chacun), \code{Time}.
#' @export
one_rep_screening <- function(rep_id, n, p, q, Beta, s, censor_rate,
                               X_type, model_type, rho, alpha, spik) {
  X <- generate_X(n, p, q, type = X_type, rho = rho, alpha = alpha, spik = spik)

  survival_times <- generate_survival_times(X, Beta, model_type)

  errc <- rnorm(n)
  c_val <- find_c_for_censoring_rate(censor_rate, survival_times, errc)
  censoring_times <- exp(c_val + errc)

  time <- pmin(survival_times, censoring_times)
  event <- as.numeric(survival_times <= censoring_times)

  t1 <- system.time({
    cwls_select <- apply_cwls(X, time, event)
  })["elapsed"]

  t2 <- system.time({
    csis_select <- apply_csis(X, time, event)
  })["elapsed"]

  m_cwls <- compute_metrics(cwls_select, s, p)
  m_csis <- compute_metrics(csis_select, s, p)

  list(
    TPR = c(m_cwls$TPR, m_csis$TPR),
    FPR = c(m_cwls$FPR, m_csis$FPR),
    FDR = c(m_cwls$FDR, m_csis$FDR),
    Time = c(t1, t2)
  )
}

#' Etude de simulation Monte-Carlo : screening CWLS vs CSIS
#'
#' Repete \code{\link{one_rep_screening}} \code{num_repeats} fois pour
#' chaque taux de censure, en parallele.
#'
#' @inheritParams one_rep_screening
#' @param censor_rates Vecteur des taux de censure a tester.
#' @param num_repeats Nombre de repetitions par taux de censure.
#' @param n_cores Nombre de coeurs (1 = sequentiel).
#' @return Une liste : \code{results} (moyennes par taux de censure),
#'   \code{TPR}, \code{FPR}, \code{FDR}, \code{Time} (matrices combinees
#'   sur tous les taux de censure), \code{Time_by_rate} (liste par taux
#'   de censure, pour \code{\link{plot_time_by_censor}}), \code{censor_rates}.
#' @export
run_simulation_screening <- function(n, p, q, Beta, s, censor_rates,
                                      num_repeats, X_type, model_type,
                                      alpha = NULL, spik = NULL, rho = NULL,
                                      n_cores = max(1, parallel::detectCores() - 1)) {
  results <- list()
  all_TPR <- list()
  all_FPR <- list()
  all_FDR <- list()
  all_Time <- list()

  for (censor_rate in censor_rates) {
    cat("-> Censoring rate:", censor_rate * 100, "%\n")

    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)

    rep_id <- NULL
    rep_results <- foreach::foreach(
      rep_id = 1:num_repeats,
      .combine = "list",
      .multicombine = TRUE,
      .packages = c("survival", "CSIRWLS")
    ) %dopar% {
      set.seed(rep_id + round(10000 * censor_rate))
      one_rep_screening(
        rep_id = rep_id, n = n, p = p, q = q, Beta = Beta, s = s,
        censor_rate = censor_rate, X_type = X_type, model_type = model_type,
        rho = rho, alpha = alpha, spik = spik
      )
    }

    parallel::stopCluster(cl)

    if (num_repeats == 1) rep_results <- list(rep_results)

    TPR <- matrix(sapply(rep_results, function(r) r$TPR), nrow = 2, ncol = num_repeats)
    FPR <- matrix(sapply(rep_results, function(r) r$FPR), nrow = 2, ncol = num_repeats)
    FDR <- matrix(sapply(rep_results, function(r) r$FDR), nrow = 2, ncol = num_repeats)
    Time <- matrix(sapply(rep_results, function(r) r$Time), nrow = 2, ncol = num_repeats)

    results[[paste0("censor_rate_", censor_rate)]] <- list(
      TPR = round(rowMeans(TPR), 2),
      FPR = round(rowMeans(FPR), 2),
      FDR = round(rowMeans(FDR), 2),
      Time = round(rowMeans(Time), 4)
    )

    all_TPR[[paste0("censor_rate_", censor_rate)]] <- TPR
    all_FPR[[paste0("censor_rate_", censor_rate)]] <- FPR
    all_FDR[[paste0("censor_rate_", censor_rate)]] <- FDR
    all_Time[[paste0("censor_rate_", censor_rate)]] <- Time
  }

  TPR_all <- do.call(cbind, all_TPR)
  FPR_all <- do.call(cbind, all_FPR)
  FDR_all <- do.call(cbind, all_FDR)
  Time_all <- do.call(cbind, all_Time)

  list(
    results = results,
    TPR = TPR_all, FPR = FPR_all, FDR = FDR_all, Time = Time_all,
    Time_by_rate = all_Time, censor_rates = censor_rates
  )
}

#' Etude de l'impact de la structure de correlation par blocs (screening)
#'
#' Repete \code{\link{run_simulation_screening}} pour plusieurs
#' parametrages \code{alpha} (structure de correlation par blocs,
#' \code{X_type = "block"}), et combine les tables de resultats.
#'
#' @param alpha_list Liste nommee de vecteurs \code{alpha} (voir
#'   \code{\link{X_rand}}).
#' @param n,p,q,Beta,censor_rates,num_repeats,model_type,n_cores Voir
#'   \code{\link{run_simulation_screening}}.
#' @return Un data.frame combinant tous les parametrages alpha.
#' @export
run_alpha_study <- function(alpha_list, n, p, q, Beta,
                             censor_rates, num_repeats, model_type,
                             n_cores = max(1, parallel::detectCores() - 1)) {
  all_tables <- list()
  method_names <- c("CWLS", "CSIS")

  for (alpha_name in names(alpha_list)) {
    cat("\n============================\n")
    cat("Alpha setting:", alpha_name, "\n")
    cat("============================\n")

    alpha <- alpha_list[[alpha_name]]

    out <- run_simulation_screening(
      n = n, p = p, q = q, Beta = Beta, s = seq_len(q),
      censor_rates = censor_rates, num_repeats = num_repeats,
      X_type = "block", model_type = model_type,
      alpha = alpha, n_cores = n_cores
    )

    tab <- data.frame()
    for (cr in censor_rates) {
      res <- out$results[[paste0("censor_rate_", cr)]]
      tab <- rbind(tab, data.frame(
        AlphaSetting = alpha_name,
        Alpha = paste(alpha, collapse = ", "),
        CensorRate = paste0(100 * cr, "%"),
        Method = method_names,
        TPR = res$TPR,
        FPR = res$FPR,
        FDR = res$FDR,
        Time_s = res$Time
      ))
    }

    all_tables[[alpha_name]] <- tab
  }

  do.call(rbind, all_tables)
}

# ============================================================================
# 3) ESTIMATION EDR, INDEX MULTIPLE (d > 1) : SIR vs CSIR vs Cox vs CWLS_SIR
#    (source : CSIRWLS-Multiple indice.R)
# ============================================================================

#' Une repetition de la comparaison EDR a index multiple (dimension fixee)
#'
#' Meme principe que \code{\link{one_rep_edr}}, generalise a d > 1
#' directions EDR (Beta est une matrice p x d_true), avec
#' \code{\link{trace_corr}}, \code{\link{subspace_dist}} et
#' \code{\link{safe_cor_multi}} comme criteres de comparaison.
#'
#' @param rep_id Identifiant de repetition.
#' @param n,p,d_true Taille d'echantillon, dimension, nombre vrai de
#'   directions EDR.
#' @param Beta Matrice p x d_true de directions EDR vraies.
#' @param censor_rate Taux de censure cible.
#' @param X_type,model_type,rho,alpha,spik Voir \code{\link{generate_X}}.
#' @param n.slice,n.slice1,n.slice0 Parametres de tranchage.
#' @return Une liste : \code{TC}, \code{SD}, \code{Cor} (ordre
#'   SIR/CSIR/COX/CWLS_SIR), \code{Time}.
#' @export
one_rep_multi <- function(rep_id, n, p, d_true, Beta, censor_rate,
                           X_type, model_type, rho, alpha, spik,
                           n.slice, n.slice1, n.slice0) {
  if (!requireNamespace("Rdimtools", quietly = TRUE)) {
    stop("Le package 'Rdimtools' est necessaire pour one_rep_multi() (methode SIR de reference).")
  }

  X <- generate_X(n, p, type = X_type, rho = rho, alpha = alpha, spik = spik)
  Y0 <- generate_survival_times(X, Beta, model_type)
  errc <- rnorm(n)
  c0 <- find_c_for_censoring_rate(censor_rate, Y0, errc)
  C <- exp(c0 + errc)
  time <- pmin(Y0, C)
  event <- as.numeric(Y0 <= C)

  beta_init <- prcomp(X, rank. = 1)$rotation[, 1]
  h <- bw.nrd(as.vector(X %*% beta_init))

  # -- SIR --
  t1 <- system.time({
    sir_out <- Rdimtools::do.sir(as.matrix(X), time, ndim = d_true)
    beta_sir <- as.matrix(sir_out$projection)
  })["elapsed"]

  # -- CSIR --
  t2 <- system.time({
    if (p > n) {
      beta_csir <- matrix(0, p, d_true)
    } else {
      ds <- double.slice(time, event, X, n.slice1, n.slice0)
      joint.edrs <- edr.n(ds, d_true)
      csir_out <- cen.sir(time, event, X, n.slice, joint.edrs, h, c = 0.05)
      beta_csir <- as.matrix(csir_out$evec[, 1:d_true, drop = FALSE])
    }
  })["elapsed"]

  # -- COX (une seule direction ; orthogonalisee pour former une base d_true) --
  t3 <- system.time({
    if (p < n) {
      fit <- survival::coxph(survival::Surv(time, event) ~ X)
      beta_cox <- matrix(0, p, d_true)
      beta_cox[, 1] <- fit$coef
      if (d_true > 1) {
        v <- rnorm(p)
        v <- v - sum(v * beta_cox[, 1]) / sum(beta_cox[, 1]^2) * beta_cox[, 1]
        beta_cox[, 2] <- v / sqrt(sum(v^2))
      }
      beta_cox <- qr.Q(qr(beta_cox))
    } else {
      beta_cox <- matrix(0, p, d_true)
    }
  })["elapsed"]

  # -- CSIR-WLS --
  t4 <- system.time({
    cwls_out <- cen.wls.sir2(
      X, time, event,
      n.slice1, n.slice0, n.slice,
      cn1 = 0.1, cn2 = 1,
      choose.dir = FALSE, c = 0.05,
      ndim = d_true
    )
    beta_cwlsir <- as.matrix(cwls_out$beta.hat)
  })["elapsed"]

  csir_ok <- !all(beta_csir == 0)
  cox_ok <- !all(beta_cox == 0)

  tc_rep <- c(
    trace_corr(Beta, beta_sir),
    ifelse(csir_ok, trace_corr(Beta, beta_csir), 0),
    ifelse(cox_ok, trace_corr(Beta, beta_cox), 0),
    trace_corr(Beta, beta_cwlsir)
  )

  sd_rep <- c(
    subspace_dist(Beta, beta_sir),
    ifelse(csir_ok, subspace_dist(Beta, beta_csir), 1),
    ifelse(cox_ok, subspace_dist(Beta, beta_cox), 1),
    subspace_dist(Beta, beta_cwlsir)
  )

  cor_rep <- c(
    safe_cor_multi(X, Beta, beta_sir),
    ifelse(csir_ok, safe_cor_multi(X, Beta, beta_csir), 0),
    ifelse(cox_ok, safe_cor_multi(X, Beta, beta_cox), 0),
    safe_cor_multi(X, Beta, beta_cwlsir)
  )

  list(TC = tc_rep, SD = sd_rep, Cor = cor_rep, Time = c(t1, t2, t3, t4))
}

#' Etude de simulation Monte-Carlo : estimation EDR a index multiple
#'
#' @inheritParams one_rep_multi
#' @param censor_rates Vecteur des taux de censure a tester.
#' @param num_repeats Nombre de repetitions par taux de censure.
#' @param n_cores Nombre de coeurs (1 = sequentiel).
#' @return Une liste : \code{results} (moyennes \code{TC}/\code{SD}/
#'   \code{Cor}/\code{Time} par taux de censure), \code{raw}.
#' @export
run_simulation_multi <- function(n, p, d_true, Beta, censor_rates,
                                  num_repeats, X_type, model_type,
                                  alpha = NULL, spik = NULL, rho = NULL,
                                  n_cores = max(1, parallel::detectCores() - 1),
                                  n.slice, n.slice1, n.slice0) {
  results <- list()
  raw_results <- list()

  for (censor_rate in censor_rates) {
    cat("-> Censoring rate:", censor_rate * 100, "%\n")

    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)

    rep_id <- NULL
    rep_results <- foreach::foreach(
      rep_id = 1:num_repeats,
      .combine = "list",
      .multicombine = TRUE,
      .packages = c("survival", "MASS", "CSIRWLS")
    ) %dopar% {
      set.seed(rep_id + round(10000 * censor_rate))
      one_rep_multi(
        rep_id = rep_id, n = n, p = p, d_true = d_true, Beta = Beta,
        censor_rate = censor_rate, X_type = X_type, model_type = model_type,
        rho = rho, alpha = alpha, spik = spik,
        n.slice = n.slice, n.slice1 = n.slice1, n.slice0 = n.slice0
      )
    }

    parallel::stopCluster(cl)

    if (num_repeats == 1) rep_results <- list(rep_results)

    TC <- matrix(sapply(rep_results, function(r) r$TC), nrow = 4, ncol = num_repeats)
    SD <- matrix(sapply(rep_results, function(r) r$SD), nrow = 4, ncol = num_repeats)
    Cor <- matrix(sapply(rep_results, function(r) r$Cor), nrow = 4, ncol = num_repeats)
    Time <- matrix(sapply(rep_results, function(r) r$Time), nrow = 4, ncol = num_repeats)

    results[[paste0("censor_rate_", censor_rate)]] <- list(
      TC = round(rowMeans(TC, na.rm = TRUE), 3),
      SD = round(rowMeans(SD, na.rm = TRUE), 3),
      Cor = round(rowMeans(Cor, na.rm = TRUE), 3),
      Time = round(rowMeans(Time, na.rm = TRUE), 4)
    )

    raw_results[[paste0("censor_rate_", censor_rate)]] <- list(TC = TC, SD = SD, Cor = Cor, Time = Time)
  }

  return(list(results = results, raw = raw_results))
}

# ============================================================================
# 4) ESTIMATION EDR, INDEX MULTIPLE, DIMENSION d AUTOMATIQUE (BIC)
#    (source : CSIRWLS-Multiple indice.R)
# ============================================================================

#' Une repetition de la comparaison EDR a index multiple, avec d automatique
#'
#' Comme \code{\link{one_rep_multi}}, mais chaque methode estime elle
#' -meme sa dimension EDR d_hat (via \code{\link{select_d_bic}} pour
#' SIR/CSIR, via \code{\link{cen_wls_sir2_auto_d}} pour CSIR-WLS ; Cox
#' est par construction limitee a d=1), puis les beta estimes sont
#' ramenes a la dimension \code{d_true} via \code{\link{pad_beta}} pour
#' rester comparables au beta vrai.
#'
#' @inheritParams one_rep_multi
#' @return Une liste : \code{TC}, \code{SD}, \code{Cor}, \code{Time},
#'   \code{d_hat} (dimension estimee par chaque methode).
#' @export
one_rep_multi_auto_d <- function(rep_id, n, p, d_true, Beta, censor_rate,
                                  X_type, model_type, rho, alpha, spik,
                                  n.slice, n.slice1, n.slice0) {
  if (!requireNamespace("Rdimtools", quietly = TRUE)) {
    stop("Le package 'Rdimtools' est necessaire pour one_rep_multi_auto_d() (methode SIR de reference).")
  }

  X <- generate_X(n, p, type = X_type, rho = rho, alpha = alpha, spik = spik)
  Y0 <- generate_survival_times(X, Beta, model_type)
  errc <- rnorm(n)
  c0 <- find_c_for_censoring_rate(censor_rate, Y0, errc)
  C <- exp(c0 + errc)
  time <- pmin(Y0, C)
  event <- as.numeric(Y0 <= C)

  beta_init <- prcomp(X, rank. = 1)$rotation[, 1]
  h <- bw.nrd(as.vector(X %*% beta_init))

  beta_sir <- matrix(0, p, d_true)
  beta_csir <- matrix(0, p, d_true)
  beta_cox <- matrix(0, p, d_true)
  beta_cwlsir <- matrix(0, p, d_true)
  d_sir <- d_csir <- d_cox <- d_cwlsir <- NA
  t1 <- t2 <- t3 <- t4 <- NA

  # -- SIR --
  tryCatch(
    {
      t1 <- system.time({
        sir_full <- Rdimtools::do.sir(as.matrix(X), time, ndim = n.slice - 1)
        sir_proj <- as.matrix(Re(sir_full$projection))

        sir_evals <- if (!is.null(sir_full$eigenvalues)) {
          as.numeric(Re(sir_full$eigenvalues))
        } else {
          apply(sir_proj, 2, function(v) sum(v^2))
        }

        d_sir <- select_d_bic(sir_evals, n)
        d_sir <- max(1L, min(as.integer(d_sir), ncol(sir_proj)))
        beta_sir <- pad_beta(sir_proj[, 1:d_sir, drop = FALSE], d_true)
      })["elapsed"]
    },
    error = function(e) cat("  SIR failed:", conditionMessage(e), "\n")
  )

  # -- CSIR --
  tryCatch(
    {
      t2 <- system.time({
        if (p > n) {
          d_csir <- NA
        } else {
          ds <- double.slice(time, event, X, n.slice1, n.slice0)
          joint.edrs <- edr.n(ds, 2)
          csir_out <- cen.sir(time, event, X, n.slice, joint.edrs, h, c = 0.05)
          csir_evals <- as.numeric(Re(csir_out$eval))
          d_csir <- select_d_bic(csir_evals, n)
          d_csir <- max(1L, min(as.integer(d_csir), ncol(csir_out$evec)))
          beta_csir <- pad_beta(Re(csir_out$evec[, 1:d_csir, drop = FALSE]), d_true)
        }
      })["elapsed"]
    },
    error = function(e) cat("  CSIR failed:", conditionMessage(e), "\n")
  )

  # -- COX (une seule direction par construction) --
  tryCatch(
    {
      t3 <- system.time({
        if (p < n) {
          fit <- survival::coxph(survival::Surv(time, event) ~ X)
          coef_cox <- as.numeric(fit$coef)
          d_cox <- 1L
          beta_cox <- pad_beta(matrix(coef_cox, ncol = 1), d_true)
        } else {
          d_cox <- NA
        }
      })["elapsed"]
    },
    error = function(e) cat("  COX failed:", conditionMessage(e), "\n")
  )

  # -- CSIR-WLS (d automatique) --
  tryCatch(
    {
      t4 <- system.time({
        cwls_out <- cen_wls_sir2_auto_d(
          X, time, event,
          n.slice1, n.slice0, n.slice,
          cn1 = 0.1, cn2 = 1, c = 0.05,
          d_max = d_true + 3
        )
        d_cwlsir <- cwls_out$d_hat
        beta_cwlsir <- pad_beta(cwls_out$beta.hat, d_true)
      })["elapsed"]
    },
    error = function(e) cat("  CSIR-WLS failed:", conditionMessage(e), "\n")
  )

  csir_ok <- !all(beta_csir == 0)
  cox_ok <- !all(beta_cox == 0)
  cwls_ok <- !all(beta_cwlsir == 0)
  sir_ok <- !all(beta_sir == 0)

  tc_rep <- c(
    ifelse(sir_ok, trace_corr(Beta, beta_sir), 0),
    ifelse(csir_ok, trace_corr(Beta, beta_csir), 0),
    ifelse(cox_ok, trace_corr(Beta, beta_cox), 0),
    ifelse(cwls_ok, trace_corr(Beta, beta_cwlsir), 0)
  )

  sd_rep <- c(
    ifelse(sir_ok, subspace_dist(Beta, beta_sir), 1),
    ifelse(csir_ok, subspace_dist(Beta, beta_csir), 1),
    ifelse(cox_ok, subspace_dist(Beta, beta_cox), 1),
    ifelse(cwls_ok, subspace_dist(Beta, beta_cwlsir), 1)
  )

  cor_rep <- c(
    ifelse(sir_ok, safe_cor_multi(X, Beta, beta_sir), 0),
    ifelse(csir_ok, safe_cor_multi(X, Beta, beta_csir), 0),
    ifelse(cox_ok, safe_cor_multi(X, Beta, beta_cox), 0),
    ifelse(cwls_ok, safe_cor_multi(X, Beta, beta_cwlsir), 0)
  )

  list(
    TC = tc_rep, SD = sd_rep, Cor = cor_rep, Time = c(t1, t2, t3, t4),
    d_hat = c(d_sir, d_csir, d_cox, d_cwlsir)
  )
}

#' Etude de simulation Monte-Carlo : EDR index multiple, d automatique
#'
#' @inheritParams run_simulation_multi
#' @return Une liste : \code{results} (avec en plus \code{d_hat} moyen
#'   par methode), \code{raw}.
#' @export
run_simulation_multi_auto_d <- function(n, p, d_true, Beta, censor_rates,
                                         num_repeats, X_type, model_type,
                                         alpha = NULL, spik = NULL, rho = NULL,
                                         n_cores = max(1, parallel::detectCores() - 1),
                                         n.slice, n.slice1, n.slice0) {
  results <- list()
  raw_results <- list()

  for (censor_rate in censor_rates) {
    cat("-> Censoring rate:", censor_rate * 100, "%\n")

    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)

    rep_id <- NULL
    rep_results <- foreach::foreach(
      rep_id = 1:num_repeats,
      .combine = "list",
      .multicombine = TRUE,
      .packages = c("survival", "MASS", "CSIRWLS")
    ) %dopar% {
      set.seed(rep_id + round(10000 * censor_rate))
      one_rep_multi_auto_d(
        rep_id = rep_id, n = n, p = p, d_true = d_true, Beta = Beta,
        censor_rate = censor_rate, X_type = X_type, model_type = model_type,
        rho = rho, alpha = alpha, spik = spik,
        n.slice = n.slice, n.slice1 = n.slice1, n.slice0 = n.slice0
      )
    }

    parallel::stopCluster(cl)
    if (num_repeats == 1) rep_results <- list(rep_results)

    TC <- matrix(sapply(rep_results, function(r) r$TC), nrow = 4, ncol = num_repeats)
    SD <- matrix(sapply(rep_results, function(r) r$SD), nrow = 4, ncol = num_repeats)
    Cor <- matrix(sapply(rep_results, function(r) r$Cor), nrow = 4, ncol = num_repeats)
    Time <- matrix(sapply(rep_results, function(r) r$Time), nrow = 4, ncol = num_repeats)
    D_hat <- matrix(sapply(rep_results, function(r) r$d_hat), nrow = 4, ncol = num_repeats)

    results[[paste0("censor_rate_", censor_rate)]] <- list(
      TC = round(rowMeans(TC, na.rm = TRUE), 3),
      SD = round(rowMeans(SD, na.rm = TRUE), 3),
      Cor = round(rowMeans(Cor, na.rm = TRUE), 3),
      Time = round(rowMeans(Time, na.rm = TRUE), 4),
      d_hat = round(rowMeans(D_hat, na.rm = TRUE), 2)
    )

    raw_results[[paste0("censor_rate_", censor_rate)]] <- list(
      TC = TC, SD = SD, Cor = Cor, Time = Time, D_hat = D_hat
    )
  }

  return(list(results = results, raw = raw_results))
}
