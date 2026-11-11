## Table S20: population linearity defect (28) of the non-elliptical designs and
## its effect on the population ranking (Proposition 4.3).
## DGP 5.1 (log y0 = beta'x + eps), AR(1) rho = 0.5, H = 5 slices of y0 at its
## population quantiles, p in {80, 300, 1500}.
##
## Because u is Gaussian and independent of (z0, W, eps), E(x | slice) lies in
## span{1, Sigma_rho beta}:
##   skew-normal  E(x|s) = delta E(|z0| | s) 1 + c E(g | s)   Sigma_rho beta / (beta' Sigma_rho beta)
##   skew-t       E(x|s) = delta E(R|z0| | s) 1 + c E(R g | s) Sigma_rho beta / (beta' Sigma_rho beta)
## with g = beta'u, c = sqrt(1 - delta^2), R = (W / nu)^(-1/2); the scalars are
## computed by Monte Carlo (M draws). Then m_l = Sigma^{-1/2}(E(x|s_l) - mu),
## theta = Sigma^{1/2} beta, xi_l = m_l - P_theta m_l, and
##   xi_rel = max_l ||xi_l|| / ||m_l||,
##   |c|    = max_l |mean_{j inactive} xi_lj|,  sd(nu) = max_l sd_{j inactive} xi_lj,
##   Delta(lin) = min_{active} omega(lin) - max_{inactive} omega(lin),
##   omega_j = sum_l p_l m_lj^2,  omega_j(lin) = sum_l p_l (P_theta m_l)_j^2.
##
## beta: the five active coefficients are drawn once from U(0.5, 1.5) with a
## fixed seed (the same vector for the three dimensions); another vector can be
## given with CSIRWLS_S20_BETA="b1,b2,b3,b4,b5".
## The defect (28) is computed for z = Sigma^{-1/2} x, Sigma the covariance of
## the design (the convention of (C.1) and Proposition 4.1).
## CSIRWLS_S20_WHITEN = "Sigma" (default), "Srho" (whiten with the AR(1)
## reference covariance) or "I" (no whitening) changes the scale.
source(file.path("reproduction", "lib", "common.R"))
M <- as.integer(Sys.getenv("CSIRWLS_S20_M", "4000000"))
bstr <- Sys.getenv("CSIRWLS_S20_BETA", "")
b5 <- if (nzchar(bstr)) as.numeric(strsplit(bstr, ",")[[1]]) else { set.seed(19L); stats::runif(5, 0.5, 1.5) }
rho <- 0.5; delta <- 0.8; nu <- 5; cc <- sqrt(1 - delta^2)
WHITEN <- Sys.getenv("CSIRWLS_S20_WHITEN", "Sigma")

msqrt <- function(S, pw) { e <- eigen(S, symmetric = TRUE); e$vectors %*% (t(e$vectors) * e$values^pw) }

defect_cell <- function(kind, p) {
  spec <- if (kind == "gaussian") spec_ar1(rho) else spec_skew(kind, rho, delta, nu)
  beta <- numeric(p); beta[1:5] <- b5
  Srho <- rho^abs(outer(1:p, 1:p, "-"))
  Sb <- as.vector(Srho %*% beta); su2 <- sum(beta * Sb); Sbeta <- sum(beta)
  set.seed(hash_seed(paste("defect", kind, p)))
  g <- stats::rnorm(M, 0, sqrt(su2)); z0 <- abs(stats::rnorm(M)); eps <- stats::rnorm(M)
  R <- if (kind == "skewt") sqrt(nu / stats::rchisq(M, nu)) else rep(1, M)
  eta <- switch(kind, gaussian = g, R * (delta * Sbeta * z0 + cc * g))
  y0 <- eta + eps
  sl <- cut(y0, stats::quantile(y0, seq(0, 1, length.out = H_SLICES + 1)), include.lowest = TRUE, labels = FALSE)
  pl <- tabulate(sl, H_SLICES) / M
  a1 <- tapply(R * z0, sl, mean); a2 <- tapply(if (kind == "gaussian") g else R * g, sl, mean)
  ER <- if (kind == "skewt") sqrt(nu / 2) * gamma((nu - 1) / 2) / gamma(nu / 2) else 1
  Ex <- sapply(seq_len(H_SLICES), function(l) switch(kind,
    gaussian = Sb * a2[l] / su2,
    delta * a1[l] + cc * Sb * a2[l] / su2))                    # p x H
  mu <- if (kind == "gaussian") numeric(p) else rep(delta * sqrt(2 / pi) * ER, p)
  Sigma <- sigma_sub(spec, N_DESIGN, p, 1:p)
  W <- switch(WHITEN, Sigma = Sigma, Srho = Srho, I = diag(p))
  m <- msqrt(W, -0.5) %*% (Ex - mu)                             # p x H, the m_l
  theta <- as.vector(msqrt(W, 0.5) %*% beta)
  Pm <- theta %*% (crossprod(theta, m) / sum(theta^2))
  xi <- m - Pm
  act <- 1:5; ina <- setdiff(1:p, act)
  om <- as.vector((m^2) %*% pl); oml <- as.vector((Pm^2) %*% pl)
  data.frame(design = kind, p = p,
             xi_rel = max(sqrt(colSums(xi^2)) / sqrt(colSums(m^2))),
             c_abs = max(abs(colMeans(xi[ina, , drop = FALSE]))),
             sd_nu = max(apply(xi[ina, , drop = FALSE], 2, stats::sd)),
             Delta_lin = min(oml[act]) - max(oml[ina]),
             Delta = min(om[act]) - max(om[ina]),
             lambda_max = max(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values),
             beta = paste(round(b5, 4), collapse = ","), whiten = WHITEN, M = M)
}
res <- do.call(rbind, lapply(c("gaussian", "skewnormal", "skewt"), function(k)
  do.call(rbind, lapply(p_grid(c(80, 300, 1500)), function(p) { message("defect ", k, " p = ", p); defect_cell(k, p) }))))
dir.create(file.path(CFG$results, "defect"), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(reps = res), file.path(CFG$results, "defect", "defect.rds"))
print(res, digits = 4)
