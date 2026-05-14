# Internal utilities: Beta/propensity parameterization, variance components.
# None of these are exported.

# --------------------------------------------------------------------------
# Numerically stable logistic and log-logistic functions
# --------------------------------------------------------------------------

.expit <- function(x) ifelse(x > 0, 1 / (1 + exp(-x)), exp(x) / (1 + exp(x)))

.log_expit <- function(x) ifelse(x > 0, -log1p(exp(-x)), x - log1p(exp(x)))


# --------------------------------------------------------------------------
# Beta / propensity score parameterization
# --------------------------------------------------------------------------

# Bhattacharyya overlap coefficient from Beta(a, b) parameters.
# Computed in log scale to avoid overflow.
.phi_from_ab <- function(a, b) {
  exp(
    lgamma(a + 0.5) - 0.5 * log(a) - lgamma(a) +
    lgamma(b + 0.5) - 0.5 * log(b) - lgamma(b)
  )
}

# Solve Beta parameters (a, b) from treatment proportion r and overlap phi.
# Uses the constraint b = a*(1-r)/r, so only a needs to be found.
# Searches on log(a) scale; uniroot bisects to tolerance 1e-12.
.solve_ab <- function(r, phi) {
  f <- function(log_a) {
    a <- exp(log_a)
    .phi_from_ab(a, a * (1 - r) / r) - phi
  }
  lo <- log(1e-8)
  hi <- 0
  while (f(hi) < 0) hi <- hi + log(2)
  log_a <- stats::uniroot(f, c(lo, hi), tol = 1e-12)$root
  a <- exp(log_a)
  list(a = a, b = a * (1 - r) / r)
}

# Empirical overlap coefficient from fitted propensity scores and treatment vector.
.phi_from_ps <- function(ps, Z) {
  r <- mean(Z)
  mean(sqrt(ps * (1 - ps))) / sqrt(r * (1 - r))
}


# --------------------------------------------------------------------------
# Variance components for continuous / binary outcomes (PS-weighted Hajek)
# --------------------------------------------------------------------------

# Solve variance components given (r, phi, rho2).
# Returns V for ATE (closed form) and for ATT, ATC, ATO (numerical integration).
# Reference: Liu, Yang and Li (2026), Annals of Statistics.
.solve_variance_ps <- function(r, phi, rho2) {

  # Step 1: solve (mu_e, sigma_sq_e) from (r, phi) via bisection.
  # The logit-normal PS distribution is approximated by Beta(alpha, beta)
  # where alpha = k*r, beta = k*(1-r), and k is found by matching phi.
  f_k <- function(k) {
    lgamma(k * r + 0.5) - 0.5 * log(k * r) - lgamma(k * r) +
      lgamma(k * (1 - r) + 0.5) - 0.5 * log(k * (1 - r)) - lgamma(k * (1 - r)) -
      log(phi)
  }
  lo <- 0.5 / min(r, 1 - r) + 1e-8
  hi <- 2 * lo
  while (f_k(hi) < 0) hi <- 2 * hi
  while (hi - lo > 1e-8) {
    mid <- (lo + hi) / 2
    if (f_k(mid) < 0) lo <- mid else hi <- mid
  }
  k_hat <- (lo + hi) / 2
  alpha_k <- k_hat * r
  beta_k  <- k_hat * (1 - r)
  mu_e       <- digamma(alpha_k) - digamma(beta_k)
  sigma_sq_e <- trigamma(alpha_k) + trigamma(beta_k)

  # Step 2: decompose outcome variance into PS-aligned and residual parts.
  a_sq    <- rho2 / sigma_sq_e   # squared slope of outcome on logit PS
  sigma_sq <- 1 - rho2            # residual variance (standardized)

  # Step 3: closed-form variance for ATE estimand.
  V_ATE <- 2 * (1 + (rho2 * sigma_sq_e + 1) *
                  exp(sigma_sq_e / 2) * cosh(abs(mu_e)))

  # Step 4: numerical variance for WATE with general tilting function h(e).
  # V_w = [a^2 * (A1 + A0) + sigma^2 * (B1 + B0)] / E[h(e)]^2
  # where Az, Bz involve integration of h(e)^2 against the arm-specific
  # PS density (logit-normal approximated by the Beta-matched normal).
  calc_V_w <- function(h) {
    dnorm_We <- function(x) stats::dnorm(x, mu_e, sqrt(sigma_sq_e))

    E_h    <- stats::integrate(function(x) h(.expit(x)) * dnorm_We(x),
                               -Inf, Inf)$value
    E_We_h <- stats::integrate(function(x) x * h(.expit(x)) * dnorm_We(x),
                               -Inf, Inf)$value
    mu_We_h <- E_We_h / E_h

    A1 <- stats::integrate(
      function(x) (x - mu_We_h)^2 * h(.expit(x))^2 *
        exp(.log_expit(-x) * (-1) + dnorm_We(x) * 0 +   # written out below
            stats::dnorm(x, mu_e, sqrt(sigma_sq_e), log = TRUE) -
            .log_expit(x)),
      -Inf, Inf)$value

    # Rewrite integrand cleanly:
    # For treated arm (z=1): h(e)^2 / e * f_We(x) = h(e)^2 * exp(-log_expit(x)) * dnorm
    # For control arm (z=0): h(e)^2 / (1-e) * f_We(x) = h(e)^2 * exp(-log_expit(-x)) * dnorm
    A1 <- stats::integrate(
      function(x) (x - mu_We_h)^2 * h(.expit(x))^2 *
        exp(-(.log_expit(x)) + stats::dnorm(x, mu_e, sqrt(sigma_sq_e), log = TRUE)),
      -Inf, Inf)$value
    A0 <- stats::integrate(
      function(x) (x - mu_We_h)^2 * h(.expit(x))^2 *
        exp(-(.log_expit(-x)) + stats::dnorm(x, mu_e, sqrt(sigma_sq_e), log = TRUE)),
      -Inf, Inf)$value
    B1 <- stats::integrate(
      function(x) h(.expit(x))^2 *
        exp(-(.log_expit(x)) + stats::dnorm(x, mu_e, sqrt(sigma_sq_e), log = TRUE)),
      -Inf, Inf)$value
    B0 <- stats::integrate(
      function(x) h(.expit(x))^2 *
        exp(-(.log_expit(-x)) + stats::dnorm(x, mu_e, sqrt(sigma_sq_e), log = TRUE)),
      -Inf, Inf)$value

    (a_sq * (A1 + A0) + sigma_sq * (B1 + B0)) / E_h^2
  }

  list(
    mu_e       = mu_e,
    sigma_sq_e = sigma_sq_e,
    V_ATE = V_ATE,
    V_ATT = calc_V_w(function(e) e),
    V_ATC = calc_V_w(function(e) 1 - e),
    V_ATO = calc_V_w(function(e) e * (1 - e)),
    calc_V_w = calc_V_w
  )
}


# --------------------------------------------------------------------------
# Variance components for survival outcomes (PS-weighted marginal Cox)
# --------------------------------------------------------------------------

# Variance of the weighted partial likelihood estimator in a randomized trial.
# Reference: Yang, Liu and Li (2026), arXiv:2605.10088.
.V_RCT <- function(tau0, r, d1, d0) {
  lam1 <- sqrt(r / (1 - r)) * exp(tau0 / 2)
  lam0 <- 1 / lam1
  d    <- r * d1 + (1 - r) * d0
  (lam1 + lam0)^2 * (r * lam0^2 * d1 + (1 - r) * lam1^2 * d0) / d^2
}

# Variance of the normalized IPW estimator in an observational study.
# Requires a > 1 and b > 1 (finite mean of inverse propensity scores).
.V_IPW <- function(tau0, r, a, b, d1, d0) {
  lam1      <- sqrt(r / (1 - r)) * exp(tau0 / 2)
  lam0      <- 1 / lam1
  d         <- r * d1 + (1 - r) * d0
  E_inv_e   <- (a + b - 1) / (a - 1)
  E_inv_1e  <- (a + b - 1) / (b - 1)
  (lam1 + lam0)^2 *
    (d1 * r^2 * lam0^2 * E_inv_e + d0 * (1 - r)^2 * lam1^2 * E_inv_1e) / d^2
}

# Schoenfeld variance (valid only for randomized trials).
.V_Schoenfeld <- function(r, d1, d0) {
  d <- r * d1 + (1 - r) * d0
  1 / (r * (1 - r) * d)
}

# Propensity score weights for a given scheme.
# Returns the vector of per-unit weights (actual, not potential).
.compute_weights <- function(Z, e, scheme) {
  if (scheme == "nipw") {
    r_hat <- mean(Z)
    w1 <- r_hat / e
    w0 <- (1 - r_hat) / (1 - e)
  } else if (scheme == "treated") {
    w1 <- rep(1, length(Z))
    w0 <- e / (1 - e)
  } else if (scheme == "ow") {
    w1 <- (1 - e) / sum(Z * (1 - e))
    w0 <- e / sum((1 - Z) * e)
  } else {
    stop("scheme must be one of: nipw, treated, ow")
  }
  Z * w1 + (1 - Z) * w0
}

# Monte Carlo estimate of the Kish design effect for a given weighting scheme.
# Samples n_mc (e_i, Z_i) pairs from Beta(a,b) and evaluates the Kish ratio.
# Used for ATO (scheme = "ow") and ATT (scheme = "treated").
.kappa_Kish <- function(r, phi, scheme, n_mc = 1e6) {
  ab <- .solve_ab(r, phi)
  e  <- stats::rbeta(n_mc, ab$a, ab$b)
  Z  <- stats::rbinom(n_mc, 1, e)
  w  <- .compute_weights(Z, e, scheme)
  n1 <- sum(Z);       n0 <- n_mc - n1
  sw1  <- sum(Z * w);       sw12 <- sum(Z * w^2)
  sw0  <- sum((1 - Z) * w); sw02 <- sum((1 - Z) * w^2)
  (n1 * n0 / n_mc) * (sw12 / sw1^2 + sw02 / sw0^2)
}


# --------------------------------------------------------------------------
# Shared sample size / power inversion
# --------------------------------------------------------------------------

# Required N from variance V, effect size delta, alpha, and power.
.N_from_V <- function(V, delta, sig_level, power, test) {
  z_alpha <- if (test == "two-sided") stats::qnorm(1 - sig_level / 2)
             else                     stats::qnorm(1 - sig_level)
  z_beta  <- stats::qnorm(power)
  ceiling((z_alpha + z_beta)^2 * V / delta^2)
}

# Achieved power from variance V, effect size delta, N, alpha.
.power_from_V <- function(V, delta, N, sig_level, test) {
  z_alpha <- if (test == "two-sided") stats::qnorm(1 - sig_level / 2)
             else                     stats::qnorm(1 - sig_level)
  stats::pnorm(sqrt(N / V) * abs(delta) - z_alpha)
}


# --------------------------------------------------------------------------
# Weight label helpers for print / summary
# --------------------------------------------------------------------------

.weight_label <- function(estimand) {
  switch(estimand,
    ATE = "inverse probability weights",
    ATT = "weights for group 1 (treated population)",
    ATC = "weights for group 0 (control population)",
    ATO = "overlap weights",
    "custom tilting function"
  )
}

# Human-readable axis/legend label for a design parameter column.
# legend = TRUE wraps two-word names with \n for right-side placement.
.col_label <- function(col, legend = FALSE) {
  br <- if (legend) "\n" else " "
  switch(col,
    effect_size = "Effect size",
    r           = paste0("Treatment", br, "proportion"),
    phi         = paste0("Overlap", br, "coefficient"),
    rho2        = paste0("Confounder", br, "coefficient"),
    d1          = paste0("Event rate", br, "(group 1)"),
    d0          = paste0("Event rate", br, "(group 0)"),
    study_type  = "Study type",
    estimand    = "Estimand",
    method      = "Method",
    sample_size = "Required sample size",
    power       = "Achieved power",
    col
  )
}
