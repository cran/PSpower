expit <- function(x) ifelse(x > 0, 1 / (1 + exp(-x)), exp(x) / (1 + exp(x)))
log_expit <- function(x) ifelse(x > 0, -log1p(exp(-x)), x - log1p(exp(x)))
logit <- function(x) log(x) - log(1 - x)

dlogitnorm <- function(x, mu, sigma) {
  1 / (sigma * sqrt(2*pi) * x * (1 - x)) * exp(-(logit(x) - mu)^2 / (2 * sigma^2))
}

overlap <- function(r, k) {
  alpha <- k * r
  beta <- k * (1 - r)
  return(list(
    alpha = alpha, beta = beta,
    mu = digamma(alpha) - digamma(beta),
    sigma = sqrt(trigamma(alpha) + trigamma(beta)),
    log_overlap = lgamma(alpha + 0.5) - log(alpha)/2 - lgamma(alpha) +
      lgamma(beta + 0.5) - log(beta)/2 - lgamma(beta)
  ))
}

solve_overlap <- function(r, phi) {
  f_overlap <- function(k) {
    overlap(r, k)$log_overlap - log(phi)
  }
  L <- 0.5 / min(r, 1 - r) + 1e-8
  R <- 2 * L
  while (f_overlap(R) < 0) {R <- 2 * R}
  while (R - L > 1e-8) {
    M <- (L + R) / 2
    if (f_overlap(M) < 0) {L <- M}
    else {R <- M}
  }
  return (overlap(r, M))
}

#' @importFrom stats dnorm
#' @importFrom stats integrate
solve_parameters <- function(r, phi, rho_sq) {

  overlap_object <- solve_overlap(r, phi)
  mu_e <- overlap_object$mu
  sigma_sq_e <- overlap_object$sigma^2

  f_int <- function(u, d, z) {
    u^d * dnorm(u, mu_e, sqrt(sigma_sq_e)) * expit((-1 + 2 * z) * u)
  }

  I01 <- integrate(f_int, lower = -Inf, upper = Inf, d = 0, z = 1)$value
  I11 <- integrate(f_int, lower = -Inf, upper = Inf, d = 1, z = 1)$value
  I21 <- integrate(f_int, lower = -Inf, upper = Inf, d = 2, z = 1)$value
  I00 <- integrate(f_int, lower = -Inf, upper = Inf, d = 0, z = 0)$value
  I10 <- integrate(f_int, lower = -Inf, upper = Inf, d = 1, z = 0)$value
  I20 <- integrate(f_int, lower = -Inf, upper = Inf, d = 2, z = 0)$value

  E1_e <- I11 / I01
  V1_e <- I21 / I01 - E1_e^2
  a_1_sq <- rho_sq / V1_e
  sigma_sq <- 1 - rho_sq

  E0_e <- I10 / I00
  V0_e <- I20 / I00 - E0_e^2
  a_0_sq <- rho_sq / V0_e

  calc_V <- function(h) {
    E_We_h <- integrate(function(x) x * h(expit(x)) * dnorm(x, mu_e, sqrt(sigma_sq_e)),
                        lower = -Inf, upper = Inf)$value
    E_h <- integrate(function(x) h(expit(x)) * dnorm(x, mu_e, sqrt(sigma_sq_e)),
                     lower = -Inf, upper = Inf)$value
    mu_We_h <- E_We_h / E_h

    term1 <- integrate(
      function(x) (x - mu_We_h)^2 * h(expit(x))^2 * exp(-log_expit(x) + dnorm(x, mu_e, sqrt(sigma_sq_e), log = T)),
      lower = -Inf, upper = Inf
    )$value
    term0 <- integrate(
      function(x) (x - mu_We_h)^2 * h(expit(x))^2 * exp(-log_expit(-x) + dnorm(x, mu_e, sqrt(sigma_sq_e), log = T)),
      lower = -Inf, upper = Inf
    )$value
    termn1 <- integrate(
      function(x) h(expit(x))^2 * exp(-log_expit(x) + dnorm(x, mu_e, sqrt(sigma_sq_e), log = T)),
      lower = -Inf, upper = Inf
    )$value
    termn0 <- integrate(
      function(x) h(expit(x))^2 * exp(-log_expit(-x) + dnorm(x, mu_e, sqrt(sigma_sq_e), log = T)),
      lower = -Inf, upper = Inf
    )$value
    denom <- integrate(
      function(x) h(expit(x))^2 * dnorm(x, mu_e, sqrt(sigma_sq_e)),
      lower = -Inf, upper = Inf
    )$value
    return ((a_1_sq * term1 + a_0_sq * term0 + sigma_sq * termn0 + sigma_sq * termn1) / denom)
  }

  return (list(
    parameters = c(
      mu_e = mu_e, sigma_sq_e = sigma_sq_e,
      a_1_sq = a_1_sq, a_0_sq = a_0_sq,
      sigma_sq = sigma_sq
    ),
    V_ATE = calc_V(function(x) 1),
    V_ATT = calc_V(function(x) x),
    V_ATC = calc_V(function(x) 1 - x),
    V_ATO = calc_V(function(x) x * (1-x)),
    calc_V = calc_V
  ))
}
