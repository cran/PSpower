## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 7,
  fig.height = 4.5
)
library(PSpower)

## ----ps-single----------------------------------------------------------------
res <- power_ps(
  effect_size = 0.2,
  r           = 0.5,
  phi         = 0.9,
  estimand    = "ATE",
  power       = 0.80
)
res

## ----ps-r---------------------------------------------------------------------
res_r <- power_ps(
  effect_size = 0.2,
  r           = c(0.30, 0.40, 0.50, 0.60, 0.70),
  phi         = 0.90,
  estimand    = c("ATE", "ATT", "ATO"),
  power       = 0.80
)
plot(res_r)

## ----ps-phi-------------------------------------------------------------------
res_phi <- power_ps(
  effect_size = 0.2,
  r           = 0.5,
  phi         = c(0.85, 0.90, 0.95),
  estimand    = c("ATE", "ATT", "ATO"),
  power       = 0.80
)
plot(res_phi)

## ----ps-phi-summary-----------------------------------------------------------
summary(res_phi)

## ----ps-rho2------------------------------------------------------------------
power_ps(
  effect_size = 0.2,
  r           = 0.5,
  phi         = 0.90,
  rho2        = c(0, 0.01, 0.02, 0.05),
  estimand    = "ATE",
  power       = 0.80
)

## ----cox-rct-single-----------------------------------------------------------
res_rct <- power_cox(
  effect_size = log(0.75),   # HR = 0.75
  r           = 0.5,
  d1          = 0.40,
  study_type  = "rct",
  power       = 0.80
)
res_rct

## ----cox-schoenfeld-----------------------------------------------------------
power_cox(
  effect_size = log(0.75),
  r           = 0.5,
  d1          = 0.40,
  study_type  = "rct",
  method      = "schoenfeld",
  power       = 0.80
)

## ----cox-rct-multi------------------------------------------------------------
res_rct_grid <- power_cox(
  effect_size = c(log(0.65), log(0.75), log(0.85)),
  r           = c(0.40, 0.50, 0.60),
  d1          = 0.40,
  study_type  = "rct",
  power       = 0.80
)
plot(res_rct_grid)

## ----cox-rct-summary----------------------------------------------------------
summary(res_rct_grid)

## ----cox-obs-single-----------------------------------------------------------
res_obs <- power_cox(
  effect_size = log(0.75),
  r           = 0.5,
  d1          = 0.40,
  phi         = 0.90,
  study_type  = "obs",
  estimand    = "ATE",
  power       = 0.80
)
res_obs

## ----cox-obs-multi------------------------------------------------------------
res_obs_grid <- power_cox(
  effect_size = log(0.75),
  r           = 0.5,
  d1          = 0.40,
  phi         = c(0.85, 0.90, 0.95),
  study_type  = "obs",
  estimand    = c("ATE", "ATO"),
  power       = 0.80,
  n_mc        = 5e4    # use the default 1e6 in practice
)
plot(res_obs_grid)

## ----cox-obs-summary----------------------------------------------------------
summary(res_obs_grid)

## ----overlap-empirical--------------------------------------------------------
set.seed(2026)
n  <- 500
X  <- rnorm(n)
ps <- plogis(0.5 * X)   # fitted propensity scores from a PS model
Z  <- rbinom(n, 1, ps)

overlap_coef(ps = ps, Z = Z)

## ----overlap-analytical-------------------------------------------------------
overlap_coef(a = 3, b = 2)   # r = a/(a+b) = 0.60

