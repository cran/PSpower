# PSpower 2.0.0

This is a major revision with a fully rewritten codebase, new theoretical
foundations, and substantially expanded functionality. The API is not
backward-compatible with version 0.1.1.

## New features

* `power_ps()` — sample size and power for the PS-weighted Hájek estimator
  with continuous or binary outcomes. Supports four estimands (ATE, ATT, ATC,
  ATO) via closed-form (ATE) or numerical integration (ATT, ATC, ATO, custom
  tilting functions). Accounts for the confounder coefficient ρ² and the
  Bhattacharyya overlap coefficient φ.

* `power_cox()` — sample size and power for the PS-weighted partial likelihood
  estimator in a Cox proportional hazards model with time-to-event outcomes.
  Supports randomized trials (robust sandwich variance or Schoenfeld formula)
  and observational studies (ATE via IPW; ATO and ATT via Monte Carlo
  design-effect adjustment).

* `overlap_coef()` — estimates the Bhattacharyya overlap coefficient φ from
  fitted propensity scores and a treatment indicator, or analytically from
  Beta distribution parameters.

* S3 `print()`, `summary()`, and `plot()` methods for both `power_ps` and
  `power_cox` result objects. Scalar inputs produce a formatted single-scenario
  summary; vector inputs produce a multi-scenario grid with a five-number
  distribution summary and a `ggplot2`-based sensitivity plot.

## Theoretical basis

* Continuous/binary outcomes: Liu, Yang and Li (2026), *Annals of Statistics*,
  forthcoming.
* Survival outcomes: Yang, Liu and Li (2026), arXiv:2605.10088.

## Changes from version 0.1.1

* The previous single `PSpower()` function has been replaced by `power_ps()`
  and `power_cox()`, covering a broader set of estimands and outcome types.
* Survival outcomes (time-to-event) are now supported.
* The overlap coefficient φ replaces the previous parameterization.
* The four WATE estimands (ATE, ATT, ATC, ATO) are all supported for
  continuous/binary outcomes; ATE, ATT, and ATO for survival.
