#' Sample size and power for PS-weighted average treatment effect estimators
#'
#' Computes the required sample size or the achieved power for the
#' propensity-score-weighted Hajek estimator of a weighted average treatment
#' effect (WATE) with continuous or binary outcomes.
#'
#' @description
#' The required sample size is
#' \deqn{N = \bar{V}\,(z_{1-\alpha/2} + z_{\beta})^2 \,/\, \tilde{\tau}^2,}
#' where \eqn{\tilde{\tau}} is the standardized effect size and
#' \eqn{\bar{V}} is the asymptotic variance of the Hajek estimator.
#'
#' For the **ATE** estimand, \eqn{\bar{V}} has the closed form
#' \deqn{\bar{V} = 2\!\left\{1 + \bigl(\rho^2\sigma_e^2+1\bigr)
#'   \exp\!\bigl(\sigma_e^2/2\bigr)\cosh(\mu_e)\right\},}
#' where \eqn{(\mu_e,\,\sigma_e^2)} are uniquely determined by \eqn{(r,\phi)}.
#'
#' For the **ATT**, **ATC**, and **ATO** estimands, \eqn{\bar{V}} is computed
#' by numerical integration of the same variance expression with the
#' corresponding tilting function \eqn{h(e)}.  A custom tilting function may
#' also be supplied.
#'
#' For binary outcomes, the estimand is the risk difference; the same formula
#' applies with \eqn{S^2 = \mathrm{Var}(Y(0))} estimated from a linear
#' probability model.
#'
#' @param effect_size Standardized effect size \eqn{\tilde{\tau} = \tau / S},
#'   where \eqn{\tau} is the treatment effect and
#'   \eqn{S = \sqrt{\mathrm{Var}(Y(0))}}. Scalar or vector.
#' @param r Treatment proportion \eqn{r = \Pr(Z = 1)}, in \eqn{(0, 1)}.
#'   Scalar or vector.
#' @param phi Overlap coefficient \eqn{\phi \in (0, 1)}, measuring covariate
#'   similarity between groups.  Rule of thumb: \eqn{< 0.80} very poor,
#'   \eqn{[0.80, 0.90)} poor, \eqn{[0.90, 0.95)} moderate,
#'   \eqn{\ge 0.95} good.  Scalar or vector.
#' @param rho2 Confounding coefficient \eqn{\rho^2 \in [0, 1)}, the squared
#'   correlation between the potential outcome and the propensity score linear
#'   predictor. Bounded above by the \eqn{R^2} of regressing the outcome on
#'   covariates. Sensitivity analysis over \eqn{\rho^2 \in [0, 0.05)} is
#'   recommended. Default \code{0}. Scalar or vector.
#' @param estimand Target estimand.  One of \code{"ATE"} (average treatment
#'   effect), \code{"ATT"} (average treatment effect for group 1),
#'   \code{"ATC"} (average treatment effect for group 0), \code{"ATO"}
#'   (average treatment effect for the overlap population), or a custom
#'   tilting function \code{h(e)} (must be a scalar function, not vectorized
#'   with other parameters). Scalar or character vector.
#' @param sig_level Significance level \eqn{\alpha}. Default \code{0.05}.
#' @param power Target power \eqn{\beta}. Provide for sample size calculation;
#'   mutually exclusive with \code{sample_size}.
#' @param sample_size Total sample size \eqn{N}. Provide for power
#'   calculation; mutually exclusive with \code{power}.
#' @param test \code{"two-sided"} (default) or \code{"one-sided"}.
#'
#' @return An object of class \code{"power_ps"}, a list containing:
#'   \describe{
#'     \item{\code{call}}{The matched call.}
#'     \item{\code{calculation}}{\code{"sample_size"} or \code{"power"}.}
#'     \item{\code{result}}{A data frame with one row per scenario (all
#'       combinations of vector inputs) and columns for every design parameter
#'       plus the computed \code{sample_size} or \code{power}.}
#'     \item{\code{settings}}{A list with \code{sig_level}, \code{power},
#'       \code{sample_size}, and \code{test}.}
#'     \item{\code{n_scenarios}}{Number of rows in \code{result}.}
#'     \item{\code{rho2_is_default}}{Logical; \code{TRUE} when \code{rho2}
#'       was left at its default value of \code{0}.}
#'   }
#'
#' @references
#' Bo Liu, Chengxin Yang, and Fan Li. Sample size and power calculations for
#' causal inference in observational studies.
#' \emph{Annals of Statistics} (2026), forthcoming.
#'
#' @seealso \code{\link{power_cox}}, \code{\link{overlap_coef}}
#'
#' @examples
#' # Sample size for ATE, scalar inputs
#' power_ps(effect_size = 0.2, r = 0.5, phi = 0.9, power = 0.8)
#'
#' # Power at a fixed N
#' power_ps(effect_size = 0.2, r = 0.5, phi = 0.9, sample_size = 250)
#'
#' # Sensitivity over r and estimand (vector inputs)
#' power_ps(effect_size = 0.2, r = c(0.3, 0.5, 0.7), phi = 0.9,
#'          estimand = c("ATE", "ATO"), power = 0.8)
#'
#' @export
power_ps <- function(effect_size,
                     r,
                     phi,
                     rho2        = 0,
                     estimand    = "ATE",
                     sig_level   = 0.05,
                     power       = NULL,
                     sample_size = NULL,
                     test        = "two-sided") {

  # ---- Input validation ----
  if (is.null(power) && is.null(sample_size))
    stop("Exactly one of 'power' or 'sample_size' must be provided.")
  if (!is.null(power) && !is.null(sample_size))
    stop("Provide 'power' or 'sample_size', not both.")
  if (!is.null(power) && (length(power) != 1 || power <= 0 || power >= 1))
    stop("'power' must be a single number in (0, 1).")
  if (!is.null(sample_size) && (length(sample_size) != 1 || sample_size <= 0))
    stop("'sample_size' must be a single positive number.")
  if (length(sig_level) != 1 || sig_level <= 0 || sig_level >= 1)
    stop("'sig_level' must be a single number in (0, 1).")
  if (!test %in% c("two-sided", "one-sided"))
    stop("'test' must be \"two-sided\" or \"one-sided\".")
  if (any(r <= 0 | r >= 1))
    stop("All values of 'r' must be in (0, 1).")
  if (any(phi <= 0 | phi >= 1))
    stop("All values of 'phi' must be in (0, 1).")
  if (any(rho2 < 0 | rho2 >= 1))
    stop("All values of 'rho2' must be in [0, 1).")

  rho2_is_default <- isTRUE(all.equal(rho2, 0)) && missing(rho2)
  calculation     <- if (!is.null(power)) "sample_size" else "power"

  # ---- Handle custom tilting function ----
  custom_h   <- NULL
  char_estim <- NULL
  if (is.function(estimand)) {
    custom_h   <- estimand
    char_estim <- "custom"
    estimand   <- "custom"
  } else {
    valid <- c("ATE", "ATT", "ATC", "ATO")
    if (!all(estimand %in% valid))
      stop("'estimand' must be one of \"ATE\", \"ATT\", \"ATC\", \"ATO\", ",
           "or a tilting function. Custom functions cannot be vectorized.")
    char_estim <- estimand
  }

  # ---- Build scenario grid ----
  grid <- expand.grid(
    effect_size = effect_size,
    r           = r,
    phi         = phi,
    rho2        = rho2,
    estimand    = char_estim,
    stringsAsFactors = FALSE
  )

  # ---- Compute V and result for each scenario ----
  result_col <- vapply(seq_len(nrow(grid)), function(i) {
    es  <- grid$effect_size[i]
    ri  <- grid$r[i]
    pi  <- grid$phi[i]
    rh  <- grid$rho2[i]
    est <- grid$estimand[i]

    params <- .solve_variance_ps(ri, pi, rh)

    V <- if (est == "custom") {
      params$calc_V_w(custom_h)
    } else {
      switch(est,
        ATE = params$V_ATE,
        ATT = params$V_ATT,
        ATC = params$V_ATC,
        ATO = params$V_ATO
      )
    }

    if (calculation == "sample_size") {
      .N_from_V(V, es, sig_level, power, test)
    } else {
      .power_from_V(V, es, sample_size, sig_level, test)
    }
  }, numeric(1))

  grid[[calculation]] <- result_col

  structure(
    list(
      call           = match.call(),
      calculation    = calculation,
      result         = grid,
      settings       = list(
        sig_level   = sig_level,
        power       = power,
        sample_size = sample_size,
        test        = test
      ),
      n_scenarios    = nrow(grid),
      rho2_is_default = rho2_is_default
    ),
    class = "power_ps"
  )
}


# --------------------------------------------------------------------------
# S3 methods
# --------------------------------------------------------------------------

#' @export
print.power_ps <- function(x, ...) {
  cat("PSpower: power_ps\n")
  cat(strrep("-", 50), "\n")

  res    <- x$result
  s      <- x$settings
  single <- x$n_scenarios == 1L

  if (single) {
    est      <- res$estimand[1]
    rho2_tag <- if (x$rho2_is_default) "  [default]" else ""

    cat(sprintf("%-26s %s\n",    "Estimand:",             est))
    cat(sprintf("%-26s %s\n",    "Weights:",              .weight_label(est)))
    cat(sprintf("%-26s %.3f\n",  "Effect size:",          res$effect_size[1]))
    cat(sprintf("%-26s %.3f\n",  "Treatment proportion:", res$r[1]))
    cat(sprintf("%-26s %.3f\n",  "Overlap (phi):",        res$phi[1]))
    cat(sprintf("%-26s %.3f%s\n","Confounding (rho2):",   res$rho2[1], rho2_tag))
    cat(strrep("-", 50), "\n")
    if (x$calculation == "sample_size") {
      cat(sprintf("Required sample size: %d\n", res$sample_size[1]))
      cat(sprintf("Target power: %.2f  |  sig. level: %.3f  |  %s test\n",
                  s$power, s$sig_level, s$test))
    } else {
      cat(sprintf("Achieved power: %.3f\n", res$power[1]))
      cat(sprintf("Sample size: %d  |  sig. level: %.3f  |  %s test\n",
                  as.integer(s$sample_size), s$sig_level, s$test))
    }

  } else {
    design_cols <- c("effect_size", "r", "phi", "rho2", "estimand")
    vary_cols   <- design_cols[sapply(design_cols, function(v)
      length(unique(res[[v]])) > 1)]
    fixed_cols  <- setdiff(design_cols, vary_cols)

    cat(sprintf("Scenarios: %d\n", x$n_scenarios))
    if (length(vary_cols))
      cat("Varying: ", paste(vary_cols, collapse = ", "), "\n")
    if (length(fixed_cols)) {
      fixed_str <- paste(sapply(fixed_cols, function(v) {
        val      <- unique(res[[v]])
        rho2_tag <- if (v == "rho2" && x$rho2_is_default) " [default]" else ""
        paste0(v, " = ", format(val, digits = 3), rho2_tag)
      }), collapse = ",  ")
      cat("Fixed:   ", fixed_str, "\n")
    }
    pow_str <- if (!is.null(s$power))
      sprintf("Target power: %.2f  |  sig. level: %.3f  |  %s test",
              s$power, s$sig_level, s$test)
    else
      sprintf("Sample size: %d  |  sig. level: %.3f  |  %s test",
              as.integer(s$sample_size), s$sig_level, s$test)
    cat(pow_str, "\n")
    cat(strrep("-", 50), "\n")

    show_cols <- c(vary_cols, x$calculation)
    n_show    <- min(5L, x$n_scenarios)
    print(res[seq_len(n_show), show_cols, drop = FALSE], row.names = FALSE)
    if (x$n_scenarios > n_show)
      cat(sprintf("... %d more rows; use summary() for the full distribution.\n",
                  x$n_scenarios - n_show))
  }

  invisible(x)
}


#' @export
summary.power_ps <- function(object, ...) {
  cat("PSpower: power_ps\n")
  cat(strrep("-", 55), "\n")

  res     <- object$result
  s       <- object$settings
  single  <- object$n_scenarios == 1L
  val_col <- object$calculation

  design_cols <- c("effect_size", "r", "phi", "rho2", "estimand")
  vary_cols   <- design_cols[sapply(design_cols, function(v)
    length(unique(res[[v]])) > 1)]
  fixed_cols  <- setdiff(design_cols, vary_cols)

  if (!single) cat(sprintf("Scenarios: %d\n", object$n_scenarios))
  if (length(vary_cols))
    cat("Varying: ", paste(vary_cols, collapse = ", "), "\n")
  if (length(fixed_cols)) {
    fixed_str <- paste(sapply(fixed_cols, function(v) {
      val      <- unique(res[[v]])
      rho2_tag <- if (v == "rho2" && object$rho2_is_default) " [default]" else ""
      paste0(v, " = ", format(val, digits = 3), rho2_tag)
    }), collapse = ",  ")
    cat("Fixed:   ", fixed_str, "\n")
  }
  if (!is.null(s$power))
    cat(sprintf("Target power: %.2f  |  sig. level: %.3f  |  %s test\n",
                s$power, s$sig_level, s$test))
  else
    cat(sprintf("Sample size: %d  |  sig. level: %.3f  |  %s test\n",
                as.integer(s$sample_size), s$sig_level, s$test))
  cat(strrep("-", 55), "\n")

  if (single) {
    est      <- res$estimand[1]
    rho2_tag <- if (object$rho2_is_default) "  [default]" else ""
    cat(sprintf("%-26s %s\n",    "Estimand:",             est))
    cat(sprintf("%-26s %s\n",    "Weights:",              .weight_label(est)))
    cat(sprintf("%-26s %.3f\n",  "Effect size:",          res$effect_size[1]))
    cat(sprintf("%-26s %.3f\n",  "Treatment proportion:", res$r[1]))
    cat(sprintf("%-26s %.3f\n",  "Overlap (phi):",        res$phi[1]))
    cat(sprintf("%-26s %.3f%s\n","Confounding (rho2):",   res$rho2[1], rho2_tag))
    cat(strrep("-", 55), "\n")
    if (val_col == "sample_size")
      cat(sprintf("Required sample size: %d\n", res$sample_size[1]))
    else
      cat(sprintf("Achieved power: %.3f\n", res$power[1]))
    return(invisible(object))
  }

  # Multi-scenario: 5-number distribution summary
  vals       <- res[[val_col]]
  n          <- object$n_scenarios
  sorted_idx <- order(vals)
  pos        <- pmin(c(1L, floor(n * 0.25) + 1L, floor(n * 0.5) + 1L,
                       floor(n * 0.75) + 1L, n), n)
  sel        <- sorted_idx[pos]
  qlabels    <- c("Min", "25%", "50%", "75%", "Max")

  label <- if (val_col == "sample_size") "Sample size (N, ceiling integer)" else "Power"
  cat(sprintf("\n%s distribution:\n\n", label))

  hdr_cols <- c(val_col, vary_cols)
  cat(sprintf("  %-8s", ""))
  cat(paste(sprintf("%-14s", hdr_cols), collapse = ""), "\n")

  for (k in seq_along(qlabels)) {
    i        <- sel[k]
    vals_row <- sapply(hdr_cols, function(v) {
      val <- res[[v]][i]
      if (v == val_col && val_col == "sample_size") as.character(as.integer(val))
      else if (v == val_col) sprintf("%.3f", val)
      else if (is.numeric(val)) format(round(val, 3), nsmall = 0)
      else as.character(val)
    })
    cat(sprintf("  %-8s%s\n", paste0(qlabels[k], ":"),
                paste(sprintf("%-14s", vals_row), collapse = "")))
  }

  invisible(object)
}


#' @export
plot.power_ps <- function(x, x_var = NULL, group = NULL, facet = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for plotting.")

  res     <- x$result
  val_col <- x$calculation

  design_cols <- c("effect_size", "r", "phi", "rho2", "estimand")
  vary_cols   <- design_cols[sapply(design_cols, function(v)
    length(unique(res[[v]])) > 1)]

  if (is.null(x_var)) x_var <- if (length(vary_cols) >= 1) vary_cols[1] else design_cols[1]
  if (is.null(group) && length(vary_cols) >= 2) group <- vary_cols[2]
  if (is.null(facet) && length(vary_cols) >= 3) facet  <- vary_cols[3]

  if (!x_var %in% names(res))
    stop("'x_var' must be a column in the result data frame.")

  if (!is.null(group)) res[[group]] <- factor(res[[group]])
  if (!is.null(facet)) res[[facet]] <- factor(res[[facet]])

  aes_args <- if (!is.null(group)) {
    ggplot2::aes(x = .data[[x_var]], y = .data[[val_col]],
                 color = .data[[group]], group = .data[[group]])
  } else {
    ggplot2::aes(x = .data[[x_var]], y = .data[[val_col]])
  }

  p <- ggplot2::ggplot(res, aes_args) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x     = .col_label(x_var),
      y     = .col_label(val_col),
      color = if (!is.null(group)) .col_label(group, legend = TRUE) else NULL
    ) +
    ggplot2::theme_bw()

  if (!is.null(facet))
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet)))

  p
}


# Null-coalescing operator (internal)
`%||%` <- function(a, b) if (!is.null(a)) a else b
