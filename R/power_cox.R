#' Sample size and power for PS-weighted marginal Cox model
#'
#' Computes the required sample size or the achieved power for the
#' propensity-score-weighted partial likelihood estimator of the marginal
#' hazard ratio in a Cox proportional hazards model.
#'
#' @description
#' The required sample size is
#' \deqn{N = V\,(z_{1-\alpha} + z_{\beta})^2 \,/\, \tau_0^2,}
#' where \eqn{\tau_0 = \log(\text{HR})} is the target log hazard ratio and
#' \eqn{V} is the asymptotic variance of the estimator.
#'
#' **Randomized trial — robust sandwich variance** (\code{method = "robust"}):
#' \deqn{V_{RCT} = \frac{(\lambda_1 + \lambda_0)^2
#'   \bigl[r\lambda_0^2 d_1 + (1-r)\lambda_1^2 d_0\bigr]}{d^2},}
#' where \eqn{\lambda_1 = \sqrt{r/(1-r)}\,e^{\tau_0/2}},
#' \eqn{\lambda_0 = 1/\lambda_1}, and \eqn{d = r\,d_1 + (1-r)\,d_0}.
#'
#' **Randomized trial — Schoenfeld formula** (\code{method = "schoenfeld"}):
#' \deqn{V_{Sch} = \frac{1}{r(1-r)\,d}.}
#' Note: the Schoenfeld formula is derived under a null effect and may
#' underestimate or overestimate the required sample size at non-null effects.
#'
#' **Observational study — inverse probability weights (ATE),
#' robust sandwich variance** (\code{study_type = "obs"}, \code{estimand = "ATE"}):
#' \deqn{V_{obs} = \frac{(\lambda_1+\lambda_0)^2}{d^2}
#'   \left[r^2\lambda_0^2 d_1\,\frac{a+b-1}{a-1}
#'        + (1-r)^2\lambda_1^2 d_0\,\frac{a+b-1}{b-1}\right],}
#' where \eqn{a, b > 1} are Beta distribution parameters determined by
#' \eqn{(r, \phi)}.  Requires \eqn{\min(a, b) > 1}.
#'
#' **Observational study — overlap weights (ATO) or
#' treated population weights (ATT)**
#' (\code{study_type = "obs"}, \code{estimand} in \code{"ATO"}, \code{"ATT"}):
#' \deqn{N = \kappa_{DE} \times N_{RCT},}
#' where \eqn{N_{RCT}} uses \eqn{V_{RCT}} above and \eqn{\kappa_{DE}} is a
#' design effect estimated by Monte Carlo simulation from the Beta
#' approximation of propensity scores.
#'
#' @param effect_size Log hazard ratio \eqn{\tau_0 = \log(\text{HR})}.
#'   Negative values indicate benefit (lower hazard in group 1).
#'   Scalar or vector.
#' @param r Treatment proportion \eqn{r = \Pr(Z = 1)}, in \eqn{(0, 1)}.
#'   Scalar or vector.
#' @param d1 Event rate in group 1 (treated), in \eqn{(0, 1]}.
#'   Scalar or vector.
#' @param d0 Event rate in group 0 (control), in \eqn{(0, 1]}.
#'   If \code{NULL} (default), set equal to \code{d1}.  Scalar or vector.
#' @param phi Overlap coefficient \eqn{\phi \in (0, 1)}.  Required when
#'   \code{study_type = "obs"}; ignored for \code{"rct"}.  Rule of thumb:
#'   \eqn{< 0.80} very poor, \eqn{[0.80, 0.90)} poor,
#'   \eqn{[0.90, 0.95)} moderate, \eqn{\ge 0.95} good.
#'   Scalar or vector.
#' @param study_type \code{"obs"} (observational study, default) or
#'   \code{"rct"} (randomized trial).
#' @param estimand Target estimand.  \code{"ATE"} (average treatment effect,
#'   uses inverse probability weights), \code{"ATO"} (overlap population,
#'   uses overlap weights), or \code{"ATT"} (group 1 population, uses
#'   weights for the treated). Ignored when \code{study_type = "rct"}.
#'   Scalar or character vector.
#' @param method Variance approximation method.  \code{"robust"} (default)
#'   for the robust sandwich variance; \code{"schoenfeld"} for the classical
#'   Schoenfeld formula. \code{"schoenfeld"} is only available when
#'   \code{study_type = "rct"}. Scalar or character vector.
#' @param sig_level Significance level \eqn{\alpha}. Default \code{0.05}.
#' @param power Target power \eqn{\beta}. Provide for sample size calculation;
#'   mutually exclusive with \code{sample_size}.
#' @param sample_size Total sample size \eqn{N}. Provide for power
#'   calculation; mutually exclusive with \code{power}.
#' @param test \code{"one-sided"} (default) or \code{"two-sided"}.
#' @param n_mc Number of Monte Carlo samples used to estimate the design
#'   effect for \code{estimand} in \code{"ATO"}, \code{"ATT"}.
#'   Default \code{1e6}.
#'
#' @return An object of class \code{"power_cox"}, a list containing:
#'   \describe{
#'     \item{\code{call}}{The matched call.}
#'     \item{\code{calculation}}{\code{"sample_size"} or \code{"power"}.}
#'     \item{\code{result}}{A data frame with one row per scenario and columns
#'       for every design parameter plus the computed \code{sample_size} or
#'       \code{power}.}
#'     \item{\code{settings}}{A list with \code{sig_level}, \code{power},
#'       \code{sample_size}, and \code{test}.}
#'     \item{\code{n_scenarios}}{Number of rows in \code{result}.}
#'     \item{\code{d0_set_equal}}{Logical; \code{TRUE} when \code{d0} was
#'       not specified and set equal to \code{d1}.}
#'   }
#'
#' @references
#' Chengxin Yang, Bo Liu, and Fan Li. Sample size and power calculations for
#' causal inference with time-to-event outcomes.
#' \emph{arXiv preprint arXiv:2605.10088} (2026).
#'
#' @seealso \code{\link{power_ps}}, \code{\link{overlap_coef}}
#'
#' @examples
#' # RCT sample size, robust variance
#' power_cox(effect_size = log(0.6), r = 0.5, d1 = 0.8, study_type = "rct",
#'           power = 0.8)
#'
#' # Observational study, ATE
#' power_cox(effect_size = log(0.6), r = 0.5, d1 = 0.8, phi = 0.9,
#'           study_type = "obs", estimand = "ATE", power = 0.8)
#'
#' # Sensitivity over phi and estimand
#' power_cox(effect_size = log(0.6), r = 0.5, d1 = 0.8,
#'           phi = c(0.9, 0.95), estimand = c("ATE", "ATO"),
#'           power = 0.8)
#'
#' @export
power_cox <- function(effect_size,
                      r,
                      d1,
                      d0          = NULL,
                      phi         = NULL,
                      study_type  = "obs",
                      estimand    = "ATE",
                      method      = "robust",
                      sig_level   = 0.05,
                      power       = NULL,
                      sample_size = NULL,
                      test        = "one-sided",
                      n_mc        = 1e6) {

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
  if (!test %in% c("one-sided", "two-sided"))
    stop("'test' must be \"one-sided\" or \"two-sided\".")
  if (!all(study_type %in% c("obs", "rct")))
    stop("'study_type' must be \"obs\" or \"rct\".")
  if (!all(estimand %in% c("ATE", "ATO", "ATT")))
    stop("'estimand' must be one of \"ATE\", \"ATO\", \"ATT\".")
  if (!all(method %in% c("robust", "schoenfeld")))
    stop("'method' must be \"robust\" or \"schoenfeld\".")
  if (any(r <= 0 | r >= 1))
    stop("All values of 'r' must be in (0, 1).")
  if (any(d1 <= 0 | d1 > 1))
    stop("All values of 'd1' must be in (0, 1].")

  # Handle d0 default
  d0_set_equal <- is.null(d0)
  if (d0_set_equal) d0 <- d1
  if (any(d0 <= 0 | d0 > 1))
    stop("All values of 'd0' must be in (0, 1].")

  # phi required for obs
  obs_requested <- any(study_type == "obs")
  if (obs_requested && is.null(phi))
    stop("'phi' is required when study_type = \"obs\".")
  if (!is.null(phi) && any(phi <= 0 | phi >= 1))
    stop("All values of 'phi' must be in (0, 1).")

  # schoenfeld only for RCT
  if ("schoenfeld" %in% method && any(study_type == "rct") == FALSE)
    stop("method = \"schoenfeld\" is only available when study_type = \"rct\".")

  calculation <- if (!is.null(power)) "sample_size" else "power"

  # ---- Build scenario grid ----
  # phi participates in grid only for obs; for rct scenarios phi is NA
  phi_vec <- if (!is.null(phi)) phi else NA_real_

  grid <- expand.grid(
    effect_size = effect_size,
    r           = r,
    d1          = d1,
    d0          = d0,
    phi         = phi_vec,
    study_type  = study_type,
    estimand    = estimand,
    method      = method,
    stringsAsFactors = FALSE
  )

  # Validate schoenfeld + obs combinations row-wise and drop/error
  bad_rows <- grid$method == "schoenfeld" & grid$study_type == "obs"
  if (any(bad_rows))
    stop("method = \"schoenfeld\" is not applicable for observational studies. ",
         "Found ", sum(bad_rows), " such scenario(s) in the grid.")

  # ---- Compute V and result for each scenario ----
  result_col <- vapply(seq_len(nrow(grid)), function(i) {
    es  <- grid$effect_size[i]
    ri  <- grid$r[i]
    d1i <- grid$d1[i]
    d0i <- grid$d0[i]
    pi  <- grid$phi[i]
    st  <- grid$study_type[i]
    est <- grid$estimand[i]
    mth <- grid$method[i]

    if (st == "rct") {
      V <- if (mth == "robust") {
        .V_RCT(es, ri, d1i, d0i)
      } else {
        .V_Schoenfeld(ri, d1i, d0i)
      }
    } else {
      # Observational study
      ab <- .solve_ab(ri, pi)
      a  <- ab$a;  b <- ab$b

      if (est == "ATE") {
        if (a <= 1 || b <= 1)
          stop(sprintf(
            "Beta parameters a=%.3f, b=%.3f are not both > 1 (r=%.2f, phi=%.2f). ",
            a, b, ri, pi,
            "IPW variance (ATE) requires phi large enough that min(a,b) > 1. ",
            "Consider using estimand = \"ATO\" or increasing phi."))
        V <- .V_IPW(es, ri, a, b, d1i, d0i)
      } else {
        # ATO or ATT: design-effect approach
        scheme <- if (est == "ATO") "ow" else "treated"
        kappa  <- .kappa_Kish(ri, pi, scheme, n_mc = n_mc)
        V      <- kappa * .V_RCT(es, ri, d1i, d0i)
      }
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
      call         = match.call(),
      calculation  = calculation,
      result       = grid,
      settings     = list(
        sig_level   = sig_level,
        power       = power,
        sample_size = sample_size,
        test        = test
      ),
      n_scenarios  = nrow(grid),
      d0_set_equal = d0_set_equal
    ),
    class = "power_cox"
  )
}


# --------------------------------------------------------------------------
# S3 methods
# --------------------------------------------------------------------------

#' @export
print.power_cox <- function(x, ...) {
  cat("PSpower: power_cox\n")
  cat(strrep("-", 55), "\n")

  res    <- x$result
  s      <- x$settings
  single <- x$n_scenarios == 1L

  if (single) {
    st  <- res$study_type[1]
    mth <- res$method[1]
    est <- res$estimand[1]

    cat(sprintf("%-30s %s\n", "Study type:",
                if (st == "rct") "Randomized trial" else "Observational study"))
    cat(sprintf("%-30s %s\n", "Method:",
                if (mth == "robust") "Robust sandwich variance" else "Schoenfeld formula"))
    if (st == "obs") {
      cat(sprintf("%-30s %s\n", "Estimand:", est))
      cat(sprintf("%-30s %s\n", "Weights:", .weight_label(est)))
    }
    cat(sprintf("%-30s %.3f  [hazard ratio = %.3f]\n",
                "Effect size (log hazard ratio):",
                res$effect_size[1], exp(res$effect_size[1])))
    cat(sprintf("%-30s %.3f\n", "Treatment proportion:", res$r[1]))
    cat(sprintf("%-30s %.3f\n", "Event rate (group 1):", res$d1[1]))
    d0_tag <- if (x$d0_set_equal) "  [= group 1]" else ""
    cat(sprintf("%-30s %.3f%s\n", "Event rate (group 0):", res$d0[1], d0_tag))
    if (st == "obs")
      cat(sprintf("%-30s %.3f\n", "Overlap (phi):", res$phi[1]))

    cat(strrep("-", 55), "\n")
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
    design_cols <- c("effect_size","r","d1","d0","phi","study_type","estimand","method")
    vary_cols   <- design_cols[sapply(design_cols, function(v)
      v %in% names(res) && length(unique(res[[v]])) > 1)]
    fixed_cols  <- setdiff(design_cols[design_cols %in% names(res)], vary_cols)

    cat(sprintf("Scenarios: %d\n", x$n_scenarios))
    if (length(vary_cols))
      cat("Varying: ", paste(vary_cols, collapse = ", "), "\n")
    if (length(fixed_cols)) {
      fixed_str <- paste(sapply(fixed_cols, function(v) {
        val    <- unique(res[[v]])
        d0_tag <- if (v == "d0" && x$d0_set_equal) " [= group 1]" else ""
        paste0(v, " = ", format(val, digits = 3), d0_tag)
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
    cat(strrep("-", 55), "\n")

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
summary.power_cox <- function(object, ...) {
  cat("PSpower: power_cox\n")
  cat(strrep("-", 60), "\n")

  res     <- object$result
  s       <- object$settings
  single  <- object$n_scenarios == 1L
  val_col <- object$calculation

  design_cols <- c("effect_size","r","d1","d0","phi","study_type","estimand","method")
  vary_cols   <- design_cols[sapply(design_cols, function(v)
    v %in% names(res) && length(unique(res[[v]])) > 1)]
  fixed_cols  <- setdiff(design_cols[design_cols %in% names(res)], vary_cols)

  if (!single) cat(sprintf("Scenarios: %d\n", object$n_scenarios))
  if (length(vary_cols))
    cat("Varying: ", paste(vary_cols, collapse = ", "), "\n")
  if (length(fixed_cols)) {
    fixed_str <- paste(sapply(fixed_cols, function(v) {
      val    <- unique(res[[v]])
      d0_tag <- if (v == "d0" && object$d0_set_equal) " [= group 1]" else ""
      paste0(v, " = ", format(val, digits = 3), d0_tag)
    }), collapse = ",  ")
    cat("Fixed:   ", fixed_str, "\n")
  }
  if (!is.null(s$power))
    cat(sprintf("Target power: %.2f  |  sig. level: %.3f  |  %s test\n",
                s$power, s$sig_level, s$test))
  else
    cat(sprintf("Sample size: %d  |  sig. level: %.3f  |  %s test\n",
                as.integer(s$sample_size), s$sig_level, s$test))
  cat(strrep("-", 60), "\n")

  if (single) {
    st  <- res$study_type[1]; mth <- res$method[1]; est <- res$estimand[1]
    cat(sprintf("%-30s %s\n", "Study type:",
                if (st == "rct") "Randomized trial" else "Observational study"))
    cat(sprintf("%-30s %s\n", "Method:",
                if (mth == "robust") "Robust sandwich variance" else "Schoenfeld formula"))
    if (st == "obs") {
      cat(sprintf("%-30s %s\n", "Estimand:", est))
      cat(sprintf("%-30s %s\n", "Weights:", .weight_label(est)))
    }
    cat(sprintf("%-30s %.3f  [HR = %.3f]\n",
                "Effect size (log hazard ratio):",
                res$effect_size[1], exp(res$effect_size[1])))
    cat(sprintf("%-30s %.3f\n", "Treatment proportion:", res$r[1]))
    cat(sprintf("%-30s %.3f\n", "Event rate (group 1):", res$d1[1]))
    d0_tag <- if (object$d0_set_equal) "  [= group 1]" else ""
    cat(sprintf("%-30s %.3f%s\n", "Event rate (group 0):", res$d0[1], d0_tag))
    if (st == "obs") cat(sprintf("%-30s %.3f\n", "Overlap (phi):", res$phi[1]))
    cat(strrep("-", 60), "\n")
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
plot.power_cox <- function(x, x_var = NULL, group = NULL, facet = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for plotting.")

  res     <- x$result
  val_col <- x$calculation

  design_cols <- c("effect_size","r","d1","d0","phi","study_type","estimand","method")
  vary_cols   <- design_cols[sapply(design_cols, function(v)
    v %in% names(res) && length(unique(res[[v]])) > 1)]

  if (is.null(x_var)) x_var <- if (length(vary_cols) >= 1) vary_cols[1] else "r"
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
