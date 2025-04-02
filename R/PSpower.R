#' Calculate sample size needed to achieve a prespecified power
#'
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @param tau the anticipated standardized treatment effect
#' @param sig.level the significance level, or the type-I error rate (default 0.05)
#' @param power the desired power to achieve (only specify for sample size calculation)
#' @param sample.size the total sample size (only specify for power calculation)
#' @param r the proportion of treated units
#' @param phi the overlap coefficient (usually between 0.8 and 1); use function plot_overlap(r, phi) for visual aid
#' @param rho_sq the squared correlation between propensity score and outcome; recommend treating as a sensitivity parameter: a grid of values between 0 and the R-squared statistic of predicting the outcomes with covariates.
#' @param test whether one-sided or two-sided test is considered
#' @param estimand the estimand (ATE, ATT, ATC or ATO), or a customized tilting function h(e(x))
#' @returns an object with the calculated sample size
#' @examples
#' PSpower(tau = 1/sqrt(20), sig.level = 0.05, power = 0.956, r = 0.5, phi = 0.99, rho_sq = 0.02)
#' @export
PSpower <- function(tau, sig.level = 0.05, power = NULL, sample.size = NULL, r, phi, rho_sq,
                    test = 'two-sided', estimand = 'ATE') {
  params <- solve_parameters(r, phi, rho_sq)
  if (typeof(estimand) == 'character') {
    if(!estimand %in% c('ATE', 'ATT', 'ATC', 'ATO'))
      stop('estimand should be one of ATE, ATT, ATC and ATO, or the tilting function.')
    if (estimand == 'ATE') {V <- params$V_ATE}
    else if (estimand == 'ATT') {V <- params$V_ATT}
    else if (estimand == 'ATC') {V <- params$V_ATC}
    else if (estimand == 'ATO') {V <- params$V_ATO}
  }
  else if (typeof(estimand) == 'closure') {
    V <- params$calc_V(estimand)
  }
  else {
    stop('estimand should be one of ATE, ATT, ATC and ATO, or the tilting function.')
  }

  if (!is.null(power)) {
    if (test == 'two-sided') {
      coef <- (qnorm(1 - sig.level / 2) + qnorm(power))^2
    }
    else if (test == 'one-sided') {
      coef <- (qnorm(1 - sig.level) + qnorm(power))^2
    }
    else {
      stop('test should be either one-sided or two-sided.')
    }
    sample_size <- coef * V / tau^2
    power <- power
    ss_calculation <- TRUE
  }
  else {
    if (test == 'two-sided') {
      power <- pnorm(sqrt(sample.size / V) * tau - qnorm(1 - sig.level / 2))
    }
    else if (test == 'one-sided') {
      power <- pnorm(sqrt(sample.size / V) * tau - qnorm(1 - sig.level))
    }
    else {
      stop('test should be either one-sided or two-sided.')
    }
    sample_size <- sample.size
    power <- power
    ss_calculation <- FALSE
  }

  return (structure(
    list(
      sample_size = sample_size,
      tau = tau, sig.level = sig.level, power = power, test = test, estimand = estimand,
      summaries = c(r = r, phi = phi, rho_sq = rho_sq),
      params = params,
      ss_calculation = ss_calculation
    ),
    class = 'PSpower'
  ))
}

#' Prints PSpower object
#' @param x PSpower object
#' @param ... ignored
#' @returns no return value; called for side effect to output a string
#' @export
print.PSpower <- function(x, ...) {
  if (x$ss_calculation) {
    cat('Estimated sample size to achieve', round(x$power, 3), 'power:', round(x$sample_size), '\n')
  }
  else {
    cat('Calculated power at sample size', round(x$sample_size), 'is', round(x$power, 3), '\n')
  }
}

#' Plots PSpower object
#' @param x PSpower object
#' @param power a range of powers to plot the power curve
#' @param ... ignored
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_line
#' @examples
#' obj <- PSpower(tau = 1/sqrt(20), sig.level = 0.05, power = 0.956,
#'                r = 0.5, phi = 0.99, rho_sq = 0.02)
#' plot(obj)
#' @returns an object (class ggplot) containing a figure
#' @export
plot.PSpower <- function(x, power = seq(0.6, 0.99, length.out = 100), ...) {
  size <- sapply(power,
         function(p)
           PSpower(x$tau, x$sig.level, p, x$summaries['r'], x$summaries['phi'],
                   x$summaries['rho_sq'], sample.size = NULL, x$test, x$estimand)$sample_size
  )
  ggplot2::ggplot(data.frame(power = power, size = size), ggplot2::aes(size, power)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(data = data.frame(power = x$power, size = x$sample_size))
}
