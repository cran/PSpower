#' Plot density of propensity scores given treatment probability and overlap coefficient
#'
#' @param r treatment probability
#' @param phi overlap coefficient
#' @returns a ggplot of the density of propensity scores in two treatment arms
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 theme_bw
#' @examples
#' plot_overlap(0.6, 0.9)
#' @export
plot_overlap <- function(r, phi) {
  overlap_object <- solve_overlap(r, phi)
  mu_e <- overlap_object$mu
  sigma_e <- overlap_object$sigma

  ps <- seq(0, 1, length.out = 1000)
  f0 <- ps / r * dlogitnorm(ps, mu_e, sigma_e)
  f1 <- (1 - ps) / (1 - r) * dlogitnorm(ps, mu_e, sigma_e)
  df <- data.frame(
    ps = ps, f0 = f0, f1 = f1
  )
  ggplot(df) + geom_line(aes(ps, f0, color = 'control')) +
    geom_line(aes(ps, f1, color = 'treatment')) +
    xlab('propensity score') + ylab('density') +
    theme_bw()
}
