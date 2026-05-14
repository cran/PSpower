#' Bhattacharyya overlap coefficient for propensity score distributions
#'
#' Computes the Bhattacharyya overlap coefficient \eqn{\phi}, a scalar
#' measure of propensity score overlap between the treatment and control
#' groups.  Values close to 1 indicate near-complete overlap (little
#' confounding); values well below 1 indicate poor overlap.
#'
#' @description
#' Two calling conventions are provided:
#'
#' \describe{
#'   \item{From propensity scores (empirical):}{
#'     Supply \code{ps} and \code{Z}.  The empirical formula is
#'     \deqn{\hat\phi = \frac{\mathrm{E}[\sqrt{e(1-e)}]}{\sqrt{r(1-r)}},}
#'     where \eqn{r = \mathrm{E}[Z]}.  This is the sample mean of
#'     \eqn{\sqrt{e_i(1-e_i)}} divided by \eqn{\sqrt{\hat r(1-\hat r)}}.
#'   }
#'   \item{From Beta parameters (analytical):}{
#'     Supply \code{a} and \code{b}.  Under the
#'     \eqn{\mathrm{Beta}(a,b)} approximation,
#'     \deqn{\phi = \exp\!\Bigl[
#'       \log\Gamma(a+\tfrac12) - \tfrac12\log a - \log\Gamma(a)
#'       + \log\Gamma(b+\tfrac12) - \tfrac12\log b - \log\Gamma(b)
#'     \Bigr].}
#'   }
#' }
#'
#' @param ps Numeric vector of estimated propensity scores
#'   \eqn{e_i = \Pr(Z_i = 1 \mid X_i)}, all in \eqn{(0, 1)}.
#'   Required when \code{a} and \code{b} are not supplied.
#' @param Z Integer or numeric vector of treatment indicators
#'   (\eqn{Z_i \in \{0, 1\}}), the same length as \code{ps}.
#'   Required when \code{ps} is supplied.
#' @param a Shape parameter \eqn{a > 0} of the Beta distribution.
#'   Supply together with \code{b} to use the analytical formula.
#' @param b Shape parameter \eqn{b > 0} of the Beta distribution.
#'   Supply together with \code{a} to use the analytical formula.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{phi}}{The overlap coefficient \eqn{\hat\phi}.}
#'     \item{\code{r}}{Treatment proportion: \code{mean(Z)} (empirical) or
#'       \code{a / (a + b)} (analytical).}
#'   }
#'
#' @references
#' Chengxin Yang, Bo Liu, and Fan Li. Sample size and power calculations for
#' causal inference with time-to-event outcomes.
#' \emph{arXiv preprint arXiv:2605.10088} (2026).
#'
#' Bo Liu, Chengxin Yang, and Fan Li. Sample size and power calculations for
#' causal inference with continuous and binary outcomes.
#' \emph{Annals of Statistics} (2026).
#'
#' @seealso \code{\link{power_ps}}, \code{\link{power_cox}}
#'
#' @examples
#' # From propensity scores
#' set.seed(1)
#' n  <- 500
#' X  <- rnorm(n)
#' ps <- plogis(0.5 * X)
#' Z  <- rbinom(n, 1, ps)
#' overlap_coef(ps = ps, Z = Z)
#'
#' # From Beta parameters
#' overlap_coef(a = 2, b = 3)
#'
#' @export
overlap_coef <- function(ps = NULL, Z = NULL, a = NULL, b = NULL) {

  empirical  <- !is.null(ps) || !is.null(Z)
  analytical <- !is.null(a)  || !is.null(b)

  if (empirical && analytical)
    stop("Supply either (ps, Z) or (a, b), not both.")
  if (!empirical && !analytical)
    stop("Supply either (ps, Z) for the empirical formula or (a, b) for the analytical formula.")

  if (empirical) {
    if (is.null(ps) || is.null(Z))
      stop("Both 'ps' and 'Z' are required for the empirical formula.")
    if (length(ps) != length(Z))
      stop("'ps' and 'Z' must have the same length.")
    if (any(ps <= 0 | ps >= 1))
      stop("All values of 'ps' must be in (0, 1).")
    if (!all(Z %in% c(0, 1)))
      stop("'Z' must contain only 0s and 1s.")

    r   <- mean(Z)
    phi <- mean(sqrt(ps * (1 - ps))) / sqrt(r * (1 - r))
    return(list(phi = phi, r = r))
  }

  # Analytical
  if (is.null(a) || is.null(b))
    stop("Both 'a' and 'b' are required for the analytical formula.")
  if (length(a) != 1 || length(b) != 1)
    stop("'a' and 'b' must each be a single positive number.")
  if (a <= 0 || b <= 0)
    stop("'a' and 'b' must be positive.")

  phi <- .phi_from_ab(a, b)
  r   <- a / (a + b)
  list(phi = phi, r = r)
}
