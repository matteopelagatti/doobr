#' Energy score of a predictive sample
#'
#' The multivariate energy score of an ensemble against a realised path
#' (Gneiting and Raftery, 2007, eq. 22).  It is a strictly proper scoring rule,
#' sensitive to both location and spread.
#'
#' @param y_obs Numeric vector of length \eqn{h}: the realised path.
#' @param dat Numeric matrix (\eqn{h \times m}) of simulated paths, one per
#'   column.  A plain vector is treated as a single path.
#'
#' @return A non-negative scalar; lower is better.
#' @references Gneiting, T. and Raftery, A. E. (2007). Strictly proper scoring
#'   rules, prediction, and estimation. *JASA* 102, 359--378.
#' @seealso [vs_score()]
#' @export
#' @examples
#' es_score(rnorm(5), matrix(rnorm(50), nrow = 5))
es_score <- function(y_obs, dat) {
  y_obs <- as.numeric(y_obs)
  dat   <- if (is.matrix(dat)) dat else matrix(dat, ncol = 1L)
  m  <- ncol(dat)
  t1 <- mean(sqrt(colSums((dat - y_obs)^2)))
  t2 <- sum(stats::dist(t(dat))) / m^2
  as.numeric(t1 - t2)
}

#' Variogram score of a predictive sample
#'
#' The variogram score of order `p` (Scheuerer and Hamill, 2015).  It is
#' sensitive to the cross-horizon dependence structure but not to the marginal
#' distributions, so it complements [es_score()].
#'
#' @inheritParams es_score
#' @param p Order of the variogram score.
#'
#' @return A non-negative scalar; lower is better.  Zero when \eqn{h = 1}.
#' @references Scheuerer, M. and Hamill, T. M. (2015). Variogram-based proper
#'   scoring rules for probabilistic forecasts of multivariate quantities.
#'   *Monthly Weather Review* 143, 1321--1334.
#' @seealso [es_score()]
#' @export
#' @examples
#' vs_score(rnorm(5), matrix(rnorm(50), nrow = 5))
vs_score <- function(y_obs, dat, p = 0.5) {
  y_obs <- as.numeric(y_obs)
  dat   <- if (is.matrix(dat)) dat else matrix(dat, ncol = 1L)
  h <- nrow(dat)
  if (h < 2L) return(0)
  obs_pw <- abs(outer(y_obs, y_obs, "-"))^p
  sim_pw <- matrix(0, h, h)
  for (i in seq_len(h))
    for (j in seq_len(h))
      sim_pw[i, j] <- mean(abs(dat[i, ] - dat[j, ])^p)
  sum((obs_pw - sim_pw)^2)
}
