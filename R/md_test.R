#' Martingale-difference hypothesis test
#'
#' Tests \eqn{E[y_t \mid \mathcal{F}_{t-1}] = 0} by regressing the series on
#' B-spline transforms of the lagged values mapped to the empirical copula
#' scale.  This is stronger than a test of zero autocorrelation: it has power
#' against nonlinear predictable structure, which a Ljung-Box test cannot see,
#' and the copula transform makes it invariant to monotone transformations of
#' the series.
#'
#' @param y Numeric vector; the series to test.  In the scenario-generation
#'   pipeline this is the out-of-fold residual series (step 1) or the squared
#'   residuals net of the fitted conditional variance (step 2).
#' @param lags Integer vector of lag indices.  `NULL` uses
#'   \eqn{\lfloor n^{0.3}\rfloor}.
#' @param zero_mean If `TRUE` tests the strict hypothesis; if `FALSE` tests
#'   whether the conditional mean is constant but possibly non-zero.
#' @param df Degrees of freedom per B-spline basis.
#' @param q Number of random projections when `type = "random"`.
#' @param type Design-matrix construction: `"additive"` (independent splines per
#'   lag), `"tensor"` (additive plus all pairwise products) or `"random"`
#'   (projection pursuit on `q` random directions, which captures multi-lag
#'   interactions with a bounded number of parameters).
#' @param test Test statistic: `"F"`, `"Chisq"`, `"HC3"` or `"WildBoot"`.  The
#'   wild bootstrap is recommended: it is valid in finite samples under unknown
#'   conditional heteroskedasticity, which is the usual situation here.
#' @param B Bootstrap replications for `test = "WildBoot"`.
#'
#' @return The p-value, or `NA` if the design is singular.  Small values are
#'   evidence against the martingale-difference hypothesis.
#'
#' @seealso [doob_fit()], which applies this to both estimation steps.
#' @export
#' @examples
#' set.seed(1)
#' md_test(rnorm(400), type = "additive", test = "F")        # should not reject
#' x <- as.numeric(arima.sim(list(ar = 0.5), 400))
#' md_test(x, type = "additive", test = "F")                 # should reject
md_test <- function(y, lags = NULL, zero_mean = TRUE, df = 3, q = 3,
                    type = c("additive", "tensor", "random"),
                    test = c("F", "Chisq", "HC3", "WildBoot"),
                    B = 499) {
  type <- match.arg(type)
  test <- match.arg(test)
  y    <- as.numeric(y)
  n    <- length(y)
  if (is.null(lags)) lags <- 1:max(1, floor(n^0.3))
  lags <- as.integer(lags)

  M     <- max(lags)
  n_obs <- n - M
  if (n_obs <= length(lags) + 2L) return(NA_real_)
  y_vec <- y[(M + 1):n]

  Z <- matrix(NA_real_, nrow = n_obs, ncol = length(lags))
  for (i in seq_along(lags)) Z[, i] <- y[(M + 1 - lags[i]):(n - lags[i])]

  if (type %in% c("additive", "tensor")) {
    spline_list <- vector("list", length(lags))
    for (i in seq_along(lags)) {
      u_lag <- rank(Z[, i]) / (n_obs + 1)
      spline_list[[i]] <- splines::bs(u_lag, df = df)
    }
    XX_base <- do.call(cbind, spline_list)
    if (type == "tensor" && ncol(XX_base) > 1) {
      idx <- utils::combn(ncol(XX_base), 2)
      XX  <- cbind(1, XX_base, XX_base[, idx[1, ]] * XX_base[, idx[2, ]])
    } else {
      XX <- cbind(1, XX_base)
    }
  } else {
    U_Z  <- apply(Z, 2, rank) / (n_obs + 1)
    Beta <- matrix(stats::rnorm(length(lags) * q), nrow = length(lags), ncol = q)
    Proj <- U_Z %*% Beta
    spline_list <- vector("list", q)
    for (j in seq_len(q)) {
      u_proj <- rank(Proj[, j]) / (n_obs + 1)
      spline_list[[j]] <- splines::bs(u_proj, df = df)
    }
    XX <- cbind(1, do.call(cbind, spline_list))
  }

  if (zero_mean) {
    rss0 <- sum(y_vec^2); df0 <- n_obs; e0 <- y_vec
  } else {
    rss0 <- sum((y_vec - mean(y_vec))^2); df0 <- n_obs - 1; e0 <- y_vec - mean(y_vec)
  }

  fit1   <- stats::lm.fit(XX, y_vec)
  rss1   <- sum(fit1$residuals^2)
  df1    <- fit1$df.residual
  df_num <- df0 - df1
  if (df_num <= 0 || rss1 <= 1e-10) return(NA_real_)

  if (test == "F") {
    f_stat <- ((rss0 - rss1) / df_num) / (rss1 / df1)
    return(unname(stats::pf(f_stat, df_num, df1, lower.tail = FALSE)))
  }
  if (test == "Chisq") {
    lr <- n_obs * log(rss0 / rss1)
    return(unname(stats::pchisq(lr, df_num, lower.tail = FALSE)))
  }
  if (test == "HC3") {
    keep   <- !is.na(fit1$coefficients)
    XX_eff <- XX[, keep, drop = FALSE]
    p_eff  <- sum(keep)
    hh     <- pmin(rowSums(qr.Q(fit1$qr)^2), 0.9999)
    omega  <- fit1$residuals / (1 - hh)
    meat   <- crossprod(XX_eff * omega)
    bread  <- chol2inv(fit1$qr$qr[1:p_eff, 1:p_eff, drop = FALSE])
    V      <- bread %*% meat %*% bread
    if (zero_mean) {
      Vs <- V; bs <- fit1$coefficients[keep]
    } else {
      Vs <- V[-1, -1, drop = FALSE]; bs <- fit1$coefficients[keep][-1]
    }
    wald <- as.numeric(t(bs) %*% solve(Vs, bs))
    return(unname(stats::pf(wald / df_num, df_num, df1, lower.tail = FALSE)))
  }
  ## WildBoot
  f_obs  <- ((rss0 - rss1) / df_num) / (rss1 / df1)
  boot_f <- numeric(B)
  for (b in seq_len(B)) {
    v <- sample(c(-1, 1), n_obs, replace = TRUE)
    if (zero_mean) {
      y_star <- e0 * v; rss0_s <- sum(y_star^2)
    } else {
      y_star <- mean(y_vec) + e0 * v; rss0_s <- sum((y_star - mean(y_star))^2)
    }
    rss1_s    <- sum(qr.resid(fit1$qr, y_star)^2)
    boot_f[b] <- ((rss0_s - rss1_s) / df_num) / (rss1_s / df1)
  }
  unname(mean(boot_f >= f_obs))
}
