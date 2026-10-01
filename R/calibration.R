#' Calibration tests for simulated sample paths
#'
#' Three complementary tests of whether the realised path is distinguishable
#' from the simulated ensemble, applied to the per-fold evaluation windows of a
#' [doob_fit()] object and combined across folds by Fisher's method.
#'
#' @param object A `"doob_fit"` object fitted with `eval_folds = TRUE`, or a
#'   list of matrices in the layout of `object$eval$window_sims`: one matrix per
#'   fold, column 1 the realised path and the rest the simulated ones.
#' @param vs_p Order of the variogram score used by the permutation test.
#' @param alpha Significance level for the `reject` flags.
#'
#' @return An object of class `"doob_calibration"` with components `rank`,
#'   `perm_es`, `perm_vs` and `band`, each carrying per-fold p-values where
#'   applicable and a combined p-value.
#'
#' @details
#' **Rank histogram.** The rank of each realised value among the simulated ones
#' is uniform under calibration; a chi-squared goodness-of-fit test assesses it.
#' The shape is interpretable: U-shaped means too little simulated spread,
#' hump-shaped too much, skewed means a biased conditional mean.  The test is
#' marginal and blind to temporal dependence.  Note that it is **invariant to a
#' monotone transformation** of the series, whereas the other two are not.
#'
#' **Permutation test on proper scores.** The score of the realised path is
#' compared with the scores obtained by treating each simulated path in turn as
#' the realisation.  Exact in finite samples by exchangeability under the null.
#'
#' **Band depth.** A path-level test: the centrality of the realised path within
#' the ensemble is compared with that of the simulated paths.  A p-value near
#' one means the realised path is *more* central than the simulated ones, i.e.
#' the ensemble is too dispersed at the path level.
#'
#' The pooled rank test treats ranks at adjacent steps within a fold as
#' independent when they are positively dependent, so it is mildly liberal; read
#' its p-values against a baseline obtained under a correctly specified model
#' rather than against the nominal level.
#'
#' @seealso [doob_fit()], [fisher_combine()]
#' @export
#' @examples
#' set.seed(1)
#' y <- as.numeric(arima.sim(list(ar = 0.6), 400)) + 20
#' fit <- doob_fit(y, horizon = 10, mean_model_fn = lm_model_fn, mlags = 1,
#'                 predict_mean_fn = lm_predict_fn, nsim_eval = 50,
#'                 run_md_test = FALSE)
#' calibration_tests(fit)
calibration_tests <- function(object, vs_p = 0.5, alpha = 0.05) {
  ws <- if (inherits(object, "doob_fit")) {
    if (is.null(object$eval))
      stop("this fit has no evaluation windows; refit with eval_folds = TRUE",
           call. = FALSE)
    object$eval$window_sims
  } else object
  if (!is.list(ws) || !length(ws))
    stop("expected a doob_fit object or a list of evaluation windows", call. = FALSE)

  K <- length(ws)
  m <- ncol(ws[[1L]]) - 1L

  ## ---- rank histogram ----------------------------------------------------
  ranks <- unlist(lapply(ws, function(W) {
    a <- W[, 1L]; S <- W[, -1L, drop = FALSE]
    vapply(seq_len(nrow(S)), function(t) sum(S[t, ] < a[t]) + 1L, integer(1L))
  }))
  nb     <- min(m + 1L, 10L)
  breaks <- seq(0.5, m + 1.5, length.out = nb + 1L)
  obs    <- as.integer(table(cut(ranks, breaks = breaks, include.lowest = TRUE)))
  chi2   <- suppressWarnings(stats::chisq.test(obs, p = rep(1 / nb, nb)))

  ## ---- permutation tests on proper scores --------------------------------
  perm <- function(score) {
    p <- vapply(ws, function(W) {
      a <- W[, 1L]; S <- W[, -1L, drop = FALSE]
      f <- if (score == "ES") function(o, d) es_score(o, d)
           else function(o, d) vs_score(o, d, p = vs_p)
      s_obs <- f(a, S)
      s_sim <- vapply(seq_len(ncol(S)), function(j)
        f(S[, j], cbind(S[, -j, drop = FALSE], a)), numeric(1L))
      (sum(s_sim >= s_obs) + 1) / (ncol(S) + 1)
    }, numeric(1L))
    fc <- fisher_combine(p)
    list(pvalues = p, fisher_stat = fc$stat, fisher_df = fc$df,
         fisher_pval = fc$pval, reject = fc$pval < alpha)
  }

  ## ---- band depth --------------------------------------------------------
  p_band <- vapply(ws, function(W) {
    a <- W[, 1L]; S <- W[, -1L, drop = FALSE]
    ctr   <- rowMeans(S)
    d_sim <- sqrt(colSums((S - ctr)^2))
    d_obs <- sqrt(sum((a - ctr)^2))
    (sum(d_sim >= d_obs) + 1) / (ncol(S) + 1)
  }, numeric(1L))
  fb <- fisher_combine(p_band)

  structure(list(
    rank = list(ranks = ranks, chisq = chi2, pval = unname(chi2$p.value),
                reject = chi2$p.value < alpha),
    perm_es = perm("ES"),
    perm_vs = perm("VS"),
    band = list(pvalues = p_band, fisher_stat = fb$stat, fisher_df = fb$df,
                fisher_pval = fb$pval, reject = fb$pval < alpha),
    n_folds = K, nsim = m, alpha = alpha),
    class = "doob_calibration")
}
