## Core recursion, shared by simulate.doob_fit() and evaluate_folds().
##
## Returns an (eff_h x nsim) matrix on the MODELLING scale (differenced if
## difference = TRUE).  `scale_const` is the constant paired with `pool`: the
## two must always move together, so that a leave-one-fold-out bank is used with
## its own constant.
run_recursion <- function(spec, mean_model, var_model, w_hist, sq_hist,
                          pool, eff_h, nsim, xreg_fut, scale_const) {
  mlags <- spec$mlags; vlags <- spec$vlags
  mln <- paste0("lag", mlags); vln <- paste0("lag", vlags)
  nh  <- length(w_hist)
  nq  <- if (spec$vlag_source == "sq_res") length(sq_hist) else 0L
  out <- matrix(NA_real_, nrow = eff_h, ncol = nsim)

  for (s in seq_len(nsim)) {
    wb <- c(w_hist, numeric(eff_h))
    qb <- if (spec$vlag_source == "sq_res") c(sq_hist, numeric(eff_h)) else NULL
    for (t in seq_len(eff_h)) {
      iy <- nh + t
      yl <- stats::setNames(wb[iy - mlags], mln)
      xt <- if (!is.null(xreg_fut)) xreg_fut[t, , drop = FALSE] else NULL
      mt <- as.numeric(spec$predict_mean_fn(mean_model, yl, xt))

      if (!spec$has_var) {
        v_raw <- max(spec$const_v, .EPS)
      } else {
        vl <- if (spec$vlag_source == "y") stats::setNames(wb[iy - vlags], vln)
              else stats::setNames(qb[nq + t - vlags], vln)
        g  <- as.numeric(spec$predict_variance_fn(var_model, vl, xt))
        if (!is.finite(g)) g <- 0
        v_raw <- if (spec$var_target == "squared") max(g, .EPS)
                 else exp(2 * min(max(g, -.LAM_CAP), .LAM_CAP))
      }
      ## the RAW scale is clipped to the admissible range, then multiplied by
      ## the constant paired with `pool`
      sd_t <- scale_const * min(max(sqrt(v_raw), spec$s_lo), spec$s_hi)
      wb[iy] <- mt + sd_t * sample(pool, 1L)
      if (spec$vlag_source == "sq_res") qb[nq + t] <- (wb[iy] - mt)^2
    }
    out[, s] <- wb[(nh + 1L):(nh + eff_h)]
  }
  out
}

## Integrate (if differenced) and map back to the original scale.
## `ref` is the last transformed-scale level before the first simulated step.
to_original <- function(spec, raw, ref) {
  m <- if (spec$difference) ref + apply(raw, 2L, cumsum) else raw
  if (!is.matrix(m)) m <- matrix(m, nrow = nrow(raw))
  list(model = m, original = spec$trf$inv(m))
}


#' Simulate sample paths or a scenario tree from a fitted model
#'
#' @param object A `"doob_fit"` object from [doob_fit()].
#' @param nsim Number of paths when `output = "paths"`.  Ignored for trees,
#'   whose size is set by `branching`.
#' @param seed Optional integer seed.
#' @param horizon Steps ahead.  Defaults to the horizon the model was fitted
#'   for; a different value is allowed but the calibration diagnostics in
#'   `object` only validate the horizon they were run at, and a warning is
#'   issued.
#' @param output `"paths"` for a flat ensemble, `"tree"` for a scenario tree.
#' @param branching Integer vector of length `horizon` giving the number of
#'   children per node at each stage, or a single integer recycled.  Required
#'   shape for `output = "tree"`; the tree has `prod(branching)` leaves.
#' @param newxreg Matrix of deterministic regressors for the future, with
#'   `horizon` rows.  Required when the model was fitted with regressors.
#' @param ... Unused.
#'
#' @return For `output = "paths"`, an object of class `"doob_paths"`: a list
#'   with `paths` (a `horizon` x `nsim` matrix on the original scale),
#'   `paths_model` (the modelling-scale counterpart) and `smearing`.  For
#'   `output = "tree"`, an object of class `"doob_tree"`; see [doob_tree()].
#'
#' @details
#' The simulated increments are a martingale difference sequence on the
#' MODELLING scale, which is where the Doob decomposition is applied.  When a
#' nonlinear `transform` was used, `smearing` reports the per-step ratio between
#' the ensemble mean on the original scale and the back-transform of the
#' modelling-scale mean.  This is a descriptive feature of a skewed conditional
#' law, not a bias: the ensemble already carries the correct conditional mean.
#' It is reported because `inv(mu)` is the conditional median, not the
#' conditional mean, so a decision problem formulated on the original scale
#' should take its conditional mean from the ensemble or as `inv(mu) * smearing`.
#'
#' @seealso [doob_fit()], [doob_tree()]
#' @importFrom stats simulate
#' @exportS3Method stats::simulate
#' @examples
#' set.seed(1)
#' y <- as.numeric(arima.sim(list(ar = 0.6), 300)) + 20
#' fit <- doob_fit(y, horizon = 10, mean_model_fn = lm_model_fn, mlags = 1,
#'                 predict_mean_fn = lm_predict_fn,
#'                 run_md_test = FALSE, eval_folds = FALSE)
#' p <- simulate(fit, nsim = 50)
#' dim(p$paths)
simulate.doob_fit <- function(object, nsim = 100L, seed = NULL,
                              horizon = NULL,
                              output = c("paths", "tree"),
                              branching = NULL, newxreg = NULL, ...) {
  output <- match.arg(output)
  if (!is.null(seed)) set.seed(seed)
  spec <- object$spec
  if (is.null(horizon)) horizon <- object$horizon
  horizon <- as.integer(horizon)
  stopifnot(horizon >= 1L)

  if (horizon != object$fold_length)
    warning(sprintf(paste0("simulating %d steps from a model whose evaluation folds ",
      "are %d long; the calibration diagnostics in this object validate the ",
      "%d-step horizon only"), horizon, object$fold_length, object$fold_length),
      call. = FALSE)

  if (spec$n_xreg > 0L) {
    if (is.null(newxreg))
      stop("`newxreg` is required: the model was fitted with regressors",
           call. = FALSE)
    newxreg <- as.matrix(newxreg)
    if (nrow(newxreg) < horizon)
      stop(sprintf("`newxreg` needs at least %d rows, got %d", horizon,
                   nrow(newxreg)), call. = FALSE)
    if (ncol(newxreg) != spec$n_xreg)
      stop(sprintf("`newxreg` has %d columns, the model was fitted with %d",
                   ncol(newxreg), spec$n_xreg), call. = FALSE)
    newxreg <- newxreg[seq_len(horizon), , drop = FALSE]
  } else if (!is.null(newxreg)) {
    warning("`newxreg` ignored: the model was fitted without regressors",
            call. = FALSE)
    newxreg <- NULL
  }

  if (output == "tree")
    return(doob_tree(object, horizon = horizon, branching = branching,
                     newxreg = newxreg))

  nsim <- as.integer(nsim)
  stopifnot(nsim >= 1L)
  raw <- run_recursion(spec, object$models$mean, object$models$scale,
                       spec$w, spec$sq_series, object$bank, horizon, nsim,
                       newxreg, spec$c_hat)
  both <- to_original(spec, raw, spec$y_last)
  smear <- if (spec$transform == "none") rep(1, horizon)
           else rowMeans(both$original) / spec$trf$inv(rowMeans(both$model))

  structure(list(paths = both$original, paths_model = both$model,
                 smearing = smear, horizon = horizon, nsim = nsim,
                 transform = spec$transform),
            class = "doob_paths")
}
