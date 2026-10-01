#' Fit a Doob scenario-generation model by purged cross-fitting
#'
#' Fits the conditional mean and conditional scale of a univariate time series
#' with user-supplied learners, using purged \eqn{K}-fold cross-fitting over
#' contiguous blocks whose length equals the intended simulation horizon.  The
#' returned object carries the full-sample models, an honest innovation bank,
#' and out-of-sample diagnostics; pass it to [simulate.doob_fit()] to generate
#' sample paths or a scenario tree.
#'
#' @section Why cross-fitting:
#' The residuals of a model evaluated at the observations used to fit it are
#' systematically too small and retain structure that is an artefact of the
#' fitting procedure, so resampling them injects artificial predictability into
#' the generated paths.  Each fold here is simultaneously a source of
#' out-of-fold innovations, an evaluation window over which paths are simulated
#' from models that never saw it, and a block of out-of-sample residuals for the
#' martingale-difference tests.
#'
#' @param y Numeric vector; the complete observed series on its ORIGINAL scale.
#'   Nothing is held out: every usable observation is predicted out-of-fold
#'   exactly once.
#' @param horizon Positive integer; the intended simulation horizon.  This is
#'   also the default fold length, which is what makes a fold double as an
#'   evaluation window.
#' @param mean_model_fn,predict_mean_fn Conditional-mean learner.  `mean_model_fn`
#'   has signature `function(y, y_lags, xreg) -> model`; `predict_mean_fn` has
#'   signature `function(model, lags_vec, xreg_row) -> scalar`, where `lags_vec`
#'   is a named numeric vector and `xreg_row` a 1-row matrix or `NULL`.
#' @param mlags Integer vector of lag offsets for the conditional mean.
#' @param variance_model_fn,predict_variance_fn,vlags Conditional-scale learner,
#'   same signatures.  `variance_model_fn = NULL` (the default) holds the
#'   conditional variance constant at the out-of-fold residual variance.
#' @param xreg Numeric matrix of deterministic regressors with one row per
#'   observation of `y`, or `NULL`.  Future values are supplied to
#'   [simulate.doob_fit()], not here.
#' @param transform One of `"none"`, `"log"`, `"sqrt"`.  Applied to `y` before
#'   modelling; simulated paths are mapped back automatically.
#' @param difference Logical; model the first difference of the (transformed)
#'   series and integrate simulated paths back.
#' @param var_target `"squared"` fits the conditional-scale learner to the
#'   squared residuals by least squares; `"logabs"` fits it to
#'   \eqn{\log|\hat\varepsilon|}, which guarantees positivity and needs only
#'   \eqn{\mathrm{Var}(\log|\eta|) < \infty} instead of \eqn{E\eta^4 < \infty}.
#' @param vlag_source `"y"` (default) indexes the conditional-scale model by
#'   lags of the modelling series, so the state is fully observed; `"sq_res"`
#'   reproduces the older behaviour of indexing it by lagged squared residuals.
#' @param n_folds Number of folds; `NULL` partitions the whole usable range.
#' @param fold_length Fold length; defaults to `horizon`.
#' @param purge Observations within `purge` of a fold are dropped from that
#'   fold's training set.  Defaults to the model memory.
#' @param trim_q Trimming quantile for `var_target = "logabs"`.
#' @param sigma_floor_frac,sigma_ceiling_mult Operational form of the assumption
#'   that the conditional scale is bounded away from zero and infinity.
#' @param eval_folds Logical; run the per-fold simulation and calibration
#'   diagnostics.
#' @param nsim_eval Number of paths simulated per evaluation fold.
#' @param run_md_test Logical; run the martingale-difference tests.
#' @param md_lags,md_zero_mean,md_df,md_q,md_type,md_stat,md_B Arguments passed
#'   to [md_test()].
#' @param seed Optional integer seed.
#' @param verbose Logical; report the martingale-difference p-values.
#'
#' @return An object of class `"doob_fit"`: a list with components
#'   \describe{
#'     \item{`bank`}{Normalised innovation bank, mean exactly 0 and variance
#'       exactly 1.}
#'     \item{`c_hat`}{Scale constant carried into the conditional variance; see
#'       Details.}
#'     \item{`oof_resid`,`oof_var`}{Out-of-fold residuals and conditional
#'       variances.}
#'     \item{`md`}{Martingale-difference p-values for step 1 (the residuals) and
#'       step 2 (the squared residuals net of the fitted variance).}
#'     \item{`eval`}{Per-fold evaluation windows and proper scores.}
#'     \item{`models`}{Full-sample conditional-mean and conditional-scale fits.}
#'     \item{`spec`}{Everything [simulate.doob_fit()] needs to continue the
#'       recursion.}
#'   }
#'
#' @details
#' The bank is normalised to mean exactly 0 and variance exactly 1, so that
#' simulated paths satisfy the conditional-mean and conditional-variance
#' identities exactly.  Normalisation alone would destroy the self-correction a
#' raw bank gives against a multiplicative bias in the fitted variance, so the
#' scale constant is carried into the variance rather than discarded.  For
#' `var_target = "logabs"` this is exactly the correction of the
#' \eqn{\exp(\kappa)} factor, \eqn{\kappa = E\log|\eta|}, which is about
#' \eqn{-0.635} for Gaussian innovations: omitting it understates the
#' conditional standard deviation by nearly a factor of two.
#'
#' States at which the fitted scale leaves the admissible range are excluded
#' from the bank rather than clipped into it.  They are states at which the
#' conditional-scale model failed, so their standardised residual is
#' meaningless; this removes states with a broken variance estimate, not
#' residuals with a large absolute value, so the innovation tails are untouched.
#'
#' @seealso [simulate.doob_fit()], [md_test()], [calibration_tests()]
#' @export
#' @examples
#' set.seed(1)
#' y <- as.numeric(arima.sim(list(ar = 0.6), 300)) + 20
#' fit <- doob_fit(y, horizon = 10,
#'                 mean_model_fn = lm_model_fn, mlags = 1,
#'                 predict_mean_fn = lm_predict_fn,
#'                 run_md_test = FALSE, eval_folds = FALSE)
#' fit
doob_fit <- function(y, horizon,
                     mean_model_fn, mlags, predict_mean_fn,
                     variance_model_fn = NULL, vlags = NULL,
                     predict_variance_fn = NULL,
                     xreg = NULL,
                     transform = c("none", "log", "sqrt"),
                     difference = FALSE,
                     var_target = c("squared", "logabs"),
                     vlag_source = c("y", "sq_res"),
                     n_folds = NULL, fold_length = NULL, purge = NULL,
                     trim_q = 0.01,
                     sigma_floor_frac = 0.05, sigma_ceiling_mult = 3,
                     eval_folds = TRUE, nsim_eval = 200L,
                     run_md_test = TRUE,
                     md_lags = NULL, md_zero_mean = TRUE,
                     md_df = 3L, md_q = 3L,
                     md_type = "random", md_stat = "WildBoot", md_B = 499L,
                     seed = NULL, verbose = TRUE) {

  transform   <- match.arg(transform)
  var_target  <- match.arg(var_target)
  vlag_source <- match.arg(vlag_source)
  if (!is.null(seed)) set.seed(seed)

  mlags   <- as.integer(mlags)
  horizon <- as.integer(horizon)
  stopifnot(all(mlags >= 1L), horizon >= 1L, nsim_eval >= 2L)
  M_m <- max(mlags)

  has_var <- !is.null(variance_model_fn)
  if (has_var) {
    if (is.null(vlags) || is.null(predict_variance_fn))
      stop("`vlags` and `predict_variance_fn` are required with `variance_model_fn`",
           call. = FALSE)
    vlags <- as.integer(vlags)
    stopifnot(all(vlags >= 1L))
    M_v <- max(vlags)
  } else {
    vlags <- 1L; M_v <- 1L
  }

  y     <- as.numeric(y)
  trf   <- check_transform(transform, y)
  y_t   <- trf$fwd(y)
  w     <- if (difference) diff(y_t) else y_t
  n_w   <- length(w)

  ## ---- regressors: one row per observation of y, trimmed to match w -------
  if (!is.null(xreg)) {
    xreg <- as.matrix(xreg)
    if (nrow(xreg) != length(y))
      stop(sprintf("`xreg` must have one row per observation of `y` (%d), got %d",
                   length(y), nrow(xreg)), call. = FALSE)
    xreg_hist <- if (difference) xreg[-1L, , drop = FALSE] else xreg
  } else {
    xreg_hist <- NULL
  }
  xr <- function(idx) if (is.null(xreg_hist)) NULL else xreg_hist[idx, , drop = FALSE]

  ## ---- memory, usable range, folds ---------------------------------------
  M <- if (has_var && vlag_source == "sq_res") M_m + M_v
       else if (has_var) max(M_m, M_v) else M_m
  if (is.null(purge)) purge <- M
  purge <- as.integer(purge)

  usable <- seq.int(M + 1L, n_w)
  if (length(usable) < 3L)
    stop("series too short for the requested lag structure", call. = FALSE)

  h <- as.integer(if (is.null(fold_length)) horizon else fold_length)
  K_max <- length(usable) %/% h
  if (K_max < 2L)
    stop(sprintf("only %d fold(s) of length %d fit in %d usable observations",
                 K_max, h, length(usable)), call. = FALSE)
  K <- if (is.null(n_folds)) K_max else min(as.integer(n_folds), K_max)
  if (K < 2L) stop("`n_folds` must be at least 2", call. = FALSE)
  tail_idx <- usable[(length(usable) - K * h + 1L):length(usable)]
  folds <- lapply(seq_len(K), function(k) tail_idx[((k - 1L) * h + 1L):(k * h)])

  train_idx <- function(fold) {
    lo <- fold[1L] - purge; hi <- fold[length(fold)] + purge
    usable[usable < lo | usable > hi]
  }

  ## =======================================================================
  ## pass 1 -- out-of-fold conditional-mean residuals
  ## =======================================================================
  oof_resid   <- rep(NA_real_, n_w)
  mean_models <- vector("list", K)
  for (k in seq_len(K)) {
    fold <- folds[[k]]
    tr   <- train_idx(fold)
    if (length(tr) <= length(mlags) + 1L)
      stop(sprintf("fold %d: training set too small after purging", k), call. = FALSE)
    mdl <- mean_model_fn(w[tr], build_lag_matrix(w, tr, mlags), xr(tr))
    mean_models[[k]] <- mdl
    oof_resid[fold] <- w[fold] -
      predict_rows(predict_mean_fn, mdl, build_lag_matrix(w, fold, mlags), xr(fold))
  }
  fitted_idx <- sort(unlist(folds))

  ## =======================================================================
  ## pass 2 -- out-of-fold conditional scale (raw, before the constant)
  ## =======================================================================
  oof_vraw   <- rep(NA_real_, n_w)
  var_models <- vector("list", K)
  kappa_hat  <- NA_real_
  tau        <- NA_real_
  const_v    <- NA_real_

  if (!has_var) {
    const_v <- stats::var(oof_resid[fitted_idx])
    oof_vraw[fitted_idx] <- max(const_v, .EPS)
  } else {
    if (var_target == "squared") {
      resp_all <- oof_resid^2
    } else {
      absr <- abs(oof_resid[fitted_idx]); absr <- absr[is.finite(absr) & absr > 0]
      tau  <- max(if (length(absr)) unname(stats::quantile(absr, trim_q)) else .EPS, .EPS)
      resp_all <- rep(NA_real_, n_w)
      resp_all[fitted_idx] <- log(pmax(abs(oof_resid[fitted_idx]), tau))
    }
    if (vlag_source == "y") {
      vlag_of <- function(idx) build_lag_matrix(w, idx, vlags)
      keep_ok <- function(idx) idx
    } else {
      sq <- rep(NA_real_, n_w); sq[fitted_idx] <- oof_resid[fitted_idx]^2
      vlag_of <- function(idx) build_lag_matrix(sq, idx, vlags)
      keep_ok <- function(idx)
        idx[vapply(idx, function(i) all(is.finite(sq[i - vlags])), logical(1L))]
    }
    oof_g <- rep(NA_real_, n_w)
    for (k in seq_len(K)) {
      fold <- folds[[k]]
      tr   <- keep_ok(train_idx(fold))
      tr   <- tr[is.finite(resp_all[tr])]
      Xtr  <- vlag_of(tr)
      good <- stats::complete.cases(Xtr)
      tr   <- tr[good]; Xtr <- Xtr[good, , drop = FALSE]
      if (length(tr) <= length(vlags) + 1L)
        stop(sprintf("fold %d: conditional-scale training set too small", k),
             call. = FALSE)
      vm <- variance_model_fn(resp_all[tr], Xtr, xr(tr))
      var_models[[k]] <- vm
      fo <- keep_ok(fold)
      oof_g[fo] <- predict_rows(predict_variance_fn, vm, vlag_of(fo), xr(fo))
    }
    valid <- fitted_idx[is.finite(oof_g[fitted_idx])]
    oof_vraw[valid] <- if (var_target == "squared")
      pmax(oof_g[valid], .EPS)
    else
      pmax(exp(2 * pmin(pmax(oof_g[valid], -.LAM_CAP), .LAM_CAP)), .EPS)
    fitted_idx <- valid
  }

  ## =======================================================================
  ## bank, scale constant, admissible range for the scale
  ## =======================================================================
  bank_idx <- fitted_idx[is.finite(oof_vraw[fitted_idx]) &
                           is.finite(oof_resid[fitted_idx])]
  sd_raw    <- sqrt(oof_vraw[bank_idx])
  res_scale <- robust_sd(oof_resid[bank_idx])
  s_lo <- sigma_floor_frac * res_scale
  s_hi <- sigma_ceiling_mult *
            max(unname(stats::quantile(sd_raw, 0.995)), res_scale)
  if (!is.finite(s_lo) || !is.finite(s_hi) || s_lo <= 0 || s_lo >= s_hi) {
    s_lo <- .EPS; s_hi <- Inf
  }
  bad <- sd_raw < s_lo | sd_raw > s_hi
  oof_vraw[bank_idx] <- pmin(pmax(sd_raw, s_lo), s_hi)^2
  if (all(bad))
    stop("the conditional-scale model failed at every state; check the learner ",
         "or use var_target = \"logabs\"", call. = FALSE)
  if (mean(bad) > 0.5)
    warning(sprintf(paste0("the fitted conditional scale was out of bounds at %.0f%% ",
      "of states; the bank is built from the remainder. This usually means the ",
      "conditional-scale model is badly specified (a least-squares fit to squared ",
      "residuals going negative is the common cause)."), 100 * mean(bad)),
      call. = FALSE)
  bank_idx <- bank_idx[!bad]

  z_raw <- oof_resid[bank_idx] / sqrt(oof_vraw[bank_idx])
  z0    <- z_raw - mean(z_raw)
  c_hat <- sqrt(mean(z0^2))
  if (!is.finite(c_hat) || c_hat <= 0)
    stop("degenerate scale constant: the innovation bank has no dispersion",
         call. = FALSE)
  bank <- z0 / c_hat
  if (var_target == "logabs") kappa_hat <- -log(c_hat)

  oof_var <- rep(NA_real_, n_w)
  oof_var[bank_idx] <- oof_vraw[bank_idx] * c_hat^2
  sq_series <- ifelse(is.finite(oof_resid), oof_resid^2, NA_real_)

  ## =======================================================================
  ## full-sample refit (the conditional scale is refitted on OUT-OF-FOLD
  ## responses, which target sigma^2 + MSE rather than the deflated in-sample
  ## analogue)
  ## =======================================================================
  full_mean <- mean_model_fn(w[usable], build_lag_matrix(w, usable, mlags), xr(usable))
  full_var  <- NULL
  if (has_var) {
    resp_full <- if (var_target == "squared") oof_resid[bank_idx]^2
                 else log(pmax(abs(oof_resid[bank_idx]), tau))
    Xfull <- if (vlag_source == "y") build_lag_matrix(w, bank_idx, vlags)
             else build_lag_matrix(sq_series, bank_idx, vlags)
    good  <- stats::complete.cases(Xfull) & is.finite(resp_full)
    full_var <- variance_model_fn(resp_full[good], Xfull[good, , drop = FALSE],
                                  xr(bank_idx[good]))
  }

  spec <- list(
    mlags = mlags, vlags = vlags, has_var = has_var,
    predict_mean_fn = predict_mean_fn, predict_variance_fn = predict_variance_fn,
    transform = transform, difference = difference,
    var_target = var_target, vlag_source = vlag_source,
    c_hat = c_hat, const_v = const_v, s_lo = s_lo, s_hi = s_hi,
    trf = trf, w = w, y_last = y_t[length(y_t)], sq_series = sq_series,
    n_xreg = if (is.null(xreg_hist)) 0L else ncol(xreg_hist)
  )

  obj <- structure(list(
    call        = match.call(),
    bank        = bank,
    bank_idx    = bank_idx,
    c_hat       = c_hat,
    kappa_hat   = kappa_hat,
    trim_tau    = tau,
    oof_resid   = oof_resid,
    oof_var     = oof_var,
    sigma_range = c(s_lo, s_hi),
    n_scale_failures = sum(bad),
    models      = list(mean = full_mean, scale = full_var,
                       fold_mean = mean_models, fold_scale = var_models),
    folds       = folds,
    n_folds     = K,
    fold_length = h,
    horizon     = horizon,
    spec        = spec,
    md          = list(step1 = NULL, step2 = NULL),
    eval        = NULL
  ), class = "doob_fit")

  ## =======================================================================
  ## diagnostics
  ## =======================================================================
  if (run_md_test) {
    p1 <- md_test(oof_resid[bank_idx], lags = md_lags, zero_mean = md_zero_mean,
                  df = md_df, q = md_q, type = md_type, test = md_stat, B = md_B)
    p2 <- md_test(oof_resid[bank_idx]^2 - oof_var[bank_idx], lags = md_lags,
                  zero_mean = md_zero_mean, df = md_df, q = md_q,
                  type = md_type, test = md_stat, B = md_B)
    obj$md <- list(step1 = p1, step2 = p2)
    if (verbose) {
      message(sprintf("MDH step 1  (eps_t         | F_t-1):  p = %.4f", p1))
      message(sprintf("MDH step 2  (eps_t^2 - v_t | F_t-1):  p = %.4f", p2))
    }
  }

  if (eval_folds)
    obj$eval <- evaluate_folds(obj, nsim = nsim_eval, xreg_hist = xreg_hist,
                               z_raw = z_raw, w = w, y_t = y_t)

  obj
}
