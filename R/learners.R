#' Learner adapters
#'
#' Matched `model_fn` / `predict_fn` pairs for common learners, in the form
#' [doob_fit()] expects.  Any learner can be used by supplying a pair with these
#' signatures:
#'
#' \preformatted{
#' model_fn(y, y_lags, xreg) -> model
#' predict_fn(model, lags_vec, xreg_row) -> scalar
#' }
#'
#' where `y_lags` is an \eqn{m \times p} matrix with columns named `"lag1"`,
#' `"lag12"` and so on, `lags_vec` is a single named row of it, and `xreg` /
#' `xreg_row` are a matrix and a 1-row matrix or `NULL`.
#'
#' @param y Numeric response vector.
#' @param y_lags Numeric matrix of lagged predictors.
#' @param xreg Numeric matrix of regressors, or `NULL`.
#' @param model A fitted model as returned by the matching `model_fn`.
#' @param lags_vec Named numeric vector of lags for one time point.
#' @param xreg_row A 1-row matrix of regressors, or `NULL`.
#' @param seed Optional seed, for learners that are stochastic.  Supplying one
#'   makes the fit reproducible, which matters: an unseeded random forest gives
#'   different diagnostics on every run.
#' @param kernel,cost,epsilon Passed to [e1071::svm()].
#' @param nrounds,params Passed to `xgboost::xgboost()`.
#' @param ... Further arguments passed to the underlying fitting function.
#'
#' @return `*_model_fn` returns a fitted model; `*_predict_fn` returns a scalar.
#' @name learners
#' @examples
#' set.seed(1)
#' y <- rnorm(50); L <- cbind(lag1 = rnorm(50))
#' m <- lm_model_fn(y, L, NULL)
#' lm_predict_fn(m, c(lag1 = 0.3), NULL)
NULL

#' @rdname learners
#' @export
lm_model_fn <- function(y, y_lags, xreg = NULL) {
  X <- cbind(1, y_lags)
  if (!is.null(xreg)) X <- cbind(X, as.matrix(xreg))
  stats::lm.fit(X, y)
}

#' @rdname learners
#' @export
lm_predict_fn <- function(model, lags_vec, xreg_row = NULL) {
  b <- model$coefficients
  b[is.na(b)] <- 0
  x <- if (is.null(xreg_row)) c(1, lags_vec)
       else c(1, lags_vec, as.numeric(xreg_row))
  as.numeric(b %*% x)
}

#' @rdname learners
#' @export
rf_model_fn <- function(y, y_lags, xreg = NULL, seed = NULL, ...) {
  if (!requireNamespace("randomForest", quietly = TRUE))
    stop("package 'randomForest' is needed for rf_model_fn", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  df <- as.data.frame(y_lags)
  if (!is.null(xreg)) df <- cbind(df, as.data.frame(xreg))
  randomForest::randomForest(x = df, y = y, ...)
}

#' @rdname learners
#' @export
rf_predict_fn <- function(model, lags_vec, xreg_row = NULL) {
  nd <- as.data.frame(t(lags_vec))
  if (!is.null(xreg_row)) nd <- cbind(nd, as.data.frame(xreg_row))
  unname(stats::predict(model, newdata = nd))
}

#' @rdname learners
#' @export
svm_model_fn <- function(y, y_lags, xreg = NULL, kernel = "radial",
                         cost = 1, epsilon = 0.1, ...) {
  if (!requireNamespace("e1071", quietly = TRUE))
    stop("package 'e1071' is needed for svm_model_fn", call. = FALSE)
  df <- as.data.frame(y_lags)
  if (!is.null(xreg)) df <- cbind(df, as.data.frame(xreg))
  e1071::svm(x = df, y = y, kernel = kernel, cost = cost, epsilon = epsilon, ...)
}

#' @rdname learners
#' @export
svm_predict_fn <- function(model, lags_vec, xreg_row = NULL) {
  nd <- as.data.frame(t(lags_vec))
  if (!is.null(xreg_row)) nd <- cbind(nd, as.data.frame(xreg_row))
  unname(stats::predict(model, newdata = nd))
}

#' @rdname learners
#' @export
xgb_model_fn <- function(y, y_lags, xreg = NULL, nrounds = 100L,
                         params = list(), ...) {
  if (!requireNamespace("xgboost", quietly = TRUE))
    stop("package 'xgboost' is needed for xgb_model_fn", call. = FALSE)
  X <- as.matrix(y_lags)
  if (!is.null(xreg)) X <- cbind(X, as.matrix(xreg))
  xgboost::xgboost(data = X, label = y, nrounds = nrounds, params = params,
                   verbose = 0, ...)
}

#' @rdname learners
#' @export
xgb_predict_fn <- function(model, lags_vec, xreg_row = NULL) {
  x <- matrix(as.numeric(lags_vec), nrow = 1L)
  if (!is.null(xreg_row)) x <- cbind(x, as.matrix(xreg_row))
  as.numeric(stats::predict(model, newdata = x))
}
