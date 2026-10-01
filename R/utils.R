## Internal helpers.  Nothing here is exported.

.EPS     <- .Machine$double.eps
.LAM_CAP <- 50   # hard cap on |log sigma| so exp() stays finite

## Variance-stabilising transformations -------------------------------------
## The series is supplied on its ORIGINAL scale; doob_fit() transforms it,
## models on the transformed scale, and simulate() maps the paths back.  No
## correction is involved and none is needed: simulation is a statement about
## the conditional law, and laws transform by change of variables.  See
## vignette material in DESIGN.md.
.transforms <- list(
  none = list(fwd = function(x) x, inv = function(x) x),
  log  = list(fwd = log,           inv = exp),
  sqrt = list(fwd = sqrt,          inv = function(x) x^2)
)

check_transform <- function(name, y) {
  if (!name %in% names(.transforms))
    stop("`transform` must be one of: ", paste(names(.transforms), collapse = ", "),
         call. = FALSE)
  if (name == "log" && any(y <= 0, na.rm = TRUE))
    stop(sprintf("transform = \"log\" needs a strictly positive series; %d value(s) are <= 0",
                 sum(y <= 0, na.rm = TRUE)), call. = FALSE)
  if (name == "sqrt" && any(y < 0, na.rm = TRUE))
    stop("transform = \"sqrt\" needs a non-negative series", call. = FALSE)
  .transforms[[name]]
}

## MAD-based scale, with a root-mean-square fallback -------------------------
robust_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  s <- 1.4826 * stats::median(abs(x - stats::median(x)))
  if (!is.finite(s) || s <= 0) s <- sqrt(mean((x - mean(x))^2))
  s
}

## (length(idx), length(lags)) matrix with column j equal to x[idx - lags[j]]
build_lag_matrix <- function(x, idx, lags) {
  mat <- matrix(NA_real_, nrow = length(idx), ncol = length(lags))
  for (j in seq_along(lags)) mat[, j] <- x[idx - lags[j]]
  colnames(mat) <- paste0("lag", lags)
  mat
}

## Row-wise prediction: preserves the single-row contract of the predict
## functions, which take a NAMED numeric vector of lags and a 1-row xreg.
predict_rows <- function(pred_fn, model, lag_mat, xreg_mat) {
  vapply(seq_len(nrow(lag_mat)), function(i)
    as.numeric(pred_fn(model,
                       stats::setNames(lag_mat[i, ], colnames(lag_mat)),
                       if (is.null(xreg_mat)) NULL else xreg_mat[i, , drop = FALSE])),
    numeric(1L))
}

#' Combine independent p-values by Fisher's method
#'
#' @param pvals Numeric vector of p-values.  Non-finite entries are dropped and
#'   the remainder are clipped away from zero so that a single exact zero does
#'   not make the statistic infinite.
#'
#' @return A list with components `stat`, `df` and `pval`.
#' @export
#' @examples
#' fisher_combine(c(0.4, 0.02, 0.31, 0.77))
fisher_combine <- function(pvals) {
  p <- pvals[is.finite(pvals)]
  if (!length(p)) return(list(stat = NA_real_, df = 0L, pval = NA_real_))
  p    <- pmin(pmax(p, 1e-300), 1)
  stat <- -2 * sum(log(p))
  df   <- 2L * length(p)
  list(stat = stat, df = df, pval = stats::pchisq(stat, df, lower.tail = FALSE))
}
