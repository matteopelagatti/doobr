## Quantise the innovation bank into `b` atoms with probabilities, then impose
## the discrete martingale-difference condition exactly.
##
## Atoms are the within-stratum means of `b` equiprobable strata of the bank
## (the stratified conditional-mean quantiser, which minimises within-stratum
## variance for equal probabilities).  They are then recentred so that
## sum(q * e) = 0 and rescaled so that sum(q * e^2) = 1.  Both corrections are
## O(n^{-1/2}) and buy the exact finite-sample conditions.
quantise_bank <- function(bank, b) {
  b <- as.integer(b)
  if (b < 2L) stop("each stage needs at least 2 branches", call. = FALSE)
  z <- sort(bank)
  n <- length(z)
  if (n < b)
    stop(sprintf("innovation bank has %d atoms, cannot form %d strata", n, b),
         call. = FALSE)
  cuts  <- floor(seq(0, n, length.out = b + 1L))
  atoms <- vapply(seq_len(b), function(i) mean(z[(cuts[i] + 1L):cuts[i + 1L]]),
                  numeric(1L))
  probs <- diff(cuts) / n
  atoms <- atoms - sum(probs * atoms)                 # sum q e  = 0
  s2    <- sum(probs * atoms^2)
  if (s2 <= 0) stop("degenerate quantisation of the innovation bank", call. = FALSE)
  atoms <- atoms / sqrt(s2)                           # sum q e^2 = 1
  list(atoms = atoms, probs = probs)
}


#' Build a scenario tree from a fitted model
#'
#' Branches on realisations of the martingale component rather than of the raw
#' process, so that the tree inherits the Doob decomposition node by node.  At
#' every non-leaf node the children satisfy \eqn{\sum_c q_c e_c = 0} and
#' \eqn{\sum_c q_c e_c^2 = 1} exactly, whence the conditional expectation of the
#' next-stage value on the tree equals the model forecast exactly.
#'
#' @param object A `"doob_fit"` object.
#' @param horizon Number of stages.  Defaults to the fitted horizon.
#' @param branching Integer vector of length `horizon`, or a single integer
#'   recycled, giving the number of children per node at each stage.  The tree
#'   has `prod(branching)` leaves.
#' @param newxreg Matrix of deterministic regressors with `horizon` rows, or
#'   `NULL`.
#' @param max_nodes Safety limit; an informative error is raised if the
#'   requested tree would exceed it.
#'
#' @return An object of class `"doob_tree"`: a list with
#'   \describe{
#'     \item{`nodes`}{A data frame with one row per node and columns `id`,
#'       `parent`, `stage`, `value` (original scale), `value_model`,
#'       `cond_prob` and `prob` (unconditional).}
#'     \item{`branching`,`horizon`,`n_leaves`}{Shape of the tree.}
#'   }
#'
#' @seealso [simulate.doob_fit()]
#' @export
#' @examples
#' set.seed(1)
#' y <- as.numeric(arima.sim(list(ar = 0.6), 300)) + 20
#' fit <- doob_fit(y, horizon = 3, mean_model_fn = lm_model_fn, mlags = 1,
#'                 predict_mean_fn = lm_predict_fn,
#'                 run_md_test = FALSE, eval_folds = FALSE)
#' tr <- doob_tree(fit, horizon = 3, branching = c(4, 3, 2))
#' tr
doob_tree <- function(object, horizon = NULL, branching, newxreg = NULL,
                      max_nodes = 1e5) {
  stopifnot(inherits(object, "doob_fit"))
  spec <- object$spec
  if (is.null(horizon)) horizon <- object$horizon
  horizon <- as.integer(horizon)
  if (missing(branching) || is.null(branching))
    stop("`branching` is required for a scenario tree", call. = FALSE)
  branching <- as.integer(branching)
  if (length(branching) == 1L) branching <- rep(branching, horizon)
  if (length(branching) != horizon)
    stop(sprintf("`branching` must have length 1 or %d, got %d",
                 horizon, length(branching)), call. = FALSE)

  n_leaves <- prod(as.numeric(branching))
  n_nodes  <- sum(cumprod(as.numeric(branching)))
  if (n_nodes > max_nodes)
    stop(sprintf(paste0("the requested tree has %.0f nodes (%.0f leaves), above ",
      "max_nodes = %.0f; reduce `branching` or raise `max_nodes`"),
      n_nodes, n_leaves, max_nodes), call. = FALSE)

  mlags <- spec$mlags; vlags <- spec$vlags
  mln <- paste0("lag", mlags); vln <- paste0("lag", vlags)
  Mw  <- max(c(mlags, if (spec$vlag_source == "y") vlags else 0L))
  Mq  <- if (spec$vlag_source == "sq_res") max(vlags) else 0L

  ## node state: the tail of the modelling series and of the squared residuals
  root <- list(w = utils::tail(spec$w, Mw),
               q = if (Mq > 0L) utils::tail(spec$sq_series, Mq) else numeric(0),
               level = spec$y_last, prob = 1)

  nodes <- data.frame(id = 1L, parent = NA_integer_, stage = 0L,
                      value = spec$trf$inv(spec$y_last),
                      value_model = spec$y_last,
                      cond_prob = 1, prob = 1)
  states  <- list(root)
  current <- 1L
  next_id <- 2L

  for (t in seq_len(horizon)) {
    qb  <- quantise_bank(object$bank, branching[t])
    xt  <- if (!is.null(newxreg)) newxreg[t, , drop = FALSE] else NULL
    new_states <- list(); new_rows <- list()

    for (ci in seq_along(current)) {
      st <- states[[ci]]
      yl <- stats::setNames(st$w[length(st$w) - mlags + 1L], mln)
      mt <- as.numeric(spec$predict_mean_fn(object$models$mean, yl, xt))

      if (!spec$has_var) {
        v_raw <- max(spec$const_v, .EPS)
      } else {
        vl <- if (spec$vlag_source == "y")
                stats::setNames(st$w[length(st$w) - vlags + 1L], vln)
              else stats::setNames(st$q[length(st$q) - vlags + 1L], vln)
        g  <- as.numeric(spec$predict_variance_fn(object$models$scale, vl, xt))
        if (!is.finite(g)) g <- 0
        v_raw <- if (spec$var_target == "squared") max(g, .EPS)
                 else exp(2 * min(max(g, -.LAM_CAP), .LAM_CAP))
      }
      sd_t <- spec$c_hat * min(max(sqrt(v_raw), spec$s_lo), spec$s_hi)

      for (j in seq_along(qb$atoms)) {
        inc   <- mt + sd_t * qb$atoms[j]          # modelling-scale increment
        level <- if (spec$difference) st$level + inc else inc
        new_rows[[length(new_rows) + 1L]] <- data.frame(
          id = next_id, parent = nodes$id[nodes$stage == t - 1L][ci],
          stage = t, value = spec$trf$inv(level), value_model = level,
          cond_prob = qb$probs[j], prob = st$prob * qb$probs[j])
        new_states[[length(new_states) + 1L]] <- list(
          w = c(st$w[-1L], inc),
          q = if (Mq > 0L) c(st$q[-1L], (inc - mt)^2) else numeric(0),
          level = level, prob = st$prob * qb$probs[j])
        next_id <- next_id + 1L
      }
    }
    nodes   <- rbind(nodes, do.call(rbind, new_rows))
    states  <- new_states
    current <- seq_along(new_states)
  }

  structure(list(nodes = nodes, branching = branching, horizon = horizon,
                 n_leaves = as.integer(n_leaves),
                 transform = spec$transform),
            class = "doob_tree")
}
