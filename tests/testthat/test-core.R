ar1 <- function(n = 400, phi = 0.6, seed = 1) {
  set.seed(seed)
  as.numeric(stats::arima.sim(list(ar = phi), n)) + 20
}

basic_fit <- function(y = ar1(), ...) {
  doob_fit(y, horizon = 10, mean_model_fn = lm_model_fn, mlags = 1,
           predict_mean_fn = lm_predict_fn, run_md_test = FALSE,
           eval_folds = FALSE, seed = 2, verbose = FALSE, ...)
}

test_that("the innovation bank has mean 0 and variance 1 exactly", {
  f <- basic_fit()
  expect_lt(abs(mean(f$bank)), 1e-12)
  expect_lt(abs(mean(f$bank^2) - 1), 1e-12)
})

test_that("a constant conditional scale reproduces the out-of-fold variance", {
  f <- basic_fit()
  v <- f$oof_var[f$bank_idx]
  expect_true(all(abs(v - v[1]) < 1e-12))              # constant across states
  r <- f$oof_resid[f$bank_idx]
  expect_equal(unname(v[1]), mean((r - mean(r))^2), tolerance = 1e-10)
  expect_equal(f$n_scale_failures, 0L)
})

test_that("c_hat is sqrt((n-1)/n) in the homoskedastic case", {
  f <- basic_fit()
  n <- length(f$bank)
  expect_equal(f$c_hat, sqrt((n - 1) / n), tolerance = 1e-12)
})

test_that("simulate() returns the requested shape", {
  f <- basic_fit()
  p <- simulate(f, nsim = 37, seed = 3)
  expect_s3_class(p, "doob_paths")
  expect_equal(dim(p$paths), c(10L, 37L))
  expect_true(all(is.finite(p$paths)))
})

test_that("a horizon different from the fold length warns", {
  f <- basic_fit()
  expect_warning(simulate(f, nsim = 5, horizon = 4, seed = 3),
                 "calibration diagnostics")
})

test_that("the tree satisfies the conditional-mean and variance identities", {
  for (dif in c(FALSE, TRUE)) {
    f  <- doob_fit(ar1(), horizon = 3, mean_model_fn = lm_model_fn, mlags = 1,
                   predict_mean_fn = lm_predict_fn, difference = dif,
                   run_md_test = FALSE, eval_folds = FALSE, seed = 2,
                   verbose = FALSE)
    tr <- doob_tree(f, horizon = 3, branching = c(5, 4, 3))
    s  <- f$spec
    yl <- stats::setNames(utils::tail(s$w, max(s$mlags))[max(s$mlags) - s$mlags + 1L],
                          paste0("lag", s$mlags))
    mu <- as.numeric(s$predict_mean_fn(f$models$mean, yl, NULL))
    ch <- tr$nodes[tr$nodes$stage == 1, ]
    lhs <- sum(ch$cond_prob * ch$value_model)
    rhs <- if (dif) s$y_last + mu else mu
    expect_equal(lhs, rhs, tolerance = 1e-10)

    v_emp <- sum(ch$cond_prob * (ch$value_model - lhs)^2)
    sd_t  <- s$c_hat * min(max(sqrt(max(s$const_v, .Machine$double.eps)),
                               s$s_lo), s$s_hi)
    expect_equal(v_emp, sd_t^2, tolerance = 1e-10)
  }
})

test_that("tree probabilities sum to one at every stage", {
  f  <- basic_fit()
  tr <- doob_tree(f, horizon = 3, branching = 3)
  for (t in 0:3)
    expect_equal(sum(tr$nodes$prob[tr$nodes$stage == t]), 1, tolerance = 1e-12)
  expect_equal(tr$n_leaves, 27L)
})

test_that("quantise_bank enforces the discrete moment conditions", {
  set.seed(4)
  qb <- doobr:::quantise_bank(stats::rnorm(500), 7)
  expect_equal(sum(qb$probs), 1, tolerance = 1e-12)
  expect_lt(abs(sum(qb$probs * qb$atoms)), 1e-12)
  expect_lt(abs(sum(qb$probs * qb$atoms^2) - 1), 1e-12)
})

test_that("an oversized tree is refused with a useful message", {
  f <- basic_fit()
  expect_error(doob_tree(f, horizon = 10, branching = 10, max_nodes = 1e4),
               "max_nodes")
})
