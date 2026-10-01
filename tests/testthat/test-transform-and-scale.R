pos_series <- function(n = 400, seed = 11) {
  set.seed(seed)
  exp(cumsum(stats::rnorm(n, 0, 0.02)) + 5)
}

test_that("transform='log' equals logging by hand", {
  y <- pos_series()
  common <- list(horizon = 10, mean_model_fn = lm_model_fn, mlags = 1,
                 predict_mean_fn = lm_predict_fn, difference = TRUE,
                 run_md_test = FALSE, eval_folds = FALSE, seed = 2,
                 verbose = FALSE)
  a <- do.call(doob_fit, c(list(y = y, transform = "log"), common))
  b <- do.call(doob_fit, c(list(y = log(y), transform = "none"), common))
  expect_equal(a$bank, b$bank, tolerance = 1e-12)
  expect_equal(a$c_hat, b$c_hat, tolerance = 1e-12)

  pa <- simulate(a, nsim = 50, seed = 5)
  pb <- simulate(b, nsim = 50, seed = 5)
  expect_equal(pa$paths_model, pb$paths, tolerance = 1e-12)
  expect_equal(pa$paths, exp(pb$paths), tolerance = 1e-12)
})

test_that("the smearing factor matches exp(s^2/2) and is 1 without a transform", {
  y <- pos_series()
  f <- doob_fit(y, horizon = 10, mean_model_fn = lm_model_fn, mlags = 1,
                predict_mean_fn = lm_predict_fn, transform = "log",
                difference = TRUE, run_md_test = FALSE, eval_folds = FALSE,
                seed = 2, verbose = FALSE)
  p <- simulate(f, nsim = 4000, seed = 7)
  expect_equal(p$smearing, exp(apply(p$paths_model, 1, var) / 2),
               tolerance = 0.01)
  ## the ensemble already carries the correct conditional mean
  expect_equal(rowMeans(p$paths),
               exp(rowMeans(p$paths_model)) * p$smearing, tolerance = 1e-10)

  g <- doob_fit(as.numeric(y), horizon = 10, mean_model_fn = lm_model_fn,
                mlags = 1, predict_mean_fn = lm_predict_fn,
                run_md_test = FALSE, eval_folds = FALSE, seed = 2,
                verbose = FALSE)
  expect_equal(simulate(g, nsim = 20, seed = 7)$smearing, rep(1, 10))
})

test_that("log of a non-positive series is refused", {
  expect_error(doob_fit(c(1, 2, -1, 4), horizon = 2, mean_model_fn = lm_model_fn,
                        mlags = 1, predict_mean_fn = lm_predict_fn,
                        transform = "log"),
               "strictly positive")
})

test_that("states where the conditional scale fails are excluded, not clipped", {
  ## a least-squares fit to squared residuals readily goes out of range
  set.seed(3)
  y <- as.numeric(stats::arima.sim(list(ar = 0.5), 600)) * 1000
  f <- doob_fit(y, horizon = 20, mean_model_fn = lm_model_fn, mlags = 1,
                predict_mean_fn = lm_predict_fn,
                variance_model_fn = lm_model_fn, vlags = 1,
                predict_variance_fn = lm_predict_fn,
                var_target = "squared", run_md_test = FALSE,
                eval_folds = FALSE, seed = 2, verbose = FALSE)
  expect_true(f$n_scale_failures >= 0L)
  expect_equal(length(f$bank), length(f$bank_idx))
  expect_lt(abs(mean(f$bank)), 1e-12)
  expect_lt(abs(mean(f$bank^2) - 1), 1e-12)
  expect_true(all(is.finite(f$oof_var[f$bank_idx])))
})

test_that("xreg must have one row per observation and newxreg is required", {
  set.seed(6)
  y <- as.numeric(stats::arima.sim(list(ar = 0.4), 300)) + 10
  X <- cbind(a = stats::rnorm(length(y)))
  expect_error(doob_fit(y, horizon = 5, mean_model_fn = lm_model_fn, mlags = 1,
                        predict_mean_fn = lm_predict_fn, xreg = X[-1, , drop = FALSE]),
               "one row per observation")
  f <- doob_fit(y, horizon = 5, mean_model_fn = lm_model_fn, mlags = 1,
                predict_mean_fn = lm_predict_fn, xreg = X,
                run_md_test = FALSE, eval_folds = FALSE, seed = 2, verbose = FALSE)
  expect_error(simulate(f, nsim = 5), "newxreg")
  p <- simulate(f, nsim = 5, newxreg = cbind(a = stats::rnorm(5)), seed = 1)
  expect_equal(dim(p$paths), c(5L, 5L))
})
