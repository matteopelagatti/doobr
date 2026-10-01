test_that("md_test does not reject white noise and does reject an AR(1)", {
  set.seed(1)
  expect_gt(md_test(stats::rnorm(500), type = "additive", test = "F"), 0.05)
  x <- as.numeric(stats::arima.sim(list(ar = 0.5), 500))
  expect_lt(md_test(x, type = "additive", test = "F"), 0.01)
})

test_that("md_test returns NA rather than failing on a degenerate input", {
  expect_true(is.na(md_test(rep(1, 5), lags = 1)))
})

test_that("fisher_combine behaves", {
  expect_equal(fisher_combine(c(0.5, 0.5))$df, 4L)
  expect_true(fisher_combine(c(1e-300, 0.5))$pval < 0.01)
  expect_true(is.na(fisher_combine(numeric(0))$pval))
})

test_that("calibration_tests runs and returns p-values in range", {
  set.seed(1)
  y <- as.numeric(stats::arima.sim(list(ar = 0.6), 400)) + 20
  f <- doob_fit(y, horizon = 10, mean_model_fn = lm_model_fn, mlags = 1,
                predict_mean_fn = lm_predict_fn, nsim_eval = 40,
                run_md_test = FALSE, seed = 2, verbose = FALSE)
  cal <- calibration_tests(f)
  expect_s3_class(cal, "doob_calibration")
  for (p in c(cal$rank$pval, cal$perm_es$fisher_pval,
              cal$perm_vs$fisher_pval, cal$band$fisher_pval)) {
    expect_gte(p, 0); expect_lte(p, 1)
  }
  expect_error(calibration_tests(
    doob_fit(y, horizon = 10, mean_model_fn = lm_model_fn, mlags = 1,
             predict_mean_fn = lm_predict_fn, eval_folds = FALSE,
             run_md_test = FALSE, seed = 2, verbose = FALSE)),
    "eval_folds")
})

test_that("scores are non-negative and vs_score is 0 at horizon 1", {
  set.seed(1)
  expect_gte(es_score(stats::rnorm(5), matrix(stats::rnorm(100), 5)), 0)
  expect_equal(vs_score(1, matrix(stats::rnorm(20), 1)), 0)
})
