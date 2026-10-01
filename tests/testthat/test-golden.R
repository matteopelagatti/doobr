## Cross-language reference fixture.  doobpy must reproduce the same numbers
## from the same input; see DESIGN.md for the tolerances and data-raw/ for the
## generator.
test_that("the golden fixture is reproduced", {
  inp <- utils::read.csv(system.file("extdata", "golden_input.csv",
                                     package = "doobr"))
  exp_df <- utils::read.csv(system.file("extdata", "golden_expected.csv",
                                        package = "doobr"))
  skip_if(nrow(inp) == 0L || nrow(exp_df) == 0L, "fixture not installed")

  fit <- doob_fit(inp$y, horizon = 20,
                  mean_model_fn = lm_model_fn, mlags = c(1, 2),
                  predict_mean_fn = lm_predict_fn,
                  variance_model_fn = lm_model_fn, vlags = 1,
                  predict_variance_fn = lm_predict_fn,
                  transform = "log", difference = TRUE,
                  var_target = "logabs", vlag_source = "y",
                  n_folds = 10, run_md_test = FALSE, eval_folds = FALSE,
                  seed = 42, verbose = FALSE)

  got <- c(c_hat = fit$c_hat, kappa_hat = fit$kappa_hat,
           bank_n = length(fit$bank), bank_mean = mean(fit$bank),
           bank_var = mean(fit$bank^2),
           n_scale_failures = fit$n_scale_failures,
           oof_resid_mean = mean(fit$oof_resid[fit$bank_idx]),
           oof_resid_sd = stats::sd(fit$oof_resid[fit$bank_idx]),
           oof_var_mean = mean(fit$oof_var[fit$bank_idx]),
           oof_var_sd = stats::sd(fit$oof_var[fit$bank_idx]),
           sigma_lo = fit$sigma_range[1], sigma_hi = fit$sigma_range[2])
  qs <- stats::quantile(fit$bank, c(.01, .05, .25, .5, .75, .95, .99))
  names(qs) <- paste0("bank_q", c(1, 5, 25, 50, 75, 95, 99))
  got <- c(got, qs)

  for (q in exp_df$quantity) {
    expect_equal(unname(got[[q]]), exp_df$value[exp_df$quantity == q],
                 tolerance = 1e-8, info = q)
  }
})
