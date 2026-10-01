## Generates the cross-language reference fixture.
##
## Both doobr and doobpy must reproduce inst/extdata/golden_expected.csv from
## inst/extdata/golden_input.csv to the tolerances stated in DESIGN.md.  Run
## this ONLY when the reference is deliberately changed, and record why in
## DECISIONS.md -- a silent regeneration defeats the purpose.
library(doobr)

set.seed(20240101)
n <- 600
w <- numeric(n); d <- 0
for (t in 2:n) {
  s <- sqrt(0.02^2 + 0.15^2 * d^2)
  d <- 0.4 * d + s * rnorm(1)
  w[t] <- w[t - 1] + d
}
y <- exp(w + 5)
utils::write.csv(data.frame(t = seq_len(n), y = y),
                 "inst/extdata/golden_input.csv", row.names = FALSE)

fit <- doob_fit(y, horizon = 20,
                mean_model_fn = lm_model_fn, mlags = c(1, 2),
                predict_mean_fn = lm_predict_fn,
                variance_model_fn = lm_model_fn, vlags = 1,
                predict_variance_fn = lm_predict_fn,
                transform = "log", difference = TRUE,
                var_target = "logabs", vlag_source = "y",
                n_folds = 10, run_md_test = FALSE, eval_folds = FALSE,
                seed = 42, verbose = FALSE)

out <- data.frame(
  quantity = c("c_hat", "kappa_hat", "bank_n", "bank_mean", "bank_var",
               "n_scale_failures", "oof_resid_mean", "oof_resid_sd",
               "oof_var_mean", "oof_var_sd", "sigma_lo", "sigma_hi",
               paste0("bank_q", c(1, 5, 25, 50, 75, 95, 99))),
  value = c(fit$c_hat, fit$kappa_hat, length(fit$bank),
            mean(fit$bank), mean(fit$bank^2), fit$n_scale_failures,
            mean(fit$oof_resid[fit$bank_idx]), sd(fit$oof_resid[fit$bank_idx]),
            mean(fit$oof_var[fit$bank_idx]), sd(fit$oof_var[fit$bank_idx]),
            fit$sigma_range[1], fit$sigma_range[2],
            unname(quantile(fit$bank, c(.01, .05, .25, .5, .75, .95, .99)))))
utils::write.csv(out, "inst/extdata/golden_expected.csv", row.names = FALSE)
cat("golden fixture written\n"); print(out, digits = 10)
