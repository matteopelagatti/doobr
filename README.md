# doobr

Nonparametric simulation of future sample paths and scenario trees for a
univariate time series, without assuming a parametric data generating process.

Any learner may supply the conditional mean and the conditional scale. The
innovations are resampled from an honest bank of **out-of-fold** standardised
residuals obtained by purged K-fold cross-fitting, so the simulated increments
form a martingale difference sequence rather than inheriting the optimism of
in-sample residuals — which is what makes naive residual resampling inject
artificial predictability into generated scenarios.

Companion Python package: [doobpy](https://github.com/matteopelagatti/doobpy).

## Install

```r
# install.packages("remotes")
remotes::install_github("matteopelagatti/doobr", build_vignettes = TRUE)
```

## Use

```r
library(doobr)

fit <- doob_fit(
  y, horizon = 12,
  mean_model_fn  = lm_model_fn, mlags = c(1, 12, 13),
  predict_mean_fn = lm_predict_fn,
  variance_model_fn = lm_model_fn, vlags = c(1, 12),
  predict_variance_fn = lm_predict_fn,
  transform = "log", difference = TRUE,
  var_target = "logabs"
)

fit                                   # fit, bank, MDH p-values, mean scores
summary(fit)                          # adds the calibration tests

paths <- simulate(fit, nsim = 500)    # 12 x 500, on the original scale
tree  <- simulate(fit, output = "tree", branching = c(10, 5, rep(2, 10)))
```

Any learner works by supplying a matched pair:

```r
model_fn(y, y_lags, xreg) -> model
predict_fn(model, lags_vec, xreg_row) -> scalar
```

Adapters for linear models, random forests, SVMs and XGBoost are included.

## What it checks

- **Martingale-difference tests on both estimation steps** — on the residuals,
  and on the squared residuals net of the fitted conditional variance. A
  pipeline that passes the first and fails the second produces paths with the
  right drift and the wrong spread.
- **Calibration of the simulated paths** — rank histogram, permutation tests on
  the energy and variogram scores, and a band-depth test, combined across folds
  by Fisher's method.

All diagnostics run on out-of-fold quantities, so none is inflated by in-sample
optimism.

## Data

`itaload` is Italian daily electricity load from 2023-01-01 to 2025-09-30 (1004
days) with daily mean temperature and the Italian calendar — national holidays,
bridge days, weekday and month. It is the series used in the applied example.

```r
data(itaload)
plot(itaload$date, itaload$load / 1e3, type = "l", ylab = "GWh")
```

## Documentation

```r
vignette("doobr")        # the API, and the two details that are easy to get wrong
vignette("electricity")  # a worked application on itaload
```

The electricity vignette is worth reading even if you only want the API, because
it shows the two things the diagnostics are for. A calendar specification that
improves the residual standard deviation by 27% and the energy score by 14% is
caught by the step-1 test as a gross violation of the martingale property
(p from 2.6e-03 to 9.5e-32), and the fix improves both. Then the three
conditional-scale targets, which the energy score separates by 1%, are separated
by the step-2 test by eighteen orders of magnitude.

`DESIGN.md` specifies the algorithm (shared with doobpy); `DECISIONS.md` records
the non-obvious choices and the failures that prompted them.

## Licence

GPL-3.
