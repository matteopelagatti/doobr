# doobr — working notes for Claude

R package for nonparametric simulation of future sample paths and scenario
trees. Companion Python package: **doobpy** (separate repo, same algorithm,
same golden fixture). The method is described in `DESIGN.md`; the non-obvious
choices and why they were made are in `DECISIONS.md`. **Read both before
changing anything in `R/doob_fit.R` or `R/tree.R`.**

## Repositories

- **doobr** (this package): <https://github.com/matteopelagatti/doobr>
- **doobpy** (companion): <https://github.com/matteopelagatti/doobpy>

The GitHub account is **matteopelagatti**, not `mpelagatti`. An earlier draft of
`DESCRIPTION`, `README.md` and doobpy's `pyproject.toml` had the wrong name and
shipped dead `URL:`/`BugReports:` links; if you see `mpelagatti` anywhere, it is
a bug.

## Environment setup

R is not preinstalled in the cloud container, and **CRAN is blocked** by the
egress proxy. Install R and the toolchain from the Ubuntu archive instead:

```bash
apt-get update -qq
apt-get install -y -qq r-base-core r-cran-testthat r-cran-roxygen2 \
                      r-cran-randomforest r-cran-e1071
```

Takes about two minutes. `xgboost` is not available this way, so the
`xgb_*` learners cannot be exercised here; `R CMD check` reports them as a NOTE
("Package suggested but not available"), which is expected and not a defect.

The user's own machine, reachable through `device_bash`, runs an isolated Linux
VM with **no R**. Do not try to run R there. Edit files there if the repo is a
connected folder; run R here.

## Build, document, test

```bash
cd doobr
Rscript -e 'roxygen2::roxygenise(".")'       # after changing any roxygen block
R CMD INSTALL --no-multiarch .
Rscript -e 'library(testthat); library(doobr); test_dir("tests/testthat", reporter="summary")'

cd .. && R CMD build doobr
LC_ALL=C.UTF-8 _R_CHECK_FORCE_SUGGESTS_=0 R CMD check --no-manual doobr_*.tar.gz
```

Do **not** use `--as-cran` here: it needs CRAN for the incoming checks and will
fail with a connection error that has nothing to do with the package. Set
`LC_ALL=C.UTF-8`, otherwise the check emits a spurious locale WARNING.

Expected clean result: `Status: 1 NOTE` (the xgboost suggestion).

## Architecture

| File | Contents |
|---|---|
| `R/doob_fit.R` | `doob_fit()` — cross-fitting, the bank, the scale constant, diagnostics |
| `R/simulate.R` | `run_recursion()` (shared engine), `simulate.doob_fit()` |
| `R/tree.R` | `quantise_bank()`, `doob_tree()` |
| `R/evaluate.R` | `evaluate_folds()` — per-fold simulation for the diagnostics |
| `R/md_test.R` | `md_test()` — martingale-difference test |
| `R/calibration.R` | `calibration_tests()` — rank, permutation, band depth, Fisher |
| `R/scores.R` | `es_score()`, `vs_score()` |
| `R/learners.R` | model/predict adapter pairs |
| `R/methods.R` | print / summary / plot |
| `R/utils.R` | transforms, `robust_sd()`, lag matrices, `fisher_combine()` |

`doob_fit()` returns everything `simulate()` needs in `$spec`. The two are
deliberately separate: fitting costs `2K+2` model fits plus the diagnostics,
after which paths and trees are cheap at any size.

## Invariants — break these and the method is wrong

These are enforced by `tests/testthat/` and must stay true:

1. **The bank has mean exactly 0 and variance exactly 1.**
2. **The scale constant `c_hat` is carried into the conditional variance, never
   discarded.** Normalising the bank alone destroys the self-correction against
   a multiplicative bias in `v_hat`. See `DECISIONS.md` §1 — this was a real bug
   that inflated the simulated scale 75-fold on one series.
3. **States where the fitted scale leaves the admissible range are excluded from
   the bank, not clipped into it.** They are states where the scale model
   failed; their standardised residual is meaningless. This removes broken
   variance estimates, *not* large residuals — the innovation tails must survive.
4. **At every tree node, `sum(q*e) == 0` and `sum(q*e^2) == 1` exactly**, so the
   conditional mean on the tree equals the model forecast to machine precision.
5. **The martingale-difference tests run on the modelling scale**; the
   calibration tests and the returned paths are on the original scale.

## Traps found by running on real data

- **`n_folds` sets the step-2 sample size.** The conditional-scale model trains
  on out-of-fold residuals, so its training set is the bank
  (`n_folds * horizon`), not the series. `n_folds = 8` with `horizon = 7` on
  `itaload` leaves 49 observations for a 13-parameter scale model and drives
  `c_hat` to 7-9. Keep `n_folds = NULL` unless there is a reason. See
  `DECISIONS.md` section 10.
- **A specification can improve every score while breaking the martingale
  property.** On `itaload`, adding same-day holiday dummies cuts the residual sd
  27% and the energy score 14% while taking the step-1 p-value from 2.6e-03 to
  9.5e-32 (residual lag-1 acf 0.004 -> 0.399), because Italian holiday effects
  spill onto adjacent days. Allowing the spillover restores it. This is the
  worked example in `vignettes/electricity.Rmd`; do not "simplify" that
  three-calendar comparison away, it is the point of the vignette.
- `itaload` leaves a lag-7 residual autocorrelation of ~0.23 that annual
  harmonics, dow-by-season interactions, a trend, lags to 28, and heating and
  cooling degree days all fail to remove. Unexplained; stated as such in the
  vignette. Do not claim it is fixed.

## Vignettes

Two, both run live at build time: `vignettes/doobr.Rmd` (~13 s) and
`vignettes/electricity.Rmd` (~70 s). `R CMD check` re-builds them, so keep
`nsim_eval` at 100, the horizon at 7, and the trees small (`prod(branching)`
leaves; the node loop builds one `data.frame` per node). `xgboost` is
unavailable here, so any xgboost chunk must be `eval = FALSE` or guarded by
`requireNamespace()`.

## Conventions

- `xreg` has **one row per observation of `y`**; the future rows go to
  `simulate(newxreg=)`. Do not reintroduce the old pre-concatenated layout.
- Learner adapters keep the **single-row predict contract**:
  `predict_fn(model, lags_vec, xreg_row)` where `lags_vec` is a *named* vector.
  Do not vectorise it without providing a compatibility path — users have their
  own adapters.
- Stochastic learners take a `seed` argument. An unseeded random forest gives
  different diagnostics every run; that bit us.
- Keep `Imports` to base R. CRAN is unreachable here, so a new hard dependency
  cannot be tested.

## Cross-language parity

`inst/extdata/golden_input.csv` → `inst/extdata/golden_expected.csv`, checked by
`tests/testthat/test-golden.R`. **doobpy must reproduce the same numbers from
the same input to 1e-8.** Regenerate only with `data-raw/make_golden.R`, only
deliberately, and record why in `DECISIONS.md`. A silent regeneration defeats
the purpose of the fixture.
