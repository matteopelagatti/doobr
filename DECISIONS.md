# Decisions, and why

Non-obvious choices, each with the reasoning and — where there was one — the
failure that prompted it. Several of these were found by running the code on
real data; they are not design preferences and should not be reverted without
reproducing the failure first.

---

## 1. The scale constant is carried into `v_hat`, not discarded

**Decision.** `c_hat = sd(z_raw)`; the bank is `z_raw / c_hat` and the
conditional variance is `v_raw * c_hat^2`.

**Why.** Normalising the bank to unit variance *destroys* a self-correction. If
`v_hat` is $c^2$ times the truth, the raw standardised residuals have variance
$c^{-2}$ and the product $\sqrt{v}\,z$ is correctly scaled; after normalising the
bank, the entire scale rests on `v_hat`. The constant must move from one to the
other, not vanish.

**Failure that prompted it.** First implemented with the bank normalised and
`v_hat` untouched. The simulated ensembles were systematically over-dispersed by
the factor by which `v_hat` was biased. Caught on a synthetic design where the
true variance was known.

For `var_target = "logabs"` this correction *is* the $e^{\kappa}$ factor,
$\kappa = \mathbb E\log|\eta| \approx -0.635$ for Gaussian innovations. Omitting
it understates the conditional standard deviation by about a factor of two —
while every diagnostic computed on the normalised bank still passes. On the
airline series `c_hat = 2.15`, so the conditional variance would be wrong by
$4.6\times$.

---

## 2. Scale failures are excluded from the bank, not clipped into it

**Decision.** States where the fitted scale leaves
$[s_{\mathrm{lo}}, s_{\mathrm{hi}}]$ are dropped from `bank_idx`.

**Why.** `c_hat` is a root-mean-square and is not robust. A state where `v_hat`
is near zero against an ordinary residual produces an enormous standardised
residual; a handful of them dominate `c_hat` and inflate the scale of every
simulated path.

**Failure that prompted it.** On Italian daily electricity load, a least-squares
fit to squared residuals went out of range at 14% of states. Clipping them into
the bank gave `c_hat = 8.67` — the simulated scale inflated 75-fold — and the
mean energy score rose from 230,247 to 753,751. Excluding them gives
`c_hat = 1.06` and ES 227,590.

**Rejected alternative: winsorising `z_raw`.** Tried, and wrong: it trims the
*tails of the innovation distribution*, which is exactly what scenario
generation must reproduce. The contamination comes from broken variance
estimates, not from large residuals, so the fix must target the former. A
constant-variance sanity case gave `c_hat = 0.9637` under winsorising instead of
the correct `1.000`, which is how the error was found.

**Consequence to preserve.** The floor is anchored to a robust estimate of the
**residual** scale, not to quantiles of the fitted scale — anchoring it to the
fitted values is useless precisely when those are pathological.

---

## 3. The state excludes lagged squared innovations (`vlag_source = "y"`)

**Decision.** Default to indexing the conditional-scale model by lags of the
modelling series. `"sq_res"` is retained for comparison.

**Why.** With lagged squared innovations in the state, the design of the second
regression is itself estimated, so the two-step nesting has a generated
*regressor* as well as a generated *regressand*. The regressor error is
first-order, $O_P(a_n)$, and `sigma_hat` stops being adaptive to a smooth
`sigma_0`. With lags of $y$ the state is fully observed and the contamination is
second-order, $O_P(a_n^2)$. Nothing is lost for ARCH: $\varepsilon_{t-j}$ is
itself a function of lagged observations.

**Cost, accepted knowingly.** Genuinely infinite memory (GARCH) is no longer
representable; truncating at $p$ lags leaves $O(\beta^p)$. Without structural
restrictions the optimal $p$ grows like $\sqrt{\log n}$ and the rate degrades to
$\exp(-c\sqrt{\log n})$ — consistent but sub-polynomial. With additive or
compositional structure the low-dimensional rate is restored.

---

## 4. The transformation is handled inside the package

**Decision.** `transform = "none" | "log" | "sqrt"`; the series is supplied on
its original scale.

**Why.** Not because the back-transformation needs correcting — it does not.
Simulation is a statement about the conditional law, and laws transform by
change of variables, so exponentiating a correct log-scale simulation is exactly
right. Verified against a known DGP: mean ratio 1.0003, median 1.0007, 1% and
99% quantiles within 0.3%.

The reasons are practical: it validates positivity before taking logs, it
returns paths on the scale the user judges them on, and it makes the
**diagnostics scale explicit**. The last is a real issue that was otherwise
invisible — the rank histogram is invariant to a monotone transformation but the
energy, variogram and band-depth statistics are not, so those three were
previously being reported in log units with no label.

**A claim that was wrong and is recorded here so it is not reintroduced.** An
earlier draft asserted that the log-scale martingale property "biases" the level
scale and needs correcting. It does not. If $\log X$ is a martingale then $X$ is
genuinely a submartingale, and the simulator reproducing that is correct
behaviour. What survives is a *labelling* point — $\exp(\hat\mu)$ is the
conditional median, not the mean — and the `smearing` output exists to serve
that, as a descriptive quantity. The one case where the scale genuinely matters
is a hard martingale constraint on the original scale (risk-neutral pricing),
which is a modelling choice, not a correction.

---

## 5. Fold length equals the simulation horizon

**Decision.** `fold_length` defaults to `horizon`.

**Why.** One partition then serves three purposes at once: honest innovations,
evaluation windows of the right length, and out-of-sample residuals for the
martingale tests. Purging by the model memory is **not optional** — with lagged
values in the state, omitting the fold alone does not break the dependence.

**Consequence.** Separating fit from simulate means someone will fit with
28-step folds and simulate 90 steps. Legitimate for the paths, but the
calibration diagnostics only validate the horizon they ran at, so
`simulate()` warns when asked for a different one.

---

## 6. `xreg` has one row per observation of `y`

**Decision.** Historical regressors go to `doob_fit(xreg=)` aligned with `y`;
future regressors go to `simulate(newxreg=)`.

**Why.** The predecessor required a single pre-concatenated matrix covering
history *and* future, with an off-by-one when differencing that was easy to get
wrong silently. Splitting them makes the contract checkable, and both lengths
are validated with explicit error messages.

---

## 7. Learner adapters keep the single-row predict contract

**Decision.** `predict_fn(model, lags_vec, xreg_row)` with a *named* vector.

**Why.** Users have their own adapters written against this contract.
Vectorising would be faster but breaks them. If vectorisation is added later it
must be opt-in with a compatibility path.

Stochastic learners take a `seed`. An unseeded random forest gives different
diagnostics on every run — the step-2 p-value ranged over 0.042–0.170 across 20
forest seeds on the airline series, straddling the conventional threshold, while
the energy score was stable to three decimals.

---

## 8. `c_hat` is not exactly 1 under constant variance

**Observed, not decided.** With `variance_model_fn = NULL`, `c_hat` comes out as
$\sqrt{(n-1)/n}$, because the constant uses `var()` (denominator $n-1$) while the
bank variance is a plain mean of squares. Their product is the $n$-denominator
variance of the out-of-fold residuals, which is self-consistent: the bank has
variance exactly 1 and the scale matches it.

Recorded because it looks like a bug and is not. Changing `var()` to the
$n$-denominator form in that branch would make `c_hat` exactly 1 with no other
effect; left as is because the current form is internally consistent.

---

## 9. The constant variance comes from out-of-fold residuals

**Decision.** `const_v = var(oof_resid)`, not the in-sample analogue.

**Why.** The homoskedastic branch was never exempt from the in-sample optimism
the whole method exists to correct — it was just invisible there. On the airline
series the out-of-fold variance is 2.9% larger, and the simulated dispersion
3.1% wider.

---

## 10. The number of folds sets the step-2 sample size

**Observed, and a trap.** The conditional-scale model is trained on *out-of-fold*
residuals, so its training sample is the innovation bank — `n_folds * horizon`
observations — not the length of the series. Reducing `n_folds` to save time
silently starves step 2.

**How it shows up.** On the packaged `itaload` series with `horizon = 7`,
`n_folds = 8` leaves 49 observations to fit a 13-parameter conditional-scale
model. The fitted scales are wild, `c_hat` comes out between 7 and 9 (against
about 2.1 with the default exhaustive partition), and the simulated ensembles are
over-dispersed by the same factor. The symptom looks exactly like a bug in the
scale constant and is not one.

**Consequence.** The default `n_folds = NULL` (partition the series
exhaustively) is the right default and should stay. Anything that reduces it is
a step-2 sample-size decision, and the documentation says so. If `c_hat` is far
from 1 for `var_target = "squared"` or far from `e^{-\kappa}` for `"logabs"`,
check `length(fit$bank)` before suspecting the estimator.

---

## 11. Pending

- **V3 (Gaussian QML) is specified in `DESIGN.md` but not implemented.** It
  needs a custom objective, so it cannot use a stock least-squares learner; the
  adapter interface would need a loss argument.
- **The pooled rank test is mildly liberal.** Under a correctly specified model
  it rejects 12–14% of the time at the nominal 5%, because ranks at adjacent
  steps within a fold are treated as independent when they are positively
  dependent. Either document the baseline or use one rank per fold.
- **Multivariate extension**, which is where stochastic regressors belong.
