# The algorithm

This document is the shared specification for **doobr** (R) and **doobpy**
(Python). It describes what the two packages must do; it is not an API
reference. Where the two languages differ in style they must not differ in
numbers — see §8.

## 1. Setting

Let $X_t$ be a univariate series and $W_t$ deterministic regressors known over
the whole simulation span. Write $w$ for the modelling series: $X$ itself, or a
variance-stabilising transform of it, differenced if requested. The state is

$$Z_{t-1} = (w_{t-1}, \dots, w_{t-p},\ W_t),$$

and the model is location–scale,

$$w_t = \mu_0(Z_{t-1}) + \sigma_0(Z_{t-1})\,\eta_t,
\qquad \eta_t \ \text{i.i.d.},\ \mathbb E\eta=0,\ \mathbb E\eta^2=1,$$

with $\eta_t$ independent of the past. The Doob decomposition splits $w$ into a
predictable part, which any learner may estimate, and a martingale, whose
increments are what scenario generation must sample.

**The state carries no lagged squared innovations.** It does not need them:
$\varepsilon_{t-j} = w_{t-j} - \mu_0(Z_{t-j-1})$ is a function of lagged
observations, so an ARCH($r$) conditional variance is exactly representable with
$p = p_\mu + r$ lags. What this excludes is genuinely infinite memory (GARCH),
at a cost of $O(\beta^p)$. The payoff is that $Z_{t-1}$ is fully observed, so
the second estimation step has a generated *regressand* only, never a generated
*regressor*.

**The location–scale assumption is not cosmetic.** A pooled innovation bank is
valid only if the conditional law has a shape that does not depend on the state.
A martingale difference with unit conditional variance may still have
state-dependent skewness or tail index. The martingale-difference property plus
a correct conditional variance is *necessary but not sufficient*.

## 2. Cross-fitting

Partition the usable indices into $K$ contiguous folds of length $h$, the
intended simulation horizon. When fitting for fold $k$, delete the $M$
observations on each side of it (*purging*); with lagged values in the state,
omitting the fold alone does not break the dependence between training and
evaluation samples.

Fold length $=$ horizon is what makes one partition serve three purposes: each
fold yields $h$ out-of-fold residuals for the bank, is an evaluation window of
exactly the length the method simulates, and contributes to the
martingale-difference tests — every pooled residual being out-of-sample.

Two passes:

1. Fit $\hat\mu^{(-k)}$, predict on $B_k$, store $\hat\varepsilon_t$.
2. Fit the conditional-scale learner on the out-of-fold responses outside $B_k$,
   predict on $B_k$, store the raw scale.

Then refit both on the full sample for forward simulation — the scale model on
the **out-of-fold** responses, which target $\sigma^2 + \mathrm{MSE}$ rather
than the deflated in-sample analogue.

## 3. The conditional-scale step

Three losses, none dominant.

| | V1 squared | V2 log-abs | V3 Gaussian QML |
|---|---|---|---|
| response / loss | $(\hat\varepsilon^2 - g)^2$ | $(\log\max(\lvert\hat\varepsilon\rvert,\tau) - g)^2$ | $g + \hat\varepsilon^2 e^{-g}$ |
| link | $\sqrt{g}$ | $\hat c\,e^{g}$ | $e^{g/2}$ |
| target | $\sigma_0^2$ | $\log\sigma_0 + \kappa$ | $\log\sigma_0^2$ |
| moment condition | $\mathbb E\eta^4<\infty$ | $\operatorname{Var}\log\lvert\eta\rvert<\infty$ | $\mathbb E\eta^4<\infty$ |
| positivity | clamp | automatic | automatic |
| two-step contamination | $O(a_n^2)$ | $O(a_n^{2/3})$ | $O(a_n^2)$ |
| at $\varepsilon\approx 0$ | harmless | fatal (trim) | graceful |
| any LS learner | yes | yes | no |

V1 and V2 are implemented. **V3 is not yet implemented** and is the natural next
addition: it keeps V1's target and robustness while gaining positivity and a
homoskedastic score, at the cost of a custom objective.

For V2, least squares estimates $\log\sigma_0 + \kappa$ with
$\kappa = \mathbb E\log\lvert\eta\rvert \approx -0.635$ for Gaussian
innovations, so exponentiating the fit is biased by $e^{\kappa}\approx 0.53$.
The correction is §4 and requires no knowledge of $F_0$.

## 4. The bank and the scale constant

With $v^{\mathrm{raw}}$ the raw fitted variance,

$$z^{\mathrm{raw}}_t = \frac{\hat\varepsilon_t}{\sqrt{v^{\mathrm{raw}}_t}},
\qquad
\hat c = \operatorname{sd}(z^{\mathrm{raw}}),
\qquad
z_t = \frac{z^{\mathrm{raw}}_t - \bar z}{\hat c},
\qquad
\hat v_t = \hat c^{\,2} v^{\mathrm{raw}}_t .$$

The bank has mean exactly 0 and variance exactly 1, so the simulated paths
satisfy the conditional-mean and conditional-variance identities exactly. **The
constant must be carried into $\hat v$, not discarded** — normalising the bank
alone would destroy the self-correction a raw bank gives against a multiplicative
bias in $\hat v$. For V2 this *is* the $e^{\kappa}$ correction, with
$\hat c = e^{-\hat\kappa}$ estimated implicitly at the parametric rate.

Admissible range for the scale: $[s_{\mathrm{lo}}, s_{\mathrm{hi}}]$ with
$s_{\mathrm{lo}} = 0.05\cdot\operatorname{MADsd}(\hat\varepsilon)$ and
$s_{\mathrm{hi}} = 3\max(q_{0.995}(\hat\sigma^{\mathrm{raw}}),
\operatorname{MADsd}(\hat\varepsilon))$. The floor is anchored to the **residual**
scale, not the fitted scale, which is useless exactly when the fitted values are
pathological. States outside the range are **excluded from the bank**, not
clipped into it: $\hat c$ is a root-mean-square and a handful of states with a
near-zero $\hat v$ would otherwise dominate it. This drops broken variance
estimates, not large residuals.

## 5. Simulation

$$X^{(s)}_{n+\tau} = \hat\mu(Z^{(s)}_{n+\tau-1})
  + \hat c\,\operatorname{clip}\!\big(\hat\sigma^{\mathrm{raw}}(Z^{(s)}_{n+\tau-1})\big)\, z^{(s)}_\tau,
  \qquad z^{(s)}_\tau \sim \text{bank}.$$

Per-fold evaluation uses the fold-$k$ models and a **leave-one-fold-out bank
with its own constant** — bank and constant always move together.

## 6. Trees

Branch on the martingale component, so the tree inherits the decomposition node
by node. At each node the bank is quantised into $b$ atoms: within-stratum means
of $b$ equiprobable strata (the stratified conditional-mean quantiser), then
recentred and rescaled so that

$$\sum_c q_c e_c = 0, \qquad \sum_c q_c e_c^2 = 1$$

**exactly**. Both corrections are $O(n^{-1/2})$ and buy exact finite-sample
conditions; hence $\mathbb E_q[X_c \mid \mathcal F_n] = \hat\mu(\mathcal F_n)$ to
machine precision, which the test suite checks.

## 7. Transformations

The series is supplied on its original scale; the package transforms, models,
simulates, and maps back. **No correction is involved and none is needed:**
simulation is a statement about the conditional law, and laws transform by
change of variables, so exponentiating a correct log-scale simulation gives a
correct simulation of the level.

Two consequences, neither a defect. First, $\exp(\hat\mu)$ is the conditional
*median*, not the mean; the mean is $\exp(\hat\mu)\varsigma$ with $\varsigma$ the
Duan (1983) smearing factor, which the ensemble already reproduces and which is
reported per step. Second, the rank histogram is invariant to a monotone
transformation while the energy, variogram and band-depth statistics are not, so
the scale on which the calibration tests are reported must be stated.

## 8. Cross-language parity

`golden_input.csv` → `golden_expected.csv`, reproduced by both packages to
**1e-8** absolute/relative. Quantities compared: `c_hat`, `kappa_hat`, bank size,
mean and variance, scale-failure count, out-of-fold residual and variance
moments, the admissible range, and seven bank quantiles.

Deterministic given the input: everything above. Not deterministic across
languages: anything depending on the RNG (simulated paths, wild-bootstrap
p-values). Parity is therefore asserted on the *fit*, not on simulated draws.

## 9. Asymptotics, in one paragraph

With $a_n, b_n$ the $L^2$ rates of the conditional mean and scale learners,
$\mathcal W_1(\hat F_n, F_0) = O_P(a_n + b_n + n^{-1/2})$, and the one-step
kernel error is of the same order. Under a stability condition
$\sum_j \ell_j < 1$ the $h$-step error is bounded by $\varepsilon_n/(1-L)$
**uniformly in $h$**; without it, the admissible horizon grows as
$h_n = o(\varepsilon_n^{-1})$. Because the coupling is built one step at a time
it is bicausal, so the same bound holds in the nested distance, hence optimal
values of the associated multistage stochastic program converge. Cross-fitting
is what lets the learner enter through its excess risk alone, with no entropy
condition on the hypothesis class.
