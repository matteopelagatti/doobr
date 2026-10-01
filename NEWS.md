# doobr 0.1.0

First release.

* `doob_fit()` fits the conditional mean and scale by purged K-fold
  cross-fitting, with fold length equal to the simulation horizon, and returns
  an honest innovation bank together with out-of-sample diagnostics.
* `simulate()` generates either a flat ensemble of sample paths or a scenario
  tree, from the same fitted object.
* `md_test()` tests the martingale-difference hypothesis; `doob_fit()` applies
  it to both estimation steps.
* `calibration_tests()` assesses whether the realised path is distinguishable
  from the simulated ensemble.
* Variance targets `"squared"` and `"logabs"`; transformations `"log"` and
  `"sqrt"` handled internally, with the smearing factor reported.
* Learner adapters for linear models, random forests, SVMs and XGBoost; any
  learner can be supplied as a `model_fn` / `predict_fn` pair.
* Two vignettes: `vignette("doobr")` for the API and the two details that are
  easy to get wrong, and `vignette("electricity")` for a worked application.
* `itaload`: Italian daily electricity load with temperature and the Italian
  calendar, 2023-01-01 to 2025-09-30, as used in the applied example.
