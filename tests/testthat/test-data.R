## The packaged dataset.
##
## The reference values here are duplicated in doobpy's tests/test_data.py.
## They are a cheap cross-language check that the two packages ship the SAME
## series: if one copy is regenerated and the other is not, both suites fail and
## say so.

N         <- 1004L
LOAD_SUM  <- 428071007.164
TEMP_SUM  <- 14989.709963768115
N_HOLIDAY <- 25L
N_BRIDGE  <- 6L

test_that("itaload has the documented shape", {
  data(itaload, package = "doobr", envir = environment())
  expect_s3_class(itaload, "data.frame")
  expect_identical(nrow(itaload), N)
  expect_identical(names(itaload),
                   c("date", "load", "temperature", "holiday", "bridge",
                     "dow", "month"))
  expect_s3_class(itaload$date, "Date")
  expect_type(itaload$load, "double")
  expect_type(itaload$holiday, "logical")
  expect_type(itaload$bridge, "logical")
  expect_s3_class(itaload$dow, "factor")
  expect_s3_class(itaload$month, "factor")
})

test_that("itaload matches the reference values shared with doobpy", {
  data(itaload, package = "doobr", envir = environment())
  expect_equal(sum(itaload$load), LOAD_SUM, tolerance = 1e-12)
  expect_equal(sum(itaload$temperature), TEMP_SUM, tolerance = 1e-12)
  expect_identical(sum(itaload$holiday), N_HOLIDAY)
  expect_identical(sum(itaload$bridge), N_BRIDGE)
  expect_equal(itaload$load[1], 269120.357)
  expect_equal(itaload$load[N], 459096.5)
})

test_that("the dates are consecutive days and the factors agree with them", {
  data(itaload, package = "doobr", envir = environment())
  expect_identical(format(itaload$date[1]), "2023-01-01")
  expect_identical(format(itaload$date[N]), "2025-09-30")
  expect_true(all(diff(itaload$date) == 1))

  expect_identical(levels(itaload$dow),
                   c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"))
  expect_identical(levels(itaload$month),
                   c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                     "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

  ## %u is the ISO weekday, Monday = 1, which is the factor's level order
  expect_identical(as.integer(itaload$dow),
                   as.integer(format(itaload$date, "%u")))
  expect_identical(as.integer(itaload$month),
                   as.integer(format(itaload$date, "%m")))
})

test_that("flagged holidays never fall at a weekend", {
  ## By construction: the flag marks national holidays that cost a working day.
  data(itaload, package = "doobr", envir = environment())
  expect_false(any(itaload$holiday & itaload$dow %in% c("Sat", "Sun")))
})

test_that("itaload runs end to end at minimal settings", {
  data(itaload, package = "doobr", envir = environment())
  X <- cbind(stats::model.matrix(~ dow, itaload)[, -1, drop = FALSE],
             hol = as.numeric(itaload$holiday))
  fit <- doob_fit(itaload$load, horizon = 7L, xreg = X,
                  mean_model_fn = lm_model_fn, mlags = c(1L, 7L),
                  predict_mean_fn = lm_predict_fn,
                  n_folds = 10L, nsim_eval = 20L,
                  md_type = "additive", md_stat = "F",
                  seed = 1, verbose = FALSE)
  expect_identical(length(fit$bank), 70L)
  expect_equal(mean(fit$bank), 0, tolerance = 1e-12)
  expect_equal(mean(fit$bank^2), 1, tolerance = 1e-12)

  paths <- simulate(fit, nsim = 50L, seed = 1, newxreg = X[1:7, , drop = FALSE])
  expect_identical(dim(paths$paths), c(7L, 50L))
  expect_true(all(paths$paths > 0))
})
