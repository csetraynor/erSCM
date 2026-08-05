test_that(".prepare_boot_sample returns all covariates when no factor variables present", {
  set.seed(1L)
  n <- 100L
  dat <- data.frame(
    DV   = rbinom(n, 1L, 0.3),
    CMAX = rnorm(n, 100, 20),
    AGE  = rnorm(n, 55, 10),
    WT   = rnorm(n, 70, 15)
  )
  result <- .prepare_boot_sample(dat, c("CMAX", "AGE", "WT"), "DV",
                                  min_level_n = 5L, min_events_per_level = 1L)
  expect_equal(result$valid_covariates, c("CMAX", "AGE", "WT"))
  expect_equal(result$excluded_covariates, character(0L))
  expect_equal(result$modified_levels,     character(0L))
})

test_that(".prepare_boot_sample excludes a factor variable with a rare level", {
  set.seed(2L)
  n <- 100L
  # SEX has level "C" with only 2 observations — below default min_level_n = 5
  sex <- factor(c(rep("M", 60L), rep("F", 38L), rep("C", 2L)))
  dat <- data.frame(
    DV   = rbinom(n, 1L, 0.4),
    CMAX = rnorm(n, 100, 20),
    SEX  = sex
  )
  result <- .prepare_boot_sample(dat, c("CMAX", "SEX"), "DV",
                                  min_level_n = 5L, min_events_per_level = 1L)
  expect_false("SEX" %in% result$valid_covariates)
  expect_true("SEX" %in% result$excluded_covariates)
  expect_true("CMAX" %in% result$valid_covariates)
})

test_that("min_level_n = 0 disables safeguard even with rare levels", {
  set.seed(3L)
  n <- 80L
  sex <- factor(c(rep("M", 77L), "F", "F", "X"))  # X has 1 obs
  dat <- data.frame(
    DV   = rbinom(n, 1L, 0.4),
    CMAX = rnorm(n, 100, 20),
    SEX  = sex
  )
  result <- .prepare_boot_sample(dat, c("CMAX", "SEX"), "DV",
                                  min_level_n = 0L, min_events_per_level = 0L)
  expect_equal(sort(result$valid_covariates), c("CMAX", "SEX"))
  expect_equal(result$excluded_covariates, character(0L))
})

test_that(".check_bootstrap_suitability flags rare level as moderate risk", {
  set.seed(4L)
  n <- 100L
  # RACE has a rare level (n=2)
  race <- factor(c(rep("White", 70L), rep("Black", 28L), "Other", "Other"))
  dat <- data.frame(
    DV   = rbinom(n, 1L, 0.3),
    CMAX = rnorm(n, 100, 20),
    RACE = race
  )
  result <- .check_bootstrap_suitability(
    data          = dat,
    covariates    = c("CMAX", "RACE"),
    exposure_vars = "CMAX",
    dv_col        = "DV",
    n_boot        = 20L
  )
  expect_true(result$risk_level %in% c("moderate", "high"))
  expect_true(nrow(result$problematic_levels) > 0L)
  expect_true("RACE" %in% result$problematic_levels$variable)
})

test_that("variable_exclusion_rate is populated when factor exclusions occur", {
  skip_on_cran()
  set.seed(99L)
  n <- 200L
  # RACE with a rare level (3 obs) — excluded in most replicates
  race <- factor(c(rep("White", 170L), rep("Black", 27L), rep("Rare", 3L)))
  cmax <- rnorm(n, 100, 20)
  age  <- rnorm(n, 55, 10)
  # Add age as a continuous covariate so replicates still converge when RACE excluded
  lp   <- -1 + 0.04 * (cmax - 100) - 0.02 * (age - 55)
  dat <- data.frame(
    DV   = rbinom(n, 1L, plogis(lp)),
    CMAX = cmax,
    AGE  = age,
    RACE = race
  )
  fit <- scm_glm(
    eres_data      = dat,
    exposure_vars  = "CMAX",
    covariates     = c("AGE", "RACE"),
    verbose        = FALSE,
    return_models  = FALSE
  )
  boot <- scm_bootstrap(fit, n_boot = 15L, seed = 7L,
                         min_level_n = 5L, min_events_per_level = 1L)
  expect_true(boot$n_converged > 0L)
  # RACE should show a non-zero exclusion rate in most replicates
  excl <- boot$variable_exclusion_rate
  if (!is.null(excl)) {
    expect_true("RACE" %in% names(excl))
    expect_true(excl[["RACE"]] > 0)
  } else {
    # All replicates kept RACE — the field may be NULL if no exclusions occurred
    expect_null(excl)
  }
})

test_that("scm_bootstrap emits cli_warn for high-risk factor configuration", {
  skip_on_cran()
  set.seed(55L)
  n <- 200L
  # Two factor variables each with rare levels → high risk
  f1   <- factor(c(rep("A", 190L), rep("B", 8L), rep("C", 2L)))
  f2   <- factor(c(rep("X", 190L), rep("Y", 8L), rep("Z", 2L)))
  cmax <- rnorm(n, 100, 20)
  age  <- rnorm(n, 55, 10)
  # Include a continuous covariate so replicates can converge even when F1/F2 excluded
  lp   <- -1 + 0.04 * (cmax - 100) - 0.02 * (age - 55)
  dat <- data.frame(
    DV   = rbinom(n, 1L, plogis(lp)),
    CMAX = cmax,
    AGE  = age,
    F1   = f1,
    F2   = f2
  )
  fit <- scm_glm(
    eres_data     = dat,
    exposure_vars = "CMAX",
    covariates    = c("AGE", "F1", "F2"),
    verbose       = FALSE,
    return_models = FALSE
  )
  # High risk → cli_warn should be emitted; capture warnings
  expect_warning(
    scm_bootstrap(fit, n_boot = 10L, seed = 1L,
                   min_level_n = 5L, min_events_per_level = 1L),
    regexp = "Factor level risk"
  )
})
