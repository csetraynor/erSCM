# Tests for Phase 5: start_covariates and direction parameters
# in scm_glm() and scm_coxph().

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

.make_glm_data <- function(n = 300, seed = 42) {
  set.seed(seed)
  x   <- rnorm(n, 100, 20)
  age <- rnorm(n, 55, 12)
  wt  <- rnorm(n, 75, 15)
  lp  <- -1 + 0.05 * (x - 100) - 0.02 * (age - 55)
  data.frame(
    DV   = rbinom(n, 1L, plogis(lp)),
    CMAX = x,
    AGE  = age,
    WT   = wt
  )
}

.make_coxph_data <- function(n = 300, seed = 42) {
  set.seed(seed)
  x   <- rnorm(n, 100, 20)
  age <- rnorm(n, 55, 12)
  wt  <- rnorm(n, 75, 15)
  lp  <- 0.04 * (x - 100) - 0.02 * (age - 55)
  data.frame(
    TIME = rexp(n, rate = exp(lp) * 0.05),
    DV   = rbinom(n, 1L, 0.6),
    CMAX = x,
    AGE  = age,
    WT   = wt
  )
}

# ---------------------------------------------------------------------------
# Validation: start_covariates
# ---------------------------------------------------------------------------

test_that("start_covariates not in covariates raises error", {
  dat <- .make_glm_data()
  expect_error(
    scm_glm(
      eres_data          = dat,
      exposure_vars      = "CMAX",
      covariates         = c("AGE", "WT"),
      significance_level = c(0.05, 0.01),
      dv_col             = "DV",
      family             = "binomial",
      start_covariates   = "NONEXISTENT",
      center = FALSE, scale_to_median = FALSE, verbose = FALSE
    ),
    "start_covariates.*not in.*covariates"
  )
})

test_that("start_covariates overlapping adjust_for raises error", {
  dat <- .make_glm_data()
  expect_error(
    scm_glm(
      eres_data          = dat,
      exposure_vars      = "CMAX",
      covariates         = c("AGE", "WT"),
      significance_level = c(0.05, 0.01),
      dv_col             = "DV",
      family             = "binomial",
      adjust_for         = "AGE",
      start_covariates   = "AGE",
      center = FALSE, scale_to_median = FALSE, verbose = FALSE
    ),
    "start_covariates.*adjust_for"
  )
})

test_that("direction='backward' without start_covariates raises error", {
  dat <- .make_glm_data()
  expect_error(
    scm_glm(
      eres_data          = dat,
      exposure_vars      = "CMAX",
      covariates         = c("AGE", "WT"),
      significance_level = c(0.05, 0.01),
      dv_col             = "DV",
      family             = "binomial",
      direction          = "backward",
      center = FALSE, scale_to_median = FALSE, verbose = FALSE
    ),
    "direction.*backward.*start_covariates"
  )
})

test_that("invalid direction raises error", {
  dat <- .make_glm_data()
  expect_error(
    scm_glm(
      eres_data          = dat,
      exposure_vars      = "CMAX",
      covariates         = c("AGE", "WT"),
      significance_level = c(0.05, 0.01),
      dv_col             = "DV",
      family             = "binomial",
      direction          = "sideways",
      center = FALSE, scale_to_median = FALSE, verbose = FALSE
    )
  )
})

# ---------------------------------------------------------------------------
# start_covariates: warm start
# ---------------------------------------------------------------------------

test_that("start_covariates enters model before forward selection (GLM)", {
  skip_on_cran()
  dat <- .make_glm_data()
  fit <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT"),
    significance_level = c(0.05, 0.01),
    dv_col             = "DV",
    family             = "binomial",
    start_covariates   = "AGE",
    center = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  expect_s3_class(fit, "eres_scm_glm")
  expect_equal(fit$metadata$start_covariates, "AGE")
})

test_that("start_covariates stored in metadata", {
  skip_on_cran()
  dat <- .make_glm_data()
  fit <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT"),
    significance_level = c(0.05, 0.01),
    dv_col             = "DV",
    family             = "binomial",
    start_covariates   = c("AGE", "WT"),
    center = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  expect_equal(fit$metadata$start_covariates, c("AGE", "WT"))
})

test_that("start_covariates are not re-added in forward step", {
  skip_on_cran()
  dat <- .make_glm_data()
  # With AGE as start covariate and direction=forward (no backward removal),
  # AGE should appear exactly once — not duplicated by the forward addition loop.
  fit <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT"),
    significance_level = c(0.05, 0.01),
    dv_col             = "DV",
    family             = "binomial",
    start_covariates   = "AGE",
    direction          = "forward",
    center = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  # AGE should appear exactly once in finalRelation, not duplicated
  age_count <- sum(fit$metadata$finalRelation == "AGE")
  expect_equal(age_count, 1L)
})

# ---------------------------------------------------------------------------
# direction = "forward" (no backward step)
# ---------------------------------------------------------------------------

test_that("direction='forward' returns eres_scm_glm", {
  skip_on_cran()
  dat <- .make_glm_data()
  fit <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT"),
    significance_level = c(0.3, 0.1),
    dv_col             = "DV",
    family             = "binomial",
    direction          = "forward",
    center = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  expect_s3_class(fit, "eres_scm_glm")
  expect_equal(fit$metadata$direction, "forward")
})

test_that("direction='forward' stores direction in metadata", {
  skip_on_cran()
  dat <- .make_glm_data()
  fit <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT"),
    significance_level = c(0.3, 0.1),
    dv_col             = "DV",
    family             = "binomial",
    direction          = "forward",
    center = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  expect_equal(fit$metadata$direction, "forward")
})

# ---------------------------------------------------------------------------
# direction = "backward" (skip forward addition)
# ---------------------------------------------------------------------------

test_that("direction='backward' with start_covariates returns eres_scm_glm", {
  skip_on_cran()
  dat <- .make_glm_data()
  fit <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT"),
    significance_level = c(0.05, 0.01),
    dv_col             = "DV",
    family             = "binomial",
    direction          = "backward",
    start_covariates   = c("AGE", "WT"),
    center = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  expect_s3_class(fit, "eres_scm_glm")
  expect_equal(fit$metadata$direction, "backward")
})

test_that("direction='backward' result contains only start+exposure terms (no extras)", {
  skip_on_cran()
  dat <- .make_glm_data()
  fit <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT"),
    significance_level = c(0.05, 0.01),
    dv_col             = "DV",
    family             = "binomial",
    direction          = "backward",
    start_covariates   = c("AGE", "WT"),
    center = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  # finalRelation should only contain terms from start_covariates + exposure
  # (backward may have removed some, but cannot add new ones)
  allowed <- c("CMAX", "AGE", "WT")
  expect_true(all(fit$metadata$finalRelation %in% allowed))
})

# ---------------------------------------------------------------------------
# CoxPH validation
# ---------------------------------------------------------------------------

test_that("start_covariates validation works for scm_coxph", {
  dat <- .make_coxph_data()
  expect_error(
    scm_coxph(
      eres_data          = dat,
      exposure_vars      = "CMAX",
      covariates         = c("AGE", "WT"),
      significance_level = c(0.05, 0.01),
      dv_col             = "DV", time_col = "TIME",
      start_covariates   = "NONEXISTENT",
      center = FALSE, scale_to_median = FALSE, verbose = FALSE
    ),
    "start_covariates.*not in.*covariates"
  )
})

test_that("direction='backward' without start_covariates raises error in scm_coxph", {
  dat <- .make_coxph_data()
  expect_error(
    scm_coxph(
      eres_data          = dat,
      exposure_vars      = "CMAX",
      covariates         = c("AGE", "WT"),
      significance_level = c(0.05, 0.01),
      dv_col             = "DV", time_col = "TIME",
      direction          = "backward",
      center = FALSE, scale_to_median = FALSE, verbose = FALSE
    ),
    "direction.*backward.*start_covariates"
  )
})

test_that("direction='backward' with start_covariates works for scm_coxph", {
  skip_on_cran()
  dat <- .make_coxph_data()
  fit <- suppressMessages(scm_coxph(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT"),
    significance_level = c(0.05, 0.01),
    dv_col             = "DV", time_col = "TIME",
    direction          = "backward",
    start_covariates   = c("AGE", "WT"),
    center = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  expect_s3_class(fit, "eres_scm_coxph")
  expect_equal(fit$metadata$direction, "backward")
})

# ---------------------------------------------------------------------------
# Bootstrap replay
# ---------------------------------------------------------------------------

test_that("scm_bootstrap replays start_covariates correctly", {
  skip_on_cran()
  dat <- .make_glm_data()
  fit <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT"),
    significance_level = c(0.05, 0.01),
    dv_col             = "DV",
    family             = "binomial",
    start_covariates   = "AGE",
    center = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  s <- .extract_scm_settings(fit)
  expect_equal(s$start_covariates, "AGE")
  expect_equal(s$direction, "both")
})

test_that("scm_bootstrap replays direction='backward' correctly", {
  skip_on_cran()
  dat <- .make_glm_data()
  fit <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT"),
    significance_level = c(0.05, 0.01),
    dv_col             = "DV",
    family             = "binomial",
    direction          = "backward",
    start_covariates   = c("AGE", "WT"),
    center = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  s <- .extract_scm_settings(fit)
  expect_equal(s$direction, "backward")
  expect_equal(s$start_covariates, c("AGE", "WT"))
})
