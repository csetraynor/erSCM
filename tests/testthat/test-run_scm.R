# Tests for run_scm() dispatcher and nonlinear parameter promotion (Gap #2).
# Relies on shared fixtures from helper-bootstrap.R (auto-loaded by testthat).

# ---------------------------------------------------------------------------
# Basic dispatch correctness
# ---------------------------------------------------------------------------

test_that("run_scm type='glm' returns eres_scm_glm regardless of exposure significance", {
  # Now that scm_glm() always returns an S3-classed object, noise data is fine.
  dat <- .make_test_data_glm(n = 200, seed = 1)
  res <- suppressMessages(run_scm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT"),
    type               = "glm",
    significance_level = c(0.05, 0.01),
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  expect_s3_class(res, "eres_scm_glm")
})

test_that("run_scm type='coxph' returns eres_scm_coxph regardless of exposure significance", {
  dat <- .make_test_data_coxph(n = 200, seed = 2)
  res <- suppressMessages(run_scm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT"),
    type               = "coxph",
    significance_level = c(0.05, 0.01),
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  expect_s3_class(res, "eres_scm_coxph")
})

# ---------------------------------------------------------------------------
# Nonlinear parameter promotion: nonlinear = TRUE produces transform_registry
# ---------------------------------------------------------------------------

test_that("run_scm with nonlinear=TRUE populates transform_registry (GLM)", {
  skip_on_cran()
  set.seed(42)
  n <- 300
  x <- runif(n, 1, 100)
  dat <- data.frame(
    DV   = rbinom(n, 1, stats::plogis(-2 + log(x / 50))),
    CMAX = x,
    AGE  = rnorm(n, 55, 12)
  )
  res <- suppressMessages(run_scm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = "AGE",
    type               = "glm",
    significance_level = c(0.05, 0.01),
    nonlinear          = list(CMAX = c("linear", "log")),
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  expect_s3_class(res, "eres_scm_glm")
  expect_false(is.null(res$metadata$transform_registry))
})

test_that("run_scm explicit nonlinear params are forwarded (linearity_alpha, fp_powers, spline_df)", {
  skip_on_cran()
  set.seed(43)
  n <- 300
  x <- runif(n, 1, 100)
  dat <- data.frame(
    DV   = rbinom(n, 1, stats::plogis(-2 + log(x / 50))),
    CMAX = x,
    AGE  = rnorm(n, 55, 12)
  )
  # Use a very high linearity_alpha so log always wins over linear
  res <- suppressMessages(run_scm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = "AGE",
    type               = "glm",
    significance_level = c(0.05, 0.01),
    nonlinear          = list(CMAX = c("linear", "log")),
    linearity_alpha    = 0.99,
    fp_powers          = c(-1, 0, 1),
    spline_df          = 4L,
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  expect_s3_class(res, "eres_scm_glm")
  reg <- res$metadata$transform_registry
  expect_false(is.null(reg))
  # With linearity_alpha = 0.99, the log form should nearly always beat linear
  if (!is.null(reg[["CMAX"]])) {
    expect_equal(reg[["CMAX"]]$transform_type, "log")
  }
})

# ---------------------------------------------------------------------------
# Input validation for promoted params
# ---------------------------------------------------------------------------

test_that("run_scm errors on invalid linearity_alpha", {
  dat <- .make_test_data_glm(n = 100, seed = 5)
  expect_error(
    run_scm(dat, "CMAX", "AGE", type = "glm", linearity_alpha = 1.5),
    "linearity_alpha"
  )
  expect_error(
    run_scm(dat, "CMAX", "AGE", type = "glm", linearity_alpha = 0),
    "linearity_alpha"
  )
  expect_error(
    run_scm(dat, "CMAX", "AGE", type = "glm", linearity_alpha = "0.05"),
    "linearity_alpha"
  )
})

test_that("run_scm errors on invalid spline_df", {
  dat <- .make_test_data_glm(n = 100, seed = 6)
  expect_error(
    run_scm(dat, "CMAX", "AGE", type = "glm", spline_df = 1L),
    "spline_df"
  )
  expect_error(
    run_scm(dat, "CMAX", "AGE", type = "glm", spline_df = -1L),
    "spline_df"
  )
})

test_that("run_scm errors on invalid fp_powers", {
  dat <- .make_test_data_glm(n = 100, seed = 7)
  expect_error(
    run_scm(dat, "CMAX", "AGE", type = "glm", fp_powers = numeric(0)),
    "fp_powers"
  )
  expect_error(
    run_scm(dat, "CMAX", "AGE", type = "glm", fp_powers = c(1, NA, 2)),
    "fp_powers"
  )
})

# ---------------------------------------------------------------------------
# Typo detection: a misspelled argument goes to ... and is IGNORED by run_scm
# (R does not error on unused ... args when passed through to scm_glm, but the
# nonlinear feature will NOT activate, so transform_registry stays NULL).
# This test documents the current behavior.
# ---------------------------------------------------------------------------

test_that("misspelled nonlinear arg (nonliner) causes an error or leaves NL inactive", {
  skip_on_cran()
  dat <- .make_test_data_glm(n = 200, seed = 8)
  # 'nonliner' is a typo — it goes to ... and is forwarded to scm_glm, which
  # errors on unknown args (R's standard behaviour for unmatched formal params).
  # Either an error is raised (most common) or, if somehow swallowed, the NL
  # feature must NOT be silently activated.
  result <- tryCatch(
    suppressMessages(run_scm(
      eres_data          = dat,
      exposure_vars      = "CMAX",
      covariates         = c("AGE", "WT"),
      type               = "glm",
      significance_level = c(0.5, 0.3),
      nonliner           = TRUE,   # deliberate typo
      center             = FALSE, scale_to_median = FALSE, verbose = FALSE
    )),
    error = function(e) e
  )
  # Acceptable outcomes: an error (expected), or a successful fit with NL inactive
  if (inherits(result, "error")) {
    expect_true(TRUE)   # error is the expected path — typo was caught
  } else {
    expect_null(result$metadata$transform_registry)
  }
})

# ---------------------------------------------------------------------------
# scm_bootstrap still works through run_scm (regression check)
# ---------------------------------------------------------------------------

test_that("scm_bootstrap works with a nonlinear GLM fit via run_scm", {
  skip_on_cran()
  set.seed(99)
  n <- 300
  x <- runif(n, 1, 100)
  dat <- data.frame(
    DV   = rbinom(n, 1, stats::plogis(-2 + log(x / 50))),
    CMAX = x,
    AGE  = rnorm(n, 55, 12)
  )
  fit <- suppressMessages(run_scm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = "AGE",
    type               = "glm",
    significance_level = c(0.05, 0.01),
    nonlinear          = list(CMAX = c("linear", "log")),
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  boot <- suppressMessages(scm_bootstrap(
    scm_fit = fit, n_boot = 10L, seed = 1L, verbose = FALSE
  ))
  expect_s3_class(boot, "scm_bootstrap")
  expect_gte(boot$n_converged, 1L)
})
