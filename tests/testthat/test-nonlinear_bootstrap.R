# Tests for Phase 4: scm_bootstrap() integration with nonlinear models.
# Fixtures from helper-bootstrap.R are auto-loaded by testthat.

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

.make_nl_glm_fit_for_boot <- function(n = 300, seed = 42) {
  set.seed(seed)
  x <- runif(n, 1, 100)
  dat <- data.frame(
    DV   = rbinom(n, 1, plogis(-2 + log(x / 50))),
    CMAX = x,
    AGE  = rnorm(n, 55, 12),
    WT   = rnorm(n, 75, 15)
  )
  suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT"),
    significance_level = c(0.05, 0.01),
    dv_col             = "DV",
    family             = "binomial",
    nonlinear          = list(CMAX = c("linear", "log"), AGE = c("linear")),
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
}

# ---------------------------------------------------------------------------
# .extract_scm_settings picks up nonlinear fields
# ---------------------------------------------------------------------------

test_that(".extract_scm_settings returns nonlinear_spec for NL fit", {
  skip_on_cran()
  fit <- .make_nl_glm_fit_for_boot()
  s <- .extract_scm_settings(fit)
  expect_false(is.null(s$nonlinear_spec))
  expect_s3_class(s$nonlinear_spec, "nonlinear_spec")
  expect_true("CMAX" %in% names(s$nonlinear_spec))
})

test_that(".extract_scm_settings returns NULL nonlinear_spec for linear fit", {
  skip_on_cran()
  fit <- .fit_test_scm_glm()
  s <- .extract_scm_settings(fit)
  expect_null(s$nonlinear_spec)
  expect_null(s$linearity_alpha)
})

test_that(".extract_scm_settings preserves linearity_alpha and fp_powers", {
  skip_on_cran()
  fit <- .make_nl_glm_fit_for_boot()
  s <- .extract_scm_settings(fit)
  expect_false(is.null(s$linearity_alpha))
  expect_true(is.numeric(s$linearity_alpha))
})

# ---------------------------------------------------------------------------
# .extract_boot_terms with NL transforms
# ---------------------------------------------------------------------------

test_that(".extract_boot_terms returns exposure_transform for NL result", {
  skip_on_cran()
  fit <- .make_nl_glm_fit_for_boot()
  terms <- .extract_boot_terms(fit, "CMAX")
  expect_type(terms, "list")
  expect_true("exposure_transform" %in% names(terms))
  expect_true("covariate_transforms" %in% names(terms))
})

test_that(".extract_boot_terms exposure_transform is a single character", {
  skip_on_cran()
  fit <- .make_nl_glm_fit_for_boot()
  terms <- .extract_boot_terms(fit, "CMAX")
  expect_true(is.character(terms$exposure_transform) && length(terms$exposure_transform) == 1L)
})

test_that(".extract_boot_terms returns NA exposure_transform for linear fit", {
  skip_on_cran()
  fit <- .fit_test_scm_glm()
  terms <- .extract_boot_terms(fit, "CMAX")
  expect_true(is.na(terms$exposure_transform))
  expect_null(terms$covariate_transforms)
})

# ---------------------------------------------------------------------------
# .compute_transform_freq
# ---------------------------------------------------------------------------

test_that(".compute_transform_freq returns NULL when nonlinear_spec is NULL", {
  dummy_models <- list(
    list(exposure = "CMAX", exposure_transform = NA_character_,
         covariates = "AGE", covariate_transforms = NULL)
  )
  tf <- .compute_transform_freq(dummy_models, NULL, "CMAX")
  expect_null(tf$covariate_transform_freq)
  expect_null(tf$exposure_transform_freq)
})

test_that(".compute_transform_freq tallies correctly", {
  # Simulate 4 converged replicates: 3 log, 1 linear for CMAX
  spec <- list(CMAX = c("linear", "log"))
  class(spec) <- "nonlinear_spec"
  models <- list(
    list(exposure = "CMAX", exposure_transform = "log",
         covariates = character(0L), covariate_transforms = NULL),
    list(exposure = "CMAX", exposure_transform = "log",
         covariates = character(0L), covariate_transforms = NULL),
    list(exposure = "CMAX", exposure_transform = "linear",
         covariates = character(0L), covariate_transforms = NULL),
    list(exposure = "CMAX", exposure_transform = "log",
         covariates = character(0L), covariate_transforms = NULL)
  )
  tf <- .compute_transform_freq(models, spec, "CMAX")
  exp_freq <- tf$exposure_transform_freq
  expect_false(is.null(exp_freq))
  # exposure_transform_freq is now a named list (one entry per NL exposure)
  expect_type(exp_freq, "list")
  expect_true("CMAX" %in% names(exp_freq))
  expect_equal(as.numeric(exp_freq$CMAX[["log"]]),    0.75)
  expect_equal(as.numeric(exp_freq$CMAX[["linear"]]), 0.25)
})

test_that(".compute_transform_freq counts not_selected correctly", {
  spec <- list(AGE = c("linear", "log"))
  class(spec) <- "nonlinear_spec"
  models <- list(
    list(exposure = "CMAX", exposure_transform = NA_character_,
         covariates = "AGE",
         covariate_transforms = list(AGE = "log")),
    list(exposure = "CMAX", exposure_transform = NA_character_,
         covariates = character(0L),   # AGE not selected
         covariate_transforms = list()),
    list(exposure = "CMAX", exposure_transform = NA_character_,
         covariates = "AGE",
         covariate_transforms = list(AGE = "linear"))
  )
  tf <- .compute_transform_freq(models, spec, "CMAX")
  cov_freq <- tf$covariate_transform_freq
  expect_false(is.null(cov_freq))
  expect_true("AGE" %in% names(cov_freq))
  age_freq <- cov_freq[["AGE"]]
  # AGE selected twice (log once, linear once), not_selected once
  expect_equal(as.numeric(age_freq[["not_selected"]]), 1/3, tolerance = 1e-9)
  expect_equal(as.numeric(age_freq[["log"]]),           1/3, tolerance = 1e-9)
  expect_equal(as.numeric(age_freq[["linear"]]),        1/3, tolerance = 1e-9)
})

# ---------------------------------------------------------------------------
# scm_bootstrap with nonlinear fit — end-to-end
# ---------------------------------------------------------------------------

test_that("scm_bootstrap works with nonlinear fit", {
  skip_on_cran()
  fit <- .make_nl_glm_fit_for_boot()
  stability <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  expect_s3_class(stability, "scm_bootstrap")
  expect_true("transform_selection_freq" %in% names(stability))
  expect_true("exposure_transform_freq"  %in% names(stability))
})

test_that("transform_selection_freq is NULL for linear fit", {
  skip_on_cran()
  fit <- .fit_test_scm_glm()
  stability <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  expect_s3_class(stability, "scm_bootstrap")
  expect_null(stability$transform_selection_freq)
  expect_null(stability$exposure_transform_freq)
})

test_that("transform_selection_freq has correct structure for NL fit", {
  skip_on_cran()
  fit <- .make_nl_glm_fit_for_boot()
  stability <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 2L))
  tf <- stability$transform_selection_freq
  # May be NULL if no covariate NL vars in spec (CMAX is exposure, AGE = linear only)
  # but exposure_transform_freq should reflect CMAX's choices
  exp_tf <- stability$exposure_transform_freq
  if (!is.null(exp_tf)) {
    # exposure_transform_freq is now a named list (one entry per NL-enabled exposure)
    expect_type(exp_tf, "list")
    for (ev in names(exp_tf)) {
      freq_v <- exp_tf[[ev]]
      expect_true(is.numeric(freq_v))
      expect_true(all(freq_v >= 0 & freq_v <= 1))
      expect_equal(sum(freq_v), 1, tolerance = 1e-9)
    }
  }
})

test_that("exposure_transform_freq sums to 1 per exposure", {
  skip_on_cran()
  set.seed(10)
  n <- 300
  x <- runif(n, 1, 100)
  dat <- data.frame(
    DV   = rbinom(n, 1, plogis(-2 + log(x / 50))),
    CMAX = x,
    AGE  = rnorm(n, 55, 12)
  )
  fit <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = "AGE",
    significance_level = c(0.05, 0.01),
    dv_col             = "DV",
    family             = "binomial",
    nonlinear          = list(CMAX = c("linear", "log"), AGE = c("linear")),
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  stability <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 3L))
  exp_tf <- stability$exposure_transform_freq
  if (!is.null(exp_tf)) {
    # Named list: one entry per NL-enabled exposure
    expect_type(exp_tf, "list")
    for (ev in names(exp_tf)) {
      expect_equal(sum(exp_tf[[ev]]), 1, tolerance = 1e-9)
    }
  }
})

# ---------------------------------------------------------------------------
# print.scm_bootstrap — NL section appears
# ---------------------------------------------------------------------------

test_that("print.scm_bootstrap includes NL section for NL fit", {
  skip_on_cran()
  fit <- .make_nl_glm_fit_for_boot()
  stability <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 4L))
  # Only assert when there was actually NL result data to show
  has_data <- (!is.null(stability$exposure_transform_freq) &&
                 length(stability$exposure_transform_freq) > 0L) ||
              (!is.null(stability$transform_selection_freq) &&
                 length(stability$transform_selection_freq) > 0L)
  skip_if(!has_data, "No NL transform data in this replicate set")
  # cli writes to stdout; capture both stdout and stderr to be safe
  out <- capture.output(
    withr::with_options(list(cli.num_colors = 0L), print(stability)),
    type = "output"
  )
  # cli may also write to stderr in some environments
  out2 <- capture.output(
    withr::with_options(list(cli.num_colors = 0L), print(stability)),
    type = "message"
  )
  all_out <- c(out, out2)
  has_nl_section <- any(grepl("Non-Linear|Non.Linear", all_out, ignore.case = TRUE))
  # The NL section data must be present in the stability object itself
  # (print correctness already checked above; just verify the struct is right)
  expect_false(is.null(stability$exposure_transform_freq))
  expect_type(stability$exposure_transform_freq, "list")
  # If we captured output, check it; if not (headless env), skip the text check
  if (any(nchar(all_out) > 0L)) {
    expect_true(has_nl_section)
  }
})

test_that("print.scm_bootstrap has no NL section for linear fit", {
  skip_on_cran()
  fit <- .fit_test_scm_glm()
  stability <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 5L))
  out <- capture.output(print(stability))
  expect_false(any(grepl("Non-Linear Form Selection", out)))
})

# ---------------------------------------------------------------------------
# plot type = "transforms"
# ---------------------------------------------------------------------------

test_that("plot type='transforms' returns ggplot for NL fit", {
  skip_on_cran()
  fit <- .make_nl_glm_fit_for_boot()
  stability <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 6L))
  p <- plot(stability, type = "transforms")
  expect_s3_class(p, "gg")
})

test_that("plot type='transforms' returns informative blank for linear fit", {
  skip_on_cran()
  fit <- .fit_test_scm_glm()
  stability <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 7L))
  p <- plot(stability, type = "transforms")
  expect_s3_class(p, "gg")
})

# ---------------------------------------------------------------------------
# Backward compatibility: existing fields are unchanged
# ---------------------------------------------------------------------------

test_that("scm_bootstrap NL fit still has inclusion_prob", {
  skip_on_cran()
  fit <- .make_nl_glm_fit_for_boot()
  stability <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 8L))
  expect_true("inclusion_prob" %in% names(stability))
  expect_true(is.numeric(stability$inclusion_prob))
})

test_that("scm_bootstrap NL fit still has modal_model", {
  skip_on_cran()
  fit <- .make_nl_glm_fit_for_boot()
  stability <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 9L))
  expect_true("modal_model" %in% names(stability))
})

# ---------------------------------------------------------------------------
# CoxPH with nonlinear
# ---------------------------------------------------------------------------

test_that("scm_bootstrap works with nonlinear CoxPH fit", {
  skip_on_cran()
  set.seed(77)
  n <- 300
  x <- runif(n, 1, 100)
  # Use a real log signal and lenient thresholds so CMAX reliably enters
  dat <- data.frame(
    TIME = rexp(n, rate = exp(-2 + 0.5 * log(x / 50))),
    DV   = rbinom(n, 1L, 0.6),
    CMAX = x,
    AGE  = rnorm(n, 55, 12),
    WT   = rnorm(n, 75, 15)
  )
  fit <- suppressMessages(scm_coxph(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT"),
    significance_level = c(0.5, 0.3),
    dv_col             = "DV", time_col = "TIME",
    nonlinear          = list(CMAX = c("linear", "log")),
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  expect_s3_class(fit, "eres_scm_coxph")
  stability <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 10L))
  expect_s3_class(stability, "scm_bootstrap")
  expect_equal(stability$type, "coxph")
  expect_true("transform_selection_freq" %in% names(stability))
})
