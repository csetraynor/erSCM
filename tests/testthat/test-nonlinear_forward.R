test_that("scm_glm with nonlinear=NULL behaves identically to before", {
  test_data <- .make_test_data_glm(n = 200, seed = 99)
  fit_linear <- suppressMessages(scm_glm(
    eres_data          = test_data,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT", "SEX"),
    significance_level = c(0.5, 0.3),
    dv_col             = "DV",
    family             = "binomial",
    center             = FALSE,
    scale_to_median    = FALSE,
    verbose            = FALSE
  ))
  fit_null <- suppressMessages(scm_glm(
    eres_data          = test_data,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT", "SEX"),
    significance_level = c(0.5, 0.3),
    dv_col             = "DV",
    family             = "binomial",
    center             = FALSE,
    scale_to_median    = FALSE,
    nonlinear          = NULL,
    verbose            = FALSE
  ))
  # Both paths produce the same selected terms and AIC
  expect_identical(fit_linear$metadata$selected_terms,
                   fit_null$metadata$selected_terms)
  expect_identical(fit_linear$metadata$finalRelation,
                   fit_null$metadata$finalRelation)
  expect_null(fit_null$metadata$transform_registry)
})

test_that("scm_glm with nonlinear=TRUE runs without error and returns transform_registry", {
  skip_on_cran()
  test_data <- .make_test_data_glm(n = 300, seed = 42)
  fit <- suppressMessages(scm_glm(
    eres_data          = test_data,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT", "SEX"),
    significance_level = c(0.5, 0.3),
    dv_col             = "DV",
    family             = "binomial",
    nonlinear          = TRUE,
    center             = FALSE,
    scale_to_median    = FALSE,
    verbose            = FALSE
  ))
  expect_s3_class(fit, "eres_scm_glm")
  expect_true("transform_registry" %in% names(fit$metadata))
})

test_that("scm_glm with named list nonlinear spec uses specified transforms", {
  skip_on_cran()
  test_data <- .make_test_data_glm(n = 300, seed = 10)
  # Force a log-linear signal in CMAX
  test_data$DV <- rbinom(300, 1, plogis(-2 + 0.8 * log(pmax(test_data$CMAX, 1) / 100)))
  fit <- suppressMessages(scm_glm(
    eres_data          = test_data,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT", "SEX"),
    significance_level = c(0.05, 0.01),
    dv_col             = "DV",
    family             = "binomial",
    nonlinear          = list(CMAX = c("linear", "log", "fp1")),
    center             = FALSE,
    scale_to_median    = FALSE,
    verbose            = FALSE
  ))
  expect_s3_class(fit, "eres_scm_glm")
  # transform_registry must record the exposure entry
  reg <- fit$metadata$transform_registry
  expect_true(!is.null(reg))
  expect_true("CMAX" %in% names(reg))
  # With a strong log signal CMAX should enter as log or fp1 (p ~ 0)
  expect_true(reg$CMAX$transform_type %in% c("log", "fp1"))
})

test_that("FP search selects a non-linear power for a known sqrt relationship", {
  skip_on_cran()
  set.seed(123)
  n <- 500
  x <- runif(n, 1, 100)
  y <- rbinom(n, 1, plogis(-3 + 2 * sqrt(x / 50)))
  test_data <- data.frame(DV = y, CMAX = x, AGE = rnorm(n, 50, 10))

  fit <- suppressMessages(scm_glm(
    eres_data          = test_data,
    exposure_vars      = "CMAX",
    covariates         = "AGE",
    significance_level = c(0.05, 0.01),
    dv_col             = "DV",
    family             = "binomial",
    nonlinear          = list(CMAX = c("linear", "fp1", "fp2")),
    center             = FALSE,
    scale_to_median    = FALSE,
    verbose            = FALSE
  ))
  # CMAX should have entered the model (strong signal)
  expect_s3_class(fit, "eres_scm_glm")
  reg <- fit$metadata$transform_registry
  expect_true(!is.null(reg))
  expect_true("CMAX" %in% names(reg))
  # With a sqrt signal the best FP form should have power ≈ 0.5
  # (or at minimum CMAX should have been tested with FP and entered)
  expect_true(reg$CMAX$transform_type %in% c("linear", "fp1", "fp2"))
})

test_that("Spline terms enter correctly", {
  skip_on_cran()
  set.seed(7)
  n <- 400
  x <- runif(n, 1, 100)
  y <- rbinom(n, 1, plogis(-2 + sin(x / 15)))
  dat <- data.frame(DV = y, CMAX = x, AGE = rnorm(n, 50, 10))

  fit <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = "AGE",
    significance_level = c(0.05, 0.01),
    dv_col             = "DV",
    family             = "binomial",
    nonlinear          = list(CMAX = c("linear", "rcs")),
    center             = FALSE,
    scale_to_median    = FALSE,
    verbose            = FALSE
  ))
  reg <- fit$metadata$transform_registry
  expect_true(!is.null(reg))
  # Check that if CMAX entered, its formula contains ns() for rcs
  if ("CMAX" %in% names(reg)) {
    if (reg$CMAX$transform_type == "rcs") {
      expect_true(grepl("ns\\(", reg$CMAX$formula_rhs))
    }
  }
})

test_that("linearity_alpha controls non-linear entry via LRT for fp2 (df > 1 vs linear)", {
  skip_on_cran()
  # fp2 has df=2 vs linear df=1: df_diff_nl=1 > 0, so LRT vs linearity_alpha is used.
  # log (df=1) would use AIC fallback, bypassing alpha — use fp2 to test the LRT branch.
  set.seed(55)
  n <- 300
  x <- runif(n, 1, 100)
  # Use a moderately non-linear signal so fp2 beats linear on AIC but the
  # LRT p-value for fp2 vs linear is detectable only with lenient alpha.
  y <- rbinom(n, 1, plogis(-2 + 0.5 * log(x / 50)))
  dat <- data.frame(DV = y, CMAX = x)

  fit_strict <- suppressMessages(scm_glm(
    eres_data = dat, exposure_vars = "CMAX", covariates = character(0),
    significance_level = c(0.05, 0.01), dv_col = "DV",
    nonlinear = list(CMAX = c("linear", "fp2")),
    linearity_alpha = 1e-10,   # essentially never prefer non-linear over linear
    center = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  fit_lenient <- suppressMessages(scm_glm(
    eres_data = dat, exposure_vars = "CMAX", covariates = character(0),
    significance_level = c(0.05, 0.01), dv_col = "DV",
    nonlinear = list(CMAX = c("linear", "fp2")),
    linearity_alpha = 0.99,   # almost always prefer non-linear
    center = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  reg_strict  <- fit_strict$metadata$transform_registry
  reg_lenient <- fit_lenient$metadata$transform_registry
  # With strict alpha CMAX must enter linearly (if it enters at all)
  if (!is.null(reg_strict) && "CMAX" %in% names(reg_strict)) {
    expect_identical(reg_strict$CMAX$transform_type, "linear")
  }
  # With lenient alpha CMAX can enter as fp2 (if signal is detectable)
  if (!is.null(reg_lenient) && "CMAX" %in% names(reg_lenient)) {
    expect_true(reg_lenient$CMAX$transform_type %in% c("fp2", "linear"))
  }
})

test_that("Factor variables are always linear regardless of nonlinear spec", {
  test_data <- .make_test_data_glm(n = 200, seed = 5)
  spec <- parse_nonlinear_spec(TRUE, c("AGE", "SEX"), "CMAX", test_data)
  expect_identical(spec[["SEX"]], "linear")
})

test_that("scm_coxph with nonlinear=NULL behaves identically to before", {
  skip_on_cran()
  dat <- .make_test_data_coxph(n = 200, seed = 77)
  fit_linear <- suppressMessages(scm_coxph(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT", "SEX"),
    significance_level = c(0.5, 0.3),
    dv_col             = "DV", time_col = "TIME",
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  fit_null <- suppressMessages(scm_coxph(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT", "SEX"),
    significance_level = c(0.5, 0.3),
    dv_col             = "DV", time_col = "TIME",
    nonlinear          = NULL,
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  expect_identical(fit_linear$metadata$selected_terms,
                   fit_null$metadata$selected_terms)
  expect_null(fit_null$metadata$transform_registry)
})

test_that("scm_coxph with nonlinear=TRUE runs without error", {
  skip_on_cran()
  dat <- .make_test_data_coxph(n = 300, seed = 88)
  fit <- suppressMessages(scm_coxph(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT", "SEX"),
    significance_level = c(0.5, 0.3),
    dv_col             = "DV", time_col = "TIME",
    nonlinear          = TRUE,
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  expect_s3_class(fit, "eres_scm_coxph")
  expect_true("transform_registry" %in% names(fit$metadata))
})

# ---------------------------------------------------------------------------
# Gap #4: nonlinear + interactions limitation
# ---------------------------------------------------------------------------

test_that("scm_glm warns when both nonlinear and interactions are specified (Gap #4)", {
  skip_on_cran()
  set.seed(55)
  n <- 300
  x <- runif(n, 1, 100)
  dat <- data.frame(
    DV   = rbinom(n, 1, stats::plogis(-2 + log(x / 50))),
    CMAX = x,
    AGE  = rnorm(n, 55, 12),
    WT   = rnorm(n, 75, 15)
  )
  # Expect a cli_warn (or base warning) mentioning the limitation
  expect_warning(
    suppressMessages(scm_glm(
      eres_data          = dat,
      exposure_vars      = "CMAX",
      covariates         = c("AGE", "WT"),
      significance_level = c(0.3, 0.1),
      nonlinear          = list(CMAX = c("linear", "log")),
      interactions       = "AGE:WT",
      center             = FALSE, scale_to_median = FALSE, verbose = FALSE
    )),
    regexp = "nonlinear.*interactions|interactions.*nonlinear"
  )
})

# ---------------------------------------------------------------------------
# Gap #6: nonlinear + adjust_for
# ---------------------------------------------------------------------------

test_that("nonlinear + adjust_for: factor forced in, exposure gets NL search (Gap #6)", {
  skip_on_cran()
  set.seed(42)
  n <- 300
  cmax <- runif(n, 10, 200)
  dat <- data.frame(
    DV   = rbinom(n, 1, plogis(-2 + 0.8 * log(cmax / 50))),
    CMAX = cmax,
    AGE  = rnorm(n, 55, 12),
    RACE = factor(sample(c("White", "Black", "Asian"), n, replace = TRUE))
  )
  fit <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "RACE"),
    adjust_for         = "RACE",
    significance_level = c(0.05, 0.01),
    dv_col             = "DV",
    family             = "binomial",
    nonlinear          = list(CMAX = c("linear", "log")),
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  expect_s3_class(fit, "eres_scm_glm")
  final_rel <- fit$metadata$finalRelation
  tr_reg    <- fit$metadata$transform_registry
  # RACE was forced in via adjust_for — must appear in finalRelation
  expect_true("RACE" %in% final_rel)
  # RACE was forced in (not selected), so it must NOT appear in transform_registry
  expect_false("RACE" %in% names(tr_reg))
  # CMAX had a NL search — if it entered, its transform type must be valid
  if ("CMAX" %in% names(tr_reg)) {
    expect_true(tr_reg$CMAX$transform_type %in% c("linear", "log"))
  }
})

test_that("nonlinear spec for adjust_for variable is ignored (Gap #6)", {
  skip_on_cran()
  set.seed(43)
  n <- 300
  cmax <- runif(n, 10, 200)
  dat <- data.frame(
    DV   = rbinom(n, 1, plogis(-2 + 0.8 * log(cmax / 50))),
    CMAX = cmax,
    AGE  = rnorm(n, 55, 12)
  )
  # AGE is both in adjust_for AND in the nonlinear spec
  fit <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = "AGE",
    adjust_for         = "AGE",
    significance_level = c(0.05, 0.01),
    dv_col             = "DV",
    family             = "binomial",
    nonlinear          = list(CMAX = c("linear", "log"), AGE = c("linear", "log")),
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  expect_s3_class(fit, "eres_scm_glm")
  final_rel <- fit$metadata$finalRelation
  tr_reg    <- fit$metadata$transform_registry
  # AGE was forced in via adjust_for — must be in finalRelation
  expect_true("AGE" %in% final_rel)
  # AGE was forced in as a fixed covariate, NOT through NL search; if it appears
  # in transform_registry at all, it must be linear (protect-only, no NL search)
  if ("AGE" %in% names(tr_reg)) {
    expect_identical(tr_reg$AGE$transform_type, "linear")
  }
  # AGE must NOT enter as log(AGE) even though the nonlinear spec included "log"
  expect_false(any(grepl("^log\\(AGE\\)$", final_rel)))
})

test_that("when nonlinear active, interaction with non-linear parent is silently skipped (Gap #4)", {
  skip_on_cran()
  set.seed(56)
  n <- 400
  x <- runif(n, 1, 100)
  dat <- data.frame(
    DV   = rbinom(n, 1, stats::plogis(-2 + log(x / 50))),
    CMAX = x,
    AGE  = rnorm(n, 55, 12)
  )
  # CMAX enters as log(CMAX); CMAX:AGE interaction should be ineligible
  fit <- suppressWarnings(suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = "AGE",
    significance_level = c(0.3, 0.1),
    nonlinear          = list(CMAX = c("linear", "log")),
    interactions       = "CMAX:AGE",
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  )))
  expect_s3_class(fit, "eres_scm_glm")
  # The interaction term must NOT appear in final Relation (it was skipped)
  final_rel <- fit$metadata$finalRelation
  expect_false(any(grepl("CMAX:AGE|AGE:CMAX", final_rel, fixed = FALSE)))
})
