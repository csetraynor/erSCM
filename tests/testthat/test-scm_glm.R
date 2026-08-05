# ---------------------------------------------------------------------------
# Early-exit (no significant exposure) always returns eres_scm_glm S3 object
# ---------------------------------------------------------------------------

test_that("scm_glm returns eres_scm_glm even when no exposure is significant", {
  set.seed(7)
  n <- 100
  # Pure noise: DV independent of CMAX
  dat <- data.frame(
    DV   = rbinom(n, 1L, 0.3),
    CMAX = rnorm(n, 100, 1),   # near-zero variance to ensure no signal
    AGE  = rnorm(n, 55, 12)
  )
  fit <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = "AGE",
    significance_level = c(0.001, 0.0001),  # very strict: guarantees early exit
    center             = FALSE,
    verbose            = FALSE
  ))
  expect_s3_class(fit, "eres_scm_glm")
  expect_true(inherits(fit, "eres_scm"))
  expect_null(fit$final_model)
  expect_true(is.na(fit$metadata$selected_exposure))
  expect_identical(fit$metadata$finalrun, FALSE)
  expect_length(fit$metadata$finalRelation, 0L)
  expect_length(fit$metadata$selected_terms, 0L)
  # Standard fields are present
  expect_true(!is.null(fit$model_data))
  expect_true(!is.null(fit$metadata))
})

test_that("print.eres_scm_glm handles no-model case without error", {
  set.seed(8)
  n <- 100
  dat <- data.frame(
    DV   = rbinom(n, 1L, 0.3),
    CMAX = rnorm(n, 100, 1),
    AGE  = rnorm(n, 55, 12)
  )
  fit <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = "AGE",
    significance_level = c(0.001, 0.0001),
    center             = FALSE
  ))
  out <- capture.output(print(fit))
  expect_true(any(grepl("No exposure", out, ignore.case = TRUE)))
  expect_true(any(grepl("not fit|was not fit", out, ignore.case = TRUE)))
})

test_that("plot_nonlinear_effect errors clearly when no model was fit", {
  set.seed(9)
  n <- 100
  dat <- data.frame(
    DV   = rbinom(n, 1L, 0.3),
    CMAX = rnorm(n, 100, 1),
    AGE  = rnorm(n, 55, 12)
  )
  fit <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = "AGE",
    significance_level = c(0.001, 0.0001),
    center             = FALSE
  ))
  expect_error(plot_nonlinear_effect(fit, "CMAX"), "Cannot plot")
})

test_that("scm_bootstrap errors clearly when original fit has no model", {
  set.seed(10)
  n <- 100
  dat <- data.frame(
    DV   = rbinom(n, 1L, 0.3),
    CMAX = rnorm(n, 100, 1),
    AGE  = rnorm(n, 55, 12)
  )
  fit <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = "AGE",
    significance_level = c(0.001, 0.0001),
    center             = FALSE
  ))
  expect_error(scm_bootstrap(fit, n_boot = 10L), "Cannot bootstrap")
})

# ---------------------------------------------------------------------------
# Normal fit: basic structure check
# ---------------------------------------------------------------------------

test_that("scm_glm runs and returns expected structure", {
  skip_on_cran()

  data("demo_eres", package = "erSCM")

  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX))

  result <- scm_glm(
    eres_data = demo_eres,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    significance_level = c(0.05, 0.01)
  )

  expect_type(result, "list")
  expect_named(result, c("final_model", "path", "summary", "model_data", "metadata", "var_spec", "naming_df"),
               ignore.order = TRUE)
  expect_true(inherits(result$final_model, "glm"))
  expect_s3_class(result$path, "data.frame")
  expect_s3_class(result$summary, "data.frame")
})
