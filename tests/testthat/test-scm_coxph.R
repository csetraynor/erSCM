test_that("scm_coxph runs and returns expected structure", {
 skip_on_cran()

  data("demo_eres", package = "erSCM")

  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) %>%
    dplyr::mutate(SEX = as.factor(SEX)) %>%
    dplyr::mutate(DVID = "ALL")

  result <- scm_coxph(
    eres_data = demo_eres,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    time_col = "sim_time",
    significance_level = c(0.9, 0.9),
    dvid = "ALL",
    dvid_col = "DVID"
  )

  expect_type(result, "list")
  expect_named(result, c("final_model", "path", "summary", "model_data", "metadata", "var_spec", "naming_df"),
               ignore.order = TRUE)
  expect_true(inherits(result$final_model, "coxph"))
  expect_s3_class(result$path, "data.frame")
  expect_s3_class(result$summary, "data.frame")
})

test_that("scm_coxph path has non-NA nobs (regression: broom::glance uses 'n' not 'nobs' for coxph)", {
  skip_on_cran()
  set.seed(42)
  n <- 300
  # Two covariates so backward step actually fires
  dat <- data.frame(
    TIME = rexp(n, rate = 0.1),
    DV   = rbinom(n, 1, 0.4),
    EXP  = rnorm(n),
    AGE  = rnorm(n)
  )
  result <- suppressMessages(scm_coxph(
    eres_data          = dat,
    covariates         = "AGE",
    exposure_vars      = "EXP",
    significance_level = c(0.9, 0.9),
    dv_col             = "DV",
    time_col           = "TIME",
    center             = FALSE,
    verbose            = FALSE
  ))
  # The path table contains rows from all backward steps; none should have NA nobs
  expect_false(any(is.na(result$path$nobs)),
               info = "nobs should not be NA in coxph path (was caused by $nobs vs $n mismatch)")
})
