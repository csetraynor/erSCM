test_that("scm_coxph returns eres_scm_coxph class", {
  skip_on_cran()
  data("demo_eres", package = "erSCM")

  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX)) |>
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

  expect_s3_class(result, "eres_scm_coxph")
  expect_s3_class(result, "eres_scm")
  expect_true(inherits(result$final_model, "coxph"))
})

test_that("scm_coxph verbose produces messages", {
  skip_on_cran()
  data("demo_eres", package = "erSCM")

  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX)) |>
    dplyr::mutate(DVID = "ALL")

  expect_message(
    scm_coxph(
      eres_data = demo_eres,
      covariates = all_covariates,
      exposure_vars = all_exposure,
      dv_col = "sim_status",
      time_col = "sim_time",
      significance_level = c(0.9, 0.9),
      dvid = "ALL",
      dvid_col = "DVID",
      verbose = TRUE
    ),
    "SCM Forward Selection"
  )
})

test_that("print.eres_scm_coxph runs and returns invisibly", {
  skip_on_cran()
  data("demo_eres", package = "erSCM")

  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX)) |>
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

  out <- capture.output(ret <- print(result))
  expect_identical(ret, result)
  expect_true(length(out) > 0)
})

test_that("plot.eres_scm_coxph runs without error", {
  skip_on_cran()
  data("demo_eres", package = "erSCM")

  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX)) |>
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

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(result))
})
