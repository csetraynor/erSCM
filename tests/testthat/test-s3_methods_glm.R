test_that("scm_glm returns eres_scm_glm class", {
  skip_on_cran()
  data("demo_eres", package = "erSCM")

  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
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

  expect_s3_class(result, "eres_scm_glm")
  expect_s3_class(result, "eres_scm")
  # list access still works

  expect_true(inherits(result$final_model, "glm"))
})

test_that("scm_glm verbose produces messages", {
  skip_on_cran()
  data("demo_eres", package = "erSCM")

  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX))

  expect_message(
    scm_glm(
      eres_data = demo_eres,
      covariates = all_covariates,
      exposure_vars = all_exposure,
      dv_col = "sim_status",
      significance_level = c(0.05, 0.01),
      verbose = TRUE
    ),
    "SCM Forward Selection"
  )
})

test_that("print.eres_scm_glm runs and returns invisibly", {
  skip_on_cran()
  data("demo_eres", package = "erSCM")

  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
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

  out <- capture.output(ret <- print(result))
  expect_identical(ret, result)
  expect_true(length(out) > 0)
})

test_that("plot.eres_scm_glm runs without error", {
  skip_on_cran()
  data("demo_eres", package = "erSCM")

  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
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

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(result))
})
