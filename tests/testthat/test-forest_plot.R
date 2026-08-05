test_that("make_forest_plot does not error on glm result", {
  skip_on_cran()

  data("demo_eres", package = "erSCM")

  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX))

  lr_search <- scm_glm(
    eres_data = demo_eres,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    significance_level = c(0.05, 0.01)
  )

  # Test that forest plot returns data when plot = FALSE
  fp_data <- make_forest_plot(lr_search, plot = FALSE)
  expect_s3_class(fp_data, "data.frame")
  expect_true(nrow(fp_data) > 0)

  # Test that forest plot returns a ggplot/patchwork when plot = TRUE
  p <- make_forest_plot(lr_search, plot = TRUE)
  expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
})
