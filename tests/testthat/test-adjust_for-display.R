## Tests for the display/classification feature tied to adjust_for in
## get_sum_tbl() and make_forest_plot().

# ---- get_sum_tbl() term_type column ----

test_that("get_sum_tbl adds term_type column when adjust_for is provided", {
  set.seed(1)
  d <- data.frame(y = rbinom(100, 1, 0.3), age = rnorm(100), wt = rnorm(100), dose = rnorm(100))
  mod <- glm(y ~ age + wt + dose, data = d, family = binomial)

  result <- get_sum_tbl(mod, adjust_for = c("age"))
  expect_true("term_type" %in% colnames(result))
})

test_that("get_sum_tbl term_type classifies intercept as base", {
  set.seed(1)
  d <- data.frame(y = rbinom(100, 1, 0.3), age = rnorm(100), wt = rnorm(100))
  mod <- glm(y ~ age + wt, data = d, family = binomial)

  result <- get_sum_tbl(mod, adjust_for = c("age"))
  intercept_row <- result[result$term == "(Intercept)", ]
  if (nrow(intercept_row) > 0) {
    expect_equal(intercept_row$term_type, "base")
  }
})

test_that("get_sum_tbl term_type classifies adjust_for terms as adjusted", {
  set.seed(1)
  d <- data.frame(y = rbinom(100, 1, 0.3), age = rnorm(100), wt = rnorm(100))
  mod <- glm(y ~ age + wt, data = d, family = binomial)

  result <- get_sum_tbl(mod, adjust_for = c("age"))
  age_row <- result[result$term == "age", ]
  expect_equal(age_row$term_type, "adjusted")
})

test_that("get_sum_tbl term_type classifies non-adjust_for terms as selected", {
  set.seed(1)
  d <- data.frame(y = rbinom(100, 1, 0.3), age = rnorm(100), wt = rnorm(100))
  mod <- glm(y ~ age + wt, data = d, family = binomial)

  result <- get_sum_tbl(mod, adjust_for = c("age"))
  wt_row <- result[result$term == "wt", ]
  expect_equal(wt_row$term_type, "selected")
})

test_that("get_sum_tbl show_terms = 'selected' keeps only selected rows", {
  set.seed(1)
  d <- data.frame(y = rbinom(100, 1, 0.3), age = rnorm(100), wt = rnorm(100))
  mod <- glm(y ~ age + wt, data = d, family = binomial)

  result <- get_sum_tbl(mod, adjust_for = c("age"), show_terms = "selected")
  expect_true(all(result$term_type == "selected"))
  expect_false("age" %in% result$term)
})

test_that("get_sum_tbl show_terms = 'adjusted' keeps only adjusted rows", {
  set.seed(1)
  d <- data.frame(y = rbinom(100, 1, 0.3), age = rnorm(100), wt = rnorm(100))
  mod <- glm(y ~ age + wt, data = d, family = binomial)

  result <- get_sum_tbl(mod, adjust_for = c("age"), show_terms = "adjusted")
  expect_true(all(result$term_type == "adjusted"))
  expect_true("age" %in% result$term)
})

test_that("get_sum_tbl without adjust_for has no term_type column", {
  set.seed(1)
  d <- data.frame(y = rbinom(100, 1, 0.3), age = rnorm(100), wt = rnorm(100))
  mod <- glm(y ~ age + wt, data = d, family = binomial)

  result <- get_sum_tbl(mod)
  expect_false("term_type" %in% colnames(result))
})

# ---- make_forest_plot() term_type column ----

test_that("make_forest_plot returns term_type column as 'main' when adjust_for = NULL", {
  skip_on_cran()
  data("demo_eres", package = "erSCM")
  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")
  d <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)),
                                ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX))

  result <- scm_glm(eres_data = d, covariates = c("WT", "AGE"),
                    exposure_vars = all_exposure, dv_col = "sim_status",
                    significance_level = c(0.9, 0.9), center = FALSE)

  fp <- make_forest_plot(result, plot = FALSE)
  expect_true("term_type" %in% colnames(fp))
  expect_true(all(fp$term_type %in% c("main", "interaction")))
})

test_that("make_forest_plot term_type contains 'adjusted' when adjust_for provided", {
  skip_on_cran()
  data("demo_eres", package = "erSCM")
  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")
  d <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)),
                                ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX))

  result <- scm_glm(eres_data = d, covariates = c("WT", "AGE"),
                    exposure_vars = all_exposure, dv_col = "sim_status",
                    significance_level = c(0.9, 0.9), adjust_for = "AGE", center = FALSE)

  fp <- make_forest_plot(result, plot = FALSE, adjust_for = "AGE")
  expect_true("adjusted" %in% fp$term_type)
})

test_that("make_forest_plot show_terms = 'adjusted' drops selected panels", {
  skip_on_cran()
  data("demo_eres", package = "erSCM")
  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")
  d <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)),
                                ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX))

  result <- scm_glm(eres_data = d, covariates = c("WT", "AGE"),
                    exposure_vars = all_exposure, dv_col = "sim_status",
                    significance_level = c(0.9, 0.9), adjust_for = "AGE", center = FALSE)

  fp <- make_forest_plot(result, plot = FALSE, adjust_for = "AGE", show_terms = "adjusted")
  expect_false("selected" %in% fp$term_type)
})
