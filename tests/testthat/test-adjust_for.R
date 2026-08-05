## Tests for the adjust_for parameter in scm_glm() and scm_coxph()
## These tests use skip_on_cran() because they fit real models.

# ---- Shared test data setup helpers ----

make_glm_data <- function() {
  data("demo_eres", package = "erSCM")
  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")
  demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)),
                                ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX))
}

make_cox_data <- function() {
  data("demo_eres", package = "erSCM")
  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")
  demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)),
                                ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX)) |>
    dplyr::mutate(DVID = "ALL")
}

all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
all_exposure   <- c("CMAXss", "CMINss", "AUCss")

# ---- Input validation tests (fast, no model fitting) ----

test_that("adjust_for must be a character vector", {
  d <- make_glm_data()
  expect_error(
    scm_glm(eres_data = d, covariates = all_covariates,
            exposure_vars = all_exposure, dv_col = "sim_status",
            adjust_for = 123),
    "must be a character vector"
  )
})

test_that("adjust_for rejects non-existent columns", {
  d <- make_glm_data()
  expect_error(
    scm_glm(eres_data = d, covariates = all_covariates,
            exposure_vars = all_exposure, dv_col = "sim_status",
            adjust_for = c("NONEXISTENT_COL")),
    "not found in the data"
  )
})

test_that("adjust_for = NULL gives same result as omitting it (scm_glm)", {
  skip_on_cran()
  d <- make_glm_data()
  set.seed(42)
  r1 <- scm_glm(eres_data = d, covariates = all_covariates,
                exposure_vars = all_exposure, dv_col = "sim_status",
                significance_level = c(0.05, 0.01), adjust_for = NULL,
                center = FALSE)
  set.seed(42)
  r2 <- scm_glm(eres_data = d, covariates = all_covariates,
                exposure_vars = all_exposure, dv_col = "sim_status",
                significance_level = c(0.05, 0.01), center = FALSE)
  expect_equal(r1$metadata$finalRelation, r2$metadata$finalRelation)
})

# ---- Functional tests: GLM ----

test_that("adjust_for terms appear in final model regardless of significance (scm_glm)", {
  skip_on_cran()
  d <- make_glm_data()
  result <- scm_glm(
    eres_data = d,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    significance_level = c(0.05, 0.01),
    adjust_for = c("AGE", "SEX"),
    center = FALSE,
    verbose = FALSE
  )
  final_terms <- result$metadata$finalRelation
  expect_true("AGE" %in% final_terms,
              info = "AGE must be in final model as an adjust_for term")
  expect_true("SEX" %in% final_terms,
              info = "SEX must be in final model as an adjust_for term")
})

test_that("adjust_for terms are not listed in selected_terms (scm_glm)", {
  skip_on_cran()
  d <- make_glm_data()
  result <- scm_glm(
    eres_data = d,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    significance_level = c(0.05, 0.01),
    adjust_for = c("AGE"),
    center = FALSE,
    verbose = FALSE
  )
  expect_false("AGE" %in% result$metadata$selected_terms,
               info = "AGE should not appear in selected_terms (it is a forced adjustment)")
})

test_that("adjust_for is stored in metadata (scm_glm)", {
  skip_on_cran()
  d <- make_glm_data()
  result <- scm_glm(
    eres_data = d,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    significance_level = c(0.05, 0.01),
    adjust_for = c("AGE", "SEX"),
    center = FALSE,
    verbose = FALSE
  )
  expect_equal(result$metadata$adjust_for, c("AGE", "SEX"))
})

# ---- Functional tests: Cox PH ----

test_that("adjust_for must be a character vector (scm_coxph)", {
  d <- make_cox_data()
  expect_error(
    scm_coxph(eres_data = d, covariates = all_covariates,
              exposure_vars = all_exposure, dv_col = "sim_status",
              time_col = "sim_time", dvid = "ALL", dvid_col = "DVID",
              adjust_for = 123),
    "must be a character vector"
  )
})

test_that("adjust_for rejects non-existent columns (scm_coxph)", {
  d <- make_cox_data()
  expect_error(
    scm_coxph(eres_data = d, covariates = all_covariates,
              exposure_vars = all_exposure, dv_col = "sim_status",
              time_col = "sim_time", dvid = "ALL", dvid_col = "DVID",
              adjust_for = "NO_SUCH_COL"),
    "not found in the data"
  )
})

test_that("adjust_for terms appear in final model (scm_coxph)", {
  skip_on_cran()
  d <- make_cox_data()
  result <- scm_coxph(
    eres_data = d,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    time_col = "sim_time",
    significance_level = c(0.9, 0.9),
    dvid = "ALL",
    dvid_col = "DVID",
    adjust_for = c("AGE"),
    center = FALSE,
    verbose = FALSE
  )
  expect_true("AGE" %in% result$metadata$finalRelation)
})

test_that("adjust_for terms not in selected_terms (scm_coxph)", {
  skip_on_cran()
  d <- make_cox_data()
  result <- scm_coxph(
    eres_data = d,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    time_col = "sim_time",
    significance_level = c(0.9, 0.9),
    dvid = "ALL",
    dvid_col = "DVID",
    adjust_for = c("AGE"),
    center = FALSE,
    verbose = FALSE
  )
  expect_false("AGE" %in% result$metadata$selected_terms)
  expect_equal(result$metadata$adjust_for, c("AGE"))
})
