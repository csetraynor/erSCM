# helper-bootstrap.R
# Shared fixtures for bootstrap tests.
# testthat auto-loads files named helper-*.R before running any test file.

# Small synthetic dataset for GLM testing (n = 200)
.make_test_data_glm <- function(n = 200, seed = 123) {
  set.seed(seed)
  data.frame(
    DV     = rbinom(n, 1L, 0.3),
    CMAX   = rnorm(n, 100, 20),
    AUC    = rnorm(n, 500, 100),
    AGE    = rnorm(n, 55, 12),
    WT     = rnorm(n, 75, 15),
    SEX    = factor(sample(c("M", "F"), n, replace = TRUE)),
    RACE   = factor(sample(c("White", "Black", "Asian"), n, replace = TRUE,
                           prob = c(0.6, 0.2, 0.2))),
    SUBJID = paste0("SUBJ", seq_len(n)),
    stringsAsFactors = FALSE
  )
}

# Small synthetic dataset for CoxPH testing (n = 200)
.make_test_data_coxph <- function(n = 200, seed = 456) {
  set.seed(seed)
  data.frame(
    TIME   = rexp(n, rate = 0.05),
    DV     = rbinom(n, 1L, 0.4),
    CMAX   = rnorm(n, 100, 20),
    AUC    = rnorm(n, 500, 100),
    AGE    = rnorm(n, 55, 12),
    WT     = rnorm(n, 75, 15),
    SEX    = factor(sample(c("M", "F"), n, replace = TRUE)),
    SUBJID = paste0("SUBJ", seq_len(n)),
    stringsAsFactors = FALSE
  )
}

# Fit a GLM SCM model for testing (quick, small data).
# Use relaxed significance_level so at least one covariate enters.
.fit_test_scm_glm <- function() {
  test_data <- .make_test_data_glm()
  suppressMessages(
    scm_glm(
      eres_data         = test_data,
      exposure_vars     = "CMAX",
      covariates        = c("AGE", "WT", "SEX"),
      significance_level = c(0.5, 0.3),   # relaxed to ensure entry with small n
      dv_col            = "DV",
      family            = "binomial",
      center            = FALSE,
      scale_to_median   = FALSE,
      verbose           = FALSE
    )
  )
}

# Fit a GLM SCM model with multiple exposure candidates.
.fit_test_scm_glm_multi_exposure <- function() {
  test_data <- .make_test_data_glm()
  suppressMessages(
    scm_glm(
      eres_data         = test_data,
      exposure_vars     = c("CMAX", "AUC"),
      covariates        = c("AGE", "WT", "SEX"),
      significance_level = c(0.5, 0.3),   # relaxed to ensure entry with small n
      dv_col            = "DV",
      family            = "binomial",
      center            = FALSE,
      scale_to_median   = FALSE,
      verbose           = FALSE
    )
  )
}

# Fit a CoxPH SCM model for testing.
.fit_test_scm_coxph <- function() {
  test_data <- .make_test_data_coxph()
  suppressMessages(
    scm_coxph(
      eres_data         = test_data,
      exposure_vars     = "CMAX",
      covariates        = c("AGE", "WT", "SEX"),
      significance_level = c(0.5, 0.3),   # relaxed to ensure entry with small n
      dv_col            = "DV",
      time_col          = "TIME",
      center            = FALSE,
      scale_to_median   = FALSE,
      verbose           = FALSE
    )
  )
}
