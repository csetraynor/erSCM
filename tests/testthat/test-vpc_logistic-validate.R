## test-vpc_logistic-validate.R
## vpc_logistic() input validation aborts (Phase E cli::cli_abort style).

.fit_vpc_valid_scm_glm <- function() {
  set.seed(303L)
  n <- 50L
  CMAX <- rnorm(n, 100, 20)
  AGE  <- rnorm(n, 55, 12)
  DV   <- rbinom(n, 1L, plogis(-6 + 0.06 * CMAX))
  test_data <- data.frame(DV = DV, CMAX = CMAX, AGE = AGE)
  suppressMessages(
    scm_glm(
      eres_data          = test_data,
      exposure_vars      = "CMAX",
      covariates         = "AGE",
      significance_level = c(0.5, 0.5),
      dv_col             = "DV",
      center             = FALSE,
      verbose            = FALSE
    )
  )
}

test_that("vpc_logistic() rejects a non-eres_scm_glm object", {
  skip_on_cran()
  expect_error(vpc_logistic(list(1, 2), "AGE"), "eres_scm_glm")
})

test_that("vpc_logistic() rejects a Cox fit", {
  skip_on_cran()
  data(demo_eres)
  demo_eres$SEX <- as.factor(demo_eres$SEX)
  cox <- suppressMessages(
    scm_coxph(
      eres_data = demo_eres, exposure_vars = "CMAXss", covariates = c("AGE", "WT"),
      dv_col = "sim_status", time_col = "sim_time",
      significance_level = c(0.3, 0.3), center = FALSE
    )
  )
  expect_error(vpc_logistic(cox, "WT"), "eres_scm_glm")
})

test_that("vpc_logistic() rejects interaction terms", {
  skip_on_cran()
  fit <- .fit_vpc_valid_scm_glm()
  expect_error(vpc_logistic(fit, "CMAX:AGE"), "interaction")
})

test_that("vpc_logistic() rejects an unknown covariate", {
  skip_on_cran()
  fit <- .fit_vpc_valid_scm_glm()
  expect_error(vpc_logistic(fit, "NOPE"), "predictor")
})

test_that("vpc_logistic() rejects invalid n_bins, n_sim, and probs", {
  skip_on_cran()
  fit <- .fit_vpc_valid_scm_glm()
  expect_error(vpc_logistic(fit, "AGE", n_bins = 1L), "n_bins")
  expect_error(vpc_logistic(fit, "AGE", n_sim = 0L), "n_sim")
  expect_error(vpc_logistic(fit, "AGE", n_sim = 2.5), "n_sim")
  expect_error(vpc_logistic(fit, "AGE", probs = c(0.5)), "probs")
  expect_error(vpc_logistic(fit, "AGE", probs = c(0, 1)), "probs")
  expect_error(vpc_logistic(fit, "AGE", probs = c(0.1, 0.2, 0.3, 0.4)), "probs")
})
