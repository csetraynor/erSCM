## test-vpc_coxph_sim-shape.R
## vpc_coxph_sim() return shape, time grid handling, and strata ordering.

.fit_vpc_sim_shape <- function() {
  set.seed(1404L)
  n <- 240L
  ARM <- factor(sample(c("Placebo", "Active"), n, replace = TRUE),
                levels = c("Placebo", "Active"), ordered = TRUE)
  SEX <- factor(sample(c("F", "M"), n, replace = TRUE),
                levels = c("F", "M"), ordered = TRUE)
  lp <- -0.2 + 0.8 * (ARM == "Active") + 0.6 * (SEX == "M")
  dat <- data.frame(
    TIME = rexp(n, rate = exp(lp - 3.0)),
    DV = rbinom(n, 1L, plogis(lp - 0.7)),
    ARM = ARM,
    SEX = SEX
  )

  suppressWarnings(
    scm_coxph(
      eres_data = dat,
      exposure_vars = "ARM",
      covariates = "SEX",
      dv_col = "DV",
      time_col = "TIME",
      significance_level = c(0.95, 0.95),
      center = FALSE,
      verbose = FALSE
    )
  )
}

test_that("vpc_coxph_sim(plot = FALSE) returns required columns", {
  skip_on_cran()
  fit <- .fit_vpc_sim_shape()
  d <- vpc_coxph_sim(fit, "ARM", n_sim = 60L, plot = FALSE)

  expect_s3_class(d, "data.frame")
  expect_true(all(c("strata", "time", "sim_lo", "sim_mid", "sim_hi", "obs_est", "obs_lo", "obs_hi") %in% names(d)))
  expect_s3_class(d$strata, "factor")
  expect_true(is.ordered(d$strata))
  expect_true(all(d$time >= 0))
})

test_that("vpc_coxph_sim() clamps user time_grid to observed event range", {
  skip_on_cran()
  fit <- .fit_vpc_sim_shape()
  max_evt <- max(fit$model_data$TIME[fit$model_data$DV == 1])

  d <- vpc_coxph_sim(
    fit,
    "ARM",
    n_sim = 30L,
    time_grid = c(-10, 0, 0.5, 2, 1e6),
    plot = FALSE
  )

  expect_true(all(d$time >= 0))
  expect_true(all(d$time <= max_evt + 1e-12))
})

test_that("vpc_coxph_sim() returns ggplot when plot = TRUE", {
  skip_on_cran()
  fit <- .fit_vpc_sim_shape()
  p <- vpc_coxph_sim(fit, "ARM", n_sim = 40L, plot = TRUE)
  expect_s3_class(p, "ggplot")
})
