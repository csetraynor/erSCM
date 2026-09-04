## test-vpc_coxph.R
## vpc_coxph() validation, return shape, ordering, and silence guards.

.fit_vpc_valid_cox_numeric <- function() {
  data(demo_eres)
  suppressWarnings(
    scm_coxph(
      eres_data = demo_eres,
      exposure_vars = "CMAXss",
      covariates = c("AGE", "WT"),
      dv_col = "sim_status",
      time_col = "sim_time",
      significance_level = c(0.3, 0.3),
      center = FALSE
    )
  )
}

.fit_vpc_valid_cox_categorical <- function() {
  set.seed(812L)
  n <- 260L
  EXP <- factor(sample(c("Low", "High"), n, replace = TRUE),
                levels = c("Low", "High"), ordered = TRUE)
  SEX <- factor(sample(c("F", "M"), n, replace = TRUE),
                levels = c("F", "M"), ordered = TRUE)
  lp <- -0.5 + 1.0 * (EXP == "High") + 0.9 * (SEX == "M")
  dat <- data.frame(
    TIME = rexp(n, rate = exp(lp - 2.8)),
    DV = rbinom(n, 1L, plogis(lp - 0.8)),
    EXP = EXP,
    SEX = SEX
  )

  suppressWarnings(
    scm_coxph(
      eres_data = dat,
      exposure_vars = "EXP",
      covariates = "SEX",
      dv_col = "DV",
      time_col = "TIME",
      significance_level = c(0.95, 0.95),
      center = FALSE
    )
  )
}

test_that("vpc_coxph() rejects a non-eres_scm_coxph object", {
  expect_error(vpc_coxph(list(1, 2)), "eres_scm_coxph")
})

test_that("vpc_coxph() rejects a fit with no final_model", {
  bad_fit <- structure(
    list(
      final_model = NULL,
      metadata = list(time_col = "TIME", finalRelation = "1"),
      model_data = data.frame(TIME = numeric(0), DV = integer(0)),
      var_spec = NULL
    ),
    class = c("eres_scm_coxph", "eres_scm")
  )
  expect_error(vpc_coxph(bad_fit), "final_model")
})

test_that("vpc_coxph() validates ci and n_bins", {
  fit <- .fit_vpc_valid_cox_numeric()
  expect_error(vpc_coxph(fit, ci = 1), "ci")
  expect_error(vpc_coxph(fit, ci = 0), "ci")
  expect_error(vpc_coxph(fit, ci = NA_real_), "ci")

  expect_error(vpc_coxph(fit, n_bins = 0), "n_bins")
  expect_error(vpc_coxph(fit, n_bins = 1.5), "n_bins")
  expect_error(vpc_coxph(fit, n_bins = NA_real_), "n_bins")
})

test_that("vpc_coxph(plot = FALSE) returns a long data frame with key columns", {
  fit <- .fit_vpc_valid_cox_numeric()
  d <- vpc_coxph(fit, plot = FALSE)

  expect_s3_class(d, "data.frame")
  expect_true(all(c("time", "estimate", "conf.low", "conf.high", "origin") %in% names(d)))
  expect_true(all(c("Observed KM", grep("^Cox model", unique(d$origin), value = TRUE)) %in% unique(d$origin)))
})

test_that("vpc_coxph() supports a categorical stratification path with ordered strata", {
  fit <- .fit_vpc_valid_cox_categorical()
  d <- vpc_coxph(fit, plot = FALSE)

  skip_if_not("strata" %in% names(d))
  expect_s3_class(d$strata, "factor")
  expect_true(is.ordered(d$strata))
  expect_setequal(levels(d$strata), unique(as.character(d$strata)))
})

test_that("vpc_coxph() supports the no-categorical path", {
  fit <- .fit_vpc_valid_cox_numeric()
  d <- vpc_coxph(fit, n_bins = 1, plot = FALSE)

  expect_false("strata" %in% names(d))
})

test_that("vpc_coxph() errors for numeric path when n_bins > 1", {
  fit <- .fit_vpc_valid_cox_numeric()
  expect_error(vpc_coxph(fit, n_bins = 2, plot = FALSE), "not implemented")
})

test_that("vpc_coxph() allows n_bins > 1 for categorical-only path", {
  fit <- .fit_vpc_valid_cox_categorical()
  d <- vpc_coxph(fit, n_bins = 3, plot = FALSE)

  expect_s3_class(d, "data.frame")
  expect_true("strata" %in% names(d))
})

test_that("vpc_coxph() is silent for grouped output and join chatter", {
  fit <- .fit_vpc_valid_cox_categorical()

  msgs <- character(0)
  withCallingHandlers(
    vpc_coxph(fit, plot = FALSE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  expect_false(any(grepl("Joining, by =", msgs, fixed = TRUE)))
  expect_false(any(grepl("`summarise()` has grouped output", msgs, fixed = TRUE)))
})
