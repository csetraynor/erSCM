## test-vpc_coxph_sim-validate.R
## vpc_coxph_sim() input validation and unsupported-model abort paths.

.fit_vpc_sim_valid_categorical <- function() {
  set.seed(912L)
  n <- 220L
  EXP <- factor(sample(c("Low", "High"), n, replace = TRUE),
                levels = c("Low", "High"), ordered = TRUE)
  SEX <- factor(sample(c("F", "M"), n, replace = TRUE),
                levels = c("F", "M"), ordered = TRUE)
  lp <- -0.6 + 1.0 * (EXP == "High") + 0.7 * (SEX == "M")
  dat <- data.frame(
    TIME = rexp(n, rate = exp(lp - 2.9)),
    DV = rbinom(n, 1L, plogis(lp - 0.9)),
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
      center = FALSE,
      verbose = FALSE
    )
  )
}

.fit_vpc_sim_valid_numeric <- function() {
  data(demo_eres)
  suppressWarnings(
    scm_coxph(
      eres_data = demo_eres,
      exposure_vars = "CMAXss",
      covariates = c("AGE", "WT"),
      dv_col = "sim_status",
      time_col = "sim_time",
      significance_level = c(0.3, 0.3),
      center = FALSE,
      verbose = FALSE
    )
  )
}

.as_eres_cox <- function(mod, data) {
  structure(
    list(
      final_model = mod,
      model_data = data,
      metadata = list(time_col = "TIME", dv_col = "DV", finalRelation = attr(mod$terms, "term.labels")),
      var_spec = NULL
    ),
    class = c("eres_scm_coxph", "eres_scm", "list")
  )
}

test_that("vpc_coxph_sim() rejects a non-eres_scm_coxph object", {
  expect_error(vpc_coxph_sim(list(1, 2), "EXP"), "eres_scm_coxph")
})

test_that("vpc_coxph_sim() rejects fit with no final_model", {
  bad_fit <- structure(
    list(
      final_model = NULL,
      metadata = list(time_col = "TIME", finalRelation = "1"),
      model_data = data.frame(TIME = numeric(0), DV = integer(0)),
      var_spec = NULL
    ),
    class = c("eres_scm_coxph", "eres_scm")
  )
  expect_error(vpc_coxph_sim(bad_fit, "EXP"), "final_model")
})

test_that("vpc_coxph_sim() rejects interaction and unknown covariate terms", {
  fit <- .fit_vpc_sim_valid_categorical()
  expect_error(vpc_coxph_sim(fit, "EXP:SEX"), "interaction")
  expect_error(vpc_coxph_sim(fit, "NOPE"), "predictor")
})

test_that("vpc_coxph_sim() validates n_bins, n_sim, probs, and time_grid", {
  fit <- .fit_vpc_sim_valid_categorical()

  expect_error(vpc_coxph_sim(fit, "EXP", n_bins = 0), "n_bins")
  expect_error(vpc_coxph_sim(fit, "EXP", n_sim = 0), "n_sim")
  expect_error(vpc_coxph_sim(fit, "EXP", n_sim = 2.5), "n_sim")
  expect_error(vpc_coxph_sim(fit, "EXP", probs = c(0.5)), "probs")
  expect_error(vpc_coxph_sim(fit, "EXP", probs = c(0, 1)), "probs")
  expect_error(vpc_coxph_sim(fit, "EXP", probs = c(0.1, 0.2, 0.3, 0.4)), "probs")
  expect_error(vpc_coxph_sim(fit, "EXP", time_grid = c(NA_real_, 1)), "time_grid")
})

test_that("vpc_coxph_sim() aborts numeric n_bins > 1 in v2", {
  fit <- .fit_vpc_sim_valid_numeric()
  expect_error(vpc_coxph_sim(fit, "WT", n_bins = 2, plot = FALSE), "S-6")
})

test_that("vpc_coxph_sim() aborts unsupported formula specials", {
  set.seed(1101L)
  dat <- data.frame(
    TIME = rexp(120, rate = 0.08),
    DV = rbinom(120, 1L, 0.6),
    x = rnorm(120),
    grp = factor(sample(c("A", "B"), 120, TRUE))
  )
  mod_strata <- survival::coxph(survival::Surv(TIME, DV) ~ x + strata(grp), data = dat)
  fit_strata <- .as_eres_cox(mod_strata, dat)
  expect_error(vpc_coxph_sim(fit_strata, "grp"), "does not support")

  dat_ss <- transform(dat, START = pmax(TIME - runif(nrow(dat), 0, 0.5), 0), STOP = TIME)
  mod_ss <- survival::coxph(survival::Surv(START, STOP, DV) ~ x, data = dat_ss)
  fit_ss <- .as_eres_cox(mod_ss, dat_ss)
  expect_error(vpc_coxph_sim(fit_ss, "x"), "start-stop")
})
