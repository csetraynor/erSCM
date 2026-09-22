## test-vpc_logistic-e2e.R
## vpc_logistic() end-to-end on a small simulated binary dataset, plus the
## length-2 probs (no median line) path.

test_that("vpc_logistic() produces a well-formed band end-to-end", {
  skip_on_cran()
  set.seed(404L)
  n <- 120L
  x  <- rnorm(n, 0, 1)
  DV <- rbinom(n, 1L, plogis(-0.2 + 1.1 * x))
  dat <- data.frame(DV = DV, EXP = x, COV = rnorm(n))

  fit <- suppressMessages(
    scm_glm(
      eres_data          = dat,
      exposure_vars      = "EXP",
      covariates         = "COV",
      dv_col             = "DV",
      significance_level = c(0.5, 0.5),
      center             = FALSE,
      verbose            = FALSE
    )
  )

  d <- vpc_logistic(fit, "EXP", n_bins = 4L, n_sim = 60L, plot = FALSE)

  expect_equal(nrow(d), 4L)
  expect_true(all(c("sim_lo", "sim_mid", "sim_hi", "obs_rate") %in% names(d)))
  # Band is well-ordered and within [0, 1].
  expect_true(all(d$sim_lo <= d$sim_mid + 1e-9))
  expect_true(all(d$sim_mid <= d$sim_hi + 1e-9))
  expect_true(all(d$sim_lo >= 0 & d$sim_hi <= 1))
  expect_true(all(d$obs_rate >= 0 & d$obs_rate <= 1))
  expect_true(all(d$n > 0))
})

test_that("vpc_logistic() with length-2 probs omits the median (sim_mid = NA)", {
  skip_on_cran()
  set.seed(505L)
  n <- 100L
  x  <- rnorm(n)
  DV <- rbinom(n, 1L, plogis(0.5 * x))
  dat <- data.frame(DV = DV, EXP = x, COV = rnorm(n))
  fit <- suppressMessages(
    scm_glm(
      eres_data = dat, exposure_vars = "EXP", covariates = "COV",
      dv_col = "DV", significance_level = c(0.5, 0.5),
      center = FALSE, verbose = FALSE
    )
  )

  d <- vpc_logistic(fit, "EXP", n_bins = 3L, n_sim = 40L,
                    probs = c(0.1, 0.9), plot = FALSE)
  expect_true(all(is.na(d$sim_mid)))
  expect_true(all(d$sim_lo <= d$sim_hi + 1e-9))

  p <- vpc_logistic(fit, "EXP", n_bins = 3L, n_sim = 40L, probs = c(0.1, 0.9))
  expect_s3_class(p, "ggplot")
})
