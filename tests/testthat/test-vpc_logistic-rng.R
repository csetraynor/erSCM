## test-vpc_logistic-rng.R
## vpc_logistic() RNG parity and hygiene (seed-parity sibling of scm_bootstrap).
## The parametric-bootstrap band must be identical for the same seed regardless
## of parallel = FALSE vs parallel = TRUE, deterministic under a fixed seed, and
## must not leak the internal L'Ecuyer-CMRG switch into the caller's RNG state.

.fit_vpc_test_scm_glm <- function() {
  set.seed(202L)
  n <- 80L
  CMAX <- rnorm(n, 100, 20)
  AGE  <- rnorm(n, 55, 12)
  WT   <- rnorm(n, 75, 15)
  # DV depends on AGE so that AGE reliably enters the model and is a valid
  # covariate_term for the VPC.
  lp   <- -4 + 0.10 * (AGE - 55)
  DV   <- rbinom(n, 1L, plogis(lp))
  test_data <- data.frame(
    DV = DV, CMAX = CMAX, AGE = AGE, WT = WT,
    stringsAsFactors = FALSE
  )
  suppressMessages(
    scm_glm(
      eres_data          = test_data,
      exposure_vars      = "CMAX",
      covariates         = c("AGE", "WT"),
      significance_level = c(0.8, 0.8),
      dv_col             = "DV",
      family             = "binomial",
      center             = FALSE,
      scale_to_median    = FALSE,
      verbose            = FALSE
    )
  )
}

test_that("vpc_logistic() band is identical for parallel = FALSE vs TRUE", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

  fit <- .fit_vpc_test_scm_glm()

  seq_run <- vpc_logistic(fit, "AGE", n_sim = 40L, seed = 20260831L,
                          parallel = FALSE, plot = FALSE)
  par_run <- vpc_logistic(fit, "AGE", n_sim = 40L, seed = 20260831L,
                          parallel = TRUE, plot = FALSE)

  expect_equal(seq_run$sim_lo, par_run$sim_lo)
  expect_equal(seq_run$sim_mid, par_run$sim_mid)
  expect_equal(seq_run$sim_hi, par_run$sim_hi)
})

test_that("vpc_logistic() is deterministic for a fixed seed", {
  fit <- .fit_vpc_test_scm_glm()

  a <- vpc_logistic(fit, "AGE", n_sim = 40L, seed = 99L, plot = FALSE)
  b <- vpc_logistic(fit, "AGE", n_sim = 40L, seed = 99L, plot = FALSE)
  expect_identical(a, b)
})

test_that("vpc_logistic() does not leak its RNG kind or seed into the caller", {
  fit <- .fit_vpc_test_scm_glm()

  RNGkind("Mersenne-Twister")
  set.seed(12345L)
  before_kind <- RNGkind()
  before_seed <- .Random.seed

  invisible(vpc_logistic(fit, "AGE", n_sim = 30L, seed = 7L, plot = FALSE))

  expect_identical(RNGkind(), before_kind)
  expect_identical(.Random.seed, before_seed)
})
