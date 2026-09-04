## test-vpc_coxph_sim-rng.R
## vpc_coxph_sim() reproducibility parity and RNG hygiene.

.fit_vpc_sim_rng <- function() {
  set.seed(777L)
  n <- 180L
  EXP <- factor(sample(c("Low", "High"), n, replace = TRUE),
                levels = c("Low", "High"), ordered = TRUE)
  SEX <- factor(sample(c("F", "M"), n, replace = TRUE),
                levels = c("F", "M"), ordered = TRUE)
  lp <- -0.4 + 0.9 * (EXP == "High") + 0.5 * (SEX == "M")
  dat <- data.frame(
    TIME = rexp(n, rate = exp(lp - 3.1)),
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
      center = FALSE,
      verbose = FALSE
    )
  )
}

test_that("vpc_coxph_sim() is deterministic for a fixed seed", {
  fit <- .fit_vpc_sim_rng()
  a <- vpc_coxph_sim(fit, "EXP", n_sim = 40L, seed = 99L, plot = FALSE)
  b <- vpc_coxph_sim(fit, "EXP", n_sim = 40L, seed = 99L, plot = FALSE)
  expect_identical(a, b)
})

test_that("vpc_coxph_sim() band is identical for parallel false vs true", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

  fit <- .fit_vpc_sim_rng()
  seq_run <- vpc_coxph_sim(fit, "EXP", n_sim = 40L, seed = 20260901L,
                           parallel = FALSE, plot = FALSE)
  par_run <- vpc_coxph_sim(fit, "EXP", n_sim = 40L, seed = 20260901L,
                           parallel = TRUE, plot = FALSE)

  expect_equal(seq_run$sim_lo, par_run$sim_lo)
  expect_equal(seq_run$sim_mid, par_run$sim_mid)
  expect_equal(seq_run$sim_hi, par_run$sim_hi)
})

test_that("vpc_coxph_sim() does not leak RNG kind or seed", {
  fit <- .fit_vpc_sim_rng()

  RNGkind("Mersenne-Twister")
  set.seed(12345L)
  before_kind <- RNGkind()
  before_seed <- .Random.seed

  invisible(vpc_coxph_sim(fit, "EXP", n_sim = 30L, seed = 7L, plot = FALSE))

  expect_identical(RNGkind(), before_kind)
  expect_identical(.Random.seed, before_seed)
})
