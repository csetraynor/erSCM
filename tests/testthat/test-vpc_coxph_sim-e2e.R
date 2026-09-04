## test-vpc_coxph_sim-e2e.R
## End-to-end checks, including the approved median-vs-analytic tolerance spot check.

.fit_vpc_sim_e2e <- function() {
  set.seed(2205L)
  n <- 260L
  TRT <- factor(sample(c("Control", "Drug"), n, replace = TRUE),
                levels = c("Control", "Drug"), ordered = TRUE)
  AGEG <- factor(sample(c("<65", ">=65"), n, replace = TRUE),
                 levels = c("<65", ">=65"), ordered = TRUE)
  lp <- -0.3 + 0.8 * (TRT == "Drug") + 0.5 * (AGEG == ">=65")
  dat <- data.frame(
    TIME = rexp(n, rate = exp(lp - 3.2)),
    DV = rbinom(n, 1L, plogis(lp - 0.8)),
    TRT = TRT,
    AGEG = AGEG
  )

  suppressWarnings(
    scm_coxph(
      eres_data = dat,
      exposure_vars = "TRT",
      covariates = "AGEG",
      dv_col = "DV",
      time_col = "TIME",
      significance_level = c(0.95, 0.95),
      center = FALSE,
      verbose = FALSE
    )
  )
}

test_that("vpc_coxph_sim() produces ordered, bounded predictive bands", {
  fit <- .fit_vpc_sim_e2e()
  d <- vpc_coxph_sim(fit, "TRT", n_sim = 80L, plot = FALSE)

  expect_true(all(d$sim_lo <= d$sim_hi + 1e-9))
  expect_true(all(d$sim_lo >= 0 & d$sim_hi <= 1))
  expect_true(all(d$obs_est >= 0 & d$obs_est <= 1))
})

test_that("vpc_coxph_sim median tracks analytic reference within tolerance", {
  skip_on_cran()

  set.seed(3306L)
  n <- 500L
  # Near-constant design and no censoring make the reference-profile comparison stable.
  EXP <- rnorm(n, mean = 100, sd = 0.005)
  TIME <- rexp(n, rate = 0.08 * exp(0.01 * (EXP - 100)))
  DV <- rep(1L, n)
  dat <- data.frame(TIME = TIME, DV = DV, EXP = EXP)

  mod <- survival::coxph(survival::Surv(TIME, DV) ~ EXP, data = dat)
  fit <- structure(
    list(
      final_model = mod,
      model_data = dat,
      metadata = list(time_col = "TIME", dv_col = "DV", finalRelation = "EXP"),
      var_spec = NULL
    ),
    class = c("eres_scm_coxph", "eres_scm", "list")
  )

  sim_df <- vpc_coxph_sim(fit, "EXP", n_bins = 1L, n_sim = 2000L,
                          probs = c(0.05, 0.5, 0.95), seed = 401L, plot = FALSE)

  analytic <- survival::survfit(fit$final_model,
                                newdata = data.frame(EXP = mean(EXP)))
  analytic_est <- {
    t <- c(0, analytic$time)
    s <- c(1, analytic$surv)
    idx <- findInterval(sim_df$time, t)
    idx[idx < 1L] <- 1L
    s[idx]
  }

  # Acceptance criterion approved for S-5: max pointwise abs diff <= 0.05.
  max_abs_diff <- max(abs(sim_df$sim_mid - analytic_est), na.rm = TRUE)
  expect_lte(max_abs_diff, 0.05)
})
