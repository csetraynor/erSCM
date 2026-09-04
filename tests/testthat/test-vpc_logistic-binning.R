## test-vpc_logistic-binning.R
## vpc_logistic() bin construction, ordering, and identity-join guard.

.fit_vpc_bin_scm_glm <- function(center = FALSE, scale_to_median = FALSE,
                                 var_spec = NULL) {
  data(demo_eres)
  demo_eres$SEX <- as.factor(demo_eres$SEX)
  suppressMessages(
    scm_glm(
      eres_data          = demo_eres,
      exposure_vars      = "CMAXss",
      covariates         = c("AGE", "WT"),
      dv_col             = "sim_status",
      significance_level = c(0.3, 0.3),
      var_spec           = var_spec,
      center             = center,
      scale_to_median    = scale_to_median,
      verbose            = FALSE
    )
  )
}

test_that("vpc_logistic() returns one row per requested bin, numerically ordered", {
  fit <- .fit_vpc_bin_scm_glm()
  skip_if_not("WT" %in% fit$metadata$finalRelation)

  d <- vpc_logistic(fit, "WT", n_bins = 6L, n_sim = 40L, plot = FALSE)

  expect_equal(nrow(d), 6L)
  # mean_x must increase monotonically across the ordered bins (never lexical).
  expect_false(is.unsorted(d$mean_x))
  expect_true(is.ordered(d$bin))
})

test_that("vpc_logistic() honours var_spec back-transformed labels under scaling", {
  spec <- parse_var_spec(c(
    CMAXss = "Cmax,ss (ng/mL)",
    WT     = "Body Weight (kg)",
    AGE    = "Age (years)"
  ))
  fit <- .fit_vpc_bin_scm_glm(center = FALSE, scale_to_median = TRUE,
                              var_spec = spec)
  skip_if_not("WT" %in% fit$metadata$finalRelation)

  # Plot builds without error and carries the var_spec x label.
  p <- vpc_logistic(fit, "WT", n_bins = 4L, n_sim = 40L)
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$x, "Body Weight")
})

test_that("vpc_logistic() ordering survives n_bins >= 10", {
  fit <- .fit_vpc_bin_scm_glm()
  skip_if_not("WT" %in% fit$metadata$finalRelation)

  d <- vpc_logistic(fit, "WT", n_bins = 10L, n_sim = 30L, plot = FALSE)
  expect_equal(nrow(d), 10L)
  expect_false(is.unsorted(d$mean_x))
})

test_that("vpc_logistic() observed rate for each bin comes from that bin's own rows", {
  # Identity-join regression guard (Phase A analogue): the observed proportion
  # attached to a bin must equal the mean outcome over exactly the rows that
  # fall in that bin -- not a lexically/positionally mismatched set.
  fit <- .fit_vpc_bin_scm_glm()
  skip_if_not("WT" %in% fit$metadata$finalRelation)

  mf  <- stats::model.frame(fit$final_model)
  yv  <- as.numeric(stats::model.response(mf))
  wt  <- mf[["WT"]]
  n_bins <- 5L
  bin <- ggplot2::cut_number(wt, n = n_bins)

  ref <- tapply(yv, bin, mean)
  ref <- ref[order(tapply(wt, bin, mean))]

  d <- vpc_logistic(fit, "WT", n_bins = n_bins, n_sim = 20L, plot = FALSE)

  expect_equal(unname(d$obs_rate), unname(as.numeric(ref)))
  # And the observed CI must bracket the observed rate it was computed from.
  expect_true(all(d$obs_lo <= d$obs_rate + 1e-9))
  expect_true(all(d$obs_hi >= d$obs_rate - 1e-9))
})

test_that("vpc_logistic() supports a categorical covariate with ordered levels", {
  # Synthetic data where SEX strongly predicts the outcome, so it reliably
  # enters the model and can be used as the VPC stratifier.
  set.seed(707L)
  n   <- 120L
  SEX <- factor(sample(c("F", "M"), n, replace = TRUE))
  EXP <- rnorm(n, 100, 20)
  DV  <- rbinom(n, 1L, plogis(-1.5 + 3 * (SEX == "M")))
  dat <- data.frame(DV = DV, EXP = EXP, SEX = SEX)
  fit <- suppressMessages(
    scm_glm(
      eres_data          = dat,
      exposure_vars      = "EXP",
      covariates         = "SEX",
      dv_col             = "DV",
      significance_level = c(0.8, 0.8),
      center             = FALSE,
      verbose            = FALSE
    )
  )
  skip_if_not("SEX" %in% fit$metadata$finalRelation)

  d <- vpc_logistic(fit, "SEX", n_sim = 30L, plot = FALSE)
  expect_true(is.ordered(d$bin))
  expect_equal(nrow(d), nlevels(SEX))
})
