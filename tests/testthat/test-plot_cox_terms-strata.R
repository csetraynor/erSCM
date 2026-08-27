# Regression test for plot_cox_terms predicted-curve / strata alignment.
# Guards against the Finding #1 regression: predicted Cox curves were assigned
# to strata by lexical position rather than by bin/level identity, causing
# silent mislabelling when bin labels did not sort in numeric order (e.g. with
# var_spec back-transformed labels or n_bins >= 10).

test_that("plot_cox_terms attaches predicted curves to the correct strata (var_spec labels)", {
  skip_if_not_installed("survival")

  data(demo_eres)
  spec <- parse_var_spec(c(
    CMAXss = "Cmax,ss (ng/mL)",
    WT     = "Body Weight (kg)",
    AGE    = "Age (years)"
  ))
  fit <- scm_coxph(
    eres_data = demo_eres, exposure_vars = "CMAXss", covariates = c("AGE", "WT"),
    dv_col = "sim_status", time_col = "sim_time",
    significance_level = c(0.3, 0.3),
    var_spec = spec, scale_to_median = TRUE, center = FALSE
  )
  skip_if_not("WT" %in% fit$metadata$finalRelation)

  d <- plot_cox_terms(fit, "WT", n_bins = 4, plot = FALSE)

  # Within each facet, the predicted (Cox) curve and the observed KM curve must
  # describe the SAME weight bin. We verify this by an ordering invariant:
  # higher weight bins must have lower predicted survival at the last time point.
  ci_lab <- grep("^Cox model", unique(d$origin), value = TRUE)
  pred <- d[d$origin == ci_lab, ]
  strata_u <- unique(pred$strata)
  # Extract the first numeric edge of each bin from the strata label. Strip the
  # legend prefix up to the newline, then drop any leading non-numeric symbols
  # ("<=", ">", spaces) before the first number so labels like "<= 62" and
  # "> 80" parse robustly.
  bin_txt    <- sub(".*\n", "", strata_u)
  lower_edge <- suppressWarnings(as.numeric(
    sub("^[^0-9.]*([0-9.]+).*", "\\1", bin_txt)
  ))
  # Last-time predicted survival per stratum, kept in the SAME order as
  # strata_u / lower_edge (do not use split(), which reorders by sorted name).
  last_surv <- vapply(strata_u, function(st) {
    s <- pred[pred$strata == st, ]
    s$estimate[which.max(s$time)]
  }, numeric(1))
  ord <- order(lower_edge)
  # Monotone non-increasing survival with increasing weight bin (allow tie slack)
  expect_true(
    all(diff(last_surv[ord]) <= 1e-6),
    info = "Predicted survival is not monotone across ordered weight bins -- predicted curves are misaligned to strata."
  )
})
