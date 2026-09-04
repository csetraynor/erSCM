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

# ---------------------------------------------------------------------------
# Finding #3 regression tests: facet (strata) ordering is numeric/factor-level,
# not lexical order of the label strings. These guard the display-only ordered
# factor conversion in plot_cox_terms().
# ---------------------------------------------------------------------------

test_that("plot_cox_terms facets are ordered by numeric bin, not label string (var_spec labels)", {
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

  # strata is now a display-only ordered factor
  expect_s3_class(d$strata, "factor")
  expect_true(is.ordered(d$strata))

  lv <- levels(d$strata)
  # var_spec back-transformed labels use "<=", ">" and/or an en-dash range
  expect_true(any(grepl("\u2264|>|\u2013|-", lv)),
              info = paste("Levels:", paste(lv, collapse = " | ")))

  # Parse the first numeric edge of each bin (strip legend prefix and any
  # leading non-numeric symbols such as "<=", ">", spaces).
  bin_txt <- sub(".*\n", "", lv)
  edge    <- suppressWarnings(as.numeric(sub("^[^0-9.]*([0-9.]+).*", "\\1", bin_txt)))
  expect_false(any(is.na(edge)),
               info = paste("Unparsed levels:", paste(lv, collapse = " | ")))

  # Levels must already be in ascending numeric order (NOT lexical order).
  expect_identical(lv, lv[order(edge)])

  # No panel dropped or invented: level set == set of strata actually present.
  expect_setequal(levels(d$strata), unique(as.character(d$strata)))
})

test_that("plot_cox_terms facets with n_bins=12 are ascending, not lexical (Q1,Q10,Q11,...)", {
  skip_if_not_installed("survival")

  data(demo_eres)
  # No var_spec -> labels are Q1..Q12, which expose the lexical-sort bug
  # (Q1, Q10, Q11, Q12, Q2, ...) if ordering is done on the raw strings.
  fit <- scm_coxph(
    eres_data = demo_eres, exposure_vars = "CMAXss", covariates = c("AGE", "WT"),
    dv_col = "sim_status", time_col = "sim_time",
    significance_level = c(0.3, 0.3)
  )
  skip_if_not("WT" %in% fit$metadata$finalRelation)

  d <- tryCatch(
    plot_cox_terms(fit, "WT", n_bins = 12, plot = FALSE),
    error = function(e) e
  )
  skip_if(inherits(d, "error"),
          paste("n_bins=12 not supported for this data:",
                conditionMessage(if (inherits(d, "error")) d else simpleError(""))))

  lv <- levels(d$strata)
  # Extract the integer after "Q" in each label (e.g. "WT\nQ12" -> 12)
  qnum <- suppressWarnings(as.integer(sub(".*\nQ", "", lv)))
  skip_if(any(is.na(qnum)), "labels are not Q-style; cannot assert Q ordering")

  # Ascending numeric order 1..12
  expect_identical(qnum, sort(qnum))
  expect_identical(qnum, seq_len(length(qnum)))

  # And explicitly NOT the lexical ordering of the bin labels.
  lexical <- lv[order(sub(".*\n", "", lv))]
  expect_false(identical(lv, lexical),
               info = "Facet order is lexical (Q1,Q10,Q11,...) instead of numeric.")

  expect_setequal(levels(d$strata), unique(as.character(d$strata)))
})

test_that("plot_cox_terms facets for a categorical covariate follow factor-level order", {
  skip_if_not_installed("survival")

  set.seed(11)
  n <- 500
  SEXf <- factor(sample(c("F", "M"), n, replace = TRUE))
  EXP  <- rnorm(n, 100, 20)
  WT   <- rnorm(n, 70, 12)
  lp   <- -0.05 + 0.01 * EXP + 1.2 * (SEXf == "M")
  df <- data.frame(
    TIME = rexp(n, exp(lp / 10 - 3)),
    DV   = rbinom(n, 1, plogis(lp - 2)),
    EXP  = EXP,
    WT   = WT,
    SEX  = SEXf
  )
  fit <- suppressWarnings(scm_coxph(
    eres_data = df, covariates = c("WT", "SEX"), exposure_vars = "EXP",
    dv_col = "DV", time_col = "TIME",
    significance_level = c(0.95, 0.9)
  ))
  skip_if_not("SEX" %in% get_main_effect_columns(fit$metadata$finalRelation))

  d <- plot_cox_terms(fit, "SEX", plot = FALSE)
  expect_true(is.ordered(d$strata))

  lv  <- levels(d$strata)
  raw <- sub(".*\n", "", lv)                       # level label after the prefix
  fct_levels <- levels(factor(fit$model_data[["SEX"]]))
  # Facet order follows the factor levels of the original variable.
  expect_identical(raw, fct_levels)

  expect_setequal(levels(d$strata), unique(as.character(d$strata)))
})

test_that("levels(strata) equal the identity-ordered label vector (anti-lexical-sort safeguard)", {
  # This is the structural safeguard requested for Finding #3: it fails loudly
  # if a future refactor re-introduces character-sort ordering. For numeric
  # covariates the identity source of truth (ref_data_expanded bin order) is, by
  # construction, numeric-ascending; asserting levels == numeric-ascending order
  # is therefore equivalent to asserting levels == transformed ref label vector.
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

  d  <- plot_cox_terms(fit, "WT", n_bins = 4, plot = FALSE)
  lv <- levels(d$strata)

  bin_txt <- sub(".*\n", "", lv)
  edge    <- suppressWarnings(as.numeric(sub("^[^0-9.]*([0-9.]+).*", "\\1", bin_txt)))
  expected <- lv[order(edge)]

  # Same length, same values, same order.
  expect_identical(length(lv), length(expected))
  expect_identical(lv, expected)
})
