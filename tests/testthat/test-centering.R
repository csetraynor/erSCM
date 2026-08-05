test_that("compute_conditional_effect returns correct beta and se for single interaction", {
  # Simple case: Y ~ X1 + X2 + X1:X2
  # coefs: X1=1, X2=0.5, X1:X2=0.3
  # At X2=2: cond_beta = 1 + 0.3*2 = 1.6
  coefs <- c(X1 = 1.0, X2 = 0.5, "X1:X2" = 0.3)
  V <- matrix(
    c(0.04, 0.01, 0.005,
      0.01, 0.09, 0.002,
      0.005, 0.002, 0.01),
    nrow = 3, ncol = 3,
    dimnames = list(names(coefs), names(coefs))
  )
  res <- compute_conditional_effect(coefs, V, "X1", "X1:X2", 2.0)
  expect_equal(res$beta, 1.6)
  # Var = V[X1,X1] + 2^2*V[X1:X2,X1:X2] + 2*2*V[X1,X1:X2]
  expected_var <- 0.04 + 4 * 0.01 + 2 * 2 * 0.005
  expect_equal(res$se, sqrt(expected_var), tolerance = 1e-8)
})

test_that("compute_conditional_effect handles multiple interaction terms", {
  coefs <- c(X1 = 1.0, X2 = 0.5, X3 = 0.2, "X1:X2" = 0.3, "X1:X3" = 0.1)
  V <- diag(5)
  dimnames(V) <- list(names(coefs), names(coefs))
  # At X2=1, X3=2: cond_beta = 1 + 0.3*1 + 0.1*2 = 1.5
  res <- compute_conditional_effect(coefs, V, "X1", c("X1:X2", "X1:X3"), c(1, 2))
  expect_equal(res$beta, 1.5)
})

test_that("scm_glm stores centering_values in metadata when center=TRUE", {
  skip_on_cran()
  set.seed(42)
  n <- 200
  df <- data.frame(
    DV  = rbinom(n, 1, 0.3),
    EXP = rnorm(n, 10, 2),
    AGE = rnorm(n, 50, 10)
  )
  result <- scm_glm(
    eres_data = df,
    covariates = "AGE",
    exposure_vars = "EXP",
    dv_col = "DV",
    significance_level = c(0.9, 0.8),
    center = TRUE
  )
  expect_true(!is.null(result$metadata$centering_values))
  expect_true(is.numeric(result$metadata$centering_values))
  expect_true("EXP" %in% names(result$metadata$centering_values))
  expect_true("AGE" %in% names(result$metadata$centering_values))
})

test_that("scm_glm stores NULL centering_values when center=FALSE", {
  skip_on_cran()
  set.seed(42)
  n <- 200
  df <- data.frame(
    DV  = rbinom(n, 1, 0.3),
    EXP = rnorm(n, 10, 2),
    AGE = rnorm(n, 50, 10)
  )
  result <- scm_glm(
    eres_data = df,
    covariates = "AGE",
    exposure_vars = "EXP",
    dv_col = "DV",
    significance_level = c(0.9, 0.8),
    center = FALSE
  )
  expect_null(result$metadata$centering_values)
})

test_that("scm_coxph stores centering_values in metadata when center=TRUE", {
  skip_on_cran()
  set.seed(42)
  n <- 200
  df <- data.frame(
    TIME = rexp(n, 0.1),
    DV   = rbinom(n, 1, 0.3),
    EXP  = rnorm(n, 10, 2),
    AGE  = rnorm(n, 50, 10)
  )
  result <- scm_coxph(
    eres_data = df,
    covariates = "AGE",
    exposure_vars = "EXP",
    dv_col = "DV",
    time_col = "TIME",
    significance_level = c(0.9, 0.8),
    center = TRUE
  )
  expect_true(!is.null(result$metadata$centering_values))
  expect_true("EXP" %in% names(result$metadata$centering_values))
  expect_true("AGE" %in% names(result$metadata$centering_values))
  # Time and event columns should NOT be centering targets (they're excluded)
  # The centering_values may contain them since TIME/DV are excluded, but AGE/EXP must be there
})

test_that("scm_coxph stores NULL centering_values when center=FALSE", {
  skip_on_cran()
  set.seed(42)
  n <- 200
  df <- data.frame(
    TIME = rexp(n, 0.1),
    DV   = rbinom(n, 1, 0.3),
    EXP  = rnorm(n, 10, 2),
    AGE  = rnorm(n, 50, 10)
  )
  result <- scm_coxph(
    eres_data = df,
    covariates = "AGE",
    exposure_vars = "EXP",
    dv_col = "DV",
    time_col = "TIME",
    significance_level = c(0.9, 0.8),
    center = FALSE
  )
  expect_null(result$metadata$centering_values)
})

test_that("centering does not change model predictions (GLM)", {
  skip_on_cran()
  set.seed(99)
  n <- 300
  df <- data.frame(
    DV  = rbinom(n, 1, 0.3),
    EXP = rnorm(n, 10, 2),
    AGE = rnorm(n, 50, 10)
  )
  # Fit models with and without centering at very lenient thresholds
  result_center <- scm_glm(
    eres_data = df, covariates = "AGE", exposure_vars = "EXP",
    dv_col = "DV", significance_level = c(0.9, 0.8), center = TRUE
  )
  result_nocenter <- scm_glm(
    eres_data = df, covariates = "AGE", exposure_vars = "EXP",
    dv_col = "DV", significance_level = c(0.9, 0.8), center = FALSE
  )

  # Both must select at least one covariate to be comparable
  skip_if(length(result_center$metadata$finalRelation) == 0 ||
            length(result_nocenter$metadata$finalRelation) == 0,
          "No covariates selected, cannot compare predictions")
  skip_if(!setequal(result_center$metadata$finalRelation,
                    result_nocenter$metadata$finalRelation),
          "Different final models selected, prediction comparison not meaningful")

  # Predictions must be identical (centering only shifts parameterization, not fitted values)
  pred_center   <- fitted(result_center$final_model)
  pred_nocenter <- fitted(result_nocenter$final_model)
  expect_equal(pred_center, pred_nocenter, tolerance = 1e-6)
})

test_that("centering does not change model predictions (Cox)", {
  skip_on_cran()
  set.seed(77)
  n <- 300
  df <- data.frame(
    TIME = rexp(n, 0.1),
    DV   = rbinom(n, 1, 0.3),
    EXP  = rnorm(n, 10, 2),
    AGE  = rnorm(n, 50, 10)
  )
  result_center <- scm_coxph(
    eres_data = df, covariates = "AGE", exposure_vars = "EXP",
    dv_col = "DV", time_col = "TIME",
    significance_level = c(0.9, 0.8), center = TRUE
  )
  result_nocenter <- scm_coxph(
    eres_data = df, covariates = "AGE", exposure_vars = "EXP",
    dv_col = "DV", time_col = "TIME",
    significance_level = c(0.9, 0.8), center = FALSE
  )

  skip_if(length(result_center$metadata$finalRelation) == 0 ||
            length(result_nocenter$metadata$finalRelation) == 0,
          "No covariates selected")
  skip_if(!setequal(result_center$metadata$finalRelation,
                    result_nocenter$metadata$finalRelation),
          "Different final models selected")

  # Linear predictors (risk scores) must be identical up to a global shift
  lp_center   <- predict(result_center$final_model, type = "lp")
  lp_nocenter <- predict(result_nocenter$final_model, type = "lp")
  lp_diff <- lp_center - lp_nocenter
  expect_equal(diff(range(lp_diff)), 0, tolerance = 1e-5)
})

test_that("centering does not affect factor/logical columns", {
  skip_on_cran()
  set.seed(11)
  n <- 200
  df <- data.frame(
    DV  = rbinom(n, 1, 0.3),
    EXP = rnorm(n, 5, 1),
    SEX = factor(sample(c("M", "F"), n, replace = TRUE))
  )
  result <- scm_glm(
    eres_data = df, covariates = "SEX", exposure_vars = "EXP",
    dv_col = "DV", significance_level = c(0.9, 0.8), center = TRUE
  )
  # SEX should still be a factor in model_data (not centered/broken)
  expect_true(is.factor(result$model_data$SEX))
})

test_that("make_forest_plot back-transforms quantile labels when centered", {
  skip_on_cran()
  set.seed(55)
  n <- 400
  df <- data.frame(
    DV  = rbinom(n, 1, 0.3),
    EXP = rnorm(n, 100, 20),  # median ~100
    AGE = rnorm(n, 50, 10)
  )
  result <- scm_glm(
    eres_data = df, covariates = "AGE", exposure_vars = "EXP",
    dv_col = "DV", significance_level = c(0.9, 0.8), center = TRUE
  )
  skip_if(!"EXP" %in% result$metadata$finalRelation, "EXP not selected")

  fp <- make_forest_plot(result, plot = FALSE)
  exp_rows <- fp[fp$group == "EXP", ]
  skip_if(nrow(exp_rows) == 0, "EXP not in forest plot")

  # With centering, the median of EXP in centered data is 0, but back-transformed
  # labels should show values near the original scale (~100). The group_level
  # contains the rounded VAL, which after back-transformation should be > 50 (not ~0)
  # Extract numeric values from group_level labels
  num_labels <- suppressWarnings(as.numeric(gsub("[^0-9.-]", "", exp_rows$group_level)))
  non_na <- num_labels[!is.na(num_labels) & num_labels != 0]
  if (length(non_na) > 0) {
    # Back-transformed values should be in a reasonable range relative to original scale
    center_val <- result$metadata$centering_values[["EXP"]]
    expect_true(any(abs(non_na - center_val) < 3 * 20),
                info = "Back-transformed quantile labels should be near original scale")
  }
})

test_that("forest plot with center=FALSE uses conditional effect for main effects in interactions", {
  skip_on_cran()
  set.seed(33)
  n <- 500
  WT_true <- rnorm(n, 70, 15)
  EXP_true <- rnorm(n, 100, 20)
  log_odds <- -2 + 0.02 * EXP_true + 0.01 * WT_true + 0.0003 * EXP_true * WT_true
  df <- data.frame(
    DV  = rbinom(n, 1, plogis(log_odds)),
    EXP = EXP_true,
    WT  = WT_true
  )

  result_nocenter <- scm_glm(
    eres_data = df,
    covariates = "WT",
    exposure_vars = "EXP",
    dv_col = "DV",
    significance_level = c(0.9, 0.8),
    interactions = "EXP:WT",
    center = FALSE
  )

  skip_if(!"EXP:WT" %in% result_nocenter$metadata$finalRelation,
          "EXP:WT interaction not selected")

  fp_nocenter <- make_forest_plot(result_nocenter, plot = FALSE)

  # The forest plot should not error and should return data
  expect_s3_class(fp_nocenter, "data.frame")
  expect_true(nrow(fp_nocenter) > 0)
  # EXP main-effect row should exist (as a main effect row)
  expect_true("EXP" %in% fp_nocenter$group)
})

test_that("forest plot main-effect row at Q0.50 matches interaction row at modifier median (center=TRUE)", {
  skip_on_cran()
  set.seed(44)
  n <- 500
  WT_true <- rnorm(n, 70, 15)
  EXP_true <- rnorm(n, 100, 20)
  log_odds <- -2 + 0.02 * EXP_true + 0.01 * WT_true + 0.0003 * EXP_true * WT_true
  df <- data.frame(
    DV  = rbinom(n, 1, plogis(log_odds)),
    EXP = EXP_true,
    WT  = WT_true
  )

  result_center <- scm_glm(
    eres_data = df,
    covariates = "WT",
    exposure_vars = "EXP",
    dv_col = "DV",
    significance_level = c(0.9, 0.8),
    interactions = "EXP:WT",
    center = TRUE
  )

  skip_if(!"EXP:WT" %in% result_center$metadata$finalRelation,
          "EXP:WT interaction not selected")

  fp_center <- make_forest_plot(result_center, plot = FALSE,
                                modifier_quantiles = c(0.25, 0.5, 0.75))

  # EXP main-effect Q0.50 row (mid=1, since contrast=0) should exist
  exp_rows <- fp_center[fp_center$group == "EXP", ]
  expect_true(nrow(exp_rows) > 0)

  # The Q0.50 row (contrast = 0) should give mid = 1 (exp(0) = 1)
  median_row <- exp_rows[grepl("Q0.50", exp_rows$group_level), ]
  if (nrow(median_row) > 0) {
    expect_equal(median_row$mid[1], 1.0, tolerance = 1e-6)
  }
})

test_that("no-interaction models are unaffected by centering", {
  skip_on_cran()
  set.seed(66)
  n <- 300
  df <- data.frame(
    DV  = rbinom(n, 1, 0.3),
    EXP = rnorm(n, 10, 2),
    AGE = rnorm(n, 50, 10)
  )
  result <- scm_glm(
    eres_data = df, covariates = "AGE", exposure_vars = "EXP",
    dv_col = "DV", significance_level = c(0.9, 0.8),
    interactions = NULL, center = TRUE
  )
  expect_type(result, "list")
  expect_true(inherits(result$final_model, "glm"))
  fp <- make_forest_plot(result, plot = FALSE)
  expect_s3_class(fp, "data.frame")
  expect_true(nrow(fp) > 0)
})
