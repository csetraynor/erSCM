## Tests for the apply_median_scaling() helper and normalize_inputs refactor

test_that("apply_median_scaling scales numeric columns to median=1", {
  df <- data.frame(X = c(2, 4, 6, 8, 10), Y = c(10, 20, 30, 40, 50), Z = c(1, 1, 1, 1, 1))
  result <- apply_median_scaling(df, c("X", "Y", "Z"))
  # median(X) = 6, so scaled X should have median 1
  expect_equal(stats::median(result$data$X), 1.0, tolerance = 1e-8)
  expect_equal(stats::median(result$data$Y), 1.0, tolerance = 1e-8)
  # scaling_factors
  expect_equal(unname(result$scaling_factors["X"]), 6.0)
  expect_equal(unname(result$scaling_factors["Y"]), 30.0)
})

test_that("apply_median_scaling skips non-numeric columns silently", {
  df <- data.frame(X = c(2, 4, 6), G = c("a", "b", "c"), stringsAsFactors = FALSE)
  result <- apply_median_scaling(df, c("X", "G"))
  # G unchanged
  expect_equal(result$data$G, c("a", "b", "c"))
  # Only X in naming_df
  expect_equal(result$naming_df$covname, "X")
})

test_that("apply_median_scaling returns NULL naming_df when no numeric cols in cols", {
  df <- data.frame(G = c("a", "b", "c"), stringsAsFactors = FALSE)
  result <- apply_median_scaling(df, "G")
  expect_null(result$naming_df)
  expect_length(result$scaling_factors, 0)
})

test_that("apply_median_scaling naming_df has transform_type = 'median_scaled' by default", {
  df <- data.frame(X = 1:10)
  result <- apply_median_scaling(df, "X")
  expect_equal(result$naming_df$transform_type, "median_scaled")
})

test_that("apply_median_scaling naming_df has transform_type = 'median_scaled_and_centered' when is_centered=TRUE", {
  df <- data.frame(X = 1:10)
  result <- apply_median_scaling(df, "X", is_centered = TRUE)
  expect_equal(result$naming_df$transform_type, "median_scaled_and_centered")
})

test_that("apply_median_scaling attaches units from units_df", {
  df <- data.frame(EXP = c(1, 2, 3, 4, 5), WT = c(50, 60, 70, 80, 90))
  udf <- data.frame(covname = c("EXP", "WT"), units = c("mg/L", "kg"), stringsAsFactors = FALSE)
  result <- apply_median_scaling(df, c("EXP", "WT"), units_df = udf)
  expect_equal(result$naming_df$units[result$naming_df$covname == "EXP"], "mg/L")
  expect_equal(result$naming_df$units[result$naming_df$covname == "WT"], "kg")
})

test_that("apply_median_scaling defaults to units='U' when no units_df", {
  df <- data.frame(X = 1:5)
  result <- apply_median_scaling(df, "X")
  expect_equal(result$naming_df$units, "U")
})

test_that("apply_median_scaling does not change column order", {
  df <- data.frame(A = 1:5, B = letters[1:5], C = 6:10, stringsAsFactors = FALSE)
  result <- apply_median_scaling(df, c("A", "C"))
  expect_equal(colnames(result$data), c("A", "B", "C"))
})

test_that("apply_median_scaling only scales cols argument, not all numeric columns", {
  df <- data.frame(X = c(2, 4, 6), Y = c(10, 20, 30), Z_unreferenced = c(100, 200, 300))
  result <- apply_median_scaling(df, c("X", "Y"))
  # Z_unreferenced should be untouched
  expect_equal(result$data$Z_unreferenced, c(100, 200, 300))
  # Only X and Y in naming_df
  expect_setequal(result$naming_df$covname, c("X", "Y"))
})

test_that("scm_glm with normalize_inputs produces naming_df with correct transform_type", {
  skip_on_cran()
  set.seed(42)
  n <- 200
  df <- data.frame(
    DV  = rbinom(n, 1, 0.3),
    EXP = rnorm(n, 10, 2),
    AGE = rnorm(n, 50, 10)
  )
  result <- scm_glm(
    eres_data = df, covariates = "AGE", exposure_vars = "EXP",
    dv_col = "DV", significance_level = c(0.9, 0.8),
    normalize_inputs = TRUE, center = FALSE
  )
  expect_false(is.null(result$naming_df))
  expect_true("transform_type" %in% colnames(result$naming_df))
  expect_true(all(result$naming_df$transform_type == "median_scaled"))
})

test_that("scm_glm with normalize_inputs=TRUE and center=TRUE produces warning and correct transform_type", {
  skip_on_cran()
  set.seed(42)
  n <- 200
  df <- data.frame(
    DV  = rbinom(n, 1, 0.3),
    EXP = rnorm(n, 10, 2),
    AGE = rnorm(n, 50, 10)
  )
  # normalize_inputs=TRUE triggers a deprecation warning; capture all warnings
  warns <- character(0)
  withCallingHandlers(
    result <- scm_glm(
      eres_data = df, covariates = "AGE", exposure_vars = "EXP",
      dv_col = "DV", significance_level = c(0.9, 0.8),
      normalize_inputs = TRUE, center = TRUE
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("deprecated", warns)),
              info = paste("Warnings:", paste(warns, collapse = "; ")))
  expect_true("transform_type" %in% colnames(result$naming_df))
  expect_true(all(result$naming_df$transform_type == "median_scaled_and_centered"))
})

test_that("scm_coxph with normalize_inputs produces naming_df with correct transform_type", {
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
    eres_data = df, covariates = "AGE", exposure_vars = "EXP",
    dv_col = "DV", time_col = "TIME",
    significance_level = c(0.9, 0.8),
    normalize_inputs = TRUE, center = FALSE
  )
  expect_false(is.null(result$naming_df))
  expect_true("transform_type" %in% colnames(result$naming_df))
  expect_true(all(result$naming_df$transform_type == "median_scaled"))
})

test_that("scm_glm stores scaling_method in metadata", {
  skip_on_cran()
  set.seed(1)
  n <- 200
  df <- data.frame(DV = rbinom(n, 1, 0.3), EXP = rnorm(n, 5, 1))

  r1 <- scm_glm(eres_data = df, covariates = character(0), exposure_vars = "EXP",
                dv_col = "DV", significance_level = c(0.9, 0.8),
                normalize_inputs = FALSE, center = FALSE)
  expect_equal(r1$metadata$scaling_method, "none")

  r2 <- scm_glm(eres_data = df, covariates = character(0), exposure_vars = "EXP",
                dv_col = "DV", significance_level = c(0.9, 0.8),
                normalize_inputs = FALSE, center = TRUE)
  expect_equal(r2$metadata$scaling_method, "centered")

  r3 <- scm_glm(eres_data = df, covariates = character(0), exposure_vars = "EXP",
                dv_col = "DV", significance_level = c(0.9, 0.8),
                normalize_inputs = TRUE, center = FALSE)
  expect_equal(r3$metadata$scaling_method, "median_scaled")
})

test_that("forest plot handles naming_df with transform_type dispatch correctly", {
  skip_on_cran()
  set.seed(22)
  n <- 300
  df <- data.frame(
    DV  = rbinom(n, 1, 0.3),
    EXP = rnorm(n, 100, 20),
    AGE = rnorm(n, 50, 10)
  )
  result <- scm_glm(
    eres_data = df, covariates = "AGE", exposure_vars = "EXP",
    dv_col = "DV", significance_level = c(0.9, 0.8),
    normalize_inputs = TRUE, center = FALSE
  )
  skip_if(length(result$metadata$finalRelation) == 0, "No covariates selected")
  fp <- make_forest_plot(result, plot = FALSE, namings = result$naming_df)
  expect_s3_class(fp, "data.frame")
  expect_true(nrow(fp) > 0)
})

test_that("forest plot backward compat: old naming_df without transform_type column (median_scaled branch)", {
  skip_on_cran()
  set.seed(22)
  n <- 300
  df <- data.frame(
    DV  = rbinom(n, 1, 0.3),
    EXP = rnorm(n, 100, 20)
  )
  result <- scm_glm(
    eres_data = df, covariates = character(0), exposure_vars = "EXP",
    dv_col = "DV", significance_level = c(0.9, 0.8),
    normalize_inputs = TRUE, center = FALSE
  )
  skip_if(length(result$metadata$finalRelation) == 0, "No covariates selected")

  # Simulate legacy naming_df without transform_type
  legacy_naming <- result$naming_df[, c("covname", "median_point", "units")]

  fp <- make_forest_plot(result, plot = FALSE, namings = legacy_naming)
  expect_s3_class(fp, "data.frame")
  expect_true(nrow(fp) > 0)
})

test_that("forest plot backward compat: old naming_df without transform_type (log-transformed branch)", {
  skip_on_cran()
  set.seed(22)
  n <- 300
  df <- data.frame(
    DV  = rbinom(n, 1, 0.3),
    EXP = rnorm(n, 100, 20)
  )
  result <- scm_glm(
    eres_data = df, covariates = character(0), exposure_vars = "EXP",
    dv_col = "DV", significance_level = c(0.9, 0.8),
    normalize_inputs = FALSE, center = FALSE
  )
  skip_if(length(result$metadata$finalRelation) == 0, "No covariates selected")

  # Simulate legacy log-transformed naming_df (median_point as character)
  legacy_naming <- data.frame(covname = "EXP", median_point = "log", units = "U",
                               stringsAsFactors = FALSE)

  fp <- make_forest_plot(result, plot = FALSE, namings = legacy_naming)
  expect_s3_class(fp, "data.frame")
  expect_true(nrow(fp) > 0)
})

test_that("model fit equivalence: normalize_inputs only scales input, final predictions proportional", {
  skip_on_cran()
  set.seed(99)
  n <- 300
  df <- data.frame(
    DV  = rbinom(n, 1, 0.3),
    EXP = rnorm(n, 10, 2)
  )
  r_raw <- scm_glm(eres_data = df, covariates = character(0), exposure_vars = "EXP",
                   dv_col = "DV", significance_level = c(0.9, 0.8),
                   normalize_inputs = FALSE, center = FALSE)
  r_norm <- scm_glm(eres_data = df, covariates = character(0), exposure_vars = "EXP",
                    dv_col = "DV", significance_level = c(0.9, 0.8),
                    normalize_inputs = TRUE, center = FALSE)

  skip_if(length(r_raw$metadata$finalRelation) == 0 || length(r_norm$metadata$finalRelation) == 0,
          "No covariates selected")
  skip_if(!setequal(r_raw$metadata$finalRelation, r_norm$metadata$finalRelation),
          "Different models selected")

  # Fitted probabilities must be identical (scaling doesn't change fitted values)
  expect_equal(fitted(r_raw$final_model), fitted(r_norm$final_model), tolerance = 1e-6)
})
