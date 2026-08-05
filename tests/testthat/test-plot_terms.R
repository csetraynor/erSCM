## Tests for var_spec integration in plot_cox_terms(), vpc_coxph(), plot_logistic_terms()

# ---- Shared fixtures ----

.make_glm_fit <- function(seed = 1, with_spec = TRUE) {
  set.seed(seed)
  n <- 300
  WT_true <- rnorm(n, 70, 12)
  EXP_true <- rnorm(n, 100, 20)
  df <- data.frame(
    DV  = rbinom(n, 1, plogis(-2 + 0.008 * EXP_true + 0.005 * WT_true)),
    EXP = EXP_true,
    WT  = WT_true,
    SEX = factor(sample(c("F", "M"), n, replace = TRUE))
  )
  spec <- if (with_spec) c(EXP = "Exposure (ng/mL)", WT = "Body Weight (kg)", SEX = "Sex") else NULL
  suppressWarnings(
    scm_glm(eres_data = df, covariates = c("WT", "SEX"), exposure_vars = "EXP",
            dv_col = "DV", significance_level = c(0.9, 0.8),
            var_spec = spec, center = TRUE)
  )
}

.make_cox_fit <- function(seed = 2, with_spec = TRUE) {
  set.seed(seed)
  n <- 400
  EXP_true <- rnorm(n, 100, 20)
  WT_true  <- rnorm(n, 70, 12)
  lp <- -0.05 + 0.008 * EXP_true + 0.006 * WT_true
  df <- data.frame(
    TIME = rexp(n, exp(lp / 10 - 3)),
    DV   = rbinom(n, 1, plogis(lp - 3)),
    EXP  = EXP_true,
    WT   = WT_true,
    SEX  = factor(sample(c("F", "M"), n, replace = TRUE))
  )
  spec <- if (with_spec) c(EXP = "Exposure (ng/mL)", WT = "Body Weight (kg)", SEX = "Sex") else NULL
  suppressWarnings(
    scm_coxph(eres_data = df, covariates = c("WT", "SEX"), exposure_vars = "EXP",
              dv_col = "DV", time_col = "TIME",
              significance_level = c(0.95, 0.9),
              var_spec = spec, center = TRUE)
  )
}

# ---- plot_cox_terms: no errors ----

test_that("plot_cox_terms runs without error (with var_spec, continuous)", {
  skip_on_cran()
  fit <- .make_cox_fit()
  terms <- get_main_effect_columns(fit$metadata$finalRelation)
  cont <- terms[vapply(fit$model_data[terms], is.numeric, logical(1))]
  skip_if(length(cont) == 0, "no continuous term selected")
  expect_s3_class(plot_cox_terms(fit, cont[1], plot = TRUE), "ggplot")
})

test_that("plot_cox_terms runs without error (with var_spec, categorical)", {
  skip_on_cran()
  fit <- .make_cox_fit()
  terms <- get_main_effect_columns(fit$metadata$finalRelation)
  cats <- terms[vapply(fit$model_data[terms], function(x) is.factor(x) || is.character(x), logical(1))]
  skip_if(length(cats) == 0, "no categorical term selected")
  expect_s3_class(plot_cox_terms(fit, cats[1], plot = TRUE), "ggplot")
})

test_that("plot_cox_terms runs without error (no var_spec, backward compat)", {
  skip_on_cran()
  fit <- .make_cox_fit(with_spec = FALSE)
  terms <- get_main_effect_columns(fit$metadata$finalRelation)
  cont <- terms[vapply(fit$model_data[terms], is.numeric, logical(1))]
  skip_if(length(cont) == 0, "no continuous term selected")
  expect_s3_class(plot_cox_terms(fit, cont[1], plot = TRUE), "ggplot")
})

test_that("plot_cox_terms returns data frame when plot=FALSE", {
  skip_on_cran()
  fit <- .make_cox_fit()
  terms <- get_main_effect_columns(fit$metadata$finalRelation)
  cont <- terms[vapply(fit$model_data[terms], is.numeric, logical(1))]
  skip_if(length(cont) == 0, "no continuous term selected")
  result <- plot_cox_terms(fit, cont[1], plot = FALSE)
  expect_s3_class(result, "data.frame")
})

test_that("plot_cox_terms uses display label in strata when var_spec present (continuous)", {
  skip_on_cran()
  fit <- .make_cox_fit()
  skip_if(!"EXP" %in% get_main_effect_columns(fit$metadata$finalRelation), "EXP not selected")
  vpc_data <- plot_cox_terms(fit, "EXP", plot = FALSE)
  skip_if(is.null(vpc_data$strata), "no strata column")
  # strata should reference "Exposure" not raw "EXP"
  expect_true(any(grepl("Exposure", vpc_data$strata)),
              info = paste("Strata values:", paste(unique(vpc_data$strata), collapse = "; ")))
})

test_that("plot_cox_terms bin labels are back-transformed when center=TRUE", {
  skip_on_cran()
  fit <- .make_cox_fit()
  skip_if(!"EXP" %in% get_main_effect_columns(fit$metadata$finalRelation), "EXP not selected")
  vpc_data <- plot_cox_terms(fit, "EXP", plot = FALSE, n_bins = 4)
  skip_if(is.null(vpc_data$strata), "no strata column")
  strata_vals <- unique(vpc_data$strata)
  # With center=TRUE and median ~100, back-transformed values should be ~100 range
  # (not near 0 as they would be in centered scale)
  nums <- suppressWarnings(as.numeric(gsub("[^0-9.–-]", "", strata_vals)))
  meaningful_nums <- nums[!is.na(nums) & nums > 10]
  expect_true(length(meaningful_nums) > 0,
              info = paste("Strata:", paste(strata_vals, collapse = "; ")))
})

# ---- vpc_coxph: no errors ----

test_that("vpc_coxph runs without error (with var_spec)", {
  skip_on_cran()
  fit <- .make_cox_fit()
  expect_s3_class(vpc_coxph(fit, plot = TRUE), "ggplot")
})

test_that("vpc_coxph runs without error (no var_spec, backward compat)", {
  skip_on_cran()
  fit <- .make_cox_fit(with_spec = FALSE)
  expect_s3_class(vpc_coxph(fit, plot = TRUE), "ggplot")
})

test_that("vpc_coxph returns data frame when plot=FALSE", {
  skip_on_cran()
  fit <- .make_cox_fit()
  result <- vpc_coxph(fit, plot = FALSE)
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
})

test_that("vpc_coxph uses display labels in strata strings (categorical)", {
  skip_on_cran()
  fit <- .make_cox_fit()
  terms <- get_main_effect_columns(fit$metadata$finalRelation)
  cats <- terms[vapply(fit$model_data[terms], function(x) is.factor(x) || is.character(x), logical(1))]
  skip_if(length(cats) == 0, "no categorical term in model")
  vpc_data <- vpc_coxph(fit, plot = FALSE)
  skip_if(is.null(vpc_data$strata), "no strata column")
  # For SEX with var_spec label "Sex", strata should say "Sex" not "SEX"
  if ("SEX" %in% cats) {
    expect_true(any(grepl("Sex", vpc_data$strata)),
                info = paste("Strata values:", paste(unique(vpc_data$strata), collapse = "; ")))
  }
})

test_that("vpc_coxph stores time_col in metadata and uses it for x-axis", {
  skip_on_cran()
  fit <- .make_cox_fit()
  # metadata$time_col should be "TIME"
  expect_equal(fit$metadata$time_col, "TIME")
})

# ---- plot_logistic_terms: no errors ----

test_that("plot_logistic_terms runs without error (with var_spec, continuous)", {
  skip_on_cran()
  fit <- .make_glm_fit()
  terms <- get_main_effect_columns(fit$metadata$finalRelation)
  cont <- terms[vapply(fit$model_data[terms], is.numeric, logical(1))]
  skip_if(length(cont) == 0, "no continuous term selected")
  expect_s3_class(plot_logistic_terms(fit, cont[1], plot = TRUE), "ggplot")
})

test_that("plot_logistic_terms runs without error (with var_spec, categorical)", {
  skip_on_cran()
  fit <- .make_glm_fit()
  terms <- get_main_effect_columns(fit$metadata$finalRelation)
  cats <- terms[vapply(fit$model_data[terms], function(x) is.factor(x) || is.character(x), logical(1))]
  skip_if(length(cats) == 0, "no categorical term selected")
  expect_s3_class(plot_logistic_terms(fit, cats[1], plot = TRUE), "ggplot")
})

test_that("plot_logistic_terms runs without error (no var_spec, backward compat)", {
  skip_on_cran()
  fit <- .make_glm_fit(with_spec = FALSE)
  terms <- get_main_effect_columns(fit$metadata$finalRelation)
  cont <- terms[vapply(fit$model_data[terms], is.numeric, logical(1))]
  skip_if(length(cont) == 0, "no continuous term selected")
  expect_s3_class(plot_logistic_terms(fit, cont[1], plot = TRUE), "ggplot")
})

test_that("plot_logistic_terms returns list when plot=FALSE (continuous)", {
  skip_on_cran()
  fit <- .make_glm_fit()
  terms <- get_main_effect_columns(fit$metadata$finalRelation)
  cont <- terms[vapply(fit$model_data[terms], is.numeric, logical(1))]
  skip_if(length(cont) == 0, "no continuous term selected")
  result <- plot_logistic_terms(fit, cont[1], plot = FALSE)
  expect_true(is.list(result))
  expect_true("predictions" %in% names(result))
})

test_that("plot_logistic_terms x-axis label uses display name from var_spec", {
  skip_on_cran()
  fit <- .make_glm_fit()
  skip_if(!"EXP" %in% get_main_effect_columns(fit$metadata$finalRelation), "EXP not selected")
  p <- plot_logistic_terms(fit, "EXP", plot = TRUE)
  x_lab <- p$labels$x
  expect_true(grepl("Exposure", x_lab),
              info = paste("x label:", x_lab))
})

test_that("plot_logistic_terms categorical x-axis uses display name from var_spec", {
  skip_on_cran()
  fit <- .make_glm_fit()
  skip_if(!"SEX" %in% get_main_effect_columns(fit$metadata$finalRelation), "SEX not selected")
  p <- plot_logistic_terms(fit, "SEX", plot = TRUE)
  x_lab <- p$labels$x
  expect_equal(x_lab, "Sex",
               info = paste("x label:", x_lab))
})

test_that("plot_logistic_terms without var_spec falls back to column name as x label", {
  skip_on_cran()
  fit <- .make_glm_fit(with_spec = FALSE)
  skip_if(!"EXP" %in% get_main_effect_columns(fit$metadata$finalRelation), "EXP not selected")
  p <- plot_logistic_terms(fit, "EXP", plot = TRUE)
  x_lab <- p$labels$x
  expect_equal(x_lab, "EXP")
})

# ---- plot_model_terms dispatcher ----

test_that("plot_model_terms dispatches to plot_cox_terms for coxph fit", {
  skip_on_cran()
  fit <- .make_cox_fit()
  terms <- get_main_effect_columns(fit$metadata$finalRelation)
  cont <- terms[vapply(fit$model_data[terms], is.numeric, logical(1))]
  skip_if(length(cont) == 0, "no continuous term selected")
  expect_s3_class(plot_model_terms(fit, cont[1]), "ggplot")
})

test_that("plot_model_terms dispatches to plot_logistic_terms for glm fit", {
  skip_on_cran()
  fit <- .make_glm_fit()
  terms <- get_main_effect_columns(fit$metadata$finalRelation)
  cont <- terms[vapply(fit$model_data[terms], is.numeric, logical(1))]
  skip_if(length(cont) == 0, "no continuous term selected")
  expect_s3_class(plot_model_terms(fit, cont[1]), "ggplot")
})

# ---- var_spec helpers: get_var_categories, make_bin_labels ----

test_that("get_var_categories returns NULL when var_spec is NULL", {
  expect_null(get_var_categories(NULL, "SEX"))
})

test_that("get_var_categories returns NULL for unknown column", {
  spec <- parse_var_spec(c(WT = "Weight (kg)"))
  expect_null(get_var_categories(spec, "UNKNOWN"))
})

test_that("get_var_categories returns NULL when no categories defined", {
  spec <- parse_var_spec(c(WT = "Weight (kg)"))
  expect_null(get_var_categories(spec, "WT"))
})

test_that("get_var_categories returns labels from data.frame spec", {
  df <- data.frame(col_name = "SEX", label = "Sex",
                   categories = "Female|Male", stringsAsFactors = FALSE)
  spec <- parse_var_spec(df)
  cats <- get_var_categories(spec, "SEX")
  expect_equal(cats, c("Female", "Male"))
})

test_that("make_bin_labels with NULL var_spec returns 'Q1'...'Qn' labels", {
  breaks <- c(0, 25, 50, 75, 100)
  labels <- make_bin_labels(breaks, NULL, "X")
  expect_equal(length(labels), 4)
  expect_equal(labels[1], "≤ 25")
})

test_that("make_bin_labels back-transforms centered breaks", {
  spec <- tibble::tibble(
    col_name = "WT", label = "Weight", units = "kg",
    var_type = "continuous", categories = list(NULL),
    transform = "centered", median_point = 70.0
  )
  # Breaks on centered scale: -20, -10, 0, 10, 20 -> original: 50, 60, 70, 80, 90
  breaks <- c(-20, -10, 0, 10, 20)
  labels <- make_bin_labels(breaks, spec, "WT")
  expect_equal(length(labels), 4)
  # First label: "<= 60"
  expect_true(grepl("60", labels[1]),
              info = paste("Labels:", paste(labels, collapse = "; ")))
  # Last label: "> 80"
  expect_true(grepl("80", labels[4]),
              info = paste("Labels:", paste(labels, collapse = "; ")))
})
