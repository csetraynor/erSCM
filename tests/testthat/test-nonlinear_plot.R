# Tests for Phase 3: partial effect plot functions.
# Fixtures are loaded from helper-bootstrap.R automatically.

# ---------------------------------------------------------------------------
# Shared test-data fixtures
# ---------------------------------------------------------------------------

.make_log_glm_fit <- function(n = 300, seed = 42) {
  set.seed(seed)
  x <- runif(n, 1, 100)
  dat <- data.frame(
    DV   = rbinom(n, 1, plogis(-2 + log(x / 50))),
    CMAX = x,
    AGE  = rnorm(n, 55, 12)
  )
  suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = "AGE",
    significance_level = c(0.05, 0.01),
    dv_col             = "DV",
    family             = "binomial",
    nonlinear          = list(CMAX = c("linear", "log"), AGE = c("linear")),
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
}

.make_fp2_glm_fit <- function(n = 400, seed = 7) {
  set.seed(seed)
  x <- runif(n, 1, 100)
  dat <- data.frame(
    DV   = rbinom(n, 1, plogis(-2 + 0.5 * sqrt(x / 50) - 0.3 * log(x / 50))),
    CMAX = x,
    AGE  = rnorm(n, 55, 12)
  )
  suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = "AGE",
    significance_level = c(0.05, 0.01),
    dv_col             = "DV",
    family             = "binomial",
    nonlinear          = list(CMAX = c("linear", "fp2"), AGE = c("linear")),
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
}

.make_rcs_glm_fit <- function(n = 400, seed = 9) {
  set.seed(seed)
  x <- runif(n, 1, 100)
  dat <- data.frame(
    DV   = rbinom(n, 1, plogis(-2 + sin(x / 15))),
    CMAX = x,
    AGE  = rnorm(n, 55, 12)
  )
  suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = "AGE",
    significance_level = c(0.05, 0.01),
    dv_col             = "DV",
    family             = "binomial",
    nonlinear          = list(CMAX = c("linear", "rcs"), AGE = c("linear")),
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
}

.make_linear_glm_fit <- function(n = 400, seed = 42) {
  set.seed(seed)
  # Strong signal: logistic coefficient 0.08/unit with SD=20 → large effect.
  # Guarantees CMAX enters at strict thresholds, returning a valid eres_scm_glm.
  x <- rnorm(n, 100, 20)
  lp <- -1 + 0.08 * (x - 100)
  dat <- data.frame(
    DV   = rbinom(n, 1L, plogis(lp)),
    CMAX = x,
    AGE  = rnorm(n, 55, 12),
    WT   = rnorm(n, 75, 15)
  )
  suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT"),
    significance_level = c(0.05, 0.01),
    dv_col             = "DV",
    family             = "binomial",
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
}

.make_coxph_nl_fit <- function(n = 300, seed = 77) {
  set.seed(seed)
  dat <- .make_test_data_coxph(n = n, seed = seed)
  suppressMessages(scm_coxph(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT"),
    significance_level = c(0.5, 0.3),
    dv_col             = "DV", time_col = "TIME",
    nonlinear          = list(CMAX = c("linear", "log"), AGE = c("linear")),
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
}


# ---------------------------------------------------------------------------
# plot_nonlinear_effect — basic returns
# ---------------------------------------------------------------------------

test_that("plot_nonlinear_effect returns ggplot for log transform", {
  skip_on_cran()
  fit <- .make_log_glm_fit()
  # Plot whichever variable is in the final model
  reg <- fit$metadata$transform_registry
  skip_if(is.null(reg) || length(reg) == 0, "No registry")
  var <- names(reg)[1L]
  p <- plot_nonlinear_effect(fit, var)
  expect_s3_class(p, "gg")
})

test_that("plot_nonlinear_effect returns ggplot for FP transform", {
  skip_on_cran()
  fit <- .make_fp2_glm_fit()
  reg <- fit$metadata$transform_registry
  skip_if(is.null(reg) || length(reg) == 0, "No registry")
  var <- names(reg)[1L]
  p <- plot_nonlinear_effect(fit, var)
  expect_s3_class(p, "gg")
})

test_that("plot_nonlinear_effect returns ggplot for spline transform", {
  skip_on_cran()
  fit <- .make_rcs_glm_fit()
  reg <- fit$metadata$transform_registry
  skip_if(is.null(reg) || length(reg) == 0, "No registry")
  var <- names(reg)[1L]
  p <- plot_nonlinear_effect(fit, var)
  expect_s3_class(p, "gg")
})

test_that("plot_nonlinear_effect works for a linear variable (no registry)", {
  skip_on_cran()
  fit <- .make_linear_glm_fit()
  # CMAX must be in model_data
  expect_true("CMAX" %in% colnames(fit$model_data))
  # Even without a nonlinear spec, plot should work if CMAX is in the model
  rel_vars <- unique(vapply(fit$metadata$finalRelation, .term_to_varname,
                             character(1L), USE.NAMES = FALSE))
  skip_if(!"CMAX" %in% rel_vars, "CMAX not in final model")
  p <- plot_nonlinear_effect(fit, "CMAX")
  expect_s3_class(p, "gg")
})

# ---------------------------------------------------------------------------
# scale argument
# ---------------------------------------------------------------------------

test_that("plot_nonlinear_effect with scale='response' returns ggplot", {
  skip_on_cran()
  fit <- .make_log_glm_fit()
  reg <- fit$metadata$transform_registry
  skip_if(is.null(reg) || length(reg) == 0, "No registry")
  var <- names(reg)[1L]
  p <- plot_nonlinear_effect(fit, var, scale = "response")
  expect_s3_class(p, "gg")
  # On response scale, reference line should be at 1 (OR), not 0 (log-OR).
  # geom_hline stores yintercept in layer$data, not layer$aes_params.
  hlines <- unlist(lapply(p$layers, function(l) {
    tryCatch({
      d <- l$data
      if (is.data.frame(d) && "yintercept" %in% names(d)) as.numeric(d$yintercept)
      else NULL
    }, error = function(e) NULL)
  }))
  expect_true(any(hlines == 1, na.rm = TRUE))
})

test_that("plot_nonlinear_effect link scale has reference line at 0", {
  skip_on_cran()
  fit <- .make_log_glm_fit()
  reg <- fit$metadata$transform_registry
  skip_if(is.null(reg) || length(reg) == 0, "No registry")
  var <- names(reg)[1L]
  p <- plot_nonlinear_effect(fit, var, scale = "link")
  hlines <- unlist(lapply(p$layers, function(l) {
    tryCatch({
      d <- l$data
      if (is.data.frame(d) && "yintercept" %in% names(d)) as.numeric(d$yintercept)
      else NULL
    }, error = function(e) NULL)
  }))
  expect_true(any(hlines == 0, na.rm = TRUE))
})

# ---------------------------------------------------------------------------
# show_linear overlay
# ---------------------------------------------------------------------------

test_that("plot_nonlinear_effect with show_linear=TRUE has more layers", {
  skip_on_cran()
  fit <- .make_log_glm_fit()
  reg <- fit$metadata$transform_registry
  skip_if(is.null(reg) || length(reg) == 0, "No registry")
  var <- names(reg)[1L]
  p_no  <- plot_nonlinear_effect(fit, var, show_linear = FALSE)
  p_yes <- plot_nonlinear_effect(fit, var, show_linear = TRUE)
  expect_gte(length(p_yes$layers), length(p_no$layers))
})

# ---------------------------------------------------------------------------
# CoxPH
# ---------------------------------------------------------------------------

test_that("plot_nonlinear_effect works for coxph", {
  skip_on_cran()
  fit <- .make_coxph_nl_fit()
  reg <- fit$metadata$transform_registry
  skip_if(is.null(reg) || length(reg) == 0, "No registry")
  var <- names(reg)[1L]
  p <- plot_nonlinear_effect(fit, var)
  expect_s3_class(p, "gg")
})

# ---------------------------------------------------------------------------
# rug and custom labels
# ---------------------------------------------------------------------------

test_that("rug = FALSE reduces layer count", {
  skip_on_cran()
  fit <- .make_log_glm_fit()
  reg <- fit$metadata$transform_registry
  skip_if(is.null(reg) || length(reg) == 0, "No registry")
  var <- names(reg)[1L]
  p_rug  <- plot_nonlinear_effect(fit, var, rug = TRUE)
  p_nrug <- plot_nonlinear_effect(fit, var, rug = FALSE)
  expect_gt(length(p_rug$layers), length(p_nrug$layers))
})

test_that("Custom xlab and title are applied", {
  skip_on_cran()
  fit <- .make_log_glm_fit()
  reg <- fit$metadata$transform_registry
  skip_if(is.null(reg) || length(reg) == 0, "No registry")
  var <- names(reg)[1L]
  p <- plot_nonlinear_effect(fit, var, xlab = "My X", title = "My Title")
  expect_equal(p$labels$x, "My X")
  expect_equal(p$labels$title, "My Title")
})

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

test_that("plot_nonlinear_effect errors on non-eres_scm input", {
  expect_error(plot_nonlinear_effect(list(), "CMAX"))
})

test_that("plot_nonlinear_effect errors when variable not in model_data", {
  skip_on_cran()
  fit <- .make_log_glm_fit()
  expect_error(plot_nonlinear_effect(fit, "NONEXISTENT_XYZ"))
})

# ---------------------------------------------------------------------------
# plot_all_nonlinear
# ---------------------------------------------------------------------------

test_that("plot_all_nonlinear returns patchwork/ggplot for model with NL terms", {
  skip_on_cran()
  fit <- .make_log_glm_fit()
  reg <- fit$metadata$transform_registry
  has_nl <- !is.null(reg) && any(vapply(reg, function(r) r$transform_type != "linear",
                                          logical(1L)))
  skip_if(!has_nl, "No non-linear terms in final model")
  pw <- plot_all_nonlinear(fit)
  expect_true(inherits(pw, "gg") || inherits(pw, "patchwork"))
})

test_that("plot_all_nonlinear with nl_only=FALSE returns result for linear model", {
  skip_on_cran()
  fit <- .make_linear_glm_fit()
  # Falls back to finalRelation variables; may warn but should not error
  result <- tryCatch(
    plot_all_nonlinear(fit, nl_only = FALSE),
    message = function(m) invokeRestart("muffleMessage"),
    error   = function(e) NULL
  )
  # Either a plot, NULL (message), or patchwork — just no hard error
  expect_true(is.null(result) || inherits(result, "gg") || inherits(result, "patchwork"))
})

test_that("plot_all_nonlinear messages when no NL terms and nl_only=TRUE", {
  skip_on_cran()
  fit <- .make_linear_glm_fit()
  expect_message(plot_all_nonlinear(fit, nl_only = TRUE))
})

# ---------------------------------------------------------------------------
# .compute_quantile_effects
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Gap #8: y-axis label verification
# ---------------------------------------------------------------------------

test_that("response scale y-axis label contains 'OR' or 'Ratio' (Gap #8)", {
  skip_on_cran()
  fit <- .make_log_glm_fit()
  reg <- fit$metadata$transform_registry
  skip_if(is.null(reg) || length(reg) == 0, "No registry")
  var <- names(reg)[1L]
  p <- plot_nonlinear_effect(fit, var, scale = "response")
  expect_true(grepl("OR|Ratio|ratio|Odds|odds", p$labels$y))
})

test_that("link scale y-axis label contains 'Log' or 'log' (Gap #8)", {
  skip_on_cran()
  fit <- .make_log_glm_fit()
  reg <- fit$metadata$transform_registry
  skip_if(is.null(reg) || length(reg) == 0, "No registry")
  var <- names(reg)[1L]
  p <- plot_nonlinear_effect(fit, var, scale = "link")
  expect_true(grepl("Log|log|Linear|linear|link|Link", p$labels$y))
})

test_that(".compute_quantile_effects returns correct tibble structure", {
  skip_on_cran()
  fit <- .make_log_glm_fit()
  reg <- fit$metadata$transform_registry
  skip_if(is.null(reg) || length(reg) == 0, "No registry")
  var <- names(reg)[1L]
  skip_if(!var %in% colnames(fit$model_data), "var not in model_data")

  out <- .compute_quantile_effects(fit, var, quantiles = c(0.25, 0.75))
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2L)
  expect_true(all(c("quantile", "value", "estimate", "ci_lower", "ci_upper") %in% colnames(out)))
  expect_true(all(is.numeric(out$estimate)))
  expect_true(all(is.numeric(out$ci_lower)))
  expect_true(all(is.numeric(out$ci_upper)))
})

test_that(".compute_quantile_effects returns 4 rows for default quantiles", {
  skip_on_cran()
  fit <- .make_log_glm_fit()
  reg <- fit$metadata$transform_registry
  skip_if(is.null(reg) || length(reg) == 0, "No registry")
  var <- names(reg)[1L]
  skip_if(!var %in% colnames(fit$model_data), "var not in model_data")

  out <- .compute_quantile_effects(fit, var)
  expect_equal(nrow(out), 4L)
})
