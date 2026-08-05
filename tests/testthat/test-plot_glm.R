## Tests for plot_glm_terms(), plot_glm_continuous(), plot_glm_categorical(),
## and the internal helpers bin_continuous_covariate() / predict_glm_curve() /
## predict_glm_categorical().

# ---- Shared fixtures --------------------------------------------------------

make_cont_data <- function(n = 200L) {
  set.seed(99)
  data.frame(x = rnorm(n), y = rbinom(n, 1L, 0.35))
}

make_cat_data <- function(n = 300L) {
  set.seed(99)
  data.frame(
    g = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
    y = rbinom(n, 1L, 0.30)
  )
}

make_mixed_data <- function(n = 300L) {
  set.seed(99)
  data.frame(
    cont = rnorm(n),
    cat  = factor(sample(c("X", "Y"), n, replace = TRUE)),
    y    = rbinom(n, 1L, 0.40)
  )
}

# ---- bin_continuous_covariate -----------------------------------------------

test_that("bin_continuous_covariate returns the correct number of bins", {
  d <- make_cont_data()
  res <- bin_continuous_covariate(d, "x", "y", n_bins = 4L)
  expect_equal(nrow(res), 4L)
})

test_that("bin_continuous_covariate has the required columns", {
  d <- make_cont_data()
  res <- bin_continuous_covariate(d, "x", "y", n_bins = 4L)
  expect_true(all(c("bin_midpoint", "n_events", "n_total",
                    "observed_prob", "ci_lower", "ci_upper") %in% names(res)))
})

test_that("bin_continuous_covariate: CIs bracket the observed proportion", {
  d <- make_cont_data()
  res <- bin_continuous_covariate(d, "x", "y", n_bins = 4L)
  expect_true(all(res$ci_lower <= res$observed_prob))
  expect_true(all(res$ci_upper >= res$observed_prob))
})

test_that("bin_continuous_covariate: all-events bins have observed_prob = 1", {
  d <- data.frame(x = 1:20, y = rep(1L, 20L))
  res <- bin_continuous_covariate(d, "x", "y", n_bins = 4L)
  expect_true(all(res$observed_prob == 1))
})

test_that("bin_continuous_covariate: no-events bins have observed_prob = 0", {
  d <- data.frame(x = 1:20, y = rep(0L, 20L))
  res <- bin_continuous_covariate(d, "x", "y", n_bins = 4L)
  expect_true(all(res$observed_prob == 0))
})

test_that("bin_continuous_covariate handles group_var", {
  set.seed(1)
  d <- data.frame(
    x   = rnorm(200),
    y   = rbinom(200, 1L, 0.30),
    grp = rep(c("A", "B"), each = 100L)
  )
  res <- bin_continuous_covariate(d, "x", "y", n_bins = 4L, group_var = "grp")
  expect_true("grp" %in% names(res))
  expect_equal(length(unique(res$grp)), 2L)
  expect_equal(nrow(res), 8L)
})

test_that("bin_continuous_covariate stops on missing term column", {
  d <- make_cont_data()
  expect_error(
    bin_continuous_covariate(d, "NOPE", "y"),
    "not found"
  )
})

test_that("bin_continuous_covariate stops on missing response column", {
  d <- make_cont_data()
  expect_error(
    bin_continuous_covariate(d, "x", "NOPE"),
    "not found"
  )
})

# ---- predict_glm_curve ------------------------------------------------------

test_that("predict_glm_curve returns the requested number of rows", {
  d <- make_cont_data()
  m <- stats::glm(y ~ x, family = stats::binomial, data = d)
  res <- predict_glm_curve(m, "x", d, n_points = 50L)
  expect_equal(nrow(res), 50L)
})

test_that("predict_glm_curve: predicted probabilities are in [0, 1]", {
  d <- make_cont_data()
  m <- stats::glm(y ~ x, family = stats::binomial, data = d)
  res <- predict_glm_curve(m, "x", d)
  expect_true(all(res$predicted_prob >= 0 & res$predicted_prob <= 1))
})

test_that("predict_glm_curve: CIs bracket the prediction", {
  d <- make_cont_data()
  m <- stats::glm(y ~ x, family = stats::binomial, data = d)
  res <- predict_glm_curve(m, "x", d)
  expect_true(all(res$ci_lower <= res$predicted_prob))
  expect_true(all(res$ci_upper >= res$predicted_prob))
})

test_that("predict_glm_curve has required columns", {
  d <- make_cont_data()
  m <- stats::glm(y ~ x, family = stats::binomial, data = d)
  res <- predict_glm_curve(m, "x", d)
  expect_true(all(c("term_value", "predicted_prob", "ci_lower", "ci_upper") %in% names(res)))
})

# ---- predict_glm_categorical ------------------------------------------------

test_that("predict_glm_categorical probability metric returns one row per level", {
  d <- make_cat_data()
  m <- stats::glm(y ~ g, family = stats::binomial, data = d)
  res <- predict_glm_categorical(m, "g", d, metric = "probability")
  expect_equal(nrow(res), nlevels(d$g))
})

test_that("predict_glm_categorical: estimates in [0,1] for probability metric", {
  d <- make_cat_data()
  m <- stats::glm(y ~ g, family = stats::binomial, data = d)
  res <- predict_glm_categorical(m, "g", d, metric = "probability")
  expect_true(all(res$estimate >= 0 & res$estimate <= 1))
})

test_that("predict_glm_categorical odds_ratio metric: reference level has OR = 1", {
  d <- make_cat_data()
  m <- stats::glm(y ~ g, family = stats::binomial, data = d)
  res <- predict_glm_categorical(m, "g", d, metric = "odds_ratio")
  ref_level <- levels(d$g)[1L]
  ref_row <- res[res$level == ref_level, ]
  if (nrow(ref_row) > 0L) {
    expect_equal(ref_row$estimate, 1, tolerance = 1e-8)
  }
})

# ---- plot_glm_continuous ----------------------------------------------------

test_that("plot_glm_continuous returns a ggplot", {
  d <- make_cont_data()
  m <- stats::glm(y ~ x, family = stats::binomial, data = d)
  p <- plot_glm_continuous(m, d, "x", "y")
  expect_s3_class(p, "gg")
})

test_that("plot_glm_continuous: show_rug = FALSE produces no GeomRug layer", {
  d <- make_cont_data()
  m <- stats::glm(y ~ x, family = stats::binomial, data = d)
  p <- plot_glm_continuous(m, d, "x", "y", show_rug = FALSE)
  layer_geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1L))
  expect_false("GeomRug" %in% layer_geoms)
})

test_that("plot_glm_continuous: show_rug = TRUE produces GeomRug layers", {
  d <- make_cont_data()
  m <- stats::glm(y ~ x, family = stats::binomial, data = d)
  p <- plot_glm_continuous(m, d, "x", "y", show_rug = TRUE)
  layer_geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1L))
  expect_true("GeomRug" %in% layer_geoms)
})

test_that("plot_glm_continuous multivariable model does not error", {
  d <- make_mixed_data()
  m <- stats::glm(y ~ cont + cat, family = stats::binomial, data = d)
  expect_no_error(plot_glm_continuous(m, d, "cont", "y"))
})

# ---- plot_glm_categorical ---------------------------------------------------

test_that("plot_glm_categorical probability metric returns a ggplot", {
  d <- make_cat_data()
  m <- stats::glm(y ~ g, family = stats::binomial, data = d)
  p <- plot_glm_categorical(m, d, "g", "y")
  expect_s3_class(p, "gg")
})

test_that("plot_glm_categorical odds_ratio metric returns a ggplot", {
  d <- make_cat_data()
  m <- stats::glm(y ~ g, family = stats::binomial, data = d)
  p <- plot_glm_categorical(m, d, "g", "y", metric = "odds_ratio")
  expect_s3_class(p, "gg")
})

test_that("plot_glm_categorical odds_ratio has a vline at x = 1", {
  d <- make_cat_data()
  m <- stats::glm(y ~ g, family = stats::binomial, data = d)
  p <- plot_glm_categorical(m, d, "g", "y", metric = "odds_ratio")
  layer_geoms <- vapply(p$layers, function(l) class(l$geom)[1L], character(1L))
  expect_true("GeomVline" %in% layer_geoms)
})

# ---- plot_glm_terms (plain glm) ---------------------------------------------

test_that("plot_glm_terms with glm + combine = TRUE returns a patchwork/gg", {
  d <- make_mixed_data()
  m <- stats::glm(y ~ cont + cat, family = stats::binomial, data = d)
  p <- plot_glm_terms(m, data = d, response_var = "y")
  expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
})

test_that("plot_glm_terms with glm + combine = FALSE returns a named list of ggplots", {
  d <- make_mixed_data()
  m <- stats::glm(y ~ cont + cat, family = stats::binomial, data = d)
  plots <- plot_glm_terms(m, data = d, response_var = "y", combine = FALSE)
  expect_type(plots, "list")
  expect_equal(length(plots), 2L)
  expect_equal(sort(names(plots)), c("cat", "cont"))
  expect_true(all(vapply(plots, function(x) inherits(x, "gg"), logical(1L))))
})

test_that("plot_glm_terms: interaction terms are skipped with a message", {
  d <- make_mixed_data()
  m <- stats::glm(y ~ cont * cat, family = stats::binomial, data = d)
  expect_message(
    plot_glm_terms(m, data = d, response_var = "y", combine = FALSE),
    "Interaction terms are not plotted individually"
  )
})

test_that("plot_glm_terms stops when no data supplied for glm input", {
  d <- make_cont_data()
  m <- stats::glm(y ~ x, family = stats::binomial, data = d)
  expect_error(plot_glm_terms(m, response_var = "y"), "`data` must be supplied")
})

test_that("plot_glm_terms: adjust_for subtitle annotation is added", {
  d <- make_mixed_data()
  m <- stats::glm(y ~ cont + cat, family = stats::binomial, data = d)
  plots <- plot_glm_terms(m, data = d, response_var = "y",
                           adjust_for = "cont", combine = FALSE)
  sub <- plots[["cont"]]$labels$subtitle
  expect_true(!is.null(sub) && grepl("Fixed adjustment", sub))
})

# ---- plot_glm_terms (scm_result) --------------------------------------------

test_that("plot_glm_terms with eres_scm_glm result does not error", {
  skip_on_cran()
  data("demo_eres", package = "erSCM")
  con_cov  <- c("WT", "AGE", "EGFR")
  exp_vars <- c("CMAXss", "CMINss", "AUCss")
  d <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_cov, exp_vars)),
                                ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX))

  result <- scm_glm(
    eres_data = d, covariates = c("WT", "AGE", "SEX"),
    exposure_vars = exp_vars, dv_col = "sim_status",
    significance_level = c(0.9, 0.9), center = FALSE, verbose = FALSE
  )

  expect_no_error(plot_glm_terms(result, combine = FALSE))
})

test_that("plot_glm_terms with eres_scm_glm + adjust_for annotates fixed terms", {
  skip_on_cran()
  data("demo_eres", package = "erSCM")
  con_cov  <- c("WT", "AGE", "EGFR")
  exp_vars <- c("CMAXss", "CMINss", "AUCss")
  d <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_cov, exp_vars)),
                                ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX))

  result <- scm_glm(
    eres_data = d, covariates = c("WT", "AGE", "SEX"),
    exposure_vars = exp_vars, dv_col = "sim_status",
    significance_level = c(0.9, 0.9), adjust_for = "AGE",
    center = FALSE, verbose = FALSE
  )

  plots <- plot_glm_terms(result, combine = FALSE)
  age_plot <- plots[["AGE"]]
  if (!is.null(age_plot)) {
    sub <- age_plot$labels$subtitle
    expect_true(!is.null(sub) && grepl("Fixed adjustment", sub))
  }
})
