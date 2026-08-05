## test-scm_bootstrap_integration.R
## End-to-end integration tests: full GLM/CoxPH fit → bootstrap → methods chain.

# ---- GLM full pipeline -------------------------------------------------------

test_that("GLM full pipeline: fit → bootstrap → print → summary → plot", {
  skip_on_cran()

  fit <- .fit_test_scm_glm()
  expect_s3_class(fit, "eres_scm_glm")

  stability <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  expect_s3_class(stability, "scm_bootstrap")

  # print runs without error
  expect_output(print(stability))

  # summary returns correct class and fields
  s <- summary(stability)
  expect_s3_class(s, "summary_scm_bootstrap")
  expect_true(all(c("inclusion_table", "stable_covariates",
                    "modal_model", "n_converged") %in% names(s)))

  # plot — inclusion
  p <- plot(stability, type = "inclusion")
  expect_s3_class(p, "ggplot")

  # stable_covariates
  sc <- stable_covariates(stability, threshold = 0.3)
  expect_type(sc, "character")
})

# ---- CoxPH full pipeline -----------------------------------------------------

test_that("CoxPH full pipeline: fit → bootstrap → print → summary → plot", {
  skip_on_cran()

  fit <- .fit_test_scm_coxph()
  expect_s3_class(fit, "eres_scm_coxph")

  stability <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 2L))
  expect_s3_class(stability, "scm_bootstrap")
  expect_equal(stability$type, "coxph")

  expect_output(print(stability))

  s <- summary(stability)
  expect_s3_class(s, "summary_scm_bootstrap")
  expect_equal(s$type, "coxph")

  p <- plot(stability, type = "inclusion")
  expect_s3_class(p, "ggplot")
})

# ---- Multi-exposure full pipeline --------------------------------------------

test_that("multi-exposure pipeline: fit → bootstrap → exposure plot → summary tibble", {
  skip_on_cran()

  fit <- .fit_test_scm_glm_multi_exposure()
  expect_s3_class(fit, "eres_scm_glm")

  stability <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 3L))
  expect_s3_class(stability, "scm_bootstrap")

  # Exposure selection fields populated
  expect_true("CMAX" %in% names(stability$exposure_selection_prob) ||
                "AUC" %in% names(stability$exposure_selection_prob))

  # Exposure plot
  p <- plot(stability, type = "exposure")
  expect_s3_class(p, "ggplot")

  # summary exposure_stability tibble
  s <- summary(stability)
  expect_true("exposure_stability" %in% names(s))
  expect_true(is.data.frame(s$exposure_stability))
  expect_true(all(c("exposure_var", "selection_prob", "is_original") %in%
                    names(s$exposure_stability)))
})

# ---- Reproducibility ---------------------------------------------------------

test_that("same seed → identical results across two bootstrap calls", {
  skip_on_cran()

  fit <- .fit_test_scm_glm()
  r1  <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 77L))
  r2  <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 77L))
  expect_identical(r1$inclusion_prob,            r2$inclusion_prob)
  expect_identical(r1$exposure_selection_prob,   r2$exposure_selection_prob)
  expect_identical(r1$boot_models,               r2$boot_models)
  expect_identical(r1$boot_converged,            r2$boot_converged)
})

# ---- Structural invariants ---------------------------------------------------

test_that("inclusion_prob values are all in [0, 1]", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 4L))
  expect_true(all(result$inclusion_prob >= 0 & result$inclusion_prob <= 1))
})

test_that("exposure_selection_prob values are all in [0, 1]", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 4L))
  expect_true(all(result$exposure_selection_prob >= 0 &
                    result$exposure_selection_prob <= 1))
})

test_that("model_variability_index = 1 - modal_model_freq (end-to-end)", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 5L))
  expect_equal(result$model_variability_index, 1 - result$modal_model_freq,
               tolerance = 1e-10)
})

test_that("co_inclusion diagonal equals marginal inclusion_prob (end-to-end)", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 6L))
  mat    <- result$co_inclusion
  probs  <- c(result$inclusion_prob, result$interaction_inclusion_prob)
  common <- intersect(rownames(mat), names(probs))
  if (length(common) > 0L) {
    expect_equal(unname(diag(mat)[common]), unname(probs[common]),
                 tolerance = 1e-10)
  }
})

test_that("n_converged <= n_boot and boot_converged length = n_boot (end-to-end)", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 7L))
  expect_true(result$n_converged <= result$n_boot)
  expect_length(result$boot_converged, result$n_boot)
})

# ---- original_selected list format -------------------------------------------

test_that("original_selected is a list with exposure and covariates (end-to-end)", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 8L))
  expect_type(result$original_selected, "list")
  expect_true("exposure"   %in% names(result$original_selected))
  expect_true("covariates" %in% names(result$original_selected))
  expect_type(result$original_selected$exposure,   "character")
  expect_type(result$original_selected$covariates, "character")
  expect_length(result$original_selected$exposure, 1L)
})

# ---- modal_model list format -------------------------------------------------

test_that("modal_model is a list with exposure and covariates (end-to-end)", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 9L))
  expect_type(result$modal_model, "list")
  expect_true("exposure"   %in% names(result$modal_model))
  expect_true("covariates" %in% names(result$modal_model))
})

# ---- stable_covariates excludes exposure by default -------------------------

test_that("stable_covariates(include_exposure=FALSE) never contains exposure_vars", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 10L))
  sc     <- stable_covariates(result, threshold = 0)
  expect_false(any(result$exposure_vars %in% sc))
})

# ---- print.summary_scm_bootstrap round-trips --------------------------------

test_that("print(summary(stability)) produces visible output without error", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 11L))
  s      <- summary(result)
  expect_output(print(s))
})

# ---- coinclusion plot --------------------------------------------------------

test_that("plot(type='coinclusion') returns a ggplot (end-to-end)", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 12L))
  p      <- plot(result, type = "coinclusion")
  expect_s3_class(p, "ggplot")
})
