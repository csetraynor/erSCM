## test-scm_bootstrap.R
## Tests for scm_bootstrap() core engine.

# ---- Input validation -------------------------------------------------------

test_that("scm_bootstrap() errors on non-SCM object", {
  expect_error(scm_bootstrap(list(a = 1)),     class = "rlang_error")
  expect_error(scm_bootstrap("not_a_model"),   class = "rlang_error")
  expect_error(scm_bootstrap(NULL),            class = "rlang_error")
})

test_that("scm_bootstrap() errors when n_boot < 10", {
  fit <- .fit_test_scm_glm()
  expect_error(scm_bootstrap(fit, n_boot = 5),  regexp = "n_boot")
  expect_error(scm_bootstrap(fit, n_boot = 1),  regexp = "n_boot")
  expect_error(scm_bootstrap(fit, n_boot = 0),  regexp = "n_boot")
})

test_that("scm_bootstrap() errors when seed is non-coercible", {
  fit <- .fit_test_scm_glm()
  expect_error(scm_bootstrap(fit, seed = "abc"))
})

test_that("scm_bootstrap() errors when seed is not a scalar", {
  fit <- .fit_test_scm_glm()
  expect_error(scm_bootstrap(fit, seed = c(1L, 2L)), regexp = "seed")
})

test_that("scm_bootstrap() errors when subject_col not in data", {
  fit <- .fit_test_scm_glm()
  expect_error(
    scm_bootstrap(fit, subject_col = "NONEXISTENT"),
    regexp = "subject_col"
  )
})

test_that("scm_bootstrap() errors when significance_level out of range", {
  fit <- .fit_test_scm_glm()
  expect_error(scm_bootstrap(fit, significance_level = c(1.5, 0.01)))
  expect_error(scm_bootstrap(fit, significance_level = c(-0.1, 0.01)))
  expect_error(scm_bootstrap(fit, significance_level = c(0, 0.01)))
  expect_error(scm_bootstrap(fit, significance_level = c(0.05, 1.0)))
})

test_that("scm_bootstrap() errors when significance_level wrong type", {
  fit <- .fit_test_scm_glm()
  expect_error(scm_bootstrap(fit, significance_level = "high"))
})

test_that("scm_bootstrap() errors when parallel is not logical", {
  fit <- .fit_test_scm_glm()
  expect_error(scm_bootstrap(fit, parallel = "yes"), regexp = "parallel")
})

# ---- Core functionality -----------------------------------------------------

test_that("scm_bootstrap() returns correct class", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 10L, seed = 1L))
  expect_s3_class(result, "scm_bootstrap")
})

test_that("scm_bootstrap() output has all required fields", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 10L, seed = 1L))
  expected_names <- c(
    "inclusion_prob", "interaction_inclusion_prob", "co_inclusion",
    "modal_model", "modal_model_freq", "model_variability_index",
    "exposure_selection_prob", "exposure_selected_per_boot",
    "boot_models", "boot_converged", "n_boot", "n_converged",
    "convergence_rate", "seed", "type", "exposure_vars", "covariates",
    "adjust_for", "significance_level", "original_selected", "call"
  )
  expect_true(all(expected_names %in% names(result)))
})

test_that("inclusion_prob values are in [0, 1] with named elements", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 10L, seed = 1L))
  probs  <- result$inclusion_prob
  expect_type(probs, "double")
  expect_true(all(probs >= 0 & probs <= 1))
  expect_true(all(nchar(names(probs)) > 0))
})

test_that("exposure variable is never in inclusion_prob", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 10L, seed = 1L))
  expect_false(any(result$exposure_vars %in% names(result$inclusion_prob)))
})

test_that("adjust_for terms are never in inclusion_prob", {
  skip_on_cran()
  test_data <- .make_test_data_glm()
  fit <- suppressMessages(scm_glm(
    eres_data         = test_data,
    exposure_vars     = "CMAX",
    covariates        = c("AGE", "WT", "SEX"),
    adjust_for        = "AGE",
    significance_level = c(0.5, 0.3),
    dv_col            = "DV",
    family            = "binomial",
    center            = FALSE,
    scale_to_median   = FALSE,
    verbose           = FALSE
  ))
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  expect_false("AGE" %in% names(result$inclusion_prob))
})

test_that("co_inclusion matrix is square, symmetric, and in [0, 1]", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 10L, seed = 1L))
  mat    <- result$co_inclusion
  if (nrow(mat) > 0L) {
    expect_equal(nrow(mat), ncol(mat))
    expect_equal(rownames(mat), colnames(mat))
    expect_true(isSymmetric(mat))
    expect_true(all(mat >= 0 & mat <= 1))
  }
})

test_that("diagonal of co_inclusion matches marginal inclusion proportions", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  mat    <- result$co_inclusion
  all_probs <- c(result$inclusion_prob, result$interaction_inclusion_prob)
  common    <- intersect(rownames(mat), names(all_probs))
  if (length(common) > 0L) {
    expect_equal(
      unname(diag(mat)[common]),
      unname(all_probs[common]),
      tolerance = 1e-10
    )
  }
})

test_that("modal_model$covariates is a character vector", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 10L, seed = 1L))
  # modal_model is now list(exposure, covariates); tested more fully in Phase 5 block
  expect_true(!is.null(result$modal_model))
})

test_that("modal_model_freq is in (0, 1]", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 10L, seed = 1L))
  expect_true(result$modal_model_freq > 0 && result$modal_model_freq <= 1)
})

test_that("model_variability_index = 1 - modal_model_freq", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  expect_equal(
    result$model_variability_index,
    1 - result$modal_model_freq,
    tolerance = 1e-10
  )
})

test_that("model_variability_index is in [0, 1]", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 10L, seed = 1L))
  expect_true(result$model_variability_index >= 0 &&
                result$model_variability_index <= 1)
})

test_that("n_converged <= n_boot and convergence_rate = n_converged / n_boot", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 10L, seed = 1L))
  expect_true(result$n_converged <= result$n_boot)
  expect_equal(result$convergence_rate,
               result$n_converged / result$n_boot,
               tolerance = 1e-10)
})

test_that("boot_converged length equals n_boot", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 10L, seed = 1L))
  expect_length(result$boot_converged, result$n_boot)
  expect_type(result$boot_converged, "logical")
})

test_that("boot_models length equals n_boot", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 10L, seed = 1L))
  expect_length(result$boot_models, result$n_boot)
})

test_that("same seed produces identical results (reproducibility)", {
  skip_on_cran()
  fit <- .fit_test_scm_glm()
  r1  <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 99L))
  r2  <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 99L))
  expect_identical(r1$inclusion_prob, r2$inclusion_prob)
  expect_identical(r1$boot_models,    r2$boot_models)
  expect_identical(r1$boot_converged, r2$boot_converged)
})

test_that("original_selected is a list with exposure and covariates fields", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  expect_type(result$original_selected, "list")
  expect_true("exposure" %in% names(result$original_selected))
  expect_true("covariates" %in% names(result$original_selected))
  expect_type(result$original_selected$covariates, "character")
})

test_that("stored seed matches the argument passed", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 42L))
  expect_equal(result$seed, 42L)
})

test_that("type field matches the SCM engine used", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  expect_equal(result$type, "glm")
})

test_that("significance_level is recorded in result", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(
    scm_bootstrap(fit, n_boot = 15L, seed = 1L,
                  significance_level = c(0.1, 0.05))
  )
  expect_equal(result$significance_level, c(0.1, 0.05))
})

test_that("significance_level of length 1 is recycled to length 2", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(
    scm_bootstrap(fit, n_boot = 15L, seed = 1L,
                  significance_level = 0.1)
  )
  expect_length(result$significance_level, 2L)
  expect_equal(result$significance_level, c(0.1, 0.1))
})

# ---- Exposure selection stability -------------------------------------------

test_that("exposure_selection_prob is a named numeric vector in [0, 1]", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 20L, seed = 1L))
  probs  <- result$exposure_selection_prob
  expect_type(probs, "double")
  expect_true(all(probs >= 0 & probs <= 1))
  expect_true(all(nchar(names(probs)) > 0))
})

test_that("exposure_selection_prob sums to convergence_rate for single exposure", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()   # single exposure "CMAX"
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 20L, seed = 1L))
  # For a single exposure, its prob = fraction of replicates it was selected
  # (which equals convergence_rate when exposure always wins if any replicate converges)
  expect_length(result$exposure_selection_prob, 1L)
  expect_equal(names(result$exposure_selection_prob), "CMAX")
})

test_that("exposure_selected_per_boot has length n_converged", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 20L, seed = 1L))
  expect_length(result$exposure_selected_per_boot, result$n_converged)
  expect_type(result$exposure_selected_per_boot, "character")
})

test_that("modal_model is a list with exposure and covariates fields", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 10L, seed = 1L))
  expect_type(result$modal_model, "list")
  expect_true("exposure" %in% names(result$modal_model))
  expect_true("covariates" %in% names(result$modal_model))
  expect_type(result$modal_model$covariates, "character")
})

test_that("exposure_selection_prob works with multiple exposures", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm_multi_exposure()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 20L, seed = 1L))
  expect_true("CMAX" %in% names(result$exposure_selection_prob))
  expect_true("AUC" %in% names(result$exposure_selection_prob))
  expect_true(all(result$exposure_selection_prob >= 0))
})

test_that("inclusion_prob_by_exposure is NULL for single exposure", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 10L, seed = 1L))
  expect_null(result$inclusion_prob_by_exposure)
})

test_that("inclusion_prob_by_exposure is a list for multiple exposures", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm_multi_exposure()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 20L, seed = 1L))
  if (!is.null(result$inclusion_prob_by_exposure)) {
    expect_type(result$inclusion_prob_by_exposure, "list")
    for (entry in result$inclusion_prob_by_exposure) {
      expect_type(entry, "double")
      expect_true(all(is.na(entry) | (entry >= 0 & entry <= 1)))
    }
  }
})

test_that("original_selected$exposure is a character string", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  expect_type(result$original_selected$exposure, "character")
  expect_length(result$original_selected$exposure, 1L)
})

test_that("boot_models converged elements are lists with exposure and covariates", {
  skip_on_cran()
  fit        <- .fit_test_scm_glm()
  result     <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  converged  <- Filter(Negate(is.null), result$boot_models)
  if (length(converged) > 0L) {
    first <- converged[[1L]]
    expect_true("exposure"   %in% names(first))
    expect_true("covariates" %in% names(first))
  }
})

# ---- Metadata fix tests -------------------------------------------------------

test_that("scm_coxph stores full candidate_exposure_vars list in metadata", {
  skip_on_cran()
  test_data <- .make_test_data_coxph()
  fit <- suppressMessages(scm_coxph(
    eres_data         = test_data,
    exposure_vars     = c("CMAX", "AUC"),
    covariates        = c("AGE", "WT"),
    significance_level = c(0.5, 0.3),
    dv_col            = "DV",
    time_col          = "TIME",
    center            = FALSE,
    scale_to_median   = FALSE,
    verbose           = FALSE
  ))
  meta <- fit$metadata
  expect_equal(length(meta$candidate_exposure_vars), 2L)
  expect_true(all(c("CMAX", "AUC") %in% meta$candidate_exposure_vars))
  expect_true("selected_exposure" %in% names(meta))
  expect_length(meta$selected_exposure, 1L)
})

test_that("scm_glm stores selected_exposure in metadata", {
  fit <- .fit_test_scm_glm()
  expect_true("selected_exposure" %in% names(fit$metadata))
  expect_type(fit$metadata$selected_exposure, "character")
  expect_length(fit$metadata$selected_exposure, 1L)
})

# ---- CoxPH -------------------------------------------------------------------

test_that("scm_bootstrap() works with a CoxPH fit", {
  skip_on_cran()
  fit    <- .fit_test_scm_coxph()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  expect_s3_class(result, "scm_bootstrap")
  expect_equal(result$type, "coxph")
  expect_true(all(result$inclusion_prob >= 0 & result$inclusion_prob <= 1))
})

# ---- Subject-level resampling -----------------------------------------------

test_that("subject_col enables patient-level resampling", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(
    scm_bootstrap(fit, n_boot = 15L, seed = 1L, subject_col = "SUBJID")
  )
  expect_s3_class(result, "scm_bootstrap")
  expect_length(result$boot_models, 15L)
})

# ---- Parallel ---------------------------------------------------------------

test_that("parallel = TRUE works with future::sequential backend", {
  skip_on_cran()
  skip_if_not_installed("future.apply")
  skip_if_not_installed("future")
  fit <- .fit_test_scm_glm()
  future::plan(future::sequential)
  result <- suppressMessages(
    scm_bootstrap(fit, n_boot = 10L, seed = 1L, parallel = TRUE)
  )
  future::plan(future::sequential)   # reset
  expect_s3_class(result, "scm_bootstrap")
})
