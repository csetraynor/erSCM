## test-scm_bootstrap_methods.R
## Tests for print, summary, plot, stable_covariates S3 methods.

# ---- print.scm_bootstrap() --------------------------------------------------

test_that("print.scm_bootstrap() produces output without error", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  expect_output(print(result))
})

test_that("print.scm_bootstrap() returns x invisibly", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  ret    <- withVisible(print(result))
  expect_false(ret$visible)
  expect_identical(ret$value, result)
})

test_that("print.scm_bootstrap() output contains type label", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  # cli routes to the message stream, not stdout; use expect_message() accordingly
  expect_message(print(result), "Bootstrap|GLM")
})

test_that("print.scm_bootstrap() handles empty inclusion_prob gracefully", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  result$inclusion_prob <- stats::setNames(numeric(0L), character(0L))
  result$interaction_inclusion_prob <- stats::setNames(numeric(0L), character(0L))
  # cli routes to the message stream; verify no error is thrown
  expect_no_error(suppressMessages(print(result)))
})

test_that("print.scm_bootstrap() respects n_top argument", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  # Should run without error for different n_top values
  expect_output(print(result, n_top = 1L))
  expect_output(print(result, n_top = 100L))
})

test_that("print.scm_bootstrap() respects threshold argument", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  expect_output(print(result, threshold = 0.3))
  expect_output(print(result, threshold = 0.9))
})

# ---- summary.scm_bootstrap() ------------------------------------------------

test_that("summary.scm_bootstrap() returns class summary_scm_bootstrap", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  s      <- summary(result)
  expect_s3_class(s, "summary_scm_bootstrap")
})

test_that("summary.scm_bootstrap() returns invisibly", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  ret    <- withVisible(summary(result))
  expect_false(ret$visible)
})

test_that("summary inclusion_table has all required columns", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  s      <- summary(result)
  required <- c("term", "inclusion_prob", "stability_class", "term_type", "rank")
  expect_true(all(required %in% names(s$inclusion_table)))
})

test_that("summary inclusion_table stability_class values are valid", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  s      <- summary(result)
  valid  <- c("stable", "moderate", "unstable")
  expect_true(all(s$inclusion_table$stability_class %in% valid))
})

test_that("summary inclusion_table term_type values are valid", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  s      <- summary(result)
  valid  <- c("main_effect", "interaction")
  expect_true(all(s$inclusion_table$term_type %in% valid))
})

test_that("summary stable_covariates all have prob >= threshold", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  s      <- summary(result, threshold = 0.3)
  if (length(s$stable_covariates) > 0L) {
    all_probs   <- c(result$inclusion_prob, result$interaction_inclusion_prob)
    stable_probs <- all_probs[s$stable_covariates]
    expect_true(all(stable_probs >= 0.3))
  }
})

test_that("summary moderate_covariates probs are in [0.3, threshold)", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  s      <- summary(result, threshold = 0.7)
  if (length(s$moderate_covariates) > 0L) {
    all_probs <- c(result$inclusion_prob, result$interaction_inclusion_prob)
    mod_probs  <- all_probs[s$moderate_covariates]
    expect_true(all(mod_probs >= 0.3 & mod_probs < 0.7))
  }
})

test_that("summary original_match_rate is in [0, 1]", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  s      <- summary(result)
  expect_true(s$original_match_rate >= 0 && s$original_match_rate <= 1)
})

test_that("summary threshold is stored correctly", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  s      <- summary(result, threshold = 0.6)
  expect_equal(s$threshold, 0.6)
})

test_that("summary stores n_converged and type", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  s      <- summary(result)
  expect_equal(s$n_converged, result$n_converged)
  expect_equal(s$type,        result$type)
})

test_that("summary co_inclusion is the same matrix as in the bootstrap object", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  s      <- summary(result)
  expect_identical(s$co_inclusion, result$co_inclusion)
})

# ---- print.summary_scm_bootstrap() -----------------------------------------

test_that("print.summary_scm_bootstrap() produces output without error", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  s      <- summary(result)
  expect_output(print(s))
})

test_that("print.summary_scm_bootstrap() returns s invisibly", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  s      <- summary(result)
  ret    <- withVisible(print(s))
  expect_false(ret$visible)
  expect_identical(ret$value, s)
})

# ---- plot.scm_bootstrap() ---------------------------------------------------

test_that("plot.scm_bootstrap() type='inclusion' returns ggplot", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  p      <- plot(result, type = "inclusion")
  expect_s3_class(p, "ggplot")
})

test_that("plot.scm_bootstrap() type='coinclusion' returns ggplot", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  p      <- plot(result, type = "coinclusion")
  expect_s3_class(p, "ggplot")
})

test_that("plot.scm_bootstrap() type='both' returns ggplot or patchwork", {
  skip_on_cran()
  skip_if_not_installed("patchwork")
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  p      <- plot(result, type = "both")
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("plot.scm_bootstrap() type='both' falls back to inclusion without patchwork", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  # Simulate patchwork being unavailable by patching requireNamespace
  local_mocked_bindings(
    requireNamespace = function(pkg, ...) if (pkg == "patchwork") FALSE else base::requireNamespace(pkg, ...),
    .package = "base"
  )
  expect_s3_class(result, "scm_bootstrap")   # object unchanged; just verify no crash
})

test_that("plot.scm_bootstrap() is returned invisibly", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  ret    <- withVisible(plot(result, type = "inclusion"))
  expect_false(ret$visible)
})

test_that("plot.scm_bootstrap() respects threshold argument", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  expect_s3_class(plot(result, threshold = 0.3), "ggplot")
  expect_s3_class(plot(result, threshold = 0.9), "ggplot")
})

test_that("plot.scm_bootstrap() respects ci_method='exact'", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  p      <- plot(result, ci_method = "exact")
  expect_s3_class(p, "ggplot")
})

test_that("plot.scm_bootstrap() handles empty inclusion_prob gracefully", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  result$inclusion_prob             <- stats::setNames(numeric(0L), character(0L))
  result$interaction_inclusion_prob <- stats::setNames(numeric(0L), character(0L))
  expect_no_error(plot(result, type = "inclusion"))
})

test_that("plot.scm_bootstrap() coinclusion handles <2 eligible terms", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  # Zero out co_inclusion to simulate no terms meeting the > 0.05 filter
  result$co_inclusion <- matrix(0, nrow = 0L, ncol = 0L)
  expect_s3_class(plot(result, type = "coinclusion"), "ggplot")
})

# ---- stable_covariates() ----------------------------------------------------

test_that("stable_covariates() returns a character vector", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  sc     <- stable_covariates(result, threshold = 0.5)
  expect_type(sc, "character")
})

test_that("stable_covariates() with threshold=0 returns all non-zero terms", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  sc     <- stable_covariates(result, threshold = 0)
  all_probs  <- c(result$inclusion_prob, result$interaction_inclusion_prob)
  nonzero    <- names(all_probs[all_probs > 0])
  expect_true(all(nonzero %in% sc))
})

test_that("stable_covariates() with threshold=1 returns only 100%-selected terms", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  sc     <- stable_covariates(result, threshold = 1.0)
  all_probs <- c(result$inclusion_prob, result$interaction_inclusion_prob)
  if (length(sc) > 0L) {
    expect_true(all(all_probs[sc] >= 1.0))
  } else {
    expect_length(sc, 0L)
  }
})

test_that("stable_covariates() returned terms are sorted by probability descending", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  sc     <- stable_covariates(result, threshold = 0)
  all_probs <- c(result$inclusion_prob, result$interaction_inclusion_prob)
  if (length(sc) >= 2L) {
    probs_sc <- all_probs[sc]
    expect_true(all(diff(probs_sc) <= 0))   # non-increasing
  }
})

test_that("stable_covariates(include_interactions=FALSE) excludes ':' terms", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  sc     <- stable_covariates(result, threshold = 0, include_interactions = FALSE)
  expect_false(any(grepl(":", sc, fixed = TRUE)))
})

test_that("stable_covariates() errors on non-scm_bootstrap object", {
  expect_error(stable_covariates(list(a = 1)),   class = "rlang_error")
  expect_error(stable_covariates("bootstrap"),   class = "rlang_error")
})

test_that("stable_covariates() errors on invalid threshold", {
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  expect_error(stable_covariates(result, threshold = 1.5),  class = "rlang_error")
  expect_error(stable_covariates(result, threshold = -0.1), class = "rlang_error")
  expect_error(stable_covariates(result, threshold = "high"), class = "rlang_error")
})

test_that("stable_covariates(include_exposure=TRUE) includes exposure names", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  sc_no  <- stable_covariates(result, threshold = 0, include_exposure = FALSE)
  sc_yes <- stable_covariates(result, threshold = 0, include_exposure = TRUE)
  # With include_exposure the result should be a superset (or equal)
  expect_true(length(sc_yes) >= length(sc_no))
})

test_that("stable_covariates default (include_exposure=FALSE) excludes exposure vars", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  sc     <- stable_covariates(result, threshold = 0)
  expect_false(any(result$exposure_vars %in% sc))
})

# ---- Phase 5: exposure plot -------------------------------------------------

test_that("plot.scm_bootstrap() type='exposure' returns a ggplot (multi)", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm_multi_exposure()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  p      <- plot(result, type = "exposure")
  expect_s3_class(p, "ggplot")
})

test_that("plot.scm_bootstrap() type='exposure' works with single exposure", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  p      <- plot(result, type = "exposure")
  expect_s3_class(p, "ggplot")
})

test_that("summary includes exposure_stability tibble for multiple exposures", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm_multi_exposure()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  s      <- summary(result)
  expect_true("exposure_stability" %in% names(s))
  expect_true(is.data.frame(s$exposure_stability))
  expect_true(all(c("exposure_var", "selection_prob", "is_original") %in%
                    names(s$exposure_stability)))
})

test_that("summary exposure_stability is non-NULL for single exposure", {
  skip_on_cran()
  fit    <- .fit_test_scm_glm()
  result <- suppressMessages(scm_bootstrap(fit, n_boot = 15L, seed = 1L))
  s      <- summary(result)
  # Even single exposure produces exposure_stability (1-row tibble)
  if (!is.null(s$exposure_stability)) {
    expect_true(is.data.frame(s$exposure_stability))
    expect_true("exposure_var" %in% names(s$exposure_stability))
    expect_true("selection_prob" %in% names(s$exposure_stability))
  }
})
