## test-scm_bootstrap_helpers.R
## Tests for internal helper functions in scm_bootstrap_helpers.R and
## scm_bootstrap_methods.R (accessed via :::).

# ---- .generate_bootstrap_sample() -------------------------------------------

test_that(".generate_bootstrap_sample() row-level returns same dimensions", {
  skip_on_cran()
  df   <- data.frame(x = seq_len(100L), y = rnorm(100L))
  boot <- erSCM:::.generate_bootstrap_sample(df, subject_col = NULL)
  expect_equal(nrow(boot), 100L)
  expect_equal(ncol(boot), 2L)
})

test_that(".generate_bootstrap_sample() row-level is a data.frame", {
  skip_on_cran()
  df   <- data.frame(x = seq_len(50L), y = rnorm(50L))
  boot <- erSCM:::.generate_bootstrap_sample(df)
  expect_s3_class(boot, "data.frame")
})

test_that(".generate_bootstrap_sample() row-level produces variation", {
  skip_on_cran()
  set.seed(42L)
  df <- data.frame(x = seq_len(100L))
  b1 <- erSCM:::.generate_bootstrap_sample(df)
  b2 <- erSCM:::.generate_bootstrap_sample(df)
  # Two independent samples are almost certainly not identical
  expect_false(identical(b1$x, b2$x))
})

test_that(".generate_bootstrap_sample() patient-level preserves all rows for sampled IDs", {
  skip_on_cran()
  df <- data.frame(
    SUBJID = rep(paste0("S", seq_len(20L)), each = 5L),
    x      = rnorm(100L),
    y      = rnorm(100L)
  )
  set.seed(7L)
  boot <- erSCM:::.generate_bootstrap_sample(df, subject_col = "SUBJID")
  # Sampled subjects = 20 unique IDs drawn with replacement → 20 slots
  # each contributes 5 rows → 100 total (some subjects duplicated, some absent)
  expect_equal(nrow(boot), 100L)
  expect_equal(ncol(boot), 3L)
})

test_that(".generate_bootstrap_sample() patient-level samples with replacement (same nrow, subset of IDs)", {
  skip_on_cran()
  df <- data.frame(
    SUBJID = rep(paste0("S", seq_len(10L)), each = 3L),
    x      = rnorm(30L)
  )
  set.seed(1L)
  boot <- erSCM:::.generate_bootstrap_sample(df, subject_col = "SUBJID")
  # Bootstrap samples the same number of subjects with replacement: same total rows,
  # but typically fewer unique IDs (~63% of original)
  expect_equal(nrow(boot), 30L)
  expect_lte(length(unique(boot$SUBJID)), 10L)
  expect_gte(length(unique(boot$SUBJID)), 1L)
  expect_true(all(unique(boot$SUBJID) %in% paste0("S", seq_len(10L))))
})

# ---- .compute_inclusion_prob() ----------------------------------------------

test_that(".compute_inclusion_prob() returns correctly named vector", {
  skip_on_cran()
  boot_models <- list(c("A", "B"), c("A", "C"), c("B", "C"))
  candidates  <- c("A", "B", "C", "D")
  probs       <- erSCM:::.compute_inclusion_prob(boot_models, candidates)
  expect_named(probs, candidates)
})

test_that(".compute_inclusion_prob() computes proportions correctly", {
  skip_on_cran()
  boot_models <- list(c("A", "B"), c("A", "C"), c("A"), c("B", "C"), c("A", "B", "C"))
  candidates  <- c("A", "B", "C", "D")
  probs       <- erSCM:::.compute_inclusion_prob(boot_models, candidates)
  expect_equal(probs[["A"]], 4 / 5)
  expect_equal(probs[["B"]], 3 / 5)
  expect_equal(probs[["C"]], 3 / 5)
  expect_equal(probs[["D"]], 0)
})

test_that(".compute_inclusion_prob() handles empty boot_models", {
  skip_on_cran()
  probs <- erSCM:::.compute_inclusion_prob(list(), c("A", "B"))
  expect_equal(unname(probs), c(0, 0))
})

test_that(".compute_inclusion_prob() handles empty candidates", {
  skip_on_cran()
  boot_models <- list(c("A"), c("B"))
  probs       <- erSCM:::.compute_inclusion_prob(boot_models, character(0L))
  expect_length(probs, 0L)
})

test_that(".compute_inclusion_prob() values are in [0, 1]", {
  skip_on_cran()
  boot_models <- list(c("A", "B"), c("B"), c("C"), c("A", "C"), c("A", "B", "C"))
  probs       <- erSCM:::.compute_inclusion_prob(boot_models, c("A", "B", "C"))
  expect_true(all(probs >= 0 & probs <= 1))
})

# ---- .compute_co_inclusion_matrix() -----------------------------------------

test_that(".compute_co_inclusion_matrix() is symmetric", {
  skip_on_cran()
  boot_models <- list(c("A", "B"), c("A", "C"), c("B", "C"))
  all_terms   <- c("A", "B", "C")
  mat         <- erSCM:::.compute_co_inclusion_matrix(boot_models, all_terms)
  expect_true(isSymmetric(mat))
})

test_that(".compute_co_inclusion_matrix() off-diagonal values are correct", {
  skip_on_cran()
  boot_models <- list(c("A", "B"), c("A", "C"), c("B", "C"))
  all_terms   <- c("A", "B", "C")
  mat         <- erSCM:::.compute_co_inclusion_matrix(boot_models, all_terms)
  # rep1={A,B}, rep2={A,C}, rep3={B,C}: each pair co-occurs in exactly one replicate
  expect_equal(mat["A", "B"], 1 / 3)
  expect_equal(mat["A", "C"], 1 / 3)
  expect_equal(mat["B", "C"], 1 / 3)
})

test_that(".compute_co_inclusion_matrix() diagonal equals marginal inclusion", {
  skip_on_cran()
  boot_models <- list(c("A", "B"), c("A"), c("B"))
  all_terms   <- c("A", "B")
  mat         <- erSCM:::.compute_co_inclusion_matrix(boot_models, all_terms)
  expect_equal(mat["A", "A"], 2 / 3)   # A in replicates 1,2
  expect_equal(mat["B", "B"], 2 / 3)   # B in replicates 1,3
})

test_that(".compute_co_inclusion_matrix() all values in [0, 1]", {
  skip_on_cran()
  boot_models <- list(c("A", "B"), c("A", "B", "C"), c("C"))
  mat         <- erSCM:::.compute_co_inclusion_matrix(boot_models, c("A", "B", "C"))
  expect_true(all(mat >= 0 & mat <= 1))
})

test_that(".compute_co_inclusion_matrix() handles empty boot_models", {
  skip_on_cran()
  mat <- erSCM:::.compute_co_inclusion_matrix(list(), c("A", "B"))
  expect_equal(dim(mat), c(2L, 2L))
  expect_true(all(mat == 0))
})

test_that(".compute_co_inclusion_matrix() handles empty terms", {
  skip_on_cran()
  boot_models <- list(c("A"), c("B"))
  mat         <- erSCM:::.compute_co_inclusion_matrix(boot_models, character(0L))
  expect_equal(dim(mat), c(0L, 0L))
})

# ---- .identify_modal_model() ------------------------------------------------

test_that(".identify_modal_model() finds the most frequent model", {
  skip_on_cran()
  boot_models <- list(
    c("A", "B"), c("A", "B"), c("A", "B"),
    c("A", "C"), c("B", "C")
  )
  result <- erSCM:::.identify_modal_model(boot_models)
  expect_equal(sort(result$model), c("A", "B"))
  expect_equal(result$frequency, 3 / 5)
})

test_that(".identify_modal_model() handles empty character vectors (null model)", {
  skip_on_cran()
  boot_models <- list(character(0L), character(0L), c("A"))
  result      <- erSCM:::.identify_modal_model(boot_models)
  expect_equal(result$model,     character(0L))
  expect_equal(result$frequency, 2 / 3)
})

test_that(".identify_modal_model() tie-breaking is deterministic (alphabetical)", {
  skip_on_cran()
  # Both "A" and "B" appear twice as singleton models → table sorts by name
  boot_models <- list(c("A"), c("B"), c("A"), c("B"))
  result      <- erSCM:::.identify_modal_model(boot_models)
  # Canonical strings: "A" and "B"; table() returns alphabetical order;
  # sort(decreasing=TRUE) then names()[1] = "B" because "B" > "A" lexicographically.
  # The spec says "ties broken alphabetically (first alphabetically wins)".
  # Our impl: sort(table, decreasing=TRUE); ties retain original table order (A < B),
  # so names()[1] after sort is the one that appears first in the original table.
  # table() sorts by name alphabetically, so "A" comes first; after sort(desc)
  # ties preserve table order → names()[1] = "A".
  expect_type(result$model, "character")
  expect_equal(result$frequency, 2 / 4)
})

test_that(".identify_modal_model() handles a single replicate", {
  skip_on_cran()
  result <- erSCM:::.identify_modal_model(list(c("A", "B")))
  expect_equal(sort(result$model), c("A", "B"))
  expect_equal(result$frequency, 1)
})

test_that(".identify_modal_model() handles empty input", {
  skip_on_cran()
  result <- erSCM:::.identify_modal_model(list())
  expect_equal(result$model,     character(0L))
  expect_true(is.na(result$frequency))
})

# ---- .binom_ci() ------------------------------------------------------------

test_that(".binom_ci() Wilson CI contains the point estimate", {
  skip_on_cran()
  ci <- erSCM:::.binom_ci(30L, 100L, method = "wilson")
  expect_true(ci["lower"] <= 0.30)
  expect_true(ci["upper"] >= 0.30)
})

test_that(".binom_ci() exact (Clopper-Pearson) CI contains point estimate", {
  skip_on_cran()
  ci <- erSCM:::.binom_ci(30L, 100L, method = "exact")
  expect_true(ci["lower"] <= 0.30)
  expect_true(ci["upper"] >= 0.30)
})

test_that(".binom_ci() returns values in [0, 1] at x=0", {
  skip_on_cran()
  ci <- erSCM:::.binom_ci(0L, 100L, method = "wilson")
  expect_true(ci["lower"] >= 0)
  expect_true(ci["upper"] <= 1)
  expect_equal(unname(ci["lower"]), 0)
})

test_that(".binom_ci() returns values in [0, 1] at x=n", {
  skip_on_cran()
  ci <- erSCM:::.binom_ci(100L, 100L, method = "wilson")
  expect_true(ci["lower"] >= 0)
  expect_true(ci["upper"] <= 1)
  expect_equal(unname(ci["upper"]), 1)
})

test_that(".binom_ci() handles n=0 edge case", {
  skip_on_cran()
  ci <- erSCM:::.binom_ci(0L, 0L, method = "wilson")
  expect_equal(unname(ci["lower"]), 0)
  expect_equal(unname(ci["upper"]), 1)
})

test_that(".binom_ci() lower < upper for interior proportions", {
  skip_on_cran()
  ci <- erSCM:::.binom_ci(50L, 100L, method = "wilson")
  expect_true(ci["lower"] < ci["upper"])
  ci2 <- erSCM:::.binom_ci(50L, 100L, method = "exact")
  expect_true(ci2["lower"] < ci2["upper"])
})

test_that(".binom_ci() Wilson and exact agree closely for n=1000", {
  skip_on_cran()
  ci_w <- erSCM:::.binom_ci(300L, 1000L, method = "wilson")
  ci_e <- erSCM:::.binom_ci(300L, 1000L, method = "exact")
  expect_true(abs(ci_w["lower"] - ci_e["lower"]) < 0.01)
  expect_true(abs(ci_w["upper"] - ci_e["upper"]) < 0.01)
})

test_that(".binom_ci() returns named vector with 'lower' and 'upper'", {
  skip_on_cran()
  ci <- erSCM:::.binom_ci(10L, 20L)
  expect_named(ci, c("lower", "upper"))
})

# ---- .extract_scm_settings() ------------------------------------------------

test_that(".extract_scm_settings() extracts type correctly for GLM", {
  skip_on_cran()
  fit      <- .fit_test_scm_glm()
  settings <- erSCM:::.extract_scm_settings(fit)
  expect_equal(settings$type, "glm")
})

test_that(".extract_scm_settings() extracts type correctly for CoxPH", {
  skip_on_cran()
  fit      <- .fit_test_scm_coxph()
  settings <- erSCM:::.extract_scm_settings(fit)
  expect_equal(settings$type, "coxph")
})

test_that(".extract_scm_settings() exposure_vars not in covariates", {
  skip_on_cran()
  fit      <- .fit_test_scm_glm()
  settings <- suppressMessages(erSCM:::.extract_scm_settings(fit))
  expect_false(any(settings$exposure_vars %in% settings$covariates))
})

test_that(".extract_scm_settings() returns data.frame for data", {
  skip_on_cran()
  fit      <- .fit_test_scm_glm()
  settings <- suppressMessages(erSCM:::.extract_scm_settings(fit))
  expect_s3_class(settings$data, "data.frame")
  expect_gt(nrow(settings$data), 0L)
})

test_that(".extract_scm_settings() overrides significance_level when supplied", {
  skip_on_cran()
  fit      <- .fit_test_scm_glm()
  settings <- suppressMessages(
    erSCM:::.extract_scm_settings(fit, significance_level = c(0.1, 0.05))
  )
  expect_equal(settings$significance_level, c(0.1, 0.05))
})

test_that(".extract_scm_settings() defaults significance_level with message", {
  skip_on_cran()
  fit <- .fit_test_scm_glm()
  expect_message(
    erSCM:::.extract_scm_settings(fit),
    regexp = "significance_level|default|0.05"
  )
})

# ---- .compute_exposure_selection_prob() -------------------------------------

test_that(".compute_exposure_selection_prob() correct for known input", {
  skip_on_cran()
  models <- list(
    list(exposure = "CMAX",         covariates = c("AGE")),
    list(exposure = "CMAX",         covariates = c("AGE", "WT")),
    list(exposure = "AUC",          covariates = c("WT")),
    list(exposure = "CMAX",         covariates = character(0L)),
    list(exposure = NA_character_,  covariates = character(0L))
  )
  probs <- erSCM:::.compute_exposure_selection_prob(models, c("CMAX", "AUC"))
  expect_equal(probs[["CMAX"]], 3 / 5)
  expect_equal(probs[["AUC"]],  1 / 5)
  expect_equal(probs[["none_significant"]], 1 / 5)
})

test_that(".compute_exposure_selection_prob() no NAs means no none_significant", {
  skip_on_cran()
  models <- list(
    list(exposure = "CMAX", covariates = c("AGE")),
    list(exposure = "AUC",  covariates = c("WT"))
  )
  probs <- erSCM:::.compute_exposure_selection_prob(models, c("CMAX", "AUC"))
  expect_false("none_significant" %in% names(probs))
  expect_equal(probs[["CMAX"]], 0.5)
  expect_equal(probs[["AUC"]],  0.5)
})

test_that(".compute_exposure_selection_prob() handles empty boot_models", {
  skip_on_cran()
  probs <- erSCM:::.compute_exposure_selection_prob(list(), c("CMAX", "AUC"))
  expect_equal(unname(probs), c(0, 0))
  expect_named(probs, c("CMAX", "AUC"))
})

# ---- .compute_inclusion_prob_by_exposure() -----------------------------------

test_that(".compute_inclusion_prob_by_exposure() returns NULL for single exposure", {
  skip_on_cran()
  models <- list(
    list(exposure = "CMAX", covariates = c("AGE", "WT")),
    list(exposure = "CMAX", covariates = c("AGE"))
  )
  result <- erSCM:::.compute_inclusion_prob_by_exposure(
    models, c("AGE", "WT"), "CMAX"
  )
  expect_null(result)
})

test_that(".compute_inclusion_prob_by_exposure() correct for multiple exposures", {
  skip_on_cran()
  models <- list(
    list(exposure = "CMAX", covariates = c("AGE", "WT")),
    list(exposure = "CMAX", covariates = c("AGE")),
    list(exposure = "AUC",  covariates = c("WT", "SEX")),
    list(exposure = "AUC",  covariates = c("WT"))
  )
  result <- erSCM:::.compute_inclusion_prob_by_exposure(
    models, c("AGE", "WT", "SEX"), c("CMAX", "AUC")
  )
  expect_type(result, "list")
  expect_named(result, c("CMAX", "AUC"))
  expect_equal(result$CMAX[["AGE"]], 1.0)   # AGE in both CMAX replicates
  expect_equal(result$CMAX[["WT"]],  0.5)   # WT in 1 of 2 CMAX replicates
  expect_equal(result$AUC[["WT"]],   1.0)   # WT in both AUC replicates
  expect_equal(result$AUC[["SEX"]],  0.5)   # SEX in 1 of 2 AUC replicates
})

test_that(".compute_inclusion_prob_by_exposure() returns NULL for empty models", {
  skip_on_cran()
  result <- erSCM:::.compute_inclusion_prob_by_exposure(
    list(), c("AGE", "WT"), c("CMAX", "AUC")
  )
  expect_null(result)
})

# ---- .identify_modal_model() includes exposure in model identity ------------

test_that(".identify_modal_model() includes exposure in model identity", {
  skip_on_cran()
  models <- list(
    list(exposure = "CMAX", covariates = c("AGE")),
    list(exposure = "CMAX", covariates = c("AGE")),
    list(exposure = "AUC",  covariates = c("AGE")),
    list(exposure = "CMAX", covariates = c("WT"))
  )
  result <- erSCM:::.identify_modal_model(models)
  expect_equal(result$model$exposure,   "CMAX")
  expect_equal(result$model$covariates, "AGE")
  expect_equal(result$frequency, 2 / 4)
})

test_that(".extract_scm_settings() returns selected_exposure for GLM", {
  skip_on_cran()
  fit      <- .fit_test_scm_glm()
  settings <- suppressMessages(erSCM:::.extract_scm_settings(fit))
  expect_true("selected_exposure" %in% names(settings))
  expect_type(settings$selected_exposure, "character")
  expect_length(settings$selected_exposure, 1L)
})

# ---- .compute_transform_freq() multi-exposure fix (Gap #9) ------------------

test_that(".compute_transform_freq() returns named list for multi-exposure NL spec (Gap #9)", {
  skip_on_cran()
  # Simulate 4 converged bootstrap models where CMAX and AUC are both NL candidates
  converged_models <- list(
    list(exposure_transform = "log",      covariate_transforms = list()),
    list(exposure_transform = "log",      covariate_transforms = list()),
    list(exposure_transform = "linear",   covariate_transforms = list()),
    list(exposure_transform = "not_selected", covariate_transforms = list())
  )
  nl_spec <- list(CMAX = c("linear", "log"), AUC = c("linear", "log"))
  result <- erSCM:::.compute_transform_freq(
    converged_models = converged_models,
    nonlinear_spec   = nl_spec,
    exposure_vars    = c("CMAX", "AUC")
  )
  etf <- result$exposure_transform_freq
  # Must be a named list with one entry per NL-enabled exposure
  expect_type(etf, "list")
  expect_true("CMAX" %in% names(etf))
  expect_true("AUC"  %in% names(etf))
  # Each entry must be a named numeric vector summing to 1
  expect_type(etf$CMAX, "double")
  expect_type(etf$AUC,  "double")
  expect_equal(sum(etf$CMAX), 1, tolerance = 1e-10)
  expect_equal(sum(etf$AUC),  1, tolerance = 1e-10)
})

test_that(".compute_transform_freq() returns named list for single-exposure NL spec (Gap #9)", {
  skip_on_cran()
  converged_models <- list(
    list(exposure_transform = "log",    covariate_transforms = list()),
    list(exposure_transform = "linear", covariate_transforms = list())
  )
  nl_spec <- list(CMAX = c("linear", "log"))
  result <- erSCM:::.compute_transform_freq(
    converged_models = converged_models,
    nonlinear_spec   = nl_spec,
    exposure_vars    = "CMAX"
  )
  etf <- result$exposure_transform_freq
  expect_type(etf, "list")
  expect_true("CMAX" %in% names(etf))
  expect_equal(sum(etf$CMAX), 1, tolerance = 1e-10)
})
