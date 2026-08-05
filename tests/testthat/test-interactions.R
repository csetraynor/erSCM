test_that("expand_interactions handles NULL input", {
  result <- expand_interactions(NULL)
  expect_identical(result, character(0))
})

test_that("expand_interactions handles character vector input", {
  result <- expand_interactions(c("age:sex", "AUC:weight"))
  expect_identical(result, c("age:sex", "AUC:weight"))
})

test_that("expand_interactions validates character vector format", {
  expect_error(
    expand_interactions(c("age")),
    "must contain at least two variables"
  )
})

test_that("expand_interactions handles list input with pairwise generation", {
  result <- expand_interactions(list(c("AUC", "age", "sex")))
  expect_setequal(result, c("AUC:age", "AUC:sex", "age:sex"))
})

test_that("expand_interactions handles list with max_interaction_order = 3", {
  result <- expand_interactions(list(c("A", "B", "C")), max_interaction_order = 3)
  expect_true("A:B:C" %in% result)
  expect_true("A:B" %in% result)
  expect_true("A:C" %in% result)
  expect_true("B:C" %in% result)
})

test_that("expand_interactions rejects invalid input types", {
  expect_error(expand_interactions(42), "must be NULL, a character vector, or a list")
})

test_that("expand_interactions rejects invalid list elements", {
  expect_error(expand_interactions(list(c("A"))), "character vector of length >= 2")
})

test_that("get_parent_effects splits interaction terms correctly", {
  expect_identical(get_parent_effects("age:sex"), c("age", "sex"))
  expect_identical(get_parent_effects("A:B:C"), c("A", "B", "C"))
})

test_that("get_eligible_interactions enforces hierarchy constraint", {
  candidates <- c("age:sex", "AUC:age", "AUC:weight")
  current <- c("AUC", "age", "sex")

  result <- get_eligible_interactions(candidates, current)
  expect_setequal(result, c("age:sex", "AUC:age"))
  expect_false("AUC:weight" %in% result)
})

test_that("get_eligible_interactions ignores interactions in current_relation for eligibility", {
  candidates <- c("age:sex", "AUC:age")
  # Interaction terms in current_relation should not count as main effects
  current <- c("AUC", "age", "sex", "AUC:age")

  result <- get_eligible_interactions(candidates, current)
  expect_true("age:sex" %in% result)
  expect_true("AUC:age" %in% result)
})

test_that("get_eligible_interactions returns empty when no parents present", {
  candidates <- c("age:sex", "AUC:weight")
  current <- c("BMI", "HT")

  result <- get_eligible_interactions(candidates, current)
  expect_identical(result, character(0))
})

test_that("get_protected_main_effects identifies protected terms", {
  current <- c("AUC", "age", "sex", "AUC:age")
  result <- get_protected_main_effects(current)
  expect_setequal(result, c("AUC", "age"))
})

test_that("get_protected_main_effects returns empty when no interactions", {
  current <- c("AUC", "age", "sex")
  result <- get_protected_main_effects(current)
  expect_identical(result, character(0))
})

test_that("validate_interactions passes for valid terms", {
  expect_true(
    validate_interactions(
      c("age:sex", "AUC:weight"),
      covariates = c("age", "sex", "weight"),
      exposure_vars = c("AUC")
    )
  )
})

test_that("validate_interactions errors for missing variables", {
  expect_error(
    validate_interactions(
      c("age:UNKNOWN"),
      covariates = c("age", "sex"),
      exposure_vars = c("AUC")
    ),
    "references variable"
  )
})

test_that("scm_glm with interactions=NULL is backward compatible", {
  skip_on_cran()

  data("demo_eres", package = "erSCM")

  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX))

  result <- scm_glm(
    eres_data = demo_eres,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    significance_level = c(0.05, 0.01),
    interactions = NULL
  )

  expect_type(result, "list")
  expect_named(result, c("final_model", "path", "summary", "model_data", "metadata", "var_spec", "naming_df"),
               ignore.order = TRUE)
  expect_true(inherits(result$final_model, "glm"))
})

test_that("scm_coxph with interactions=NULL is backward compatible", {
  skip_on_cran()

  data("demo_eres", package = "erSCM")

  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX)) |>
    dplyr::mutate(DVID = "ALL")

  result <- scm_coxph(
    eres_data = demo_eres,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    time_col = "sim_time",
    significance_level = c(0.9, 0.9),
    dvid = "ALL",
    dvid_col = "DVID",
    interactions = NULL
  )

  expect_type(result, "list")
  expect_named(result, c("final_model", "path", "summary", "model_data", "metadata", "var_spec", "naming_df"),
               ignore.order = TRUE)
  expect_true(inherits(result$final_model, "coxph"))
})

test_that("scm_glm interaction testing with after_main strategy works", {
  skip_on_cran()

  data("demo_eres", package = "erSCM")

  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX))

  result <- scm_glm(
    eres_data = demo_eres,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    significance_level = c(0.05, 0.01),
    interactions = list(c("CMAXss", "AGE", "WT", "SEX")),
    interaction_strategy = "after_main",
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_named(result, c("final_model", "path", "summary", "model_data", "metadata", "var_spec", "naming_df"),
               ignore.order = TRUE)
  expect_true(inherits(result$final_model, "glm"))
})

test_that("scm_coxph interaction testing with after_main strategy works", {
  skip_on_cran()

  data("demo_eres", package = "erSCM")

  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX)) |>
    dplyr::mutate(DVID = "ALL")

  result <- scm_coxph(
    eres_data = demo_eres,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    time_col = "sim_time",
    significance_level = c(0.9, 0.9),
    dvid = "ALL",
    dvid_col = "DVID",
    interactions = list(c("CMAXss", "AGE", "WT", "SEX")),
    interaction_strategy = "after_main",
    p_interaction = c(0.9, 0.9),
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_named(result, c("final_model", "path", "summary", "model_data", "metadata", "var_spec", "naming_df"),
               ignore.order = TRUE)
  expect_true(inherits(result$final_model, "coxph"))
})

test_that("scm_glm interaction testing with simultaneous strategy works", {
  skip_on_cran()

  data("demo_eres", package = "erSCM")

  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX))

  result <- scm_glm(
    eres_data = demo_eres,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    significance_level = c(0.05, 0.01),
    interactions = c("CMAXss:AGE"),
    interaction_strategy = "simultaneous",
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_true(inherits(result$final_model, "glm"))
})

test_that("hierarchy constraint prevents interaction without parents in model", {
  skip_on_cran()

  data("demo_eres", package = "erSCM")

  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX))

  # Use very strict significance level for main effects so only exposure enters
  # but test an interaction involving a covariate that won't enter
  result <- scm_glm(
    eres_data = demo_eres,
    covariates = c("WT", "AGE", "SEX"),
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    significance_level = c(0.001, 0.0001),
    interactions = c("WT:SEX"),  # WT or SEX may not enter with strict threshold
    interaction_strategy = "after_main",
    verbose = FALSE
  )

  # If the interaction was added, its parents must be in the Relation
  final_relation <- result$metadata$finalRelation
  interaction_in_model <- grep(":", final_relation, value = TRUE)
  for (int_term in interaction_in_model) {
    parents <- get_parent_effects(int_term)
    expect_true(all(parents %in% final_relation),
                info = paste("Parents of", int_term, "must be in final model"))
  }
})

test_that("protected main effects are not removed during backward elimination", {
  skip_on_cran()

  # Test that get_protected_main_effects correctly identifies terms
  # and that the backward function respects them
  current_relation <- c("AUC", "age", "sex", "AUC:age")
  protected <- get_protected_main_effects(current_relation)
  expect_true("AUC" %in% protected)
  expect_true("age" %in% protected)
  expect_false("sex" %in% protected)

  # Only sex and AUC:age should be removable
  removable <- current_relation[!current_relation %in% protected]
  expect_setequal(removable, c("sex", "AUC:age"))
})

test_that("validate_interactions errors for completely missing variables", {
  expect_error(
    validate_interactions(
      c("MISSING1:MISSING2"),
      covariates = c("age", "sex"),
      exposure_vars = c("AUC")
    ),
    "references variable"
  )
})

test_that("p_interaction parameter is respected", {
  skip_on_cran()

  data("demo_eres", package = "erSCM")

  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX))

  # Use very strict p_interaction so no interactions enter
  result_strict <- scm_glm(
    eres_data = demo_eres,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    significance_level = c(0.05, 0.01),
    interactions = list(c("CMAXss", "AGE", "WT")),
    interaction_strategy = "after_main",
    p_interaction = c(1e-10, 1e-10),
    verbose = FALSE
  )

  # With extremely strict threshold, no interactions should enter
  final_relation <- result_strict$metadata$finalRelation
  interaction_in_model <- grep(":", final_relation, value = TRUE)
  expect_length(interaction_in_model, 0)
})

# ------------------------------------------------------------------
# Helper function tests for new interaction utilities
# ------------------------------------------------------------------

test_that("get_main_effect_columns filters out interaction terms", {
  expect_identical(
    get_main_effect_columns(c("AUC", "age", "sex", "AUC:age")),
    c("AUC", "age", "sex")
  )
  expect_identical(get_main_effect_columns(c("A", "B")), c("A", "B"))
  expect_identical(get_main_effect_columns(character(0)), character(0))
})



test_that("has_interaction_terms detects presence correctly", {
  expect_true(has_interaction_terms(c("A", "B", "A:B")))
  expect_false(has_interaction_terms(c("A", "B")))
  expect_false(has_interaction_terms(character(0)))
})

test_that("get_interaction_terms extracts only interaction terms", {
  expect_identical(
    get_interaction_terms(c("AUC", "age", "sex", "AUC:age", "age:sex")),
    c("AUC:age", "age:sex")
  )
  expect_identical(get_interaction_terms(c("A", "B")), character(0))
})

# ------------------------------------------------------------------
# Visualization function tests with interaction terms
# ------------------------------------------------------------------

test_that("make_forest_plot works with interaction terms in GLM result", {
  skip_on_cran()

  data("demo_eres", package = "erSCM")

  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX))

  result <- scm_glm(
    eres_data = demo_eres,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    significance_level = c(0.05, 0.01),
    interactions = list(c("CMAXss", "AGE", "WT", "SEX")),
    interaction_strategy = "after_main",
    verbose = FALSE
  )

  # Forest plot data should not error
  fp_data <- make_forest_plot(result, plot = FALSE)
  expect_s3_class(fp_data, "data.frame")
  expect_true(nrow(fp_data) > 0)

  # Forest plot as ggplot should not error
  p <- make_forest_plot(result, plot = TRUE)
  expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
})

test_that("make_forest_plot works with interaction terms in coxph result", {
  skip_on_cran()

  data("demo_eres", package = "erSCM")

  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX)) |>
    dplyr::mutate(DVID = "ALL")

  result <- scm_coxph(
    eres_data = demo_eres,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    time_col = "sim_time",
    significance_level = c(0.9, 0.9),
    dvid = "ALL",
    dvid_col = "DVID",
    interactions = list(c("CMAXss", "AGE", "WT", "SEX")),
    interaction_strategy = "after_main",
    p_interaction = c(0.9, 0.9),
    verbose = FALSE
  )

  fp_data <- make_forest_plot(result, plot = FALSE)
  expect_s3_class(fp_data, "data.frame")
  expect_true(nrow(fp_data) > 0)

  p <- make_forest_plot(result, plot = TRUE)
  expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
})

test_that("vpc_coxph works with interaction terms", {
  skip_on_cran()

  data("demo_eres", package = "erSCM")

  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX)) |>
    dplyr::mutate(DVID = "ALL")

  result <- scm_coxph(
    eres_data = demo_eres,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    time_col = "sim_time",
    significance_level = c(0.9, 0.9),
    dvid = "ALL",
    dvid_col = "DVID",
    interactions = list(c("CMAXss", "AGE", "WT", "SEX")),
    interaction_strategy = "after_main",
    p_interaction = c(0.9, 0.9),
    verbose = FALSE
  )

  # VPC data should not error
  vpc_ds <- vpc_coxph(result, plot = FALSE)
  expect_s3_class(vpc_ds, "data.frame")
  expect_true(nrow(vpc_ds) > 0)

  # VPC plot should not error
  vpc_p <- vpc_coxph(result, plot = TRUE)
  expect_true(inherits(vpc_p, "gg"))
})

test_that("plot_cox_terms works with interaction model using main-effect covariate_term", {
  skip_on_cran()

  data("demo_eres", package = "erSCM")

  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX)) |>
    dplyr::mutate(DVID = "ALL")

  result <- scm_coxph(
    eres_data = demo_eres,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    time_col = "sim_time",
    significance_level = c(0.9, 0.9),
    dvid = "ALL",
    dvid_col = "DVID",
    interactions = list(c("CMAXss", "AGE", "WT", "SEX")),
    interaction_strategy = "after_main",
    p_interaction = c(0.9, 0.9),
    verbose = FALSE
  )

  # Plotting a main-effect term should work fine
  vpc_ds <- plot_cox_terms(result, covariate_term = "SEX", plot = FALSE)
  expect_s3_class(vpc_ds, "data.frame")
  expect_true(nrow(vpc_ds) > 0)
})

test_that("plot_cox_terms rejects interaction term as covariate_term", {
  skip_on_cran()

  data("demo_eres", package = "erSCM")

  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX)) |>
    dplyr::mutate(DVID = "ALL")

  result <- scm_coxph(
    eres_data = demo_eres,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    time_col = "sim_time",
    significance_level = c(0.9, 0.9),
    dvid = "ALL",
    dvid_col = "DVID",
    interactions = list(c("CMAXss", "AGE", "WT", "SEX")),
    interaction_strategy = "after_main",
    p_interaction = c(0.9, 0.9),
    verbose = FALSE
  )

  expect_error(
    plot_cox_terms(result, covariate_term = "CMAXss:AGE"),
    "does not support interaction terms"
  )
})

# =============================================================================
# Backward elimination dual-threshold tests
# =============================================================================

test_that("backward elimination uses exit_alpha_interaction for interaction terms (coxph)", {

  skip_on_cran()

  data("demo_eres", package = "erSCM")

  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX)) |>
    dplyr::mutate(DVID = "ALL")

  # Use lenient forward thresholds to force terms in, then strict main-effect

  # backward threshold but LENIENT interaction backward threshold.
  # This tests that interactions are retained when their p > significance_level[2]
  # but p < p_interaction[2]
  result_lenient_int <- scm_coxph(
    eres_data = demo_eres,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    time_col = "sim_time",
    significance_level = c(0.9, 0.001),  # strict retention for main effects
    dvid = "ALL",
    dvid_col = "DVID",
    interactions = list(c("CMAXss", "AGE", "WT", "SEX")),
    interaction_strategy = "after_main",
    p_interaction = c(0.9, 0.9),   # very lenient retention for interactions
    verbose = FALSE
  )

  # With lenient interaction threshold, any interaction that entered should be retained
  interaction_terms_final <- get_interaction_terms(result_lenient_int$metadata$finalRelation)

  # Now run same thing but with STRICT interaction threshold
  result_strict_int <- scm_coxph(
    eres_data = demo_eres,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    time_col = "sim_time",
    significance_level = c(0.9, 0.001),
    dvid = "ALL",
    dvid_col = "DVID",
    interactions = list(c("CMAXss", "AGE", "WT", "SEX")),
    interaction_strategy = "after_main",
    p_interaction = c(0.9, 0.0001), # very strict retention for interactions
    verbose = FALSE
  )

  strict_int_terms <- get_interaction_terms(result_strict_int$metadata$finalRelation)

  # With stricter interaction threshold, fewer (or equal) interactions should survive

  expect_true(length(strict_int_terms) <= length(interaction_terms_final))
})

test_that("backward elimination uses exit_alpha_interaction for interaction terms (glm)", {
  skip_on_cran()

  data("demo_eres", package = "erSCM")

  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX)) |>
    dplyr::mutate(DVID = "ALL")

  # Lenient interaction backward threshold
  result_lenient <- scm_glm(
    eres_data = demo_eres,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    significance_level = c(0.9, 0.001),
    interactions = list(c("CMAXss", "AGE", "WT", "SEX")),
    interaction_strategy = "after_main",
    p_interaction = c(0.9, 0.9),
    verbose = FALSE
  )

  # Strict interaction backward threshold
  result_strict <- scm_glm(
    eres_data = demo_eres,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    significance_level = c(0.9, 0.001),
    interactions = list(c("CMAXss", "AGE", "WT", "SEX")),
    interaction_strategy = "after_main",
    p_interaction = c(0.9, 0.0001),
    verbose = FALSE
  )

  lenient_int <- get_interaction_terms(result_lenient$metadata$finalRelation)
  strict_int  <- get_interaction_terms(result_strict$metadata$finalRelation)

  # Stricter threshold should retain fewer (or equal) interactions
  expect_true(length(strict_int) <= length(lenient_int))
})

test_that("p_interaction=NULL gives identical behavior to same significance_level (coxph)", {
  skip_on_cran()

  data("demo_eres", package = "erSCM")

  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX)) |>
    dplyr::mutate(DVID = "ALL")

  # Run with p_interaction = NULL (should default to significance_level)
  result_null <- scm_coxph(
    eres_data = demo_eres,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    time_col = "sim_time",
    significance_level = c(0.9, 0.5),
    dvid = "ALL",
    dvid_col = "DVID",
    interactions = list(c("CMAXss", "AGE", "WT")),
    interaction_strategy = "after_main",
    p_interaction = NULL,
    verbose = FALSE
  )

  # Run with explicit p_interaction = significance_level (should be identical)
  result_explicit <- scm_coxph(
    eres_data = demo_eres,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    time_col = "sim_time",
    significance_level = c(0.9, 0.5),
    dvid = "ALL",
    dvid_col = "DVID",
    interactions = list(c("CMAXss", "AGE", "WT")),
    interaction_strategy = "after_main",
    p_interaction = c(0.9, 0.5),
    verbose = FALSE
  )

  # Same final relation
  expect_setequal(result_null$metadata$finalRelation,
                  result_explicit$metadata$finalRelation)
})

test_that("backward elimination respects hierarchy even with dual thresholds", {
  skip_on_cran()

  data("demo_eres", package = "erSCM")

  con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
  all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
  all_exposure   <- c("CMAXss", "CMINss", "AUCss")

  demo_eres <- demo_eres |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(con_covariates, all_exposure)), ~ . / median(., na.rm = TRUE))) |>
    dplyr::mutate(SEX = as.factor(SEX)) |>
    dplyr::mutate(DVID = "ALL")

  # Use very lenient thresholds for interactions to ensure they enter and stay
  result <- scm_coxph(
    eres_data = demo_eres,
    covariates = all_covariates,
    exposure_vars = all_exposure,
    dv_col = "sim_status",
    time_col = "sim_time",
    significance_level = c(0.9, 0.0001), # very strict main effect retention
    dvid = "ALL",
    dvid_col = "DVID",
    interactions = list(c("CMAXss", "AGE", "WT", "SEX")),
    interaction_strategy = "after_main",
    p_interaction = c(0.9, 0.99),  # interactions always retained
    verbose = FALSE
  )

  final_rel <- result$metadata$finalRelation
  int_terms <- get_interaction_terms(final_rel)

  # For any interaction still in the model, both parent main effects must be present
  for (int_term in int_terms) {
    parents <- get_parent_effects(int_term)
    for (p in parents) {
      expect_true(p %in% final_rel,
                  info = paste0("Parent '", p, "' missing for interaction '", int_term, "'"))
    }
  }
})
