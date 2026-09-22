## Tests for the var_spec system: parse_var_spec, validate_var_spec, integration

# ---- Parse: named vector ----

test_that("named vector with units in parentheses extracts label and units", {
  spec <- c(WT = "Body Weight (kg)", AGE = "Age (years)", CMAXss = "Cmax,ss (ng/mL)")
  result <- parse_var_spec(spec)
  expect_equal(result$col_name, c("WT", "AGE", "CMAXss"))
  expect_equal(result$label,    c("Body Weight", "Age", "Cmax,ss"))
  expect_equal(result$units,    c("kg", "years", "ng/mL"))
})

test_that("named vector without units sets units to NA", {
  spec <- c(SEX = "Sex", RACE = "Race")
  result <- parse_var_spec(spec)
  expect_equal(result$label, c("Sex", "Race"))
  expect_true(all(is.na(result$units)))
})

test_that("data.frame with all columns produces correct canonical tibble", {
  df <- data.frame(
    col_name = c("WT", "AGE", "SEX"),
    label    = c("Body Weight", "Age", "Sex"),
    units    = c("kg", "years", NA),
    stringsAsFactors = FALSE
  )
  result <- parse_var_spec(df)
  expect_equal(result$col_name, c("WT", "AGE", "SEX"))
  expect_equal(result$label,    c("Body Weight", "Age", "Sex"))
  expect_equal(result$units[1], "kg")
  expect_true(is.na(result$units[3]))
  expect_equal(result$transform, rep("none", 3))
})

test_that("data.frame with only required columns sets correct defaults", {
  df <- data.frame(col_name = c("WT", "AGE"), label = c("Weight", "Age"),
                   stringsAsFactors = FALSE)
  result <- parse_var_spec(df)
  expect_true(all(is.na(result$units)))
  expect_true(all(is.na(result$var_type)))
  expect_equal(result$transform, rep("none", 2))
  expect_true(is.list(result$categories))
})

test_that("data.frame with pipe-separated categories parses them", {
  df <- data.frame(
    col_name   = c("SEX", "RACE"),
    label      = c("Sex", "Race"),
    categories = c("Female|Male", "White|Black|Asian"),
    stringsAsFactors = FALSE
  )
  result <- parse_var_spec(df)
  expect_equal(result$categories[[1]], c("Female", "Male"))
  expect_equal(result$categories[[2]], c("White", "Black", "Asian"))
})

test_that("YAML path parses continuous and categorical variables", {
  # Build a temp YAML file
  yaml_content <- "
WT:
  label: Body Weight
  units: kg
  type: continuous
SEX:
  label: Sex
  type: categorical
  categories:
    F: Female
    M: Male
"
  tmp <- tempfile(fileext = ".yaml")
  writeLines(yaml_content, tmp)
  on.exit(unlink(tmp))

  skip_if_not_installed("yaml")
  result <- parse_var_spec(tmp)
  expect_equal(result$col_name, c("WT", "SEX"))
  expect_equal(result$label,    c("Body Weight", "Sex"))
  expect_equal(result$units[1], "kg")
  expect_equal(result$var_type, c("continuous", "categorical"))
})

test_that("invalid input (unnamed numeric) gives informative error", {
  expect_error(parse_var_spec(42),
               regexp = "data.frame, named character vector, or path")
})

test_that("invalid input (unnamed character vector) gives informative error", {
  expect_error(parse_var_spec(c("WT", "AGE")),
               regexp = "data.frame, named character vector, or path")
})

test_that("non-existent YAML path gives informative error", {
  expect_error(parse_var_spec("/no/such/file.yaml"),
               regexp = "does not exist")
})

# ---- Validate ----

test_that("all col_names present in data produces no warning", {
  df <- data.frame(WT = 1:5, AGE = 11:15)
  spec <- parse_var_spec(c(WT = "Weight (kg)", AGE = "Age (years)"))
  expect_no_warning(erSCM:::validate_var_spec(spec, df, "WT", "AGE"))
})

test_that("missing col_names in data produce warning and are filtered out", {
  df <- data.frame(WT = 1:5)
  spec <- parse_var_spec(c(WT = "Weight (kg)", MISSING = "Ghost"))
  result <- expect_warning(
    erSCM:::validate_var_spec(spec, df, "WT", character(0)),
    regexp = "not found in data"
  )
  expect_false("MISSING" %in% result$col_name)
})

test_that("model variables not in spec trigger message and are auto-added", {
  df <- data.frame(WT = 1:5, AGE = 11:15)
  spec <- parse_var_spec(c(WT = "Weight (kg)"))
  msgs <- character(0)
  result <- withCallingHandlers(
    erSCM:::validate_var_spec(spec, df, "WT", "AGE"),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(any(grepl("not in var_spec", msgs)),
              info = paste("Messages:", paste(msgs, collapse = "; ")))
  expect_true("AGE" %in% result$col_name)
})

test_that("var_type auto-detection: numeric -> continuous, factor -> categorical", {
  df <- data.frame(WT = 1:5, SEX = factor(c("M", "F", "M", "F", "M")))
  spec <- tibble::tibble(
    col_name = c("WT", "SEX"), label = c("Weight", "Sex"),
    units = NA_character_, var_type = NA_character_,
    categories = list(NULL, NULL), transform = "none"
  )
  result <- erSCM:::validate_var_spec(spec, df, c("WT", "SEX"), character(0))
  expect_equal(result$var_type[result$col_name == "WT"],  "continuous")
  expect_equal(result$var_type[result$col_name == "SEX"], "categorical")
})

# ---- Integration: scm_glm ----

test_that("scm_glm with var_spec named vector stores canonical tibble in fit$var_spec", {
  skip_on_cran()
  set.seed(1)
  n <- 200
  df <- data.frame(DV = rbinom(n, 1, 0.3), EXP = rnorm(n, 10, 2), AGE = rnorm(n, 50, 10))
  spec <- c(EXP = "Exposure (ng/mL)", AGE = "Age (years)")
  result <- scm_glm(eres_data = df, covariates = "AGE", exposure_vars = "EXP",
                    dv_col = "DV", significance_level = c(0.9, 0.8),
                    var_spec = spec, center = FALSE)
  expect_false(is.null(result$var_spec))
  expect_true(is.data.frame(result$var_spec))
  expect_true("col_name" %in% colnames(result$var_spec))
  expect_true("EXP" %in% result$var_spec$col_name)
  expect_equal(result$var_spec$label[result$var_spec$col_name == "EXP"], "Exposure")
  expect_equal(result$var_spec$units[result$var_spec$col_name == "EXP"], "ng/mL")
})

test_that("scm_glm with var_spec data.frame stores canonical tibble", {
  skip_on_cran()
  set.seed(1)
  n <- 200
  df <- data.frame(DV = rbinom(n, 1, 0.3), EXP = rnorm(n, 10, 2))
  spec_df <- data.frame(col_name = "EXP", label = "Exposure", units = "ng/mL",
                        stringsAsFactors = FALSE)
  result <- scm_glm(eres_data = df, covariates = character(0), exposure_vars = "EXP",
                    dv_col = "DV", significance_level = c(0.9, 0.8),
                    var_spec = spec_df, center = FALSE)
  expect_false(is.null(result$var_spec))
  expect_equal(result$var_spec$units[result$var_spec$col_name == "EXP"], "ng/mL")
})

test_that("scm_glm with var_spec = NULL stores NULL in fit$var_spec (backward compat)", {
  skip_on_cran()
  set.seed(3)
  n <- 200
  df <- data.frame(DV = rbinom(n, 1, 0.3), EXP = rnorm(n, 10, 2))
  result <- scm_glm(eres_data = df, covariates = character(0), exposure_vars = "EXP",
                    dv_col = "DV", significance_level = c(0.9, 0.8))
  expect_null(result$var_spec)
})

test_that("scm_coxph with var_spec stores canonical tibble", {
  skip_on_cran()
  set.seed(4)
  n <- 200
  df <- data.frame(TIME = rexp(n, 0.1), DV = rbinom(n, 1, 0.3),
                   EXP = rnorm(n, 10, 2), WT = rnorm(n, 70, 10))
  spec <- c(EXP = "Exposure (ng/mL)", WT = "Body Weight (kg)")
  result <- scm_coxph(eres_data = df, covariates = "WT", exposure_vars = "EXP",
                      dv_col = "DV", time_col = "TIME",
                      significance_level = c(0.9, 0.8),
                      var_spec = spec, center = FALSE)
  expect_false(is.null(result$var_spec))
  expect_equal(result$var_spec$label[result$var_spec$col_name == "EXP"], "Exposure")
})

test_that("forest plot labels use display names from var_spec", {
  skip_on_cran()
  set.seed(5)
  n <- 300
  df <- data.frame(DV = rbinom(n, 1, 0.3), EXP = rnorm(n, 100, 20))
  spec <- c(EXP = "Exposure (ng/mL)")
  result <- scm_glm(eres_data = df, covariates = character(0), exposure_vars = "EXP",
                    dv_col = "DV", significance_level = c(0.9, 0.8),
                    var_spec = spec, center = FALSE)
  skip_if(!"EXP" %in% result$metadata$finalRelation, "EXP not selected")
  fp <- make_forest_plot(result, plot = FALSE)
  # group labels should use the display name "Exposure", not raw "EXP"
  expect_true("Exposure" %in% fp$group,
              info = paste("Groups found:", paste(unique(fp$group), collapse = ", ")))
})

test_that("forest plot labels show correct units from var_spec", {
  skip_on_cran()
  set.seed(6)
  n <- 300
  df <- data.frame(DV = rbinom(n, 1, 0.3), EXP = rnorm(n, 100, 20))
  spec <- c(EXP = "Exposure (ng/mL)")
  result <- scm_glm(eres_data = df, covariates = character(0), exposure_vars = "EXP",
                    dv_col = "DV", significance_level = c(0.9, 0.8),
                    var_spec = spec, center = FALSE)
  skip_if(!"EXP" %in% result$metadata$finalRelation, "EXP not selected")
  fp <- make_forest_plot(result, plot = FALSE)
  exp_rows <- fp[fp$group == "Exposure", ]
  skip_if(nrow(exp_rows) == 0, "Exposure not in forest plot")
  # group_level should contain "ng/mL" for numeric rows
  numeric_levels <- exp_rows$group_level[grepl("[0-9]", exp_rows$group_level)]
  if (length(numeric_levels) > 0) {
    expect_true(any(grepl("ng/mL", numeric_levels)),
                info = paste("group_level values:", paste(numeric_levels, collapse = "; ")))
  }
})

test_that("interaction labels use var_spec display names", {
  skip_on_cran()
  set.seed(7)
  n <- 500
  WT_true <- rnorm(n, 70, 15)
  EXP_true <- rnorm(n, 100, 20)
  log_odds <- -2 + 0.02 * EXP_true + 0.001 * WT_true + 0.0003 * EXP_true * WT_true
  df <- data.frame(DV = rbinom(n, 1, plogis(log_odds)), EXP = EXP_true, WT = WT_true)
  spec <- c(EXP = "Exposure (ng/mL)", WT = "Body Weight (kg)")
  result <- scm_glm(eres_data = df, covariates = "WT", exposure_vars = "EXP",
                    dv_col = "DV", significance_level = c(0.9, 0.8),
                    interactions = "EXP:WT", center = TRUE, var_spec = spec)
  skip_if(!"EXP:WT" %in% result$metadata$finalRelation, "EXP:WT not selected")
  fp <- make_forest_plot(result, plot = FALSE)
  # Interaction group should use display labels
  int_groups <- fp$group[grepl("×", fp$group)]
  expect_true(length(int_groups) > 0)
  expect_true(any(grepl("Exposure|Body Weight", int_groups)),
              info = paste("Interaction groups:", paste(unique(int_groups), collapse = ", ")))
})

test_that("print method uses var_spec display labels", {
  skip_on_cran()
  set.seed(8)
  n <- 200
  df <- data.frame(DV = rbinom(n, 1, 0.3), EXP = rnorm(n, 10, 2))
  spec <- c(EXP = "Exposure (ng/mL)")
  result <- scm_glm(eres_data = df, covariates = character(0), exposure_vars = "EXP",
                    dv_col = "DV", significance_level = c(0.9, 0.8),
                    var_spec = spec, center = FALSE)
  skip_if(!"EXP" %in% result$metadata$finalRelation, "EXP not selected")
  out <- capture.output(print(result))
  # Should mention "Exposure" not raw "EXP" in the covariates section
  expect_true(any(grepl("Exposure", out)),
              info = paste("Print output:", paste(out, collapse = "\n")))
})

# ---- Backward compatibility ----

test_that("scm_glm with normalize_inputs=TRUE emits deprecation warning", {
  skip_on_cran()
  set.seed(9)
  n <- 200
  df <- data.frame(DV = rbinom(n, 1, 0.3), EXP = rnorm(n, 10, 2))
  expect_warning(
    scm_glm(eres_data = df, covariates = character(0), exposure_vars = "EXP",
            dv_col = "DV", significance_level = c(0.9, 0.8),
            normalize_inputs = TRUE, center = FALSE),
    regexp = "deprecated"
  )
})

test_that("scm_coxph with normalize_inputs=TRUE emits deprecation warning", {
  skip_on_cran()
  set.seed(10)
  n <- 200
  df <- data.frame(TIME = rexp(n, 0.1), DV = rbinom(n, 1, 0.3), EXP = rnorm(n, 10, 2))
  expect_warning(
    scm_coxph(eres_data = df, covariates = character(0), exposure_vars = "EXP",
              dv_col = "DV", time_col = "TIME", significance_level = c(0.9, 0.8),
              normalize_inputs = TRUE, center = FALSE),
    regexp = "deprecated"
  )
})

test_that("make_forest_plot with explicit namings overrides fit$var_spec", {
  skip_on_cran()
  set.seed(11)
  n <- 300
  df <- data.frame(DV = rbinom(n, 1, 0.3), EXP = rnorm(n, 100, 20))
  spec <- c(EXP = "Exposure (ng/mL)")
  result <- scm_glm(eres_data = df, covariates = character(0), exposure_vars = "EXP",
                    dv_col = "DV", significance_level = c(0.9, 0.8),
                    var_spec = spec, center = FALSE)
  skip_if(!"EXP" %in% result$metadata$finalRelation, "EXP not selected")

  # Pass an old-style namings override with different units
  old_namings <- data.frame(covname = "EXP", median_point = median(df$EXP),
                             units = "override_units", stringsAsFactors = FALSE)
  fp <- make_forest_plot(result, plot = FALSE, namings = old_namings)
  # The override units should appear (not ng/mL from var_spec)
  numeric_levels <- fp$group_level[grepl("[0-9]", fp$group_level)]
  if (length(numeric_levels) > 0) {
    expect_true(any(grepl("override_units", numeric_levels)),
                info = paste("Levels:", paste(numeric_levels[1:min(3, length(numeric_levels))],
                                              collapse = "; ")))
  }
})

test_that("old fit objects without var_spec fall back to naming_df gracefully", {
  skip_on_cran()
  set.seed(12)
  n <- 300
  df <- data.frame(DV = rbinom(n, 1, 0.3), EXP = rnorm(n, 100, 20))
  result <- suppressWarnings(
    scm_glm(eres_data = df, covariates = character(0), exposure_vars = "EXP",
            dv_col = "DV", significance_level = c(0.9, 0.8),
            normalize_inputs = TRUE, center = FALSE)
  )
  skip_if(!"EXP" %in% result$metadata$finalRelation, "EXP not selected")

  # Simulate old object: remove var_spec
  result_old <- result
  result_old$var_spec <- NULL

  fp <- make_forest_plot(result_old, plot = FALSE)
  expect_s3_class(fp, "data.frame")
  expect_true(nrow(fp) > 0)
})

# ---- Scale-to-median + center combination ----

test_that("scale_to_median=TRUE and center=TRUE together issue a warning", {
  skip_on_cran()
  set.seed(13)
  n <- 200
  df <- data.frame(DV = rbinom(n, 1, 0.3), EXP = rnorm(n, 10, 2))
  expect_warning(
    scm_glm(eres_data = df, covariates = character(0), exposure_vars = "EXP",
            dv_col = "DV", significance_level = c(0.9, 0.8),
            scale_to_median = TRUE, center = TRUE),
    regexp = "scale_to_median"
  )
})

test_that("apply_var_labels covname stripping is element-wise (vectorised label regression)", {
  # Regression: gsub(label, "", covname, fixed=TRUE) inside mutate uses label[1] only.
  # With mapply the replacement must be per-row, not first-element.
  fp <- data.frame(
    covname = c("SEXFemale", "RACEWhite", "RACEBlack"),
    label   = c("Female",    "White",     "Black"),
    VAL     = c(NA_real_,    NA_real_,    NA_real_),
    stringsAsFactors = FALSE
  )
  out <- erSCM:::apply_var_labels(fp, var_meta = NULL)
  expect_equal(out$covname, c("SEX", "RACE", "RACE"))
})

test_that("var_spec transform = 'centered' correctly back-transforms values in forest plot", {
  skip_on_cran()
  set.seed(14)
  n <- 300
  EXP_vals <- rnorm(n, 100, 20)
  df <- data.frame(DV = rbinom(n, 1, 0.3), EXP = EXP_vals)
  spec <- c(EXP = "Exposure (ng/mL)")
  result <- scm_glm(eres_data = df, covariates = character(0), exposure_vars = "EXP",
                    dv_col = "DV", significance_level = c(0.9, 0.8),
                    var_spec = spec, center = TRUE)
  skip_if(!"EXP" %in% result$metadata$finalRelation, "EXP not selected")

  fp <- make_forest_plot(result, plot = FALSE)
  exp_rows <- fp[fp$group == "Exposure", ]
  skip_if(nrow(exp_rows) == 0, "Exposure not in forest plot")

  # Back-transformed quantile labels should be in original scale (~100 ± 3*20)
  num_labels <- suppressWarnings(as.numeric(gsub("[^0-9.-]", "", exp_rows$group_level)))
  non_na <- num_labels[!is.na(num_labels) & num_labels != 0]
  skip_if(length(non_na) == 0, "no numeric labels in forest plot")
  expect_true(any(abs(non_na - 100) < 4 * 20),
              info = paste("Label values:", paste(round(non_na, 1), collapse = ", ")))
})
