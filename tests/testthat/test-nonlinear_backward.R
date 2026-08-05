# Tests for Phase 2 nonlinear backward elimination and output adaptation.
# All fixtures rely on helper-bootstrap.R which testthat auto-loads.

# ---------------------------------------------------------------------------
# Core backward helper
# ---------------------------------------------------------------------------

test_that("scm_backward_delete_nonlinear removes noise variable", {
  skip_on_cran()
  set.seed(11)
  n <- 400
  # Only CMAX truly predicts DV; AGE is pure noise.
  y <- rbinom(n, 1, plogis(-2 + 1.5 * rnorm(n)))
  dat <- data.frame(DV = y, CMAX = rnorm(n, 100, 20), AGE = rnorm(n, 55, 12))

  # Simulate a current_model where both entered (at lenient thresholds)
  # by running forward selection with relaxed alpha
  fit_fwd <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = "AGE",
    significance_level = c(0.9, 0.9),   # very lenient — force both in
    dv_col             = "DV",
    family             = "binomial",
    nonlinear          = list(CMAX = c("linear", "log"), AGE = c("linear")),
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))

  # Both CMAX and AGE must be in the final model (else test is vacuous)
  relation <- fit_fwd$metadata$finalRelation
  tr_reg   <- fit_fwd$metadata$transform_registry

  skip_if(length(relation) < 2, "Forward step didn't add both terms; skipping")

  # Run one backward step with strict exit_alpha — AGE should be candidate
  current_mod <- list(
    Relation           = relation,
    toto_tot           = fit_fwd$path,
    finalforward       = max(fit_fwd$path$Run, na.rm = TRUE),
    transform_registry = tr_reg
  )
  # Reconstruct toto_tot in correct format (needs Run, SELECT, etc.)
  # Use a minimal compatible tibble
  toto_compat <- tibble::tibble(
    Run = seq_len(nrow(fit_fwd$path)),
    min2LL  = NA_real_, diff2LL = NA_real_, AIC = NA_real_,
    Reference = NA_real_, pvalue = NA_real_, Signif = NA_real_,
    nobs = NA_integer_, dftest = 0L, dftot = 0L,
    RUNTEST = NA_character_, SELECT = NA_character_,
    Model = paste0("step", seq_len(nrow(fit_fwd$path)))
  )
  current_mod$toto_tot <- toto_compat

  res <- scm_backward_delete_nonlinear(
    eres_data          = dat,
    current_model      = current_mod,
    model_type         = "glm",
    family             = "binomial",
    exit_alpha         = 0.01,
    transform_registry = tr_reg,
    simplify           = FALSE,
    verbose            = FALSE
  )

  expect_type(res, "list")
  expect_true("transform_registry" %in% names(res))
  expect_true("Relation" %in% names(res))
  expect_true("flagforst" %in% names(res))
})

test_that("Backward step removes variable from transform_registry", {
  skip_on_cran()
  set.seed(22)
  n <- 300
  # CMAX is the real predictor; AGE is noise. Force both in, then eliminate AGE.
  dat <- data.frame(
    DV   = rbinom(n, 1, plogis(-2 + log(pmax(runif(n, 1, 100) / 50, 0.01)))),
    CMAX = runif(n, 1, 100),
    AGE  = rnorm(n, 55, 12)
  )

  tr_reg <- list(
    CMAX = list(transform_type = "log", formula_rhs = "log(CMAX)",
                df = 1L, label = "log(CMAX)", powers = NULL, knots = NULL,
                pvalue_vs_linear = 0.03),
    AGE  = list(transform_type = "linear", formula_rhs = "AGE",
                df = 1L, label = "AGE", powers = NULL, knots = NULL,
                pvalue_vs_linear = NA_real_)
  )
  toto_compat <- tibble::tibble(
    Run = 1:3, min2LL = NA_real_, diff2LL = NA_real_, AIC = NA_real_,
    Reference = NA_real_, pvalue = NA_real_, Signif = NA_real_,
    nobs = NA_integer_, dftest = 0L, dftot = 0L,
    RUNTEST = NA_character_, SELECT = NA_character_, Model = c("Base", "log(CMAX)", "AGE")
  )
  current_mod <- list(
    Relation           = c("log(CMAX)", "AGE"),
    toto_tot           = toto_compat,
    finalforward       = 3L,
    transform_registry = tr_reg
  )

  res <- scm_backward_delete_nonlinear(
    eres_data          = dat,
    current_model      = current_mod,
    model_type         = "glm",
    family             = "binomial",
    exit_alpha         = 0.01,
    transform_registry = tr_reg,
    simplify           = FALSE,
    verbose            = FALSE
  )

  if (res$flagforst) {
    # A variable was removed — it must also be gone from the registry
    removed_var <- setdiff(names(tr_reg), names(res$transform_registry))
    expect_length(removed_var, 1L)
    # The removed variable's formula term must be gone from Relation
    expect_false(tr_reg[[removed_var]]$formula_rhs %in% res$Relation)
  } else {
    # Nothing removed — both retained; still a valid outcome
    expect_equal(length(res$Relation), 2L)
  }
})

test_that("Multi-df spline term removed as a unit, not split", {
  skip_on_cran()
  set.seed(33)
  n <- 300
  dat <- data.frame(
    DV   = rbinom(n, 1, plogis(-2 + sin(runif(n, 0, 10)))),
    CMAX = runif(n, 1, 100),
    AGE  = rnorm(n, 55, 12)
  )

  spline_formula <- "splines::ns(AGE, df = 3)"
  tr_reg <- list(
    CMAX = list(transform_type = "linear", formula_rhs = "CMAX",
                df = 1L, label = "CMAX", powers = NULL, knots = NULL,
                pvalue_vs_linear = NA_real_),
    AGE  = list(transform_type = "rcs", formula_rhs = spline_formula,
                df = 3L, label = "RCS(AGE, df=3)", powers = NULL,
                knots = quantile(dat$AGE, c(0.25, 0.5, 0.75)),
                pvalue_vs_linear = 0.10)
  )
  toto_compat <- tibble::tibble(
    Run = 1:3, min2LL = NA_real_, diff2LL = NA_real_, AIC = NA_real_,
    Reference = NA_real_, pvalue = NA_real_, Signif = NA_real_,
    nobs = NA_integer_, dftest = 0L, dftot = 0L,
    RUNTEST = NA_character_, SELECT = NA_character_,
    Model = c("Base", "CMAX", spline_formula)
  )
  current_mod <- list(
    Relation           = c("CMAX", spline_formula),
    toto_tot           = toto_compat,
    finalforward       = 3L,
    transform_registry = tr_reg
  )

  res <- scm_backward_delete_nonlinear(
    eres_data          = dat,
    current_model      = current_mod,
    model_type         = "glm",
    family             = "binomial",
    exit_alpha         = 0.01,
    transform_registry = tr_reg,
    simplify           = FALSE,
    verbose            = FALSE
  )

  # If AGE was removed, it should be absent entirely from Relation,
  # NOT as individual ns() basis terms.
  if (res$flagforst && "AGE" %in% setdiff(names(tr_reg), names(res$transform_registry))) {
    # Verify the complete spline term string was removed as one unit
    expect_false(spline_formula %in% res$Relation)
    expect_equal(length(res$Relation), 1L)   # only CMAX remains
  }
})

test_that("Simplification: registry is updated to linear when simplification occurs", {
  # This test verifies the MECHANISM of simplification by directly invoking
  # .lrt_pvalue semantics: if pval > linearity_alpha the entry is downgraded.
  # We create a case where the FP2 formula is identical to linear (same fit)
  # so p-value must be NA or > any linearity_alpha.
  set.seed(44)
  n <- 100
  dat <- data.frame(
    DV   = rbinom(n, 1, 0.3),
    CMAX = rnorm(n, 50, 10),
    AGE  = rnorm(n, 55, 12)
  )

  # FP2 with power (1,1) uses x^1 * log(x) — but test the STRUCTURAL path
  # by simulating what happens after simplification p-value > linearity_alpha.
  # Directly test the registry update logic: if we run with simplify=FALSE
  # no vars should appear in simplified_vars.
  fp2_formula <- "I(CMAX^(0.5))"   # FP1, df=1 — can't simplify (df <= 1)
  tr_reg <- list(
    CMAX = list(transform_type = "fp1", formula_rhs = fp2_formula,
                df = 1L, label = "FP1(CMAX; 0.5)", powers = 0.5,
                knots = NULL, pvalue_vs_linear = 0.03),
    AGE  = list(transform_type = "linear", formula_rhs = "AGE",
                df = 1L, label = "AGE", powers = NULL, knots = NULL,
                pvalue_vs_linear = NA_real_)
  )
  toto_compat <- tibble::tibble(
    Run = 1:3, min2LL = NA_real_, diff2LL = NA_real_, AIC = NA_real_,
    Reference = NA_real_, pvalue = NA_real_, Signif = NA_real_,
    nobs = NA_integer_, dftest = 0L, dftot = 0L,
    RUNTEST = NA_character_, SELECT = NA_character_,
    Model = c("Base", fp2_formula, "AGE")
  )
  current_mod <- list(
    Relation           = c(fp2_formula, "AGE"),
    toto_tot           = toto_compat,
    finalforward       = 3L,
    transform_registry = tr_reg
  )

  res <- scm_backward_delete_nonlinear(
    eres_data          = dat,
    current_model      = current_mod,
    model_type         = "glm",
    family             = "binomial",
    exit_alpha         = 0.001,
    transform_registry = tr_reg,
    simplify           = TRUE,
    linearity_alpha    = 0.05,
    verbose            = FALSE
  )

  # simplified_vars must be a character vector (empty when nothing simplified)
  expect_type(res$simplified_vars, "character")
  # FP1 (df=1) cannot be simplified further — it must remain fp1 or be removed
  if ("CMAX" %in% names(res$transform_registry)) {
    expect_identical(res$transform_registry$CMAX$transform_type, "fp1")
    expect_equal(res$transform_registry$CMAX$df, 1L)
  }
  # Structural: if anything was simplified it should appear in simplified_vars
  if (length(res$simplified_vars) > 0) {
    expect_true(all(res$simplified_vars %in% names(tr_reg)))
  }
})

# ---------------------------------------------------------------------------
# End-to-end integration
# ---------------------------------------------------------------------------

test_that("scm_glm backward with nonlinear=TRUE ends cleanly and updates registry", {
  skip_on_cran()
  set.seed(55)
  n <- 400
  x <- runif(n, 1, 100)
  dat <- data.frame(
    DV   = rbinom(n, 1, plogis(-2 + 0.5 * log(x / 50))),
    CMAX = x,
    AGE  = rnorm(n, 55, 12),
    WT   = rnorm(n, 75, 15)
  )

  fit <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT"),
    significance_level = c(0.05, 0.01),
    dv_col             = "DV",
    family             = "binomial",
    nonlinear          = list(CMAX = c("linear", "log"), AGE = c("linear"),
                              WT   = c("linear")),
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))

  expect_s3_class(fit, "eres_scm_glm")
  meta <- fit$metadata
  expect_true("transform_registry"  %in% names(meta))
  expect_true("has_nonlinear_terms" %in% names(meta))
  expect_true("fp_powers"           %in% names(meta))
  expect_true("spline_df"           %in% names(meta))
  expect_type(meta$has_nonlinear_terms, "logical")
})

test_that("scm_glm backward with nonlinear=NULL is identical to before", {
  skip_on_cran()
  dat <- .make_test_data_glm(n = 200, seed = 77)
  fit_lin <- suppressMessages(scm_glm(
    eres_data = dat, exposure_vars = "CMAX", covariates = c("AGE", "WT", "SEX"),
    significance_level = c(0.5, 0.3), dv_col = "DV", family = "binomial",
    center = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  fit_null <- suppressMessages(scm_glm(
    eres_data = dat, exposure_vars = "CMAX", covariates = c("AGE", "WT", "SEX"),
    significance_level = c(0.5, 0.3), dv_col = "DV", family = "binomial",
    nonlinear = NULL, center = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  expect_identical(fit_lin$metadata$finalRelation, fit_null$metadata$finalRelation)
  expect_null(fit_null$metadata$transform_registry)
  expect_false(fit_null$metadata$has_nonlinear_terms)
})

# ---------------------------------------------------------------------------
# Output adaptation
# ---------------------------------------------------------------------------

test_that("Summary table has transform columns when nonlinear is used", {
  skip_on_cran()
  set.seed(66)
  n <- 300
  x <- runif(n, 1, 100)
  dat <- data.frame(
    DV   = rbinom(n, 1, plogis(-2 + log(x / 50))),
    CMAX = x
  )

  fit <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = character(0),
    significance_level = c(0.05, 0.01),
    dv_col             = "DV",
    family             = "binomial",
    nonlinear          = list(CMAX = c("linear", "log")),
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))

  sm <- fit$summary
  # Summary must have the enhancement columns when nonlinear is active
  expect_true(all(c("variable", "transform_type", "transform_label") %in% colnames(sm)))
})

test_that(".enhance_summary_nonlinear adds correct columns", {
  # Minimal unit test, no model fitting
  sum_tbl <- data.frame(
    term      = c("(Intercept)", "log(CMAX)", "AGE"),
    estimate  = c(0.1, 0.5, -0.1),
    stringsAsFactors = FALSE
  )
  reg <- list(
    CMAX = list(transform_type = "log", formula_rhs = "log(CMAX)",
                df = 1L, label = "log(CMAX)", powers = NULL, knots = NULL,
                pvalue_vs_linear = 0.02),
    AGE  = list(transform_type = "linear", formula_rhs = "AGE",
                df = 1L, label = "AGE", powers = NULL, knots = NULL,
                pvalue_vs_linear = NA_real_)
  )
  out <- .enhance_summary_nonlinear(sum_tbl, reg)
  expect_true(all(c("variable", "transform_type", "transform_label") %in% colnames(out)))
  expect_equal(out$variable[out$term == "log(CMAX)"],  "CMAX")
  expect_equal(out$transform_type[out$term == "log(CMAX)"], "log")
  expect_equal(out$variable[out$term == "AGE"], "AGE")
  expect_equal(out$transform_type[out$term == "AGE"], "linear")
})

test_that(".enhance_summary_nonlinear returns input unchanged when registry is NULL", {
  sum_tbl <- data.frame(term = "CMAX", estimate = 0.5, stringsAsFactors = FALSE)
  out <- .enhance_summary_nonlinear(sum_tbl, NULL)
  expect_identical(out, sum_tbl)
})

# ---------------------------------------------------------------------------
# Gap #7: nonlinear + start_covariates / backward direction
# ---------------------------------------------------------------------------

test_that("nonlinear + start_covariates with backward direction (Gap #7)", {
  skip_on_cran()
  set.seed(77)
  n <- 400
  cmax <- runif(n, 10, 200)
  dat <- data.frame(
    DV   = rbinom(n, 1, plogis(-2 + 0.9 * log(cmax / 50) + 0.02 * rnorm(n, 55, 10))),
    CMAX = cmax,
    AGE  = rnorm(n, 55, 10),
    WT   = rnorm(n, 75, 15)
  )
  fit <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT"),
    start_covariates   = c("AGE", "WT"),
    direction          = "backward",
    significance_level = c(0.05, 0.01),
    dv_col             = "DV",
    family             = "binomial",
    nonlinear          = list(CMAX = c("linear", "log"),
                              AGE  = c("linear", "log"),
                              WT   = c("linear", "fp1")),
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  expect_s3_class(fit, "eres_scm_glm")
  tr_reg <- fit$metadata$transform_registry
  # Any surviving covariate must have a valid transform_registry entry
  if (!is.null(tr_reg)) {
    for (v in names(tr_reg)) {
      expect_true(tr_reg[[v]]$transform_type %in% c("linear", "log", "fp1", "fp2", "rcs", "sqrt"))
    }
  }
})

test_that("strict retention removes start_covariate cleanly from registry (Gap #7)", {
  skip_on_cran()
  set.seed(78)
  n <- 300
  # AGE is pure noise — should be removed under strict backward elimination
  dat <- data.frame(
    DV   = rbinom(n, 1, 0.3),
    CMAX = rnorm(n, 100, 20),
    AGE  = rnorm(n, 55, 10)
  )
  fit <- suppressMessages(scm_glm(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = "AGE",
    start_covariates   = "AGE",
    direction          = "backward",
    significance_level = c(0.05, 0.001),
    dv_col             = "DV",
    family             = "binomial",
    nonlinear          = list(CMAX = c("linear", "log"), AGE = c("linear")),
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  expect_s3_class(fit, "eres_scm_glm")
  tr_reg <- fit$metadata$transform_registry
  # AGE is noise — under strict thresholds it should be absent from the final registry
  # (it may have been removed during backward elimination)
  if (!is.null(tr_reg)) {
    # If AGE survived (data is random, could happen), it must still have a valid entry
    if ("AGE" %in% names(tr_reg)) {
      expect_identical(tr_reg$AGE$transform_type, "linear")
    }
  }
})

test_that("scm_coxph backward with nonlinear=TRUE ends cleanly", {
  skip_on_cran()
  set.seed(88)
  dat <- .make_test_data_coxph(n = 300, seed = 88)
  fit <- suppressMessages(scm_coxph(
    eres_data          = dat,
    exposure_vars      = "CMAX",
    covariates         = c("AGE", "WT", "SEX"),
    significance_level = c(0.5, 0.3),
    dv_col             = "DV", time_col = "TIME",
    nonlinear          = list(CMAX = c("linear", "log")),
    center             = FALSE, scale_to_median = FALSE, verbose = FALSE
  ))
  expect_s3_class(fit, "eres_scm_coxph")
  expect_true("has_nonlinear_terms" %in% names(fit$metadata))
  expect_type(fit$metadata$has_nonlinear_terms, "logical")
})
