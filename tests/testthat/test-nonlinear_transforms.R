test_that(".fp_term generates correct formula for various powers", {
  expect_identical(.fp_term("X", 0),    "log(X)")
  expect_identical(.fp_term("X", 1),    "X")
  expect_identical(.fp_term("X", 2),    "I(X^(2))")
  expect_identical(.fp_term("X", -1),   "I(X^(-1))")
  expect_identical(.fp_term("X", 0.5),  "I(X^(0.5))")
  expect_identical(.fp_term("X", -2),   "I(X^(-2))")
})

test_that(".fp_term aborts on non-finite power (NA / Inf / -Inf)", {
  # A valid finite power still succeeds.
  expect_identical(.fp_term("X", 3), "I(X^(3))")
  # Non-finite powers previously produced an invalid term or `if (NA)` error.
  expect_error(.fp_term("X", NA_real_), "finite numeric")
  expect_error(.fp_term("X", Inf),      "finite numeric")
  expect_error(.fp_term("X", -Inf),     "finite numeric")
})

test_that(".find_best_fp rejects fp_powers outside the allowed FP power set", {
  dat <- .make_test_data_glm(n = 50, seed = 1)
  expect_error(
    .find_best_fp("AGE", degree = 1L, current_rhs = "1", data = dat,
                  response_lhs = "DV", fp_powers = c(0.5, 4)),
    "allowed FP power set"
  )
})

test_that(".fp2_terms handles distinct powers", {
  terms <- .fp2_terms("AGE", -1, 2)
  expect_length(terms, 2L)
  expect_identical(terms[1L], "I(AGE^(-1))")
  expect_identical(terms[2L], "I(AGE^(2))")
})

test_that(".fp2_terms handles repeated powers (p1 == p2)", {
  # p1 == p2 != 0: second term is x^p * log(x)
  terms <- .fp2_terms("WT", 2, 2)
  expect_length(terms, 2L)
  expect_true(grepl("log", terms[2L]))
  # p1 == p2 == 0: second term is (log x)^2
  terms0 <- .fp2_terms("WT", 0, 0)
  expect_true(grepl("log.*\\^2", terms0[2L]))
})

test_that(".make_simple_transform produces correct formula fragments", {
  lin  <- .make_simple_transform("CMAX", "linear")
  expect_identical(lin$formula_rhs, "CMAX")
  expect_equal(lin$df, 1L)

  lg   <- .make_simple_transform("CMAX", "log")
  expect_identical(lg$formula_rhs, "log(CMAX)")
  expect_equal(lg$df, 1L)

  sq   <- .make_simple_transform("CMAX", "sqrt")
  expect_identical(sq$formula_rhs, "sqrt(CMAX)")
  expect_equal(sq$df, 1L)
})

test_that(".make_spline_transform produces valid ns() terms", {
  dat <- .make_test_data_glm(n = 200)
  sp  <- .make_spline_transform("AGE", df = 3L, data = dat)
  expect_true(grepl("splines::ns", sp$formula_rhs))
  expect_equal(sp$df, 3L)
  expect_length(sp$knots, 3L)
  expect_true(all(is.numeric(sp$knots)))
})

test_that(".make_spline_transform accepts rcs4 via .rcs_df", {
  expect_equal(.rcs_df("rcs4"), 4L)
  expect_equal(.rcs_df("rcs"),  3L)    # default
  expect_equal(.rcs_df("rcs5"), 5L)
})

test_that(".find_best_fp returns sensible FP1 result", {
  skip_on_cran()
  dat  <- .make_test_data_glm(n = 200, seed = 1)
  best <- .find_best_fp("AGE", degree = 1L, current_rhs = "1",
                         data = dat, response_lhs = "DV",
                         model_type = "glm", family = "binomial")
  expect_type(best, "list")
  expect_equal(best$df, 1L)
  expect_length(best$powers, 1L)
  expect_true(is.finite(best$deviance))
  expect_match(best$label, "FP1")
})

test_that(".find_best_fp returns sensible FP2 result", {
  skip_on_cran()
  dat  <- .make_test_data_glm(n = 200, seed = 2)
  best <- .find_best_fp("AGE", degree = 2L, current_rhs = "1",
                         data = dat, response_lhs = "DV",
                         model_type = "glm", family = "binomial")
  expect_type(best, "list")
  expect_equal(best$df, 2L)
  expect_length(best$powers, 2L)
  expect_match(best$label, "FP2")
})

test_that(".find_best_fp works for coxph", {
  skip_on_cran()
  dat  <- .make_test_data_coxph(n = 200, seed = 3)
  best <- .find_best_fp("AGE", degree = 1L, current_rhs = "1",
                         data = dat, response_lhs = "survival::Surv(TIME, DV)",
                         model_type = "coxph")
  expect_type(best, "list")
  expect_equal(best$df, 1L)
})

test_that(".fit_model returns NULL on bad formula", {
  dat <- .make_test_data_glm()
  fit <- .fit_model("DV", "NONEXISTENT_COL", dat, "glm", "binomial")
  expect_null(fit)
})

test_that(".fit_model fits glm successfully", {
  dat <- .make_test_data_glm()
  fit <- .fit_model("DV", "AGE", dat, "glm", "binomial")
  expect_s3_class(fit, "glm")
})

test_that(".fit_model fits coxph successfully", {
  skip_on_cran()
  dat <- .make_test_data_coxph()
  fit <- .fit_model("survival::Surv(TIME, DV)", "AGE", dat, "coxph")
  expect_s3_class(fit, "coxph")
})
