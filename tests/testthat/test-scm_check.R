test_that(".check_inputs validates correctly", {
  d <- data.frame(DV = 1:5, x1 = rnorm(5), x2 = rnorm(5))

  # Should not error with valid inputs
 expect_message(
    .check_inputs(modds = d, dvid_col = NULL, dvid = NULL, dv_col = "DV", covariates = c("x1", "x2"))
  )

  # Should error when covariates not in data
  expect_error(
    .check_inputs(modds = d, dvid_col = NULL, dvid = NULL, dv_col = "DV", covariates = c("x1", "MISSING"))
  )

  # Should error when dvid provided without dvid_col
  expect_error(
    .check_inputs(modds = d, dvid_col = NULL, dvid = "A", dv_col = "DV", covariates = c("x1"))
  )
})

test_that(".check_forward returns a scalar finalforward for a single selected row", {
  run_summary <- data.frame(
    Run = c(1, 2, 3),
    Model = c("a", "b", "SEL"),
    SELECT = "SEL",
    stringsAsFactors = FALSE
  )

  res <- .check_forward(run_summary, relation = "base")

  expect_type(res$finalforward, "double")
  expect_length(res$finalforward, 1L)
  expect_identical(res$finalforward, 3)
  expect_true(res$flagforst)
  expect_true("SEL" %in% res$Relation)
})

test_that(".check_forward returns the first match (scalar) when several rows tie on SELECT", {
  run_summary <- data.frame(
    Run = c(5, 7),
    Model = c("SEL", "SEL"),
    SELECT = "SEL",
    stringsAsFactors = FALSE
  )

  res <- .check_forward(run_summary, relation = "base")

  # Regression guard for the pre-fix bug where a length > 1 vector was returned,
  # which cascaded into downstream length-1 validators. Must be a scalar.
  expect_length(res$finalforward, 1L)
  expect_identical(res$finalforward, 5)
  expect_true(res$flagforst)
})

test_that(".check_forward returns scalar NA_real_ when nothing is selected", {
  run_summary <- data.frame(
    Run = c(1, 2),
    Model = c("a", "b"),
    SELECT = "SEL",
    stringsAsFactors = FALSE
  )

  expect_message(
    res <- .check_forward(run_summary, relation = "base"),
    "not significant"
  )

  expect_type(res$finalforward, "double")
  expect_length(res$finalforward, 1L)
  expect_true(is.na(res$finalforward))
  expect_false(res$flagforst)
})
