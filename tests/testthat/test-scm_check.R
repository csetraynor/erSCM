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
