test_that("parse_nonlinear_spec returns NULL for NULL input", {
  dat <- .make_test_data_glm()
  expect_null(parse_nonlinear_spec(NULL, c("AGE", "WT"), "CMAX", dat))
})

test_that("parse_nonlinear_spec returns NULL for FALSE input", {
  dat <- .make_test_data_glm()
  expect_null(parse_nonlinear_spec(FALSE, c("AGE", "WT"), "CMAX", dat))
})

test_that("parse_nonlinear_spec handles TRUE / auto", {
  dat <- .make_test_data_glm()
  spec <- parse_nonlinear_spec(TRUE, c("AGE", "WT", "SEX"), "CMAX", dat)
  expect_s3_class(spec, "nonlinear_spec")
  # Numeric variables get the default transform set
  expect_true("fp1" %in% spec[["AGE"]])
  expect_true("log" %in% spec[["CMAX"]])
  # Factor variable is forced to linear
  expect_identical(spec[["SEX"]], "linear")

  # "auto" string works the same way
  spec2 <- parse_nonlinear_spec("auto", c("AGE", "WT", "SEX"), "CMAX", dat)
  expect_identical(spec, spec2)
})

test_that("parse_nonlinear_spec handles character vector", {
  dat <- .make_test_data_glm()
  spec <- parse_nonlinear_spec(c("linear", "log", "fp1"),
                               c("AGE", "WT", "SEX"), "CMAX", dat)
  expect_s3_class(spec, "nonlinear_spec")
  expect_true("log" %in% spec[["AGE"]])
  expect_true("fp1" %in% spec[["WT"]])
  # Factor still forced to linear
  expect_identical(spec[["SEX"]], "linear")
})

test_that("parse_nonlinear_spec handles named list", {
  dat <- .make_test_data_glm()
  spec <- parse_nonlinear_spec(
    list(CMAX = c("linear", "log", "fp2"), AGE = c("linear", "log")),
    covariates = c("AGE", "WT", "SEX"),
    exposure_vars = "CMAX",
    data = dat
  )
  expect_s3_class(spec, "nonlinear_spec")
  expect_identical(spec[["CMAX"]], c("linear", "log", "fp2"))
  expect_identical(spec[["AGE"]],  c("linear", "log"))
  # Unlisted variables default to "linear"
  expect_identical(spec[["WT"]],  "linear")
})

test_that("parse_nonlinear_spec warns for non-numeric variables with explicit nonlinear spec", {
  dat <- .make_test_data_glm()
  expect_warning(
    parse_nonlinear_spec(
      list(SEX = c("linear", "log")),
      covariates = c("AGE", "SEX"),
      exposure_vars = "CMAX",
      data = dat
    ),
    regexp = "not numeric"
  )
})

test_that("parse_nonlinear_spec errors on unknown transform names", {
  dat <- .make_test_data_glm()
  expect_error(
    parse_nonlinear_spec(c("linear", "cubic"), c("AGE"), "CMAX", dat),
    regexp = "Unknown transformation"
  )
})

test_that("parse_nonlinear_spec accepts rcsN names", {
  dat <- .make_test_data_glm()
  spec <- parse_nonlinear_spec(c("linear", "rcs4"), c("AGE"), "CMAX", dat)
  expect_true("rcs4" %in% spec[["AGE"]])
})

test_that("validate_nonlinear_spec warns for log of non-positive values", {
  dat <- .make_test_data_glm()
  dat$AGE[1:5] <- -1   # introduce negatives
  spec <- parse_nonlinear_spec(list(AGE = c("linear", "log")),
                               "AGE", "CMAX", dat)
  expect_warning(
    validate_nonlinear_spec(spec, dat),
    regexp = "non-positive"
  )
})

test_that("validate_nonlinear_spec warns for sqrt of negative values", {
  dat <- .make_test_data_glm()
  dat$WT[1:3] <- -5
  spec <- parse_nonlinear_spec(list(WT = c("linear", "sqrt")),
                               "WT", "CMAX", dat)
  expect_warning(
    validate_nonlinear_spec(spec, dat),
    regexp = "negative"
  )
})

test_that("validate_nonlinear_spec warns for variables with too few unique values", {
  dat <- .make_test_data_glm(n = 50)
  dat$SPARSE <- rep(1:3, length.out = nrow(dat))
  spec <- parse_nonlinear_spec(list(SPARSE = c("linear", "fp1")),
                               "SPARSE", "CMAX", dat)
  expect_warning(
    validate_nonlinear_spec(spec, dat),
    regexp = "unique"
  )
})
