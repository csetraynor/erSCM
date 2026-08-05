test_that("get_sum_tbl returns expected structure", {
  set.seed(42)
  d <- data.frame(y = rbinom(100, 1, 0.3), x1 = rnorm(100), x2 = rnorm(100))
  mod <- glm(y ~ x1 + x2, data = d, family = binomial)

  result <- get_sum_tbl(mod, exponentiate = TRUE)
  expect_s3_class(result, "data.frame")
  expect_true("term" %in% colnames(result))
  expect_true(nrow(result) >= 3)
})

test_that("get_sum_tbl_coxph returns expected structure", {
  skip_if_not_installed("survival")
  set.seed(42)
  d <- data.frame(time = rexp(100), status = rbinom(100, 1, 0.5), x1 = rnorm(100))
  mod <- survival::coxph(survival::Surv(time, status) ~ x1, data = d)

  result <- get_sum_tbl_coxph(mod, exponentiate = TRUE, include_conf_int = TRUE)
  expect_s3_class(result, "data.frame")
  expect_true("term" %in% colnames(result))
})

test_that("run_sum_tbl fits a glm", {
  set.seed(42)
  d <- data.frame(DV = rbinom(100, 1, 0.3), x1 = rnorm(100))
  mod <- run_sum_tbl(d, covariates = "x1", family = "binomial")
  expect_true(inherits(mod, "glm"))
})
