test_that(".anova_pvalue extracts the p-value from the glm anova layout", {
  # glm anova(test = "Chisq") uses the column name "Pr(>Chi)".
  tab <- data.frame(
    check.names = FALSE,
    "Resid. Df" = c(99, 98),
    "Resid. Dev" = c(120.5, 110.2),
    "Df" = c(NA, 1),
    "Deviance" = c(NA, 10.3),
    "Pr(>Chi)" = c(NA, 0.0013)
  )
  expect_equal(.anova_pvalue(tab), 0.0013)
})

test_that(".anova_pvalue extracts the p-value from the coxph anova layout", {
  # survival::coxph anova(test = "Chisq") uses the column name "P(>|Chi|)".
  tab <- data.frame(
    check.names = FALSE,
    "loglik" = c(-200.1, -195.4),
    "Chisq" = c(NA, 9.4),
    "Df" = c(NA, 1),
    "P(>|Chi|)" = c(NA, 0.0021)
  )
  expect_equal(.anova_pvalue(tab), 0.0021)
})

test_that(".anova_pvalue also accepts the 'Pr(>|Chi|)' spelling", {
  tab <- data.frame(
    check.names = FALSE,
    "loglik" = c(-200.1, -195.4),
    "Chisq" = c(NA, 9.4),
    "Df" = c(NA, 1),
    "Pr(>|Chi|)" = c(NA, 0.0047)
  )
  expect_equal(.anova_pvalue(tab), 0.0047)
})

test_that(".anova_pvalue aborts when no recognized p-value column is present", {
  tab <- data.frame(
    check.names = FALSE,
    "Df" = c(NA, 1),
    "Deviance" = c(NA, 10.3),
    "P.value" = c(NA, 0.01)
  )
  expect_error(.anova_pvalue(tab), "p-value column")
})
