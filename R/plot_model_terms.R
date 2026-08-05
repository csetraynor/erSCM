#' Plot covariate effects -- dispatcher for coxph and logistic models
#'
#' Routes to \code{\link{plot_cox_terms}} for \code{coxph} fits or
#' \code{\link{plot_logistic_terms}} for binomial \code{glm} fits.
#' All extra arguments are forwarded to the underlying function.
#'
#' @param fit A fit object returned by \code{scm_coxph()} or \code{scm_glm()}.
#' @param covariate_term Character; name of the main-effect covariate to plot.
#' @param ... Additional arguments passed to \code{plot_cox_terms()} or
#'   \code{plot_logistic_terms()}.
#'
#' @return A \code{ggplot2} object (or data, when \code{plot = FALSE} is passed
#'   in \code{...}).
#'
#' @seealso \code{\link{plot_cox_terms}}, \code{\link{plot_logistic_terms}}
#'
#' @examples
#' data(demo_eres)
#' demo_eres$SEX <- as.factor(demo_eres$SEX)
#' lr_fit <- scm_glm(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", covariates = c("AGE", "WT"),
#'   dv_col = "sim_status", significance_level = c(0.2, 0.2)
#' )
#' plot_model_terms(lr_fit, "WT")
#'
#' @export
plot_model_terms <- function(fit, covariate_term, ...) {
  model <- fit$final_model
  if (inherits(model, "coxph")) {
    plot_cox_terms(fit, covariate_term, ...)
  } else if (inherits(model, "glm") && model$family$family == "binomial") {
    plot_logistic_terms(fit, covariate_term, ...)
  } else {
    cli::cli_abort("Model class {.cls {class(model)[1]}} is not yet supported by {.fun plot_model_terms}.")
  }
}
