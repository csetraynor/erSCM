#' Print method for SCM GLM results
#'
#' Displays a concise summary of the stepwise covariate modelling result
#' from [scm_glm()].
#'
#' @param x An object of class `eres_scm_glm`.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#' @examples
#' data(demo_eres)
#' fit <- scm_glm(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", covariates = c("AGE", "WT"),
#'   dv_col = "sim_status", significance_level = c(0.2, 0.2)
#' )
#' print(fit)
#' @export
print.eres_scm_glm <- function(x, ...) {
  cli_header <- "Stepwise Covariate Model (GLM)"
  cat(cli_header, "\n")
  cat(strrep("-", nchar(cli_header)), "\n")

  if (is.null(x$final_model)) {
    cat("No exposure variable reached the entry significance threshold.\n")
    cat("Model was not fit.\n")
    n_tested <- length(x$metadata$candidate_exposure_vars)
    cat("Exposure candidates tested: ", n_tested, "\n")
    if (!is.null(x$metadata$dvid)) {
      cat("DVID:                ", x$metadata$dvid, "\n")
    }
    return(invisible(x))
  }

  n_exp_tested <- length(x$metadata$candidate_exposure_vars)
  n_cov_tested <- length(x$metadata$candidate_covariates)
  n_selected <- length(x$metadata$finalRelation)
  fam <- x$final_model$family$family
  aic_val <- stats::AIC(x$final_model)

  cat("Family:              ", fam, "\n")
  cat("Exposure candidates: ", n_exp_tested, "\n")
  cat("Covariate candidates:", n_cov_tested, "\n")
  cat("Terms in final model:", n_selected, "\n")

  cat("Final covariates:\n")
  for (rel in x$metadata$finalRelation) {
    lbl <- get_var_label(x$var_spec, rel)
    unt <- get_var_units(x$var_spec, rel)
    display <- if (!is.na(unt) && !is.null(unt)) {
      paste0("  ", lbl, " (", unt, ")")
    } else {
      paste0("  ", lbl)
    }
    cat(display, "\n")
  }

  cat("AIC:                 ", round(aic_val, 2), "\n")
  if (!is.null(x$metadata$dvid)) {
    cat("DVID:                ", x$metadata$dvid, "\n")
  }

  invisible(x)
}

#' Print method for SCM Cox PH results
#'
#' Displays a concise summary of the stepwise covariate modelling result
#' from [scm_coxph()].
#'
#' @param x An object of class `eres_scm_coxph`.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#' @examples
#' data(demo_eres)
#' fit <- scm_coxph(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", covariates = c("AGE", "WT"),
#'   dv_col = "sim_status", time_col = "sim_time", significance_level = c(0.2, 0.2)
#' )
#' print(fit)
#' @export
print.eres_scm_coxph <- function(x, ...) {
  cli_header <- "Stepwise Covariate Model (Cox PH)"
  cat(cli_header, "\n")
  cat(strrep("-", nchar(cli_header)), "\n")

  if (is.null(x$final_model)) {
    cat("No exposure variable reached the entry significance threshold.\n")
    cat("Model was not fit.\n")
    n_tested <- length(x$metadata$candidate_exposure_vars)
    cat("Exposure candidates tested: ", n_tested, "\n")
    if (!is.null(x$metadata$dvid_name)) {
      cat("DVID:                ", x$metadata$dvid_name, "\n")
    }
    return(invisible(x))
  }

  n_exp_tested <- length(x$metadata$candidate_exposure_vars)
  n_cov_tested <- length(x$metadata$candidate_covariates)
  n_selected <- length(x$metadata$finalRelation)
  aic_val <- stats::AIC(x$final_model)

  cat("Model:               Cox Proportional Hazards\n")
  cat("Exposure candidates: ", n_exp_tested, "\n")
  cat("Covariate candidates:", n_cov_tested, "\n")
  cat("Terms in final model:", n_selected, "\n")

  cat("Final covariates:\n")
  for (rel in x$metadata$finalRelation) {
    lbl <- get_var_label(x$var_spec, rel)
    unt <- get_var_units(x$var_spec, rel)
    display <- if (!is.na(unt) && !is.null(unt)) {
      paste0("  ", lbl, " (", unt, ")")
    } else {
      paste0("  ", lbl)
    }
    cat(display, "\n")
  }

  cat("AIC:                 ", round(aic_val, 2), "\n")
  if (!is.null(x$metadata$dvid_name)) {
    cat("DVID:                ", x$metadata$dvid_name, "\n")
  }

  invisible(x)
}

#' Plot method for SCM GLM results
#'
#' Produces a forest plot of covariate effects from a [scm_glm()] result.
#'
#' @param x An object of class `eres_scm_glm`.
#' @param type Character. Type of plot to produce. Currently only `"forest"`
#'   is supported (the default).
#' @param ... Additional arguments passed to [make_forest_plot()].
#'
#' @return A `ggplot`/`patchwork` object (invisibly).
#' @examples
#' data(demo_eres)
#' fit <- scm_glm(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", covariates = c("AGE", "WT"),
#'   dv_col = "sim_status", significance_level = c(0.2, 0.2)
#' )
#' plot(fit)
#' @export
plot.eres_scm_glm <- function(x, type = "forest", ...) {
  type <- match.arg(type, choices = c("forest"))
  p <- make_forest_plot(fit = x, ...)
  print(p)
  invisible(p)
}

#' Plot method for SCM Cox PH results
#'
#' Produces a forest plot of covariate effects from a [scm_coxph()] result.
#'
#' @param x An object of class `eres_scm_coxph`.
#' @param type Character. Type of plot to produce. Currently only `"forest"`
#'   is supported (the default).
#' @param ... Additional arguments passed to [make_forest_plot()].
#'
#' @return A `ggplot`/`patchwork` object (invisibly).
#' @examples
#' data(demo_eres)
#' fit <- scm_coxph(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", covariates = c("AGE", "WT"),
#'   dv_col = "sim_status", time_col = "sim_time", significance_level = c(0.2, 0.2)
#' )
#' plot(fit)
#' @export
plot.eres_scm_coxph <- function(x, type = "forest", ...) {
  type <- match.arg(type, choices = c("forest"))
  p <- make_forest_plot(fit = x, ...)
  print(p)
  invisible(p)
}
