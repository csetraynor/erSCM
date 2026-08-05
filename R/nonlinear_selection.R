# Nonlinear selection logic
# Determines the best functional form for a candidate variable at each SCM
# forward step, following the MFP closed-testing approach.

#' Compute LRT p-value between two fitted models
#'
#' @param fit_full Fitted model (larger).
#' @param fit_reduced Fitted model (smaller / null). May be NULL for
#'   deviance-only comparisons.
#' @param df_diff Integer, degrees-of-freedom difference. When NULL, computed
#'   from model df automatically.
#' @return Numeric p-value, or NA_real_ if either model is NULL or deviances
#'   are non-finite.
#' @keywords internal
.lrt_pvalue <- function(fit_full, fit_reduced, df_diff = NULL) {
  if (is.null(fit_full) || is.null(fit_reduced)) return(NA_real_)
  dev_full    <- .model_deviance(fit_full)
  dev_reduced <- .model_deviance(fit_reduced)
  if (!is.finite(dev_full) || !is.finite(dev_reduced)) return(NA_real_)
  chi2 <- dev_reduced - dev_full
  if (chi2 < 0) return(NA_real_)
  if (is.null(df_diff)) {
    df_full    <- tryCatch(as.integer(attr(stats::logLik(fit_full),    "df")), error = function(e) NA_integer_)
    df_reduced <- tryCatch(as.integer(attr(stats::logLik(fit_reduced), "df")), error = function(e) NA_integer_)
    df_diff <- df_full - df_reduced
  }
  if (is.na(df_diff) || df_diff <= 0L) return(NA_real_)
  stats::pchisq(chi2, df = df_diff, lower.tail = FALSE)
}


#' Select Best Transformation for a Variable at a Forward Step
#'
#' Given the current model, a candidate variable, and its allowed
#' transformations, determines the best functional form and returns the LRT
#' test result vs the current model.
#'
#' @param var_name Character.
#' @param transforms Character vector of allowed transformations for this
#'   variable.
#' @param current_rhs Character, the current model's RHS formula (terms
#'   already in the model, not including the candidate variable). Use
#'   \code{"1"} for intercept-only baseline.
#' @param response_lhs Character, LHS formula string (e.g., \code{"DV"} or
#'   \code{"survival::Surv(TIME, DV)"}).
#' @param data data.frame.
#' @param model_type Character, \code{"glm"} or \code{"coxph"}.
#' @param family Family for GLM.
#' @param fp_powers FP power set.
#' @param spline_df Spline degrees of freedom (used when \code{"rcs"} or
#'   \code{"rcsN"} appears without an embedded number).
#' @param linearity_alpha Numeric. Significance level for preferring a
#'   non-linear form over linear. Default \code{0.05}.
#'
#' @return A list with:
#'   \describe{
#'     \item{var_name}{The variable.}
#'     \item{best_transform}{Character, the winning transformation type.}
#'     \item{formula_rhs}{Character, the formula terms to add.}
#'     \item{df}{Integer, degrees of freedom of the best form.}
#'     \item{label}{Human-readable label.}
#'     \item{pvalue_vs_null}{Numeric, LRT p-value of best form vs current model.}
#'     \item{pvalue_vs_linear}{Numeric, LRT p-value of best form vs linear form.
#'       \code{NA} if best is linear.}
#'     \item{aic_best}{AIC of the model with best form.}
#'     \item{aic_linear}{AIC of the model with linear form.}
#'     \item{powers}{FP powers if applicable, else \code{NULL}.}
#'     \item{knots}{Spline knots if applicable, else \code{NULL}.}
#'   }
#'
#' @details
#' The selection procedure follows the MFP closed testing approach:
#' \enumerate{
#'   \item Fit the current model + variable in ALL candidate forms.
#'   \item Among non-linear forms: select the one with lowest deviance.
#'   \item Compare best non-linear vs linear:
#'     \itemize{
#'       \item Non-linear is preferred only when LRT p < \code{linearity_alpha}.
#'       \item Otherwise: linear is used.
#'     }
#'   \item Compare the chosen form vs the current model (no variable):
#'     this p-value feeds into the standard SCM forward selection criterion.
#' }
#'
#' @keywords internal
.select_best_transform <- function(var_name, transforms,
                                    current_rhs, response_lhs = "DV",
                                    data,
                                    model_type = "glm", family = "binomial",
                                    fp_powers = .FP_POWERS,
                                    spline_df = 3L,
                                    linearity_alpha = 0.05) {

  # Ensure "linear" is always among tested forms (needed for comparison)
  transforms <- unique(c("linear", transforms))

  # --- Fit the baseline (current model, no candidate variable) ---
  base_rhs <- if (nchar(trimws(current_rhs)) == 0L || current_rhs == "1") {
    "1"
  } else {
    current_rhs
  }
  fit_null <- .fit_model(response_lhs, base_rhs, data, model_type, family)

  # --- Fit model + variable in linear form ---
  linear_rhs <- if (base_rhs == "1") var_name else paste(base_rhs, "+", var_name)
  fit_linear <- .fit_model(response_lhs, linear_rhs, data, model_type, family)
  aic_linear <- if (!is.null(fit_linear)) tryCatch(stats::AIC(fit_linear), error = function(e) Inf) else Inf

  # --- Fit all non-linear forms and track the best ---
  # result_tbl: list of candidates with their fit info
  # We store each as a named list: transform, formula_rhs, df, label, deviance, powers, knots, fit
  candidates <- list()

  # Always add the linear candidate
  candidates[["linear"]] <- list(
    transform   = "linear",
    formula_rhs = var_name,
    df          = 1L,
    label       = var_name,
    deviance    = .model_deviance(fit_linear),
    powers      = NULL,
    knots       = NULL,
    fit         = fit_linear
  )

  for (tr in transforms) {
    if (tr == "linear") next  # already handled

    result <- tryCatch(
      .fit_one_transform(
        transform_type = tr,
        var_name       = var_name,
        base_rhs       = base_rhs,
        response_lhs   = response_lhs,
        data           = data,
        model_type     = model_type,
        family         = family,
        fp_powers      = fp_powers,
        spline_df      = spline_df
      ),
      error = function(e) NULL
    )

    if (!is.null(result)) {
      candidates[[tr]] <- result
    }
  }

  # --- Among non-linear forms, find the best (lowest deviance) ---
  nonlinear_keys <- setdiff(names(candidates), "linear")

  best_nonlinear <- NULL
  best_nonlinear_dev <- Inf

  for (k in nonlinear_keys) {
    dev_k <- candidates[[k]]$deviance
    if (is.finite(dev_k) && dev_k < best_nonlinear_dev) {
      best_nonlinear_dev <- dev_k
      best_nonlinear     <- candidates[[k]]
    }
  }

  # --- Compare best non-linear vs linear using LRT ---
  pvalue_vs_linear <- NA_real_
  use_nonlinear    <- FALSE

  if (!is.null(best_nonlinear) && !is.null(fit_linear)) {
    df_diff_nl <- best_nonlinear$df - 1L  # additional df over linear
    if (df_diff_nl > 0L) {
      # Fit the actual model for the best non-linear form to compute LRT properly
      fit_best_nl <- best_nonlinear$fit
      pvalue_vs_linear <- .lrt_pvalue(fit_best_nl, fit_linear, df_diff = df_diff_nl)
      if (!is.na(pvalue_vs_linear) && pvalue_vs_linear < linearity_alpha) {
        use_nonlinear <- TRUE
      }
    } else if (df_diff_nl == 0L) {
      # FP1 with 1 df: compare AIC instead (same df as linear, so LRT meaningless)
      aic_nl <- tryCatch(
        if (!is.null(best_nonlinear$fit)) stats::AIC(best_nonlinear$fit) else Inf,
        error = function(e) Inf
      )
      if (is.finite(aic_nl) && aic_nl < aic_linear) {
        use_nonlinear <- TRUE
      }
    }
  }

  # --- Select the winning form ---
  if (use_nonlinear && !is.null(best_nonlinear)) {
    winner <- best_nonlinear
  } else {
    winner <- candidates[["linear"]]
  }

  # --- LRT: winner vs null (current model without this variable) ---
  pvalue_vs_null <- .lrt_pvalue(winner$fit, fit_null, df_diff = winner$df)

  # --- AIC for the winner ---
  aic_best <- if (!is.null(winner$fit)) {
    tryCatch(stats::AIC(winner$fit), error = function(e) NA_real_)
  } else {
    NA_real_
  }

  list(
    var_name           = var_name,
    best_transform     = winner$transform,
    formula_rhs        = winner$formula_rhs,
    df                 = winner$df,
    label              = winner$label,
    pvalue_vs_null     = pvalue_vs_null,
    pvalue_vs_linear   = if (use_nonlinear) pvalue_vs_linear else NA_real_,
    aic_best           = aic_best,
    aic_linear         = aic_linear,
    powers             = winner$powers,
    knots              = winner$knots
  )
}


#' Fit a Single Non-Linear Transformation Form
#'
#' Internal helper called by \code{\link{.select_best_transform}} for each
#' non-linear transform type. Returns a candidate record.
#'
#' @param transform_type Character.
#' @param var_name Character.
#' @param base_rhs Character, current model RHS (no candidate variable).
#' @param response_lhs Character, LHS formula string.
#' @param data data.frame.
#' @param model_type Character.
#' @param family Family.
#' @param fp_powers Numeric vector.
#' @param spline_df Integer.
#'
#' @return A named list (transform, formula_rhs, df, label, deviance, powers,
#'   knots, fit), or \code{NULL} on failure.
#' @keywords internal
.fit_one_transform <- function(transform_type, var_name, base_rhs,
                                response_lhs, data,
                                model_type, family,
                                fp_powers, spline_df) {

  # Dispatch to the appropriate builder
  if (transform_type %in% c("log", "sqrt")) {
    info <- .make_simple_transform(var_name, transform_type)
    rhs  <- if (base_rhs == "1") info$formula_rhs else paste(base_rhs, "+", info$formula_rhs)
    fit  <- .fit_model(response_lhs, rhs, data, model_type, family)
    return(list(
      transform   = transform_type,
      formula_rhs = info$formula_rhs,
      df          = info$df,
      label       = info$label,
      deviance    = .model_deviance(fit),
      powers      = NULL,
      knots       = NULL,
      fit         = fit
    ))
  }

  if (transform_type == "fp1") {
    best <- .find_best_fp(
      var_name     = var_name,
      degree       = 1L,
      current_rhs  = base_rhs,
      data         = data,
      response_lhs = response_lhs,
      model_type   = model_type,
      family       = family,
      fp_powers    = fp_powers
    )
    if (is.null(best)) return(NULL)
    rhs <- if (base_rhs == "1") best$formula_rhs else paste(base_rhs, "+", best$formula_rhs)
    fit <- .fit_model(response_lhs, rhs, data, model_type, family)
    return(list(
      transform   = "fp1",
      formula_rhs = best$formula_rhs,
      df          = best$df,
      label       = best$label,
      deviance    = .model_deviance(fit),
      powers      = best$powers,
      knots       = NULL,
      fit         = fit
    ))
  }

  if (transform_type == "fp2") {
    best <- .find_best_fp(
      var_name     = var_name,
      degree       = 2L,
      current_rhs  = base_rhs,
      data         = data,
      response_lhs = response_lhs,
      model_type   = model_type,
      family       = family,
      fp_powers    = fp_powers
    )
    if (is.null(best)) return(NULL)
    rhs <- if (base_rhs == "1") best$formula_rhs else paste(base_rhs, "+", best$formula_rhs)
    fit <- .fit_model(response_lhs, rhs, data, model_type, family)
    return(list(
      transform   = "fp2",
      formula_rhs = best$formula_rhs,
      df          = best$df,
      label       = best$label,
      deviance    = .model_deviance(fit),
      powers      = best$powers,
      knots       = NULL,
      fit         = fit
    ))
  }

  if (grepl("^rcs", transform_type)) {
    df_rcs <- .rcs_df(transform_type, default_df = spline_df)
    info   <- .make_spline_transform(var_name, df = df_rcs, data = data)
    rhs    <- if (base_rhs == "1") info$formula_rhs else paste(base_rhs, "+", info$formula_rhs)
    fit    <- .fit_model(response_lhs, rhs, data, model_type, family)
    return(list(
      transform   = transform_type,
      formula_rhs = info$formula_rhs,
      df          = info$df,
      label       = info$label,
      deviance    = .model_deviance(fit),
      powers      = NULL,
      knots       = info$knots,
      fit         = fit
    ))
  }

  # Unknown type -- return NULL silently (validated upstream)
  NULL
}
