# Partial effect plots for non-linear terms in fitted SCM models.
# Works for both eres_scm_glm and eres_scm_coxph objects.

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Build a prediction grid for a single focal variable
#'
#' Creates a data frame where \code{var_name} varies over its observed range
#' and all other model variables are held at reference values (median for
#' continuous, reference/first level for factors).
#'
#' @param fit  eres_scm object.
#' @param var_name Character. The focal variable.
#' @param n_points Integer. Grid resolution.
#' @param ref_values Optional named list overriding default reference values.
#' @return data.frame ready for \code{predict()}.
#' @keywords internal
.build_prediction_grid <- function(fit, var_name, n_points = 100L,
                                    ref_values = NULL) {
  data  <- fit$model_data
  model <- fit$final_model

  if (is.null(data)) {
    cli::cli_abort("{.code fit$model_data} is NULL -- cannot build prediction grid.")
  }
  if (!var_name %in% colnames(data)) {
    cli::cli_abort("Variable {.val {var_name}} not found in {.code fit$model_data}.")
  }

  x_obs   <- data[[var_name]]
  x_range <- range(x_obs, na.rm = TRUE)
  x_seq   <- seq(x_range[1L], x_range[2L], length.out = n_points)

  # Collect predictor variable names from the model frame.
  # model.frame() returns named columns: the first is the response (DV or
  # Surv()), the rest are predictor expressions like "log(CMAX)", "AGE".
  # We need the raw *data* column names so we can look up medians / levels.
  other_vars <- tryCatch({
    mf        <- stats::model.frame(model)
    resp_name <- names(mf)[attr(stats::terms(model), "response")]
    pred_cols <- setdiff(names(mf), resp_name)
    # Strip any residual Surv() column (coxph)
    pred_cols <- pred_cols[!grepl("^Surv\\(|^survival::", pred_cols)]
    # Map formula expressions -> raw column names
    raw_names <- unique(vapply(pred_cols, .term_to_varname, character(1L),
                               USE.NAMES = FALSE))
    # Keep only those actually present in data, and exclude the focal variable
    intersect(setdiff(raw_names, var_name), colnames(data))
  }, error = function(e) {
    # Fallback: all data columns except response and focal variable
    resp <- if (inherits(model, "coxph")) c("DV", "TIME") else "DV"
    setdiff(colnames(data), c(resp, var_name))
  })

  # Build reference row (scalar per variable)
  ref_row <- stats::setNames(
    lapply(other_vars, function(v) {
      col <- data[[v]]
      if (!is.null(ref_values) && v %in% names(ref_values)) {
        ref_values[[v]]
      } else if (is.numeric(col)) {
        stats::median(col, na.rm = TRUE)
      } else if (is.factor(col)) {
        factor(levels(col)[1L], levels = levels(col))
      } else {
        col[!is.na(col)][1L]
      }
    }),
    other_vars
  )

  # Expand ref_row to n_points rows
  if (length(ref_row) > 0L) {
    grid <- as.data.frame(
      lapply(ref_row, function(v) rep(v, n_points)),
      stringsAsFactors = FALSE
    )
    # Restore factor levels (as.data.frame drops them on replication)
    for (v in other_vars) {
      if (is.factor(ref_row[[v]])) {
        grid[[v]] <- factor(grid[[v]], levels = levels(ref_row[[v]]))
      }
    }
  } else {
    # Intercept-only model or single-variable model: no other vars
    grid <- data.frame(row.names = seq_len(n_points))
  }

  # Set focal variable
  grid[[var_name]] <- x_seq
  grid
}


#' Compute partial effect predictions (link scale, centred at reference)
#'
#' @param fit     eres_scm object.
#' @param grid    Prediction grid from \code{.build_prediction_grid()}.
#' @param var_name Character. Focal variable name.
#' @param ci_level Numeric. Confidence level.
#' @return data.frame with columns \code{x_val}, \code{predicted},
#'   \code{ci_lower_nl}, \code{ci_upper_nl}, \code{se}.
#' @keywords internal
.predict_partial_effect <- function(fit, grid, var_name, ci_level = 0.95) {
  model    <- fit$final_model
  is_surv  <- inherits(model, "coxph")
  z        <- stats::qnorm((1 + ci_level) / 2)

  preds <- tryCatch({
    if (is_surv) {
      # coxph: use linear predictor + vcov for SE
      lp_res <- stats::predict(model, newdata = grid, type = "lp", se.fit = TRUE)
      list(fit = lp_res$fit, se.fit = lp_res$se.fit)
    } else {
      # glm: predict on link scale
      stats::predict(model, newdata = grid, type = "link", se.fit = TRUE)
    }
  }, error = function(e) {
    # Fallback: manual X %*% beta
    tryCatch({
      mm  <- stats::model.matrix(stats::terms(model), data = grid)
      cf  <- stats::coef(model)
      # Align columns
      common <- intersect(colnames(mm), names(cf))
      mm2    <- mm[, common, drop = FALSE]
      cf2    <- cf[common]
      lp_fit <- as.numeric(mm2 %*% cf2)
      vc     <- stats::vcov(model)[common, common, drop = FALSE]
      se_fit <- sqrt(pmax(0, rowSums((mm2 %*% vc) * mm2)))
      list(fit = lp_fit, se.fit = se_fit)
    }, error = function(e2) NULL)
  })

  if (is.null(preds)) {
    return(data.frame(
      x_val = grid[[var_name]], predicted = NA_real_,
      ci_lower_nl = NA_real_, ci_upper_nl = NA_real_, se = NA_real_
    ))
  }

  fit_vals <- preds$fit
  se_vals  <- preds$se.fit

  # Centre: subtract predicted value at the reference (median) row
  ref_idx  <- which.min(abs(grid[[var_name]] -
                              stats::median(fit$model_data[[var_name]], na.rm = TRUE)))
  ref_pred <- fit_vals[ref_idx]

  centred  <- fit_vals - ref_pred

  data.frame(
    x_val       = grid[[var_name]],
    predicted   = centred,
    ci_lower_nl = centred - z * se_vals,
    ci_upper_nl = centred + z * se_vals,
    se          = se_vals
  )
}


#' Subtitle annotation for the transformation type
#' @keywords internal
.transform_subtitle <- function(reg_entry, var_name) {
  if (is.null(reg_entry)) {
    return(paste("Variable:", var_name, "(linear)"))
  }
  tt <- reg_entry$transform_type
  base <- paste0("Variable: ", var_name, "  |  Transform: ", reg_entry$label)
  if (!is.null(reg_entry$pvalue_vs_linear) && !is.na(reg_entry$pvalue_vs_linear)) {
    pv <- reg_entry$pvalue_vs_linear
    pv_str <- if (pv < 0.001) "< 0.001" else sprintf("%.3f", pv)
    base <- paste0(base, "  (p vs linear = ", pv_str, ")")
  }
  base
}


# ---------------------------------------------------------------------------
# Main exported function
# ---------------------------------------------------------------------------

#' Plot Non-Linear Exposure-Response or Covariate-Response Curve
#'
#' Generates a partial effect plot for a variable in a fitted \code{eres_scm_glm}
#' or \code{eres_scm_coxph} model. The curve shows the predicted outcome on the
#' link scale (log-odds or log-hazard) or response scale, centred at the median
#' so the y-axis reads as "change from median".
#'
#' @param fit An object of class \code{"eres_scm_glm"} or \code{"eres_scm_coxph"}.
#' @param var_name Character. The variable to plot. Must appear in the final model.
#' @param scale Character. One of \code{"link"} (default, log-odds / log-HR) or
#'   \code{"response"} (probability for GLM, hazard ratio for CoxPH).
#' @param n_points Integer. Grid resolution for the prediction curve. Default 100.
#' @param ref_values Named list. Override reference values for other covariates.
#'   If \code{NULL} (default), medians / first factor levels are used.
#' @param ci_level Numeric. Confidence interval level. Default \code{0.95}.
#' @param rug Logical. Show a rug of observed values at the bottom. Default \code{TRUE}.
#' @param show_linear Logical. Overlay the linear fit as a dashed grey line.
#'   Default \code{FALSE}.
#' @param xlab Character or NULL. Custom x-axis label (overrides \code{var_spec}).
#' @param title Character or NULL. Custom plot title.
#' @param ... Not used; kept for forward compatibility.
#'
#' @return A \code{ggplot2} object (returned invisibly).
#'
#' @details
#' The partial effect is computed over the observed range of \code{var_name},
#' with all other variables fixed at reference values. Predictions are centred
#' by subtracting the value at the median of \code{var_name}, so the y-axis
#' represents the log-OR or log-HR \emph{relative to the median}.
#'
#' For \code{scale = "response"}: GLM predictions are inverse-logit transformed
#' (probabilities); CoxPH predictions are exponentiated (hazard ratios relative
#' to median = 1).
#'
#' @seealso \code{\link{plot_all_nonlinear}}, \code{\link{scm_glm}},
#'   \code{\link{scm_coxph}}
#'
#' @examples
#' data(demo_eres)
#' demo_eres$SEX <- as.factor(demo_eres$SEX)
#' fit <- scm_glm(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", covariates = c("AGE", "WT"),
#'   dv_col = "sim_status", significance_level = c(0.2, 0.2),
#'   nonlinear = list(CMAXss = c("linear", "fp1"), WT = c("linear", "fp1"), AGE = "linear")
#' )
#' plot_nonlinear_effect(fit, "CMAXss")
#' plot_nonlinear_effect(fit, "CMAXss", scale = "response", show_linear = TRUE)
#'
#' @export
plot_nonlinear_effect <- function(fit, var_name,
                                   scale        = c("link", "response"),
                                   n_points     = 100L,
                                   ref_values   = NULL,
                                   ci_level     = 0.95,
                                   rug          = TRUE,
                                   show_linear  = FALSE,
                                   xlab         = NULL,
                                   title        = NULL,
                                   ...) {

  # --- Input validation ---
  if (!inherits(fit, "eres_scm")) {
    cli::cli_abort("{.arg fit} must be an {.cls eres_scm_glm} or {.cls eres_scm_coxph} object.")
  }
  if (is.null(fit$final_model)) {
    cli::cli_abort(
      c(
        "Cannot plot: no model was fit.",
        "i" = "The SCM did not find a significant exposure variable.",
        "i" = "Check {.code fit$metadata$selected_exposure} for details."
      )
    )
  }
  if (!is.character(var_name) || length(var_name) != 1L) {
    cli::cli_abort("{.arg var_name} must be a single character string.")
  }
  scale <- match.arg(scale)

  meta   <- fit$metadata
  tr_reg <- meta$transform_registry

  # Check var_name is accessible: either in registry, or in finalRelation
  in_registry <- !is.null(tr_reg) && var_name %in% names(tr_reg)
  in_relation <- var_name %in% unlist(lapply(meta$finalRelation, .term_to_varname))
  in_data     <- !is.null(fit$model_data) && var_name %in% colnames(fit$model_data)

  if (!in_data) {
    cli::cli_abort("{.val {var_name}} is not found in {.code fit$model_data}.")
  }
  if (!in_registry && !in_relation) {
    cli::cli_warn("{.val {var_name}} does not appear in the final model -- plotting its partial effect anyway.")
  }

  reg_entry <- if (in_registry) tr_reg[[var_name]] else NULL
  is_surv   <- inherits(fit, "eres_scm_coxph")

  # --- Build prediction grid ---
  grid <- tryCatch(
    .build_prediction_grid(fit, var_name, n_points = n_points,
                            ref_values = ref_values),
    error = function(e) {
      cli::cli_abort("Failed to build prediction grid: {conditionMessage(e)}")
    }
  )

  # --- Compute partial effect ---
  pred_df <- .predict_partial_effect(fit, grid, var_name, ci_level = ci_level)

  # --- Optionally compute linear overlay ---
  linear_df <- NULL
  if (show_linear) {
    # Build a model with the variable in linear form
    lin_reg_entry <- list(transform_type = "linear", formula_rhs = var_name,
                          df = 1L, label = var_name, powers = NULL, knots = NULL)
    linear_df <- tryCatch({
      lin_fit  <- .fit_model(
        response_formula = if (is_surv) "survival::Surv(TIME, DV)" else "DV",
        rhs_terms        = paste(c(
          setdiff(meta$finalRelation,
                  if (!is.null(reg_entry)) reg_entry$formula_rhs else var_name),
          var_name
        ), collapse = " + "),
        data       = fit$model_data,
        model_type = if (is_surv) "coxph" else "glm",
        family     = if (!is_surv && !is.null(fit$final_model$family))
          fit$final_model$family$family else "binomial"
      )
      if (!is.null(lin_fit)) {
        tmp_fit <- fit
        tmp_fit$final_model <- lin_fit
        reg_copy <- tr_reg
        reg_copy[[var_name]] <- lin_reg_entry
        tmp_fit$metadata$transform_registry <- reg_copy
        tmp_fit$metadata$finalRelation <- setdiff(meta$finalRelation,
          if (!is.null(reg_entry)) reg_entry$formula_rhs else var_name)
        tmp_fit$metadata$finalRelation <- c(tmp_fit$metadata$finalRelation, var_name)
        .predict_partial_effect(tmp_fit, grid, var_name, ci_level = ci_level)
      } else NULL
    }, error = function(e) NULL)
  }

  # --- Transform to response scale if requested ---
  if (scale == "response") {
    if (is_surv) {
      pred_df$predicted   <- exp(pred_df$predicted)
      pred_df$ci_lower_nl <- exp(pred_df$ci_lower_nl)
      pred_df$ci_upper_nl <- exp(pred_df$ci_upper_nl)
    } else {
      # Response scale shows OR relative to the reference value (median), not absolute probability.
      pred_df$predicted   <- exp(pred_df$predicted)
      pred_df$ci_lower_nl <- exp(pred_df$ci_lower_nl)
      pred_df$ci_upper_nl <- exp(pred_df$ci_upper_nl)
    }
    if (!is.null(linear_df)) {
      linear_df$predicted   <- exp(linear_df$predicted)
      linear_df$ci_lower_nl <- exp(linear_df$ci_lower_nl)
      linear_df$ci_upper_nl <- exp(linear_df$ci_upper_nl)
    }
  }

  # --- Axis labels ---
  y_label <- if (scale == "link") {
    if (is_surv) "Log-HR (relative to median)" else "Log-OR (relative to median)"
  } else {
    if (is_surv) "HR (relative to median)" else "OR (relative to median)"
  }

  if (is.null(xlab)) {
    # Try var_spec first
    vs <- fit$var_spec
    if (!is.null(vs) && var_name %in% vs$col_name) {
      row_vs <- vs[vs$col_name == var_name, , drop = FALSE]
      lbl <- row_vs$label[1L]
      un  <- if ("units" %in% colnames(row_vs) && !is.na(row_vs$units[1L]))
        row_vs$units[1L] else NULL
      xlab <- if (!is.null(un)) paste0(lbl, " (", un, ")") else lbl
    } else {
      xlab <- var_name
    }
  }

  plot_title <- if (!is.null(title)) title else {
    if (is_surv) "Partial Effect (CoxPH)" else "Partial Effect (GLM)"
  }
  subtitle_txt <- .transform_subtitle(reg_entry, var_name)

  # Reference line value
  ref_line <- if (scale == "link") 0 else 1

  # --- Build ggplot ---
  p <- ggplot2::ggplot(pred_df,
                       ggplot2::aes(x = x_val, y = predicted)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_lower_nl, ymax = ci_upper_nl),
      fill  = "#4393C3",
      alpha = 0.20
    ) +
    ggplot2::geom_line(colour = "#2166AC", linewidth = 0.9) +
    ggplot2::geom_hline(yintercept = ref_line,
                        linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    ggplot2::labs(
      x        = xlab,
      y        = y_label,
      title    = plot_title,
      subtitle = subtitle_txt
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "grey30"),
      panel.grid.minor = ggplot2::element_blank()
    )

  # Rug plot
  if (rug) {
    obs_df <- data.frame(observed = fit$model_data[[var_name]])
    obs_df <- obs_df[!is.na(obs_df$observed), , drop = FALSE]
    p <- p + ggplot2::geom_rug(
      data   = obs_df,
      ggplot2::aes(x = observed),
      inherit.aes = FALSE,
      sides  = "b",
      alpha  = 0.25,
      colour = "grey50"
    )
  }

  # Linear overlay
  if (!is.null(linear_df)) {
    p <- p + ggplot2::geom_line(
      data    = linear_df,
      ggplot2::aes(x = x_val, y = predicted),
      inherit.aes = FALSE,
      linetype = "dashed",
      colour   = "grey40",
      linewidth = 0.7
    )
  }

  invisible(p)
}


# ---------------------------------------------------------------------------
# Plot all non-linear effects
# ---------------------------------------------------------------------------

#' Plot All Non-Linear Effects in a Fitted SCM Model
#'
#' Generates a partial effect plot for every variable in the model. When
#' \code{nonlinear} was used, variables in non-linear form are highlighted;
#' linear terms are still plotted but labelled accordingly.
#' Panels are arranged with \pkg{patchwork}.
#'
#' @param fit A fitted \code{eres_scm_glm} or \code{eres_scm_coxph} object.
#' @param nl_only Logical. If \code{TRUE} (default), only plot variables that
#'   entered in a non-linear form. If \code{FALSE}, plot all terms.
#' @param ncol Integer or NULL. Number of columns in the panel layout.
#'   \code{NULL} lets patchwork choose.
#' @param ... Further arguments passed to \code{\link{plot_nonlinear_effect}}.
#'
#' @return A \code{patchwork} object (or a single \code{ggplot} if only one
#'   variable), returned invisibly. If no plottable variables are found, a
#'   message is issued and \code{NULL} is returned invisibly.
#'
#' @seealso \code{\link{plot_nonlinear_effect}}
#'
#' @examples
#' data(demo_eres)
#' demo_eres$SEX <- as.factor(demo_eres$SEX)
#' fit <- scm_glm(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", covariates = c("AGE", "WT"),
#'   dv_col = "sim_status", significance_level = c(0.2, 0.2),
#'   nonlinear = list(CMAXss = c("linear", "fp1"), WT = c("linear", "fp1"), AGE = "linear")
#' )
#' plot_all_nonlinear(fit)
#' plot_all_nonlinear(fit, nl_only = FALSE)
#'
#' @export
plot_all_nonlinear <- function(fit, nl_only = TRUE, ncol = NULL, ...) {

  if (!inherits(fit, "eres_scm")) {
    cli::cli_abort("{.arg fit} must be an {.cls eres_scm_glm} or {.cls eres_scm_coxph} object.")
  }

  meta   <- fit$metadata
  tr_reg <- meta$transform_registry

  # Determine variables to plot
  if (!is.null(tr_reg) && length(tr_reg) > 0L) {
    if (nl_only) {
      plot_vars <- names(tr_reg)[vapply(tr_reg, function(r) {
        r$transform_type != "linear"
      }, logical(1L))]
    } else {
      plot_vars <- names(tr_reg)
    }
  } else {
    # No registry: fall back to finalRelation variable names
    if (nl_only) {
      message("No transform_registry found. Set nl_only = FALSE to plot all variables.")
      return(invisible(NULL))
    }
    plot_vars <- unique(vapply(meta$finalRelation, .term_to_varname,
                                character(1L), USE.NAMES = FALSE))
    # Keep only those present in model_data
    if (!is.null(fit$model_data)) {
      plot_vars <- plot_vars[plot_vars %in% colnames(fit$model_data)]
    }
  }

  # Exclude adjust_for from plotting (they are confounders, not of primary interest)
  if (!is.null(meta$adjust_for)) {
    plot_vars <- setdiff(plot_vars, meta$adjust_for)
  }

  if (length(plot_vars) == 0L) {
    message("No non-linear terms found. Use nl_only = FALSE to plot all terms.")
    return(invisible(NULL))
  }

  plots <- lapply(plot_vars, function(v) {
    tryCatch(
      plot_nonlinear_effect(fit, v, ...),
      error = function(e) {
        message("  Could not plot '", v, "': ", conditionMessage(e))
        NULL
      }
    )
  })
  plots <- Filter(Negate(is.null), plots)

  if (length(plots) == 0L) {
    message("No plots could be generated.")
    return(invisible(NULL))
  }
  if (length(plots) == 1L) {
    return(invisible(plots[[1L]]))
  }

  pw <- if (!is.null(ncol)) {
    patchwork::wrap_plots(plots, ncol = ncol)
  } else {
    patchwork::wrap_plots(plots)
  }
  invisible(pw)
}


# ---------------------------------------------------------------------------
# Quantile effects helper (for forest plot integration)
# ---------------------------------------------------------------------------

#' Compute Effects at Quantiles for a Non-Linear Term
#'
#' For a variable that entered the model in non-linear form, computes the
#' predicted effect (OR or HR) at specified quantiles relative to the median.
#'
#' @param fit  Fitted \code{eres_scm} object.
#' @param var_name Character. The variable name.
#' @param quantiles Numeric vector of probabilities. Default
#'   \code{c(0.10, 0.25, 0.75, 0.90)}.
#' @param transform_info Transform registry entry for this variable (optional;
#'   taken from \code{fit$metadata$transform_registry} if NULL).
#'
#' @return A tibble with columns \code{quantile}, \code{value}, \code{estimate},
#'   \code{ci_lower}, \code{ci_upper}.
#'
#' @keywords internal
.compute_quantile_effects <- function(fit, var_name,
                                       quantiles      = c(0.10, 0.25, 0.75, 0.90),
                                       transform_info = NULL) {
  if (!inherits(fit, "eres_scm")) {
    cli::cli_abort("{.arg fit} must be an {.cls eres_scm} object.")
  }

  data <- fit$model_data
  if (is.null(data) || !var_name %in% colnames(data)) {
    cli::cli_abort("{.val {var_name}} not found in {.code fit$model_data}.")
  }

  if (is.null(transform_info)) {
    tr_reg <- fit$metadata$transform_registry
    transform_info <- if (!is.null(tr_reg)) tr_reg[[var_name]] else NULL
  }

  x_obs  <- data[[var_name]][!is.na(data[[var_name]])]
  q_vals <- as.numeric(stats::quantile(x_obs, probs = quantiles))
  q_med  <- stats::median(x_obs)

  # Build a grid with the focal variable at each quantile value + the median
  all_x  <- c(q_med, q_vals)
  grid   <- .build_prediction_grid(fit, var_name, n_points = length(all_x))
  grid[[var_name]] <- all_x

  pred_df <- .predict_partial_effect(fit, grid, var_name, ci_level = 0.95)

  # Row 1 is the median (reference = 0 on link scale)
  is_surv <- inherits(fit, "eres_scm_coxph")

  rows <- lapply(seq_along(quantiles), function(i) {
    pr <- pred_df[i + 1L, ]
    est <- if (is_surv) exp(pr$predicted) else exp(pr$predicted)
    cil <- if (is_surv) exp(pr$ci_lower_nl) else exp(pr$ci_lower_nl)
    ciu <- if (is_surv) exp(pr$ci_upper_nl) else exp(pr$ci_upper_nl)
    tibble::tibble(
      quantile  = quantiles[i],
      value     = q_vals[i],
      estimate  = est,
      ci_lower  = cil,
      ci_upper  = ciu
    )
  })
  dplyr::bind_rows(rows)
}
