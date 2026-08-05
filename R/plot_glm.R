# ---- Internal helpers -------------------------------------------------------

#' Bin a continuous covariate into quantile groups and compute observed proportions
#'
#' @param data Data frame.
#' @param term Character; continuous covariate column name.
#' @param response_var Character; binary response column name (0/1).
#' @param n_bins Integer; number of quantile bins. Default \code{4L}.
#' @param conf_level Numeric; confidence level for exact binomial CIs. Default \code{0.95}.
#' @param group_var Character; optional column to stratify binning by. Default \code{NULL}.
#'
#' @return A data frame with columns: \code{bin_midpoint}, \code{n_events},
#'   \code{n_total}, \code{observed_prob}, \code{ci_lower}, \code{ci_upper},
#'   and optionally the grouping column.
#' @keywords internal
#' @noRd
bin_continuous_covariate <- function(data, term, response_var,
                                     n_bins = 4L, conf_level = 0.95,
                                     group_var = NULL) {
  if (!term %in% names(data)) {
    cli::cli_abort("Column {.val {term}} not found in {.arg data}.")
  }
  if (!response_var %in% names(data)) {
    cli::cli_abort("Column {.val {response_var}} not found in {.arg data}.")
  }

  .bin_single <- function(df) {
    breaks <- stats::quantile(df[[term]], probs = seq(0, 1, length.out = n_bins + 1L),
                              na.rm = TRUE)
    # Ensure unique breaks (handles degenerate distributions)
    breaks <- unique(breaks)
    if (length(breaks) < 2L) {
      return(data.frame(
        bin_midpoint = NA_real_, n_events = NA_integer_,
        n_total = nrow(df), observed_prob = NA_real_,
        ci_lower = NA_real_, ci_upper = NA_real_
      ))
    }
    df$bin <- cut(df[[term]], breaks = breaks, include.lowest = TRUE)
    result <- df[!is.na(df$bin), ]
    result <- split(result, result$bin)
    rows <- lapply(result, function(grp) {
      ne  <- sum(grp[[response_var]] == 1L, na.rm = TRUE)
      nt  <- nrow(grp)
      obs <- if (nt > 0L) ne / nt else NA_real_
      ci  <- if (nt > 0L && ne >= 0L) {
        stats::binom.test(ne, nt, conf.level = conf_level)$conf.int
      } else {
        c(NA_real_, NA_real_)
      }
      data.frame(
        bin_midpoint = mean(grp[[term]], na.rm = TRUE),
        n_events     = ne,
        n_total      = nt,
        observed_prob = obs,
        ci_lower     = ci[1L],
        ci_upper     = ci[2L]
      )
    })
    do.call(rbind, rows)
  }

  if (!is.null(group_var)) {
    groups <- unique(data[[group_var]])
    parts <- lapply(groups, function(g) {
      sub <- data[data[[group_var]] == g, , drop = FALSE]
      res <- .bin_single(sub)
      res[[group_var]] <- g
      res
    })
    return(do.call(rbind, parts))
  }

  .bin_single(data)
}


#' Generate prediction curve data for a single continuous covariate from a GLM
#'
#' @param model Fitted \code{glm} object.
#' @param term Character; covariate name to vary.
#' @param data Data frame used for model fitting (supplies variable ranges and
#'   reference values for other covariates).
#' @param n_points Integer; number of prediction points. Default \code{500L}.
#' @param conf_level Numeric; confidence level. Default \code{0.95}.
#'
#' @return Data frame with columns: \code{term_value}, \code{predicted_prob},
#'   \code{ci_lower}, \code{ci_upper}.
#' @keywords internal
#' @noRd
predict_glm_curve <- function(model, term, data, n_points = 500L, conf_level = 0.95) {
  z <- stats::qnorm((1 + conf_level) / 2)

  term_range <- range(data[[term]], na.rm = TRUE)
  newdata <- data.frame(
    x__ = seq(term_range[1L], term_range[2L], length.out = n_points)
  )
  names(newdata)[1L] <- term

  # All variables visible to the model (excludes response)
  mf_vars <- setdiff(names(stats::model.frame(model)), names(model$model)[1L])
  other_vars <- setdiff(mf_vars, term)

  for (v in other_vars) {
    col <- data[[v]]
    if (is.null(col)) next
    if (is.numeric(col)) {
      newdata[[v]] <- stats::median(col, na.rm = TRUE)
    } else {
      tab <- sort(table(col), decreasing = TRUE)
      mode_val <- names(tab)[1L]
      newdata[[v]] <- if (is.factor(col)) {
        factor(mode_val, levels = levels(col))
      } else {
        mode_val
      }
    }
  }

  pred <- stats::predict(model, newdata = newdata, type = "link", se.fit = TRUE)

  data.frame(
    term_value    = newdata[[term]],
    predicted_prob = stats::plogis(pred$fit),
    ci_lower      = stats::plogis(pred$fit - z * pred$se.fit),
    ci_upper      = stats::plogis(pred$fit + z * pred$se.fit)
  )
}


#' Generate predictions for each level of a categorical covariate
#'
#' @param model Fitted \code{glm} object.
#' @param term Character; categorical covariate name.
#' @param data Data frame.
#' @param metric Character; \code{"probability"} or \code{"odds_ratio"}.
#' @param conf_level Numeric; confidence level. Default \code{0.95}.
#'
#' @return Data frame with columns: \code{level}, \code{estimate},
#'   \code{ci_lower}, \code{ci_upper}, \code{n_obs}, \code{n_events}.
#' @keywords internal
#' @noRd
predict_glm_categorical <- function(model, term, data,
                                    metric = "probability",
                                    conf_level = 0.95) {
  z <- stats::qnorm((1 + conf_level) / 2)
  col_data   <- data[[term]]
  cat_levels <- if (is.factor(col_data)) levels(col_data) else sort(unique(col_data))

  # Variables visible to the model minus response and the focal term
  mf_vars    <- setdiff(names(stats::model.frame(model)), names(model$model)[1L])
  other_vars <- setdiff(mf_vars, term)

  # Reference values: median for numeric, mode for categorical
  ref_vals <- list()
  for (v in other_vars) {
    col <- data[[v]]
    if (is.null(col)) next
    if (is.numeric(col)) {
      ref_vals[[v]] <- stats::median(col, na.rm = TRUE)
    } else {
      tab <- sort(table(col), decreasing = TRUE)
      mode_val <- names(tab)[1L]
      ref_vals[[v]] <- if (is.factor(col)) factor(mode_val, levels = levels(col)) else mode_val
    }
  }

  if (metric == "probability") {
    rows <- lapply(cat_levels, function(lev) {
      row <- data.frame(x__ = lev, stringsAsFactors = FALSE)
      names(row)[1L] <- term
      if (is.factor(col_data)) row[[term]] <- factor(lev, levels = levels(col_data))

      for (v in other_vars) {
        if (!is.null(ref_vals[[v]])) row[[v]] <- ref_vals[[v]]
      }

      lp <- stats::predict(model, newdata = row, type = "link", se.fit = TRUE)
      data.frame(
        level    = as.character(lev),
        estimate = stats::plogis(lp$fit),
        ci_lower = stats::plogis(lp$fit - z * lp$se.fit),
        ci_upper = stats::plogis(lp$fit + z * lp$se.fit)
      )
    })
  } else {
    # Odds ratio: extract from model coefficients
    coef_mat  <- summary(model)$coefficients
    ci_mat    <- tryCatch(
      stats::confint.default(model, level = conf_level),
      error = function(e) NULL
    )

    rows <- lapply(cat_levels, function(lev) {
      # Reference level: OR = 1
      coef_name <- paste0(term, lev)
      if (!coef_name %in% rownames(coef_mat)) {
        return(data.frame(level = as.character(lev), estimate = 1,
                          ci_lower = 1, ci_upper = 1))
      }
      beta <- coef_mat[coef_name, 1L]
      if (!is.null(ci_mat) && coef_name %in% rownames(ci_mat)) {
        lo <- ci_mat[coef_name, 1L]
        hi <- ci_mat[coef_name, 2L]
      } else {
        se <- coef_mat[coef_name, 2L]
        lo <- beta - z * se
        hi <- beta + z * se
      }
      data.frame(
        level    = as.character(lev),
        estimate = exp(beta),
        ci_lower = exp(lo),
        ci_upper = exp(hi)
      )
    })
  }

  result <- do.call(rbind, rows)

  # Merge observed counts
  resp_col <- names(model$model)[1L]
  dv <- model$model[[resp_col]]
  grp <- as.character(col_data)
  obs <- tapply(dv == 1L, grp, sum, na.rm = TRUE)
  tot <- tapply(rep(1L, length(grp)), grp, sum)
  result$n_events <- as.integer(obs[result$level])
  result$n_obs    <- as.integer(tot[result$level])

  result
}


# ---- Exported lower-level functions -----------------------------------------

#' Plot predicted probability curve for a continuous covariate (GLM)
#'
#' For a single continuous covariate from a binomial GLM, produces:
#' \itemize{
#'   \item A predicted probability curve with confidence ribbon.
#'   \item Binned observed event proportions with exact binomial CIs (pointrange).
#'   \item Optional rug plots: events on top, non-events on bottom.
#'   \item Optional vertical reference lines at the median and 5th/95th percentiles.
#'   \item Optional p-value and N annotations.
#' }
#' Predictions beyond \code{extrapolation_pct} of the covariate's observed
#' range are shown with a dashed line (or hidden) to signal extrapolation.
#'
#' @param model A fitted \code{glm} object with \code{family = binomial}.
#' @param data Data frame containing the original observations.
#' @param term Character string; name of the continuous covariate column.
#' @param response_var Character string; name of the binary response column.
#' @param n_bins Integer; number of quantile bins for observed data overlay.
#'   Default \code{4L}.
#' @param n_pred Integer; number of points for the prediction curve.
#'   Default \code{500L}.
#' @param conf_level Numeric; confidence level for all intervals. Default \code{0.95}.
#' @param show_rug Logical; show rug tick marks. Default \code{TRUE}.
#' @param show_annotations Logical; annotate with p-value and N. Default \code{TRUE}.
#' @param show_reference_lines Logical; show vertical lines at the median and
#'   5th/95th percentiles. Default \code{FALSE}.
#' @param extrapolation_pct Numeric in \code{(0, 1]}; percentile beyond which
#'   predictions are marked as extrapolation. Default \code{0.99}.
#' @param extrapolation_style Character; \code{"dashed"} (default) or \code{"hide"}.
#' @param colour Character; colour for the prediction line and ribbon.
#'   Default \code{"steelblue"}.
#' @param fill_alpha Numeric; transparency for the confidence ribbon.
#'   Default \code{0.15}.
#' @param point_colour Character; colour for binned observed data points.
#'   Default \code{"grey30"}.
#' @param point_size Numeric; size of the observed data pointrange markers.
#'   Default \code{3}.
#' @param xlab Character; x-axis label. Default is \code{term}.
#' @param ylab Character; y-axis label. Default \code{"Predicted Probability"}.
#' @param title Character or \code{NULL}; plot title. Default \code{NULL}.
#' @param rug_colour_var Character or \code{NULL}; optional column name to
#'   colour rug marks by (e.g., a dose group). Default \code{NULL}.
#'
#' @return A \code{ggplot} object.
#' @seealso \code{\link{plot_glm_categorical}}, \code{\link{plot_glm_terms}}
#' @export
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(200), y = rbinom(200, 1, 0.4))
#' m <- stats::glm(y ~ x, family = stats::binomial, data = d)
#' plot_glm_continuous(m, d, "x", "y")
plot_glm_continuous <- function(model, data, term, response_var,
                                 n_bins = 4L, n_pred = 500L,
                                 conf_level = 0.95,
                                 show_rug = TRUE,
                                 show_annotations = TRUE,
                                 show_reference_lines = FALSE,
                                 extrapolation_pct = 0.99,
                                 extrapolation_style = c("dashed", "hide"),
                                 colour = "steelblue",
                                 fill_alpha = 0.15,
                                 point_colour = "grey30",
                                 point_size = 3,
                                 xlab = term,
                                 ylab = "Predicted Probability",
                                 title = NULL,
                                 rug_colour_var = NULL) {

  extrapolation_style <- match.arg(extrapolation_style)

  pred_data <- predict_glm_curve(model, term, data, n_points = n_pred,
                                  conf_level = conf_level)
  bin_data  <- bin_continuous_covariate(data, term, response_var,
                                         n_bins = n_bins, conf_level = conf_level)

  extrap_boundary <- stats::quantile(data[[term]], probs = extrapolation_pct, na.rm = TRUE)
  pred_interp <- pred_data[pred_data$term_value <= extrap_boundary, ]
  pred_extrap <- pred_data[pred_data$term_value > extrap_boundary, ]

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = pred_interp,
      ggplot2::aes(x = .data$term_value, ymin = .data$ci_lower, ymax = .data$ci_upper),
      fill = colour, alpha = fill_alpha
    ) +
    ggplot2::geom_line(
      data = pred_interp,
      ggplot2::aes(x = .data$term_value, y = .data$predicted_prob),
      colour = colour, linewidth = 1
    )

  if (nrow(pred_extrap) > 0L && extrapolation_style == "dashed") {
    p <- p +
      ggplot2::geom_line(
        data = pred_extrap,
        ggplot2::aes(x = .data$term_value, y = .data$predicted_prob),
        colour = colour, linewidth = 1, linetype = "dashed"
      )
  }

  p <- p +
    ggplot2::geom_pointrange(
      data = bin_data,
      ggplot2::aes(x = .data$bin_midpoint, y = .data$observed_prob,
                   ymin = .data$ci_lower, ymax = .data$ci_upper),
      colour = point_colour, size = point_size / 3, linewidth = 0.8,
      inherit.aes = FALSE
    )

  if (show_rug) {
    events     <- data[data[[response_var]] == 1L, , drop = FALSE]
    non_events <- data[data[[response_var]] == 0L, , drop = FALSE]

    if (!is.null(rug_colour_var)) {
      p <- p +
        ggplot2::geom_rug(
          data = events,
          ggplot2::aes(x = .data[[term]], colour = .data[[rug_colour_var]]),
          sides = "t", alpha = 0.5, inherit.aes = FALSE
        ) +
        ggplot2::geom_rug(
          data = non_events,
          ggplot2::aes(x = .data[[term]], colour = .data[[rug_colour_var]]),
          sides = "b", alpha = 0.5, inherit.aes = FALSE
        )
    } else {
      p <- p +
        ggplot2::geom_rug(
          data = events,
          ggplot2::aes(x = .data[[term]]),
          sides = "t", alpha = 0.3, colour = "grey40", inherit.aes = FALSE
        ) +
        ggplot2::geom_rug(
          data = non_events,
          ggplot2::aes(x = .data[[term]]),
          sides = "b", alpha = 0.3, colour = "grey60", inherit.aes = FALSE
        )
    }
  }

  if (show_reference_lines) {
    med <- stats::median(data[[term]], na.rm = TRUE)
    pct <- stats::quantile(data[[term]], probs = c(0.05, 0.95), na.rm = TRUE)
    p <- p +
      ggplot2::geom_vline(xintercept = med, linetype = "solid",
                          colour = "grey50", linewidth = 0.5) +
      ggplot2::geom_vline(xintercept = pct, linetype = "dashed",
                          colour = "grey50", linewidth = 0.4)
  }

  if (show_annotations) {
    coef_tbl <- tryCatch(summary(model)$coefficients, error = function(e) NULL)
    if (!is.null(coef_tbl)) {
      match_idx <- grep(paste0("^", term), rownames(coef_tbl))
      if (length(match_idx) > 0L) {
        pval <- coef_tbl[match_idx[1L], "Pr(>|z|)"]
        pval_label <- if (pval < 0.001) "p < 0.001" else
          paste0("p = ", formatC(pval, format = "f", digits = 3))
        n_label <- paste0("N = ", nrow(data))
        x_pos <- stats::quantile(data[[term]], 0.95, na.rm = TRUE)
        p <- p +
          ggplot2::annotate("text", x = x_pos, y = 0.95, label = pval_label,
                            hjust = 1, size = 3.5) +
          ggplot2::annotate("text", x = x_pos, y = 0.87, label = n_label,
                            hjust = 1, size = 3.5)
      }
    }
  }

  p +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::labs(x = xlab, y = ylab, title = title) +
    ggplot2::theme_bw() +
    ggpubr::labs_pubr()
}


#' Plot predicted probability or odds ratio by category (GLM)
#'
#' For a single categorical covariate from a binomial GLM, produces a
#' pointrange plot showing either predicted probabilities or odds ratios for
#' each category level, optionally overlaid with observed event proportions.
#'
#' @param model A fitted \code{glm} object with \code{family = binomial}.
#' @param data Data frame containing the original observations.
#' @param term Character string; name of the categorical covariate column.
#' @param response_var Character string; name of the binary response column.
#' @param metric Character; \code{"probability"} (default) or
#'   \code{"odds_ratio"}. \code{"probability"} shows predicted P(event) for
#'   each level (other covariates held at reference/median).
#'   \code{"odds_ratio"} shows OR relative to the model's reference level.
#' @param conf_level Numeric; confidence level for all intervals. Default \code{0.95}.
#' @param show_observed Logical; overlay observed proportions per level.
#'   Default \code{TRUE}.
#' @param show_n Logical; annotate each level with its sample size.
#'   Default \code{TRUE}.
#' @param colour Character; colour for model estimate markers.
#'   Default \code{"steelblue"}.
#' @param observed_colour Character; colour for observed proportion markers.
#'   Default \code{"grey50"}.
#' @param xlab Character or \code{NULL}; x-axis label.
#' @param ylab Character or \code{NULL}; y-axis label. Default is \code{term}.
#' @param title Character or \code{NULL}; plot title. Default \code{NULL}.
#'
#' @return A \code{ggplot} object.
#' @seealso \code{\link{plot_glm_continuous}}, \code{\link{plot_glm_terms}}
#' @export
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(g = factor(sample(c("A","B","C"), 200, replace = TRUE)),
#'                 y = rbinom(200, 1, 0.3))
#' m <- stats::glm(y ~ g, family = stats::binomial, data = d)
#' plot_glm_categorical(m, d, "g", "y")
plot_glm_categorical <- function(model, data, term, response_var,
                                  metric = c("probability", "odds_ratio"),
                                  conf_level = 0.95,
                                  show_observed = TRUE,
                                  show_n = TRUE,
                                  colour = "steelblue",
                                  observed_colour = "grey50",
                                  xlab = NULL,
                                  ylab = term,
                                  title = NULL) {

  metric <- match.arg(metric)
  pred_data <- predict_glm_categorical(model, term, data,
                                        metric = metric,
                                        conf_level = conf_level)

  if (metric == "probability") {
    if (is.null(xlab)) xlab <- "Predicted Probability"

    obs_data <- data.frame(
      level = as.character(data[[term]]),
      dv    = data[[response_var]]
    )
    obs_sum <- stats::aggregate(dv ~ level, data = obs_data,
                                FUN = function(x) c(mean = mean(x == 1L), n = length(x)))
    obs_df <- data.frame(
      level         = obs_sum$level,
      observed_prob = obs_sum$dv[, "mean"],
      n_obs         = as.integer(obs_sum$dv[, "n"])
    )

    plot_data <- merge(pred_data, obs_df, by = "level", all.x = TRUE)

    p <- ggplot2::ggplot(plot_data,
                         ggplot2::aes(x = .data$estimate, y = .data$level)) +
      ggplot2::geom_pointrange(
        ggplot2::aes(xmin = .data$ci_lower, xmax = .data$ci_upper),
        colour = colour, size = 0.8, linewidth = 0.7
      ) +
      ggplot2::coord_cartesian(xlim = c(0, 1)) +
      ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
      ggplot2::theme_bw() +
      ggpubr::labs_pubr()

    if (show_observed) {
      p <- p +
        ggplot2::geom_point(
          ggplot2::aes(x = .data$observed_prob),
          colour = observed_colour, shape = 4L, size = 3
        )
    }

    if (show_n) {
      x_right <- max(plot_data$ci_upper, na.rm = TRUE) + 0.03
      p <- p +
        ggplot2::geom_text(
          ggplot2::aes(x = x_right, label = paste0("n=", .data$n_obs.x)),
          hjust = 0, size = 3
        )
    }

  } else {
    # Odds ratio
    if (is.null(xlab)) xlab <- "Odds Ratio"

    p <- ggplot2::ggplot(pred_data,
                         ggplot2::aes(x = .data$estimate, y = .data$level)) +
      ggplot2::geom_pointrange(
        ggplot2::aes(xmin = .data$ci_lower, xmax = .data$ci_upper),
        colour = colour, size = 0.8, linewidth = 0.7
      ) +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
      ggplot2::scale_x_log10() +
      ggplot2::theme_bw() +
      ggpubr::labs_pubr()

    if (show_n) {
      x_right <- max(pred_data$ci_upper, na.rm = TRUE) * 1.1
      p <- p +
        ggplot2::geom_text(
          ggplot2::aes(x = x_right, label = paste0("n=", .data$n_obs)),
          hjust = 0, size = 3
        )
    }
  }

  p + ggplot2::labs(x = xlab, y = ylab, title = title)
}


# ---- High-level wrapper -----------------------------------------------------

#' Plot covariate effects from a fitted GLM or SCM result
#'
#' Produces diagnostic/presentation plots for each covariate in a binomial GLM
#' or in the final model of a \code{\link{scm_glm}()} result. Continuous
#' covariates receive predicted probability curves with binned observed data
#' overlay; categorical covariates receive pointrange plots showing predicted
#' probability (or odds ratio) by category level.
#'
#' @param scm_result A \code{\link{scm_glm}()} result object
#'   (class \code{eres_scm_glm}) \emph{or} a fitted \code{glm} object.
#' @param data Data frame used for model fitting. Required when
#'   \code{scm_result} is a plain \code{glm}. When \code{scm_result} is an
#'   \code{eres_scm_glm} object, \code{model_data} is extracted automatically
#'   and this argument is optional.
#' @param terms Character vector of covariate names to plot. If \code{NULL}
#'   (default), all model terms are plotted (excluding interaction terms and
#'   the intercept).
#' @param type Character; how to classify covariate types. One of
#'   \code{"auto"} (default; numeric columns become continuous, all others
#'   become categorical), \code{"continuous"}, or \code{"categorical"} (force
#'   all terms to one type).
#' @param response_var Character string naming the binary response column.
#'   Extracted automatically from the model when \code{NULL} (default).
#' @param n_bins Integer; number of quantile bins for the observed data overlay
#'   on continuous plots. Default \code{4L}.
#' @param conf_level Numeric; confidence level for all intervals. Default \code{0.95}.
#' @param show_rug Logical; show rug tick marks on continuous plots.
#'   Default \code{TRUE}.
#' @param show_annotations Logical; annotate continuous plots with p-value and
#'   N. Default \code{TRUE}.
#' @param show_reference_lines Logical; show vertical lines at the median and
#'   5th/95th percentiles for continuous covariates. Default \code{FALSE}.
#' @param extrapolation_style Character; \code{"dashed"} (default) or
#'   \code{"hide"} for predictions beyond the 99th percentile of observed data.
#' @param colour Character; primary colour for prediction lines/ribbons.
#'   Default \code{"steelblue"}.
#' @param fill_alpha Numeric; transparency for confidence ribbons.
#'   Default \code{0.15}.
#' @param point_colour Character; colour for binned observed data markers.
#'   Default \code{"grey30"}.
#' @param combine Logical; if \code{TRUE} (default), combine individual plots
#'   into a multi-panel figure using \pkg{patchwork}. If \code{FALSE}, return
#'   a named list of individual \code{ggplot} objects.
#' @param ncol Integer or \code{NULL}; number of columns when combining plots.
#'   \code{NULL} (default) uses \code{min(2, number_of_terms)}.
#' @param adjust_for Character vector of terms that are fixed adjustments.
#'   Extracted automatically from \code{scm_result$metadata$adjust_for} when
#'   \code{NULL} (default). Fixed-adjustment panels receive a
#'   a "Fixed adjustment" subtitle annotation (dagger symbol).
#' @param ... Additional arguments passed to \code{\link{plot_glm_continuous}}
#'   (for continuous terms) or \code{\link{plot_glm_categorical}} (for
#'   categorical terms).
#'
#' @return If \code{combine = TRUE} and there are multiple terms, a
#'   \pkg{patchwork} object. If \code{combine = TRUE} and there is one term,
#'   a single \code{ggplot} object. If \code{combine = FALSE}, a named list
#'   of \code{ggplot} objects.
#'
#' @seealso \code{\link{plot_glm_continuous}}, \code{\link{plot_glm_categorical}},
#'   \code{\link{scm_glm}}
#'
#' @export
#'
#' @examples
#' data(demo_eres)
#' demo_eres$SEX <- as.factor(demo_eres$SEX)
#' result <- scm_glm(
#'   eres_data = demo_eres,
#'   covariates = c("AGE", "WT"),
#'   exposure_vars = "CMAXss",
#'   dv_col = "sim_status",
#'   significance_level = c(0.2, 0.2)
#' )
#' plot_glm_terms(result, terms = "WT")
#' plot_glm_terms(result, terms = "WT", combine = FALSE)
plot_glm_terms <- function(scm_result, data = NULL,
                            terms = NULL,
                            type = c("auto", "continuous", "categorical"),
                            response_var = NULL,
                            n_bins = 4L,
                            conf_level = 0.95,
                            show_rug = TRUE,
                            show_annotations = TRUE,
                            show_reference_lines = FALSE,
                            extrapolation_style = c("dashed", "hide"),
                            colour = "steelblue",
                            fill_alpha = 0.15,
                            point_colour = "grey30",
                            combine = TRUE,
                            ncol = NULL,
                            adjust_for = NULL,
                            ...) {

  type               <- match.arg(type)
  extrapolation_style <- match.arg(extrapolation_style)

  # ---- Unpack scm_result or plain glm ----
  if (inherits(scm_result, "eres_scm_glm")) {
    model <- scm_result$final_model
    if (is.null(data)) data <- scm_result$model_data
    if (is.null(response_var)) {
      response_var <- as.character(stats::formula(model)[[2L]])
    }
    if (is.null(adjust_for)) adjust_for <- scm_result$metadata$adjust_for
  } else if (inherits(scm_result, "glm")) {
    model <- scm_result
    if (is.null(data)) {
      cli::cli_abort("{.arg data} must be supplied when {.arg scm_result} is a {.cls glm} object.")
    }
    if (is.null(response_var)) {
      response_var <- as.character(stats::formula(model)[[2L]])
    }
  } else {
    cli::cli_abort("{.arg scm_result} must be a {.cls glm} object or an {.cls eres_scm_glm} result.")
  }

  if (is.null(data)) {
    cli::cli_abort("{.arg data} could not be determined automatically. Please supply it explicitly.")
  }

  # ---- Determine terms ----
  if (is.null(terms)) {
    terms <- attr(stats::terms(model), "term.labels")
  }

  # Skip interaction terms with a message
  interaction_mask <- grepl(":", terms, fixed = TRUE)
  if (any(interaction_mask)) {
    cli::cli_inform(c(
      "!" = "Interaction terms are not plotted individually; skipping:",
      " " = "{.val {terms[interaction_mask]}}"
    ))
    terms <- terms[!interaction_mask]
  }

  if (length(terms) == 0L) {
    cli::cli_abort("No plottable terms found.")
  }

  # ---- Classify each term ----
  term_types <- vapply(terms, function(t) {
    if (type != "auto") return(type)
    if (t %in% names(data) && is.numeric(data[[t]])) "continuous" else "categorical"
  }, character(1L))

  # ---- Build individual plots ----
  plots <- vector("list", length(terms))
  names(plots) <- terms

  for (i in seq_along(terms)) {
    t  <- terms[i]
    tt <- term_types[i]

    if (tt == "continuous") {
      p <- plot_glm_continuous(
        model = model, data = data, term = t,
        response_var = response_var,
        n_bins = n_bins, conf_level = conf_level,
        show_rug = show_rug,
        show_annotations = show_annotations,
        show_reference_lines = show_reference_lines,
        extrapolation_style = extrapolation_style,
        colour = colour, fill_alpha = fill_alpha,
        point_colour = point_colour,
        ...
      )
    } else {
      p <- plot_glm_categorical(
        model = model, data = data, term = t,
        response_var = response_var,
        conf_level = conf_level,
        colour = colour,
        ...
      )
    }

    if (!is.null(adjust_for) && t %in% adjust_for) {
      p <- p + ggplot2::labs(subtitle = "\u2020 Fixed adjustment")
    }

    plots[[i]] <- p
  }

  # ---- Combine or return ----
  if (combine && length(plots) > 1L) {
    if (is.null(ncol)) ncol <- min(2L, length(plots))
    return(patchwork::wrap_plots(plots, ncol = ncol))
  }

  if (combine && length(plots) == 1L) return(plots[[1L]])

  plots
}
