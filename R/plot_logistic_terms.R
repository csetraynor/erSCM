#' Plot covariate effects for a fitted logistic regression model
#'
#' Produces a covariate effect plot for a single predictor from a binomial GLM
#' fit object returned by \code{scm_glm()} (or compatible). For continuous
#' covariates the plot shows a smooth predicted probability curve with confidence
#' ribbon, a rug of observed events/non-events, and optionally binned observed
#' event rates. For categorical covariates it shows predicted vs. observed
#' probabilities per level with dodged points and error bars.
#'
#' @param fit A list from \code{scm_glm()} containing \code{fit$final_model} (a
#'   \code{glm} object), \code{fit$model_data}, and \code{fit$metadata}.
#' @param covariate_term Character; name of the covariate to plot. Must be a
#'   main-effect variable -- interaction terms (containing \code{:}) are not
#'   supported.
#' @param ci Numeric; confidence level for intervals. Default \code{0.95}.
#' @param n_bins Integer; number of quantile bins used to compute binned
#'   observed event rates for continuous covariates. Default \code{4}.
#' @param show_rug Logical; show rug ticks for events (top) and non-events
#'   (bottom) on continuous covariate plots. Default \code{TRUE}.
#' @param show_observed Logical; overlay binned observed event rates
#'   (continuous) or per-level observed rates (categorical). Default \code{TRUE}.
#' @param plot Logical; if \code{TRUE} (default) return a \code{ggplot} object,
#'   if \code{FALSE} return a list with \code{predictions} and \code{observed}
#'   tibbles.
#' @param n_grid Integer; number of points in the prediction grid for continuous
#'   covariates. Default \code{300}.
#'
#' @return If \code{plot = TRUE}, a \code{ggplot2} object. If \code{plot = FALSE},
#'   a named list with elements \code{predictions} and \code{observed}.
#'
#' @seealso \code{\link{plot_cox_terms}}, \code{\link{plot_model_terms}}
#'
#' @examples
#' data(demo_eres)
#' demo_eres$SEX <- as.factor(demo_eres$SEX)
#' lr_fit <- scm_glm(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", covariates = c("AGE", "WT"),
#'   dv_col = "sim_status", significance_level = c(0.2, 0.2)
#' )
#' plot_logistic_terms(lr_fit, "WT")
#' plot_logistic_terms(lr_fit, "WT", plot = FALSE)
#'
#' @export
plot_logistic_terms <- function(fit, covariate_term, ci = 0.95, n_bins = 4,
                                show_rug = TRUE, show_observed = TRUE,
                                plot = TRUE, n_grid = 300) {

  if (grepl(":", covariate_term, fixed = TRUE)) {
    main_effects <- get_main_effect_columns(fit$metadata$finalRelation)
    cli::cli_abort(c(
      "{.fun plot_logistic_terms} does not support interaction terms as {.arg covariate_term}.",
      "i" = "Use a main-effect variable name, e.g. one of: {.val {main_effects}}."
    ))
  }

  z <- stats::qnorm(1 - (1 - ci) / 2)

  main_effect_cols <- get_main_effect_columns(fit$metadata$finalRelation)

  ori_data <- fit$model_data %>%
    dplyr::select(dplyr::all_of(main_effect_cols), DV)

  tmp_data <- ori_data %>% dplyr::select(-DV)

  finalCharacters <- dplyr::bind_cols(
    dplyr::select(tmp_data, dplyr::where(is.factor)),
    dplyr::select(tmp_data, dplyr::where(is.character))
  ) %>%
    dplyr::bind_cols(dplyr::select(tmp_data, dplyr::where(is.logical))) %>%
    colnames()

  finalNumeric <- dplyr::select(tmp_data, dplyr::where(is.numeric)) %>% colnames()

  is_covariate_numeric <- covariate_term %in% finalNumeric

  other_cols <- setdiff(main_effect_cols, covariate_term)

  # Reference values: median for numeric, mode for categorical
  ref_vals <- list()
  for (col in other_cols) {
    if (col %in% finalNumeric) {
      ref_vals[[col]] <- stats::median(tmp_data[[col]], na.rm = TRUE)
    } else {
      tbl <- sort(table(tmp_data[[col]]), decreasing = TRUE)
      ref_vals[[col]] <- names(tbl)[1]
    }
  }

  dv_label <- if (!is.null(fit$metadata$dvid_name)) {
    fit$metadata$dvid_name
  } else if (!is.null(fit$metadata$dvid)) {
    fit$metadata$dvid
  } else {
    "Event"
  }

  # Resolve display label and units from var_spec
  var_spec_obj  <- fit$var_spec
  display_label <- get_var_label(var_spec_obj, covariate_term)
  display_units <- get_var_units(var_spec_obj, covariate_term)
  x_label <- if (!is.na(display_units) && !is.null(display_units)) {
    paste0(display_label, " (", display_units, ")")
  } else {
    display_label
  }

  if (is_covariate_numeric) {
    .plot_logistic_continuous(
      fit = fit, covariate_term = covariate_term,
      ori_data = ori_data, tmp_data = tmp_data,
      finalCharacters = finalCharacters, finalNumeric = finalNumeric,
      other_cols = other_cols, ref_vals = ref_vals,
      z = z, ci = ci, n_bins = n_bins,
      show_rug = show_rug, show_observed = show_observed,
      plot = plot, n_grid = n_grid, dv_label = dv_label,
      x_label = x_label, var_spec_obj = var_spec_obj
    )
  } else {
    # Resolve category display labels for the x-axis
    cat_labels <- get_var_categories(var_spec_obj, covariate_term)
    .plot_logistic_categorical(
      fit = fit, covariate_term = covariate_term,
      ori_data = ori_data, tmp_data = tmp_data,
      finalCharacters = finalCharacters,
      other_cols = other_cols, ref_vals = ref_vals,
      z = z, plot = plot, dv_label = dv_label,
      x_label = x_label, cat_labels = cat_labels
    )
  }
}


# Internal: continuous covariate path
#' @keywords internal
#' @noRd
.plot_logistic_continuous <- function(fit, covariate_term, ori_data, tmp_data,
                                      finalCharacters, finalNumeric,
                                      other_cols, ref_vals,
                                      z, ci, n_bins, show_rug, show_observed,
                                      plot, n_grid, dv_label,
                                      x_label = covariate_term,
                                      var_spec_obj = NULL) {

  x_range <- range(tmp_data[[covariate_term]], na.rm = TRUE)

  # Prediction grid
  grid <- tibble::tibble(..x.. = seq(x_range[1], x_range[2], length.out = n_grid))
  names(grid)[1] <- covariate_term

  for (col in other_cols) {
    grid[[col]] <- ref_vals[[col]]
    if (col %in% finalCharacters && is.factor(tmp_data[[col]])) {
      grid[[col]] <- factor(grid[[col]], levels = levels(tmp_data[[col]]))
    }
  }

  pred <- stats::predict(fit$final_model, newdata = grid, type = "link", se.fit = TRUE)
  grid$PredictedProb <- stats::plogis(pred$fit)
  grid$LL            <- stats::plogis(pred$fit - z * pred$se.fit)
  grid$UL            <- stats::plogis(pred$fit + z * pred$se.fit)

  # Binned observed rates (Wilson CI)
  obs_rates <- NULL
  if (show_observed) {
    obs_rates <- ori_data %>%
      dplyr::mutate(bin = ggplot2::cut_number(.data[[covariate_term]], n = n_bins)) %>%
      dplyr::group_by(bin) %>%
      dplyr::summarise(
        n        = dplyr::n(),
        events   = sum(DV, na.rm = TRUE),
        rate     = events / n,
        mean_x   = mean(.data[[covariate_term]], na.rm = TRUE),
        ci_low   = (rate + z^2 / (2 * n) - z * sqrt((rate * (1 - rate) + z^2 / (4 * n)) / n)) / (1 + z^2 / n),
        ci_high  = (rate + z^2 / (2 * n) + z * sqrt((rate * (1 - rate) + z^2 / (4 * n)) / n)) / (1 + z^2 / n),
        .groups  = "drop"
      )
  }

  if (!plot) {
    return(list(predictions = grid, observed = obs_rates))
  }

  n_total  <- nrow(ori_data)
  n_events <- sum(ori_data$DV, na.rm = TRUE)
  subtitle_text <- paste0(
    "N=", n_total,
    "  Events=", n_events,
    "  Non-events=", n_total - n_events
  )

  events_data    <- ori_data %>% dplyr::filter(DV == 1)
  nonevents_data <- ori_data %>% dplyr::filter(DV == 0)

  p <- ggplot2::ggplot(grid, ggplot2::aes(x = .data[[covariate_term]])) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$LL, ymax = .data$UL),
      alpha = 0.2, fill = "#2166AC"
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$PredictedProb),
      color = "black", linewidth = 0.8
    )

  if (show_observed && !is.null(obs_rates)) {
    err_width <- diff(x_range) / (n_bins * 4)
    p <- p +
      ggplot2::geom_errorbar(
        data = obs_rates,
        ggplot2::aes(x = mean_x, ymin = ci_low, ymax = ci_high),
        color = "grey40", width = err_width, inherit.aes = FALSE
      ) +
      ggplot2::geom_point(
        data = obs_rates,
        ggplot2::aes(x = mean_x, y = rate),
        color = "grey40", size = 4, inherit.aes = FALSE
      )
  }

  if (show_rug) {
    p <- p +
      ggplot2::geom_rug(
        data = events_data,
        ggplot2::aes(x = .data[[covariate_term]]),
        sides = "t", color = "grey50", alpha = 0.4, inherit.aes = FALSE
      ) +
      ggplot2::geom_rug(
        data = nonevents_data,
        ggplot2::aes(x = .data[[covariate_term]]),
        sides = "b", color = "grey70", alpha = 0.4, inherit.aes = FALSE
      )
  }

  # Back-transform x-axis tick labels when data is on a transformed scale
  back_fn <- .get_back_transform_fn(var_spec_obj, covariate_term)
  x_scale <- if (identical(back_fn, identity)) {
    ggplot2::scale_x_continuous()
  } else {
    ggplot2::scale_x_continuous(labels = function(x) round(back_fn(x), 2))
  }

  p +
    x_scale +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::labs(
      x = x_label,
      y = paste0("P(", dv_label, ")"),
      subtitle = subtitle_text
    ) +
    ggplot2::theme_bw() +
    ggpubr::labs_pubr()
}


# Internal: categorical covariate path
#' @keywords internal
#' @noRd
.plot_logistic_categorical <- function(fit, covariate_term, ori_data, tmp_data,
                                       finalCharacters, other_cols, ref_vals,
                                       z, plot, dv_label,
                                       x_label = covariate_term,
                                       cat_labels = NULL) {

  col_data   <- tmp_data[[covariate_term]]
  cat_levels <- if (is.factor(col_data)) levels(col_data) else sort(unique(col_data))

  # Predicted probability per level
  pred_rows <- lapply(cat_levels, function(lev) {
    row <- tibble::tibble(!!covariate_term := lev)
    if (is.factor(col_data)) {
      row[[covariate_term]] <- factor(lev, levels = levels(col_data))
    }
    for (col in other_cols) {
      row[[col]] <- ref_vals[[col]]
      if (col %in% finalCharacters && is.factor(tmp_data[[col]])) {
        row[[col]] <- factor(row[[col]], levels = levels(tmp_data[[col]]))
      }
    }
    lp <- stats::predict(fit$final_model, newdata = row, type = "link", se.fit = TRUE)
    tibble::tibble(
      level    = as.character(lev),
      pred_prob = stats::plogis(lp$fit),
      pred_lo  = stats::plogis(lp$fit - z * lp$se.fit),
      pred_hi  = stats::plogis(lp$fit + z * lp$se.fit)
    )
  })
  pred_ds <- dplyr::bind_rows(pred_rows)

  # Observed event rate per level (Wilson CI)
  obs_ds <- ori_data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(covariate_term))) %>%
    dplyr::summarise(
      n       = dplyr::n(),
      events  = sum(DV, na.rm = TRUE),
      rate    = events / n,
      ci_low  = (rate + z^2 / (2 * n) - z * sqrt((rate * (1 - rate) + z^2 / (4 * n)) / n)) / (1 + z^2 / n),
      ci_high = (rate + z^2 / (2 * n) + z * sqrt((rate * (1 - rate) + z^2 / (4 * n)) / n)) / (1 + z^2 / n),
      .groups = "drop"
    ) %>%
    dplyr::mutate(level = as.character(.data[[covariate_term]])) %>%
    dplyr::select(-dplyr::all_of(covariate_term))

  if (!plot) {
    return(list(predictions = pred_ds, observed = obs_ds))
  }

  plot_data <- dplyr::bind_rows(
    pred_ds %>%
      dplyr::select(level, mid = pred_prob, lo = pred_lo, hi = pred_hi) %>%
      dplyr::mutate(source = "Model predicted"),
    obs_ds %>%
      dplyr::select(level, mid = rate, lo = ci_low, hi = ci_high) %>%
      dplyr::mutate(source = "Observed")
  )
  # Apply display labels to factor levels if provided via var_spec
  display_levels <- if (!is.null(cat_labels) && length(cat_labels) == length(cat_levels)) {
    cat_labels
  } else {
    as.character(cat_levels)
  }
  plot_data$level <- factor(plot_data$level, levels = as.character(cat_levels),
                            labels = display_levels)

  n_labels <- paste(obs_ds$level, paste0("N=", obs_ds$n), sep = ": ", collapse = ";  ")

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = level, y = .data$mid,
      color = .data$source, shape = .data$source
    )
  ) +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.4), size = 3
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
      position = ggplot2::position_dodge(width = 0.4),
      width = 0.15
    ) +
    ggplot2::scale_color_manual(
      values = c("Model predicted" = "#2166AC", "Observed" = "grey40")
    ) +
    ggplot2::scale_shape_manual(
      values = c("Model predicted" = 16L, "Observed" = 1L)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::labs(
      x = x_label,
      y = paste0("P(", dv_label, ")"),
      color = NULL, shape = NULL,
      subtitle = n_labels
    ) +
    ggplot2::theme_bw() +
    ggpubr::labs_pubr() +
    ggplot2::theme(legend.position = "bottom")
}
