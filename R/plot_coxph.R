#' Plot Observed vs. Predicted Survival Stratified by a Covariate Term
#'
#' Visualizes model fit for a single covariate term from a fitted Cox proportional
#' hazards SCM object, by stratifying the model data on \code{covariate_term}
#' (binning numeric covariates into \code{n_bins} groups) and overlaying, per
#' stratum, the observed Kaplan-Meier survival curve against the model-predicted
#' survival curve with a confidence ribbon.
#'
#' @param fit An object of class \code{eres_scm_coxph} (as returned by
#'   \code{\link{scm_coxph}} or \code{\link{run_scm}}).
#' @param covariate_term Character; name of a main-effect covariate term in
#'   \code{fit} to stratify by. Interaction terms (containing \code{":"}) are
#'   not supported -- use one of the parent main effects instead.
#' @param ci Numeric; confidence level for the predicted survival ribbon. Default 0.95.
#' @param n_bins Integer; number of bins used to discretize a numeric
#'   \code{covariate_term} into strata. Ignored for factor/character covariates.
#'   Default 4.
#' @param plot Logical; if \code{TRUE} (default), return a \code{ggplot} object.
#'   If \code{FALSE}, return the underlying data instead.
#'
#' @return If \code{plot = TRUE}, a \code{ggplot2} object with one panel per
#'   stratum of \code{covariate_term}, showing the observed Kaplan-Meier step
#'   curve and the model-predicted survival curve with a shaded confidence
#'   ribbon. If \code{plot = FALSE}, a data.frame with columns \code{time},
#'   \code{estimate}, \code{conf.low}, \code{conf.high}, \code{origin}
#'   (\code{"Observed KM"} or the predicted-curve label), and \code{strata}.
#'
#' @examples
#' data(demo_eres)
#' demo_eres$SEX <- as.factor(demo_eres$SEX)
#' fit <- scm_coxph(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", covariates = c("AGE", "WT"),
#'   dv_col = "sim_status", time_col = "sim_time", significance_level = c(0.2, 0.2)
#' )
#' plot_cox_terms(fit, "WT")
#'
#' @seealso \code{\link[survival]{coxph}}, \code{\link[survival]{Surv}},
#'   \code{\link{plot_logistic_terms}}, \code{\link{plot_model_terms}}
#'
#' @importFrom survival coxph Surv basehaz survfit
#' @importFrom stats runif model.matrix
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon facet_wrap labs theme_minimal
#' @importFrom tidyselect where
#'
#' @references
#' - Kalbfleisch JD, Prentice RL. The Statistical Analysis of Failure Time Data.
#' - Cox DR. Regression Models and Life-Tables (1972).
#'
#' @export
plot_cox_terms <- function(fit, covariate_term, ci = 0.95, n_bins = 4, plot = TRUE) {
  
  # Interaction terms cannot be used as covariate_term (no single data column to stratify by)
  if (grepl(":", covariate_term, fixed = TRUE)) {
    main_effects <- get_main_effect_columns(fit$metadata$finalRelation)
    cli::cli_abort(c(
      "{.fun plot_cox_terms} does not support interaction terms as {.arg covariate_term}.",
      "i" = "Use a main-effect variable name, e.g. one of: {.val {main_effects}}."
    ))
  }

  # Select only main-effect columns (interaction terms are formula constructs, not data columns)
  main_effect_cols <- get_main_effect_columns(fit$metadata$finalRelation)

  ori_data <- fit$model_data %>%
    dplyr::select(dplyr::all_of(main_effect_cols), DV, TIME)
  
  tmp_data <- ori_data %>% dplyr::select(-TIME, -DV)
  
  finalCharacters <- dplyr::bind_cols(
    dplyr::select(tmp_data, dplyr::where(is.factor)) ,
    dplyr::select(tmp_data, dplyr::where(is.character)) ) %>%
    dplyr::bind_cols(dplyr::select(tmp_data, dplyr::where(is.logical)))  %>%
    colnames()
  
  finalNumeric <- dplyr::select(tmp_data, dplyr::where(is.numeric)) %>% colnames()
  
  with_numeric <- length(finalNumeric) > 0
  with_bins <- with_numeric & n_bins > 1 
  
  model_with_factors <- length(finalCharacters) > 0 
  
  model_with_categories <- model_with_factors | with_bins
  
  ref_data <-   tmp_data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(finalCharacters))) %>%
    dplyr::summarise( dplyr::across(dplyr::all_of(finalNumeric), ~ mean(.x, na.rm = TRUE)), n = dplyr::n(), .groups = "drop" ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(finalCharacters), ~ factor(as.character(.x), levels = sort(unique(as.character(.x)))))) %>%
    dplyr::arrange(dplyr::across(dplyr::all_of(finalCharacters))) %>%
    dplyr::mutate(ID = dplyr::row_number())

  if(model_with_factors) {

    ref_data <-   ref_data %>%
      tidyr::pivot_longer(dplyr::all_of(finalCharacters)) %>%
      dplyr::group_by(name) %>%
      dplyr::filter(n > stats::quantile(n, .8)) %>%
      dplyr::filter(dplyr::row_number() == 1) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_wider(names_from = "name", values_from = "value")

    ori_data <- ori_data %>%
      dplyr::inner_join( ref_data %>% dplyr::select(dplyr::all_of(finalCharacters)), by = finalCharacters)

    ref_data <-   ori_data %>%
      dplyr::select(-TIME, -DV) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(finalCharacters))) %>%
      dplyr::summarise( dplyr::across(dplyr::all_of(finalNumeric), ~ mean(.x, na.rm = TRUE)), n = dplyr::n(), .groups = "drop" ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(ID = dplyr::row_number())

  }
  
  covariate_data <- tmp_data %>% dplyr::select( dplyr::all_of(covariate_term) )
  
  is_covariate_numeric <- is.numeric( dplyr::pull(covariate_data, covariate_term) )
  
  # Resolve display label and units from var_spec
  var_spec_obj   <- fit$var_spec
  display_label  <- get_var_label(var_spec_obj, covariate_term)
  display_units  <- get_var_units(var_spec_obj, covariate_term)
  legend_title   <- if (!is.na(display_units) && !is.null(display_units)) {
    paste0(display_label, " (", display_units, ")")
  } else {
    display_label
  }

  if(is_covariate_numeric) {

    # Compute break points for back-transformation to original scale
    cov_vals <- tmp_data[[covariate_term]]
    bin_breaks <- stats::quantile(cov_vals, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)

    # Build back-transformed range labels when var_spec is available
    if (!is.null(var_spec_obj)) {
      bin_display_labels <- make_bin_labels(bin_breaks, var_spec_obj, covariate_term)
    } else {
      bin_display_labels <- paste0("Q", seq_len(n_bins))
    }

    ori_data <- ori_data %>%
      dplyr::mutate( quartile_cat = ggplot2::cut_number(.data[[covariate_term]], n = n_bins,
                                                         labels = bin_display_labels),
              quartile_mean = stats::ave(.data[[covariate_term]], quartile_cat, FUN = function(z) base::mean(z, na.rm = TRUE)) )

    ref_data_expanded <- ori_data %>%
      dplyr::select(-TIME, -DV) %>%
      dplyr::group_by(.data$quartile_cat) %>%
      dplyr::summarise(dplyr::across(tidyselect::where(is.numeric), stats::median), .groups = "drop")

    # Broadcast the single modal reference profile for any categorical covariate
    # in the model across every numeric bin. Assert exactly one such profile row
    # rather than relying on `bind_cols()`'s silent size-1 recycling, which
    # errors cryptically (or mis-recycles) if the upstream modal-profile
    # reduction did not collapse to one row.
    ref_factor_profile <- ref_data %>% dplyr::select(dplyr::all_of(finalCharacters))
    if (ncol(ref_factor_profile) > 0) {
      if (nrow(ref_factor_profile) != 1) {
        cli::cli_abort(c(
          "Expected exactly one modal reference profile for the categorical covariate(s).",
          "x" = "Got {nrow(ref_factor_profile)} row{?s} for {.val {finalCharacters}}.",
          "i" = "The upstream modal-profile reduction did not collapse to a single row."
        ))
      }
      ref_data_expanded <- dplyr::bind_cols(ref_data_expanded, ref_factor_profile)
    }

    survformula <- "survival::Surv(TIME, DV) ~ quartile_cat"

  } else {
    
    ref_data_categorical <- ref_data %>% 
      dplyr::select(dplyr::all_of(finalCharacters), -tidyselect::all_of(covariate_term) )
    
    if(ncol(ref_data_categorical) > 0) {
      # Join key = the categorical covariates other than the one being stratified
      # on (the right table is `ref_data` restricted to exactly those columns).
      cat_join_keys <- base::setdiff(finalCharacters, covariate_term)
      ref_data_expanded <-  fit$model_data %>%
        dplyr::select(dplyr::all_of(main_effect_cols)) %>%
        dplyr::inner_join(ref_data %>% dplyr::select(dplyr::all_of(finalCharacters), -dplyr::all_of(covariate_term) ), by = cat_join_keys ) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(finalCharacters))) %>%
        dplyr::summarise(dplyr::across(tidyselect::where(is.numeric), stats::median), .groups = "drop")
      
      ori_data <- fit$model_data %>%
        dplyr::select(TIME, DV, dplyr::all_of(main_effect_cols)) %>%
        dplyr::inner_join(ref_data %>% dplyr::select(dplyr::all_of(finalCharacters), -dplyr::all_of(covariate_term) ), by = cat_join_keys )
      
      
    } else {
      ref_data_expanded <-  fit$model_data %>%
        dplyr::select(dplyr::all_of(main_effect_cols))  %>%
        dplyr::mutate(dplyr::across(dplyr::all_of(finalCharacters), ~ factor(as.character(.x), levels = sort(unique(as.character(.x)))))) %>%
        dplyr::arrange(dplyr::across(dplyr::all_of(finalCharacters))) %>% 
        dplyr::group_by(dplyr::across(dplyr::all_of(covariate_term))) %>%
        dplyr::summarise(dplyr::across(tidyselect::where(is.numeric), stats::median), .groups = "drop")
      
      ori_data <-  fit$model_data 
    }
    
    survformula <- paste0("survival::Surv(TIME, DV) ~ ", covariate_term)
    
  }
  
  km_fit <- survival::survfit(as.formula(survformula), data = ori_data )
  
  obs <- broom::tidy(km_fit) %>%
    dplyr::mutate(origin = "Observed KM")
  
  if("strata" %in% colnames(obs)) {
    obs <- dplyr::bind_rows(
      obs,
      obs %>% 
        dplyr::group_by(strata, origin) %>%
        dplyr::filter(dplyr::row_number() == 1) %>% 
        dplyr::mutate(n.event = 0, time = 0,
                      estimate = 1, conf.high = 1, conf.low = 1) ) %>%
      dplyr::arrange(strata, time, dplyr::desc(estimate))
  }

  ci_label <- paste0("Cox model (",round(ci*100),"% CI)")
  
  # Build the survfit-style strata label for EACH reference profile directly
  # from its own key column so that predicted curves are matched to the
  # CORRECT stratum by identity rather than by lexical position. The previous
  # implementation named the predicted curves positionally from
  # `distinct(strata)` (lexically ordered by dplyr::arrange), which misaligns
  # whenever the bin/level labels do not sort lexically in numeric order
  # (e.g. var_spec back-transformed labels using "<=", ">", en-dash, or
  # n_bins >= 10 where "Q10" < "Q2"). Under dplyr's C-locale arrange() default
  # (>= 1.0.0) this produced predicted curves attached to the wrong facet.
  strata_var   <- if (is_covariate_numeric) "quartile_cat" else covariate_term
  ref_keys     <- as.character(
    dplyr::pull(dplyr::ungroup(ref_data_expanded), !!rlang::sym(strata_var))
  )
  strata_names <- paste0(strata_var, "=", ref_keys)

  # Predicted survival function per reference profile (row of ref_data_expanded)
  base_surv <- seq_len(nrow(ref_data_expanded)) %>% purrr::map(function(i){
    survival::survfit(fit$final_model, newdata = ref_data_expanded %>% dplyr::ungroup() %>% dplyr::slice(i) )
  })
  
  cox_ph <- base_surv %>% purrr::map(function(x) broom::tidy(x, conf.level = ci) %>% dplyr::mutate(origin = ci_label) ) %>%
    rlang::set_names(strata_names) %>% dplyr::bind_rows(.id = "strata")

  # Display-only label transformation. Encapsulated in a single function so the
  # exact same steps are applied to (a) the plotted strata values and (b) the
  # derived facet-level order below, guaranteeing they always agree. This does
  # NOT alter label *content* and does NOT touch the identity-based strata join
  # built above (Phase A).
  transform_strata_labels <- function(s) {
    s <- base::gsub("=", "\n", s)
    if (is_covariate_numeric) {
      s <- base::gsub("quartile_cat", legend_title, s)
    } else {
      # For categorical, replace raw column name with display label in strata strings
      s <- base::gsub(paste0("^", covariate_term, "="), paste0(legend_title, "="), s)
      # Relabel factor levels in strata using var_spec categories if available
      cats <- get_var_categories(var_spec_obj, covariate_term)
      if (!is.null(cats)) {
        for (raw_lev in names(cats)) {
          s <- base::gsub(paste0("=\n", raw_lev, "$"), paste0("=\n", cats[[raw_lev]]), s)
        }
      }
    }
    s
  }

  vpc_ds <- dplyr::bind_rows(obs, cox_ph) %>%
    dplyr::mutate(strata = transform_strata_labels(strata))

  # Finding #3 (display-only): render facets in the bin/level order taken from
  # `ref_data_expanded` (the identity source of truth used by the Phase A fix)
  # instead of the lexical order of the transformed label strings. `strata_names`
  # is built row-wise from `ref_data_expanded`, so it is already in
  # numeric-ascending (numeric) / factor-level (categorical) order; applying the
  # SAME label transform makes the level vector match the displayed labels
  # exactly.
  ordered_strata_levels <- unique(transform_strata_labels(strata_names))

  # Never silently drop panels: every plotted stratum must map to a known level.
  unmatched_strata <- base::setdiff(
    unique(as.character(vpc_ds$strata)), ordered_strata_levels
  )
  if (length(unmatched_strata) > 0) {
    cli::cli_abort(c(
      "Internal error: plotted strata not found in the derived facet ordering.",
      "x" = "Unmatched strata label{?s}: {.val {unmatched_strata}}.",
      "i" = "This means {.var vpc_ds$strata} diverged from {.var ref_data_expanded} keys."
    ))
  }

  vpc_ds <- vpc_ds %>%
    dplyr::mutate(strata = factor(strata, levels = ordered_strata_levels, ordered = TRUE))


  if(!plot) {
    
    return(vpc_ds)
    
  } else {
    
    if("strata" %in% colnames(vpc_ds)) {
      
      vpc_plot <- ggplot2::ggplot(data = vpc_ds, ggplot2::aes( x = time, y = estimate)) +
        ggplot2::geom_step(data = vpc_ds %>% dplyr::filter(origin == "Observed KM"),
                           ggplot2::aes(linetype = origin, col = origin)) +
        ggplot2::geom_ribbon(data = vpc_ds %>% dplyr::filter(origin == ci_label),
                             ggplot2::aes(ymin = conf.low, ymax = conf.high, fill = origin),
                             alpha = 0.3, col = "transparent", show.legend = FALSE) +
        ggplot2::geom_ribbon(data = vpc_ds %>% dplyr::filter(origin == "Observed KM"),
                             ggplot2::aes(ymin = conf.low, ymax = conf.high, fill = origin),
                             alpha = 0.3, col = "transparent", show.legend = FALSE) +
        ggplot2::geom_line(data = vpc_ds %>% dplyr::filter(origin == ci_label),
                           ggplot2::aes(linetype = origin, col = origin)) +
        ggplot2::scale_fill_brewer(type = "qual", palette = "Set1") +
        ggplot2::scale_colour_brewer(type = "qual", palette = "Set1") +
        ggplot2::scale_linetype_manual(values = c("solid", "solid")) +
        ggplot2::theme_bw() +
        ggplot2::facet_wrap( ~ strata) +
        ggplot2::labs(col = "", fill = "", linetype = "", y = "Survival Probability",
                      x = "Time", title = legend_title) +
        ggplot2::scale_y_continuous(limits = c(0,1))

    } else {

      vpc_plot <- ggplot2::ggplot(data = vpc_ds, ggplot2::aes( x = time, y = estimate)) +
        ggplot2::geom_step(data = vpc_ds %>% dplyr::filter(origin == "Observed KM"),
                           ggplot2::aes(linetype = origin)) +
        ggplot2::geom_ribbon(data = vpc_ds %>% dplyr::filter(origin == ci_label),
                             ggplot2::aes(ymin = conf.low, ymax = conf.high, fill = origin),
                             alpha = 0.3, col = "transparent") +
        ggplot2::geom_line(data = vpc_ds %>% dplyr::filter(origin == ci_label),
                           ggplot2::aes(col = origin)) +
        ggplot2::scale_fill_brewer(type = "qual", palette = "Set1") +
        ggplot2::theme_bw() +
        ggplot2::labs(col = "", fill = "", linetype = "", y = "Survival Probability",
                      x = "Time", title = legend_title) +
        ggplot2::scale_y_continuous(limits = c(0,1))
      
    }
    
    return(vpc_plot)
    
  }
  
}
