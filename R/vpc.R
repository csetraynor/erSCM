#' Visual Predictive Check (VPC) for Cox Proportional Hazards Models
#'
#' Produces an analytical overlay for a fitted Cox proportional hazards model
#' from \code{\link{scm_coxph}} by contrasting observed Kaplan-Meier curves with
#' model-implied survival curves and their confidence band at reference profiles.
#' This is a confidence-band diagnostic (not a simulation-based predictive check).
#'
#' @param fit An object of class \code{eres_scm_coxph} returned by
#'   \code{\link{scm_coxph}}.
#' @param ci Numeric scalar in \code{(0, 1)}; confidence level for the model
#'   survival band. Default 0.95.
#' @param n_bins Integer scalar \code{>= 1}. Default 1.
#' @param plot Logical; if \code{TRUE} (default), return a \code{ggplot} object.
#'   If \code{FALSE}, return the computed summaries only.
#'
#' @details
#' The function builds one or more reference profiles from the fitted-model data
#' (grouped by categorical covariates and averaged over numeric covariates), then
#' obtains model-implied survival curves with \code{survival::survfit()} at those
#' profiles. The band shown is the pointwise \code{ci} confidence interval from
#' \code{survfit()} for each model-implied curve.
#'
#' Observed curves are unadjusted Kaplan-Meier estimates from
#' \code{survival::survfit(survival::Surv(TIME, DV) ~ ...)} on
#' \code{fit$model_data}. The overlay is therefore an analytical
#' confidence-band diagnostic over time, not a parametric-bootstrap predictive
#' band. A simulation-based predictive-check path is deferred to
#' \code{reports/vpc_coxph_v2_plan.md}.
#'
#' @return If \code{plot = TRUE}, a \code{ggplot2} object showing:
#' - Observed Kaplan-Meier curve(s)
#' - Model-implied survival curve(s)
#' - Shaded \code{ci} confidence band around the model-implied curve(s)
#' Optionally faceted by stratum. If \code{plot = FALSE}, returns the bare long
#' data frame used for plotting.
#'
#' @section Simulation Approach:
#' The current implementation is analytical and does not simulate event times.
#' A simulation-based predictive VPC for Cox models is deferred to
#' \code{reports/vpc_coxph_v2_plan.md}.
#'
#' @examples
#' data(demo_eres)
#' demo_eres$SEX <- as.factor(demo_eres$SEX)
#' fit <- scm_coxph(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", covariates = c("AGE", "WT"),
#'   dv_col = "sim_status", time_col = "sim_time", significance_level = c(0.2, 0.2)
#' )
#' vpc_coxph(fit)
#'
#' @seealso \code{\link{vpc_coxph_sim}}, \code{\link{vpc_logistic}}, \code{\link{plot_cox_terms}},
#'   \code{\link{scm_coxph}}, \code{\link[survival]{coxph}},
#'   \code{\link[survival]{Surv}}
#'
#' @importFrom survival coxph Surv basehaz survfit
#' @importFrom stats runif model.matrix
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon facet_wrap labs theme_minimal
#'
#' @references
#' - Kalbfleisch JD, Prentice RL. The Statistical Analysis of Failure Time Data.
#' - Cox DR. Regression Models and Life-Tables (1972).
#'
#' @export
vpc_coxph <- function(fit, ci = 0.95, n_bins = 1, plot = TRUE) {

  if (!inherits(fit, "eres_scm_coxph")) {
    cli::cli_abort(c(
      "{.arg fit} must be an {.cls eres_scm_coxph} object from {.fun scm_coxph}.",
      "x" = "Got an object of class {.cls {class(fit)}}."
    ))
  }
  if (is.null(fit$final_model)) {
    cli::cli_abort(c(
      "{.arg fit} has no {.field final_model}.",
      "i" = "No model was selected, so there is nothing to check."
    ))
  }
  if (!is.numeric(ci) || length(ci) != 1L || is.na(ci) || ci <= 0 || ci >= 1) {
    cli::cli_abort("{.arg ci} must be a numeric scalar in (0, 1).")
  }
  if (!is.numeric(n_bins) || length(n_bins) != 1L || is.na(n_bins) ||
      n_bins != as.integer(n_bins) || n_bins < 1) {
    cli::cli_abort("{.arg n_bins} must be a single integer >= 1.")
  }
  n_bins <- as.integer(n_bins)

  var_spec_obj <- fit$var_spec

  # Build time-axis label from var_spec if available
  time_col_name <- fit$metadata$time_col
  if (!is.null(var_spec_obj) && !is.null(time_col_name)) {
    time_display <- get_var_label(var_spec_obj, time_col_name)
    time_units   <- get_var_units(var_spec_obj, time_col_name)
    x_label <- if (time_display != time_col_name && !is.na(time_units) && !is.null(time_units)) {
      paste0(time_display, " (", time_units, ")")
    } else if (time_display != time_col_name) {
      time_display
    } else {
      "Time"
    }
  } else {
    x_label <- "Time"
  }

  # Select only main-effect columns (interaction terms are formula constructs, not data columns)
  main_effect_cols <- get_main_effect_columns(fit$metadata$finalRelation)

  tmp_data <- fit$model_data %>%
    dplyr::select(dplyr::all_of(main_effect_cols))

  finalCharacters <- dplyr::bind_cols(
    dplyr::select(tmp_data, dplyr::where(is.factor)) ,
    dplyr::select(tmp_data, dplyr::where(is.character)) ) %>%
    dplyr::bind_cols(dplyr::select(tmp_data, dplyr::where(is.logical)))  %>%
    colnames()
  
  finalNumeric <- dplyr::select(tmp_data, dplyr::where(is.numeric)) %>% colnames()

  if (n_bins > 1L && length(finalNumeric) > 0L) {
    cli::cli_abort(c(
      "{.arg n_bins} > 1 is not implemented for numeric covariates in {.fun vpc_coxph}.",
      "x" = "Numeric terms in the fitted main effects: {.val {finalNumeric}}.",
      "i" = "This numeric-binning path is deferred; see {.file reports/vpc_coxph_v2_plan.md}."
    ))
  }
  
  with_numeric <- length(finalNumeric) > 0
  with_bins <- with_numeric & n_bins > 1 
  
  model_with_categories <- length(finalCharacters) > 0 | with_bins
  
  ref_data <-   tmp_data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(finalCharacters))) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(finalNumeric), ~ mean(.x, na.rm = TRUE)),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(finalCharacters), ~ factor(as.character(.x), levels = sort(unique(as.character(.x)))))) %>%
    dplyr::arrange(dplyr::across(dplyr::all_of(finalCharacters))) %>% 
    dplyr::mutate(ID = dplyr::row_number()) 


  # Predicted baseline survival function (depends on chosen ref profile via lp=0 vs newdata)
  base_surv <- seq_len(nrow(ref_data)) %>% purrr::map(function(i){
    survival::survfit(fit$final_model, newdata = ref_data %>% dplyr::slice(i) )
  })

  # Observed KM (unadjusted)
  if(model_with_categories) {
    # TIME/DV are hard-coded here because scm_coxph() normalizes model_data by
    # renaming user-specified time_col/dv_col to TIME/DV internally.
    survformula <- paste0("survival::Surv(TIME, DV) ~", paste(finalCharacters, collapse = "+") )
  } else {
    survformula <- "survival::Surv(TIME, DV) ~ 1" 
  }

  km_fit <- survival::survfit(stats::as.formula(survformula), data = fit$model_data)
  
  obs <- broom::tidy(km_fit) %>%
    dplyr::mutate(origin = "Observed KM")
  
  ci_label <- paste0("Cox model (",round(ci*100),"% CI)")

  if(model_with_categories) {
    
    obs <- dplyr::bind_rows(
      obs,
      obs %>% 
        dplyr::group_by(strata, origin) %>%
        dplyr::filter(dplyr::row_number() == 1) %>% 
        dplyr::mutate(n.event = 0, time = 0,
                      estimate = 1, conf.high = 1, conf.low = 1) ) %>%
        dplyr::arrange(strata, time, dplyr::desc(estimate))
    
    stratas <- obs %>%
      dplyr::distinct(strata) %>% dplyr::pull()
    
    cox_ph <- base_surv %>% purrr::map(function(x) broom::tidy(x, conf.level = ci) %>% dplyr::mutate(origin = ci_label) ) %>%
      rlang::set_names(stratas ) %>% dplyr::bind_rows(.id = "strata")
    
    vpc_ds <- dplyr::bind_rows(obs, cox_ph) %>%
      dplyr::mutate(strata = gsub("=", "\n", strata))

    # Replace raw column names and factor levels in strata with display labels
    if (!is.null(var_spec_obj)) {
      for (cv in finalCharacters) {
        cv_label <- get_var_label(var_spec_obj, cv)
        vpc_ds <- vpc_ds %>%
          dplyr::mutate(strata = base::gsub(
            paste0("^", cv, "\n"), paste0(cv_label, "\n"), strata))
        cats <- get_var_categories(var_spec_obj, cv)
        if (!is.null(cats)) {
          for (j in seq_along(cats)) {
            raw_lev <- names(cats)[j]
            vpc_ds <- vpc_ds %>%
              dplyr::mutate(strata = base::gsub(
                paste0("\n", raw_lev, "$"), paste0("\n", cats[[j]]), strata))
          }
        }
      }
    }

    # Preserve identity-derived strata order for facets (avoid lexical sorting).
    strata_levels <- unique(vpc_ds$strata)
    vpc_ds$strata <- factor(vpc_ds$strata, levels = strata_levels, ordered = TRUE)
    
  } else {
    cox_ph <- base_surv %>% purrr::map(function(x) broom::tidy(x, conf.level = ci) %>% dplyr::mutate(origin = ci_label) ) %>%
      dplyr::bind_rows() 
    
    vpc_ds <- dplyr::bind_rows(obs, cox_ph) 
  }


  if(!plot) {

    return(vpc_ds)

  } else {
    
    if(model_with_categories) {
      
      vpc_plot <- ggplot2::ggplot(data = vpc_ds, ggplot2::aes( x = time, y = estimate)) +
        ggplot2::geom_step(data = vpc_ds %>% dplyr::filter(origin == "Observed KM"),
                           ggplot2::aes(linetype = origin)) +
        ggplot2::geom_ribbon(data = vpc_ds %>% dplyr::filter(origin == ci_label),
                             ggplot2::aes(ymin = conf.low, ymax = conf.high, fill = origin),
                             # stat = "stepribbon", 
                             alpha = 0.3, col = "transparent") +
        ggplot2::geom_line(data = vpc_ds %>% dplyr::filter(origin == ci_label),
                           ggplot2::aes(col = origin)) +
        ggplot2::scale_fill_brewer(type = "qual", palette = "Set1") +
        ggplot2::theme_bw() +
        ggplot2::facet_wrap( ~ strata) +
        #labs_pubr()  +
        ggplot2::labs(col = "", fill = "", linetype = "", y = "Survival Probability", x = x_label) +
        ggplot2::scale_y_continuous(limits = c(0,1))
      
    } else {
      
      vpc_plot <- ggplot2::ggplot(data = vpc_ds, ggplot2::aes( x = time, y = estimate)) +
        ggplot2::geom_step(data = vpc_ds %>% dplyr::filter(origin == "Observed KM"),
                           ggplot2::aes(linetype = origin)) +
        ggplot2::geom_ribbon(data = vpc_ds %>% dplyr::filter(origin == ci_label),
                             ggplot2::aes(ymin = conf.low, ymax = conf.high, fill = origin),
                             # stat = "stepribbon", 
                             alpha = 0.3, col = "transparent") +
        ggplot2::geom_line(data = vpc_ds %>% dplyr::filter(origin == ci_label),
                           ggplot2::aes(col = origin)) +
        ggplot2::scale_fill_brewer(type = "qual", palette = "Set1") +
        ggplot2::theme_bw() +
        #labs_pubr()  +
        ggplot2::labs(col = "", fill = "", linetype = "", y = "Survival Probability", x = x_label) +
        ggplot2::scale_y_continuous(limits = c(0,1))
      
    }

    return(vpc_plot)

  }

}
