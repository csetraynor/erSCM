#' Visual Predictive Check (VPC) for Cox Proportional Hazards Models
#'
#' Performs a visual predictive check (VPC) for a fitted Cox proportional hazards
#' model (`coxph`), by simulating survival data under the fitted model, computing
#' summary statistics (e.g., survival probabilities or cumulative hazards) over
#' time, and comparing simulated prediction intervals with the observed Kaplan-Meier
#' estimates. This helps assess model calibration and the adequacy of the
#' proportional hazards assumption in terms of predicted survival patterns.
#'
#' @param fit A list from \code{scm_coxph()} containing model results.
#' @param ci Numeric; confidence level for the prediction interval. Default 0.95.
#' @param n_bins Integer; number of bins (currently unused). Default 1.
#' @param plot Logical; if \code{TRUE} (default), return a \code{ggplot} object.
#'   If \code{FALSE}, return the computed summaries only.
#'
#' @details
#' The VPC procedure simulates survival data from the fitted Cox model by:
#' 1) computing linear predictors for each subject, 2) sampling event times
#' under the estimated baseline hazard and covariate effects, and optionally
#' 3) sampling censoring times to reflect the observed censoring pattern.
#' For each simulation, survival curves (or cumulative hazards) are summarized
#' at \code{times}. Across simulations, pointwise prediction intervals are
#' computed and contrasted with the observed Kaplan-Meier or Nelson-Aalen
#' estimates.
#'
#' Assumptions include adequate model fit, correct specification of the baseline
#' hazard estimation and covariate effects, and--when \code{censoring = TRUE}--
#' that the censoring mechanism is non-informative and adequately approximated.
#' Stratified VPCs enable visual checks within covariate-defined subgroups.
#'
#' @return If \code{plot = TRUE}, a \code{ggplot2} object showing:
#' - Observed curve (Kaplan-Meier survival or Nelson-Aalen cumulative hazard)
#' - Simulated median curve
#' - Shaded prediction interval ribbon (defined by \code{ci})
#' Optionally faceted by \code{stratify}. If \code{plot = FALSE}, returns a
#' list with components:
#' - \code{observed}: data frame of observed summary statistics by time and stratum
#' - \code{sim_q}: data frame of simulated quantiles by time and stratum
#' - \code{sim_median}: data frame of simulated medians by time and stratum
#' - \code{times}: evaluation times used
#' - \code{meta}: settings and fit details
#'
#' @section Simulation Approach:
#' Event times can be simulated via inversion using the estimated baseline
#' cumulative hazard \eqn{\hat{H}_0(t)} and linear predictor \eqn{\eta = x^\top \hat{\beta}},
#' with \eqn{H(t|x) = \hat{H}_0(t) \exp(\eta)}. Draw \eqn{U \sim \mathrm{Uniform}(0,1)},
#' set \eqn{T = H^{-1}(-\log(U))}, and optionally apply censoring. Baseline hazards
#' may be approximated using Breslow's estimator, spline-smoothed hazards, or
#' piecewise-constant segments depending on implementation.
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
#' @seealso \code{\link[survival]{coxph}}, \code{\link[survival]{Surv}}
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
  
  with_numeric <- length(finalNumeric) > 0
  with_bins <- with_numeric & n_bins > 1 
  
  model_with_categories <- length(finalCharacters) > 0 | with_bins
  
  ref_data <-   tmp_data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(finalCharacters))) %>%
    dplyr::summarise( dplyr::across(dplyr::all_of(finalNumeric), ~ mean(.x, na.rm = TRUE)), n = dplyr::n() ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(finalCharacters), ~ factor(as.character(.x), levels = sort(unique(as.character(.x)))))) %>%
    dplyr::arrange(dplyr::across(dplyr::all_of(finalCharacters))) %>% 
    dplyr::mutate(ID = dplyr::row_number()) 


  # Predicted baseline survival function (depends on chosen ref profile via lp=0 vs newdata)
  base_surv <- seq_len(nrow(ref_data)) %>% purrr::map(function(i){
    survival::survfit(fit$final_model, newdata = ref_data %>% dplyr::slice(i) )
  })

  # Observed KM (unadjusted)
  if(model_with_categories) {
    survformula <- paste0("survival::Surv(TIME, DV) ~", paste(finalCharacters, collapse = "+") )
  } else {
    survformula <- "survival::Surv(TIME, DV) ~ 1" 
  }

  km_fit <- survival::survfit(as.formula(survformula), data = fit$model_data)
  
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
        dplyr::arrange(strata, time, desc(estimate))
    
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
