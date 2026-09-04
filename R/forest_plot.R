#' Create a forest plot from a fitted model
#'
#' @description
#' Generates a forest plot of effect estimates (e.g., odds ratios, hazard ratios, or
#' mean differences) from a fitted model object. Confidence intervals are computed
#' at the specified level, and optional label mappings and level additions can be applied.
#' The function can either return the plot object or just the processed data, depending
#' on the `plot` argument.
#'
#' For interaction terms, conditional effects are computed at meaningful quantile
#' values using the model's variance-covariance matrix for proper confidence intervals.
#' Continuous x continuous interactions show the effect of one variable at quantiles of
#' the other. Continuous x categorical interactions show the continuous effect stratified
#' by category levels.
#'
#' @param fit An object containing model results to be plotted. 
#' @param ci Numeric scalar specifying the confidence level for intervals (between 0 and 1).
#'   Defaults to `0.95` for 95% confidence intervals.
#' @param plot Logical indicating whether to return a plot (`TRUE`) or only the processed
#'   data structure (`FALSE`). Defaults to `TRUE`.
#' @param namings Optional dataframe with named character vector used to relabel variables or levels
#'   in the plot (e.g., `data.frame(covname = c("CMAXss", "WT"),  units = c("ug*day/mL", "Kg") )`). Optionally, include  `median_point` to inverse standardisation.
#'   Defaults to `NULL` (no relabeling).
#' @param n_thres Integer threshold for filtering variables by sample size or count (e.g.,
#'   exclude strata with counts below `n_thres`). Set to `0` to disable filtering. Defaults to `0`.
#' @param log_x Logical; if TRUE, use log10 x-axis scale. Defaults to FALSE.
#' @param add_lev Optional vector specifying additional levels or categories to append to
#'   the plotting data (e.g., reference levels). Defaults to `NULL`.
#' @param quantiles_interaction Numeric vector of quantile probabilities at which to display
#'   conditional effects for the "effect variable" in interaction terms. Defaults to
#'   `c(0.05, 0.25, 0.75, 0.95)`. Values should be between 0 and 1.
#' @param modifier_quantiles Numeric vector of quantile probabilities at which the
#'   "modifier variable" is evaluated for continuous x continuous interactions.
#'   Defaults to `c(0.25, 0.5, 0.75)`. Ignored for continuous x categorical interactions
#'   where category levels are used instead.
#' @param adjust_for Optional character vector of covariate names that were forced into
#'   the model (i.e., passed as \code{adjust_for} in \code{\link{scm_glm}} or
#'   \code{\link{scm_coxph}}). When supplied, forced covariates receive a distinct visual
#'   style (colour and point shape) in the plot. When \code{NULL} (the default), all main
#'   effects are styled identically, preserving the existing appearance.
#' @param show_terms Character string controlling which covariate panels to include.
#'   One of \code{"all"} (default), \code{"selected"} (SCM-selected terms only), or
#'   \code{"adjusted"} (forced-adjustment terms only). Has no effect when
#'   \code{adjust_for = NULL}.
#' @param adjust_colour Colour used for forced-adjustment covariate markers.
#'   Defaults to \code{"grey50"}.
#' @param selected_colour Colour used for SCM-selected covariate markers.
#'   Defaults to \code{"steelblue"}.
#' @param adjust_shape Integer point shape for forced-adjustment markers.
#'   Defaults to \code{17L} (filled triangle).
#' @param selected_shape Integer point shape for SCM-selected markers.
#'   Defaults to \code{16L} (filled circle).
#'
#' @details
#' The function extracts effect estimates and confidence intervals from search `fit` at the
#' specified confidence level. If `namings` is provided, variable names or levels in the
#' plot are mapped accordingly. If `n_thres` is greater than zero, strata or groups with
#' counts below the threshold are removed prior to plotting. `add_lev` can be used to
#' include additional labels or categories (such as reference groups) for clarity.
#'
#' For interaction terms, the conditional effect computation uses the variance-covariance
#' matrix of the fitted model to derive proper confidence intervals. The conditional
#' log-HR (or log-OR) for a change in X1 from median to quantile q, at X2 = x2, is:
#' \deqn{\delta = (\beta_1 + \beta_{12} \cdot x_2) \cdot (q - \text{median}_{X1})}{
#' delta = (beta_1 + beta_12 * x2) * (q - median_X1)}
#' with variance:
#' \deqn{Var(\delta) = (q - m)^2 [Var(\beta_1) + x_2^2 Var(\beta_{12}) + 2 x_2 Cov(\beta_1, \beta_{12})]}{
#' Var(delta) = (q-m)^2 * [Var(beta_1) + x2^2*Var(beta_12) + 2*x2*Cov(beta_1, beta_12)]}
#'
#' The y-axis typically lists variables or factor levels, while the x-axis shows the
#' effect size on an appropriate scale (often log scale for ratios). A vertical reference
#' line (e.g., at 1 for ratios or 0 for differences) may be included depending on the
#' implementation.
#'
#' @return
#' If `plot = TRUE`, returns a plot object (commonly a `ggplot` object) representing the
#' forest plot. If `plot = FALSE`, returns a data frame (or tibble) with the processed
#' estimates and confidence intervals used to build the plot.
#'
#' @examples
#' data(demo_eres)
#' demo_eres$SEX <- as.factor(demo_eres$SEX)
#' lr_fit <- scm_glm(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", covariates = c("AGE", "WT"),
#'   dv_col = "sim_status", significance_level = c(0.2, 0.2)
#' )
#'
#' # Return processed data only
#' dat <- make_forest_plot(lr_fit, plot = FALSE)
#'
#' # Custom labels and filtering
#' labels <- c(CMAXss = "Cmax,ss (ng/mL)", WT = "Body Weight (kg)")
#' p2 <- make_forest_plot(lr_fit, ci = 0.90, namings = labels)
#'
#' @seealso
#' - `broom::tidy()` for tidying model outputs
#' - `ggplot2` for plotting
#'
#' @export
make_forest_plot <- function(fit,
                             ci = 0.95,
                             plot = TRUE,
                             log_x = FALSE,
                             namings = NULL,
                             n_thres = 0,
                             add_lev = NULL,
                             quantiles_interaction = c(0.05, 0.25, 0.75, 0.95),
                             modifier_quantiles = c(0.25, 0.5, 0.75),
                             adjust_for = NULL,
                             show_terms = c("all", "selected", "adjusted"),
                             adjust_colour = "grey50",
                             selected_colour = "steelblue",
                             adjust_shape = 17L,
                             selected_shape = 16L) {

  show_terms <- match.arg(show_terms)

  # Separate main-effect columns from interaction terms
  all_relation <- fit$metadata$finalRelation
  main_effect_cols <- get_main_effect_columns(all_relation)
  interaction_terms <- get_interaction_terms(all_relation)

  tmp_data <- fit$model_data %>%
    dplyr::select(dplyr::all_of(main_effect_cols)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(id = 1:dplyr::n())

  finalCharacters <- dplyr::bind_cols(
    dplyr::select(tmp_data, dplyr::where(is.factor)) ,
    dplyr::select(tmp_data, dplyr::where(is.character)) ) %>%
    dplyr::bind_cols(dplyr::select(tmp_data, dplyr::where(is.logical)))  %>%
    colnames()
  finalNumeric <- dplyr::select(tmp_data, dplyr::where(is.numeric)) %>% colnames()


  covariates <- main_effect_cols
  exposure_metrics <- fit$metadata$candidate_exposure_vars

  exposurecovars <- covariates[covariates %in% exposure_metrics]
  contcovars <- covariates[covariates %in% finalNumeric & covariates %!in% exposure_metrics]
  catcovars <- covariates[covariates %in% finalCharacters]
  
  # Retrieve centering values for label back-transformation
  centering_values <- fit$metadata$centering_values

  if(length(exposurecovars) > 0 ) {
    tmp_exposuredata <- list()

    for(i in seq_along(exposurecovars)) {
      vname <- exposurecovars[i]
      df_x <- tmp_data %>% dplyr::rename(x := !!vname)
      med_x <- stats::median(df_x$x)
      # Back-transform center offset for display labels
      center_offset <- if (!is.null(centering_values) && vname %in% names(centering_values)) centering_values[[vname]] else 0

      tmp_exposuredata[[i]] <- df_x %>%
        dplyr::mutate(Q =  cut(x,
                        breaks=quantile(unique(df_x[,c('x','id')])$x,c(0,.25,.5,.75,1), na.rm =TRUE),
                        include.lowest=TRUE,
                        labels=paste0("Q",c(1,2,3,4)))) %>%
        dplyr::summarise( med = stats::median(x,na.rm=TRUE),
                   qlow =quantile(x,na.rm=TRUE,p=0.05),
                   q25=quantile(x,na.rm=TRUE,p=0.25),
                   q75=quantile(x,na.rm=TRUE,p=0.75),
                   qhigh = quantile(x,na.rm=TRUE,p=0.95)) %>%
        tidyr::pivot_longer(med:qhigh, names_to = "COV", values_to = "VAL") %>%
        dplyr::mutate(term = vname) %>%
        dplyr::mutate(medcov = VAL[COV == "med"]) %>%
        dplyr::mutate(VAL2 = VAL - medcov) %>%
        dplyr::mutate(medcov = med_x) %>%
        dplyr::mutate(VAL = VAL + center_offset)
    }

    tmp_exposuredata <- dplyr::bind_rows(tmp_exposuredata)
  } else {
    tmp_exposuredata <- NULL
  }


  if(length(contcovars) > 0) {

    tmp_numdata <- list()

    for(i in seq_along(contcovars)) {
      vname <- contcovars[i]
      df_x <- tmp_data %>% dplyr::rename(x := !!vname)
      med_x <- stats::median(df_x$x)
      center_offset <- if (!is.null(centering_values) && vname %in% names(centering_values)) centering_values[[vname]] else 0

      tmp_numdata[[i]] <- df_x %>%
        dplyr::mutate(Q =  cut(x,
                        breaks=quantile(unique(df_x[,c('x','id')])$x,c(0,.25,.5,.75,1), na.rm =TRUE),
                        include.lowest=TRUE,
                        labels=paste0("Q",c(1,2,3,4)))) %>%
        dplyr::summarise( med = stats::median(x,na.rm=TRUE),
                   qlow =quantile(x,na.rm=TRUE,p=0.05),
                   q25=quantile(x,na.rm=TRUE,p=0.25),
                   q75=quantile(x,na.rm=TRUE,p=0.75),
                   qhigh = quantile(x,na.rm=TRUE,p=0.95)) %>%
        tidyr::pivot_longer(med:qhigh, names_to = "COV", values_to = "VAL") %>%
        dplyr::mutate(term = vname) %>%
        dplyr::mutate(medcov = VAL[COV == "med"]) %>%
        dplyr::mutate(VAL2 = VAL - medcov) %>%
        dplyr::mutate(medcov = med_x) %>%
        dplyr::mutate(VAL = VAL + center_offset)
    }

    tmp_numdata <- dplyr::bind_rows(tmp_numdata)
  } else {
    tmp_numdata <- NULL
  }

  if(length(catcovars) > 0) {
    tmp_catdata <- tmp_data %>%
      dplyr::select(dplyr::all_of(catcovars)) %>%
      dplyr::distinct() %>%
      tidyr::pivot_longer(dplyr::all_of(catcovars),
                          names_to = "term",
                          values_to = "COV") %>%
      dplyr::distinct() %>%
      dplyr::mutate(term = paste0(term, COV))
  } else {
    tmp_catdata <- NULL
  }

  tmp_dat <- dplyr::bind_rows(list(tmp_exposuredata, tmp_numdata, tmp_catdata))

  # Add placeholder entries for interaction terms so their broom::tidy() rows
  # are included in the forest plot as a separate "Interactions" group.
  # For models with interactions AND accessible vcov, compute proper conditional effects.
  # Resolve variable metadata (var_spec > naming_df > centering_values > NULL)
  var_meta <- resolve_var_metadata(fit, namings)

  interaction_fp_data <- NULL
  if (length(interaction_terms) > 0) {
    # Try to compute quantile-based conditional effects using vcov
    int_result <- tryCatch(
      compute_interaction_forest_data(
        fit = fit,
        interaction_terms = interaction_terms,
        quantiles_int = quantiles_interaction,
        modifier_quantiles = modifier_quantiles,
        ci_level = ci,
        var_spec = if (!is.null(fit$var_spec)) fit$var_spec else NULL
      ),
      error = function(e) NULL
    )

    if (!is.null(int_result) && nrow(int_result) > 0) {
      # Store properly computed interaction effects for later
      interaction_fp_data <- int_result
    } else {
      # Fallback: use raw exp(beta) as placeholder
      tmp_intdata <- tibble::tibble(
        term = interaction_terms,
        COV  = interaction_terms,
        VAL  = NA_real_,
        VAL2 = NA_real_,
        medcov = NA_real_
      )
      tmp_dat <- dplyr::bind_rows(tmp_dat, tmp_intdata)
    }
  }

  if (nrow(tmp_dat) < 1 && is.null(interaction_fp_data)) {
    cli::cli_abort("Cannot build forest plot: no covariates found in the fit.")
  }

  if(all(class(fit$final_model)  %in% "coxph")) {
    parname='Hazard Ratio'
    xlab_name <- paste(parname, fit$metadata$dvid_name,"and 95% CI", sep = " ")
  } else {
    parname='Odds Ratio'
    xlab_name <- paste(parname, fit$metadata$dvid_name,"and 95% CI", sep = " ")
  }

  # Process main-effect data through broom::tidy pipeline
  if (nrow(tmp_dat) > 0) {
    mod1_raw <- broom::tidy(fit$final_model)

    # Part B: when centering was NOT applied, fix main-effect rows for variables
    # that participate in interactions by computing conditional effects at modifier median.
    if (is.null(centering_values) && length(interaction_terms) > 0) {
      coefs_all <- tryCatch(stats::coef(fit$final_model), error = function(e) NULL)
      V_all <- tryCatch(stats::vcov(fit$final_model), error = function(e) NULL)
      data_all <- fit$model_data

      if (!is.null(coefs_all) && !is.null(V_all)) {
        coef_names_all <- names(coefs_all)
        z_val <- stats::qnorm(1 - (1 - ci) / 2)

        for (int_t in interaction_terms) {
          parents <- get_parent_effects(int_t)
          # Only act on continuous parents (categorical handled separately)
          cont_parents <- parents[vapply(parents, function(p) {
            col <- data_all[[p]]
            !is.factor(col) && !is.character(col) && !is.logical(col)
          }, logical(1))]

          for (eff_var in cont_parents) {
            modifier_vars <- setdiff(cont_parents, eff_var)
            if (length(modifier_vars) == 0) next

            # Build interaction term name(s) (R may reverse order)
            int_coef_candidates <- vapply(modifier_vars, function(mv) {
              c1 <- paste0(eff_var, ":", mv)
              c2 <- paste0(mv, ":", eff_var)
              if (c1 %in% coef_names_all) c1
              else if (c2 %in% coef_names_all) c2
              else NA_character_
            }, character(1))
            int_coef_names <- int_coef_candidates[!is.na(int_coef_candidates)]
            if (length(int_coef_names) == 0) next

            mod_medians <- vapply(modifier_vars[!is.na(int_coef_candidates)], function(mv) {
              stats::median(data_all[[mv]], na.rm = TRUE)
            }, numeric(1))

            cond <- tryCatch(
              compute_conditional_effect(coefs_all, V_all, eff_var, int_coef_names, mod_medians),
              error = function(e) NULL
            )
            if (is.null(cond)) next

            # Replace the estimate and std.error for this term in mod1_raw
            idx <- which(mod1_raw$term == eff_var)
            if (length(idx) == 1) {
              mod1_raw$estimate[idx] <- cond$beta
              mod1_raw$std.error[idx] <- cond$se
            }
          }
        }
      }
    }

    mod1 <- mod1_raw %>%
      dplyr::mutate(CI95low=estimate+qnorm(.025)*std.error,
             CI95hig=estimate+qnorm(.975)*std.error) %>%
      tidyr::pivot_longer(cols=c(CI95low,CI95hig,estimate),
                          names_to='CI',values_to='CI_VAL')

    fp_data <- dplyr::left_join(mod1,tmp_dat,by=c('term'),relationship =
                           "many-to-many")
  
  ## in no numeric variables scenario
  if("VAL" %!in% colnames(fp_data)) fp_data$VAL <- rep(NA_real_, nrow(fp_data))
  if("VAL2" %!in% colnames(fp_data)) fp_data$VAL2 <- rep(NA_real_, nrow(fp_data))
  if("medcov" %!in% colnames(fp_data)) fp_data$medcov <- rep(NA_real_, nrow(fp_data))
  
  fp_data <- fp_data %>%
    dplyr::mutate(CIFIN=dplyr::case_when(!is.na(VAL2 )~ exp(CI_VAL*VAL2),
                           TRUE ~ exp(CI_VAL))) %>%
    dplyr::select(term,COV,CI,CIFIN,VAL,medcov) %>%
    tidyr::pivot_wider(names_from = c(CI),values_from=CIFIN)%>%
    dplyr::mutate(paramname=parname) %>%
    dplyr::rename(covname=term,label=COV,mid=estimate,lower=CI95low,upper=CI95hig)%>%
    dplyr::mutate(lower2=dplyr::case_when(lower>upper~upper,
                            TRUE~lower))%>%
    dplyr::mutate(upper2=dplyr::case_when(upper<lower~lower,
                            TRUE~upper))%>%
    dplyr::select(-upper,-lower)
  
  # Apply variable labels and back-transformation using resolved metadata
  fp_data <- apply_var_labels(fp_data, var_meta)
  
  fp_data <- fp_data %>% dplyr::filter(!is.na(covname)) ## drop intercept from glm model
  
  features <- fp_data %>% dplyr::select(covname) %>% dplyr::distinct() %>% dplyr::pull()
  # Filter features to actual data columns only (exclude interaction term names)
  features_cols <- features[features %in% colnames(tmp_data)]

  fp_data <- fp_data %>% dplyr::select(group = display_covname, covname, group_level, mid, lo = lower2, hi = upper2)

  if(length(catcovars) > 0) {

    feature_ds <- tmp_data %>%
      dplyr::select(dplyr::all_of(features_cols))%>%
      dplyr::select(dplyr::all_of(catcovars)) %>%
      dplyr::mutate(id = 1:dplyr::n()) %>%
      tidyr::pivot_longer(-id, names_to = "covname", values_to = "group_level") %>%
      dplyr::group_by(covname, group_level) %>%
      dplyr::count()

    categorical_fp_data <- fp_data %>% dplyr::filter(covname %in% catcovars)

    ref_ds <- feature_ds %>%
      dplyr::anti_join(categorical_fp_data, by = c("covname", "group_level")) %>%
      dplyr::mutate(mid = 1, lo = 1, hi = 1) %>%
      dplyr::left_join(
        dplyr::distinct(fp_data %>% dplyr::filter(covname %in% catcovars), covname, group),
        by = "covname"
      ) %>%
      dplyr::mutate(group = dplyr::coalesce(group, covname)) %>%
      dplyr::select(group, covname, group_level, mid, lo, hi, n)

    cat_ds <- dplyr::bind_rows(ref_ds,
                        categorical_fp_data %>% dplyr::inner_join(feature_ds, by = c("covname", "group_level")) %>%
                          dplyr::select(group, covname, group_level, mid, lo, hi, n)) %>%
      dplyr::arrange(covname, mid)

    if(n_thres > 0) {
      cat_ds <- cat_ds %>% dplyr::filter(n > n_thres)
    }

    if(!is.null(add_lev)){
      cat_ds <- cat_ds %>% dplyr::mutate(group_level = paste(group_level, add_lev))
    }

    cat_ds <- cat_ds %>%
      dplyr::ungroup() %>%
      dplyr::mutate(group_level = paste0(group_level, "\nN=", n))

    lhs <- dplyr::distinct(cat_ds, covname)

    fp_data <- dplyr::anti_join(fp_data, lhs, by = "covname") %>% dplyr::bind_rows(cat_ds)

  }

  if(length(exposurecovars) > 0) {
    exposure_ds <- fp_data %>%
      dplyr::filter(covname %in% exposurecovars)

    exposure_ds <- exposure_ds %>% dplyr::mutate(group_level = paste0(c("Q0.50", "Q0.05", "Q0.25",  "Q0.75", "Q0.95"), "; ", group_level )) %>%
      dplyr::mutate(n = c(0.5, 0.05, 0.25, 0.75, 0.95) ) %>%
      dplyr::distinct(covname, group_level, .keep_all = TRUE)
    fp_data <- fp_data %>% dplyr::anti_join(dplyr::distinct(exposure_ds, covname), by = "covname") %>% dplyr::bind_rows(exposure_ds)
  }


  if(length(contcovars) > 0) {
    cont_ds <- fp_data %>%
      dplyr::filter(covname %in% contcovars)

    cont_ds <- cont_ds %>%
      dplyr::mutate(group_level =  paste0(c("Q0.50", "Q0.05", "Q0.25",  "Q0.75", "Q0.95"), "; ", group_level) ) %>%
      dplyr::mutate(n = rep( c(0.5, 0.05, 0.25, 0.75, 0.95), length(contcovars) )  ) %>%
      dplyr::distinct(covname, group_level, .keep_all = TRUE)
    fp_data <- fp_data %>% dplyr::filter(!(covname %in% contcovars)) %>% dplyr::bind_rows(cont_ds)
  }

  } else {
    # No main-effect terms -- fp_data starts empty, interactions will be appended below
    fp_data <- tibble::tibble(group = character(), covname = character(),
                              group_level = character(),
                              mid = numeric(), lo = numeric(), hi = numeric())
  }

  # Append properly computed interaction effects (replacing any placeholder rows)
  if (!is.null(interaction_fp_data)) {
    # Remove any fallback placeholder interaction rows (group name contains ":" from broom)
    fp_data <- fp_data %>% dplyr::filter(!grepl(":", group, fixed = TRUE))
    # Add the computed conditional effects with sort order
    interaction_fp_data <- interaction_fp_data %>% dplyr::mutate(n = dplyr::row_number() / (dplyr::n() + 1))
    fp_data <- dplyr::bind_rows(fp_data, interaction_fp_data)
  }

  # Tag term type for visual styling downstream
  # Interaction groups contain x in the group name (from compute_interaction_forest_data)
  fp_data <- fp_data %>%
    dplyr::mutate(term_type = dplyr::if_else(
      grepl("\u00d7", group, fixed = TRUE) | grepl(":", group, fixed = TRUE),
      "interaction", "main"
    ))

  # When adjust_for is provided, subdivide "main" rows into "adjusted" and "selected"
  if (!is.null(adjust_for) && length(adjust_for) > 0) {
    has_covname <- "covname" %in% colnames(fp_data)
    fp_data <- fp_data %>%
      dplyr::mutate(term_type = dplyr::case_when(
        term_type == "interaction" ~ "interaction",
        has_covname & (.data$covname %in% adjust_for) ~ "adjusted",
        !has_covname & (.data$group %in% adjust_for) ~ "adjusted",
        term_type == "main" ~ "selected",
        TRUE ~ term_type
      ))

    if (show_terms != "all") {
      fp_data <- fp_data[fp_data$term_type == show_terms | fp_data$term_type == "interaction", ]
    }
  }

  # Ensure n column exists for sorting (fallback interaction rows may lack it)
  if (!"n" %in% colnames(fp_data)) fp_data$n <- NA_real_
  fp_data <- fp_data %>% dplyr::mutate(n = dplyr::if_else(is.na(n), 0.5, n))

  # Order: adjusted first, then selected, then interaction (then legacy "main" as fallback)
  type_order <- c("adjusted", "selected", "interaction", "main")
  fp_data <- fp_data %>%
    dplyr::mutate(.sort_type = match(term_type, type_order)) %>%
    dplyr::arrange(.data$.sort_type, group, n) %>%
    dplyr::select(-".sort_type", -n, -dplyr::any_of("covname"))

  if(!plot) {
    return(fp_data)
  }

  fp1 <- plot_grouped_ratios(fp_data, xlab = xlab_name, log_x = log_x,
                             adjust_colour = adjust_colour,
                             selected_colour = selected_colour,
                             adjust_shape = adjust_shape,
                             selected_shape = selected_shape)

  return(fp1)

}


#' Compute quantile-based conditional effects for interaction terms
#'
#' Internal helper that computes hazard ratios (or odds ratios) for interaction
#' terms at specified quantile values, using the variance-covariance matrix for
#' proper confidence intervals.
#'
#' @param fit The fitted SCM result object.
#' @param interaction_terms Character vector of interaction terms (e.g., "CMAXss:AGE").
#' @param quantiles_int Numeric vector of quantile probabilities for the "effect variable".
#' @param modifier_quantiles Numeric vector of quantile probabilities for the "modifier variable"
#'   (continuous x continuous case only).
#' @param ci_level Confidence level (e.g., 0.95).
#'
#' @return A tibble with columns: group, group_level, mid, lo, hi.
#'
#' @keywords internal
compute_interaction_forest_data <- function(fit, interaction_terms, quantiles_int,
                                            modifier_quantiles, ci_level, var_spec = NULL) {
  model <- fit$final_model
  data <- fit$model_data
  V <- stats::vcov(model)
  coefs <- stats::coef(model)
  z <- stats::qnorm(1 - (1 - ci_level) / 2)
  coef_names <- names(coefs)
  exposure_vars <- fit$metadata$candidate_exposure_vars

  results <- list()


  for (int_term in interaction_terms) {
    parents <- get_parent_effects(int_term)

    # Classify parent variables as continuous or categorical
    parent_types <- vapply(parents, function(p) {
      col <- data[[p]]
      if (is.factor(col) || is.character(col) || is.logical(col)) "categorical" else "continuous"
    }, character(1))

    n_cat <- sum(parent_types == "categorical")
    n_cont <- sum(parent_types == "continuous")

    if (n_cont == 2 && n_cat == 0) {
      # ---- Case A: Continuous x Continuous ----
      res <- .fp_int_cont_cont(
        parents = parents, int_term = int_term, coefs = coefs,
        V = V, data = data, exposure_vars = exposure_vars,
        quantiles_int = quantiles_int, modifier_quantiles = modifier_quantiles, z = z,
        var_spec = var_spec
      )
      results <- c(results, list(res))

    } else if (n_cont == 1 && n_cat == 1) {
      # ---- Case B: Continuous x Categorical ----
      cont_var <- parents[parent_types == "continuous"]
      cat_var <- parents[parent_types == "categorical"]
      res <- .fp_int_cont_cat(
        cont_var = cont_var, cat_var = cat_var, int_term = int_term,
        coefs = coefs, coef_names = coef_names, V = V, data = data,
        quantiles_int = quantiles_int, z = z, var_spec = var_spec
      )
      results <- c(results, list(res))

    } else {
      # ---- Case C: Categorical x Categorical (or other) ----
      # Fall back to raw exp(beta) with CI from vcov
      res <- .fp_int_cat_cat(
        int_term = int_term, coefs = coefs, coef_names = coef_names, V = V, z = z,
        var_spec = var_spec
      )
      results <- c(results, list(res))
    }
  }

  out <- dplyr::bind_rows(results)
  if (nrow(out) == 0) return(NULL)
  out
}


#' Conditional effects for continuous x continuous interaction
#' @keywords internal
#' @noRd
.fp_int_cont_cont <- function(parents, int_term, coefs, V, data,
                              exposure_vars, quantiles_int, modifier_quantiles, z,
                              var_spec = NULL) {
  # Choose effect variable: prefer exposure variable if present
  effect_var <- intersect(parents, exposure_vars)[1]
  if (is.na(effect_var)) effect_var <- parents[1]
  modifier_var <- setdiff(parents, effect_var)
  effect_label   <- get_var_label(var_spec, effect_var)
  modifier_label <- get_var_label(var_spec, modifier_var)

  beta_effect_name <- effect_var
  beta_int_name <- int_term

  # Check coefficient order -- R may store as "A:B" or "B:A"
  if (!beta_int_name %in% names(coefs)) {
    beta_int_name <- paste(rev(get_parent_effects(int_term)), collapse = ":")
  }
  if (!beta_int_name %in% names(coefs)) return(tibble::tibble())

  # Quantile values from data
  eff_quantiles <- stats::quantile(data[[effect_var]], probs = quantiles_int, na.rm = TRUE)
  eff_median <- stats::median(data[[effect_var]], na.rm = TRUE)
  mod_values <- stats::quantile(data[[modifier_var]], probs = modifier_quantiles, na.rm = TRUE)

  rows <- list()
  for (mi in seq_along(mod_values)) {
    x2_val <- mod_values[mi]
    mod_label <- paste0("Q", format(modifier_quantiles[mi], nsmall = 2))

    # Conditional coefficient: beta_effect + beta_int * x2_val
    cond_beta <- coefs[beta_effect_name] + coefs[beta_int_name] * x2_val
    # Variance of conditional coefficient
    var_cond <- V[beta_effect_name, beta_effect_name] +
      x2_val^2 * V[beta_int_name, beta_int_name] +
      2 * x2_val * V[beta_effect_name, beta_int_name]
    se_cond <- sqrt(max(var_cond, 0))

    for (qi in seq_along(eff_quantiles)) {
      q_val <- eff_quantiles[qi]
      q_label <- paste0("Q", format(quantiles_int[qi], nsmall = 2))
      contrast <- q_val - eff_median

      log_effect <- unname(cond_beta * contrast)
      se_effect <- abs(contrast) * se_cond

      rows[[length(rows) + 1]] <- tibble::tibble(
        group = paste0(effect_label, " \u00d7 ", modifier_label),
        group_level = paste0(q_label, " | ", modifier_label, "=", mod_label),
        mid = exp(log_effect),
        lo = exp(log_effect - z * se_effect),
        hi = exp(log_effect + z * se_effect)
      )
    }
  }
  dplyr::bind_rows(rows)
}


#' Conditional effects for continuous x categorical interaction
#' @keywords internal
#' @noRd
.fp_int_cont_cat <- function(cont_var, cat_var, int_term, coefs, coef_names,
                             V, data, quantiles_int, z, var_spec = NULL) {
  beta_cont_name <- cont_var
  cont_label <- get_var_label(var_spec, cont_var)
  cat_label  <- get_var_label(var_spec, cat_var)

  # Get factor levels
  col <- data[[cat_var]]
  if (is.factor(col)) {
    all_levels <- levels(col)
  } else {
    all_levels <- sort(unique(col))
  }
  ref_level <- all_levels[1]
  other_levels <- all_levels[-1]

  # Quantile values for the continuous variable
  eff_quantiles <- stats::quantile(data[[cont_var]], probs = quantiles_int, na.rm = TRUE)
  eff_median <- stats::median(data[[cont_var]], na.rm = TRUE)

  rows <- list()

  # --- At reference level: effect is just beta_cont ---
  for (qi in seq_along(eff_quantiles)) {
    q_val <- eff_quantiles[qi]
    q_label <- paste0("Q", format(quantiles_int[qi], nsmall = 2))
    contrast <- q_val - eff_median

    log_effect <- unname(coefs[beta_cont_name] * contrast)
    se_effect <- abs(contrast) * sqrt(max(V[beta_cont_name, beta_cont_name], 0))

    rows[[length(rows) + 1]] <- tibble::tibble(
      group = paste0(cont_label, " \u00d7 ", cat_label),
      group_level = paste0(q_label, " | ", cat_label, "=", ref_level, " [ref]"),
      mid = exp(log_effect),
      lo = exp(log_effect - z * se_effect),
      hi = exp(log_effect + z * se_effect)
    )
  }

  # --- At non-reference levels ---
  for (lev in other_levels) {
    # Find the interaction coefficient name (R may encode as "cont:catLevel" or "catLevel:cont")
    beta_int_name <- .find_interaction_coef(cont_var, cat_var, lev, coef_names)
    if (is.null(beta_int_name)) next

    # Conditional coefficient: beta_cont + beta_int
    cond_beta <- coefs[beta_cont_name] + coefs[beta_int_name]
    var_cond <- V[beta_cont_name, beta_cont_name] +
      V[beta_int_name, beta_int_name] +
      2 * V[beta_cont_name, beta_int_name]
    se_cond <- sqrt(max(var_cond, 0))

    for (qi in seq_along(eff_quantiles)) {
      q_val <- eff_quantiles[qi]
      q_label <- paste0("Q", format(quantiles_int[qi], nsmall = 2))
      contrast <- q_val - eff_median

      log_effect <- unname(cond_beta * contrast)
      se_effect <- abs(contrast) * se_cond

      rows[[length(rows) + 1]] <- tibble::tibble(
        group = paste0(cont_label, " \u00d7 ", cat_label),
        group_level = paste0(q_label, " | ", cat_label, "=", lev),
        mid = exp(log_effect),
        lo = exp(log_effect - z * se_effect),
        hi = exp(log_effect + z * se_effect)
      )
    }
  }
  dplyr::bind_rows(rows)
}


#' Raw exponentiated coefficients for categorical x categorical interactions
#' @keywords internal
#' @noRd
.fp_int_cat_cat <- function(int_term, coefs, coef_names, V, z, var_spec = NULL) {
  # Find all coefficient names matching this interaction
  parents <- get_parent_effects(int_term)
  int_coef_names <- .find_all_interaction_coefs(parents, coef_names)

  if (length(int_coef_names) == 0) return(tibble::tibble())

  p1_label <- get_var_label(var_spec, parents[1])
  p2_label <- get_var_label(var_spec, parents[2])

  rows <- list()
  for (cname in int_coef_names) {
    beta_val <- coefs[cname]
    se_val <- sqrt(max(V[cname, cname], 0))

    rows[[length(rows) + 1]] <- tibble::tibble(
      group = paste0(p1_label, " \u00d7 ", p2_label),
      group_level = cname,
      mid = exp(unname(beta_val)),
      lo = exp(unname(beta_val) - z * se_val),
      hi = exp(unname(beta_val) + z * se_val)
    )
  }
  dplyr::bind_rows(rows)
}


#' Find the coefficient name for a continuous x categorical interaction level
#' @keywords internal
#' @noRd
.find_interaction_coef <- function(cont_var, cat_var, level, coef_names) {
  # Try both orderings: "cont:catLevel" and "catLevel:cont"
  candidate1 <- paste0(cont_var, ":", cat_var, level)
  candidate2 <- paste0(cat_var, level, ":", cont_var)

  if (candidate1 %in% coef_names) return(candidate1)
  if (candidate2 %in% coef_names) return(candidate2)

  # Fuzzy match: look for coefficient containing both variable references with ":"
  matches <- coef_names[grepl(":", coef_names, fixed = TRUE)]
  for (m in matches) {
    parts <- strsplit(m, ":")[[1]]
    if (length(parts) == 2) {
      if ((startsWith(parts[1], cont_var) && startsWith(parts[2], cat_var)) ||
          (startsWith(parts[1], cat_var) && startsWith(parts[2], cont_var))) {
        # Check if the level suffix matches
        if (grepl(paste0(cat_var, level), m, fixed = TRUE) ||
            grepl(paste0(level), m, fixed = TRUE)) {
          return(m)
        }
      }
    }
  }
  NULL
}


#' Find all coefficient names for an interaction between two variables
#' @keywords internal
#' @noRd
.find_all_interaction_coefs <- function(parents, coef_names) {
  matches <- coef_names[grepl(":", coef_names, fixed = TRUE)]
  result <- character(0)
  for (m in matches) {
    parts <- strsplit(m, ":")[[1]]
    if (length(parts) == 2) {
      if (all(vapply(parents, function(p) any(startsWith(parts, p)), logical(1)))) {
        result <- c(result, m)
      }
    }
  }
  result
}


#' Plot grouped ratios with CIs and an adjacent table
#'
#' Creates a forest plot with per-group panels (one panel per covariate or
#' interaction term). Interaction panels are visually distinguished with diamond
#' markers, italic titles, and independent x-axis scaling.
#'
#' @param df Data frame with columns: group, group_level, mid, lo, hi.
#'   Optionally includes a `term_type` column with values "main", "interaction",
#'   "adjusted", or "selected".
#' @param xlab X-axis label for ratio (default "Ratio")
#' @param title Optional overall title (character or NULL)
#' @param xlim Optional numeric length-2 vector for x-axis limits; if NULL, auto per group
#' @param ref_line Reference ratio to highlight (default = 1)
#' @param point_size Size of points (default 2.5)
#' @param ci_width Line width for CIs (default 0.7)
#' @param theme_base A ggplot2 theme function (default theme_bw)
#' @param log_x Logical; if TRUE, use log10 x-axis scale. Default FALSE.
#' @param table_digits Number of digits for mid/lo/hi in table (default 3)
#' @param table_width_ratio Numeric vector of length 2 giving relative widths of forest and table panels.
#' @param adjust_colour Fill colour for "adjusted" term markers. Defaults to \code{"grey50"}.
#' @param selected_colour Fill colour for "selected" term markers. Defaults to \code{"steelblue"}.
#' @param adjust_shape Point shape for "adjusted" term markers. Defaults to \code{17L}.
#' @param selected_shape Point shape for "selected" term markers. Defaults to \code{16L}.
#' @return A patchwork object
#' @examples
#' data(demo_eres)
#' demo_eres$SEX <- as.factor(demo_eres$SEX)
#' lr_fit <- scm_glm(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", covariates = c("AGE", "WT"),
#'   dv_col = "sim_status", significance_level = c(0.2, 0.2)
#' )
#' dat <- make_forest_plot(lr_fit, plot = FALSE)
#' plot_grouped_ratios(dat)
#' @export
plot_grouped_ratios <- function(
    df,
    xlab = "Ratio",
    title = NULL,
    xlim = NULL,
    ref_line = 1,
    point_size = 2.5,
    ci_width = 0.7,
    log_x = FALSE,
    theme_base = ggplot2::theme_bw,
    table_digits = 3,
    table_width_ratio = c(2.5, 1),
    adjust_colour = "grey50",
    selected_colour = "steelblue",
    adjust_shape = 17L,
    selected_shape = 16L
) {

  # Validate columns
  required_cols <- c("group", "group_level", "mid", "lo", "hi")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    cli::cli_abort("Missing required columns: {.val {missing_cols}}.")
  }

  # Ensure term_type column exists (backward compat)
  if (!"term_type" %in% names(df)) {
    df$term_type <- "main"
  }

  # Coerce and standardize
  df <- df %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      group = as.character(.data$group),
      group_level = as.character(.data$group_level),
      mid = as.numeric(.data$mid),
      lo  = as.numeric(.data$lo),
      hi  = as.numeric(.data$hi),
      term_type = as.character(.data$term_type)
    )

  # group_level_f is assigned per-group in the blocks loop below (single-group subsets
  # avoid vctrs level-merging that silently reorders factors when groups are combined)

  # Number formatting for table
  fmt_num <- function(x) {
    base::formatC(x, format = "f", digits = table_digits, drop0trailing = FALSE)
  }

  # Per-group x-axis limits (avoids wide interaction CIs compressing main effects)
  compute_xlim <- function(dsub) {
    rng <- range(c(dsub$lo, dsub$hi), na.rm = TRUE)
    pad <- diff(rng) * 0.10
    if (!is.finite(pad) || pad == 0) pad <- 0.1
    c(max(0, rng[1] - pad), rng[2] + pad)
  }

  # Forest plot panel for a group
  build_forest <- function(dsub, gname, show_x = TRUE, log_x = FALSE, group_xlim = NULL,
                           is_interaction = FALSE, panel_term_type = "main") {

    n_rows <- nrow(dsub)
    # Dynamically scale y-axis text size; cap at 12-row panel baseline (7.25) to avoid oversized text in small panels
    y_text_size <- min(7.25, max(6.5, 9 - max(0, n_rows - 5) * 0.25))

    # Marker styling: interactions get diamonds; adjusted/selected get per-type colour/shape
    pt_shape <- if (is_interaction) 23 else if (panel_term_type == "adjusted") adjust_shape else if (panel_term_type == "selected") selected_shape else 21
    pt_fill  <- if (is_interaction) "#2166AC" else if (panel_term_type == "adjusted") adjust_colour else if (panel_term_type == "selected") selected_colour else "black"
    ci_color <- if (is_interaction) "#4393C3" else "grey30"
    title_face <- if (is_interaction) "bold.italic" else "bold"

    p <- ggplot2::ggplot(dsub, ggplot2::aes(x = .data$mid, y = .data$group_level_f)) +
      ggplot2::geom_vline(xintercept = ref_line, linetype = "dashed", color = "grey40") +
      ggplot2::geom_errorbarh(
        ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
        #height = 0, 
        #width = 0, 
        linewidth = ci_width, color = ci_color
      ) +
      ggplot2::geom_point(shape = pt_shape, size = point_size, stroke = 0.4,
                          fill = pt_fill, color = "black") +
      ggplot2::labs(x = if (show_x) xlab else NULL, y = NULL, title = gname) +
      theme_base() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0, face = title_face, size = 10,
                                           colour = if (is_interaction) "#2166AC" else "black",
                                           margin = ggplot2::margin(b = 6)),
        axis.text.y = ggplot2::element_text(size = y_text_size),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5)
      ) +
      ggpubr::labs_pubr()

    if (log_x) {
      p <- p +
        ggplot2::scale_x_log10(breaks = scales::trans_breaks('log10', function(x) 10^x),
                               labels = fmt_dcimals(decimals = 2))
    } else {
      p <- p + ggplot2::scale_x_continuous(limits = group_xlim)
    }

    return(p)
  }

  # Table panel for a group -- uses the SAME discrete y-axis as the forest panel
  build_table <- function(dsub, show_header = TRUE) {
    n_rows <- nrow(dsub)
    # Cap at 12-row panel baseline (2.6) to prevent oversized text in small panels
    tbl_size <- min(2.6, max(2.4, 3.2 - max(0, n_rows - 6) * 0.1))
    hdr_size <- tbl_size + 0.3

    tbl_data <- dsub %>% dplyr::mutate(
      ratio_label = fmt_num(.data$mid),
      ci_label    = paste0(fmt_num(.data$lo), " \u2013 ", fmt_num(.data$hi))
    )

    p <- ggplot2::ggplot(tbl_data, ggplot2::aes(y = .data$group_level_f)) +
      ggplot2::geom_text(ggplot2::aes(x = 0.22, label = .data$ratio_label),
                         size = tbl_size, hjust = 0.5, color = "grey20") +
      ggplot2::geom_text(ggplot2::aes(x = 0.72, label = .data$ci_label),
                         size = tbl_size, hjust = 0.5, color = "grey20")

    if (show_header) {
      p <- p +
        ggplot2::annotate("text", x = 0.22, y = Inf, label = "Ratio",
                          fontface = "bold", size = hdr_size, hjust = 0.5, vjust = -0.8) +
        ggplot2::annotate("text", x = 0.72, y = Inf, label = "95% CI",
                          fontface = "bold", size = hdr_size, hjust = 0.5, vjust = -0.8)
    }

    p <- p +
      ggplot2::scale_x_continuous(limits = c(0, 1), breaks = NULL, expand = c(0, 0)) +
      ggplot2::scale_y_discrete(breaks = NULL) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 0)
      )

    return(p)
  }

  # Build per-group blocks
  groups <- df %>% dplyr::pull(.data$group) %>% unique()

  # Compute row counts for proportional panel heights
  row_counts <- purrr::map_int(groups, function(g) sum(df$group == g))
  # Ensure minimum height for small panels (at least equivalent of 2 rows)
  panel_heights <- pmax(row_counts, 2)

  blocks <- groups %>%
    purrr::map(function(g) {
      dsub <- df %>% dplyr::filter(.data$group == g)
      # Build factor levels here, on a single-group slice, so unique() is unambiguous
      dsub <- dsub %>% dplyr::mutate(
        group_level_f = factor(.data$group_level, levels = rev(unique(.data$group_level)))
      )
      is_int <- dsub$term_type[1] == "interaction"
      panel_tt <- dsub$term_type[1]
      group_xlim <- if (!is.null(xlim)) xlim else compute_xlim(dsub)

      p_forest <- build_forest(dsub, gname = g,
                               show_x = (g == dplyr::last(groups)),
                               log_x = log_x,
                               group_xlim = group_xlim,
                               is_interaction = is_int,
                               panel_term_type = panel_tt)
      p_table  <- build_table(dsub, show_header = TRUE)
      patchwork::wrap_plots(p_forest, p_table) +
        patchwork::plot_layout(widths = table_width_ratio)
    })

  # Combine all groups vertically with proportional heights
  combined <- purrr::reduce(blocks, `/`) + patchwork::plot_layout(
    nrow = length(blocks),
    heights = panel_heights,
    guides = "collect"
  )

  # Add optional global title
  if (!is.null(title)) {
    combined <- combined + patchwork::plot_annotation(
      title = title,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
    )
  }

  combined
}