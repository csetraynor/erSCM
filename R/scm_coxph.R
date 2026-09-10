#' Stepwise Covariate Modeling (SCM) for Cox Proportional Hazards Models
#'
#' Perform stepwise covariate model building for Cox proportional hazards models
#' using event-level data. The procedure evaluates candidate covariates
#' and exposure variables against specified significance thresholds, optionally
#' adjusts for prespecified factors, and returns the final model and summary
#' outputs.
#'
#' @param eres_data A data frame containing event/survival data and all candidate
#'   covariates and exposure variables. Must include columns for time-to-event
#'   and event indicator.
#' @param covariates A character vector of candidate covariate names to be
#'   considered in the SCM procedure (e.g., baseline factors, demographics).
#' @param exposure_vars A character vector of exposure-related predictor names
#'   to be evaluated (e.g., continuous exposure metrics, categories).
#' @param significance_level Numeric vector of length 2 giving the entry and
#'   retention p-value thresholds used in SCM. The first value is the *entry*
#'   criterion (e.g., 0.01) and the second is the *retention* criterion
#'   (e.g., 0.001).
#' @param dv_col Character string naming the event indicator column (typically
#'   0/1), default \code{"DV"}.
#' @param time_col Character string naming the time-to-event column, default
#'   \code{"TIME"}.
#' @param dvid Optional character to subset the analysis by \code{dvid}. 
#' @param dvid_col Optional character string naming the column in \code{eres_data}
#'   that contains the adverse event identifier used with \code{dvid}. Required
#'   if \code{dvid} is provided.
#' @param var_spec Variable metadata specification. Defines display labels, units, and categories
#'   for each variable so that metadata flows automatically to \code{make_forest_plot()},
#'   \code{print()}, and other downstream functions. Accepts three formats:
#'   \describe{
#'     \item{data.frame}{Must have columns \code{col_name} and \code{label}. Optional:
#'       \code{units}, \code{var_type}, \code{categories}.}
#'     \item{named character vector}{Names are column names; values are display labels with
#'       optional parenthetical units, e.g. \code{c(WT = "Body Weight (kg)")}.}
#'     \item{character path}{Path to a YAML file. Requires the \pkg{yaml} package.}
#'   }
#'   When \code{var_spec} is provided, \code{units_df} is ignored (with a message).
#'   Parsed spec is stored as \code{fit$var_spec}. Default \code{NULL}.
#' @param scale_to_median Logical; if TRUE, continuous exposure and covariate columns are
#'   divided by their median so that 1.0 represents the median value. Scaling factors are
#'   stored in \code{fit$var_spec} / \code{fit$naming_df} for back-transformation in
#'   \code{make_forest_plot()}. A warning is issued if both \code{scale_to_median} and
#'   \code{center} are \code{TRUE}. Default \code{FALSE}.
#' @param normalize_inputs Deprecated. Use \code{scale_to_median} instead.
#' @param units_df Deprecated. Use \code{var_spec} instead for variable metadata.
#' @param center Logical; if TRUE (default), continuous covariates and exposure variables are
#'   mean-centered at their medians before model fitting. This ensures that main-effect
#'   coefficients in interaction models represent the effect at the modifier's median rather
#'   than at modifier = 0, producing interpretable p-values and forest plot CIs. Centering
#'   does not change model predictions or likelihood ratio tests -- only the parameterization.
#'   Set to FALSE to reproduce pre-centering behavior. Centering values are stored in
#'   \code{fit$metadata$centering_values} for downstream use (e.g., back-transforming labels).
#'   Factor, character, and logical columns are never centered.
#' @param adjust_for Character vector of column names to include as fixed adjustment
#'   covariates in the base model. These terms are:
#'   \itemize{
#'     \item Included in the model from the start (before forward selection begins)
#'     \item Protected from backward elimination (cannot be removed)
#'     \item Excluded from the candidate pool for forward selection
#'   }
#'   Use this for known confounders or protocol-mandated adjustments. Default \code{NULL}
#'   (no forced adjustments; SCM starts from an intercept-only base model).
#' @param start_covariates Character vector of covariates (a subset of
#'   \code{covariates}) to pre-load into the model before forward selection
#'   begins, instead of starting from an intercept-only base. Required when
#'   \code{direction = "backward"}. Default \code{NULL} (start empty).
#' @param direction Character string controlling the stepwise search
#'   direction: one of \code{"both"} (forward selection followed by backward
#'   elimination, the default), \code{"forward"} (forward selection only), or
#'   \code{"backward"} (backward elimination only, starting from
#'   \code{start_covariates}).
#' @param nonlinear Non-linear transformation specification. Controls whether
#'   and how non-linear functional forms are tested during forward selection.
#'   Accepted values:
#'   \describe{
#'     \item{\code{NULL} (default)}{All variables enter linearly -- existing
#'       behavior is preserved exactly.}
#'     \item{\code{TRUE} or \code{"auto"}}{Apply default transformation search
#'       (\code{c("linear", "log", "fp1", "fp2")}) to all continuous variables.}
#'     \item{Character vector}{Same transform set for every continuous variable,
#'       e.g. \code{c("linear", "log", "rcs")}.}
#'     \item{Named list}{Per-variable specification,
#'       e.g. \code{list(CMAX = c("linear", "log", "fp2"), AGE = c("linear", "log"))}.}
#'   }
#' @param linearity_alpha Numeric. Significance threshold for preferring a
#'   non-linear form over linear within the closed-testing step. Default \code{0.05}.
#' @param fp_powers Numeric vector. Candidate powers for the Fractional
#'   Polynomial search. Default \code{c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)}.
#' @param spline_df Integer. Degrees of freedom for restricted cubic splines.
#'   Default \code{3}.
#' @param verbose Logical; if \code{TRUE}, prints progress messages and
#'   intermediate model summaries during SCM. Default \code{FALSE}.
#' @param return_models Logical; if \code{TRUE}, returns fitted model objects
#'   for the intermediate steps in addition to summary outputs.
#'   Default \code{FALSE}.
#' @param interactions Character vector or list specifying interaction terms to test.
#'   - If a character vector: each element is an interaction term in R formula notation
#'     using ":" (e.g., \code{c("age:sex", "AUC:weight")}). These are tested as-is.
#'   - If a list of character vectors: each vector defines a set of variables from which
#'     all pairwise interactions are generated (e.g., \code{list(c("AUC", "age", "sex"))} generates
#'     "AUC:age", "AUC:sex", "age:sex").
#'   - If NULL (default): no interaction testing is performed (backward compatible).
#' @param interaction_strategy Character string controlling when interactions are tested
#'   during forward selection. One of:
#'   - \code{"after_main"} (default): interactions are tested in a separate forward loop AFTER
#'     no more main effects can enter the model.
#'   - \code{"simultaneous"}: interactions are tested alongside main effects at each forward step,
#'     subject to the hierarchy constraint.
#' @param p_interaction Numeric vector of length 1 or 2 specifying significance thresholds
#'   for interaction terms during BOTH forward inclusion AND backward exclusion. If length 1,
#'   the same value is used for both inclusion and exclusion. If length 2, the first is for
#'   inclusion (forward selection) and the second for exclusion (backward elimination).
#'   If NULL (default), uses the same \code{significance_level} values as main effects.
#' @param max_interaction_order Integer. Maximum order of interactions to generate when
#'   \code{interactions} is provided as a list. Default 2 (pairwise only).
#' @param parallel Logical; if \code{TRUE}, evaluates candidate terms in
#'   parallel where applicable. Default \code{FALSE}.
#' @param ... Additional arguments passed to \code{survival::coxph()} (e.g.,
#'   \code{ties}, \code{weights}, \code{cluster}, \code{robust}).
#'
#' @details
#' The SCM algorithm typically proceeds by:
#' \enumerate{
#'   \item Optional subsetting based on \code{dvid} and \code{dvid_col}.
#'   \item Fitting a base Cox model with \code{adjust_for} forced in.
#'   \item Testing candidate \code{covariates} and \code{exposure_vars} for
#'         model entry using the first \code{significance_level} threshold.
#'   \item Refining the model by removing non-significant terms using the second
#'         \code{significance_level} threshold for retention.
#' }
#' P-values are typically derived from Wald or likelihood ratio tests as
#' implemented by \code{survival::coxph}. Users should verify proportional
#' hazards assumptions and model diagnostics outside of the SCM procedure.
#'
#' @return An object of class \code{c("eres_scm_coxph", "eres_scm")}: a list containing:
#' \itemize{
#'   \item \code{final_model}: the selected \code{coxph} model object, or \code{NULL}
#'     if no exposure variable reached the entry significance threshold.
#'   \item \code{path}: the selection path (ordered steps with added/dropped covariates)
#'   \item \code{summary}: a tidy summary table of parameter estimates and fit criteria
#'   \item \code{model_data}: the (possibly centered/scaled) data.frame used for fitting
#'   \item \code{metadata}: run metadata including \code{dvid_name}, \code{centering_values},
#'     and \code{scaling_method} (\code{"none"}, \code{"centered"}, \code{"median_scaled"}, or
#'     \code{"median_scaled_and_centered"})
#'   \item \code{var_spec}: canonical variable metadata tibble (from \code{var_spec} argument);
#'     \code{NULL} when not provided
#'   \item \code{naming_df}: legacy scaling tibble (for backward compat); use \code{var_spec} instead
#' }
#'
#' @section Requirements:
#' - Columns \code{time_col} and \code{dv_col} must exist in \code{eres_data}.
#' - Variables in \code{covariates}, \code{exposure_vars}, and
#'   \code{adjust_for} must be present in \code{eres_data}.
#' - If \code{dvid} is provided, \code{dvid_col} must also be provided.
#'
#' @section Notes:
#' When \code{nonlinear} is non-\code{NULL} AND \code{interactions} is also specified,
#' interactions whose parent variable entered in a non-linear form are silently
#' skipped. The eligibility check compares bare variable names against
#' \code{Relation}, which stores formula fragments (e.g., \code{log(CMAX)}) rather
#' than bare names (\code{CMAX}), so the interaction is never found eligible.
#' A warning is emitted when both arguments are non-\code{NULL}.
#'
#' \strong{Workaround:} Use \code{nonlinear = NULL} when interaction testing is
#' required, or pre-transform the column and set \code{nonlinear = NULL}.
#'
#' @seealso \code{survival::coxph}, \code{survival::Surv}
#'
#' @examples
#' data(demo_eres)
#' demo_eres$SEX <- as.factor(demo_eres$SEX)
#' \donttest{
#' res <- scm_coxph(
#'   eres_data = demo_eres,
#'   covariates = c("AGE", "WT"),
#'   exposure_vars = "CMAXss",
#'   significance_level = c(0.05, 0.01),
#'   dv_col = "sim_status",
#'   time_col = "sim_time"
#' )
#' res$summary
#' }
#'
#' @export
scm_coxph <- function(
    eres_data,
    covariates,
    exposure_vars,
    significance_level = c(0.01, 0.001),
    dvid = NULL,
    dvid_col = NULL,
    dv_col = "DV",
    time_col = "TIME",
    var_spec = NULL,
    scale_to_median = FALSE,
    normalize_inputs = FALSE,
    units_df = NULL,
    center = TRUE,
    adjust_for = NULL,
    start_covariates = NULL,
    direction = "both",
    interactions = NULL,
    interaction_strategy = "after_main",
    p_interaction = NULL,
    max_interaction_order = 2,
    nonlinear = NULL,
    linearity_alpha = 0.05,
    fp_powers = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3),
    spline_df = 3L,
    verbose = FALSE,
    return_models = FALSE,
    parallel = FALSE,
    ...)
{
  # Validate adjust_for
  if (!is.null(adjust_for)) {
    if (!is.character(adjust_for)) {
      cli::cli_abort("{.arg adjust_for} must be a character vector of column names or NULL.")
    }
    missing_cols <- setdiff(adjust_for, colnames(eres_data))
    if (length(missing_cols) > 0) {
      cli::cli_abort(c(
        "The following {.arg adjust_for} variables are not found in the data:",
        "x" = "{.val {missing_cols}}"
      ))
    }
  }

  # Validate start_covariates
  if (!is.null(start_covariates)) {
    if (!is.character(start_covariates)) {
      cli::cli_abort("{.arg start_covariates} must be a character vector or NULL.")
    }
    missing_sc <- setdiff(start_covariates, covariates)
    if (length(missing_sc) > 0) {
      cli::cli_abort(c(
        "The following {.arg start_covariates} are not in {.arg covariates}:",
        "x" = "{.val {missing_sc}}"
      ))
    }
    overlap_sc <- intersect(start_covariates, adjust_for)
    if (length(overlap_sc) > 0) {
      cli::cli_abort(c(
        "The following variables appear in both {.arg start_covariates} and {.arg adjust_for}:",
        "x" = "{.val {overlap_sc}}"
      ))
    }
  }

  # Validate direction
  direction <- match.arg(direction, c("both", "forward", "backward"))
  if (direction == "backward" && is.null(start_covariates)) {
    cli::cli_abort("{.code direction = \"backward\"} requires {.arg start_covariates} to be provided.")
  }

  # Handle deprecated parameters
  if (!missing(normalize_inputs) && normalize_inputs) {
    warning("`normalize_inputs` is deprecated. Use `scale_to_median` instead.", call. = FALSE)
    scale_to_median <- normalize_inputs
  }
  if (!missing(units_df) && !is.null(units_df)) {
    if (!is.null(var_spec)) {
      message("`units_df` is ignored when `var_spec` is provided. ",
              "Variable units should be specified in `var_spec`.")
      units_df <- NULL
    } else {
      warning("`units_df` is deprecated. Use `var_spec` instead for variable metadata.",
              call. = FALSE)
      var_spec <- data.frame(
        col_name = units_df$covname,
        label    = units_df$covname,
        units    = units_df$units
      )
    }
  }

  .check_inputs_surv(modds = eres_data, dvid_col = dvid_col, dvid = dvid, dv_col = dv_col,
                time_col = time_col, covariates = c(covariates, exposure_vars))

  # Parse and validate variable specification
  var_spec_parsed <- parse_var_spec(var_spec)
  if (!is.null(var_spec_parsed)) {
    var_spec_parsed <- validate_var_spec(var_spec_parsed, eres_data, covariates, exposure_vars)
  }

  ## Expand and validate interaction terms
  candidate_interactions <- expand_interactions(interactions, max_interaction_order)
  if (length(candidate_interactions) > 0) {
    validate_interactions(candidate_interactions, covariates, exposure_vars)
    interaction_strategy <- match.arg(interaction_strategy, c("after_main", "simultaneous"))
  }
  # Resolve interaction p-value thresholds
  p_incl_int <- if (!is.null(p_interaction)) p_interaction[1] else significance_level[1]
  p_excl_int <- if (!is.null(p_interaction) && length(p_interaction) >= 2) p_interaction[2] else significance_level[2]

  # --- Parse nonlinear specification (on raw data before scaling/centering) ---
  nl_spec       <- parse_nonlinear_spec(nonlinear, covariates, exposure_vars, eres_data)
  use_nonlinear <- !is.null(nl_spec)

  # Warn when nonlinear + interactions are both requested: interactions whose parent
  # variable entered in a non-linear form (e.g. log(CMAX)) are silently ineligible
  # because the eligibility check matches bare names against formula fragments.
  if (use_nonlinear && length(candidate_interactions) > 0L) {
    cli::cli_warn(
      c(
        "!" = "Both {.arg nonlinear} and {.arg interactions} are specified.",
        "i" = "Interaction terms whose parent variable enters in a non-linear form
               (e.g. {.code log(CMAX)}) will be silently skipped by the eligibility
               check, because {.code Relation} stores formula fragments rather than
               bare variable names.",
        "i" = "See {.code ?scm_coxph} section {.emph Notes} for details and workarounds."
      )
    )
  }

  if (scale_to_median && center) {
    warning(
      "Both `scale_to_median` and `center` are TRUE. The combined transformation is:\n",
      "  x_final = x/median(x) - 1\n",
      "Coefficients represent the log-HR per unit change on the fractional-deviation ",
      "scale (i.e., per 100% increase relative to median). Consider using only one:\n",
      "- `scale_to_median = TRUE, center = FALSE`: coeff = effect per fold-change\n",
      "- `scale_to_median = FALSE, center = TRUE`: coeff = per-unit effect at modifier median (recommended for interaction models)",
      call. = FALSE
    )
  }

  if (scale_to_median) {
    scale_cols <- c(covariates, exposure_vars)
    scaling_result <- apply_median_scaling(eres_data, scale_cols, units_df = NULL, is_centered = center)
    eres_data  <- scaling_result$data
    naming_ds  <- scaling_result$naming_df
    scaling_method <- if (center) "median_scaled_and_centered" else "median_scaled"
  } else {
    naming_ds <- NULL
    scaling_method <- if (center) "centered" else "none"
  }

  # Median-centering of continuous covariates
  if (center) {
    numeric_cols <- names(eres_data)[vapply(eres_data, is.numeric, logical(1))]
    cols_to_center <- setdiff(numeric_cols, c(dv_col, time_col))
    centering_values <- vapply(eres_data[cols_to_center], stats::median, numeric(1), na.rm = TRUE)
    for (col in cols_to_center) {
      eres_data[[col]] <- eres_data[[col]] - centering_values[col]
    }
  } else {
    centering_values <- NULL
  }

  # Update var_spec with transform info from scaling/centering steps
  if (!is.null(var_spec_parsed)) {
    if (scale_to_median) {
      scaled_cols <- intersect(var_spec_parsed$col_name, names(scaling_result$scaling_factors))
      sf <- scaling_result$scaling_factors
      var_spec_parsed <- var_spec_parsed %>%
        dplyr::mutate(
          median_point = dplyr::if_else(
            col_name %in% scaled_cols, sf[col_name], NA_real_),
          transform = dplyr::if_else(
            col_name %in% scaled_cols,
            if (center) "median_scaled_and_centered" else "median_scaled",
            transform)
        )
    } else if (center && !is.null(centering_values)) {
      centered_cols <- intersect(var_spec_parsed$col_name, names(centering_values))
      cv <- centering_values
      var_spec_parsed <- var_spec_parsed %>%
        dplyr::mutate(
          median_point = dplyr::if_else(
            col_name %in% centered_cols, cv[col_name], NA_real_),
          transform = dplyr::if_else(
            col_name %in% centered_cols, "centered", transform)
        )
    }
  }

  if(!is.null(dvid_col)) {
    modds <- dplyr::rename(eres_data, `:=`(DVID_name, !!dvid_col))
    modds <- dplyr::filter(modds, DVID_name == dvid)
  } else {
    modds <- eres_data
  }
  
  if(!is.null(dv_col)) {
    modds <- dplyr::rename(modds, `:=`(DV, !!dv_col))
  }
  if(!is.null(time_col)) {
    modds <- dplyr::rename(modds, `:=`(TIME, !!time_col))
  }
  ## Validate nonlinear spec against the processed dataset
  if (use_nonlinear) validate_nonlinear_spec(nl_spec, modds)

  ## ------------------------------
  ## Initialise variables
  ## ------------------------------
  Relation <- c()
  transform_registry <- list()

  # Exclude adjust_for from forward candidate covariates
  covariates_scm <- if (!is.null(adjust_for)) setdiff(covariates, adjust_for) else covariates

  # start_covariates: pre-loaded into Relation before forward selection begins.
  # They are NOT protected -- backward can remove them.
  # Forward addition candidates exclude them (they are already in the model).
  if (!is.null(start_covariates)) {
    covariates_scm_forward <- setdiff(covariates_scm, start_covariates)
  } else {
    covariates_scm_forward <- covariates_scm
  }

  ## ------------------------------
  ## SCM first step forward
  ## ------------------------------
  if (verbose) message("--- SCM Forward Selection: Step 1 (exposure screening) ---")

  if (use_nonlinear) {
    models_step1 <- scm_forward_first_step_nonlinear(
      eres_data          = modds,
      exposure_vars      = exposure_vars,
      nonlinear_spec     = nl_spec,
      model_type         = "coxph",
      family             = "binomial",   # ignored for coxph
      adjust_for         = adjust_for,
      significance_level = significance_level[1],
      linearity_alpha    = linearity_alpha,
      fp_powers          = fp_powers,
      spline_df          = spline_df,
      verbose            = verbose
    )
    checks_step1 <- list(
      toto_tot     = models_step1$toto_tot,
      Relation     = models_step1$Relation,
      finalforward = models_step1$finalforward,
      flagforst    = models_step1$flagforst
    )
    transform_registry <- models_step1$transform_registry
  } else {
    models_step1 <- scm_forward_first_step_coxph(
      eres_data          = modds,
      exposure_vars      = exposure_vars,
      adjust_for         = adjust_for,
      dvid               = dvid, dvid_col = dvid_col,
      dv_col             = dv_col, time_col = time_col,
      significance_level = significance_level[1],
      verbose            = verbose,
      return_models      = return_models
    )
    checks_step1 <- .check_forward(run_summary = models_step1$step1_sum, relation = Relation)
  }

  # For the LINEAR path only: prepend adjust_for (nonlinear path already includes it).
  if (!use_nonlinear && !is.null(adjust_for)) {
    checks_step1$Relation <- c(adjust_for, checks_step1$Relation)
  }

  if (!checks_step1$flagforst) {
    if (verbose) {
      message("No significant exposure variable found. SCM complete.")
      message("End of the SCM.")
    }
    # Return a fully-formed S3 object so all downstream methods work correctly.
    early_metadata <- list(
      dvid_name               = dvid,
      time_col                = time_col,
      dv_col                  = dv_col,
      candidate_covariates    = c(exposure_vars, covariates),
      candidate_exposure_vars = exposure_vars,
      selected_exposure       = NA_character_,
      finalRelation           = character(0L),
      selected_terms          = character(0L),
      adjust_for              = adjust_for,
      start_covariates        = start_covariates,
      direction               = direction,
      finalrun                = FALSE,
      centering_values        = centering_values,
      scaling_method          = scaling_method,
      transform_registry      = if (use_nonlinear) list() else NULL,
      nonlinear_spec          = nl_spec,
      linearity_alpha         = if (use_nonlinear) linearity_alpha else NULL,
      fp_powers               = if (use_nonlinear) fp_powers else NULL,
      spline_df               = if (use_nonlinear) spline_df else NULL,
      has_nonlinear_terms     = FALSE
    )
    out <- list(
      final_model = NULL,
      path        = checks_step1$toto_tot,
      summary     = tibble::tibble(),
      model_data  = modds,
      metadata    = early_metadata,
      var_spec    = var_spec_parsed,
      naming_df   = naming_ds
    )
    class(out) <- c("eres_scm_coxph", "eres_scm", "list")
    return(out)
  } else {
    exposure_var <- utils::tail(checks_step1$Relation, 1)
    if (verbose) message("  Exposure variable selected: ", exposure_var)
  }

  forward_count = 2
  flagforst     <- checks_step1$flagforst
  finalforward  <- checks_step1$finalforward
  Relation      <- checks_step1$Relation

  # Inject start_covariates into Relation before forward selection begins.
  # They are pre-loaded but not re-tested in forward addition and not protected
  # from backward elimination.
  if (!is.null(start_covariates) && length(start_covariates) > 0L) {
    Relation <- c(Relation, start_covariates)
    checks_step1$Relation <- Relation
    if (verbose) message("  Pre-loaded start_covariates: ", paste(start_covariates, collapse = ", "))
  }

  current_model <- checks_step1
  current_model$Relation <- Relation
  if (use_nonlinear) current_model$transform_registry <- transform_registry

  # Forward addition loop -- skipped when direction == "backward"
  if (direction != "backward") {
    while (flagforst == TRUE && length(Relation) <= length(covariates_scm)) {
      if (verbose) message("--- SCM Forward Selection: Step ", forward_count, " ---")
      if (use_nonlinear) {
        tmp_list <- scm_forward_add_nonlinear(
          eres_data          = modds,
          current_model      = current_model,
          covariates         = covariates_scm_forward,
          nonlinear_spec     = nl_spec,
          model_type         = "coxph",
          family             = "binomial",   # ignored for coxph
          significance_level = significance_level[1],
          linearity_alpha    = linearity_alpha,
          fp_powers          = fp_powers,
          spline_df          = spline_df,
          verbose            = verbose
        )
        current_model <- list(
          toto_tot           = tmp_list$toto_tot,
          Relation           = tmp_list$Relation,
          finalforward       = tmp_list$finalforward,
          flagforst          = tmp_list$flagforst,
          transform_registry = tmp_list$transform_registry
        )
        transform_registry <- tmp_list$transform_registry
      } else {
        tmp_list <- scm_forward_add_coxph(
          eres_data          = modds,
          current_model      = current_model,
          covariates         = covariates_scm_forward,
          dvid_col           = dvid_col,
          dvid               = dvid,
          time_col           = time_col,
          significance_level = significance_level[1],
          verbose            = verbose
        )
        current_model <- list(
          toto_tot     = tmp_list$toto_tot,
          Relation     = tmp_list$Relation,
          finalforward = tmp_list$finalforward,
          flagforst    = tmp_list$flagforst
        )
      }
      flagforst <- current_model$flagforst
      Relation  <- current_model$Relation
      if (verbose && flagforst) message("  Covariate added: ", utils::tail(Relation, 1))
      if (verbose && !flagforst) message("  No additional significant covariate found.")
      forward_count = forward_count + 1
    }
    if (verbose) message("Forward selection complete -- ", length(Relation), " covariate(s) selected.")
  } else {
    if (verbose) message("Skipping forward addition (direction = 'backward').")
  }

  ## ------------------------------

  ## SCM Forward Interaction Step (if requested)
  ## Skipped when direction == "backward"
  ## ------------------------------
  if (direction != "backward" && length(candidate_interactions) > 0 && interaction_strategy == "after_main") {
    # Get eligible interactions (those whose parents are in current Relation)
    eligible <- get_eligible_interactions(candidate_interactions, current_model$Relation)
    # Remove interactions already in model
    eligible <- eligible[!eligible %in% current_model$Relation]

    if (length(eligible) > 0 && verbose) {
      message("--- SCM Forward Selection: Interaction Terms ---")
      message("  Eligible interactions: ", paste(eligible, collapse = ", "))
    }

    flag_int <- length(eligible) > 0

    while (flag_int) {
      # Re-check eligibility (Relation may have changed)
      eligible <- get_eligible_interactions(candidate_interactions, current_model$Relation)
      eligible <- eligible[!eligible %in% current_model$Relation]

      if (length(eligible) == 0) {
        flag_int <- FALSE
        break
      }

      tmp_list <- scm_forward_add_interaction_coxph(
        eres_data = modds,
        current_model = current_model,
        candidate_interactions = eligible,
        significance_level = p_incl_int,
        verbose = verbose,
        return_models = return_models
      )

      current_model <- list(
        toto_tot = tmp_list$toto_tot,
        Relation = tmp_list$Relation,
        finalforward = tmp_list$finalforward,
        flagforst = tmp_list$flagforst
      )

      flag_int <- tmp_list$flagforst
      Relation <- current_model$Relation
      if (verbose && flag_int) message("  Interaction added: ", utils::tail(current_model$Relation, 1))
      if (verbose && !flag_int) message("  No additional significant interaction found.")
    }
    if (verbose) message("Interaction forward selection complete.")
  }

  if (direction != "backward" && length(candidate_interactions) > 0 && interaction_strategy == "simultaneous") {
    # For simultaneous strategy, re-run forward loop including interactions
    flag_sim <- TRUE
    while (flag_sim) {
      # Get eligible interactions
      eligible <- get_eligible_interactions(candidate_interactions, current_model$Relation)
      eligible <- eligible[!eligible %in% current_model$Relation]

      if (length(eligible) > 0) {
        # Test interactions
        tmp_int <- scm_forward_add_interaction_coxph(
          eres_data = modds,
          current_model = current_model,
          candidate_interactions = eligible,
          significance_level = p_incl_int,
          verbose = verbose,
          return_models = return_models
        )
        if (tmp_int$flagforst) {
          current_model <- list(
            toto_tot = tmp_int$toto_tot,
            Relation = tmp_int$Relation,
            finalforward = tmp_int$finalforward,
            flagforst = tmp_int$flagforst
          )
          Relation <- current_model$Relation
          if (verbose) message("  Interaction added (simultaneous): ", utils::tail(current_model$Relation, 1))
          next
        }
      }
      flag_sim <- FALSE
    }
  }

  finalrun<-finalforward
  backward_count = 1
  if (direction != "forward" && length(current_model$Relation) > 1){
    if (verbose) message("--- SCM Backward Elimination ---")
    flagforst = TRUE # flag for loop
    while(flagforst == TRUE && length(current_model$Relation) > 1){
      adj_backward_count = backward_count + forward_count - 1; # used for continuity for appending to lists with forward addition data
      # call backward elimination function until no significant covariate is identified
      if (use_nonlinear) {
        tmp_list <- scm_backward_delete_nonlinear(
          eres_data              = modds,
          current_model          = current_model,
          model_type             = "coxph",
          family                 = "binomial",   # ignored for coxph
          exit_alpha             = significance_level[2],
          exit_alpha_interaction = p_excl_int,
          adjust_for             = adjust_for,
          transform_registry     = transform_registry,
          simplify               = TRUE,
          linearity_alpha        = linearity_alpha,
          verbose                = verbose
        )
        transform_registry <- tmp_list$transform_registry
        current_model <- list(
          toto_tot           = tmp_list$toto_tot,
          Relation           = tmp_list$Relation,
          finalforward       = tmp_list$finalforward,
          flagforst          = tmp_list$flagforst,
          transform_registry = transform_registry
        )
      } else {
        tmp_list = scm_backward_delete_coxph(eres_data = modds, current_model = current_model,
                                       exit_alpha = significance_level[2],
                                       exit_alpha_interaction = p_excl_int,
                                       adjust_for = adjust_for,
                                       verbose = verbose)
        current_model <- list(
          toto_tot = tmp_list$toto_tot,
          Relation = tmp_list$Relation,
          finalforward = tmp_list$finalforward,
          flagforst = tmp_list$flagforst
        )
      }
      Relation<-current_model$Relation
      flagforst<-current_model$flagforst
      backward_count = backward_count + 1
    }
  } # only perform backward elimination if there is multiple significant covariates identified in forward addition
  ## ------------------------------
  ## check if the only covariate is significant against null
  ## ------------------------------
  models_step_last <- current_model$toto_tot
  # for saving
  models_step_last <- models_step_last %>%
    dplyr::select(-SELECT,-Model) %>%
    dplyr::rename(Significant=Signif,Model=RUNTEST)%>%
    dplyr::mutate(Significant=ifelse(Significant==0,'No','Yes')) %>%
    dplyr::mutate(Model=gsub('reg', '', Model))%>%
    dplyr::mutate(pvalue=ifelse(pvalue>0.001,
                                signif_pad_like(pvalue,digits=4),'<0.001' ) )
  current_model$models_step_last <- models_step_last
  final_model <- run_sum_tbl_coxph(data = modds, covariates = current_model$Relation )
  path <- current_model$models_step_last
  summary <- get_sum_tbl_coxph(final_model, exponentiate = TRUE, include_conf_int = TRUE)
  if (use_nonlinear) {
    summary <- .enhance_summary_nonlinear(summary, transform_registry)
  }
  selected_terms <- if (!is.null(adjust_for)) setdiff(current_model$Relation, adjust_for) else current_model$Relation
  has_nonlinear_terms <- if (use_nonlinear && length(transform_registry) > 0L) {
    any(vapply(transform_registry, function(x) x$transform_type != "linear", logical(1)))
  } else {
    FALSE
  }
  metadata <- list(
    dvid_name               = dvid,
    time_col                = time_col,
    dv_col                  = dv_col,
    candidate_covariates    = c(exposure_vars, covariates),
    candidate_exposure_vars = exposure_vars,
    selected_exposure       = exposure_var,
    finalRelation           = current_model$Relation,
    selected_terms          = selected_terms,
    adjust_for              = adjust_for,
    start_covariates        = start_covariates,
    direction               = direction,
    finalrun                = current_model$finalforward,
    centering_values        = centering_values,
    scaling_method          = scaling_method,
    transform_registry      = if (use_nonlinear) transform_registry else NULL,
    nonlinear_spec          = nl_spec,
    linearity_alpha         = if (use_nonlinear) linearity_alpha else NULL,
    fp_powers               = if (use_nonlinear) fp_powers else NULL,
    spline_df               = if (use_nonlinear) spline_df else NULL,
    has_nonlinear_terms     = has_nonlinear_terms
  )
  out <- list(final_model = final_model, path = path, summary = summary, model_data = modds, metadata = metadata, var_spec = var_spec_parsed, naming_df = naming_ds )
  class(out) <- c("eres_scm_coxph", "eres_scm", "list")

  if (verbose) {
    message("--- SCM Complete ---")
    message("  Final covariates: ", paste(current_model$Relation, collapse = ", "))
  }

  return(out)
}
#' Perform a backward deletion step in SCM
#'
#' Executes a single backward-deletion iteration in the Stepwise Covariate Modeling
#' workflow: evaluates removable terms currently in the model, tests each for removal
#' using the chosen criterion, deletes the least-contributing qualifying term, and
#' returns updated artifacts.
#'
#' @param current_model model or formula/spec. The current fitted model (or a model
#'   specification) before backward deletion. For formulas, the function will fit as needed.
#' @param eres_data data.frame or tibble. Dataset used for fitting, containing the outcome,
#'   exposure metrics, and covariates.
#' @param removable_terms character. Vector of term names that are eligible for deletion
#' @param outcome character. Not used. Column name of the outcome variable.
#'   (typically non-fixed covariates currently present in \code{current_model}).
#' @param fixed_terms character. Optional terms that must remain in the model and should
#'   not be considered for deletion (e.g., exposure or baseline risk factors). Default \code{NULL}.
#'   Default \code{"binomial"}.
#' @param criterion character. Selection criterion guiding deletion decisions. One of
#'   \code{"AIC"}, \code{"BIC"}, or \code{"LRT"} (likelihood ratio test vs. full model).
#'   Default \code{"AIC"}.
#' @param exit_alpha numeric. Significance threshold for removal when \code{criterion = "LRT"}.
#'   Ignored otherwise. Default \code{0.01}.
#' @param exit_alpha_interaction numeric or NULL. Separate significance threshold for
#'   removing interaction terms during backward elimination. If NULL (default), interaction
#'   terms are evaluated against \code{exit_alpha} (same as main effects). When provided,
#'   interaction terms (those containing ":") are compared against this threshold while
#'   main effects continue to use \code{exit_alpha}.
#' @param min_terms integer. Optional lower bound on the number of terms allowed to remain
#'   in the model (excluding intercept). If reached, no deletion occurs. Default \code{0}.
#' @param run_tag character. Optional label stored in metadata. Default \code{""}.
#' @param return_models logical. If \code{TRUE}, include fitted candidate removal models in the output.
#'   Default \code{FALSE}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{finalforward}: row index in \code{toto_tot} of the term removed this
#'         step, or the previous \code{current_model$finalforward} if none qualified.
#'   \item \code{toto_back1}: tibble summarizing this step's removable-term candidates
#'         (same columns as \code{\link{scm_forward_first_step_coxph}}'s \code{step1_sum}).
#'   \item \code{toto_tot}: cumulative history of all step summaries so far
#'         (\code{current_model$toto_tot} with \code{toto_back1} appended).
#'   \item \code{Relation}: character vector of terms remaining in the model after
#'         this step (with the deleted term, if any, removed).
#'   \item \code{flagforst}: logical; \code{TRUE} if a term was deleted this step.
#'   \item \code{models} (optional): named list of fitted candidate-removal \code{coxph}
#'         objects if \code{return_models = TRUE}.
#' }
#' The term deleted in this step (if any) is
#' \code{setdiff(current_model$Relation, Relation)}. This return value has the same
#' shape as \code{current_model}, so it can be fed back into another
#' \code{scm_backward_delete_coxph()} call or into \code{\link{scm_forward_add_coxph}}.
#'
#' @details
#' For each term in \code{removable_terms} that is present in the model and not in
#' \code{fixed_terms}, the function fits a Cox model excluding that term and compares
#' it by likelihood-ratio test against the full model. The term with the largest
#' relative exceedance of \code{exit_alpha} is removed. \code{min_terms} places a
#' floor on how many terms may remain. Data preparation (missingness, factor
#' handling, scaling) should be done prior to calling this function.
#'
#' \code{current_model} must be a list with the same shape as the value returned by
#' this function or by \code{\link{scm_forward_add_coxph}} -- it is not a fitted
#' model object. This low-level building block is normally orchestrated internally
#' by \code{\link{scm_coxph}}; most users should call \code{scm_coxph()} directly
#' rather than chaining these step functions by hand.
#'
#' @examples
#' data(demo_eres)
#' demo_eres$DV <- demo_eres$sim_status
#' demo_eres$TIME <- demo_eres$sim_time
#' \donttest{
#' # Build a small current_model with two terms already selected (see
#' # scm_forward_add_coxph() for how this list is constructed step by step;
#' # this bookkeeping is normally done for you inside scm_coxph()).
#' best1 <- scm_forward_first_step_coxph(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", dvid = NULL,
#'   significance_level = 0.2, run_tag = "step1"
#' )
#' x <- best1$step1_sum
#' current_model <- list(
#'   toto_tot = x,
#'   Relation = x$Model[x$Model == x$SELECT][!is.na(x$Model[x$Model == x$SELECT])],
#'   finalforward = x$Run[x$Model == x$SELECT][1],
#'   flagforst = any(!is.na(x$Run[x$Model == x$SELECT]))
#' )
#' step_add <- scm_forward_add_coxph(
#'   current_model = current_model, eres_data = demo_eres,
#'   covariates = c("AGE", "WT"), dvid = NULL,
#'   significance_level = 0.2, run_tag = "fwd_add_1"
#' )
#'
#' step_b <- scm_backward_delete_coxph(
#'   current_model = step_add,
#'   eres_data = demo_eres,
#'   removable_terms = step_add$Relation,
#'   outcome = "sim_status",
#'   criterion = "BIC",
#'   exit_alpha = 0.01,
#'   min_terms = 1,
#'   run_tag = "bwd_del_1"
#' )
#' step_b$Relation
#' setdiff(step_add$Relation, step_b$Relation)
#' }
#'
#' @seealso \code{\link{scm_forward_add_coxph}}, \code{\link{scm_coxph}}, \code{\link{fit_univariate_models_coxph}}
#' @inheritParams scm_coxph
#' @export
scm_backward_delete_coxph <- function(
    current_model,
    eres_data,
    removable_terms,
    outcome,
    fixed_terms = NULL,
    criterion = "AIC",
    exit_alpha = 0.01,
    exit_alpha_interaction = NULL,
    adjust_for = NULL,
    min_terms = 0,
    run_tag = "",
    return_models = FALSE,
    verbose = FALSE
) {
  # implementation goes here
  modds <- eres_data
  tmp_ref <- current_model$finalforward
  toto_tot <- current_model$toto_tot
  Relation <- unique(current_model$Relation)
  covariates_back1 <- unique(current_model$Relation)

  # Hierarchy protection: exclude main effects that have child interactions in the model
  protected <- get_protected_main_effects(covariates_back1)
  # Also protect adjust_for terms from backward elimination
  protected <- union(protected, if (!is.null(adjust_for)) adjust_for else character(0))
  removable_back1 <- covariates_back1[!covariates_back1 %in% protected]

  tmp_base_cov <- glue::glue_collapse(covariates_back1, sep = " + ")
  covback1 <- NA
  for (i in seq_along(removable_back1)) {
    coco <- covariates_back1[covariates_back1 != removable_back1[i]]
    coco <- glue::glue_collapse(coco, sep = "+")
    covback1 <- c(covback1, coco)
  }
  covback1 <- covback1[!is.na(covback1)]
  mod_nul <- survival::coxph(build_model_formula("survival::Surv(TIME, DV)", tmp_base_cov), data = modds)
  mod_nul_tab <- summary(mod_nul)
  mod_nul_sum <- broom::glance(mod_nul)
  df_ref = mod_nul_tab$logtest[["df"]]
  glance_ref = broom::glance(mod_nul)
  # initialize lists
  mod_obj <- vector(mode = "list", length = length(covback1))
  mod_sum <- vector(mode = "list", length = length(covback1))
  mod_metric <- vector(mode = "list", length = length(covback1))
  summ.list <- vector(mode = "list", length = length(covback1))
  for (i in seq_along(covback1)) {
    mod_obj[[i]] <- survival::coxph(build_model_formula("survival::Surv(TIME, DV)", covback1[i]), data = modds)
    mod_sum[[i]] <- broom::glance(mod_obj[[i]])
    mod_metric[i] <- removable_back1[i]
    summ.list[[i]] <- summary(mod_obj[[i]])
  }
  metrics = unlist(mod_metric)
  nsub.tab = do.call(rbind, mod_sum) %>%
    dplyr::select(logLik, AIC, nobs = n) %>%
    dplyr::mutate(Metric = metrics)
  outtab = NULL
  for (i in seq_along(metrics)) {
    tabi = summ.list[[i]]
    dftot = tabi$logtest[["df"]]
    tabi = summ.list[[i]]
    tabi.conf = as.data.frame(confint(mod_obj[[i]], level = 0.95))
    tabi.coef = data.frame(tabi$coefficients)
    out.tab = dplyr::bind_cols(tabi.coef, tabi.conf)
    dftest = dftot-df_ref
    loglik = mod_sum[[i]]$logLik
    out.tab = out.tab %>%
      dplyr::mutate(Metric = metrics[i],
                    Effect = row.names(out.tab),
                    RSE = 100 *(se.coef./abs(coef)),
                    estRSE = sprintf("%s (%s)", signif(coef, 3), signif(RSE, 3)),
                    CI = sprintf("(%s, %s)", signif(`2.5 %`,3), signif(`97.5 %`, 3)),
                    pvalue = signif(Pr...z.., 3),
                    dftest =dftest, dftot = dftot, logLik =  loglik ) %>%
      dplyr::select(Effect, Metric, CI, dftest, dftot)
    nul = mod_nul
    fit = mod_obj[[i]]
    pval.fitx = round(anova(nul, fit, test = "Chisq")[2, 4], 10)
    out.tab$pvalue = pval.fitx  ###---replace the Wald test p-value with LRT p-value
    outtab = rbind(outtab, out.tab)
  }
  # only take covariates
  outtab <- outtab[outtab$Effect != "(Intercept)", ]
  outtab <- outtab %>%
    dplyr::filter(!is.na(pvalue)) %>%
    dplyr::left_join(nsub.tab, by = "Metric") %>%
    dplyr::mutate(min2LL = -2*logLik) %>%
    dplyr::mutate(diff2LL = -2*glance_ref$logLik -min2LL) %>%
    dplyr::select(Metric, min2LL, diff2LL, AIC,  pvalue, nobs, dftest, dftot)
  outtab<-outtab[!duplicated(outtab[,c('Metric')]),]
  if(verbose) {
    message(paste(utils::capture.output(print(outtab)), collapse = "\n"))
  }
  # Compute per-term effective threshold: interactions use exit_alpha_interaction if provided
  eff_exit_alpha <- ifelse(
    grepl(":", outtab$Metric, fixed = TRUE) & !is.null(exit_alpha_interaction),
    exit_alpha_interaction,
    exit_alpha
  )

  # Determine which terms exceed their respective threshold (candidates for removal)
  exceeds_threshold <- outtab$pvalue > eff_exit_alpha
  # Among those exceeding, select the one with largest relative exceedance (pvalue / threshold)
  if (any(exceeds_threshold, na.rm = TRUE)) {
    rel_exceedance <- ifelse(exceeds_threshold, outtab$pvalue / eff_exit_alpha, 0)
    select_term <- outtab$Metric[which.max(rel_exceedance)]
  } else {
    select_term <- NA_character_
  }

  toto_back1 <- dplyr::bind_rows(
    tibble::tibble(Model='Base',
                   pvalue = NA_real_,
                   Signif = NA_real_,
                   AIC = mod_nul$aic,
                   min2LL = -2*mod_nul_sum$logLik,
                   diff2LL = 0,
                   nobs = mod_nul_sum$n,
                   dftest = 0,
                   dftot = df_ref ),
    outtab %>%
      dplyr::rename(Model = Metric) %>%
      dplyr::mutate(Model = as.character(Model),
                    SELECT = select_term,
                    Signif = ifelse(pvalue > eff_exit_alpha, 1, 0),
                    Reference = tmp_ref)) %>%
    dplyr::mutate(Run = - 1 + (max(toto_tot$Run) + 1):(max(toto_tot$Run) + dplyr::n())) %>%
    dplyr::relocate(Run) %>%
    tidyr::fill(SELECT, .direction = "downup") %>%
    dplyr::mutate(RUNTEST = ifelse(Run != 1, paste0("Run", Reference, " - ",
                                                    Model), NA_character_))
  #
  toto_back1 <- toto_back1[-1, ]
  toto_tot <- dplyr::bind_rows(toto_tot, toto_back1)
  if (any(!is.na(toto_back1$Run[toto_back1$Model == toto_back1$SELECT]))) {
    finalbackward2 <- toto_back1$Run[toto_back1$Model == toto_back1$SELECT]
    finalrun <- finalbackward2
    Relation <- Relation[!Relation %in% toto_back1$Model[toto_back1$Model ==
                                                           toto_back1$SELECT]]  # remove selected model
  } else {
    finalbackward2 <- NA
    finalrun <- current_model$finalforward
    Relation <- current_model$Relation
  }
  flagbackst <- FALSE
  if (any(!is.na(finalbackward2))) {
    flagbackst <- TRUE
  }
  toto_back1 <- toto_back1 %>%
    dplyr::select(Run, min2LL, diff2LL, AIC, Reference, pvalue, Signif, nobs, dftest, dftot, RUNTEST,  SELECT, Model, AIC )
  if(return_models) {
    names(mod_obj) <- paste("-", metrics)
    res <-  list(finalforward = finalrun,
                 toto_back1 = toto_back1,
                 toto_tot = toto_tot,
                 Relation = Relation,
                 flagforst = flagbackst,
                 models = mod_obj)
  }  else {
    res <-  list(finalforward = finalrun,
                 toto_back1 = toto_back1,
                 toto_tot = toto_tot,
                 Relation = Relation,
                 flagforst = flagbackst)
  }
  # combine into list to return
  return(res)
}
#' Perform a forward addition step in SCM
#'
#' Executes a single forward-addition iteration in the Stepwise Covariate Modeling
#' workflow: evaluates candidate terms, selects the best qualifying term based on
#' the chosen criterion, updates the current model, and returns updated artifacts.
#'
#' @param current_model model or formula/spec. The current fitted model (or a model
#'   specification) before forward addition. For formulas, the function will fit as needed.
#' @param eres_data data.frame or tibble. Dataset used for fitting, containing the outcome,
#'   exposure metrics, and candidate covariates.
#'   (e.g., covariates or exposure variables not yet in \code{current_model}).
#' @param adjust_for character. Optional fixed terms to include in all models (e.g., baseline risk factors).
#'   Default \code{NULL}.
#'   Default \code{"binomial"}.
#' @param criterion character. Selection criterion to rank candidates. One of \code{"AIC"},
#'   \code{"BIC"}, or \code{"LRT"}. Default \code{"AIC"}.
#' @param max_terms integer. Optional cap on total terms allowed in the model (excluding intercept).
#'   If reached, the function will not add new terms. Default \code{Inf}.
#' @param run_tag character. Optional label stored in metadata. Default \code{""}.
#' @param return_models logical. If \code{TRUE}, include fitted candidate models in the output. Default \code{FALSE}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{tmp_cov}: character vector of candidate covariates that were tested
#'         this step (sorted, excluding terms already in the model).
#'   \item \code{finalforward}: row index in \code{toto_tot} of the term added this
#'         step, or the previous \code{current_model$finalforward} if none qualified.
#'   \item \code{toto_2}: tibble summarizing this step's candidates (same columns as
#'         \code{\link{scm_forward_first_step_coxph}}'s \code{step1_sum}).
#'   \item \code{toto_tot}: cumulative history of all forward-step summaries so far
#'         (\code{current_model$toto_tot} with \code{toto_2} appended).
#'   \item \code{Relation}: character vector of all terms currently in the model,
#'         updated to include the newly added term (if any).
#'   \item \code{flagforst}: logical; \code{TRUE} if a candidate was added this step.
#'   \item \code{models} (optional): named list of fitted candidate \code{coxph}
#'         objects if \code{return_models = TRUE}.
#' }
#' The term added in this step (if any) is
#' \code{setdiff(Relation, current_model$Relation)}. This return value has the same
#' shape as \code{current_model}, so it can be fed back into \code{scm_forward_add_coxph()}
#' or \code{\link{scm_backward_delete_coxph}} for the next step.
#'
#' @details
#' For each candidate in \code{covariates} not already present in
#' \code{current_model$Relation}, the function fits a Cox model of the form
#' \code{survival::Surv(TIME, DV) ~ adjust_for + current_terms + candidate}, compares
#' it by likelihood-ratio test against the model without the candidate, and adds the
#' best-qualifying candidate (by smallest AIC among significant candidates) if
#' \code{significance_level} is met. Data preparation (missingness, factor handling,
#' scaling) should be done prior to calling this function.
#'
#' \code{current_model} must be a list with the same shape as the value returned by
#' this function (or, for the very first call, reconstructed from
#' \code{\link{scm_forward_first_step_coxph}}'s output as shown in the example) -- it
#' is not a fitted model object. This low-level building block is normally
#' orchestrated internally by \code{\link{scm_coxph}}; most users should call
#' \code{scm_coxph()} directly rather than chaining these step functions by hand.
#'
#' @examples
#' data(demo_eres)
#' demo_eres$DV <- demo_eres$sim_status
#' demo_eres$TIME <- demo_eres$sim_time
#' \donttest{
#' # Build the current_model list expected by scm_forward_add_coxph() from the
#' # output of scm_forward_first_step_coxph() (this bookkeeping is normally
#' # done for you inside scm_coxph()).
#' best1 <- scm_forward_first_step_coxph(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", dvid = NULL,
#'   significance_level = 0.2, run_tag = "step1"
#' )
#' x <- best1$step1_sum
#' current_model <- list(
#'   toto_tot = x,
#'   Relation = x$Model[x$Model == x$SELECT][!is.na(x$Model[x$Model == x$SELECT])],
#'   finalforward = x$Run[x$Model == x$SELECT][1],
#'   flagforst = any(!is.na(x$Run[x$Model == x$SELECT]))
#' )
#'
#' step <- scm_forward_add_coxph(
#'   current_model = current_model,
#'   eres_data = demo_eres,
#'   covariates = c("AGE", "WT"),
#'   dvid = NULL,
#'   criterion = "AIC",
#'   significance_level = 0.2,
#'   run_tag = "fwd_add_1"
#' )
#' step$Relation
#' setdiff(step$Relation, current_model$Relation)
#' }
#'
#' @seealso \code{\link{scm_forward_first_step_coxph}}, \code{\link{scm_coxph}}, \code{\link{fit_univariate_models_coxph}}
#' @inheritParams scm_coxph
#' @export
scm_forward_add_coxph <- function(
    current_model,
    eres_data,
    covariates,
    dvid_col = NULL,
    time_col = NULL,
    dvid = NULL,
    adjust_for = NULL,
    criterion = "AIC",
    significance_level = 0.01,
    max_terms = Inf,
    run_tag = "",
    return_models = FALSE,
    verbose = FALSE
) {
  modds <- eres_data
  ## second step
  tmp_ref <- current_model$finalforward
  Relation <- current_model$Relation
  toto <- current_model$toto_tot
  covariates_2 <- covariates[(covariates %!in% unique(toto$SELECT)) & !(covariates %in% Relation)]  # removing all exposure coviartes  significant
  tmp_cov <- sort(covariates_2)
  tmp_cov_ref <- unique(Relation)
  tmp_base_cov <- paste(tmp_cov_ref, collapse = "+")
  tmp_cov_combo = paste(tmp_base_cov, "+", tmp_cov)
  mod_nul <- survival::coxph(build_model_formula("survival::Surv(TIME, DV)", tmp_base_cov), data = modds)
  mod_nul_sum <- broom::glance(mod_nul)
  mod_nul_tab <- summary(mod_nul)
  df_ref = mod_nul_tab$logtest[["df"]]
  glance_ref = broom::glance(mod_nul)
  # initialize lists
  mod_obj <- vector(mode = "list", length = length(tmp_cov_combo))
  mod_sum <- vector(mode = "list", length = length(tmp_cov_combo))
  mod_metric <- vector(mode = "list", length = length(tmp_cov_combo))
  summ.list <- vector(mode = "list", length = length(tmp_cov_combo))
  for (i in seq_along(tmp_cov)) {
    mod_obj[[i]] <- survival::coxph(build_model_formula("survival::Surv(TIME, DV)", tmp_cov_combo[i]), data = modds)
    mod_sum[[i]] <- broom::glance(mod_obj[[i]])
    mod_metric[i] <- tmp_cov[i]
    summ.list[[i]] <- summary(mod_obj[[i]])
  }
  metrics = unlist(mod_metric)
  nsub.tab = do.call(rbind, mod_sum) %>%
    dplyr::select(logLik, AIC,  nobs = n) %>%
    dplyr::mutate(Metric = metrics)
  if(verbose) {
    message(paste(utils::capture.output(print(nsub.tab)), collapse = "\n"))
  }
  outtab = NULL
  for (i in seq_len(length(metrics)) ) {
    tabi = summ.list[[i]]
    dftot = tabi$logtest[["df"]]
    tabi.conf = as.data.frame(confint(mod_obj[[i]], level = 0.95))
    tabi.coef = data.frame(tabi$coefficients)
    out.tab = dplyr::bind_cols(tabi.coef, tabi.conf)
    dftest = dftot-df_ref
    loglik = mod_sum[[i]]$logLik
    out.tab = out.tab %>%
      dplyr::mutate(Metric = metrics[i],
                    Effect = row.names(out.tab),
                    RSE = 100 *
                      (se.coef./abs(coef)),
                    estRSE = sprintf("%s (%s)", signif(coef,3), signif(RSE, 3)),
                    CI = sprintf("(%s, %s)", signif(`2.5 %`,3), signif(`97.5 %`, 3)),
                    pvalue = signif(Pr...z.., 3),
                    dftest =dftest, dftot = dftot, logLik =  loglik ) %>%
      dplyr::select(Effect, Metric, CI, dftest, dftot)
    nul = mod_nul
    fit = mod_obj[[i]]
    pval.fitx = round( anova(nul, fit, test = "Chisq")[2, 4], 10)
    out.tab$pvalue = pval.fitx  ###---replace the Wald test p-value with LRT p-value
    outtab = rbind(outtab, out.tab)
  }
  # only take covariates
  outtab <- outtab[outtab$Effect != "(Intercept)", ]
  outtab <- outtab %>%
    dplyr::filter(!is.na(pvalue)) %>%
    dplyr::left_join(nsub.tab, by = "Metric") %>%
    dplyr::mutate(min2LL = -2*logLik) %>%
    dplyr::mutate(diff2LL = -2*glance_ref$logLik -min2LL) %>%
    dplyr::select(Metric, min2LL, diff2LL, AIC,  pvalue, nobs, dftest, dftot)
  outtab<-outtab[!duplicated(outtab[,c('Metric')]),]
  if(verbose) {
    message(paste(utils::capture.output(print(outtab)), collapse = "\n"))
  }
  toto_2 <- dplyr::bind_rows(
    tibble::tibble(Model='Base',
                   pvalue = NA_real_,
                   Signif = NA_real_,
                   AIC = mod_nul$aic,
                   min2LL = -2*mod_nul_sum$logLik,
                   diff2LL = 0,
                   nobs = mod_nul_sum$n,
                   dftest = 0,
                   dftot = df_ref ),
    outtab %>%
      dplyr::rename(Model=Metric)%>%
      dplyr::mutate(Model = as.character(Model),
                    SELECT=ifelse(min(pvalue , na.rm=TRUE)<significance_level[1],
                                  Model[AIC==min(AIC , na.rm=TRUE) ],
                                  NA_character_),
                    Signif = ifelse(pvalue<significance_level[1],1,0),
                    Reference = tmp_ref)) %>%
    dplyr::mutate(Run= -1 + (max(toto$Run)+1) : (max(toto$Run)+dplyr::n())) %>%
    dplyr::relocate(Run) %>%
    tidyr::fill(SELECT, .direction="downup") %>%
    dplyr::mutate(RUNTEST=ifelse(Run!=1,paste0("Run",Reference," + ",Model),NA_character_))
  toto_2 <- toto_2[-1, ]
  toto_tot <- dplyr::bind_rows(current_model$toto_tot, toto_2)
  if (!all(is.na(toto_2$Run[toto_2$Model == toto_2$SELECT]))) {
    finalforward <- toto_2$Run[toto_2$Model == toto_2$SELECT]
    Relation <- c(Relation, toto_2$Model[toto_2$Model == toto_2$SELECT])
  } else {
    finalforward <- current_model$finalforward
    Relation <- current_model$Relation
  }
  flagforst <- FALSE
  if (any(!is.na(toto_2$Run[toto_2$Model == toto_2$SELECT]))) {
    flagforst <- TRUE
  }
  toto_2 <- toto_2 %>%
    dplyr::select(Run, min2LL, diff2LL, AIC,
                  Reference, pvalue, Signif,
                  nobs, dftest, dftot, RUNTEST,
                  SELECT, Model)
  if(return_models) {
    names(mod_obj) <- metrics
    res <-  list(tmp_cov = tmp_cov,
                 finalforward = finalforward,
                 toto_2 = toto_2,
                 toto_tot = toto_tot,
                 Relation = Relation,
                 flagforst = flagforst,
                 models = mod_obj)
  }  else {
    res <-  list(tmp_cov = tmp_cov,
                 finalforward = finalforward,
                 toto_2 = toto_2,
                 toto_tot = toto_tot,
                 Relation = Relation,
                 flagforst = flagforst)
  }
  # combine into list to return
  return(res)
}
#' Perform the First Forward Step Using Selected Exposure Variables
#'
#' Executes the initial forward-selection step of the SCM workflow, evaluating only
#' pre-selected exposure variables as candidates for inclusion in the model. The function
#' fits one model per exposure candidate (optionally adjusting for fixed covariates),
#' ranks them by the chosen selection criterion, and returns the best-performing exposure
#' along with fit statistics to decide whether to enter it into the model.
#'
#' @param dvid character. Event value of interest (matches values in \code{dvid_col}).
#' @param eres_data data.frame or tibble. Exposure-Safety dataset containing \code{dvid_col},
#'   outcome, exposure variables, and covariates.
#' @param exposure_vars character. Vector of exposure variable names to be tested in the
#'   first forward step.
#'   For binary event occurrence, this is typically a 0/1 indicator.
#' @param dvid_col character. Column name in \code{eres_data} containing event labels. Default \code{"DVID_name"}.
#' @param adjust_for character. Optional covariates (e.g., baseline risk factors) to include in
#'   all candidate models. Default \code{NULL}.
#'   Default \code{"binomial"}.
#' @param criterion character. Model selection criterion used to rank candidates. One of
#'   \code{"AIC"}, \code{"BIC"}, or \code{"LRT"} (likelihood ratio test vs. null/adjusted-only).
#'   Default \code{"AIC"}.
#'   Ignored otherwise. Default \code{0.01}.
#' @param run_tag character. Optional label for metadata and saved artifacts. Default \code{""}.
#'
#' @return A list with one element:
#' \itemize{
#'   \item \code{step1_sum} (or \code{sum_models}/\code{models} additionally when
#'     \code{return_models = TRUE}): a tibble with one row per exposure candidate
#'     (plus a \code{"Base"} reference row), containing \code{Run}, \code{min2LL},
#'     \code{diff2LL}, \code{AIC}, \code{Reference}, \code{pvalue} (LRT vs. the
#'     base model), \code{Signif} (0/1 flag for \code{pvalue < significance_level}),
#'     \code{nobs}, \code{dftest}, \code{dftot}, \code{RUNTEST}, \code{SELECT}
#'     (name of the best-performing exposure), and \code{Model} (the exposure
#'     tested on that row).
#' }
#' This raw output is post-processed internally by \code{\link{scm_coxph}} into the
#' \code{current_model} structure consumed by \code{\link{scm_forward_add_coxph}}.
#'
#' @details
#' The function subsets \code{eres_data} to rows matching \code{dvid} in \code{dvid_col} (if applicable),
#' then fits a set of candidate Cox models of the form:
#' \code{survival::Surv(TIME, DV) ~ adjust_for + exposure_var}.
#' Candidates are ranked by likelihood-ratio test p-value against the baseline model
#' \code{survival::Surv(TIME, DV) ~ adjust_for}, and \code{significance_level} decides
#' entry. Missing data handling and factor/scale preparation should be performed prior
#' to calling this function.
#'
#' @examples
#' data(demo_eres)
#' # NOTE: unlike scm_coxph(), this low-level step does not rename dv_col/time_col
#' # for you -- the columns must already be literally called "DV" and "TIME".
#' demo_eres$DV <- demo_eres$sim_status
#' demo_eres$TIME <- demo_eres$sim_time
#' \donttest{
#' best1 <- scm_forward_first_step_coxph(
#'   eres_data = demo_eres,
#'   exposure_vars = c("CMAXss", "AUCss"),
#'   dvid = NULL,
#'   criterion = "AIC",
#'   significance_level = 0.05,
#'   run_tag = "ERES_step1"
#' )
#' best1$step1_sum
#' }
#'
#' @seealso \code{\link{scm_coxph}}, \code{\link{fit_univariate_models_coxph}}
#' @inheritParams scm_coxph
#' @export
scm_forward_first_step_coxph <- function(
    eres_data,
    exposure_vars,
    dvid,
    dvid_col = "DVID_name",
    dv_col = "DV",
    time_col = "TIME",
    adjust_for = NULL,
    criterion = "AIC",
    significance_level = 0.01,
    return_models =FALSE,
    verbose = FALSE,
    run_tag = ""
) {
  modds <- eres_data
  ## ------------------------------
  ## Fit reference model
  ## ------------------------------
  base_rhs <- if (!is.null(adjust_for) && length(adjust_for) > 0) {
    paste(adjust_for, collapse = " + ")
  } else {
    "1"
  }
  moderef <- survival::coxph(build_model_formula("survival::Surv(TIME, DV)", base_rhs), data = modds)
  df_ref <- 0
  Relation <- c()
  tmp_ref <- 1
  sum_moderef <- summary(moderef)
  glance_ref = broom::glance(moderef)
  # initialize lists
  mod_obj <- vector(mode = "list", length = length(exposure_vars))
  mod_sum <- vector(mode = "list", length = length(exposure_vars))
  mod_metric <- vector(mode = "list", length = length(exposure_vars))
  summ.list <- vector(mode = "list", length = length(exposure_vars))
  for (i in seq_along(exposure_vars)) {
    cand_rhs <- paste(c(adjust_for, exposure_vars[i]), collapse = "+")
    mod_obj[[i]] <- survival::coxph(build_model_formula("survival::Surv(TIME, DV)", cand_rhs), data = modds)
    mod_sum[[i]] <- broom::glance(mod_obj[[i]])
    mod_metric[i] <- exposure_vars[i]
    summ.list[[i]] <- summary(mod_obj[[i]])
  }
  mod_nul <- moderef
  mod_nul_sum <- broom::glance(mod_nul)
  metrics = unlist(mod_metric)
  nsub.tab = do.call(rbind, mod_sum) %>%
    dplyr::select(logLik, AIC,  nobs = n) %>%
    dplyr::mutate(Metric = metrics)
  if(verbose) {
    message(paste(utils::capture.output(print(nsub.tab)), collapse = "\n"))
  }
  outtab = NULL
  for (i in seq_along(metrics)) {
    tabi = summ.list[[i]]
    dftot = tabi$logtest[["df"]]
    tabi.conf = as.data.frame(confint(mod_obj[[i]], level = 0.95))
    tabi.coef = data.frame(tabi$coefficients)
    out.tab = dplyr::bind_cols(tabi.coef, tabi.conf)
    dftest = dftot-df_ref
    out.tab = out.tab %>%
      dplyr::mutate(Metric = metrics[i],
                    Effect = row.names(out.tab),
                    RSE = 100 *
                      (se.coef./abs(coef)),
                    estRSE = sprintf("%s (%s)", signif(coef,3), signif(RSE, 3)),
                    CI = sprintf("(%s, %s)", signif(`2.5 %`,3), signif(`97.5 %`, 3)),
                    pvalue = signif(Pr...z.., 3),
                    dftest =dftest, dftot = dftot ) %>%
      dplyr::select(Effect, Metric, CI, dftest, dftot)
    nul = mod_nul
    fit = mod_obj[[i]]
    pval.fitx = round(anova(nul, fit, test = "Chisq")[2, 4], 10)
    out.tab$pvalue = pval.fitx  ###---replace the Wald test p-value with LRT p-value
    outtab = rbind(outtab, out.tab)
  }
  # only take covariates
  outtab <- outtab[outtab$Effect != "(Intercept)", ]
  outtab <- outtab %>%
    dplyr::left_join(nsub.tab, by = "Metric") %>%
    dplyr::mutate(min2LL = -2*logLik) %>%
    dplyr::mutate(diff2LL = -2*glance_ref$logLik -min2LL) %>%
    dplyr::select(Metric, min2LL, diff2LL, AIC,  pvalue, nobs, dftest, dftot)
  outtab<-outtab[!duplicated(outtab[,c('Metric')]),]
  if(verbose) {
    message(paste(utils::capture.output(print(outtab)), collapse = "\n"))
  }
  toto <- dplyr::bind_rows(tibble::tibble(Model = "Base",
                                          pvalue = NA_real_,
                                          Signif = NA_real_,
                                          AIC = mod_nul$aic,
                                          min2LL = -2*mod_nul_sum$logLik,
                                          diff2LL = 0,
                                          nobs = mod_nul_sum$n,
                                          dftest = 0,
                                          dftot = df_ref),
                           outtab %>%
                             dplyr::rename(Model=Metric) %>%
                             dplyr::mutate(Model = as.character(Model),
                                           SELECT = ifelse(min(pvalue,na.rm = TRUE) < significance_level[1],
                                                           Model[pvalue == min(pvalue, na.rm = TRUE)],
                                                           NA_character_),
                                           Signif = ifelse(pvalue < significance_level[1], 1, 0),
                                           Reference = tmp_ref), ) %>%
    dplyr::mutate(Run = 1:dplyr::n()) %>%
    dplyr::relocate(Run) %>%
    tidyr::fill(SELECT, .direction = "downup") %>%
    dplyr::mutate(RUNTEST = ifelse(Run != 1, paste0("Run", Reference, " + ",
                                                    Model), NA_character_))
  res <- toto %>% dplyr::select(Run,
                               min2LL,
                               diff2LL,
                               AIC,
                               Reference,
                               pvalue,
                               Signif,
                               nobs,
                               dftest,
                               dftot,
                               RUNTEST,
                               SELECT,
                               Model )
  if(return_models) {
    names(mod_obj) <- exposure_vars
    res <-  list(step1_sum = res, models = mod_obj)
  }  else {
    res <-  list(step1_sum = res)
  }
  return(res)
}
#' Fit Univariate Models for Candidate Covariates
#'
#' Fits separate univariate models for each candidate covariate to screen
#' associations with a specified adverse event in an Exposure-Safety dataset.
#' Returns estimates, test statistics, and p-values to inform subsequent SCM steps.
#'
#' @param dvid character. Event value of interest (matches values in \code{dvid_col}).
#' @param eres_data data.frame or tibble. Exposure-Safety dataset containing \code{dvid_col},
#'   exposure metrics, and candidate covariates.
#' @param covariates character. Vector of candidate covariate column names to test.
#' @param significance_level numeric. Screening threshold in (0, 1). Default \code{0.01}.
#' @param dvid_col character. Column name in \code{eres_data} containing event labels. Default \code{"DVID_name"}.
#' @param dv_col character. Column name of the dependent variable. Default \code{"DV"}.
#' @param time_col character. Column name of the time-to-event variable. Default \code{"TIME"}.
#' @param adjust_for character. Optional covariates to include in all models (fixed effects). Default \code{NULL}.
#' @param verbose logical. If TRUE, print progress. Default FALSE.
#' @param return_models logical. If \code{TRUE}, include fitted model objects in output. Default \code{FALSE}.
#' @param parallel logical or integer. Whether to run fits in parallel or number of workers. Default \code{FALSE}.
#'
#' @return A list with one element:
#' \itemize{
#'   \item \code{step1_sum} (or \code{sum_models} when \code{return_models = TRUE}):
#'     a tibble with one row per candidate covariate (plus a \code{"Base"}
#'     reference row), containing \code{Run}, \code{min2LL}, \code{diff2LL},
#'     \code{Reference}, \code{pvalue} (LRT vs. the base model), \code{Signif}
#'     (0/1 flag for \code{pvalue < significance_level}), \code{nobs},
#'     \code{dftest}, \code{dftot}, \code{RUNTEST}, \code{SELECT} (name of the
#'     best-performing covariate), and \code{Model} (the covariate tested on
#'     that row).
#' }
#' When \code{return_models = TRUE}, a second element \code{models} is also
#' returned: a named list of the fitted \code{coxph} objects, one per covariate.
#'
#' @examples
#' data(demo_eres)
#' \donttest{
#' uni <- fit_univariate_models_coxph(
#'   eres_data = demo_eres,
#'   covariates = c("AGE", "WT"),
#'   dv_col = "sim_status",
#'   time_col = "sim_time",
#'   significance_level = c(0.05, 0.01)
#' )
#' uni$step1_sum
#' }
#'
#' @seealso \code{\link[survival]{coxph}}, \code{\link{scm_coxph}}
#' @importFrom broom glance
#' @export
fit_univariate_models_coxph <- function(
    eres_data,
    covariates,
    significance_level = c(0.05, 0.01),
    dvid = NULL,
    dvid_col = NULL,
    dv_col = "DV",
    time_col = "TIME",
    adjust_for = NULL,
    verbose = FALSE,
    return_models = FALSE,
    parallel = FALSE) {
  .check_inputs(modds = eres_data, dvid_col = dvid_col, dvid = dvid, dv_col = dv_col, covariates = covariates)
  if(!is.null(dvid_col)) {
    modds <- dplyr::rename(eres_data, `:=`(DVID_name, !!dvid_col))
    modds <- dplyr::filter(modds, DVID_name == dvid)
  } else {
    modds <- eres_data
  }
  if(!is.null(dv_col)) {
    modds <- dplyr::rename(modds, `:=`(DV, !!dv_col))
  }
  if(!is.null(time_col)) {
    modds <- dplyr::rename(modds, `:=`(TIME, !!time_col))
  }
  ## ------------------------------
  ## Fit reference model
  ## ------------------------------
  moderef <- survival::coxph(survival::Surv(TIME, DV) ~ 1, data = modds)
  df_ref <- 0
  Relation <- c()
  tmp_ref <- 1
  sum_moderef <- summary(moderef)
  glance_ref = broom::glance(moderef)
  # initialize lists
  mod_obj <- vector(mode = "list", length = length(covariates))
  mod_sum <- vector(mode = "list", length = length(covariates))
  mod_metric <- vector(mode = "list", length = length(covariates))
  summ.list <- vector(mode = "list", length = length(covariates))
  for (i in seq_along(covariates)) {
    mod_obj[[i]] <- survival::coxph(build_model_formula("survival::Surv(TIME, DV)", covariates[i]), data = modds)
    mod_sum[[i]] <- broom::glance(mod_obj[[i]])
    mod_metric[i] <- covariates[i]
    summ.list[[i]] <- summary(mod_obj[[i]])
  }
  mod_nul <- moderef
  mod_nul_sum <- broom::glance(mod_nul)
  metrics = unlist(mod_metric)
  nsub.tab = do.call(rbind, mod_sum) %>%
    dplyr::select(logLik, AIC,  nobs = n) %>%
    dplyr::mutate(Metric = metrics)
  if(verbose) {
    message(paste(utils::capture.output(print(nsub.tab)), collapse = "\n"))
  }
  outtab = NULL
  for (i in seq_along(metrics)) {
    tabi = summ.list[[i]]
    dftot = tabi$logtest[["df"]]
    tabi.conf = as.data.frame(confint(mod_obj[[i]], level = 0.95))
    tabi.coef = data.frame(tabi$coefficients)
    out.tab = dplyr::bind_cols(tabi.coef, tabi.conf)
    dftest = dftot-df_ref
    out.tab = out.tab %>%
      dplyr::mutate(Metric = metrics[i],
                    Effect = row.names(out.tab),
                    RSE = 100 *(se.coef./abs(coef)),
                    estRSE = sprintf("%s (%s)", signif(coef, 3), signif(RSE, 3)),
                    CI = sprintf("(%s, %s)", signif(`2.5 %`,3), signif(`97.5 %`, 3)),
                    pvalue = signif(Pr...z.., 3),
                    dftest =dftest, dftot = dftot ) %>%
      dplyr::select(Effect, Metric, CI, dftest, dftot)
    nul = mod_nul
    fit = mod_obj[[i]]
    pval.fitx = round(anova(nul, fit, test = "Chisq")[2, 4], 10)
    out.tab$pvalue = pval.fitx  ###---replace the Wald test p-value with LRT p-value
    outtab = rbind(outtab, out.tab)
  }
  # only take covariates
  outtab <- outtab[outtab$Effect != "(Intercept)", ]
  outtab <- outtab %>%
    dplyr::left_join(nsub.tab, by = "Metric") %>%
    dplyr::mutate(min2LL = -2*logLik) %>%
    dplyr::mutate(diff2LL = -2*glance_ref$logLik -min2LL) %>%
    dplyr::select(Metric, min2LL, diff2LL, AIC,  pvalue, nobs, dftest, dftot)
  outtab<-outtab[!duplicated(outtab[,c('Metric')]),]
  if(verbose) {
    message(paste(utils::capture.output(print(outtab)), collapse = "\n"))
  }
  toto <- dplyr::bind_rows(tibble::tibble(Model = "Base",
                                          pvalue = NA_real_,
                                          Signif = NA_real_,
                                          AIC = mod_nul$aic,
                                          min2LL = -2*mod_nul_sum$logLik,
                                          diff2LL = 0,
                                          nobs = mod_nul_sum$n,
                                          dftest = 0,
                                          dftot = df_ref),
                           outtab %>%
                             dplyr::rename(Model=Metric) %>%
                             dplyr::mutate(Model = as.character(Model),
                                           SELECT = ifelse(min(pvalue,na.rm = TRUE) < significance_level[1],
                                                           Model[pvalue == min(pvalue, na.rm = TRUE)],
                                                           NA_character_),
                                           Signif = ifelse(pvalue < significance_level[1], 1, 0),
                                           Reference = tmp_ref), ) %>%
    dplyr::mutate(Run = 1:dplyr::n()) %>%
    dplyr::relocate(Run) %>%
    tidyr::fill(SELECT, .direction = "downup") %>%
    dplyr::mutate(RUNTEST = ifelse(Run != 1, paste0("Run", Reference, " + ",
                                                    Model), NA_character_))
  res <- toto %>% dplyr::select(Run,
                               min2LL,
                               diff2LL,
                               Reference,
                               pvalue,
                               Signif,
                               nobs,
                               dftest,
                               dftot,
                               RUNTEST,
                               SELECT,
                               Model )
  if(return_models) {
    names(mod_obj) <- covariates
    res <-  list(sum_models = res, models = mod_obj)
  }  else {
    res <-  list(step1_sum = res)
  }
  return(res)
}
