#' Run Stepwise Covariate Modeling (Unified Interface)
#'
#' A unified dispatcher that runs the full SCM procedure using either a GLM
#' or Cox proportional hazards framework. This function delegates to
#' \code{\link{scm_glm}} or \code{\link{scm_coxph}} based on the \code{type}
#' argument, providing a single entry point for programmatic use (e.g.,
#' bootstrap stability assessment).
#'
#' @param eres_data data.frame. The analysis dataset containing events/survival data,
#'   exposure metrics, and candidate covariates.
#' @param exposure_vars character. Vector of candidate exposure variable names.
#'   The best-performing exposure is selected in the first forward step and is
#'   protected throughout the procedure (never removed).
#' @param covariates character. Vector of candidate covariate names to screen.
#' @param type Character, one of \code{"glm"} or \code{"coxph"}. Determines which
#'   modelling framework to use.
#' @param significance_level Numeric vector of length 2. The first value is the
#'   entry threshold for forward selection; the second is the retention threshold
#'   for backward elimination. Default \code{c(0.05, 0.001)}.
#' @param dvid Optional character. Adverse event identifier value to subset on.
#' @param dvid_col Optional character. Column name containing the AE identifier.
#' @param dv_col Character. Name of the dependent variable column. Default \code{"DV"}.
#' @param time_col Character. Name of the time-to-event column (used only when
#'   \code{type = "coxph"}). Default \code{"TIME"}.
#' @param family Character or family object. Model family for GLM (used only when
#'   \code{type = "glm"}). Default \code{"binomial"}.
#' @param var_spec Variable metadata specification. See \code{\link{scm_coxph}} for
#'   accepted formats.
#' @param scale_to_median Logical. If TRUE, continuous columns are divided by their
#'   median. Default \code{FALSE}.
#' @param center Logical. If TRUE, continuous covariates are median-centered.
#'   Default \code{TRUE}.
#' @param adjust_for Character vector of column names forced into the base model
#'   and protected from removal. Default \code{NULL}.
#' @param start_covariates Character vector of covariates (a subset of
#'   \code{covariates}) to pre-load into the model before forward selection
#'   begins, instead of starting from an intercept-only base. Required when
#'   \code{direction = "backward"}. Default \code{NULL} (start empty).
#' @param direction Character string controlling the stepwise search
#'   direction: one of \code{"both"} (forward selection followed by backward
#'   elimination, the default), \code{"forward"} (forward selection only), or
#'   \code{"backward"} (backward elimination only, starting from
#'   \code{start_covariates}).
#' @param interactions Character vector or list specifying interaction terms to test.
#'   Default \code{NULL} (no interactions).
#' @param interaction_strategy Character. One of \code{"after_main"} or
#'   \code{"simultaneous"}. Default \code{"after_main"}.
#' @param p_interaction Numeric vector of length 1 or 2 for interaction thresholds.
#'   Default \code{NULL} (uses \code{significance_level}).
#' @param max_interaction_order Integer. Maximum interaction order. Default 2.
#' @param verbose Logical. Print progress messages. Default \code{FALSE}.
#' @param return_models Logical. Return intermediate fitted models. Default \code{FALSE}.
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
#'   non-linear form over linear within the closed-testing step. A non-linear
#'   transformation is adopted only when the LRT p-value vs. linear is below
#'   this level. Default \code{0.05}.
#' @param fp_powers Numeric vector. Candidate powers for the Fractional
#'   Polynomial search. Default \code{c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)}.
#' @param spline_df Integer. Degrees of freedom for restricted cubic splines
#'   when \code{"rcs"} (or \code{"rcsN"}) appears in the specification.
#'   Default \code{3L}.
#' @param ... Additional arguments passed to \code{\link{scm_glm}} or
#'   \code{\link{scm_coxph}} (e.g., \code{ties}, \code{weights} for Cox models).
#'
#' @return An object of class \code{"eres_scm_glm"} or \code{"eres_scm_coxph"},
#'   depending on \code{type}. See \code{\link{scm_glm}} and \code{\link{scm_coxph}}
#'   for full return value documentation.
#'
#' @details
#' This function performs no modelling itself -- it validates inputs, resolves the
#' \code{type} argument, and delegates to the appropriate engine. It exists to
#' provide a stable programmatic interface for tools that need to run SCM
#' repeatedly (e.g., \code{\link{scm_bootstrap}}).
#'
#' The full SCM algorithm (regardless of type) proceeds as:
#' \enumerate{
#'   \item \strong{Forward first step}: Screen exposure variables; select the best
#'     by LRT p-value. The selected exposure is protected throughout.
#'   \item \strong{Forward addition}: Iteratively add covariates that improve the
#'     model at \code{significance_level[1]}.
#'   \item \strong{Forward interaction} (if \code{interactions} provided): Test
#'     eligible interaction terms (both parents must be in the model).
#'   \item \strong{Backward elimination}: Remove terms not meeting
#'     \code{significance_level[2]}, respecting hierarchy and \code{adjust_for}
#'     protection.
#' }
#'
#' @seealso \code{\link{scm_glm}}, \code{\link{scm_coxph}}, \code{\link{scm_bootstrap}}
#'
#' @examples
#' data(demo_eres)
#' demo_eres$SEX <- as.factor(demo_eres$SEX)
#'
#' # GLM example
#' fit_glm <- run_scm(
#'   eres_data = demo_eres,
#'   exposure_vars = "CMAXss",
#'   covariates = c("AGE", "WT"),
#'   type = "glm",
#'   dv_col = "sim_status",
#'   family = "binomial",
#'   significance_level = c(0.05, 0.01)
#' )
#'
#' # Cox PH example
#' fit_cox <- run_scm(
#'   eres_data = demo_eres,
#'   exposure_vars = "CMAXss",
#'   covariates = c("AGE", "WT"),
#'   type = "coxph",
#'   dv_col = "sim_status",
#'   time_col = "sim_time",
#'   significance_level = c(0.05, 0.01)
#' )
#'
#' @family bootstrap
#' @export
run_scm <- function(eres_data,
                    exposure_vars,
                    covariates,
                    type = c("glm", "coxph"),
                    significance_level = c(0.05, 0.001),
                    dvid = NULL,
                    dvid_col = NULL,
                    dv_col = "DV",
                    time_col = "TIME",
                    family = "binomial",
                    var_spec = NULL,
                    scale_to_median = FALSE,
                    center = TRUE,
                    adjust_for = NULL,
                    start_covariates = NULL,
                    direction = "both",
                    interactions = NULL,
                    interaction_strategy = "after_main",
                    p_interaction = NULL,
                    max_interaction_order = 2L,
                    nonlinear = NULL,
                    linearity_alpha = 0.05,
                    fp_powers = c(-2, -1, -0.5, 0, 0.5, 1, 2, 3),
                    spline_df = 3L,
                    verbose = FALSE,
                    return_models = FALSE,
                    ...) {

  type <- match.arg(type)

  # --- Input validation ---

  if (!is.data.frame(eres_data) || nrow(eres_data) == 0L) {
    cli::cli_abort(
      "{.arg eres_data} must be a non-empty data.frame, not {.cls {class(eres_data)}}."
    )
  }

  if (!is.character(exposure_vars) || length(exposure_vars) == 0L) {
    cli::cli_abort("{.arg exposure_vars} must be a non-empty character vector.")
  }
  missing_exp <- setdiff(exposure_vars, colnames(eres_data))
  if (length(missing_exp) > 0L) {
    cli::cli_abort(
      "The following {.arg exposure_vars} are not columns in {.arg eres_data}: {.val {missing_exp}}."
    )
  }

  if (!is.character(covariates) || length(covariates) == 0L) {
    cli::cli_abort("{.arg covariates} must be a non-empty character vector.")
  }
  missing_cov <- setdiff(covariates, colnames(eres_data))
  if (length(missing_cov) > 0L) {
    cli::cli_abort(
      "The following {.arg covariates} are not columns in {.arg eres_data}: {.val {missing_cov}}."
    )
  }

  if (!is.numeric(significance_level) || !length(significance_level) %in% 1:2 ||
      any(significance_level <= 0) || any(significance_level >= 1)) {
    cli::cli_abort(
      "{.arg significance_level} must be a numeric vector of length 1 or 2 with values in (0, 1)."
    )
  }
  if (length(significance_level) == 1L) {
    significance_level <- rep(significance_level, 2L)
  }

  if (type == "coxph") {
    if (!is.character(time_col) || length(time_col) != 1L) {
      cli::cli_abort("{.arg time_col} must be a single character string.")
    }
    if (!time_col %in% colnames(eres_data)) {
      cli::cli_abort(
        "{.arg time_col} {.val {time_col}} is not a column in {.arg eres_data}."
      )
    }
  }

  if (type == "glm") {
    if (!is.character(family) && !inherits(family, "family")) {
      cli::cli_abort(
        "{.arg family} must be a character string (e.g. {.val \"binomial\"}) or a family object."
      )
    }
  }

  # --- Nonlinear parameter validation ---
  if (!is.numeric(linearity_alpha) || length(linearity_alpha) != 1L ||
      is.na(linearity_alpha) || linearity_alpha <= 0 || linearity_alpha >= 1) {
    cli::cli_abort(
      "{.arg linearity_alpha} must be a single numeric value in (0, 1), not {.val {linearity_alpha}}."
    )
  }
  if (!is.numeric(fp_powers) || length(fp_powers) < 1L || anyNA(fp_powers)) {
    cli::cli_abort(
      "{.arg fp_powers} must be a non-empty numeric vector without NA values."
    )
  }
  spline_df <- as.integer(spline_df)
  if (length(spline_df) != 1L || is.na(spline_df) || spline_df < 2L) {
    cli::cli_abort(
      "{.arg spline_df} must be a single integer >= 2, not {.val {spline_df}}."
    )
  }

  # --- Dispatch ---

  if (type == "glm") {
    result <- scm_glm(
      eres_data             = eres_data,
      exposure_vars         = exposure_vars,
      covariates            = covariates,
      significance_level    = significance_level,
      dvid                  = dvid,
      dvid_col              = dvid_col,
      dv_col                = dv_col,
      family                = family,
      var_spec              = var_spec,
      scale_to_median       = scale_to_median,
      center                = center,
      adjust_for            = adjust_for,
      start_covariates      = start_covariates,
      direction             = direction,
      interactions          = interactions,
      interaction_strategy  = interaction_strategy,
      p_interaction         = p_interaction,
      max_interaction_order = max_interaction_order,
      nonlinear             = nonlinear,
      linearity_alpha       = linearity_alpha,
      fp_powers             = fp_powers,
      spline_df             = spline_df,
      verbose               = verbose,
      return_models         = return_models,
      ...
    )
  } else {
    result <- scm_coxph(
      eres_data             = eres_data,
      exposure_vars         = exposure_vars,
      covariates            = covariates,
      significance_level    = significance_level,
      dvid                  = dvid,
      dvid_col              = dvid_col,
      dv_col                = dv_col,
      time_col              = time_col,
      var_spec              = var_spec,
      scale_to_median       = scale_to_median,
      center                = center,
      adjust_for            = adjust_for,
      start_covariates      = start_covariates,
      direction             = direction,
      interactions          = interactions,
      interaction_strategy  = interaction_strategy,
      p_interaction         = p_interaction,
      max_interaction_order = max_interaction_order,
      nonlinear             = nonlinear,
      linearity_alpha       = linearity_alpha,
      fp_powers             = fp_powers,
      spline_df             = spline_df,
      verbose               = verbose,
      return_models         = return_models,
      ...
    )
  }

  result
}


#' Extract selected covariates from an SCM result
#'
#' @param scm_result An object of class `"eres_scm_glm"` or `"eres_scm_coxph"`.
#' @param include_exposure Logical; if FALSE (default), exclude the exposure
#'   variable from the returned set.
#'
#' @return Character vector of selected covariate names.
#' @keywords internal
.extract_selected_terms <- function(scm_result, include_exposure = FALSE) {
  meta <- scm_result[["metadata"]]

  # Early termination: no significant exposure was found
  if (is.null(meta) || is.null(meta[["selected_terms"]])) {
    return(character(0L))
  }

  terms <- meta[["selected_terms"]]

  if (!include_exposure) {
    exposure_vars <- meta[["candidate_exposure_vars"]]
    if (!is.null(exposure_vars)) {
      terms <- setdiff(terms, exposure_vars)
    }
  }

  terms
}
