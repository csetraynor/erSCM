# Nonlinear-aware forward step functions for scm_glm() and scm_coxph().
# When nonlinear = NULL these functions are never called; the existing linear
# forward steps are used instead, preserving exact backward compatibility.

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Build the response LHS string for a model type
#' @keywords internal
.response_lhs <- function(model_type) {
  if (model_type == "coxph") "survival::Surv(TIME, DV)" else "DV"
}


#' Extract base variable name(s) from a formula term string
#'
#' Given a formula fragment such as \code{"log(CMAX)"}, \code{"I(AGE^(-1))"},
#' or \code{"splines::ns(WT, df = 3)"}, extracts the underlying variable name.
#' For simple names, returns the input unchanged.
#'
#' @param term Character, a formula fragment.
#' @return Character, the base variable name.
#' @keywords internal
.term_to_varname <- function(term) {
  # splines::ns(VAR, ...) or ns(VAR, ...)
  m <- regmatches(term, regexpr("(?:splines::)?ns\\(([^,)]+)", term, perl = TRUE))
  if (length(m) > 0L && nchar(m) > 0L) {
    return(trimws(sub("(?:splines::)?ns\\(", "", m, perl = TRUE)))
  }
  # log(VAR), sqrt(VAR)
  m2 <- regmatches(term, regexpr("(?:log|sqrt)\\(([^)]+)\\)", term, perl = TRUE))
  if (length(m2) > 0L && nchar(m2) > 0L) {
    inner <- sub("(?:log|sqrt)\\(([^)]+)\\)", "\\1", m2, perl = TRUE)
    return(trimws(inner))
  }
  # I(VAR^...) or I(VAR^... * log(VAR))
  m3 <- regmatches(term, regexpr("I\\(([A-Za-z_.][A-Za-z0-9_.]*)", term, perl = TRUE))
  if (length(m3) > 0L && nchar(m3) > 0L) {
    inner <- sub("I\\(", "", m3)
    return(trimws(inner))
  }
  # plain name
  trimws(term)
}


#' Get variables currently in the model from the transform_registry
#'
#' @param current_model List, the current model state. May have a
#'   \code{transform_registry} element.
#' @return Character vector of base variable names in the model.
#' @keywords internal
.vars_in_model <- function(current_model) {
  if (!is.null(current_model$transform_registry)) {
    names(current_model$transform_registry)
  } else {
    # Fallback: parse Relation terms
    unique(vapply(current_model$Relation, .term_to_varname,
                  character(1L), USE.NAMES = FALSE))
  }
}


#' Build nobs from a fitted model
#' @keywords internal
.model_nobs <- function(fit, model_type = "glm") {
  if (is.null(fit)) return(NA_integer_)
  tryCatch(
    if (model_type == "coxph") as.integer(fit$n) else as.integer(fit$df.null + 1L),
    error = function(e) NA_integer_
  )
}


#' Build the running step table used throughout the SCM
#'
#' Constructs one row per candidate, mirroring the format of toto_2 in the
#' existing forward functions, and appends it to \code{current_model$toto_tot}.
#'
#' @param results List of per-candidate results from .select_best_transform().
#' @param current_model The current model state list.
#' @param significance_level Numeric entry threshold.
#' @param model_type Character.
#' @param data data.frame (for nobs).
#' @param family Family (GLM only).
#' @return A list with \code{toto_2} (new rows) and \code{toto_tot} (full history).
#' @keywords internal
.build_nonlinear_step_table <- function(results, current_model,
                                         significance_level, model_type,
                                         data, family) {
  tmp_ref  <- current_model$finalforward
  toto_old <- current_model$toto_tot

  # Fit the null model (current terms, no candidate) to get baseline AIC/LL
  base_rhs <- if (length(current_model$Relation) > 0L) {
    paste(current_model$Relation, collapse = " + ")
  } else {
    "1"
  }
  lhs      <- .response_lhs(model_type)
  fit_null <- .fit_model(lhs, base_rhs, data, model_type, family)
  null_ll  <- if (!is.null(fit_null)) as.numeric(stats::logLik(fit_null)) else NA_real_
  null_aic <- if (!is.null(fit_null)) tryCatch(stats::AIC(fit_null), error = function(e) NA_real_) else NA_real_
  nobs_val <- if (!is.null(fit_null)) .model_nobs(fit_null, model_type) else NA_integer_
  df_ref   <- if (!is.null(fit_null)) .model_df(fit_null, model_type) else 0L

  # Build candidate rows
  rows <- lapply(results, function(r) {
    tibble::tibble(
      Model    = r$label,
      Variable = r$var_name,
      Transform = r$best_transform,
      pvalue   = r$pvalue_vs_null,
      AIC      = r$aic_best,
      min2LL   = if (!is.na(r$aic_best)) r$aic_best - 2 * r$df else NA_real_,
      diff2LL  = if (!is.na(null_ll) && !is.na(r$aic_best))
        -2 * null_ll - (r$aic_best - 2 * r$df)
      else NA_real_,
      nobs     = nobs_val,
      dftest   = as.integer(r$df),
      dftot    = df_ref + as.integer(r$df)
    )
  })
  outtab <- dplyr::bind_rows(rows)

  # Determine which candidate to select (lowest p-value, subject to threshold,
  # ties broken by AIC as in the existing forward step)
  sig_rows <- outtab[!is.na(outtab$pvalue) & outtab$pvalue < significance_level, ]
  if (nrow(sig_rows) > 0L) {
    select_label <- sig_rows$Model[which.min(sig_rows$AIC)]
  } else {
    select_label <- NA_character_
  }

  base_row <- tibble::tibble(
    Model     = "Base",
    Variable  = NA_character_,
    Transform = NA_character_,
    pvalue    = NA_real_,
    AIC       = null_aic,
    min2LL    = -2 * null_ll,
    diff2LL   = 0,
    nobs      = nobs_val,
    dftest    = 0L,
    dftot     = df_ref
  )

  toto_2 <- dplyr::bind_rows(
    base_row,
    outtab %>%
      dplyr::mutate(
        Model    = as.character(Model),
        SELECT   = select_label,
        Signif   = as.numeric(!is.na(pvalue) & pvalue < significance_level),
        Reference = tmp_ref
      )
  ) %>%
    dplyr::mutate(
      Run = -1L + (max(toto_old$Run, na.rm = TRUE) + 1L) :
        (max(toto_old$Run, na.rm = TRUE) + dplyr::n())
    ) %>%
    dplyr::relocate(Run) %>%
    tidyr::fill(SELECT, .direction = "downup") %>%
    dplyr::mutate(
      RUNTEST = ifelse(Run != 1L,
                       paste0("Run", Reference, " + ", Model),
                       NA_character_)
    )

  toto_2 <- toto_2[-1L, ]   # drop the Base placeholder row
  toto_tot <- dplyr::bind_rows(toto_old, toto_2)

  list(toto_2 = toto_2, toto_tot = toto_tot, select_label = select_label)
}


# ---------------------------------------------------------------------------
# First forward step -- exposure screening
# ---------------------------------------------------------------------------

#' First Forward Step with Non-Linear Transformation Search (Exposure)
#'
#' Enhanced first forward step that tests transformations of the exposure
#' variable. For each exposure candidate, all allowed transformation forms are
#' tested and the best (variable, transformation) combination is selected.
#'
#' @param eres_data data.frame. The (already prepared) modelling dataset.
#'   For GLM: must contain column \code{DV}. For CoxPH: must contain
#'   \code{TIME} and \code{DV}.
#' @param exposure_vars Character vector of candidate exposure variable names.
#' @param nonlinear_spec A parsed \code{"nonlinear_spec"} object from
#'   \code{\link{parse_nonlinear_spec}}.
#' @param model_type Character, \code{"glm"} or \code{"coxph"}.
#' @param family Character or family object (GLM only).
#' @param adjust_for Character vector of terms forced into every model.
#' @param significance_level Numeric. Entry p-value threshold.
#' @param linearity_alpha Numeric. Threshold for preferring non-linear over
#'   linear within the closed-testing step.
#' @param fp_powers Numeric vector. FP candidate powers.
#' @param spline_df Integer. Default spline degrees of freedom.
#' @param verbose Logical.
#'
#' @return A list compatible with the output of
#'   \code{\link{scm_forward_first_step}} / \code{scm_forward_first_step_coxph},
#'   plus:
#'   \describe{
#'     \item{\code{transform_registry}}{Named list tracking how each selected
#'       variable entered the model.}
#'     \item{\code{transform_info}}{The winning candidate's transform details.}
#'   }
#'
#' @keywords internal
scm_forward_first_step_nonlinear <- function(eres_data,
                                              exposure_vars,
                                              nonlinear_spec,
                                              model_type   = "glm",
                                              family       = "binomial",
                                              adjust_for   = NULL,
                                              significance_level = 0.05,
                                              linearity_alpha    = 0.05,
                                              fp_powers    = .FP_POWERS,
                                              spline_df    = 3L,
                                              verbose      = FALSE) {
  lhs <- .response_lhs(model_type)

  base_rhs <- if (!is.null(adjust_for) && length(adjust_for) > 0L) {
    paste(adjust_for, collapse = " + ")
  } else {
    "1"
  }

  # For the first step we initialise a dummy current_model with just a base row
  fit_null  <- .fit_model(lhs, base_rhs, eres_data, model_type, family)
  null_ll   <- if (!is.null(fit_null)) as.numeric(stats::logLik(fit_null)) else NA_real_
  null_aic  <- if (!is.null(fit_null)) tryCatch(stats::AIC(fit_null), error = function(e) NA_real_) else NA_real_
  nobs_val  <- .model_nobs(fit_null, model_type)
  df_ref    <- .model_df(fit_null, model_type)

  # --- Test every exposure x transform combination ---
  results <- vector("list", length(exposure_vars))
  for (i in seq_along(exposure_vars)) {
    ev <- exposure_vars[[i]]
    transforms <- if (!is.null(nonlinear_spec) && ev %in% names(nonlinear_spec)) {
      nonlinear_spec[[ev]]
    } else {
      "linear"
    }

    res <- tryCatch(
      .select_best_transform(
        var_name       = ev,
        transforms     = transforms,
        current_rhs    = base_rhs,
        response_lhs   = lhs,
        data           = eres_data,
        model_type     = model_type,
        family         = family,
        fp_powers      = fp_powers,
        spline_df      = spline_df,
        linearity_alpha = linearity_alpha
      ),
      error = function(e) NULL
    )
    if (!is.null(res)) {
      results[[i]] <- res
      if (verbose) {
        message(sprintf(
          "  Exposure %s: best=%s (pval=%.4g, AIC=%.1f)",
          ev, res$best_transform,
          ifelse(is.na(res$pvalue_vs_null), NA, res$pvalue_vs_null),
          ifelse(is.na(res$aic_best), NA, res$aic_best)
        ))
      }
    }
  }
  results <- Filter(Negate(is.null), results)

  # --- Build the step tracking table ---
  # We need a minimal current_model for .build_nonlinear_step_table
  toto_init <- tibble::tibble(
    Run = 1L, Model = "Base", Variable = NA_character_, Transform = NA_character_,
    pvalue = NA_real_, Signif = NA_real_, AIC = null_aic,
    min2LL = if (!is.na(null_ll)) -2 * null_ll else NA_real_,
    diff2LL = 0, nobs = nobs_val, dftest = 0L, dftot = df_ref,
    SELECT = NA_character_, Reference = NA_real_, RUNTEST = NA_character_
  )
  dummy_current <- list(
    Relation        = if (!is.null(adjust_for) && length(adjust_for) > 0L) adjust_for else character(0L),
    toto_tot        = toto_init,
    finalforward    = 1L,
    transform_registry = list()
  )

  step_tbl <- .build_nonlinear_step_table(
    results           = results,
    current_model     = dummy_current,
    significance_level = significance_level,
    model_type        = model_type,
    data              = eres_data,
    family            = family
  )

  # --- Determine whether a variable was selected ---
  select_label <- step_tbl$select_label
  winner       <- NULL
  if (!is.na(select_label)) {
    winner <- Filter(function(r) r$label == select_label, results)
    winner <- if (length(winner) > 0L) winner[[1L]] else NULL
  }

  selected_terms_formula <- if (!is.null(winner)) winner$formula_rhs else character(0L)
  new_relation <- c(
    if (!is.null(adjust_for) && length(adjust_for) > 0L) adjust_for else character(0L),
    selected_terms_formula
  )

  transform_registry <- list()
  if (!is.null(winner)) {
    transform_registry[[winner$var_name]] <- list(
      transform_type = winner$best_transform,
      formula_rhs    = winner$formula_rhs,
      df             = winner$df,
      label          = winner$label,
      powers         = winner$powers,
      knots          = winner$knots,
      pvalue_vs_linear = winner$pvalue_vs_linear
    )
  }

  flagforst   <- !is.null(winner)
  finalforward <- if (flagforst) {
    max(step_tbl$toto_2$Run[step_tbl$toto_2$Model == select_label], na.rm = TRUE)
  } else {
    1L
  }

  if (!flagforst) {
    message("Exposure relationship not significant at the selected level.")
  }

  list(
    step1_sum          = step_tbl$toto_2,
    toto_tot           = step_tbl$toto_tot,
    Relation           = new_relation,
    finalforward       = finalforward,
    flagforst          = flagforst,
    transform_registry = transform_registry,
    transform_info     = winner
  )
}


# ---------------------------------------------------------------------------
# Main forward addition step
# ---------------------------------------------------------------------------

#' Forward Addition Step with Non-Linear Transformation Search
#'
#' Enhanced version of the forward addition step that tests multiple functional
#' forms for each candidate variable. At each step the best transformation for
#' each candidate is determined via closed testing, and the overall best
#' candidate-in-best-form is added to the model.
#'
#' Called internally by \code{\link{scm_glm}} and \code{\link{scm_coxph}} when
#' \code{nonlinear} is non-NULL. When \code{nonlinear = NULL} the standard
#' linear step is used instead.
#'
#' @param eres_data data.frame. The modelling dataset.
#' @param current_model List. The current model state:
#'   \code{toto_tot}, \code{Relation} (formula terms), \code{finalforward},
#'   \code{flagforst}, and optionally \code{transform_registry}.
#' @param covariates Character vector. All candidate covariate base names.
#' @param nonlinear_spec A parsed \code{"nonlinear_spec"} object.
#' @param model_type Character, \code{"glm"} or \code{"coxph"}.
#' @param family Character or family object (GLM only).
#' @param significance_level Numeric. Entry p-value threshold.
#' @param linearity_alpha Numeric. Closed-testing threshold for non-linearity.
#' @param fp_powers Numeric vector. FP candidate powers.
#' @param spline_df Integer. Default spline degrees of freedom.
#' @param verbose Logical.
#'
#' @return A list with the same structure as \code{\link{scm_forward_add}} /
#'   \code{scm_forward_add_coxph}, plus:
#'   \describe{
#'     \item{\code{transform_registry}}{Updated registry mapping variable names
#'       to their selected transformation details.}
#'     \item{\code{transform_info}}{Details for the variable added this step, or
#'       \code{NULL} if nothing was added.}
#'   }
#'
#' @keywords internal
scm_forward_add_nonlinear <- function(eres_data,
                                       current_model,
                                       covariates,
                                       nonlinear_spec,
                                       model_type   = "glm",
                                       family       = "binomial",
                                       significance_level = 0.05,
                                       linearity_alpha    = 0.05,
                                       fp_powers    = .FP_POWERS,
                                       spline_df    = 3L,
                                       verbose      = FALSE) {
  lhs <- .response_lhs(model_type)

  # Current model RHS for use as baseline
  current_rhs <- if (length(current_model$Relation) > 0L) {
    paste(current_model$Relation, collapse = " + ")
  } else {
    "1"
  }

  # Determine which variables are already in the model (by base name)
  vars_in <- .vars_in_model(current_model)

  # Already-selected variable names (from toto SELECT column -- same semantics
  # as the linear forward step to avoid re-testing rejected candidates)
  already_selected_labels <- unique(stats::na.omit(current_model$toto_tot$SELECT))
  # Map labels back to variable names via transform_registry
  reg   <- current_model$transform_registry
  registry_labels <- vapply(reg, function(x) x$label, character(1L))
  selected_from_reg <- names(reg)[registry_labels %in% already_selected_labels]

  # Variables to test: not yet in model AND not already selected/rejected
  candidates <- setdiff(covariates, union(vars_in, selected_from_reg))
  # Also exclude anything whose formula terms appear in Relation (linear path)
  candidates <- candidates[
    !candidates %in% vapply(current_model$Relation, .term_to_varname,
                             character(1L), USE.NAMES = FALSE)
  ]

  if (length(candidates) == 0L) {
    # Nothing left to test
    return(list(
      toto_tot           = current_model$toto_tot,
      Relation           = current_model$Relation,
      finalforward       = current_model$finalforward,
      flagforst          = FALSE,
      transform_registry = current_model$transform_registry,
      transform_info     = NULL
    ))
  }

  # --- Test each candidate ---
  results <- vector("list", length(candidates))
  for (i in seq_along(candidates)) {
    v <- candidates[[i]]
    transforms <- if (!is.null(nonlinear_spec) && v %in% names(nonlinear_spec)) {
      nonlinear_spec[[v]]
    } else {
      "linear"
    }

    res <- tryCatch(
      .select_best_transform(
        var_name        = v,
        transforms      = transforms,
        current_rhs     = current_rhs,
        response_lhs    = lhs,
        data            = eres_data,
        model_type      = model_type,
        family          = family,
        fp_powers       = fp_powers,
        spline_df       = spline_df,
        linearity_alpha = linearity_alpha
      ),
      error = function(e) NULL
    )

    if (!is.null(res)) {
      results[[i]] <- res
      if (verbose) {
        message(sprintf(
          "  Testing %s: best=%s (pval=%.4g, AIC=%.1f)",
          v, res$best_transform,
          ifelse(is.na(res$pvalue_vs_null), NA, res$pvalue_vs_null),
          ifelse(is.na(res$aic_best), NA, res$aic_best)
        ))
      }
    }
  }
  results <- Filter(Negate(is.null), results)

  if (length(results) == 0L) {
    return(list(
      toto_tot           = current_model$toto_tot,
      Relation           = current_model$Relation,
      finalforward       = current_model$finalforward,
      flagforst          = FALSE,
      transform_registry = current_model$transform_registry,
      transform_info     = NULL
    ))
  }

  # --- Build step tracking table ---
  step_tbl <- .build_nonlinear_step_table(
    results            = results,
    current_model      = current_model,
    significance_level = significance_level,
    model_type         = model_type,
    data               = eres_data,
    family             = family
  )

  select_label <- step_tbl$select_label
  winner       <- NULL
  if (!is.na(select_label)) {
    winner <- Filter(function(r) r$label == select_label, results)
    winner <- if (length(winner) > 0L) winner[[1L]] else NULL
  }

  # --- Update Relation and transform_registry ---
  new_relation <- current_model$Relation
  new_registry <- if (!is.null(current_model$transform_registry)) {
    current_model$transform_registry
  } else {
    list()
  }

  if (!is.null(winner)) {
    new_relation <- c(new_relation, winner$formula_rhs)
    new_registry[[winner$var_name]] <- list(
      transform_type   = winner$best_transform,
      formula_rhs      = winner$formula_rhs,
      df               = winner$df,
      label            = winner$label,
      powers           = winner$powers,
      knots            = winner$knots,
      pvalue_vs_linear = winner$pvalue_vs_linear
    )
    if (verbose) {
      message(sprintf("  Selected: %s", winner$label))
    }
  }

  flagforst    <- !is.null(winner)
  finalforward <- if (flagforst) {
    max(step_tbl$toto_2$Run[step_tbl$toto_2$Model == select_label], na.rm = TRUE)
  } else {
    current_model$finalforward
  }

  list(
    toto_tot           = step_tbl$toto_tot,
    Relation           = new_relation,
    finalforward       = finalforward,
    flagforst          = flagforst,
    transform_registry = new_registry,
    transform_info     = winner
  )
}
