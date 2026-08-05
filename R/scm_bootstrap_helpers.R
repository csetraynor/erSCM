# Internal helpers for scm_bootstrap()
# These functions are not exported and are prefixed with . to signal internal use.

# ---- Resampling ---------------------------------------------------------------

#' Generate a bootstrap sample
#'
#' @param data data.frame to resample.
#' @param subject_col Character or NULL. If provided, resample at subject level.
#'
#' @return A resampled data.frame.
#' @keywords internal
.generate_bootstrap_sample <- function(data, subject_col = NULL) {
  if (is.null(subject_col)) {
    idx <- sample.int(nrow(data), replace = TRUE)
    data[idx, , drop = FALSE]
  } else {
    ids <- unique(data[[subject_col]])
    sampled_ids <- sample(ids, length(ids), replace = TRUE)
    row_list <- lapply(sampled_ids, function(id) {
      data[data[[subject_col]] == id, , drop = FALSE]
    })
    do.call(rbind, row_list)
  }
}

# ---- Settings extraction -------------------------------------------------------

#' Extract SCM settings from a fitted object for bootstrap replay
#'
#' Pulls all necessary configuration from an \code{eres_scm_glm} or
#' \code{eres_scm_coxph} object to replay the SCM procedure on a new dataset.
#'
#' @param scm_fit A fitted SCM object.
#' @param significance_level Optional override for entry/retention thresholds.
#'   If \code{NULL}, attempts to recover from the object or falls back to the
#'   package default \code{c(0.05, 0.001)}.
#'
#' @return A named list of settings ready to pass to \code{run_scm()}.
#' @keywords internal
.extract_scm_settings <- function(scm_fit, significance_level = NULL) {
  meta <- scm_fit[["metadata"]]

  type <- if (inherits(scm_fit, "eres_scm_glm")) "glm" else "coxph"

  # ---- data ----
  data <- scm_fit[["model_data"]]
  if (is.null(data)) {
    cli::cli_abort(
      "{.arg scm_fit} does not contain {.field model_data}. Re-run the SCM with {.code return_models = TRUE} or ensure the fit completed successfully."
    )
  }

  # ---- exposure_vars ----
  exposure_vars <- meta[["candidate_exposure_vars"]]
  if (is.null(exposure_vars) || length(exposure_vars) == 0L) {
    cli::cli_abort(
      "Cannot extract {.field candidate_exposure_vars} from {.arg scm_fit$metadata}. The original SCM may have terminated early with no significant exposure."
    )
  }

  # ---- covariates ----
  # For glm:  candidate_covariates = covariates only
  # For coxph: candidate_covariates = c(exposure_vars, covariates) combined
  all_candidates <- meta[["candidate_covariates"]]
  if (is.null(all_candidates)) {
    cli::cli_abort(
      "Cannot extract {.field candidate_covariates} from {.arg scm_fit$metadata}."
    )
  }
  adjust_for <- meta[["adjust_for"]]
  # Remove exposure vars and adjust_for to get the pure covariate candidate list
  covariates <- setdiff(all_candidates, c(exposure_vars, adjust_for))

  # ---- significance_level ----
  if (!is.null(significance_level)) {
    sig_level <- significance_level  # already validated and recycled by caller
  } else {
    # Not stored in metadata; use package default and inform the user
    sig_level <- c(0.05, 0.001)
    cli::cli_inform(
      c(
        "i" = "{.field significance_level} is not stored in the fitted SCM object.",
        " " = "Using the default {.val c(0.05, 0.001)} for bootstrap replicates.",
        " " = "Override via the {.arg significance_level} argument if the original fit used different thresholds."
      )
    )
  }

  # ---- family (GLM only) ----
  family <- "binomial"
  if (type == "glm") {
    fm <- scm_fit[["final_model"]]
    family <- tryCatch(
      fm[["family"]][["family"]],
      error = function(e) "binomial"
    )
    if (is.null(family) || !is.character(family)) family <- "binomial"
  }

  # ---- selected_exposure (from the original fit) ----
  # New field added in Phase 5; fall back gracefully for objects without it.
  selected_exposure <- meta[["selected_exposure"]]
  if (is.null(selected_exposure)) {
    # Legacy object: infer from selected_terms (first exposure_var present)
    orig_terms <- meta[["selected_terms"]]
    found <- intersect(exposure_vars, orig_terms)
    selected_exposure <- if (length(found) > 0L) found[[1L]] else NA_character_
  }

  # ---- original selected terms (covariates only) ----
  original_selected_covariates <- meta[["selected_terms"]]
  if (is.null(original_selected_covariates)) original_selected_covariates <- character(0L)
  # Exclude exposure from the reported covariate selection
  original_selected_covariates <- setdiff(original_selected_covariates, exposure_vars)

  # ---- nonlinear settings (Phase 4) ----
  nonlinear_spec  <- meta[["nonlinear_spec"]]
  linearity_alpha <- meta[["linearity_alpha"]]
  fp_powers       <- meta[["fp_powers"]]
  spline_df       <- meta[["spline_df"]]

  # ---- start_covariates / direction (Phase 5) ----
  start_covariates <- meta[["start_covariates"]]
  direction        <- meta[["direction"]]
  if (is.null(direction)) direction <- "both"

  list(
    type                  = type,
    data                  = data,
    exposure_vars         = exposure_vars,
    covariates            = covariates,
    adjust_for            = adjust_for,
    start_covariates      = start_covariates,
    direction             = direction,
    significance_level    = sig_level,
    family                = family,
    interactions          = NULL,       # not stored in metadata; Phase 1 default
    interaction_strategy  = "after_main",
    p_interaction         = NULL,
    max_interaction_order = 2L,
    selected_exposure     = selected_exposure,
    original_selected     = original_selected_covariates,
    # nonlinear (NULL for legacy objects without NL search)
    nonlinear_spec        = nonlinear_spec,
    linearity_alpha       = linearity_alpha,
    fp_powers             = fp_powers,
    spline_df             = spline_df
  )
}

# ---- Single-replicate term extractor -------------------------------------------

#' Extract selected terms from a single bootstrap SCM result
#'
#' @param scm_result A fitted SCM object (or NULL if failed).
#' @param exposure_vars Character vector of exposure variable name(s) to exclude
#'   from the covariates result.
#'
#' @return A list with \code{exposure} (selected exposure or \code{NA_character_}),
#'   \code{exposure_transform} (character, e.g. \code{"log"} or \code{"linear"},
#'   \code{NA} when no NL search was used), \code{covariates} (character vector of
#'   selected base-variable names), and \code{covariate_transforms} (named list of
#'   transform type strings, or \code{NULL} when no NL search was used), or
#'   \code{NULL} if the result is malformed.
#' @keywords internal
.extract_boot_terms <- function(scm_result, exposure_vars) {
  if (is.null(scm_result)) return(NULL)
  meta <- scm_result[["metadata"]]
  # Early termination: metadata absent (GLM) or selected_terms absent (CoxPH)
  if (is.null(meta) || is.null(meta[["selected_terms"]])) {
    return(list(
      exposure             = NA_character_,
      exposure_transform   = NA_character_,
      covariates           = character(0L),
      covariate_transforms = NULL
    ))
  }
  terms <- meta[["selected_terms"]]

  # When nonlinear is used, selected_terms / selected_exposure contain formula
  # RHS strings like "log(CMAX)". Normalise to base variable names using
  # .term_to_varname() for all comparisons and for the returned covariates.
  terms_base <- vapply(terms, .term_to_varname, character(1L), USE.NAMES = FALSE)

  # Extract selected exposure: prefer the stored field, fall back to intersection.
  sel_exp_raw <- meta[["selected_exposure"]]
  if (is.null(sel_exp_raw)) {
    found <- intersect(exposure_vars, terms_base)
    sel_exp_raw <- if (length(found) > 0L) found[[1L]] else NA_character_
  }
  # Normalise: always store as base variable name so downstream can match it.
  sel_exp <- if (!is.na(sel_exp_raw)) .term_to_varname(sel_exp_raw) else NA_character_

  # Covariate names: base-variable names, excluding exposure
  selected_covs <- setdiff(terms_base, exposure_vars)

  # Extract nonlinear transform information from transform_registry (Phase 4).
  # transform_registry is keyed by base variable name -> transform_type string.
  tr_reg <- meta[["transform_registry"]]
  if (!is.null(tr_reg) && length(tr_reg) > 0L) {
    # Exposure transform (e.g. "log", "fp2", "linear")
    exp_transform <- if (!is.na(sel_exp) && sel_exp %in% names(tr_reg)) {
      tr_reg[[sel_exp]]$transform_type
    } else {
      "linear"
    }
    # Covariate transforms: named list, one entry per selected base-var
    cov_transforms <- if (length(selected_covs) > 0L) {
      stats::setNames(
        lapply(selected_covs, function(v) {
          if (v %in% names(tr_reg)) tr_reg[[v]]$transform_type else "linear"
        }),
        selected_covs
      )
    } else {
      list()
    }
  } else {
    exp_transform  <- NA_character_
    cov_transforms <- NULL
  }

  list(
    exposure             = sel_exp,
    exposure_transform   = exp_transform,
    covariates           = selected_covs,
    covariate_transforms = cov_transforms
  )
}

# ---- Categorical safeguard helpers (Phase 6) ----------------------------------

#' Pre-bootstrap diagnostic for factor level risk
#'
#' Scans candidate covariates for factor levels that may cause complete
#' separation or convergence failures during bootstrap resampling. Returns a
#' structured report that can be used to emit pre-bootstrap warnings.
#'
#' @param data The modelling dataset.
#' @param covariates Character vector of candidate covariate names.
#' @param exposure_vars Character vector of exposure variable name(s).
#' @param dv_col Character, name of the dependent variable column.
#' @param n_boot Integer, number of bootstrap replicates planned (informational).
#' @param subject_col Character or NULL, subject ID column for clustered data.
#'
#' @return A named list with:
#'   \item{problematic_levels}{Tibble with columns \code{variable}, \code{level},
#'     \code{n_obs}, \code{n_events}, \code{issue} ("rare_level", "empty_level",
#'     or "no_events").}
#'   \item{recommendations}{Character vector of plain-text recommendations.}
#'   \item{risk_level}{Character, one of \code{"low"}, \code{"moderate"},
#'     \code{"high"}.}
#' @keywords internal
.check_bootstrap_suitability <- function(data, covariates, exposure_vars,
                                          dv_col, n_boot, subject_col = NULL) {
  all_vars    <- c(exposure_vars, covariates)
  has_dv      <- dv_col %in% colnames(data)

  factor_vars <- all_vars[vapply(all_vars, function(v) {
    if (!v %in% colnames(data)) return(FALSE)
    is.factor(data[[v]]) || is.character(data[[v]])
  }, logical(1L))]

  if (length(factor_vars) == 0L) {
    return(list(
      problematic_levels = tibble::tibble(
        variable = character(0L), level = character(0L),
        n_obs    = integer(0L),   n_events = integer(0L),
        issue    = character(0L)
      ),
      recommendations = character(0L),
      risk_level      = "low"
    ))
  }

  rows <- list()
  for (v in factor_vars) {
    col <- data[[v]]
    if (is.character(col)) col <- factor(col)
    lvls <- levels(col)
    dv   <- if (has_dv) data[[dv_col]] else NULL

    for (lv in lvls) {
      mask  <- !is.na(col) & col == lv
      n_obs <- sum(mask)

      n_events <- if (has_dv && !is.null(dv)) {
        sum(as.integer(dv[mask]) == 1L, na.rm = TRUE)
      } else {
        NA_integer_
      }

      issue <- if (n_obs == 0L) {
        "empty_level"
      } else if (n_obs < 5L) {
        "rare_level"
      } else if (!is.na(n_events) && n_events == 0L) {
        "no_events"
      } else {
        NA_character_
      }

      if (!is.na(issue)) {
        rows[[length(rows) + 1L]] <- list(
          variable = v,
          level    = as.character(lv),
          n_obs    = as.integer(n_obs),
          n_events = if (is.na(n_events)) NA_integer_ else as.integer(n_events),
          issue    = issue
        )
      }
    }
  }

  prob_tbl <- if (length(rows) > 0L) {
    tibble::tibble(
      variable = vapply(rows, `[[`, character(1L), "variable"),
      level    = vapply(rows, `[[`, character(1L), "level"),
      n_obs    = vapply(rows, `[[`, integer(1L),   "n_obs"),
      n_events = vapply(rows, function(r) {
        v <- r[["n_events"]]
        if (is.na(v)) NA_integer_ else as.integer(v)
      }, integer(1L)),
      issue    = vapply(rows, `[[`, character(1L), "issue")
    )
  } else {
    tibble::tibble(
      variable = character(0L), level = character(0L),
      n_obs    = integer(0L),   n_events = integer(0L),
      issue    = character(0L)
    )
  }

  n_prob_vars <- if (nrow(prob_tbl) > 0L) {
    length(unique(prob_tbl$variable[prob_tbl$issue %in% c("empty_level", "rare_level")]))
  } else {
    0L
  }

  risk_level <- if (nrow(prob_tbl) == 0L) {
    "low"
  } else if (n_prob_vars >= 2L) {
    "high"
  } else {
    "moderate"
  }

  recommendations <- character(0L)
  if (nrow(prob_tbl) > 0L) {
    prob_vars <- unique(prob_tbl$variable)
    recommendations <- paste0(
      "Consider combining or removing rare levels in: ",
      paste(prob_vars, collapse = ", ")
    )
    if (risk_level == "high") {
      recommendations <- c(
        recommendations,
        "Multiple variables have problematic levels; bootstrap results may be unreliable."
      )
    }
  }

  list(
    problematic_levels = prob_tbl,
    recommendations    = recommendations,
    risk_level         = risk_level
  )
}

#' Prepare a bootstrap sample: drop empty levels and exclude sparse-level variables
#'
#' Called on each bootstrap replicate before the SCM run. Drops unused factor
#' levels via \code{droplevels()} and removes any factor/character covariate
#' where at least one remaining level has too few observations or too few events
#' (complete separation risk).
#'
#' @param boot_data The resampled data.frame.
#' @param covariates Character vector of candidate covariate names to screen.
#' @param dv_col Character, name of the dependent variable column.
#' @param min_level_n Integer. Minimum observations per factor level; levels
#'   with fewer observations trigger exclusion of the whole variable. Default 5L.
#' @param min_events_per_level Integer. Minimum events per factor level; zero
#'   events (complete separation) trigger exclusion. Default 1L. Pass 0L to
#'   disable the events check.
#' @param drop_empty_levels Logical. Call \code{droplevels()} first. Default TRUE.
#'
#' @return A named list with:
#'   \item{data}{The (possibly \code{droplevels()}-cleaned) boot sample.}
#'   \item{valid_covariates}{Character vector -- covariates safe to use.}
#'   \item{excluded_covariates}{Character vector -- covariates removed this replicate.}
#'   \item{modified_levels}{Character vector -- factor variables whose levels were
#'     dropped but were not excluded entirely.}
#' @keywords internal
.prepare_boot_sample <- function(boot_data, covariates, dv_col,
                                  min_level_n          = 5L,
                                  min_events_per_level = 1L,
                                  drop_empty_levels    = TRUE) {
  # Safeguards disabled -- return unchanged
  if (min_level_n <= 0L && min_events_per_level <= 0L) {
    return(list(
      data                = boot_data,
      valid_covariates    = covariates,
      excluded_covariates = character(0L),
      modified_levels     = character(0L)
    ))
  }

  if (drop_empty_levels) boot_data <- droplevels(boot_data)

  has_dv   <- dv_col %in% colnames(boot_data)
  dv       <- if (has_dv) boot_data[[dv_col]] else NULL
  excluded <- character(0L)
  modified <- character(0L)

  for (v in covariates) {
    if (!v %in% colnames(boot_data)) next
    col <- boot_data[[v]]
    if (!is.factor(col) && !is.character(col)) next
    if (is.character(col)) col <- factor(col)

    lvls_orig <- levels(col)

    problem <- FALSE
    for (lv in lvls_orig) {
      mask  <- !is.na(col) & col == lv
      n_obs <- sum(mask)

      if (min_level_n > 0L && n_obs < min_level_n) {
        problem <- TRUE
        break
      }

      if (min_events_per_level > 0L && has_dv && !is.null(dv)) {
        n_ev <- sum(as.integer(dv[mask]) == 1L, na.rm = TRUE)
        if (n_ev < min_events_per_level) {
          problem <- TRUE
          break
        }
      }
    }

    if (problem) {
      excluded <- c(excluded, v)
    } else {
      lvls_now <- levels(factor(col))
      if (length(lvls_now) < length(lvls_orig)) {
        modified <- c(modified, v)
      }
    }
  }

  list(
    data                = boot_data,
    valid_covariates    = setdiff(covariates, excluded),
    excluded_covariates = excluded,
    modified_levels     = modified
  )
}

# ---- Replicate runner ----------------------------------------------------------

#' Run bootstrap replicates (sequential or parallel)
#'
#' @param data The modelling dataset (already scaled/centered/subsetted).
#' @param n_boot Number of replicates.
#' @param subject_col Subject column name or NULL.
#' @param scm_settings List from \code{.extract_scm_settings()}.
#' @param parallel Logical. If TRUE, uses \code{future.apply::future_lapply()}.
#' @param seed Integer seed for reproducibility.
#' @param min_level_n Integer. Passed to \code{.prepare_boot_sample()}. Default 5L.
#' @param min_events_per_level Integer. Passed to \code{.prepare_boot_sample()}.
#'   Default 1L.
#'
#' @return A list of length \code{n_boot}, each element a list with
#'   \code{selected_terms} (character vector or NULL), \code{converged} (logical),
#'   and \code{excluded_covariates} (character vector, possibly empty).
#' @keywords internal
.run_bootstrap_replicates <- function(data, n_boot, subject_col, scm_settings,
                                      parallel, seed,
                                      min_level_n          = 5L,
                                      min_events_per_level = 1L) {
  # Single-replicate worker: generate sample, run SCM, extract terms.
  # The function is self-contained so it can be passed to future_lapply unchanged.
  worker_fn <- function(i) {
    boot_data <- .generate_bootstrap_sample(data, subject_col)

    # --- Per-replicate categorical safeguard (Phase 6) ---
    prep <- .prepare_boot_sample(
      boot_data            = boot_data,
      covariates           = scm_settings$covariates,
      dv_col               = "DV",
      min_level_n          = min_level_n,
      min_events_per_level = min_events_per_level,
      drop_empty_levels    = TRUE
    )
    boot_data     <- prep$data
    valid_covs    <- prep$valid_covariates
    excluded_covs <- prep$excluded_covariates

    # Filter start_covariates to only include variables still in valid_covs
    start_covs_rep <- if (!is.null(scm_settings$start_covariates)) {
      sc <- intersect(scm_settings$start_covariates, valid_covs)
      if (length(sc) == 0L) NULL else sc
    } else {
      NULL
    }

    # direction="backward" requires start_covariates; fall back gracefully
    direction_rep <- scm_settings$direction
    if (isTRUE(direction_rep == "backward") && is.null(start_covs_rep)) {
      direction_rep <- "both"
    }

    # Build the base argument list, then append nonlinear settings when present.
    scm_args <- list(
      eres_data             = boot_data,
      exposure_vars         = scm_settings$exposure_vars,
      covariates            = valid_covs,
      type                  = scm_settings$type,
      significance_level    = scm_settings$significance_level,
      dvid                  = NULL,
      dvid_col              = NULL,
      dv_col                = "DV",
      time_col              = "TIME",
      family                = scm_settings$family,
      var_spec              = NULL,
      scale_to_median       = FALSE,
      center                = FALSE,
      adjust_for            = scm_settings$adjust_for,
      start_covariates      = start_covs_rep,
      direction             = direction_rep,
      interactions          = scm_settings$interactions,
      interaction_strategy  = scm_settings$interaction_strategy,
      p_interaction         = scm_settings$p_interaction,
      max_interaction_order = scm_settings$max_interaction_order,
      verbose               = FALSE,
      return_models         = FALSE
    )

    # Pass nonlinear settings when the original fit used NL search (Phase 4).
    # run_scm() forwards ... to scm_glm()/scm_coxph(), so these are picked up.
    if (!is.null(scm_settings$nonlinear_spec)) {
      scm_args$nonlinear <- scm_settings$nonlinear_spec
      if (!is.null(scm_settings$linearity_alpha)) {
        scm_args$linearity_alpha <- scm_settings$linearity_alpha
      }
      if (!is.null(scm_settings$fp_powers)) {
        scm_args$fp_powers <- scm_settings$fp_powers
      }
      if (!is.null(scm_settings$spline_df)) {
        scm_args$spline_df <- scm_settings$spline_df
      }
    }

    raw_result <- tryCatch(
      withCallingHandlers(
        suppressMessages(do.call(run_scm, scm_args)),
        warning = function(w) invokeRestart("muffleWarning")
      ),
      error = function(e) NULL
    )

    boot_entry <- tryCatch(
      .extract_boot_terms(raw_result, scm_settings$exposure_vars),
      error = function(e) NULL
    )

    list(
      selected_terms      = boot_entry,
      converged           = !is.null(boot_entry),
      excluded_covariates = excluded_covs
    )
  }

  if (parallel) {
    future.apply::future_lapply(
      seq_len(n_boot),
      worker_fn,
      future.seed = seed
    )
  } else {
    set.seed(seed)
    results <- vector("list", n_boot)
    cli::cli_progress_bar("Bootstrapping SCM", total = n_boot)
    for (i in seq_len(n_boot)) {
      results[[i]] <- worker_fn(i)
      cli::cli_progress_update()
    }
    cli::cli_progress_done()
    results
  }
}

# ---- Aggregation helpers -------------------------------------------------------

# Internal: extract the covariate character vector from a boot_models element.
# Handles both the new list(exposure, covariates) format and legacy plain char vectors.
.get_covariates <- function(m) {
  if (is.list(m) && !is.null(m[["covariates"]])) m[["covariates"]] else m
}

#' Compute inclusion probabilities from bootstrap models
#'
#' @param boot_models List of entries per converged replicate; each entry is
#'   either a \code{list(exposure, covariates)} (new format) or a plain
#'   character vector (legacy / unit-test format). NULL entries must be removed
#'   before calling.
#' @param candidate_terms Character vector of all candidate term names.
#'
#' @return Named numeric vector of inclusion proportions.
#' @keywords internal
.compute_inclusion_prob <- function(boot_models, candidate_terms) {
  n <- length(boot_models)
  if (n == 0L || length(candidate_terms) == 0L) {
    return(stats::setNames(numeric(length(candidate_terms)), candidate_terms))
  }
  counts <- vapply(candidate_terms, function(term) {
    sum(vapply(boot_models, function(m) term %in% .get_covariates(m), logical(1L)))
  }, numeric(1L))
  counts / n
}

#' Compute exposure selection probabilities from bootstrap models
#'
#' @param boot_models List of \code{list(exposure, covariates)} entries (converged
#'   replicates only).
#' @param candidate_exposures Character vector of all candidate exposure names.
#'
#' @return Named numeric vector of selection proportions, one entry per candidate
#'   exposure. When some replicates converged without selecting any exposure
#'   (early-termination models), a \code{"none_significant"} entry is appended.
#' @keywords internal
.compute_exposure_selection_prob <- function(boot_models, candidate_exposures) {
  n <- length(boot_models)
  if (n == 0L || length(candidate_exposures) == 0L) {
    return(stats::setNames(numeric(length(candidate_exposures)), candidate_exposures))
  }
  counts <- vapply(candidate_exposures, function(exp) {
    sum(vapply(boot_models, function(m) {
      sel <- if (is.list(m) && !is.null(m[["exposure"]])) m[["exposure"]] else NA_character_
      isTRUE(sel == exp)
    }, logical(1L)))
  }, numeric(1L))
  probs <- counts / n

  # For multi-exposure scenarios, also record fraction of replicates where no
  # exposure was selected (early-termination / null models). Not appended for
  # single-exposure runs because the complement is already implied by 1 - probs.
  if (length(candidate_exposures) > 1L) {
    n_none <- sum(vapply(boot_models, function(m) {
      sel <- if (is.list(m) && !is.null(m[["exposure"]])) m[["exposure"]] else NA_character_
      is.na(sel)
    }, logical(1L)))
    if (n_none > 0L) {
      probs <- c(probs, none_significant = n_none / n)
    }
  }
  probs
}

#' Compute covariate inclusion probabilities conditional on exposure selection
#'
#' @param boot_models List of \code{list(exposure, covariates)} entries (converged
#'   replicates only).
#' @param candidate_terms Character vector of candidate covariate names.
#' @param candidate_exposures Character vector of all candidate exposure names.
#'
#' @return A named list of numeric vectors (one per exposure), or \code{NULL}
#'   when there is only a single candidate exposure (the unconditional
#'   \code{inclusion_prob} already gives the answer).
#' @keywords internal
.compute_inclusion_prob_by_exposure <- function(boot_models, candidate_terms,
                                                candidate_exposures) {
  if (length(candidate_exposures) <= 1L) return(NULL)
  if (length(boot_models) == 0L || length(candidate_terms) == 0L) return(NULL)

  lapply(stats::setNames(candidate_exposures, candidate_exposures), function(exp) {
    subset_models <- Filter(function(m) {
      sel <- if (is.list(m) && !is.null(m[["exposure"]])) m[["exposure"]] else NA_character_
      isTRUE(sel == exp)
    }, boot_models)
    .compute_inclusion_prob(subset_models, candidate_terms)
  })
}

#' Compute co-inclusion matrix
#'
#' @param boot_models List of character vectors (converged replicates only).
#' @param all_terms Character vector of all terms to include in the matrix.
#'
#' @return A symmetric numeric matrix of co-inclusion proportions.
#' @keywords internal
.compute_co_inclusion_matrix <- function(boot_models, all_terms) {
  n  <- length(boot_models)
  p  <- length(all_terms)
  dn <- list(all_terms, all_terms)
  mat <- matrix(0, nrow = p, ncol = p, dimnames = dn)

  if (n == 0L || p == 0L) return(mat)

  # Build binary inclusion matrix: rows = terms, cols = replicates.
  # vapply returns a vector (not matrix) when p == 1; coerce explicitly.
  incl_mat <- vapply(boot_models,
                     function(m) all_terms %in% .get_covariates(m),
                     logical(p))
  if (!is.matrix(incl_mat)) {
    incl_mat <- matrix(incl_mat, nrow = p, ncol = n)
  }
  # tcrossprod gives pxp co-occurrence counts
  co_counts <- tcrossprod(incl_mat * 1L)
  mat <- co_counts / n
  dimnames(mat) <- dn
  mat
}

#' Compute transformation selection frequencies
#'
#' For each variable that had a non-linear search in the original fit, tallies
#' how often each transformation type was selected across converged bootstrap
#' replicates. Also counts replicates where the variable was not selected at all.
#'
#' @param converged_models List of per-replicate results (with
#'   \code{covariate_transforms} and \code{exposure_transform} fields).
#' @param nonlinear_spec A \code{nonlinear_spec} object (named list of
#'   transformation type character vectors), or \code{NULL}.
#' @param exposure_vars Character vector of exposure candidate names.
#'
#' @return Named list with:
#'   \item{covariate_transform_freq}{Named list; one named numeric vector per
#'     covariate that appeared in \code{nonlinear_spec}. Entries include
#'     \code{"not_selected"}.}
#'   \item{exposure_transform_freq}{Named list; one named numeric vector per
#'     NL-enabled exposure variable. \code{NULL} if no exposure had a NL search.
#'     Each entry is named by the exposure variable name.}
#' @keywords internal
.compute_transform_freq <- function(converged_models, nonlinear_spec,
                                    exposure_vars) {
  n <- length(converged_models)

  if (is.null(nonlinear_spec) || n == 0L) {
    return(list(covariate_transform_freq = NULL,
                exposure_transform_freq  = NULL))
  }

  nl_vars     <- names(nonlinear_spec)
  cov_nl_vars <- setdiff(nl_vars, exposure_vars)
  exp_nl_vars <- intersect(nl_vars, exposure_vars)

  # --- Covariate transform frequencies ---
  cov_freq <- if (length(cov_nl_vars) > 0L) {
    stats::setNames(
      lapply(cov_nl_vars, function(v) {
        allowed_types <- unique(c(nonlinear_spec[[v]], "not_selected"))
        counts <- stats::setNames(integer(length(allowed_types)), allowed_types)
        for (m in converged_models) {
          if (!is.list(m) || is.null(m[["covariate_transforms"]])) {
            counts["not_selected"] <- counts["not_selected"] + 1L
            next
          }
          ct <- m$covariate_transforms
          if (!is.null(ct) && v %in% names(ct)) {
            ttype <- ct[[v]]
            if (ttype %in% names(counts)) {
              counts[ttype] <- counts[ttype] + 1L
            } else {
              # Transform happened but not in allowed list (shouldn't occur)
              counts["not_selected"] <- counts["not_selected"] + 1L
            }
          } else {
            counts["not_selected"] <- counts["not_selected"] + 1L
          }
        }
        counts / n
      }),
      cov_nl_vars
    )
  } else {
    NULL
  }

  # --- Exposure transform frequencies ---
  exp_freq <- if (length(exp_nl_vars) > 0L) {
    lapply(
      stats::setNames(exp_nl_vars, exp_nl_vars),
      function(ev) {
        allowed_types <- unique(c(nonlinear_spec[[ev]], "not_selected"))
        counts <- stats::setNames(integer(length(allowed_types)), allowed_types)
        for (m in converged_models) {
          if (!is.list(m)) {
            counts["not_selected"] <- counts["not_selected"] + 1L
            next
          }
          et <- m[["exposure_transform"]]
          if (!is.null(et) && !is.na(et)) {
            if (et %in% names(counts)) {
              counts[et] <- counts[et] + 1L
            } else {
              counts["not_selected"] <- counts["not_selected"] + 1L
            }
          } else {
            counts["not_selected"] <- counts["not_selected"] + 1L
          }
        }
        counts / n
      }
    )
  } else {
    NULL
  }

  list(
    covariate_transform_freq = cov_freq,
    exposure_transform_freq  = exp_freq
  )
}

#' Identify the modal (most frequent) model
#'
#' @param boot_models List of boot model entries (converged replicates only).
#'   Each entry is either a \code{list(exposure, covariates)} (new format) or a
#'   plain character vector (legacy / unit-test format).
#'
#' @return A list with \code{model} (a \code{list(exposure, covariates)} in new
#'   format, or a character vector in legacy format when the input is plain
#'   character vectors) and \code{frequency} (numeric proportion).
#' @keywords internal
.identify_modal_model <- function(boot_models) {
  if (length(boot_models) == 0L) {
    return(list(model = character(0L), frequency = NA_real_))
  }

  # Detect format: if any element is a list with 'covariates', use new format.
  is_new_format <- any(vapply(boot_models, function(m) {
    is.list(m) && !is.null(m[["covariates"]])
  }, logical(1L)))

  if (is_new_format) {
    # Canonical string: "EXPOSURE | cov1 + cov2 + ..."
    model_strings <- vapply(boot_models, function(m) {
      exp_part <- if (!is.null(m[["exposure"]]) && !is.na(m[["exposure"]])) {
        m[["exposure"]]
      } else {
        "(none)"
      }
      cov_part <- paste(sort(.get_covariates(m)), collapse = " + ")
      paste0(exp_part, " | ", cov_part)
    }, character(1L))

    freq_table <- sort(table(model_strings), decreasing = TRUE)
    modal_string <- names(freq_table)[1L]
    modal_freq   <- as.numeric(freq_table[1L]) / length(boot_models)

    # Parse back into list(exposure, covariates)
    parts <- strsplit(modal_string, " \\| ", fixed = FALSE)[[1L]]
    modal_exposure   <- parts[[1L]]
    modal_cov_string <- if (length(parts) >= 2L) parts[[2L]] else ""
    modal_covariates <- if (nchar(modal_cov_string) == 0L) {
      character(0L)
    } else {
      strsplit(modal_cov_string, " \\+ ", fixed = FALSE)[[1L]]
    }

    list(
      model     = list(exposure = modal_exposure, covariates = modal_covariates),
      frequency = modal_freq
    )
  } else {
    # Legacy: plain character vectors (unit tests)
    model_strings <- vapply(boot_models, function(m) {
      paste(sort(m), collapse = " + ")
    }, character(1L))

    freq_table <- sort(table(model_strings), decreasing = TRUE)
    modal_string <- names(freq_table)[1L]
    modal_freq   <- as.numeric(freq_table[1L]) / length(boot_models)

    modal_terms <- if (nchar(modal_string) == 0L) {
      character(0L)
    } else {
      strsplit(modal_string, " \\+ ", fixed = FALSE)[[1L]]
    }

    list(model = modal_terms, frequency = modal_freq)
  }
}
