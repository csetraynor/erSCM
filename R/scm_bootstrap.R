#' Bootstrap Model Selection Stability for SCM
#'
#' Assesses the stability of covariate selection by resampling the dataset
#' \code{n_boot} times and re-running the full SCM procedure on each bootstrap
#' sample. Reports covariate inclusion frequencies, co-inclusion patterns, and
#' model selection variability.
#'
#' @param scm_fit An object of class \code{"eres_scm_glm"} or
#'   \code{"eres_scm_coxph"}, as returned by \code{\link{scm_glm}} or
#'   \code{\link{scm_coxph}}.
#' @param n_boot Integer, number of bootstrap replicates. Default 500.
#'   Must be \eqn{\ge} 10.
#' @param seed Integer, random seed for reproducibility. Default 42.
#' @param subject_col Optional character string naming the column that identifies
#'   unique subjects for patient-level resampling. When provided, entire patients
#'   (all their rows) are resampled together -- important for multi-row-per-patient
#'   data or time-to-event analyses. If \code{NULL} (default), row-level
#'   resampling is performed.
#' @param significance_level Optional numeric vector of length 2 to override the
#'   entry/retention thresholds used in the original fit. If \code{NULL} (default),
#'   uses the same thresholds as the original SCM run (or \code{c(0.05, 0.001)} if
#'   not stored in the fitted object).
#' @param parallel Logical. If \code{TRUE}, uses
#'   \code{future.apply::future_lapply()} for parallel execution. Requires the
#'   \pkg{future.apply} package and an active \code{future::plan()}.
#'   Default \code{FALSE}.
#' @param min_level_n Integer. Minimum number of observations required per
#'   factor level in a bootstrap sample for that variable to be included in the
#'   candidate pool for that replicate. Variables with any level below this
#'   threshold are excluded (counted as "not selected") from that replicate.
#'   Default \code{5L}. Set to \code{0L} to disable the check.
#' @param min_events_per_level Integer. Minimum number of outcome events
#'   required per factor level to avoid complete separation. Variables with any
#'   level having fewer events are excluded from that replicate. Default
#'   \code{1L}. Set to \code{0L} to disable the check.
#' @param ... Reserved for future use. Currently ignored.
#'
#' @return An object of class \code{"scm_bootstrap"} containing:
#' \describe{
#'   \item{inclusion_prob}{Named numeric vector of inclusion probabilities per
#'     candidate covariate (excluding exposure and \code{adjust_for} terms).}
#'   \item{interaction_inclusion_prob}{Named numeric vector of inclusion
#'     probabilities per interaction term (empty named vector if no interactions
#'     appeared in any replicate).}
#'   \item{co_inclusion}{Symmetric matrix of pairwise co-inclusion proportions.
#'     Entry (i,j) = proportion of converged replicates where both term i and
#'     term j were selected.}
#'   \item{modal_model}{Named list with \code{exposure} (character string) and
#'     \code{covariates} (character vector) giving the most frequently selected
#'     model (ties broken by first canonical string).}
#'   \item{modal_model_freq}{Proportion of converged replicates selecting the
#'     modal model.}
#'   \item{model_variability_index}{Numeric: \code{1 - modal_model_freq}. Values
#'     near 0 indicate stable selection; near 1 indicate high instability.}
#'   \item{exposure_selection_prob}{Named numeric vector of selection probabilities
#'     per candidate exposure variable (proportion of converged replicates in which
#'     each exposure was selected). Length 1 for single-exposure analyses.}
#'   \item{exposure_selected_per_boot}{Character vector of length \code{n_converged}
#'     recording which exposure was selected in each converged replicate
#'     (\code{NA} for replicates where no exposure was selected).}
#'   \item{inclusion_prob_by_exposure}{Named list of inclusion probability vectors,
#'     one per candidate exposure variable, showing covariate inclusion conditional
#'     on each exposure being selected. \code{NULL} when only a single candidate
#'     exposure is provided.}
#'   \item{boot_models}{List of length \code{n_boot}. Each element is a named list
#'     with \code{exposure} (character) and \code{covariates} (character vector) for
#'     that replicate, or \code{NULL} if the replicate failed to converge.}
#'   \item{boot_converged}{Logical vector of length \code{n_boot}.}
#'   \item{n_boot}{Integer, replicates attempted.}
#'   \item{n_converged}{Integer, replicates that succeeded.}
#'   \item{convergence_rate}{Numeric, \code{n_converged / n_boot}.}
#'   \item{seed}{The seed used.}
#'   \item{type}{Character, \code{"glm"} or \code{"coxph"}.}
#'   \item{exposure_vars}{Character vector, the exposure variable(s) from the
#'     original fit.}
#'   \item{covariates}{Character vector, the candidate covariates from the
#'     original fit (excluding exposure and \code{adjust_for}).}
#'   \item{adjust_for}{Character vector or NULL, forced adjustment terms.}
#'   \item{significance_level}{Numeric vector of length 2, the thresholds used.}
#'   \item{original_selected}{Named list with \code{exposure} (character string,
#'     the exposure selected in the original fit) and \code{covariates} (character
#'     vector of covariates selected in the original fit, exposure excluded).}
#'   \item{convergence_issues}{Tibble with one row per failed replicate, columns
#'     \code{replicate} (integer index) and \code{reason} (\code{"error"}).}
#'   \item{variable_exclusion_rate}{Named numeric vector. For each candidate
#'     covariate, the proportion of converged replicates in which that variable
#'     was excluded due to sparse/separated factor levels. Only populated when
#'     \code{min_level_n > 0} or \code{min_events_per_level > 0} and at least one
#'     exclusion occurred; otherwise \code{NULL}.}
#'   \item{diagnostics}{A list of class \code{"bootstrap_diagnostics"} from the
#'     pre-bootstrap suitability check (see \code{.check_bootstrap_suitability()}).
#'     Contains \code{problematic_levels}, \code{recommendations}, and
#'     \code{risk_level}.}
#'   \item{call}{The matched call.}
#' }
#'
#' @section Handling Categorical Covariates with Rare Levels:
#' Bootstrap resampling can produce samples where factor levels have very few
#' (or zero) observations, leading to complete separation or convergence
#' failures. \code{scm_bootstrap()} implements a two-layer safeguard:
#' \enumerate{
#'   \item \strong{Pre-bootstrap check} (\code{.check_bootstrap_suitability()}):
#'     scans the full dataset for factor levels with \eqn{< 5} observations or
#'     zero events and emits a \code{cli_warn} ("high" risk) or \code{cli_inform}
#'     ("moderate" risk) before the bootstrap loop begins.
#'   \item \strong{Per-replicate screening} (\code{.prepare_boot_sample()}):
#'     for each bootstrap sample, any factor covariate with a level failing the
#'     \code{min_level_n} or \code{min_events_per_level} thresholds is removed
#'     from that replicate's candidate pool. Excluded variables are treated as
#'     "not selected" for inclusion probability calculations (conservative).
#' }
#' Set both \code{min_level_n = 0L} and \code{min_events_per_level = 0L} to
#' disable all safeguards (restores pre-Phase-6 behavior).
#'
#' @details
#' The function extracts all necessary information from the fitted SCM object:
#' \itemize{
#'   \item The modelling dataset (\code{scm_fit$model_data}) -- already
#'     transformed (scaled, centered, subsetted) as in the original fit.
#'   \item Candidate covariates and exposure variables from
#'     \code{scm_fit$metadata}.
#'   \item Model type (GLM vs CoxPH) from the object's class.
#'   \item Significance thresholds (can be overridden via
#'     \code{significance_level}).
#' }
#'
#' For each bootstrap replicate:
#' \enumerate{
#'   \item Resample the dataset with replacement (row-level or patient-level).
#'   \item Run the full SCM procedure using the same settings as the original fit
#'     (\code{scale_to_median = FALSE}, \code{center = FALSE}: data already
#'     transformed).
#'   \item Record which covariates (and interactions, if any) were selected.
#' }
#'
#' The exposure variable is ALWAYS protected (enters first, never removed) and
#' is excluded from inclusion probability calculations. Similarly,
#' \code{adjust_for} terms are excluded since they are forced into every model
#' by design.
#'
#' @section Parallelisation:
#' When \code{parallel = TRUE}, the function uses
#' \code{future.apply::future_lapply()} with pre-generated per-replicate
#' L'Ecuyer-CMRG streams (see the Reproducibility section) for reproducible
#' parallel random number generation. Set up a parallel backend before calling:
#' \preformatted{
#' future::plan(future::multisession, workers = 4)
#' stability <- scm_bootstrap(fit, n_boot = 500, parallel = TRUE)
#' future::plan(future::sequential)  # reset
#' }
#'
#' @section Reproducibility:
#' \code{scm_bootstrap()} returns \strong{identical resamples for the same
#' \code{seed} regardless of the \code{parallel} argument}. Internally, a
#' per-replicate seed vector of length \code{n_boot} is derived from \code{seed}
#' using an L'Ecuyer-CMRG stream: the RNG kind is switched to
#' \code{"L'Ecuyer-CMRG"}, \code{seed} is set, and successive streams are built
#' with \code{parallel::nextRNGStream()}. The sequential path assigns stream
#' \code{i} before generating replicate \code{i}; the parallel path injects the
#' same streams via \code{future.apply::future_lapply(future.seed = )}. Because
#' both paths consume the same streams, they produce bit-identical resamples.
#' The caller's RNG kind and \code{.Random.seed} are saved on entry and restored
#' on exit, so the internal L'Ecuyer-CMRG switch does not affect the caller's
#' session. Note that results from prior package versions (which seeded the
#' sequential and parallel paths independently) will not match the current
#' implementation.
#'
#' @seealso \code{\link{scm_glm}}, \code{\link{scm_coxph}},
#'   \code{\link{print.scm_bootstrap}}, \code{\link{summary.scm_bootstrap}},
#'   \code{\link{stable_covariates}}
#'
#' @examples
#' data(demo_eres)
#' fit <- scm_glm(
#'   eres_data          = demo_eres,
#'   exposure_vars      = "CMAXss",
#'   covariates         = c("AGE", "WT"),
#'   dv_col             = "sim_status",
#'   significance_level = c(0.2, 0.2)
#' )
#'
#' # A small number of replicates keeps this fast; use n_boot >= 200 in practice.
#' stability <- scm_bootstrap(fit, n_boot = 10, seed = 42)
#' print(stability)
#' stable_covariates(stability, threshold = 0.5)
#'
#' \donttest{
#' # A realistic run (n_boot = 500) is too slow for automated checks; kept
#' # wrapped for illustration only.
#' stability_full <- scm_bootstrap(fit, n_boot = 500, seed = 42)
#'
#' # Parallel execution (requires the optional future/future.apply packages)
#' if (requireNamespace("future", quietly = TRUE) &&
#'     requireNamespace("future.apply", quietly = TRUE)) {
#'   future::plan(future::multisession, workers = 4)
#'   stability_par <- scm_bootstrap(fit, n_boot = 500, parallel = TRUE)
#'   future::plan(future::sequential)
#' }
#' }
#'
#' @family bootstrap
#' @export
scm_bootstrap <- function(scm_fit,
                           n_boot               = 500L,
                           seed                 = 42L,
                           subject_col          = NULL,
                           significance_level   = NULL,
                           parallel             = FALSE,
                           min_level_n          = 5L,
                           min_events_per_level = 1L,
                           ...) {
  mc <- match.call()

  # ---- Input validation -------------------------------------------------------

  if (!inherits(scm_fit, "eres_scm")) {
    cli::cli_abort(
      c(
        "{.arg scm_fit} must be an object of class {.cls eres_scm_glm} or {.cls eres_scm_coxph}.",
        "x" = "Got an object of class {.cls {class(scm_fit)}}.",
        "i" = "Fit a model first with {.fn scm_glm} or {.fn scm_coxph}, then pass the result here."
      )
    )
  }
  if (is.null(scm_fit$final_model)) {
    cli::cli_abort(
      c(
        "Cannot bootstrap: the original fit did not produce a model.",
        "i" = "No exposure variable reached the entry significance threshold in the original fit.",
        "i" = "Bootstrap requires a successfully fitted SCM with at least one selected exposure."
      )
    )
  }

  n_boot <- tryCatch(
    as.integer(n_boot),
    warning = function(w) {
      cli::cli_abort("{.arg n_boot} must be coercible to integer; got {.val {n_boot}}.")
    }
  )
  if (length(n_boot) != 1L || is.na(n_boot) || n_boot < 10L) {
    cli::cli_abort(
      "{.arg n_boot} must be a single integer >= 10; got {.val {n_boot}}."
    )
  }

  seed <- tryCatch(
    as.integer(seed),
    warning = function(w) {
      cli::cli_abort("{.arg seed} must be coercible to a single integer.")
    }
  )
  if (length(seed) != 1L || is.na(seed)) {
    cli::cli_abort("{.arg seed} must be a single integer.")
  }

  if (!is.null(subject_col)) {
    if (!is.character(subject_col) || length(subject_col) != 1L) {
      cli::cli_abort("{.arg subject_col} must be a single character string or NULL.")
    }
    model_data <- scm_fit[["model_data"]]
    if (!subject_col %in% colnames(model_data)) {
      cli::cli_abort(
        "{.arg subject_col} {.val {subject_col}} is not a column in {.code scm_fit$model_data}."
      )
    }
  }

  if (!is.null(significance_level)) {
    if (!is.numeric(significance_level) ||
        !length(significance_level) %in% 1:2 ||
        any(significance_level <= 0) ||
        any(significance_level >= 1)) {
      cli::cli_abort(
        "{.arg significance_level} must be a numeric vector of length 1 or 2 with values in (0, 1)."
      )
    }
    if (length(significance_level) == 1L) {
      significance_level <- rep(significance_level, 2L)
    }
  }

  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    cli::cli_abort("{.arg parallel} must be TRUE or FALSE.")
  }
  if (parallel) {
    rlang::check_installed(
      "future.apply",
      reason = "for parallel bootstrap execution."
    )
  }

  min_level_n <- tryCatch(
    as.integer(min_level_n),
    warning = function(w) {
      cli::cli_abort("{.arg min_level_n} must be coercible to a non-negative integer.")
    }
  )
  if (length(min_level_n) != 1L || is.na(min_level_n) || min_level_n < 0L) {
    cli::cli_abort("{.arg min_level_n} must be a single non-negative integer.")
  }

  min_events_per_level <- tryCatch(
    as.integer(min_events_per_level),
    warning = function(w) {
      cli::cli_abort("{.arg min_events_per_level} must be coercible to a non-negative integer.")
    }
  )
  if (length(min_events_per_level) != 1L || is.na(min_events_per_level) ||
      min_events_per_level < 0L) {
    cli::cli_abort("{.arg min_events_per_level} must be a single non-negative integer.")
  }

  # ---- Extract settings from the fitted object --------------------------------

  scm_settings <- .extract_scm_settings(scm_fit, significance_level)

  # ---- Pre-bootstrap suitability check (Phase 6) ------------------------------

  boot_diagnostics <- .check_bootstrap_suitability(
    data         = scm_settings$data,
    covariates   = scm_settings$covariates,
    exposure_vars = scm_settings$exposure_vars,
    dv_col       = "DV",
    n_boot       = n_boot,
    subject_col  = subject_col
  )
  class(boot_diagnostics) <- c("bootstrap_diagnostics", class(boot_diagnostics))

  if (boot_diagnostics$risk_level == "high") {
    cli::cli_warn(c(
      "!" = "Factor level risk detected in bootstrap candidates (risk level: HIGH).",
      " " = "Some variables have rare or empty factor levels that may cause convergence failures.",
      "i" = boot_diagnostics$recommendations
    ))
  } else if (boot_diagnostics$risk_level == "moderate") {
    cli::cli_inform(c(
      "i" = "Factor level risk detected in bootstrap candidates (risk level: MODERATE).",
      " " = boot_diagnostics$recommendations
    ))
  }

  # ---- Run bootstrap replicates -----------------------------------------------

  cli::cli_inform(c(
    "i" = "Running {n_boot} bootstrap replicates for a {.val {scm_settings$type}} SCM.",
    " " = "Significance thresholds: {.val {scm_settings$significance_level[1]}} (entry) / {.val {scm_settings$significance_level[2]}} (retention)."
  ))

  raw_results <- .run_bootstrap_replicates(
    data                 = scm_settings$data,
    n_boot               = n_boot,
    subject_col          = subject_col,
    scm_settings         = scm_settings,
    parallel             = parallel,
    seed                 = seed,
    min_level_n          = min_level_n,
    min_events_per_level = min_events_per_level
  )

  # ---- Unpack per-replicate results -------------------------------------------

  boot_converged  <- vapply(raw_results, function(r) isTRUE(r$converged), logical(1L))
  boot_models     <- lapply(raw_results, function(r) r$selected_terms)
  boot_excluded   <- lapply(raw_results, function(r) {
    if (!is.null(r$excluded_covariates)) r$excluded_covariates else character(0L)
  })

  # Per-boot exposure selection (n_boot length, NA for non-converged)
  exposure_selected_per_boot_all <- vapply(boot_models, function(m) {
    if (is.null(m)) NA_character_
    else if (is.list(m) && !is.null(m[["exposure"]])) m[["exposure"]]
    else NA_character_
  }, character(1L))

  n_converged      <- sum(boot_converged)
  convergence_rate <- n_converged / n_boot

  # ---- Convergence health checks ----------------------------------------------

  if (n_converged == 0L) {
    cli::cli_abort(
      c(
        "No bootstrap replicates converged.",
        "x" = "All {n_boot} replicates failed.",
        "i" = "The SCM procedure may not be appropriate for this dataset or configuration. Check data quality, reduce the number of candidate covariates, or increase {.arg significance_level}."
      )
    )
  }

  if (convergence_rate < 0.8) {
    cli::cli_warn(
      "Only {n_converged}/{n_boot} bootstrap replicates converged ({round(convergence_rate * 100)}%). Consider checking data quality or reducing the number of candidate covariates."
    )
  }

  # ---- Convergence issues tibble (Phase 6) ------------------------------------

  failed_idx <- which(!boot_converged)
  convergence_issues <- tibble::tibble(
    replicate = as.integer(failed_idx),
    reason    = rep("error", length(failed_idx))
  )

  # ---- Variable exclusion rate (Phase 6) --------------------------------------
  # For each candidate covariate, compute proportion of all replicates where it
  # was excluded due to sparse/separated factor levels. Only non-zero results kept.
  candidate_covs_for_excl <- scm_settings$covariates

  variable_exclusion_rate <- if (
    (min_level_n > 0L || min_events_per_level > 0L) && n_boot > 0L
  ) {
    excl_counts <- vapply(candidate_covs_for_excl, function(v) {
      sum(vapply(boot_excluded, function(ex) v %in% ex, logical(1L)))
    }, numeric(1L))
    excl_rates <- excl_counts / n_boot
    if (any(excl_rates > 0)) excl_rates else NULL
  } else {
    NULL
  }

  # ---- Aggregation ------------------------------------------------------------

  converged_models <- boot_models[boot_converged]

  # Exposure-selected vector for converged replicates only (length = n_converged)
  exposure_selected_per_boot <- exposure_selected_per_boot_all[boot_converged]

  # Inclusion probabilities for candidate covariates (exclude exposure/adjust_for)
  candidate_covs <- scm_settings$covariates
  inclusion_prob <- .compute_inclusion_prob(converged_models, candidate_covs)

  # Inclusion probabilities for interaction terms
  all_selected_covs <- unique(unlist(lapply(converged_models, .get_covariates)))
  is_inter_cov      <- grepl(":", all_selected_covs, fixed = TRUE)
  all_inter_cov     <- all_selected_covs[is_inter_cov]
  interaction_inclusion_prob <- if (length(all_inter_cov) > 0L) {
    .compute_inclusion_prob(converged_models, sort(all_inter_cov))
  } else {
    stats::setNames(numeric(0L), character(0L))
  }

  # Co-inclusion matrix over all unique covariate terms that appeared
  all_main_cov <- all_selected_covs[!is_inter_cov]
  all_terms_sorted <- sort(unique(c(all_main_cov, all_inter_cov)))
  co_inclusion <- .compute_co_inclusion_matrix(converged_models, all_terms_sorted)

  # Exposure selection probabilities
  exposure_selection_prob <- .compute_exposure_selection_prob(
    converged_models, scm_settings$exposure_vars
  )

  # Conditional covariate inclusion by exposure (NULL when only one candidate)
  inclusion_prob_by_exposure <- .compute_inclusion_prob_by_exposure(
    converged_models,
    candidate_terms  = candidate_covs,
    candidate_exposures = scm_settings$exposure_vars
  )

  # Modal model (returns list(exposure, covariates) for new format)
  modal_result            <- .identify_modal_model(converged_models)
  modal_model             <- modal_result$model
  modal_model_freq        <- modal_result$frequency
  model_variability_index <- 1 - modal_model_freq

  # Build original_selected as list(exposure, covariates) for rich comparison
  original_selected <- list(
    exposure   = scm_settings$selected_exposure,
    covariates = scm_settings$original_selected
  )

  # ---- Transform selection frequencies (Phase 4) ------------------------------
  tf <- .compute_transform_freq(
    converged_models = converged_models,
    nonlinear_spec   = scm_settings$nonlinear_spec,
    exposure_vars    = scm_settings$exposure_vars
  )
  transform_selection_freq  <- tf$covariate_transform_freq
  exposure_transform_freq   <- tf$exposure_transform_freq

  # ---- Assemble output --------------------------------------------------------

  out <- list(
    inclusion_prob             = inclusion_prob,
    interaction_inclusion_prob = interaction_inclusion_prob,
    co_inclusion               = co_inclusion,
    modal_model                = modal_model,
    modal_model_freq           = modal_model_freq,
    model_variability_index    = model_variability_index,
    exposure_selection_prob    = exposure_selection_prob,
    exposure_selected_per_boot = exposure_selected_per_boot,
    inclusion_prob_by_exposure = inclusion_prob_by_exposure,
    boot_models                = boot_models,
    boot_converged             = boot_converged,
    n_boot                     = n_boot,
    n_converged                = n_converged,
    convergence_rate           = convergence_rate,
    convergence_issues         = convergence_issues,
    variable_exclusion_rate    = variable_exclusion_rate,
    diagnostics                = boot_diagnostics,
    seed                       = seed,
    type                       = scm_settings$type,
    exposure_vars              = scm_settings$exposure_vars,
    covariates                 = scm_settings$covariates,
    adjust_for                 = scm_settings$adjust_for,
    significance_level         = scm_settings$significance_level,
    original_selected          = original_selected,
    transform_selection_freq   = transform_selection_freq,
    exposure_transform_freq    = exposure_transform_freq,
    call                       = mc
  )
  class(out) <- "scm_bootstrap"
  out
}
