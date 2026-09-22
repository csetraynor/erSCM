# S3 methods for "scm_bootstrap" objects produced by scm_bootstrap().
# print, summary, plot, stable_covariates.

# ---- Internal helpers --------------------------------------------------------

#' Compute binomial confidence interval
#'
#' @param x Integer, number of successes.
#' @param n Integer, number of trials.
#' @param method Character, \code{"wilson"} or \code{"exact"}.
#' @param conf_level Numeric, confidence level (default 0.95).
#'
#' @return Named numeric vector with elements \code{"lower"} and \code{"upper"}.
#' @keywords internal
.binom_ci <- function(x, n, method = "wilson", conf_level = 0.95) {
  if (n == 0L) return(c(lower = 0, upper = 1))
  p_hat <- x / n
  alpha <- 1 - conf_level
  z     <- stats::qnorm(1 - alpha / 2)

  if (method == "wilson") {
    denom  <- 1 + z^2 / n
    centre <- (p_hat + z^2 / (2 * n)) / denom
    margin <- (z / denom) * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))
    lower  <- max(0, centre - margin)
    upper  <- min(1, centre + margin)
  } else {
    lower <- if (x == 0L) 0 else stats::qbeta(alpha / 2,       x,     n - x + 1)
    upper <- if (x == n)  1 else stats::qbeta(1 - alpha / 2,   x + 1, n - x)
  }
  c(lower = lower, upper = upper)
}

#' Build a progress bar string (Unicode or ASCII fallback)
#' @param p Numeric in 0,1.
#' @param width Integer, total bar width in characters.
#' @keywords internal
.block_bar <- function(p, width = 20L) {
  if (cli::is_utf8_output()) {
    blocks <- c(" ", "\u258f", "\u258e", "\u258d", "\u258c",
                "\u258b", "\u258a", "\u2589", "\u2588")
    filled_eighths <- round(p * width * 8)
    full_blocks    <- filled_eighths %/% 8L
    remainder      <- filled_eighths %% 8L
    bar <- paste0(
      strrep("\u2588", full_blocks),
      if (remainder > 0L && full_blocks < width) blocks[remainder + 1L] else "",
      strrep(" ", max(0L, width - full_blocks - (remainder > 0L)))
    )
  } else {
    filled <- round(p * width)
    bar <- paste0(strrep("#", filled), strrep(".", max(0L, width - filled)))
  }
  paste0("[", bar, "]")
}

# ---- print.scm_bootstrap -----------------------------------------------------

#' Print Bootstrap Model Selection Stability Results
#'
#' Displays a formatted summary of bootstrap stability results including
#' inclusion probabilities, the modal model, and variability assessment.
#'
#' @param x An object of class \code{"scm_bootstrap"}.
#' @param n_top Integer, number of top covariates to display (default 15).
#' @param threshold Numeric, inclusion probability threshold to highlight
#'   "stable" covariates (default 0.7, i.e., selected in \eqn{\ge}70\% of
#'   replicates).
#' @param ... Additional arguments (ignored).
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' data(demo_eres)
#' fit <- scm_glm(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", covariates = c("AGE", "WT"),
#'   dv_col = "sim_status", significance_level = c(0.2, 0.2)
#' )
#' stability <- scm_bootstrap(fit, n_boot = 10, seed = 42)
#' print(stability)
#'
#' @seealso \code{\link{scm_bootstrap}}, \code{\link{summary.scm_bootstrap}},
#'   \code{\link{plot.scm_bootstrap}}
#' @family bootstrap
#' @export
print.scm_bootstrap <- function(x, n_top = 15L, threshold = 0.7, ...) {
  type_label <- if (x$type == "glm") "GLM (logistic)" else "Cox Proportional Hazards"

  # Support both list(exposure, covariates) and legacy character vector
  orig_covs <- if (is.list(x$original_selected)) {
    x$original_selected$covariates
  } else {
    x$original_selected
  }
  orig_exp <- if (is.list(x$original_selected)) x$original_selected$exposure else NA_character_

  cli::cli_h1("Bootstrap Model Selection Stability")

  # --- Metadata block ----------------------------------------------------------
  cli::cli_text("{.strong Type:}            {type_label}")
  cli::cli_text("{.strong Exposure candidates:} {paste(x$exposure_vars, collapse = ', ')}")
  adj <- if (is.null(x$adjust_for) || length(x$adjust_for) == 0L) "None"
         else paste(x$adjust_for, collapse = ", ")
  cli::cli_text("{.strong Forced adjustments:}   {adj}")
  cli::cli_text(
    "{.strong Replicates:}       {x$n_converged}/{x$n_boot} converged ({round(x$convergence_rate * 100, 1)}%)"
  )
  cli::cli_text(
    "{.strong Sig. levels:}      entry = {x$significance_level[1]}, retention = {x$significance_level[2]}"
  )

  # --- Exposure selection stability (only shown when >1 candidate) -------------
  if (!is.null(x$exposure_selection_prob) && length(x$exposure_selection_prob) > 1L) {
    cli::cli_h2("Exposure Selection Stability")
    probs_exp <- sort(x$exposure_selection_prob, decreasing = TRUE)
    for (nm in names(probs_exp)) {
      p   <- probs_exp[[nm]]
      bar <- .block_bar(p, 20L)
      is_orig <- isTRUE(!is.na(orig_exp) && nm == orig_exp)
      orig_mark <- if (is_orig) if (cli::is_utf8_output()) "\u25c6" else "*" else " "
      cat(sprintf("  %-12s  %s  %5.1f%%  %s\n", nm, bar, p * 100, orig_mark))
    }
  }

  # --- Inclusion probabilities -------------------------------------------------
  cli::cli_h2("Covariate Inclusion Probabilities")

  all_main  <- x$inclusion_prob
  all_inter <- x$interaction_inclusion_prob

  if (length(all_main) == 0L && length(all_inter) == 0L) {
    cli::cli_text("  (no candidate covariates screened)")
  } else {
    .print_inclusion_block(all_main, threshold, n_top, orig_covs)

    if (length(all_inter) > 0L && any(all_inter > 0)) {
      cli::cli_h3("Interaction Terms")
      .print_inclusion_block(all_inter, threshold, n_top, orig_covs)
    }
  }

  # --- Modal model -------------------------------------------------------------
  # modal_model is list(exposure, covariates) in new format, or char vector legacy
  modal_covs <- if (is.list(x$modal_model) && !is.null(x$modal_model$covariates)) {
    x$modal_model$covariates
  } else {
    x$modal_model
  }
  modal_exp <- if (is.list(x$modal_model) && !is.null(x$modal_model$exposure)) {
    x$modal_model$exposure
  } else {
    NA_character_
  }

  cli::cli_h2("Modal Model ({round(x$modal_model_freq * 100, 1)}% of converged replicates)")

  if (!is.na(modal_exp) && modal_exp != "(none)") {
    cli::cli_text("  Exposure: {modal_exp}")
  }
  if (length(modal_covs) == 0L) {
    cli::cli_text("  Covariates: (none selected)")
  } else {
    for (tm in sort(modal_covs)) cli::cli_bullets(c("*" = tm))
  }

  orig_sorted  <- sort(orig_covs)
  modal_sorted <- sort(modal_covs)
  orig_exp_match <- is.na(orig_exp) || is.na(modal_exp) || isTRUE(orig_exp == modal_exp)
  if (identical(orig_sorted, modal_sorted) && orig_exp_match) {
    cli::cli_alert_success("Matches the original fit")
  } else {
    cli::cli_alert_info("Differs from the original fit")
    exp_str  <- if (!is.na(orig_exp)) orig_exp else "(unknown)"
    covs_str <- if (length(orig_covs) > 0L) paste(sort(orig_covs), collapse = ", ")
                else "(none)"
    cli::cli_text("  Original: exposure={exp_str}, covariates={covs_str}")
  }

  # --- Non-linear form selection (Phase 4) -------------------------------------
  tf_freq   <- x[["transform_selection_freq"]]
  exp_tf    <- x[["exposure_transform_freq"]]
  has_tf    <- (!is.null(tf_freq) && length(tf_freq) > 0L) ||
               (!is.null(exp_tf)  && length(exp_tf)  > 0L)
  if (has_tf) {
    cli::cli_h2("Non-Linear Form Selection")
    .print_transform_freq_block(exp_tf,  x$exposure_vars, is_exposure = TRUE)
    .print_transform_freq_block(tf_freq, names(tf_freq),  is_exposure = FALSE)
  }

  # --- Factor Level Warnings (Phase 6) -----------------------------------------
  excl_rate <- x[["variable_exclusion_rate"]]
  if (!is.null(excl_rate) && length(excl_rate) > 0L) {
    high_excl <- excl_rate[excl_rate > 0.1]
    if (length(high_excl) > 0L) {
      cli::cli_h2("Factor Level Warnings")
      cli::cli_text(
        "  The following variable{?s} {?was/were} excluded from >10% of replicates due to sparse/separated factor levels:"
      )
      for (nm in names(sort(high_excl, decreasing = TRUE))) {
        cat(sprintf("  %-20s  excluded in %5.1f%% of replicates\n",
                    nm, high_excl[[nm]] * 100))
      }
      cli::cli_alert_warning(
        "Inclusion probabilities for these variables may be underestimated."
      )
    }
  }

  # --- Stability assessment ----------------------------------------------------
  cli::cli_h2("Stability Assessment")

  vi  <- x$model_variability_index
  vi_label <- if (vi < 0.3) {
    cli::cli_alert_success("Variability index: {round(vi, 3)} -- Low variability, stable selection")
  } else if (vi < 0.6) {
    cli::cli_alert_warning("Variability index: {round(vi, 3)} -- Moderate variability")
  } else {
    cli::cli_alert_danger("Variability index: {round(vi, 3)} -- High variability, unstable selection")
  }

  probs <- x$inclusion_prob
  n_stable   <- sum(probs >= threshold)
  n_moderate <- sum(probs >= 0.3 & probs < threshold)
  n_unstable <- sum(probs < 0.3)
  cli::cli_text(
    "  Stable (>={threshold}): {n_stable} | Moderate (0.3-{threshold}): {n_moderate} | Unstable (<0.3): {n_unstable}"
  )

  if (x$convergence_rate < 1.0) {
    n_failed <- x$n_boot - x$n_converged
    cli::cli_alert_warning(
      "{n_failed} replicate{?s} did not converge and {?was/were} excluded from all calculations."
    )
  }

  invisible(x)
}

#' Internal: print a block of inclusion probabilities with bar chart
#' @keywords internal
.print_inclusion_block <- function(probs, threshold, n_top, original_selected) {
  if (length(probs) == 0L) return(invisible(NULL))
  probs_sorted <- sort(probs, decreasing = TRUE)
  probs_show   <- utils::head(probs_sorted, n_top)

  max_name <- max(nchar(names(probs_show)), 12L)
  fmt      <- paste0("  %-", max_name, "s  %s  %5.1f%%  %s")

  for (nm in names(probs_show)) {
    p   <- probs_show[[nm]]
    bar <- .block_bar(p, 20L)
    marker <- if (p >= threshold) "***"
              else if (p >= 0.3)  " . "
              else                "   "
    orig_mark <- if (nm %in% original_selected) if (cli::is_utf8_output()) "\u25c6" else "*" else " "
    cat(sprintf(fmt, nm, bar, p * 100, paste0(marker, " ", orig_mark)), "\n")
  }

  if (length(probs_sorted) > n_top) {
    cli::cli_text("  ... and {length(probs_sorted) - n_top} more (use n_top to show all)")
  }
}

# ---- summary.scm_bootstrap ---------------------------------------------------

#' Summarise Bootstrap Model Selection Stability
#'
#' Creates a structured summary of bootstrap stability results, classifying
#' covariates by their selection stability.
#'
#' @param object An object of class \code{"scm_bootstrap"}.
#' @param threshold Numeric, stability threshold (default 0.7). Covariates
#'   selected in \eqn{\ge} this proportion of replicates are classified as
#'   \code{"stable"}.
#' @param ... Additional arguments (ignored).
#'
#' @return A list of class \code{"summary_scm_bootstrap"} (returned invisibly)
#'   with components:
#'   \describe{
#'     \item{inclusion_table}{Tibble with columns: \code{term},
#'       \code{inclusion_prob}, \code{stability_class}, \code{term_type},
#'       \code{rank}, \code{in_original}.}
#'     \item{stable_covariates}{Character vector of stable terms.}
#'     \item{moderate_covariates}{Character vector of moderate terms.}
#'     \item{unstable_covariates}{Character vector of unstable terms.}
#'     \item{co_inclusion}{The co-inclusion matrix.}
#'     \item{modal_model}{Character vector of the modal model.}
#'     \item{modal_model_freq}{Numeric.}
#'     \item{variability_index}{Numeric.}
#'     \item{original_selected}{Character vector of covariates from the original fit.}
#'     \item{original_exposure}{Character; the exposure selected in the original fit.}
#'     \item{exposure_stability}{Tibble of exposure selection probabilities
#'       (when multiple exposure candidates were screened).}
#'     \item{original_match_rate}{Proportion of converged replicates that
#'       selected exactly the same model as the original fit.}
#'     \item{threshold}{The threshold used.}
#'     \item{n_converged}{Integer.}
#'     \item{type}{Character.}
#'   }
#'
#' @examples
#' data(demo_eres)
#' fit <- scm_glm(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", covariates = c("AGE", "WT"),
#'   dv_col = "sim_status", significance_level = c(0.2, 0.2)
#' )
#' stability <- scm_bootstrap(fit, n_boot = 10, seed = 42)
#' summary(stability)
#'
#' @seealso \code{\link{scm_bootstrap}}, \code{\link{print.scm_bootstrap}},
#'   \code{\link{plot.scm_bootstrap}}
#' @family bootstrap
#' @export
summary.scm_bootstrap <- function(object, threshold = 0.7, ...) {
  x <- object

  # Support both list(exposure, covariates) and legacy character vector
  orig_covs <- if (is.list(x$original_selected)) {
    x$original_selected$covariates
  } else {
    x$original_selected
  }
  orig_exp <- if (is.list(x$original_selected)) x$original_selected$exposure else NA_character_

  # ---- Build combined inclusion tibble ----------------------------------------
  make_rows <- function(probs, type_label) {
    if (length(probs) == 0L) return(NULL)
    tibble::tibble(
      term           = names(probs),
      inclusion_prob = as.numeric(probs),
      term_type      = type_label,
      in_original    = names(probs) %in% orig_covs
    )
  }

  tbl <- dplyr::bind_rows(
    make_rows(x$inclusion_prob,             "main_effect"),
    make_rows(x$interaction_inclusion_prob, "interaction")
  )

  tbl <- tbl %>%
    dplyr::mutate(
      stability_class = dplyr::case_when(
        inclusion_prob >= threshold ~ "stable",
        inclusion_prob >= 0.3       ~ "moderate",
        TRUE                        ~ "unstable"
      )
    ) %>%
    dplyr::arrange(dplyr::desc(inclusion_prob)) %>%
    dplyr::mutate(rank = dplyr::row_number())

  # ---- Exposure stability tibble ----------------------------------------------
  exposure_stability <- if (!is.null(x$exposure_selection_prob) &&
                             length(x$exposure_selection_prob) > 0L) {
    tibble::tibble(
      exposure_var  = names(x$exposure_selection_prob),
      selection_prob = as.numeric(x$exposure_selection_prob),
      is_original   = names(x$exposure_selection_prob) == orig_exp
    ) %>%
      dplyr::arrange(dplyr::desc(selection_prob))
  } else {
    NULL
  }

  # ---- Original match rate (full model: exposure + covariates) ----------------
  converged_models <- x$boot_models[x$boot_converged]
  orig_covs_sorted <- sort(orig_covs)
  original_match_rate <- if (length(converged_models) == 0L) {
    NA_real_
  } else {
    mean(vapply(converged_models, function(m) {
      cov_match <- identical(sort(.get_covariates(m)), orig_covs_sorted)
      # Also compare exposure if stored
      exp_match <- if (is.list(m) && !is.null(m[["exposure"]]) && !is.na(orig_exp)) {
        isTRUE(m[["exposure"]] == orig_exp)
      } else {
        TRUE  # legacy: no exposure stored in m
      }
      cov_match && exp_match
    }, logical(1L)))
  }

  # ---- Classify terms ---------------------------------------------------------
  stable_terms   <- tbl$term[tbl$stability_class == "stable"]
  moderate_terms <- tbl$term[tbl$stability_class == "moderate"]
  unstable_terms <- tbl$term[tbl$stability_class == "unstable"]

  out <- list(
    inclusion_table      = tbl,
    stable_covariates    = stable_terms,
    moderate_covariates  = moderate_terms,
    unstable_covariates  = unstable_terms,
    co_inclusion         = x$co_inclusion,
    modal_model          = x$modal_model,
    modal_model_freq     = x$modal_model_freq,
    variability_index    = x$model_variability_index,
    original_selected    = orig_covs,
    original_exposure    = orig_exp,
    exposure_stability   = exposure_stability,
    original_match_rate  = original_match_rate,
    threshold            = threshold,
    n_converged          = x$n_converged,
    type                 = x$type
  )
  class(out) <- "summary_scm_bootstrap"
  invisible(out)
}

#' Print method for summary_scm_bootstrap objects
#'
#' @param x An object of class \code{"summary_scm_bootstrap"}.
#' @param ... Additional arguments (ignored).
#'
#' @return \code{x}, invisibly.
#' @examples
#' data(demo_eres)
#' fit <- scm_glm(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", covariates = c("AGE", "WT"),
#'   dv_col = "sim_status", significance_level = c(0.2, 0.2)
#' )
#' stability <- scm_bootstrap(fit, n_boot = 10, seed = 42)
#' print(summary(stability))
#' @export
print.summary_scm_bootstrap <- function(x, ...) {
  cli::cli_h1("Bootstrap SCM Stability Summary")
  cli::cli_text("Type: {x$type} | Converged replicates: {x$n_converged} | Threshold: {x$threshold}")

  # ---- Inclusion table --------------------------------------------------------
  cli::cli_h2("Inclusion Table")
  check_mark <- if (cli::is_utf8_output()) "+" else "x"
  tbl_print <- x$inclusion_table %>%
    dplyr::mutate(
      inclusion_prob = round(inclusion_prob, 3),
      in_original    = ifelse(in_original, check_mark, "")
    ) %>%
    dplyr::select(rank, term, inclusion_prob, stability_class, term_type, in_original)
  print(as.data.frame(tbl_print), row.names = FALSE)

  # ---- Stability classification counts ----------------------------------------
  cli::cli_h2("Stability Classification (threshold = {x$threshold})")
  cli::cli_bullets(c(
    "*" = "Stable   (>={x$threshold}): {length(x$stable_covariates)} term{?s} -- {paste(x$stable_covariates, collapse = ', ')}",
    "*" = "Moderate (0.3-{x$threshold}): {length(x$moderate_covariates)} term{?s}",
    "*" = "Unstable (<0.3): {length(x$unstable_covariates)} term{?s}"
  ))

  # ---- Exposure stability -----------------------------------------------------
  if (!is.null(x$exposure_stability) && nrow(x$exposure_stability) > 1L) {
    cli::cli_h2("Exposure Selection Stability")
    print(as.data.frame(dplyr::mutate(
      x$exposure_stability,
      selection_prob = round(selection_prob, 3),
      is_original    = ifelse(is_original, check_mark, "")
    )), row.names = FALSE)
  }

  # ---- Original model match rate ----------------------------------------------
  cli::cli_h2("Original Model Agreement")
  if (!is.na(x$original_match_rate)) {
    cli::cli_text(
      "{round(x$original_match_rate * 100, 1)}% of replicates selected exactly the same model as the original fit."
    )
  }
  if (!is.na(x$original_exposure)) {
    cli::cli_text("Original exposure: {x$original_exposure}")
  }
  if (length(x$original_selected) > 0L) {
    cli::cli_text("Original covariates: {paste(sort(x$original_selected), collapse = ', ')}")
  }

  # ---- Co-inclusion highlights ------------------------------------------------
  mat <- x$co_inclusion
  if (!is.null(mat) && nrow(mat) >= 2L) {
    cli::cli_h2("Co-inclusion Highlights")

    pairs <- utils::combn(rownames(mat), 2L, simplify = FALSE)
    high_coinc    <- character(0L)
    competing     <- character(0L)

    for (p in pairs) {
      a <- p[1L]; b <- p[2L]
      co_val <- mat[a, b]
      marg_a <- if (a %in% names(x$inclusion_table$term)) {
        x$inclusion_table$inclusion_prob[x$inclusion_table$term == a]
      } else 0
      marg_b <- if (b %in% names(x$inclusion_table$term)) {
        x$inclusion_table$inclusion_prob[x$inclusion_table$term == b]
      } else 0

      if (co_val > 0.8) {
        high_coinc <- c(high_coinc,
          sprintf("%s + %s  (co-incl = %.2f)", a, b, co_val))
      } else if (length(marg_a) > 0L && length(marg_b) > 0L &&
                 marg_a[1L] >= 0.5 && marg_b[1L] >= 0.5 && co_val < 0.3) {
        competing <- c(competing,
          sprintf("%s vs %s  (co-incl = %.2f, marg = %.2f / %.2f)",
                  a, b, co_val, marg_a[1L], marg_b[1L]))
      }
    }

    if (length(high_coinc) > 0L) {
      cli::cli_text("{.strong Always selected together (co-incl > 0.8):}")
      for (s in high_coinc) cli::cli_bullets(c("*" = s))
    }
    if (length(competing) > 0L) {
      cli::cli_text("{.strong Competing covariates (high marginal, low co-incl):}")
      for (s in competing) cli::cli_bullets(c("!" = s))
    }
    if (length(high_coinc) == 0L && length(competing) == 0L) {
      cli::cli_text("  No notable co-inclusion patterns detected.")
    }
  }

  invisible(x)
}

# ---- stable_covariates -------------------------------------------------------

#' Extract Stable Covariates from Bootstrap Results
#'
#' Returns covariate names whose bootstrap inclusion probability meets or
#' exceeds the specified threshold, indicating stable selection across
#' resampled datasets.
#'
#' @param x An object of class \code{"scm_bootstrap"}.
#' @param threshold Numeric, minimum inclusion probability to be considered
#'   stable (default 0.7).
#' @param include_interactions Logical, whether to include interaction terms
#'   in the output (default \code{TRUE}).
#' @param include_exposure Logical, whether to include the selected exposure
#'   variable when its selection probability meets the threshold (default
#'   \code{FALSE}). Only meaningful when multiple exposure candidates were
#'   screened; when a single exposure is provided its selection probability
#'   is trivially 1.0 and including it adds no diagnostic value.
#'
#' @return Character vector of term names meeting the threshold, sorted by
#'   inclusion probability (highest first). Returns \code{character(0)} if
#'   no terms meet the threshold.
#'
#' @examples
#' data(demo_eres)
#' fit <- scm_glm(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", covariates = c("AGE", "WT"),
#'   dv_col = "sim_status", significance_level = c(0.2, 0.2)
#' )
#' stability <- scm_bootstrap(fit, n_boot = 10, seed = 42)
#' stable_covariates(stability, threshold = 0.5)
#'
#' @seealso \code{\link{scm_bootstrap}}, \code{\link{summary.scm_bootstrap}}
#' @family bootstrap
#' @export
stable_covariates <- function(x, threshold = 0.7, include_interactions = TRUE,
                              include_exposure = FALSE) {
  if (!inherits(x, "scm_bootstrap")) {
    cli::cli_abort("{.arg x} must be an object of class {.cls scm_bootstrap}.")
  }
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      is.na(threshold) || threshold < 0 || threshold > 1) {
    cli::cli_abort(
      "{.arg threshold} must be a single numeric value in [0, 1]."
    )
  }

  probs <- x$inclusion_prob
  if (include_interactions && length(x$interaction_inclusion_prob) > 0L) {
    probs <- c(probs, x$interaction_inclusion_prob)
  }
  if (include_exposure && !is.null(x$exposure_selection_prob) &&
      length(x$exposure_selection_prob) > 0L) {
    probs <- c(probs, x$exposure_selection_prob)
  }

  stable <- probs[probs >= threshold]
  names(sort(stable, decreasing = TRUE))
}

# ---- plot.scm_bootstrap ------------------------------------------------------

#' Plot Bootstrap Model Selection Stability
#'
#' Creates visualisations of covariate selection stability from bootstrap
#' SCM results.
#'
#' @param x An object of class \code{"scm_bootstrap"}.
#' @param type Character, one of \code{"inclusion"} (default),
#'   \code{"coinclusion"}, \code{"both"}, \code{"exposure"}, or
#'   \code{"transforms"} (stacked barplot of transformation selection
#'   frequencies; only meaningful when the original fit used a nonlinear
#'   specification).
#' @param threshold Numeric, stability threshold line on the inclusion plot
#'   (default 0.7).
#' @param ci_method Character, method for confidence intervals on inclusion
#'   proportions. One of \code{"wilson"} (default) or \code{"exact"}
#'   (Clopper-Pearson).
#' @param n_top Integer, maximum number of terms to display (default 20).
#' @param show_original Logical, whether to mark covariates selected in the
#'   original fit (default \code{TRUE}).
#' @param ... Additional arguments (ignored).
#'
#' @return A \code{ggplot2} object (or a \code{patchwork} combined plot when
#'   \code{type = "both"}). Returned invisibly.
#'
#' @examples
#' data(demo_eres)
#' fit <- scm_glm(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", covariates = c("AGE", "WT"),
#'   dv_col = "sim_status", significance_level = c(0.2, 0.2)
#' )
#' stability <- scm_bootstrap(fit, n_boot = 10, seed = 42)
#' plot(stability)
#'
#' @seealso \code{\link{scm_bootstrap}}, \code{\link{print.scm_bootstrap}},
#'   \code{\link{summary.scm_bootstrap}}
#' @family bootstrap
#' @export
plot.scm_bootstrap <- function(x,
                                type         = c("inclusion", "coinclusion", "both",
                                                 "exposure", "transforms"),
                                threshold    = 0.7,
                                ci_method    = c("wilson", "exact"),
                                n_top        = 20L,
                                show_original = TRUE,
                                ...) {
  type      <- match.arg(type)
  ci_method <- match.arg(ci_method)

  p_incl  <- NULL
  p_coinc <- NULL

  if (type %in% c("inclusion", "both")) {
    p_incl <- .plot_inclusion(x, threshold, ci_method, n_top, show_original)
  }
  if (type %in% c("coinclusion", "both")) {
    p_coinc <- .plot_coinclusion(x)
  }
  if (type == "exposure") {
    return(invisible(.plot_exposure_stability(x, ci_method, show_original)))
  }
  if (type == "transforms") {
    return(invisible(.plot_transform_stability(x)))
  }

  if (type == "both") {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      p <- patchwork::wrap_plots(p_incl, p_coinc, ncol = 1L, heights = c(2, 3))
    } else {
      cli::cli_inform(
        "Install {.pkg patchwork} for combined plots. Showing inclusion plot only."
      )
      p <- p_incl
    }
  } else if (type == "inclusion") {
    p <- p_incl
  } else {
    p <- p_coinc
  }

  invisible(p)
}

# ---- Internal plot builders --------------------------------------------------

#' Build the inclusion probability bar plot
#' @keywords internal
.plot_inclusion <- function(x, threshold, ci_method, n_top, show_original) {

  # Combine main + interaction inclusion probs
  # Guard: NULL or zero-length means none tested
  all_main  <- x$inclusion_prob
  all_inter <- x$interaction_inclusion_prob
  .tag_probs <- function(probs, tag) {
    if (is.null(probs) || length(probs) == 0L) return(stats::setNames(numeric(0L), character(0L)))
    stats::setNames(as.numeric(probs), paste0(names(probs), "|", tag))
  }
  probs_all <- c(.tag_probs(all_main, "main_effect"), .tag_probs(all_inter, "interaction"))

  # Keep top n_top by probability
  probs_top <- utils::head(sort(probs_all, decreasing = TRUE), n_top)

  # Build plotting data.frame
  n_conv <- x$n_converged
  df <- data.frame(
    term_raw       = names(probs_top),
    inclusion_prob = as.numeric(probs_top),
    stringsAsFactors = FALSE
  )
  df$term      <- sub("\\|.*$", "", df$term_raw)
  df$term_type <- sub("^[^|]+\\|", "", df$term_raw)

  df$stability_class <- ifelse(
    df$inclusion_prob >= threshold, "stable",
    ifelse(df$inclusion_prob >= 0.3, "moderate", "unstable")
  )
  orig_covs_plot <- if (is.list(x$original_selected)) {
    x$original_selected$covariates
  } else {
    x$original_selected
  }
  df$in_original <- df$term %in% orig_covs_plot

  # Confidence intervals
  if (nrow(df) == 0L) {
    df$ci_lower <- numeric(0L)
    df$ci_upper <- numeric(0L)
  } else {
    successes <- round(df$inclusion_prob * n_conv)
    ci_mat <- mapply(.binom_ci, x = successes, n = n_conv,
                     MoreArgs = list(method = ci_method),
                     SIMPLIFY = TRUE)
    df$ci_lower <- ci_mat["lower", ]
    df$ci_upper <- ci_mat["upper", ]
  }

  # Factor for ordering (highest prob at top of horizontal bar)
  df$term_fct <- factor(df$term,
                        levels = rev(df$term[base::order(df$inclusion_prob)]))

  colour_map <- c(
    stable   = "#2E8B57",
    moderate = "#DAA520",
    unstable = "#CD5C5C"
  )

  type_label <- if (x$type == "glm") "GLM" else "CoxPH"

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x    = .data$inclusion_prob,
      y    = .data$term_fct,
      fill = .data$stability_class
    )
  ) +
    ggplot2::geom_col(width = 0.65, alpha = 0.9) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$ci_lower, xmax = .data$ci_upper),
      orientation = "y", width = 0.25, linewidth = 0.45, colour = "grey30"
    ) +
    ggplot2::geom_vline(
      xintercept = threshold,
      linetype   = "dashed",
      colour     = "grey40",
      linewidth  = 0.5
    ) +
    ggplot2::annotate(
      "text",
      x     = threshold + 0.02,
      y     = 0.6,
      label = paste0("threshold = ", threshold),
      hjust = 0, vjust = 0, size = 2.8, colour = "grey40"
    ) +
    ggplot2::scale_fill_manual(
      values = colour_map,
      name   = "Stability",
      limits = c("stable", "moderate", "unstable"),
      labels = c(
        stable   = paste0(">= ", threshold, " (stable)"),
        moderate = "0.3 - < threshold",
        unstable = "< 0.3 (unstable)"
      )
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      expand = ggplot2::expansion(mult = c(0, 0.06))
    ) +
    ggplot2::labs(
      title    = paste0("Covariate Selection Stability (B = ", x$n_converged, ")"),
      subtitle = paste0("Type: ", type_label,
                        " | Exposure: ", paste(x$exposure_vars, collapse = ", "),
                        " (protected)"),
      x        = "Inclusion Probability",
      y        = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y  = ggplot2::element_blank(),
      panel.grid.minor    = ggplot2::element_blank(),
      panel.grid.major.x  = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      axis.text.y         = ggplot2::element_text(size = 9),
      legend.position     = "bottom",
      legend.title        = ggplot2::element_text(size = 9),
      legend.text         = ggplot2::element_text(size = 8),
      plot.title          = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle       = ggplot2::element_text(size = 9, colour = "grey40")
    )

  # Mark original-fit selections
  if (show_original && any(df$in_original)) {
    df_orig <- df[df$in_original, ]
    p <- p + ggplot2::geom_point(
      data  = df_orig,
      ggplot2::aes(x = .data$inclusion_prob, y = .data$term_fct),
      shape = 18, size = 3.5, colour = "black",
      inherit.aes = FALSE,
      show.legend = FALSE
    )
  }

  # Separate main effects and interactions when both present
  if (length(unique(df$term_type)) > 1L) {
    p <- p + ggplot2::facet_wrap(
      ~ term_type,
      scales = "free_y",
      ncol   = 1L,
      labeller = ggplot2::labeller(
        term_type = c(main_effect = "Main Effects", interaction = "Interactions")
      )
    )
  }

  p
}

#' Build the exposure selection stability bar chart
#' @keywords internal
.plot_exposure_stability <- function(x, ci_method = "wilson", show_original = TRUE) {
  probs <- x$exposure_selection_prob
  if (is.null(probs) || length(probs) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text", x = 0.5, y = 0.5,
          label = "No exposure selection data available",
          size = 5, colour = "grey50"
        ) +
        ggplot2::theme_void() +
        ggplot2::labs(title = "Exposure Selection Stability")
    )
  }

  orig_exp <- if (is.list(x$original_selected)) x$original_selected$exposure else NA_character_

  n_conv <- x$n_converged
  df <- data.frame(
    exposure_var   = names(probs),
    selection_prob = as.numeric(probs),
    stringsAsFactors = FALSE
  )
  df$in_original <- !is.na(orig_exp) & df$exposure_var == orig_exp

  successes <- round(df$selection_prob * n_conv)
  ci_mat <- mapply(.binom_ci, x = successes, n = n_conv,
                   MoreArgs = list(method = ci_method),
                   SIMPLIFY = TRUE)
  df$ci_lower <- ci_mat["lower", ]
  df$ci_upper <- ci_mat["upper", ]

  df$exposure_fct <- factor(df$exposure_var,
                            levels = rev(df$exposure_var[base::order(df$selection_prob)]))

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$selection_prob,
      y = .data$exposure_fct
    )
  ) +
    ggplot2::geom_col(width = 0.55, alpha = 0.85, fill = "#4682B4") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$ci_lower, xmax = .data$ci_upper),
      orientation = "y", width = 0.2, linewidth = 0.45, colour = "grey30"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      expand = ggplot2::expansion(mult = c(0, 0.06))
    ) +
    ggplot2::labs(
      title    = paste0("Exposure Selection Stability (B = ", n_conv, ")"),
      subtitle = "Proportion of bootstrap replicates each exposure was selected",
      x        = "Selection Probability",
      y        = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.text.y        = ggplot2::element_text(size = 9),
      plot.title         = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle      = ggplot2::element_text(size = 9, colour = "grey40")
    )

  if (show_original && any(df$in_original)) {
    df_orig <- df[df$in_original, ]
    p <- p + ggplot2::geom_point(
      data = df_orig,
      ggplot2::aes(x = .data$selection_prob, y = .data$exposure_fct),
      shape = 18, size = 3.5, colour = "black",
      inherit.aes = FALSE, show.legend = FALSE
    )
  }

  p
}

#' Build the co-inclusion heatmap
#' @keywords internal
.plot_coinclusion <- function(x) {
  mat <- x$co_inclusion

  # Restrict to terms with marginal inclusion > 0.05
  marginal <- diag(mat)
  keep     <- names(marginal)[marginal > 0.05]

  if (length(keep) < 2L) {
    # Not enough terms -- return an informative placeholder
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text", x = 0.5, y = 0.5,
          label = "Fewer than 2 terms with\nmarginal inclusion > 5%",
          size = 5, colour = "grey50"
        ) +
        ggplot2::theme_void() +
        ggplot2::labs(title = "Covariate Co-inclusion Patterns")
    )
  }

  mat_sub <- mat[keep, keep, drop = FALSE]

  # Hierarchical clustering for ordering
  dist_mat <- stats::as.dist(1 - mat_sub)
  hc       <- stats::hclust(dist_mat, method = "complete")
  term_ord <- hc$labels[hc$order]

  # Long format (lower triangle only -- including diagonal for marginal probs)
  pairs <- expand.grid(row_var = term_ord, col_var = term_ord,
                       stringsAsFactors = FALSE)
  pairs$co_inclusion_value <- mapply(
    function(r, c) mat_sub[r, c],
    pairs$row_var, pairs$col_var
  )
  pairs$is_diagonal <- pairs$row_var == pairs$col_var
  pairs$row_idx <- match(pairs$row_var, term_ord)
  pairs$col_idx <- match(pairs$col_var, term_ord)
  # Keep lower triangle + diagonal
  df <- pairs[pairs$row_idx >= pairs$col_idx, ]

  df$label_text <- ifelse(
    df$co_inclusion_value >= 0.1,
    sprintf("%.2f", df$co_inclusion_value),
    ""
  )
  df$row_var <- factor(df$row_var, levels = rev(term_ord))
  df$col_var <- factor(df$col_var, levels = term_ord)

  # Diagonal gets a distinct fill (slightly offset to mark marginals)
  df$fill_val <- df$co_inclusion_value
  df$tile_col <- ifelse(df$is_diagonal, "marginal", "joint")

  ggplot2::ggplot(df, ggplot2::aes(
    x    = .data$col_var,
    y    = .data$row_var,
    fill = .data$fill_val
  )) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$label_text),
      size = 2.8, colour = "grey20"
    ) +
    ggplot2::scale_fill_gradient(
      low   = "white",
      high  = "#2C3E50",
      name  = "Co-inclusion\nProportion",
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25)
    ) +
    ggplot2::scale_x_discrete(position = "bottom") +
    ggplot2::labs(
      title = "Covariate Co-inclusion Patterns",
      x     = NULL,
      y     = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x         = ggplot2::element_text(angle = 35, hjust = 1, size = 8),
      axis.text.y         = ggplot2::element_text(size = 8),
      panel.grid          = ggplot2::element_blank(),
      legend.position     = "right",
      legend.title        = ggplot2::element_text(size = 9),
      plot.title          = ggplot2::element_text(face = "bold", size = 11)
    ) +
    ggplot2::coord_fixed()
}

# ---- Non-linear form helpers (Phase 4) ---------------------------------------

#' Print a transform-frequency block to the console
#' @param freq_list Named list of named numeric vectors (covariate freqs), or a
#'   single named numeric vector (exposure freq). NULL -> silent no-op.
#' @param var_labels Character vector of variable labels to prefix each line.
#' @param is_exposure Logical. If TRUE, prints as a single line for the
#'   exposure; if FALSE, iterates over the list.
#' @keywords internal
.print_transform_freq_block <- function(freq_list, var_labels, is_exposure = FALSE) {
  if (is.null(freq_list) || length(freq_list) == 0L) return(invisible(NULL))

  fmt_pcts <- function(freq_vec) {
    # Sort descending, skip 0-frequency entries, format as "type (XX%)"
    freq_vec <- sort(freq_vec, decreasing = TRUE)
    freq_vec <- freq_vec[freq_vec > 0]
    paste(
      vapply(names(freq_vec), function(nm) {
        sprintf("%s (%s%%)", nm, round(freq_vec[[nm]] * 100))
      }, character(1L)),
      collapse = ", "
    )
  }

  if (is_exposure) {
    if (is.list(freq_list)) {
      # Named list: one entry per NL-enabled exposure variable
      for (ev in names(freq_list)) {
        cat(sprintf("  %-12s  %s\n", ev, fmt_pcts(freq_list[[ev]])))
      }
    } else {
      lbl <- if (length(var_labels) > 0L) var_labels[[1L]] else "Exposure"
      cat(sprintf("  %-12s  %s\n", lbl, fmt_pcts(freq_list)))
    }
  } else {
    for (v in names(freq_list)) {
      cat(sprintf("  %-12s  %s\n", v, fmt_pcts(freq_list[[v]])))
    }
  }
  invisible(NULL)
}

#' Build a stacked horizontal barplot of transform selection stability
#'
#' Shows, for each variable in the nonlinear specification, how often each
#' transformation type (log, fp1, fp2, linear, not_selected, ...) was chosen
#' across bootstrap replicates.
#'
#' @param x An object of class \code{"scm_bootstrap"}.
#' @return A \code{ggplot2} object, or an informative blank plot when no
#'   transform frequency data is available.
#' @keywords internal
.plot_transform_stability <- function(x) {
  tf_freq <- x[["transform_selection_freq"]]
  exp_tf  <- x[["exposure_transform_freq"]]

  has_cov <- !is.null(tf_freq) && length(tf_freq) > 0L
  has_exp <- !is.null(exp_tf)  && length(exp_tf)  > 0L

  if (!has_cov && !has_exp) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text", x = 0.5, y = 0.5,
          label = "No non-linear specification in the original fit.\nUse nonlinear = list(...) in scm_glm() / scm_coxph().",
          size = 4.5, colour = "grey50", hjust = 0.5, vjust = 0.5
        ) +
        ggplot2::theme_void() +
        ggplot2::labs(title = "Transformation Selection Stability")
    )
  }

  # Build a long-format data frame - one row per (variable, transform)
  rows <- list()

  if (has_exp) {
    if (is.list(exp_tf)) {
      # Named list: one named numeric vector per NL-enabled exposure variable
      for (ev in names(exp_tf)) {
        freq_v <- exp_tf[[ev]]
        for (tt in names(freq_v)) {
          rows[[length(rows) + 1L]] <- data.frame(
            variable  = paste0(ev, " [exposure]"),
            transform = tt,
            prop      = as.numeric(freq_v[[tt]]),
            stringsAsFactors = FALSE
          )
        }
      }
    } else {
      exp_lbl <- if (length(x$exposure_vars) > 0L) x$exposure_vars[[1L]] else "Exposure"
      for (tt in names(exp_tf)) {
        rows[[length(rows) + 1L]] <- data.frame(
          variable  = paste0(exp_lbl, " [exposure]"),
          transform = tt,
          prop      = as.numeric(exp_tf[[tt]]),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (has_cov) {
    for (v in names(tf_freq)) {
      freq_v <- tf_freq[[v]]
      for (tt in names(freq_v)) {
        rows[[length(rows) + 1L]] <- data.frame(
          variable  = v,
          transform = tt,
          prop      = as.numeric(freq_v[[tt]]),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  df <- do.call(rbind, rows)

  # Canonical transform colour palette (robust to unknown types)
  base_palette <- c(
    "linear"       = "#AAAAAA",
    "log"          = "#4393C3",
    "sqrt"         = "#74ADD1",
    "fp1"          = "#F4A582",
    "fp2"          = "#D6604D",
    "rcs"          = "#4DAC26",
    "not_selected" = "#EEEEEE"
  )
  all_types     <- unique(df$transform)
  extra_types   <- setdiff(all_types, names(base_palette))
  extra_cols    <- grDevices::hcl.colors(max(1L, length(extra_types)), "Set2")
  full_palette  <- c(base_palette,
                     stats::setNames(extra_cols, extra_types))
  used_palette  <- full_palette[intersect(names(full_palette), all_types)]

  # Variable order: exposure first, then alphabetical covariates
  var_order <- rev(unique(df$variable))
  df$variable_fct <- factor(df$variable, levels = var_order)

  # Transform order within bars: not_selected last, else alphabetical
  tf_order <- c(
    sort(setdiff(unique(df$transform), "not_selected")),
    if ("not_selected" %in% df$transform) "not_selected" else character(0L)
  )
  df$transform_fct <- factor(df$transform, levels = tf_order)

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x    = .data$prop,
      y    = .data$variable_fct,
      fill = .data$transform_fct
    )
  ) +
    ggplot2::geom_col(width = 0.7, position = "stack") +
    ggplot2::geom_text(
      ggplot2::aes(
        label = ifelse(.data$prop >= 0.05,
                       paste0(round(.data$prop * 100), "%"), "")
      ),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 3, colour = "grey20", fontface = "plain"
    ) +
    ggplot2::scale_fill_manual(
      values = used_palette,
      name   = "Transform",
      breaks = tf_order
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      labels = scales::percent_format(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      title    = paste0("Transformation Selection Stability (B = ", x$n_converged, ")"),
      subtitle = "Proportion of converged bootstrap replicates selecting each form",
      x        = "Selection Proportion",
      y        = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y  = ggplot2::element_blank(),
      panel.grid.minor    = ggplot2::element_blank(),
      panel.grid.major.x  = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      axis.text.y         = ggplot2::element_text(size = 9),
      legend.position     = "right",
      legend.title        = ggplot2::element_text(size = 9),
      legend.text         = ggplot2::element_text(size = 8),
      plot.title          = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle       = ggplot2::element_text(size = 9, colour = "grey40")
    )
}
