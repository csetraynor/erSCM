# Nonlinear-aware backward elimination step and summary enhancement helpers.

# ---------------------------------------------------------------------------
# Backward deletion with simplification + closed testing
# ---------------------------------------------------------------------------

#' Backward Deletion Step with Non-Linear Terms
#'
#' Enhanced backward deletion that handles variables in non-linear forms.
#' A variable is removed as a WHOLE (all its terms together), not term-by-term.
#' Optionally re-tests whether a simpler form is adequate before full removal
#' (closed-testing: FP2 -> linear -> remove).
#'
#' @param eres_data data.frame.
#' @param current_model List with current model state (\code{Relation},
#'   \code{toto_tot}, \code{finalforward}, plus optionally
#'   \code{transform_registry}).
#' @param model_type Character, \code{"glm"} or \code{"coxph"}.
#' @param family Family for GLM.
#' @param exit_alpha Numeric. Retention threshold -- p > exit_alpha -> removal.
#' @param exit_alpha_interaction Numeric or NULL. Separate threshold for
#'   interaction terms.
#' @param adjust_for Character vector. Terms protected from removal.
#' @param transform_registry Named list mapping base variable names to their
#'   transformation details (from \code{forward_nonlinear}).
#' @param simplify Logical. If \code{TRUE} (default), before removing a
#'   variable entirely, test whether the linear form would be retained at
#'   \code{linearity_alpha}. Implements closed testing: complex -> linear -> remove.
#' @param linearity_alpha Numeric. LRT threshold for the simplification test.
#' @param verbose Logical.
#'
#' @return Same structure as \code{scm_backward_delete()} plus
#'   \code{transform_registry} (updated) and \code{simplified_vars} (character
#'   vector of variables simplified this step).
#'
#' @details
#' \strong{Simplification hierarchy (current implementation):}
#' When \code{simplify = TRUE}, a variable in a non-linear form (e.g., FP2 or RCS)
#' is tested for simplification to the \emph{linear} form before removal is
#' considered. If the LRT p-value for non-linear vs. linear exceeds
#' \code{linearity_alpha}, the variable is simplified to linear; the (now linear)
#' variable then enters the removal phase alongside all other variables.
#'
#' This implements a \strong{two-step} closed-testing procedure:
#' \enumerate{
#'   \item \code{FP2} (or RCS) vs. \code{linear}: simplify if p > \code{linearity_alpha}.
#'   \item Linear form vs. null (removal): remove if p > \code{exit_alpha}.
#' }
#'
#' The full Royston & Altman (1994) MFP hierarchy includes an \emph{intermediate}
#' FP1 step between FP2 and linear (FP2 -> FP1 -> linear -> remove), which this
#' implementation skips. FP2 terms are simplified directly to linear; the FP1
#' intermediate step is a known limitation documented in the package vignette.
#'
#' @keywords internal
scm_backward_delete_nonlinear <- function(eres_data,
                                           current_model,
                                           model_type         = "glm",
                                           family             = "binomial",
                                           exit_alpha         = 0.01,
                                           exit_alpha_interaction = NULL,
                                           adjust_for         = NULL,
                                           transform_registry = list(),
                                           simplify           = TRUE,
                                           linearity_alpha    = 0.05,
                                           verbose            = FALSE) {
  modds       <- eres_data
  tmp_ref     <- current_model$finalforward
  toto_tot    <- current_model$toto_tot
  Relation    <- unique(current_model$Relation)
  lhs         <- .response_lhs(model_type)

  # -------------------------------------------------------------------
  # 1. Identify removable variables
  # -------------------------------------------------------------------
  # Hierarchy: protect main effects that are parents of interactions in Relation
  protected_terms <- get_protected_main_effects(Relation)
  # Also protect adjust_for
  protected_terms <- union(protected_terms,
                           if (!is.null(adjust_for)) adjust_for else character(0))

  # Formula terms currently in Relation that belong to registry variables
  reg_formula_terms <- vapply(transform_registry,
                               function(x) x$formula_rhs, character(1L),
                               USE.NAMES = TRUE)  # names = var names

  # A variable is removable when its formula_rhs is in Relation AND
  # neither its name nor its formula term is protected
  removable_vars <- names(reg_formula_terms)[
    vapply(names(reg_formula_terms), function(v) {
      frhs <- reg_formula_terms[[v]]
      frhs %in% Relation &&
        !(v %in% protected_terms) &&
        !(frhs %in% protected_terms)
    }, logical(1L), USE.NAMES = FALSE)
  ]

  if (length(removable_vars) == 0L) {
    return(list(
      finalforward       = current_model$finalforward,
      toto_tot           = toto_tot,
      Relation           = Relation,
      flagforst          = FALSE,
      transform_registry = transform_registry,
      simplified_vars    = character(0)
    ))
  }

  simplified_vars <- character(0)

  # -------------------------------------------------------------------
  # 2. Simplification phase (closed-testing: non-linear -> linear)
  # -------------------------------------------------------------------
  if (simplify) {
    for (v in removable_vars) {
      reg_entry <- transform_registry[[v]]
      if (reg_entry$transform_type == "linear") next
      if (reg_entry$df <= 1L) next  # 1-df form can't be simplified further

      old_frhs <- reg_entry$formula_rhs
      idx <- which(Relation == old_frhs)
      if (length(idx) == 0L) next

      # TODO: full MFP hierarchy would insert an FP2 -> FP1 step
      # here before the FP2 -> linear test below. Specifically:
      #   if (reg_entry$transform_type == "fp2") {
      #     best_fp1 <- .find_best_fp(v, 1L, paste(Relation[-idx], collapse=" + "), ...)
      #     fit_fp1  <- .fit_model(lhs, paste(c(Relation[-idx], best_fp1$formula_rhs), ...), ...)
      #     pval_fp2_vs_fp1 <- .lrt_pvalue(fit_full, fit_fp1, df_diff = 1L)
      #     if (!is.na(pval_fp2_vs_fp1) && pval_fp2_vs_fp1 > linearity_alpha) {
      #       # simplify to FP1, update Relation and registry, then continue
      #       # to FP1 -> linear test (which the existing logic below covers)
      #     }
      #   }

      # Build linear version of Relation for this variable
      linear_rel  <- Relation
      linear_rel[idx] <- v   # replace with bare variable name

      full_rhs   <- paste(Relation, collapse = " + ")
      linear_rhs <- paste(linear_rel, collapse = " + ")
      fit_full   <- .fit_model(lhs, full_rhs,   modds, model_type, family)
      fit_lin    <- .fit_model(lhs, linear_rhs, modds, model_type, family)

      df_diff <- reg_entry$df - 1L
      pval_simp <- .lrt_pvalue(fit_full, fit_lin, df_diff = df_diff)

      if (!is.na(pval_simp) && pval_simp > linearity_alpha) {
        # Non-linear form not significantly better than linear -> simplify
        Relation <- linear_rel
        transform_registry[[v]] <- list(
          transform_type   = "linear",
          formula_rhs      = v,
          df               = 1L,
          label            = v,
          powers           = NULL,
          knots            = NULL,
          pvalue_vs_linear = pval_simp
        )
        simplified_vars <- c(simplified_vars, v)
        if (verbose) {
          message(sprintf("  Simplified %s -> linear (p=%.4g vs linearity_alpha=%.4g)",
                          v, pval_simp, linearity_alpha))
        }
      }
    }
    # Refresh: update reg_formula_terms after any simplifications
    reg_formula_terms <- vapply(transform_registry,
                                 function(x) x$formula_rhs, character(1L),
                                 USE.NAMES = TRUE)
  }

  # -------------------------------------------------------------------
  # 3. Removal phase
  # -------------------------------------------------------------------
  full_rhs <- paste(Relation, collapse = " + ")
  fit_full <- .fit_model(lhs, full_rhs, modds, model_type, family)

  if (is.null(fit_full)) {
    return(list(
      finalforward       = current_model$finalforward,
      toto_tot           = toto_tot,
      Relation           = Relation,
      flagforst          = FALSE,
      transform_registry = transform_registry,
      simplified_vars    = simplified_vars
    ))
  }

  full_ll   <- as.numeric(stats::logLik(fit_full))
  full_aic  <- tryCatch(stats::AIC(fit_full), error = function(e) NA_real_)
  df_full   <- .model_df(fit_full, model_type)
  nobs_val  <- .model_nobs(fit_full, model_type)

  # Build per-variable removal results
  metrics   <- character(0)
  pvals     <- numeric(0)
  aic_vals  <- numeric(0)
  min2lls   <- numeric(0)
  diff2lls  <- numeric(0)
  df_tests  <- integer(0)
  df_tots   <- integer(0)
  labels    <- character(0)
  transforms_out <- character(0)

  for (v in removable_vars) {
    frhs     <- transform_registry[[v]]$formula_rhs
    df_v     <- transform_registry[[v]]$df
    lbl      <- transform_registry[[v]]$label
    tr_type  <- transform_registry[[v]]$transform_type

    reduced_rel <- Relation[Relation != frhs]
    reduced_rhs <- if (length(reduced_rel) == 0L) "1" else
      paste(reduced_rel, collapse = " + ")

    fit_red <- .fit_model(lhs, reduced_rhs, modds, model_type, family)

    pval    <- .lrt_pvalue(fit_full, fit_red, df_diff = df_v)
    ll_r    <- if (!is.null(fit_red)) tryCatch(as.numeric(stats::logLik(fit_red)), error = function(e) NA_real_) else NA_real_
    aic_r   <- if (!is.null(fit_red)) tryCatch(stats::AIC(fit_red), error = function(e) Inf) else Inf
    df_r    <- if (!is.null(fit_red)) .model_df(fit_red, model_type) else 0L

    metrics        <- c(metrics, v)
    pvals          <- c(pvals, pval)
    aic_vals       <- c(aic_vals, aic_r)
    min2lls        <- c(min2lls, if (!is.na(ll_r)) -2 * ll_r else NA_real_)
    diff2lls       <- c(diff2lls, if (!is.na(ll_r)) -2 * full_ll - (-2 * ll_r) else NA_real_)
    df_tests       <- c(df_tests, df_full - df_r)
    df_tots        <- c(df_tots, df_r)
    labels         <- c(labels, lbl)
    transforms_out <- c(transforms_out, tr_type)
  }

  if (length(metrics) == 0L) {
    return(list(
      finalforward       = current_model$finalforward,
      toto_tot           = toto_tot,
      Relation           = Relation,
      flagforst          = FALSE,
      transform_registry = transform_registry,
      simplified_vars    = simplified_vars
    ))
  }

  outtab <- tibble::tibble(
    Metric    = metrics,
    Label     = labels,
    Transform = transforms_out,
    pvalue    = pvals,
    AIC       = aic_vals,
    min2LL    = min2lls,
    diff2LL   = diff2lls,
    nobs      = nobs_val,
    dftest    = df_tests,
    dftot     = df_tots
  )

  # Per-term effective threshold
  eff_exit_alpha <- ifelse(
    grepl(":", outtab$Metric, fixed = TRUE) & !is.null(exit_alpha_interaction),
    exit_alpha_interaction, exit_alpha
  )

  # Select term to remove: highest p-value / threshold among those exceeding threshold
  exceeds_threshold <- !is.na(outtab$pvalue) & outtab$pvalue > eff_exit_alpha
  if (any(exceeds_threshold, na.rm = TRUE)) {
    rel_exceedance <- ifelse(exceeds_threshold, outtab$pvalue / eff_exit_alpha, 0)
    select_term    <- outtab$Metric[which.max(rel_exceedance)]
  } else {
    select_term <- NA_character_
  }

  # -------------------------------------------------------------------
  # 4. Build toto_back1 (compatible format with scm_backward_delete)
  # -------------------------------------------------------------------
  mod_nul_sum_aic <- full_aic
  mod_nul_sum_ll  <- full_ll
  nobs_base       <- nobs_val

  toto_back1 <- dplyr::bind_rows(
    tibble::tibble(
      Model     = "Base",
      Variable  = NA_character_,
      Transform = NA_character_,
      pvalue    = NA_real_,
      Signif    = NA_real_,
      AIC       = mod_nul_sum_aic,
      min2LL    = -2 * mod_nul_sum_ll,
      diff2LL   = 0,
      nobs      = nobs_base,
      dftest    = 0L,
      dftot     = df_full
    ),
    outtab %>%
      dplyr::rename(Model = Metric) %>%
      dplyr::mutate(
        Model     = as.character(Model),
        Variable  = Model,
        SELECT    = select_term,
        Signif    = as.numeric(!is.na(pvalue) & pvalue > eff_exit_alpha),
        Reference = tmp_ref
      )
  ) %>%
    dplyr::mutate(
      Run = -1L + (max(toto_tot$Run, na.rm = TRUE) + 1L) :
        (max(toto_tot$Run, na.rm = TRUE) + dplyr::n())
    ) %>%
    dplyr::relocate(Run) %>%
    tidyr::fill(SELECT, .direction = "downup") %>%
    dplyr::mutate(
      RUNTEST = ifelse(Run != 1L,
                       paste0("Run", Reference, " - ", Model),
                       NA_character_)
    )

  toto_back1 <- toto_back1[-1L, ]   # drop Base placeholder
  toto_tot   <- dplyr::bind_rows(toto_tot, toto_back1)

  # -------------------------------------------------------------------
  # 5. Update Relation and registry for the removed variable
  # -------------------------------------------------------------------
  flagbackst   <- !is.na(select_term)
  new_relation <- Relation
  new_registry <- transform_registry

  if (flagbackst) {
    frhs_removed   <- transform_registry[[select_term]]$formula_rhs
    new_relation   <- new_relation[new_relation != frhs_removed]
    new_registry[[select_term]] <- NULL

    finalrun <- max(
      toto_back1$Run[!is.na(toto_back1$SELECT) & toto_back1$Model == select_term],
      na.rm = TRUE
    )
    if (!is.finite(finalrun)) finalrun <- current_model$finalforward
    if (verbose) {
      message(sprintf("  Removed %s (p=%.4g > exit_alpha=%.4g)",
                      select_term, outtab$pvalue[outtab$Metric == select_term], exit_alpha))
    }
  } else {
    finalrun <- current_model$finalforward
  }

  toto_back1 <- toto_back1 %>%
    dplyr::select(Run, min2LL, diff2LL, AIC, Reference, pvalue, Signif,
                  nobs, dftest, dftot, RUNTEST, SELECT, Model, Variable, Transform)

  list(
    finalforward       = finalrun,
    toto_back1         = toto_back1,
    toto_tot           = toto_tot,
    Relation           = new_relation,
    flagforst          = flagbackst,
    transform_registry = new_registry,
    simplified_vars    = simplified_vars
  )
}


# ---------------------------------------------------------------------------
# Summary table enhancement
# ---------------------------------------------------------------------------

#' Enhance Summary Table with Transformation Labels
#'
#' Post-processes the standard summary table from \code{get_sum_tbl()} /
#' \code{get_sum_tbl_coxph()} to add transformation information and group
#' multi-coefficient terms by their parent variable.
#'
#' @param sum_tbl The standard summary data.frame from \code{get_sum_tbl()}.
#' @param transform_registry The transformation registry from
#'   \code{fit$metadata$transform_registry}.
#'
#' @return The same data.frame with three additional columns:
#'   \describe{
#'     \item{\code{variable}}{Base variable name (e.g., \code{"CMAX"}).}
#'     \item{\code{transform_type}}{Transformation type string (e.g.,
#'       \code{"log"}, \code{"fp2"}, \code{"rcs"}, \code{"linear"}).}
#'     \item{\code{transform_label}}{Human-readable label for the
#'       selected form (e.g., \code{"log(CMAX)"}, \code{"FP2(CMAX; 0.5, 2)"}).}
#'   }
#'
#' @details
#' For multi-term variables (FP2, RCS with 3+ df), multiple coefficient rows
#' will share the same \code{variable} and \code{transform_label}. Use
#' \code{transform_label} to group them for display.
#'
#' When \code{transform_registry} is \code{NULL} or empty, the input
#' \code{sum_tbl} is returned unchanged.
#'
#' @keywords internal
.enhance_summary_nonlinear <- function(sum_tbl, transform_registry) {
  if (is.null(transform_registry) || length(transform_registry) == 0L) {
    return(sum_tbl)
  }
  if (!is.data.frame(sum_tbl) || nrow(sum_tbl) == 0L) {
    return(sum_tbl)
  }

  # The first column is always "term" (coefficient name from the model)
  term_col <- sum_tbl[[1L]]

  var_out      <- character(nrow(sum_tbl))
  tr_type_out  <- character(nrow(sum_tbl))
  tr_label_out <- character(nrow(sum_tbl))

  for (i in seq_len(nrow(sum_tbl))) {
    t  <- term_col[i]
    vn <- .term_to_varname(t)

    if (!is.na(vn) && vn %in% names(transform_registry)) {
      reg <- transform_registry[[vn]]
      var_out[i]      <- vn
      tr_type_out[i]  <- reg$transform_type
      tr_label_out[i] <- reg$label
    } else {
      var_out[i]      <- vn
      tr_type_out[i]  <- "linear"
      tr_label_out[i] <- t
    }
  }

  sum_tbl$variable        <- var_out
  sum_tbl$transform_type  <- tr_type_out
  sum_tbl$transform_label <- tr_label_out
  sum_tbl
}
