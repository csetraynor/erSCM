#' Internal: Check forward step status and update model relation
#'
#' Uses the accumulated summary of candidate runs and the current model relation
#' (e.g., edges between base and candidate models) to decide the next action in
#' the forward-selection phase. Updates the relation with the chosen candidate,
#' determines the provisional final forward model (if any), and sets a flag
#' indicating whether forward selection should continue.
#'
#' This is an internal helper and is not exported.
#'
#' @param run_summary tibble or data.frame. Aggregated summary of forward-run
#'   candidates, typically one row per candidate with columns such as \code{term},
#'   \code{estimate}, \code{se}, \code{criterion} (e.g., AIC/BIC), \code{p_value}
#'   (if LRT-based), and \code{enter} (logical).
#' @param relation list or data structure. Current model relation/graph capturing
#'   the base model, tested candidates, and accepted terms; commonly includes fields
#'   like \code{base_terms}, \code{accepted_terms}, and an edge or history log.
#' @param criterion character. Selection criterion to prioritize when choosing the
#'   next candidate. One of \code{"AIC"}, \code{"BIC"}, or \code{"LRT"}. Default \code{"AIC"}.
#' @param entry_alpha numeric. Entry threshold for \code{"LRT"} criterion; ignored
#'   for information-criteria. Default \code{0.01}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{run_summary}: possibly updated/annotated run summary (e.g., ranked order).
#'   \item \code{relation}: updated relation reflecting any accepted candidate.
#'   \item \code{final_forward_model}: model object or structured spec representing the
#'         model after the current forward check (may be \code{NULL} if none accepted).
#'   \item \code{continue}: logical flag indicating whether to continue forward selection.
#' }
#'
#' @details
#' The function ranks candidates based on \code{criterion}. For \code{"LRT"}, it
#' evaluates the smallest p-value against \code{entry_alpha}. For \code{"AIC"} or
#' \code{"BIC"}, it prefers the smallest criterion value and may require an improvement
#' over the current best model recorded in \code{relation}. If no candidate passes,
#' \code{continue} is set to \code{FALSE}. Missing values in \code{criterion}/\code{p_value}
#' are treated as non-qualifying.
#'
#' @examples
#' \dontrun{
#' out <- .check_forward(
#'   run_summary = cand_summary,
#'   relation = current_relation,
#'   criterion = "AIC",
#'   entry_alpha = 0.01
#' )
#' if (out$continue) {
#'   # proceed to next forward iteration
#' }
#' }
#'
#' @keywords internal
#' @noRd
.check_forward <- function(
    run_summary,
    relation
) {
  # implementation goes here
  x <- run_summary
  Relation <- relation
  if (!all(is.na(x$Run[x$Model == x$SELECT]))) {
    finalforward <- x$Run[x$Model == x$SELECT][1L]
    Relation <- c(Relation, x$Model[x$Model == x$SELECT])
  } else {
    finalforward <- NA_real_
    Relation <- Relation
  }
  flagforst <- FALSE
  if (any(!is.na(x$Run[x$Model == x$SELECT]))) {
    flagforst <- TRUE
  }
  
  if(!flagforst) {
    message("Exposure relationship not significant at the selected level.")
  }
  
  return(list(toto_tot = x, Relation = Relation, finalforward =finalforward, flagforst = flagforst  ) )
}


#' Check inputs for GLM SCM functions
#'
#' @param modds Data frame.
#' @param dvid_col Character or NULL.
#' @param dvid Character or NULL.
#' @param dv_col Character.
#' @param covariates Character vector.
#' @keywords internal
#' @noRd
.check_inputs <- function(modds, dvid_col, dvid, dv_col, covariates) {
  if (!is.null(dvid_col)) {
    if (!is.character(dvid_col)) {
      cli::cli_abort("{.arg dvid_col} must be NULL or a character string.")
    } else {
      if (is.null(dvid)) {
        cli::cli_abort("{.arg dvid} must be provided when {.arg dvid_col} is specified.")
      }
    }
  } else {
    if (!is.null(dvid)) {
      cli::cli_abort("{.arg dvid} is set but {.arg dvid_col} is NULL; provide {.arg dvid_col} to filter the dataset.")
    } else {
      cli::cli_inform("{.arg dvid_col} is NULL; model run without filtering the dataset by DVID.")
    }
  }

  if (is.null(dv_col)) {
    cli::cli_warn("{.arg dv_col} is NULL; assuming the outcome column is named {.val DV}.")
    if ("DV" %!in% colnames(modds)) {
      cli::cli_abort("Column {.val DV} not found in {.arg eres_data}. Please supply {.arg dv_col}.")
    }
  }

  missing_covs <- covariates[covariates %!in% colnames(modds)]
  if (length(missing_covs) > 0) {
    cli::cli_abort(c(
      "The following covariates are not in {.arg eres_data}:",
      "x" = "{.val {missing_covs}}"
    ))
  }
}



#' Check inputs for survival SCM functions
#'
#' @param modds Data frame.
#' @param dvid_col Character or NULL.
#' @param dvid Character or NULL.
#' @param dv_col Character.
#' @param time_col Character.
#' @param covariates Character vector.
#' @keywords internal
#' @noRd
.check_inputs_surv <- function(modds, dvid_col, dvid, dv_col, time_col, covariates) {
  if (!is.null(dvid_col)) {
    if (!is.character(dvid_col)) {
      cli::cli_abort("{.arg dvid_col} must be NULL or a character string.")
    } else {
      if (is.null(dvid)) {
        cli::cli_abort("{.arg dvid} must be provided when {.arg dvid_col} is specified.")
      }
    }
  } else {
    if (!is.null(dvid)) {
      cli::cli_abort("{.arg dvid} is set but {.arg dvid_col} is NULL; provide {.arg dvid_col} to filter the dataset.")
    } else {
      cli::cli_inform("{.arg dvid_col} is NULL; model run without filtering the dataset by DVID.")
    }
  }

  if (is.null(dv_col)) {
    cli::cli_warn("{.arg dv_col} is NULL; assuming the outcome column is named {.val DV}.")
    if ("DV" %!in% colnames(modds)) {
      cli::cli_abort("Column {.val DV} not found in {.arg eres_data}. Please supply {.arg dv_col}.")
    }
  }

  if (is.null(time_col)) {
    cli::cli_warn("{.arg time_col} is NULL; assuming the time column is named {.val TIME}.")
    if ("TIME" %!in% colnames(modds)) {
      cli::cli_abort("Column {.val TIME} not found in {.arg eres_data}. Please supply {.arg time_col}.")
    }
  }

  missing_covs <- covariates[covariates %!in% colnames(modds)]
  if (length(missing_covs) > 0) {
    cli::cli_abort(c(
      "The following covariates are not in {.arg eres_data}:",
      "x" = "{.val {missing_covs}}"
    ))
  }
}
