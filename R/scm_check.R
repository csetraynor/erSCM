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

#' Internal: Validate a specific argument list for forward-step helpers
#'
#' Validates that the supplied list (or environment) contains the required
#' arguments with the expected types for forward-selection internal helpers.
#' Required arguments are:
#' - \code{toto_tot}: data.frame
#' - \code{Relation}: character vector
#' - \code{finalforward}: numeric (length 1)
#' - \code{flagforst}: logical (length 1)
#'
#' This is an internal helper and is not exported.
#'
#' @param args list or environment. Container expected to include the four required
#'   entries named above.
#' @param error logical. If \code{TRUE}, stop with an informative message on validation
#'   failure. If \code{FALSE}, return a structured report with \code{ok = FALSE}.
#'   Default \code{TRUE}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{ok}: logical indicating overall validation success.
#'   \item \code{args}: the original arguments (converted to list if needed).
#'   \item \code{issues}: data.frame of validation issues (name, type, message).
#'   \item \code{summary}: character message suitable for logging/errors.
#' }
#'
#' @examples
#' \dontrun{
#' x <- list(
#'   toto_tot = data.frame(a = 1),
#'   Relation = c("base -> AUC"),
#'   finalforward = 1,
#'   flagforst = TRUE
#' )
#' v <- .validate_forward_args(x, error = FALSE)
#' if (!v$ok) warning(v$summary)
#' }
#'
#' @keywords internal
#' @noRd
.validate_forward_args <- function(args, error = TRUE) {
  # Normalize args to list
  if (is.environment(args)) args <- as.list(args, all.names = TRUE)
  if (!is.list(args)) {
    msg <- "args must be a list or environment."
    if (error) cli::cli_abort(msg) else return(list(
      ok = FALSE,
      args = args,
      issues = data.frame(name = NA_character_, type = "type", message = msg),
      summary = msg
    ))
  }
  
  issues <- list()
  
  # Presence checks
  required_names <- c("toto_tot", "Relation", "finalforward", "flagforst")
  missing <- setdiff(required_names, names(args))
  for (nm in missing) {
    issues[[length(issues) + 1]] <- list(
      name = nm, type = "missing",
      message = sprintf("Missing required argument '%s'.", nm)
    )
  }
  
  # Type and length checks only if present
  # toto_tot: data.frame
  if ("toto_tot" %in% names(args)) {
    x <- args[["toto_tot"]]
    if (!is.data.frame(x)) {
      issues[[length(issues) + 1]] <- list(
        name = "toto_tot", type = "type",
        message = "Argument 'toto_tot' must be a data.frame."
      )
    }
  }
  
  # Relation: character vector
  if ("Relation" %in% names(args)) {
    x <- args[["Relation"]]
    if (!is.character(x)) {
      issues[[length(issues) + 1]] <- list(
        name = "Relation", type = "type",
        message = "Argument 'Relation' must be a character vector."
      )
    }
    # allow empty or length >= 1; no additional constraints here
  }
  
  # finalforward: numeric length 1
  if ("finalforward" %in% names(args)) {
    x <- args[["finalforward"]]
    if (!is.numeric(x)) {
      issues[[length(issues) + 1]] <- list(
        name = "finalforward", type = "type",
        message = "Argument 'finalforward' must be numeric."
      )
    } else if (length(x) != 1 || any(is.na(x))) {
      issues[[length(issues) + 1]] <- list(
        name = "finalforward", type = "length",
        message = "Argument 'finalforward' must be a non-NA scalar numeric (length 1)."
      )
    }
  }
  
  # flagforst: logical length 1
  if ("flagforst" %in% names(args)) {
    x <- args[["flagforst"]]
    if (!is.logical(x)) {
      issues[[length(issues) + 1]] <- list(
        name = "flagforst", type = "type",
        message = "Argument 'flagforst' must be logical."
      )
    } else if (length(x) != 1 || any(is.na(x))) {
      issues[[length(issues) + 1]] <- list(
        name = "flagforst", type = "length",
        message = "Argument 'flagforst' must be a non-NA scalar logical (length 1)."
      )
    }
  }
  
  # Assemble issues table
  issues_df <- if (length(issues)) {
    do.call(rbind, lapply(issues, function(z) {
      data.frame(
        name = as.character(z$name),
        type = as.character(z$type),
        message = as.character(z$message),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame(name = character(), type = character(), message = character())
  }
  
  ok <- nrow(issues_df) == 0
  summary_msg <- if (ok) {
    "Argument validation passed."
  } else {
    paste0(
      "Argument validation failed (", nrow(issues_df), " issue(s)): ",
      paste(sprintf("[%s:%s] %s", issues_df$name, issues_df$type, issues_df$message), collapse = " | ")
    )
  }
  
  if (error && !ok) cli::cli_abort(summary_msg)
  
  list(
    ok = ok,
    args = args,
    issues = issues_df,
    summary = summary_msg
  )
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
