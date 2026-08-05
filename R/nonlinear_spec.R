# Nonlinear specification: parsing and validation
# Translates the user-facing `nonlinear` argument into a canonical per-variable
# list of transformation types used throughout the SCM forward/backward steps.

#' Allowed nonlinear transformation names
#' @keywords internal
.ALLOWED_TRANSFORMS <- c("linear", "log", "sqrt", "fp1", "fp2", "rcs")

#' Default transformation set for auto mode
#' @keywords internal
.DEFAULT_TRANSFORMS <- c("linear", "log", "fp1", "fp2")


#' Validate a vector of transformation type names
#'
#' @param transforms Character vector of transformation names to check.
#' @param var_name Optional variable name for cleaner error messages.
#' @return Invisible TRUE, or aborts on invalid names.
#' @keywords internal
.validate_transform_names <- function(transforms, var_name = NULL) {
  is_valid <- transforms %in% .ALLOWED_TRANSFORMS |
    grepl("^rcs[0-9]+$", transforms)
  invalid <- transforms[!is_valid]
  if (length(invalid) > 0L) {
    prefix <- if (!is.null(var_name)) paste0(" for variable ", var_name) else ""
    allowed <- .ALLOWED_TRANSFORMS
    cli::cli_abort(
      c(
        "Unknown transformation type(s){prefix}: {.val {invalid}}.",
        "i" = "Allowed values are {.val {allowed}} and
          {.code \"rcsN\"} where N is the knot count (e.g. {.code \"rcs4\"})."
      )
    )
  }
  invisible(TRUE)
}


#' Parse Non-Linear Transformation Specification
#'
#' Parses the user-provided \code{nonlinear} argument into a canonical internal
#' format used by the SCM forward/backward steps.
#'
#' @param nonlinear The user's non-linear specification. Accepts multiple formats:
#'   \describe{
#'     \item{NULL}{No non-linear terms (default). All variables enter linearly.
#'       Existing behavior is preserved exactly.}
#'     \item{TRUE or "auto"}{Apply default transformation search
#'       (\code{c("linear", "log", "fp1", "fp2")}) to all continuous candidate
#'       variables. Factor/character variables are always linear.}
#'     \item{Character vector}{A set of transformation names to apply to ALL
#'       continuous candidates. E.g., \code{c("linear", "log", "rcs")} means
#'       test these three forms for every continuous variable.}
#'     \item{Named list}{Per-variable specification. Names are variable names;
#'       values are character vectors of transformations to test for that variable.
#'       Variables not listed default to \code{"linear"} only.
#'       E.g., \code{list(CMAX = c("linear", "log", "fp2"), AGE = c("linear", "log"))}}
#'   }
#' @param covariates Character vector of candidate covariate names.
#' @param exposure_vars Character vector of exposure variable names.
#' @param data The modelling dataset (to determine which variables are numeric).
#'
#' @return A named list of class \code{"nonlinear_spec"} where each element is
#'   named by a variable and contains a character vector of transformation types
#'   to test. Variables not in the list are implicitly \code{"linear"} only.
#'   Returns \code{NULL} if no non-linear terms are requested.
#'
#' @details
#' Factor, character, and logical variables are ALWAYS treated as linear
#' (categorical) regardless of the specification. If a non-numeric variable
#' is explicitly given a non-linear transformation, a warning is emitted and
#' it is forced to \code{"linear"}.
#'
#' @examples
#' data(demo_eres)
#'
#' # All continuous variables get the full transform search
#' spec <- parse_nonlinear_spec(TRUE, covariates = c("AGE", "WT"),
#'   exposure_vars = "CMAXss", data = demo_eres)
#'
#' # Per-variable specification
#' spec <- parse_nonlinear_spec(
#'   list(CMAXss = c("linear", "log", "fp2"), AGE = c("linear", "log")),
#'   covariates = c("AGE", "WT"),
#'   exposure_vars = "CMAXss", data = demo_eres
#' )
#'
#' @keywords internal
#' @export
parse_nonlinear_spec <- function(nonlinear, covariates, exposure_vars, data) {

  # --- NULL / FALSE: disable non-linear search entirely ---
  if (is.null(nonlinear) || identical(nonlinear, FALSE)) {
    return(NULL)
  }

  all_vars <- unique(c(covariates, exposure_vars))

  # Identify which variables are numeric (continuous)
  is_numeric_var <- vapply(all_vars, function(v) {
    if (!v %in% colnames(data)) return(FALSE)
    is.numeric(data[[v]])
  }, logical(1L), USE.NAMES = TRUE)

  # --- Build raw spec (before numeric enforcement) ---

  if (isTRUE(nonlinear) ||
      (is.character(nonlinear) && length(nonlinear) == 1L && nonlinear == "auto")) {
    # Auto mode: default transforms for continuous, "linear" for categorical
    raw_spec <- stats::setNames(
      lapply(all_vars, function(v) {
        if (is_numeric_var[[v]]) .DEFAULT_TRANSFORMS else "linear"
      }),
      all_vars
    )

  } else if (is.character(nonlinear)) {
    # Character vector: same transform set for all continuous variables
    .validate_transform_names(nonlinear)
    raw_spec <- stats::setNames(
      lapply(all_vars, function(v) {
        if (is_numeric_var[[v]]) nonlinear else "linear"
      }),
      all_vars
    )

  } else if (is.list(nonlinear)) {
    # Named list: per-variable specification
    if (is.null(names(nonlinear))) {
      cli::cli_abort(
        "{.arg nonlinear} as a list must be named (each name is a variable name)."
      )
    }

    # Warn about variables listed that are not known candidates
    unknown_vars <- setdiff(names(nonlinear), all_vars)
    if (length(unknown_vars) > 0L) {
      cli::cli_warn(
        c(
          "The following variable(s) in {.arg nonlinear} are not in
            {.arg covariates} or {.arg exposure_vars}: {.val {unknown_vars}}.",
          "i" = "They will be ignored."
        )
      )
    }

    # Validate transform names per variable
    for (nm in intersect(names(nonlinear), all_vars)) {
      .validate_transform_names(nonlinear[[nm]], var_name = nm)
    }

    # Build spec: explicitly listed variables use their spec; others default to "linear"
    raw_spec <- stats::setNames(
      lapply(all_vars, function(v) {
        if (v %in% names(nonlinear)) nonlinear[[v]] else "linear"
      }),
      all_vars
    )

  } else {
    cli::cli_abort(
      c(
        "{.arg nonlinear} must be {.code NULL}, {.code FALSE}, {.code TRUE},
          {.code \"auto\"}, a character vector, or a named list.",
        "x" = "Got {.cls {class(nonlinear)}}."
      )
    )
  }

  # --- Enforce: non-numeric variables must be "linear" ---
  for (v in all_vars) {
    if (!is_numeric_var[[v]]) {
      has_nonlinear_requested <- any(raw_spec[[v]] != "linear")
      if (has_nonlinear_requested && is.list(nonlinear) && v %in% names(nonlinear)) {
        # Only warn when the user explicitly requested non-linear for a non-numeric var
        cli::cli_warn(
          c(
            "Variable {.val {v}} is not numeric; non-linear transformations are not applicable.",
            "i" = "Forcing to {.code \"linear\"} only."
          )
        )
      }
      raw_spec[[v]] <- "linear"
    }
  }

  class(raw_spec) <- "nonlinear_spec"
  raw_spec
}


#' Validate a Nonlinear Specification
#'
#' Checks that variables referenced in the spec exist in the data, and emits
#' warnings when data properties conflict with the requested transformations
#' (e.g., log of non-positive values).
#'
#' @param spec A parsed \code{nonlinear_spec} object (from
#'   \code{\link{parse_nonlinear_spec}}), or \code{NULL}.
#' @param data The modelling dataset.
#'
#' @return The \code{spec} invisibly, or aborts if variables are missing.
#' @keywords internal
validate_nonlinear_spec <- function(spec, data) {
  if (is.null(spec)) return(invisible(NULL))

  var_names <- names(spec)

  # Check all named variables exist in data
  missing_vars <- setdiff(var_names, colnames(data))
  if (length(missing_vars) > 0L) {
    cli::cli_abort(
      "Variable(s) in nonlinear_spec not found in {.arg data}: {.val {missing_vars}}."
    )
  }

  for (v in var_names) {
    transforms <- spec[[v]]
    x      <- data[[v]]
    x_obs  <- x[!is.na(x)]

    if (length(x_obs) == 0L) next

    # log requires strictly positive values
    if ("log" %in% transforms && any(x_obs <= 0)) {
      n_bad <- sum(x_obs <= 0)
      cli::cli_warn(
        c(
          "Variable {.val {v}} has {n_bad} non-positive value(s);
            the {.code \"log\"} transform requires x > 0.",
          "i" = "Log-transform fits for {.val {v}} may fail or produce {.code NaN}."
        )
      )
    }

    # sqrt requires non-negative values
    if ("sqrt" %in% transforms && any(x_obs < 0)) {
      n_bad <- sum(x_obs < 0)
      cli::cli_warn(
        c(
          "Variable {.val {v}} has {n_bad} negative value(s);
            the {.code \"sqrt\"} transform requires x >= 0.",
          "i" = "Square-root fits for {.val {v}} may fail or produce {.code NaN}."
        )
      )
    }

    # FP and RCS need sufficient unique values for reliable estimation
    fp_rcs_requested <- any(c("fp1", "fp2") %in% transforms) ||
      any(grepl("^rcs", transforms))
    if (fp_rcs_requested) {
      n_unique <- length(unique(x_obs))
      if (n_unique < 5L) {
        cli::cli_warn(
          c(
            "Variable {.val {v}} has only {n_unique} unique value(s).",
            "i" = "FP and RCS transformations may be unreliable with fewer than 5
              unique values."
          )
        )
      }
    }
  }

  invisible(spec)
}
