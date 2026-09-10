# Variable Specification System for erSCM
# Provides parse_var_spec(), validate_var_spec(), and accessor/integration helpers
# that allow variable metadata (labels, units, transform info) to flow from
# scm_coxph()/scm_glm() through all downstream functions.

#' Parse variable specification into canonical format
#'
#' Accepts a data.frame, named character vector, or path to a YAML file and
#' returns a standardised tibble of variable metadata for use in the erSCM
#' analysis pipeline.
#'
#' @param spec One of:
#'   \describe{
#'     \item{data.frame}{Must have columns \code{col_name} and \code{label}.
#'       Optional: \code{units}, \code{var_type}, \code{categories} (pipe-separated
#'       string or list-column).}
#'     \item{named character vector}{Names are column names; values are display
#'       labels optionally with units in the last parenthetical, e.g.
#'       \code{c(WT = "Body Weight (kg)")}.}
#'     \item{character scalar path}{Path to a YAML file in yspec-compatible format.
#'       Requires the \pkg{yaml} package.}
#'   }
#'
#' @return A tibble with columns \code{col_name}, \code{label}, \code{units},
#'   \code{var_type}, \code{categories} (list-column), and \code{transform}
#'   (initialised to \code{"none"}).  Returns \code{NULL} when \code{spec} is
#'   \code{NULL}.
#'
#' @examples
#' # Named vector
#' spec <- c(WT = "Body Weight (kg)", AGE = "Age (years)", SEX = "Sex")
#' parse_var_spec(spec)
#'
#' # Data frame
#' spec_df <- data.frame(
#'   col_name = c("WT", "AGE"),
#'   label    = c("Body Weight", "Age"),
#'   units    = c("kg", "years")
#' )
#' parse_var_spec(spec_df)
#'
#' @export
parse_var_spec <- function(spec) {
  if (is.null(spec)) return(NULL)

  # File path: un-named scalar string pointing to an existing file
  if (is.character(spec) && length(spec) == 1 && is.null(names(spec))) {
    if (file.exists(spec)) {
      return(.parse_var_spec_yaml(spec))
    }
    cli::cli_abort("{.arg var_spec} is a length-1 character string but the file does not exist: {.path {spec}}")
  }

  # Named character vector
  if (is.character(spec) && !is.null(names(spec))) {
    return(.parse_var_spec_vector(spec))
  }

  # Data frame
  if (is.data.frame(spec)) {
    return(.parse_var_spec_df(spec))
  }

  cli::cli_abort("{.arg var_spec} must be a data.frame, named character vector, or path to a YAML file.")
}


#' Parse named character vector format
#' @keywords internal
.parse_var_spec_vector <- function(spec) {
  labels_raw <- unname(spec)
  col_names  <- names(spec)

  # Extract units from last parenthetical: "Body Weight (kg)" -> label="Body Weight", units="kg"
  has_units <- grepl("\\(([^)]+)\\)\\s*$", labels_raw)
  units  <- ifelse(has_units,
                   sub(".*\\(([^)]+)\\)\\s*$", "\\1", labels_raw),
                   NA_character_)
  labels <- ifelse(has_units,
                   trimws(sub("\\([^)]+\\)\\s*$", "", labels_raw)),
                   labels_raw)

  tibble::tibble(
    col_name   = col_names,
    label      = labels,
    units      = units,
    var_type   = NA_character_,
    categories = vector("list", length(col_names)),
    transform  = "none"
  )
}


#' Parse data.frame format
#' @keywords internal
.parse_var_spec_df <- function(spec) {
  if (!"col_name" %in% colnames(spec)) {
    cli::cli_abort("{.arg var_spec} data.frame must have a {.col col_name} column.")
  }
  if (!"label" %in% colnames(spec)) {
    cli::cli_abort("{.arg var_spec} data.frame must have a {.col label} column.")
  }

  out <- tibble::tibble(
    col_name = as.character(spec$col_name),
    label    = as.character(spec$label),
    units    = if ("units" %in% colnames(spec)) as.character(spec$units) else NA_character_,
    var_type = if ("var_type" %in% colnames(spec)) as.character(spec$var_type) else NA_character_,
    transform = "none"
  )

  # Parse categories column if present
  if ("categories" %in% colnames(spec)) {
    out$categories <- lapply(spec$categories, function(x) {
      if (is.null(x) || (length(x) == 1 && is.na(x))) return(NULL)
      if (is.character(x) && length(x) == 1 && grepl("\\|", x)) {
        vals <- trimws(strsplit(x, "\\|")[[1]])
        return(vals)
      }
      x
    })
  } else {
    out$categories <- vector("list", nrow(out))
  }

  out
}


#' Parse YAML file format
#' @keywords internal
.parse_var_spec_yaml <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg yaml} is required to read YAML var_spec files.",
      "i" = "Install with: {.run install.packages('yaml')}"
    ))
  }

  raw <- yaml::read_yaml(path)

  rows <- lapply(names(raw), function(varname) {
    entry <- raw[[varname]]
    label_val <- if (!is.null(entry$label)) entry$label else varname
    units_val <- if (!is.null(entry$units)) entry$units else NA_character_
    type_val  <- if (!is.null(entry$type))  entry$type  else NA_character_
    cats      <- if (!is.null(entry$categories)) entry$categories else NULL

    tibble::tibble(
      col_name   = varname,
      label      = label_val,
      units      = units_val,
      var_type   = type_val,
      categories = list(cats),
      transform  = "none"
    )
  })

  dplyr::bind_rows(rows)
}


#' Validate var_spec against the analysis dataset
#'
#' Checks that specified column names exist in the data, auto-detects
#' \code{var_type} where not provided, and emits a message for model variables
#' not covered by the spec (adding minimal auto-entries for them).
#'
#' @param var_spec Canonical var_spec tibble from \code{\link{parse_var_spec}}.
#' @param data The analysis data frame.
#' @param covariates Character vector of covariate names used in the model.
#' @param exposure_vars Character vector of exposure variable names used in the model.
#'
#' @return Validated and completed var_spec tibble.
#' @keywords internal
validate_var_spec <- function(var_spec, data, covariates, exposure_vars) {
  all_model_vars <- unique(c(covariates, exposure_vars))

  # Drop entries that reference columns not in data
  missing_cols <- setdiff(var_spec$col_name, colnames(data))
  if (length(missing_cols) > 0) {
    cli::cli_warn("{.arg var_spec} references columns not found in data: {.val {missing_cols}}.")
    var_spec <- var_spec[var_spec$col_name %in% colnames(data), ]
  }

  # Auto-detect var_type where not specified
  for (i in seq_len(nrow(var_spec))) {
    if (is.na(var_spec$var_type[i])) {
      col <- data[[var_spec$col_name[i]]]
      var_spec$var_type[i] <- if (is.numeric(col)) "continuous" else "categorical"
    }
  }

  # Warn about model variables not in spec; add minimal auto-entries
  uncovered <- setdiff(all_model_vars, var_spec$col_name)
  uncovered <- uncovered[uncovered %in% colnames(data)]
  if (length(uncovered) > 0) {
    cli::cli_inform("Variables not in var_spec (will use column name as label): {.val {uncovered}}.")
    is_num <- vapply(data[uncovered], is.numeric, logical(1))
    auto_entries <- tibble::tibble(
      col_name   = uncovered,
      label      = uncovered,
      units      = ifelse(is_num, "U", NA_character_),
      var_type   = ifelse(is_num, "continuous", "categorical"),
      categories = vector("list", length(uncovered)),
      transform  = "none"
    )
    var_spec <- dplyr::bind_rows(var_spec, auto_entries)
  }

  var_spec
}


#' Get display label for a variable
#'
#' Looks up the human-readable label for a column from a canonical
#' \code{var_spec} tibble (as produced by \code{\link{parse_var_spec}} or
#' stored on a fitted model as \code{fit$var_spec}). Useful for building
#' custom summaries or report templates that need consistent variable labels.
#'
#' @param var_spec Canonical var_spec tibble, or \code{NULL}.
#' @param col_name Character; the column name to look up.
#' @return A length-1 character string: the display label for \code{col_name},
#'   or \code{col_name} itself if \code{var_spec} is \code{NULL} or does not
#'   contain an entry for it.
#' @export
#' @examples
#' spec <- parse_var_spec(c(WT = "Body Weight (kg)", AGE = "Age (years)"))
#' get_var_label(spec, "WT")
#' get_var_label(spec, "MISSING_COL")
get_var_label <- function(var_spec, col_name) {
  if (is.null(var_spec)) return(col_name)
  idx <- match(col_name, var_spec$col_name)
  if (is.na(idx)) return(col_name)
  var_spec$label[idx]
}


#' Get units string for a variable
#'
#' Looks up the measurement units for a column from a canonical
#' \code{var_spec} tibble (as produced by \code{\link{parse_var_spec}} or
#' stored on a fitted model as \code{fit$var_spec}). Useful for building
#' custom summaries or report templates that need consistent variable units.
#'
#' @param var_spec Canonical var_spec tibble, or \code{NULL}.
#' @param col_name Character; the column name to look up.
#' @return A length-1 character string with the units for \code{col_name}, or
#'   \code{NA_character_} if \code{var_spec} is \code{NULL} or does not
#'   contain an entry for it.
#' @export
#' @examples
#' spec <- parse_var_spec(c(WT = "Body Weight (kg)", AGE = "Age (years)"))
#' get_var_units(spec, "WT")
#' get_var_units(spec, "MISSING_COL")
get_var_units <- function(var_spec, col_name) {
  if (is.null(var_spec)) return(NA_character_)
  idx <- match(col_name, var_spec$col_name)
  if (is.na(idx)) return(NA_character_)
  var_spec$units[idx]
}


#' Resolve variable metadata for forest plot labeling
#'
#' Determines which variable metadata source to use, in priority order:
#' 1. Explicit \code{namings_override} (old-style or var_spec format),
#' 2. \code{fit$var_spec} (new system),
#' 3. \code{fit$naming_df} (legacy system),
#' 4. Centering values only (minimal auto-constructed entries).
#'
#' @param fit A fitted erSCM object.
#' @param namings_override Optional data.frame; overrides all other sources.
#' @return A canonical var_spec tibble or \code{NULL}.
#' @keywords internal
resolve_var_metadata <- function(fit, namings_override = NULL) {
  # Priority 1: explicit override
  if (!is.null(namings_override)) {
    if (is.data.frame(namings_override) && "covname" %in% colnames(namings_override)) {
      return(.convert_legacy_namings(namings_override))
    }
    if (is.data.frame(namings_override) && "col_name" %in% colnames(namings_override)) {
      return(namings_override)
    }
  }

  # Priority 2: fit$var_spec
  if (!is.null(fit$var_spec)) {
    return(fit$var_spec)
  }

  # Priority 3: fit$naming_df (legacy)
  if (!is.null(fit$naming_df)) {
    return(.convert_legacy_namings(fit$naming_df))
  }

  # Priority 4: centering_values only
  cv <- fit$metadata$centering_values
  if (!is.null(cv) && length(cv) > 0) {
    model_vars <- get_main_effect_columns(fit$metadata$finalRelation)
    cv <- cv[names(cv) %in% model_vars]
    if (length(cv) > 0) {
      return(tibble::tibble(
        col_name     = names(cv),
        label        = names(cv),
        units        = "U",
        var_type     = "continuous",
        median_point = unname(cv),
        transform    = "centered"
      ))
    }
  }

  NULL
}


#' Convert legacy naming_df to canonical var_spec format
#' @keywords internal
.convert_legacy_namings <- function(naming_df) {
  has_median <- "median_point" %in% colnames(naming_df)
  has_units  <- "units" %in% colnames(naming_df)
  has_ttype  <- "transform_type" %in% colnames(naming_df)

  mp_raw <- if (has_median) naming_df$median_point else NA_real_

  # Determine transform type
  if (has_ttype) {
    ttype <- naming_df$transform_type
    # Coerce character median_point to numeric for old log_transformed format
    mp <- if (is.character(mp_raw)) {
      ifelse(is.na(mp_raw), NA_real_, 1.0)
    } else {
      as.numeric(mp_raw)
    }
  } else {
    if (is.character(mp_raw)) {
      ttype <- "log_transformed"
      mp    <- ifelse(is.na(mp_raw), NA_real_, 1.0)
    } else {
      ttype <- ifelse(!is.na(mp_raw), "median_scaled", "none")
      mp    <- as.numeric(mp_raw)
    }
  }

  tibble::tibble(
    col_name     = naming_df$covname,
    label        = naming_df$covname,
    units        = if (has_units) naming_df$units else "U",
    var_type     = "continuous",
    median_point = mp,
    transform    = ttype
  )
}


#' Apply variable labels and back-transformation to forest plot data
#'
#' Uses a resolved var_spec (canonical tibble) to back-transform quantile
#' values to the original variable scale and attach display labels.
#'
#' @param fp_data Forest plot tibble with at minimum columns
#'   \code{covname}, \code{label}, \code{VAL}.
#' @param var_meta Canonical var_spec tibble from
#'   \code{\link{resolve_var_metadata}}, or \code{NULL}.
#' @return The \code{fp_data} tibble with added/modified columns
#'   \code{units}, \code{group_level}, and \code{covname} (cleaned).
#' @keywords internal
apply_var_labels <- function(fp_data, var_meta) {
  if (is.null(var_meta)) {
    return(fp_data %>%
      dplyr::mutate(
        units           = ifelse(is.na(VAL), NA_character_, "U"),
        group_level     = ifelse(is.na(VAL), label, paste(round(VAL, 3), "U")),
        covname         = as.character(mapply(gsub, label, "", covname, MoreArgs = list(fixed = TRUE), USE.NAMES = FALSE)),
        display_covname = covname
      ))
  }

  # Ensure median_point column exists in var_meta
  if (!"median_point" %in% colnames(var_meta)) {
    var_meta <- var_meta %>% dplyr::mutate(median_point = NA_real_)
  }

  fp_data %>%
    dplyr::left_join(
      var_meta %>%
        dplyr::select(dplyr::any_of(c("col_name", "label", "units",
                                       "median_point", "transform"))) %>%
        dplyr::rename(
          covname_key   = col_name,
          display_label = label
        ),
      by = c("covname" = "covname_key")
    ) %>%
    dplyr::mutate(
      # Back-transform VAL to original scale for display labels
      VAL_display = dplyr::case_when(
        is.na(VAL)                                                          ~ NA_real_,
        is.na(transform) | transform == "none"                              ~ VAL,
        transform == "median_scaled"                                        ~ VAL * median_point,
        transform == "centered"                                             ~ VAL + median_point,
        transform == "median_scaled_and_centered"                           ~ (VAL + 1.0) * median_point,
        transform == "log_transformed"                                      ~ exp(VAL),
        TRUE                                                                ~ VAL
      ),
      units = dplyr::coalesce(units, "U"),
      group_level = dplyr::case_when(
        is.na(VAL_display) ~ label,
        is.na(units) | units == "U" ~ paste(round(VAL_display, 3), "U"),
        TRUE ~ paste(round(VAL_display, 3), units)
      ),
      # Clean categorical covname suffixes (e.g. "SEXFemale" -> "SEX")
      covname = as.character(mapply(gsub, label, "", covname, MoreArgs = list(fixed = TRUE), USE.NAMES = FALSE)),
      # display_covname: use display label for continuous, keep raw for categorical
      display_covname = dplyr::case_when(
        !is.na(display_label) & !is.na(VAL_display) ~ display_label,
        TRUE ~ covname
      )
    ) %>%
    dplyr::select(-dplyr::any_of(c("display_label", "median_point", "transform",
                                    "VAL_display")))
}


#' Get category display labels for a variable
#'
#' @param var_spec Canonical var_spec tibble, or \code{NULL}.
#' @param col_name Character; the column name to look up.
#' @return Character vector of category labels (in order), or \code{NULL} if not
#'   found or no categories are defined.
#' @keywords internal
get_var_categories <- function(var_spec, col_name) {
  if (is.null(var_spec)) return(NULL)
  idx <- match(col_name, var_spec$col_name)
  if (is.na(idx)) return(NULL)
  cats <- var_spec$categories[[idx]]
  if (is.null(cats) || length(cats) == 0) return(NULL)
  unname(cats)
}


#' Build display-friendly bin labels using var_spec back-transformation
#'
#' Given quantile break points on the (potentially transformed) analysis scale,
#' back-transforms them to the original variable scale and formats them as
#' range strings (e.g., "<= 95.2", "95.2-45.3", "> 198.7").
#'
#' @param breaks Numeric vector of break points including the outer boundaries
#'   (length = n_bins + 1).
#' @param var_spec Canonical var_spec tibble, or \code{NULL}.
#' @param col_name Character; the column name whose transform to look up.
#' @return Character vector of length \code{length(breaks) - 1}.
#' @keywords internal
make_bin_labels <- function(breaks, var_spec, col_name) {
  back_fn <- .get_back_transform_fn(var_spec, col_name)
  original_breaks <- back_fn(breaks)

  n_bins <- length(original_breaks) - 1
  labels <- character(n_bins)
  for (i in seq_len(n_bins)) {
    lo <- round(original_breaks[i], 1)
    hi <- round(original_breaks[i + 1], 1)
    if (i == 1) {
      labels[i] <- paste0("\u2264 ", hi)
    } else if (i == n_bins) {
      labels[i] <- paste0("> ", lo)
    } else {
      labels[i] <- paste0(lo, "\u2013", hi)
    }
  }
  labels
}


#' Build a back-transformation function for a variable
#' @keywords internal
.get_back_transform_fn <- function(var_spec, col_name) {
  if (is.null(var_spec)) return(identity)
  idx <- match(col_name, var_spec$col_name)
  if (is.na(idx)) return(identity)

  transform <- if ("transform" %in% colnames(var_spec)) var_spec$transform[idx] else "none"
  mp <- if ("median_point" %in% colnames(var_spec)) var_spec$median_point[idx] else NA_real_

  if (is.na(transform) || transform == "none") return(identity)
  if (transform == "centered")    return(function(x) x + mp)
  if (transform == "median_scaled") return(function(x) x * mp)
  if (transform == "median_scaled_and_centered") return(function(x) (x + 1.0) * mp)
  if (transform == "log_transformed") return(exp)

  # Unrecognized transform: previously fell through to identity(), silently
  # producing wrong-scale axis/bin labels with no diagnostic. The set of valid
  # transforms is closed (set internally by scm_glm/scm_coxph and parse_var_spec)
  # so an unknown value here indicates a typo or an unsupported spec, which we
  # surface loudly rather than mislabel.
  known_transforms <- c("none", "centered", "median_scaled",
                        "median_scaled_and_centered", "log_transformed")
  cli::cli_abort(c(
    "Unrecognized {.field transform} {.val {transform}} for column {.val {col_name}}.",
    "i" = "Known transforms: {.val {known_transforms}}."
  ))
}
