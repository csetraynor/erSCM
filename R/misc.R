#' Not-in operator
#'
#' Returns a logical vector indicating which elements of \code{x}
#' are not in \code{y}.
#'
#' @param x Vector of values to check.
#' @param y Vector of values to check against.
#' @return A logical vector.
#' @name grapes-not-in-grapes
#' @export
'%!in%' <- function(x,y)!('%in%'(x,y))



#' Pad numbers to a fixed number of significant digits (without table1::signif_pad)
#'
#' Formats numeric values to a specified number of significant digits while
#' preserving trailing zeros, so the printed width reflects precision
#' (e.g., \code{1.20} instead of \code{1.2}). Uses fixed notation within a
#' configurable range and scientific notation outside that range. Handles
#' \code{NA}, \code{Inf}, \code{-Inf}, and \code{0} explicitly.
#'
#' @param x Numeric vector to format.
#' @param digits Integer. Number of significant digits to keep (default \code{3}).
#' @param sci_thresh Named numeric vector of length 2 giving thresholds for
#'   scientific notation; values with \code{|x| < sci_thresh["low"]} or
#'   \code{|x| >= sci_thresh["high"]} are printed in scientific notation.
#'   Default is \code{c(low = 1e-4, high = 1e6)}.
#' @param trim Logical. If \code{TRUE}, leading/trailing whitespace is removed
#'   from the formatted output (default \code{FALSE}).
#'
#' @return A character vector the same length as \code{x} with numbers formatted
#'   to \code{digits} significant digits and padded with trailing zeros as needed.
#'
#' @details
#' For values printed in fixed notation, the number of decimal places is computed
#' as \code{max(digits - 1 - floor(log10(abs(x))), 0)} and the value is rounded
#' accordingly, then formatted with \code{sprintf} to preserve trailing zeros.
#' For scientific notation, the mantissa has \code{digits - 1} decimal places
#' (e.g., for \code{digits = 3}, \code{1.23e+04}). Zero is shown as
#' \code{"0."} followed by \code{digits - 1} zeros. \code{NA}, \code{Inf}, and
#' \code{-Inf} are returned as strings \code{"NA"}, \code{"Inf"}, and \code{"-Inf"}.
#'
#' @examples
#' x <- c(12345, 123.456, 1.23456, 0.0123456, 0.0000123456, 0, NA, Inf, -Inf, -0.004321)
#' signif_pad_like(x, digits = 3)
#' signif_pad_like(x, digits = 4)
#'
#' # Adjust thresholds for scientific notation
#' signif_pad_like(c(1e-6, 1e-5, 1e6, 1e7), digits = 3,
#'                 sci_thresh = c(low = 1e-6, high = 1e6))
#'
#' # Trim whitespace if any
#' signif_pad_like(x, digits = 3, trim = TRUE)
#'
#' @export
signif_pad_like <- function(x, digits = 3, sci_thresh = c(low = 1e-4, high = 1e6),
                            trim = FALSE) {
  # x: numeric vector
  # digits: number of significant digits to keep
  # sci_thresh: range where fixed notation is used; outside uses scientific notation
  # trim: if TRUE, remove leading spaces (from sprintf/scientific formats)
  
  # Helper: format one number with padded significant digits
  format_one <- function(val) {
    if (is.na(val)) return("NA")
    if (!is.finite(val)) return(ifelse(val > 0, "Inf", "-Inf"))
    if (val == 0) {
      # zero is special: pad with digits-1 zeros
      # e.g., digits=3 => "0.00"
      dec <- max(digits - 1, 0)
      return(sprintf(paste0("%.", dec, "f"), 0))
    }
    
    s <- abs(val)
    # Decide notation: fixed vs scientific
    use_scientific <- (s < sci_thresh["low"]) || (s >= sci_thresh["high"])
    
    if (!use_scientific) {
      # Fixed notation: compute decimal places needed for desired significant digits
      # decimal_places = max(digits - 1 - floor(log10(|x|)), 0)
      dp <- max(digits - 1 - floor(log10(s)), 0)
      # Round to dp decimals and format with zero padding
      out <- sprintf(paste0("%.", dp, "f"), round(val, dp))
    } else {
      # Scientific notation: format with (digits-1) decimals in mantissa
      # Example: digits=3 => mantissa 1.23e+04
      dp <- max(digits - 1, 0)
      # Use sprintf with %e and then normalize exponent formatting
      out <- sprintf(paste0("%.", dp, "e"), val)
      # Make exponent always sign-padded to at least 2 digits (e.g., e+01, e-06)
      # Base R already uses that, but ensure no leading spaces:
      out <- sub("e\\+0", "e+0", out, fixed = FALSE)  # no-op mostly, retained for clarity
    }
    
    if (trim) out <- trimws(out)
    out
  }
  
  # Vectorized
  vapply(x, format_one, FUN.VALUE = character(1))
}


#' Format axis labels with significant digits
#'
#' Returns a function that formats numeric values to a given number
#' of significant digits, suitable for use as an axis label formatter.
#'
#' @param decimals Integer; number of significant digits.
#' @return A function that takes a numeric vector and returns a character
#'   vector of the values formatted to \code{decimals} significant digits.
#' @examples
#' f <- fmt_dcimals(2)
#' f(c(1.23456, 200.567, 0.001234))
#' @export
fmt_dcimals <- function(decimals = 0){
  # return a function responpsible for formatting the
  # axis labels with a given number of decimals
  function(x) as.character(signif(x, decimals))
}


#' Compute confidence intervals quietly
#'
#' Some `confint()` methods emit progress messages (for example,
#' "Waiting for profiling to be done...") when computing confidence intervals.
#' `quiet_confint()` returns the usual `confint()` result while suppressing
#' warnings and muffling messages produced during the call.
#'
#' @param fit A fitted model object with a `confint()` method.
#' @param ... Additional arguments passed to \code{stats::confint()}.
#'
#' @return The result of `stats::confint(fit)` (typically a matrix of confidence
#'   intervals).
#'
#' @keywords internal
#' @noRd
quiet_confint <- function(fit, ...) {
  suppressMessages(suppressWarnings(
    withCallingHandlers(
      stats::confint(fit, ...),
      message = function(m) invokeRestart("muffleMessage")
    )
  ))
}


#' Apply median scaling to selected columns of a data frame
#'
#' Divides each numeric column in \code{cols} by its median, so that 1.0
#' represents the median value. Non-numeric columns in \code{cols} are silently
#' skipped. Returns updated data, the per-column scaling factors, and a
#' \code{naming_df} suitable for downstream back-transformation in
#' \code{make_forest_plot()}.
#'
#' @param data data.frame. Input data (rows = observations).
#' @param cols character. Names of columns to scale (typically
#'   \code{c(covariates, exposure_vars)}).
#' @param units_df Optional data.frame with columns \code{covname} and
#'   \code{units} to attach unit labels. Joined onto \code{naming_df} when
#'   provided.
#' @param is_centered Logical (default \code{FALSE}). If \code{TRUE}, the
#'   \code{transform_type} column in \code{naming_df} is set to
#'   \code{"median_scaled_and_centered"} instead of \code{"median_scaled"}.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{\code{data}}{The modified data frame with scaled columns.}
#'   \item{\code{scaling_factors}}{Named numeric vector of medians used to
#'     scale each column.}
#'   \item{\code{naming_df}}{A tibble with columns \code{covname},
#'     \code{median_point}, \code{transform_type}, and \code{units}, or
#'     \code{NULL} if no numeric columns were found in \code{cols}.}
#' }
#'
#' @keywords internal
apply_median_scaling <- function(data, cols, units_df = NULL, is_centered = FALSE) {
  numeric_mask <- vapply(data[cols], is.numeric, logical(1))
  cols_to_scale <- cols[numeric_mask]

  if (length(cols_to_scale) == 0) {
    return(list(data = data, scaling_factors = numeric(0), naming_df = NULL))
  }

  scaling_factors <- vapply(data[cols_to_scale], stats::median, numeric(1), na.rm = TRUE)

  data <- data %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(cols_to_scale), ~ . / stats::median(., na.rm = TRUE)))

  ttype <- if (is_centered) "median_scaled_and_centered" else "median_scaled"

  naming_df <- tibble::tibble(
    covname       = cols_to_scale,
    median_point  = unname(scaling_factors),
    transform_type = ttype
  )

  if (!is.null(units_df)) {
    naming_df <- dplyr::left_join(naming_df, units_df, by = "covname")
  }
  if (!"units" %in% colnames(naming_df)) {
    naming_df <- naming_df %>% dplyr::mutate(units = "U")
  }

  list(data = data, scaling_factors = scaling_factors, naming_df = naming_df)
}
#' Build a model formula from a response string and RHS terms
#'
#' Constructs a \code{\link[stats]{formula}} from character strings, avoiding
#' \code{eval(parse())} patterns. The right-hand side may be a single
#' pre-assembled string (e.g., \code{"AGE + WT"}) or a character vector whose
#' elements are joined with \code{" + "}.
#'
#' @param response Character string for the left-hand side
#'   (e.g., \code{"DV"} or \code{"survival::Surv(TIME, DV)"}).
#' @param rhs Character string or vector for the right-hand side terms.
#'   A length-1 string is used as-is; a vector is collapsed with \code{" + "}.
#'
#' @return A \code{\link[stats]{formula}} object.
#' @keywords internal
#' @noRd
build_model_formula <- function(response, rhs) {
  rhs_str <- if (length(rhs) == 1L) rhs else paste(rhs, collapse = " + ")
  stats::as.formula(paste(response, "~", rhs_str))
}
