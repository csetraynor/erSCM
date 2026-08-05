# Nonlinear transformation engine
# Provides formula-fragment generators and model-fitting helpers used during
# the SCM forward step to evaluate non-linear functional forms.

#' Standard FP power set
#' @keywords internal
.FP_POWERS <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

#' Generate FP Term Formula Fragment
#'
#' Given a variable name and a power, produces the correct formula term.
#'
#' @param var_name Character.
#' @param power Numeric, the FP power.
#'
#' @return Character, formula fragment (e.g., \code{"I(CMAX^(-0.5))"} or
#'   \code{"log(CMAX)"}).
#' @keywords internal
.fp_term <- function(var_name, power) {
  if (power == 0) {
    paste0("log(", var_name, ")")
  } else if (power == 1) {
    var_name
  } else {
    paste0("I(", var_name, "^(", power, "))")
  }
}


#' Generate FP2 Pair Terms
#'
#' @param var_name Character.
#' @param p1 Numeric, first power.
#' @param p2 Numeric, second power (p2 >= p1).
#'
#' @return Character vector of formula terms (length 2).
#' @keywords internal
.fp2_terms <- function(var_name, p1, p2) {
  term1 <- .fp_term(var_name, p1)
  if (p1 == p2) {
    # Repeated power: second term is x^p * log(x), per Royston & Altman
    if (p1 == 0) {
      term2 <- paste0("I(log(", var_name, ")^2)")
    } else {
      term2 <- paste0("I(", var_name, "^(", p1, ") * log(", var_name, "))")
    }
  } else {
    term2 <- .fp_term(var_name, p2)
  }
  c(term1, term2)
}


# ---------------------------------------------------------------------------
# Uniform model-fitting helper
# ---------------------------------------------------------------------------

#' Fit a Model with Specified Formula Terms
#'
#' Utility to fit a GLM or CoxPH model given formula components.
#' Handles both model types uniformly.
#'
#' @param response_formula Character, the LHS
#'   (e.g., \code{"DV"} or \code{"survival::Surv(TIME, DV)"}).
#' @param rhs_terms Character, the RHS terms joined by \code{" + "}.
#' @param data data.frame.
#' @param model_type \code{"glm"} or \code{"coxph"}.
#' @param family Family for GLM (ignored for coxph).
#'
#' @return A fitted model object, or \code{NULL} if fitting fails.
#' @keywords internal
.fit_model <- function(response_formula, rhs_terms, data,
                       model_type = "glm", family = "binomial") {
  fml <- build_model_formula(response_formula, rhs_terms)
  tryCatch(
    {
      if (model_type == "coxph") {
        survival::coxph(fml, data = data)
      } else {
        stats::glm(fml, data = data, family = family)
      }
    },
    error   = function(e) NULL,
    warning = function(w) {
      # attempt to proceed; only return NULL on hard failures
      tryCatch(
        if (model_type == "coxph") {
          suppressWarnings(survival::coxph(fml, data = data))
        } else {
          suppressWarnings(stats::glm(fml, data = data, family = family))
        },
        error = function(e2) NULL
      )
    }
  )
}


#' Extract deviance (-2 log-likelihood) from a fitted model
#'
#' @param fit A fitted glm or coxph object, or NULL.
#' @return Numeric scalar, or Inf if fit is NULL or logLik fails.
#' @keywords internal
.model_deviance <- function(fit) {
  if (is.null(fit)) return(Inf)
  tryCatch(
    -2 * as.numeric(stats::logLik(fit)),
    error = function(e) Inf
  )
}


# ---------------------------------------------------------------------------
# Simple (fixed) transform helpers
# ---------------------------------------------------------------------------

#' Generate Simple (Fixed) Transformation Terms
#'
#' For transformations with a known, fixed formula (linear, log, sqrt),
#' returns the formula fragment immediately.
#'
#' @param var_name Character.
#' @param transform_type Character, one of \code{"linear"}, \code{"log"},
#'   \code{"sqrt"}.
#'
#' @return A list with: \code{formula_rhs} (character), \code{df} (integer),
#'   \code{label} (character).
#' @keywords internal
.make_simple_transform <- function(var_name, transform_type) {
  switch(transform_type,
    linear = list(
      formula_rhs = var_name,
      df          = 1L,
      label       = var_name
    ),
    log = list(
      formula_rhs = paste0("log(", var_name, ")"),
      df          = 1L,
      label       = paste0(var_name, " [log]")
    ),
    sqrt = list(
      formula_rhs = paste0("sqrt(", var_name, ")"),
      df          = 1L,
      label       = paste0(var_name, " [sqrt]")
    ),
    cli::cli_abort(
      "{.fn .make_simple_transform}: unrecognised transform {.val {transform_type}}."
    )
  )
}


# ---------------------------------------------------------------------------
# Spline helper
# ---------------------------------------------------------------------------

#' Generate Restricted Cubic Spline Terms
#'
#' @param var_name Character.
#' @param df Integer, degrees of freedom (default 3). This equals number of
#'   knots minus 1.
#' @param data data.frame (for computing knot positions at quantiles). Unused
#'   when only the formula fragment is needed, but kept for API consistency.
#'
#' @return A list with: \code{formula_rhs} (character), \code{df} (integer),
#'   \code{knots} (numeric vector or NULL), \code{label} (character).
#' @keywords internal
.make_spline_transform <- function(var_name, df = 3L, data = NULL) {
  df <- as.integer(df)

  # Compute knot positions from data quantiles if data is provided
  knots <- NULL
  if (!is.null(data) && var_name %in% colnames(data)) {
    x_obs <- data[[var_name]]
    x_obs <- x_obs[is.finite(x_obs)]
    if (length(x_obs) >= (df + 1L)) {
      probs <- seq(0, 1, length.out = df + 2L)[c(-1L, -(df + 2L))]
      knots <- unname(stats::quantile(x_obs, probs = probs))
    }
  }

  # Build formula fragment: splines::ns() with df argument
  formula_rhs <- paste0("splines::ns(", var_name, ", df = ", df, ")")
  label        <- paste0(var_name, " [RCS df=", df, "]")

  list(
    formula_rhs = formula_rhs,
    df          = df,
    knots       = knots,
    label       = label
  )
}


# ---------------------------------------------------------------------------
# FP search
# ---------------------------------------------------------------------------

#' Find Best Fractional Polynomial for a Variable
#'
#' Fits all FP models of the specified degree for a given variable within
#' the context of the current model, and returns the best one (by deviance).
#'
#' @param var_name Character, variable to transform.
#' @param degree Integer, 1 or 2. FP degree.
#' @param current_rhs Character, the current model's RHS (terms already in
#'   the model, not including the candidate variable). Use \code{"1"} for
#'   intercept-only baseline.
#' @param data data.frame.
#' @param response_lhs Character, LHS formula string (e.g., \code{"DV"} or
#'   \code{"survival::Surv(TIME, DV)"}).
#' @param model_type Character, \code{"glm"} or \code{"coxph"}.
#' @param family Family for GLM (ignored for coxph).
#' @param fp_powers Numeric vector of candidate powers.
#'
#' @return A list with:
#'   \describe{
#'     \item{formula_rhs}{Character, the formula terms to add.}
#'     \item{powers}{Numeric vector of selected powers.}
#'     \item{df}{Integer, degrees of freedom (1 for FP1, 2 for FP2).}
#'     \item{deviance}{Numeric, deviance of the best FP model.}
#'     \item{label}{Character.}
#'   }
#'   Returns \code{NULL} if no FP model converged.
#'
#' @details
#' For FP1 (degree 1): tests each power p in \code{fp_powers}. The model
#' term is: p = 0 \eqn{\Rightarrow} \code{log(x)}; p = 1 \eqn{\Rightarrow}
#' \code{x}; otherwise \code{I(x^p)}.
#'
#' For FP2 (degree 2): tests all pairs (p1, p2) with p1 <= p2. When
#' p1 == p2, the second term is \code{x^p * log(x)} per Royston & Altman.
#' Power 0 means \code{log(x)} in both cases.
#'
#' @keywords internal
.find_best_fp <- function(var_name, degree, current_rhs, data,
                          response_lhs = "DV",
                          model_type = "glm", family = "binomial",
                          fp_powers = .FP_POWERS) {
  degree <- as.integer(degree)
  if (!degree %in% 1L:2L) {
    cli::cli_abort("{.arg degree} must be 1 or 2, got {.val {degree}}.")
  }

  # Build list of candidate term sets
  if (degree == 1L) {
    candidates <- lapply(fp_powers, function(p) {
      list(terms = .fp_term(var_name, p), powers = p, df = 1L)
    })
  } else {
    # All pairs (p1, p2) with p1 <= p2
    pairs <- list()
    for (i in seq_along(fp_powers)) {
      for (j in seq_along(fp_powers)) {
        if (fp_powers[[j]] >= fp_powers[[i]]) {
          pairs <- c(pairs, list(list(p1 = fp_powers[[i]], p2 = fp_powers[[j]])))
        }
      }
    }
    candidates <- lapply(pairs, function(pr) {
      terms <- .fp2_terms(var_name, pr$p1, pr$p2)
      list(terms = terms, powers = c(pr$p1, pr$p2), df = 2L)
    })
  }

  # Base RHS (current model, no candidate variable)
  base_rhs <- if (nchar(trimws(current_rhs)) == 0L || current_rhs == "1") {
    "1"
  } else {
    current_rhs
  }

  best_dev   <- Inf
  best_cand  <- NULL

  for (cand in candidates) {
    rhs <- if (base_rhs == "1") {
      paste(cand$terms, collapse = " + ")
    } else {
      paste(c(base_rhs, cand$terms), collapse = " + ")
    }

    fit <- .fit_model(response_lhs, rhs, data, model_type, family)
    dev <- .model_deviance(fit)

    if (is.finite(dev) && dev < best_dev) {
      best_dev  <- dev
      best_cand <- cand
    }
  }

  if (is.null(best_cand)) return(NULL)

  powers_str <- paste(best_cand$powers, collapse = ", ")
  label <- paste0(var_name, " [FP", degree, ": p=(",  powers_str, ")]")

  list(
    formula_rhs = paste(best_cand$terms, collapse = " + "),
    terms       = best_cand$terms,
    powers      = best_cand$powers,
    df          = best_cand$df,
    deviance    = best_dev,
    label       = label
  )
}


# ---------------------------------------------------------------------------
# RCS degree extraction
# ---------------------------------------------------------------------------

#' Extract total model degrees of freedom from a fitted model
#'
#' @param fit A fitted glm or coxph object, or NULL.
#' @param model_type Character, \code{"glm"} or \code{"coxph"}.
#' @return Integer, total parameter df used by the model (excluding intercept for GLM).
#' @keywords internal
.model_df <- function(fit, model_type = "glm") {
  if (is.null(fit)) return(0L)
  if (model_type == "coxph") {
    tryCatch(
      as.integer(round(summary(fit)$logtest[["df"]])),
      error = function(e) as.integer(length(stats::coef(fit)))
    )
  } else {
    as.integer(fit$df.null - fit$df.residual)
  }
}


#' Extract degrees of freedom from an rcsN transform name
#'
#' @param transform_type Character, e.g. \code{"rcs"} or \code{"rcs4"}.
#' @param default_df Integer, fallback df when no number is embedded.
#' @return Integer df.
#' @keywords internal
.rcs_df <- function(transform_type, default_df = 3L) {
  if (transform_type == "rcs") return(as.integer(default_df))
  m <- regmatches(transform_type, regexpr("[0-9]+", transform_type))
  if (length(m) == 0L) return(as.integer(default_df))
  as.integer(m)
}
