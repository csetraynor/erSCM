#' Tidy summary table for glm models
#'
#' @description
#' Create a tidy summary table of parameter estimates and model fit criteria
#' from a fitted \code{glm} object. The output combines coefficient-level
#' information (estimate, standard error, test statistic, p-value) with
#' model-level metrics (AIC, BIC, deviance, null deviance, residual degrees of
#' freedom, number of observations, and pseudo-\eqn{R^2} where appropriate).
#'
#' @param model A fitted \code{glm} object.
#' @param p_min Numeric giving the minimum p-value display threshold.
#'   Defaults to \code{0.001}.
#' @param include_conf_int Logical; if \code{TRUE}, adds lower and upper
#'   confidence interval columns for coefficients. Defaults to \code{FALSE}.
#' @param exponentiate Logical: if exponentiated estimate should be included.
#' @param conf_level Confidence level for intervals (default \code{0.95}).
#' @param conf_type Character string indicating the type of interval.
#'   One of \code{"wald"} or \code{"profile"}. Defaults to \code{"wald"}.
#' @param adjust_for Optional character vector of covariate names that were
#'   forced into the model (i.e., passed as \code{adjust_for} in
#'   \code{\link{scm_glm}}). When supplied, a \code{term_type} column is added
#'   to the output classifying each row as \code{"adjusted"} (forced covariate),
#'   \code{"selected"} (SCM-selected main effect), or \code{"base"}
#'   (intercept / other model terms). Defaults to \code{NULL}.
#' @param show_terms Character string controlling which rows to return.
#'   One of \code{"all"} (default), \code{"selected"} (only SCM-selected
#'   terms), or \code{"adjusted"} (only forced-adjustment terms). Has no effect
#'   when \code{adjust_for = NULL}.
#'
#' @details
#' The coefficient table includes:
#' \itemize{
#'   \item \code{term}: Parameter name
#'   \item \code{estimate}: Coefficient estimate
#'   \item \code{std_error}: Standard error
#'   \item \code{exp(estimate)}: Exponentiated estimate (only when \code{exponentiate = TRUE})
#'   \item \code{statistic}: Wald z-statistic
#'   \item \code{p-value}: Two-sided p-value (values below \code{p_min} are displayed
#'     as \code{"<p_min"})
#' }
#'
#' @return A \code{data.frame} with one row per coefficient, containing the
#' columns above (numeric columns coerced to character for display). When
#' \code{adjust_for} is supplied, a \code{term_type} column is added
#' (\code{"adjusted"}, \code{"selected"}, or \code{"base"}) and rows are
#' reordered accordingly; \code{show_terms} can then subset to only one type.
#'
#' @examples
#' data(demo_eres)
#' m <- glm(sim_status ~ AGE + WT, data = demo_eres, family = binomial())
#' get_sum_tbl(m)
#'
#' @importFrom tidyselect where
#' @export
get_sum_tbl <- function(model, p_min = 0.001, include_conf_int = FALSE, exponentiate = FALSE,
                        conf_level = 0.95, conf_type = "wald",
                        adjust_for = NULL,
                        show_terms = c("all", "selected", "adjusted")) {

  show_terms <- match.arg(show_terms)

  if (!inherits(model, "glm")) {
    cli::cli_abort("{.arg model} must be a fitted {.cls glm} object.")
  }

  x <- model

  summary_table <- summary(x)$coefficients
  if(exponentiate){
    result <- data.frame(
      Estimate = summary_table[, 1],
      StdError = summary_table[, 2],
      expEstimate = exp(summary_table[, 1]),
      ZValue = summary_table[, 3],
      Pr = summary_table[, 4]
    )

    p_min_char <- paste0("<", p_min)

    # result <- result %>% slice(-1)
    result <- result   %>% as.data.frame() %>%
      tibble::rownames_to_column(var = "strata") %>%
      dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ as.character(round(., 6)))) %>%
      dplyr::mutate(dplyr::across("Pr", ~ ifelse(. < p_min, p_min_char, as.character(.))))
    colnames(result) <- c("term", "estimate", "std_error", "exp(estimate)", "statistic", "p-value")
  } else {
    result <- data.frame(
      Estimate = summary_table[, 1],
      StdError = summary_table[, 2],
      ZValue = summary_table[, 3],
      Pr = summary_table[, 4]
    )
    p_min_char <- paste0("<", p_min)
    # result <- result %>% slice(-1)
    result <- result   %>% as.data.frame() %>%
      tibble::rownames_to_column(var = "strata") %>%
      dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ as.character(round(., 6)))) %>%
      dplyr::mutate(dplyr::across("Pr", ~ ifelse(. < p_min, p_min_char, as.character(.))))
    colnames(result) <- c("term", "estimate", "std_error", "statistic", "p-value")
  }

  if (!is.null(adjust_for)) {
    result <- result %>%
      dplyr::mutate(term_type = dplyr::case_when(
        term == "(Intercept)" ~ "base",
        vapply(term, function(t) {
          any(vapply(adjust_for, function(a) startsWith(t, a), logical(1)))
        }, logical(1)) ~ "adjusted",
        TRUE ~ "selected"
      ))

    # Logical ordering: adjusted first, then selected, then base
    type_order <- c("adjusted", "selected", "base")
    result <- result[order(match(result$term_type, type_order)), ]

    if (show_terms != "all") {
      result <- result[result$term_type == show_terms, ]
    }
  }

  return(result)
}

#' Fit glm and return tidy summary table
#'
#' @description
#' Fits a \code{glm} on a dataset using a user-specified outcome, covariates, and family,
#' optionally standardizing specified continuous covariates by their median, and returns
#' a tidy summary table via \code{get_sum_tbl()}.
#'
#' @param data A \code{data.frame} containing the analysis variables.
#' @param covariates Character vector; names of covariates to include in the model.
#' @param family A GLM family object (e.g., \code{binomial()}, \code{gaussian()}, \code{poisson()}).
#' @param standardize Logical; if \code{TRUE}, standardizes variables listed in
#'   \code{continuous_covariates} by dividing by their median (after centering at 0 if median is used).
#'   Defaults to \code{FALSE}.
#' @param continuous_covariates Character vector; names of continuous covariates to standardize
#'   when \code{standardize = TRUE}. Must be a subset of \code{covariates}. Defaults to \code{NULL}.
#' @param ... Additional arguments passed to \code{glm()} (e.g., \code{weights}, \code{subset}).
#'
#' @details
#' The model formula is constructed as \code{DV ~ covariates}. When \code{standardize = TRUE},
#' each variable in \code{continuous_covariates} is transformed as \code{x / median(x, na.rm = TRUE)}
#' (observations with zero or missing median are left unchanged).
#'
#' @return A \code{data.frame} as produced by \code{get_sum_tbl()}, containing coefficient-level
#' statistics and model-level metrics for the fitted \code{glm}.
#'
#' @examples
#' # Binary outcome example
#' set.seed(123)
#' d <- data.frame(DV = rbinom(200, 1, 0.4),
#'                 x1 = rnorm(200),
#'                 x2 = rnorm(200),
#'                 g  = factor(sample(letters[1:3], 200, TRUE)))
#' run_sum_tbl(d, covariates = c("x1", "x2", "g"),
#'             family = binomial())
#'
#' # With standardization of continuous covariates by median
#' run_sum_tbl(d, covariates = c("x1", "x2", "g"),
#'             family = binomial(), standardize = TRUE,
#'             continuous_covariates = c("x1", "x2"))
#'
#' @seealso \code{\link[stats]{glm}}, \code{\link{get_sum_tbl}}
#' @export
run_sum_tbl <- function(data, covariates, family, standardize = FALSE, continuous_covariates, ...) {

  df <- data

  if(standardize){
    df <- df %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(continuous_covariates[continuous_covariates %in% covariates]), ~ . / median(.)))
  }
  formul <- as.formula(paste("DV ~", paste0(covariates, collapse = " + ") ))

  model <- glm(formul, data = df  , family = family)

  return(model)
}


#' Tidy summary table for coxph models
#'
#' @description
#' Create a tidy summary table of parameter estimates and model fit criteria
#' from a fitted \code{coxph} object. The output combines coefficient-level
#' information (estimate, standard error, test statistic, p-value) with
#' model-level metrics (AIC, BIC, deviance, null deviance, residual degrees of
#' freedom, number of observations, and pseudo-\eqn{R^2} where appropriate).
#'
#' @param model A fitted \code{coxph} object.
#' @param p_min Numeric giving the minimum p-value display threshold.
#'   Defaults to \code{0.001}.
#' @param include_conf_int Logical; if \code{TRUE}, adds lower and upper
#'   confidence interval columns for coefficients. Defaults to \code{FALSE}.
#' @param exponentiate Logical: if exponentiated estimate should be included.
#' @param conf_level Numeric value in \code{(0, 1)} giving the confidence level
#'   for coefficient intervals. Defaults to \code{0.95}.
#'
#' @details
#' The coefficient table includes \code{term}, \code{hazard Ratio} (the
#' coefficient estimate, exponentiated if \code{exponentiate = TRUE}),
#' \code{std.error}, \code{lower CI}/\code{upper CI}, \code{statistic}, and
#' \code{p-value} (values below \code{p_min} are displayed as \code{"< p_min"}).
#'
#' Note: \code{include_conf_int} must currently be set to \code{TRUE} -- the
#' confidence interval columns are always selected internally, so the default
#' \code{FALSE} will error.
#'
#' @return A \code{data.frame} (tibble) with one row per coefficient, containing
#' the columns described above.
#'
#' @examples
#' data(demo_eres)
#' m <- survival::coxph(survival::Surv(sim_time, sim_status) ~ AGE + WT, data = demo_eres)
#' get_sum_tbl_coxph(m, include_conf_int = TRUE)
#'
#' @export
get_sum_tbl_coxph <- function(model, p_min = 0.001, include_conf_int = FALSE, exponentiate = FALSE, conf_level = 0.95) {

  if (!inherits(model, "coxph")) {
    cli::cli_abort("{.arg model} must be a fitted {.cls coxph} object.")
  }

  tidy_cox <- broom::tidy(model, exponentiate = exponentiate, conf.int = include_conf_int, conf.level = conf_level)

  summary_table <- tidy_cox %>%
    dplyr::select(term, estimate, std.error, conf.low, conf.high, statistic, p.value) %>%
    dplyr::mutate(
      term = stringr::str_replace_all(term, "_", " "), # Replaces underscores with spaces in variable names
      estimate = round(estimate, 3),
      std.error = round(std.error, 3),
      conf.low = round(conf.low, 3),
      conf.high = round(conf.high, 3),
      statistic  = round(statistic, 3),
      p.value = ifelse(p.value < 0.005, "< 0.005", as.character( round(p.value, 3)) )  #round the p-values
    ) %>%
    dplyr::rename(
      "hazard Ratio" = estimate,
      "lower CI" = conf.low,
      "upper CI" = conf.high,
      "p-value" = p.value
    )

  #print(result)
  return(summary_table)
}

#' Fit glm and return tidy summary table
#'
#' @description
#' Fits a \code{glm} on a dataset using a user-specified outcome, covariates, and family,
#' optionally standardizing specified continuous covariates by their median, and returns
#' a tidy summary table via \code{get_sum_tbl()}.
#'
#' @param data A \code{data.frame} containing the analysis variables.
#' @param covariates Character vector; names of covariates to include in the model.
#' @param standardize Logical; if \code{TRUE}, standardizes variables listed in
#'   \code{continuous_covariates} by dividing by their median.
#'   Defaults to \code{FALSE}.
#' @param continuous_covariates Character vector; names of continuous covariates to standardize
#'   when \code{standardize = TRUE}. Must be a subset of \code{covariates}. Defaults to \code{NULL}.
#' @param ... Additional arguments passed to \code{survival::coxph()}.
#'
#' @details
#' The model formula is constructed as \code{Surv(TIME, DV) ~ covariates}. When
#' \code{standardize = TRUE}, each variable in \code{continuous_covariates} is
#' transformed as \code{x / median(x, na.rm = TRUE)}.
#'
#' @return A fitted \code{coxph} object.
#'
#' @examples
#' data(demo_eres)
#' # NOTE: the outcome/time columns must already be literally called "DV"/"TIME".
#' demo_eres$DV <- demo_eres$sim_status
#' demo_eres$TIME <- demo_eres$sim_time
#' run_sum_tbl_coxph(demo_eres, covariates = c("AGE", "WT"))
#'
#' @seealso \code{\link[survival]{coxph}}, \code{\link{get_sum_tbl_coxph}}
#' @export
run_sum_tbl_coxph <- function(data, covariates, standardize = FALSE, continuous_covariates, ...) {

  df <- data

  if(standardize){
    df <- df %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(continuous_covariates[continuous_covariates %in% covariates]), ~ . / median(.)))
  }
  formul <- as.formula(paste("survival::Surv( TIME, DV) ~", paste0(covariates, collapse = " + ") ))

  model <- survival::coxph(formul, data = df)

  return(model)
}




#' Summarize and evaluate a Cox model equation
#'
#' @param fit A survival::coxph model object
#' @param newdata Optional data.frame to compute linear predictor (eta) and HR
#' @param digits Number of digits to print for coefficients
#' @param return_values If TRUE, returns a list with coefficients, lp, and hr
#'
#' @examples
#' data(demo_eres)
#' m <- survival::coxph(survival::Surv(sim_time, sim_status) ~ AGE + WT, data = demo_eres)
#' summarize_cox_equation(m)
#' summarize_cox_equation(m, newdata = demo_eres[1:3, ])
#'
#' @return If \code{return_values = FALSE}, called only for its side effect of
#'   printing the linear-predictor equation and hazard function; returns
#'   \code{NULL} invisibly. If \code{return_values = TRUE} (the default), a
#'   list with \code{coefficients} (named numeric vector), and, when
#'   \code{newdata} is supplied, additionally \code{lp} (linear predictor per
#'   row of \code{newdata}) and \code{hr} (\code{exp(lp)}, the hazard ratio
#'   relative to \code{lp = 0}).
#' @export
summarize_cox_equation <- function(fit, newdata = NULL, digits = 6, return_values = TRUE) {
  if (!inherits(fit, "coxph")) cli::cli_abort("{.arg fit} must be a {.cls coxph} object.")
  
  coefs <- stats::coef(fit)
  term_names <- names(coefs)
  
  # Build human-readable linear predictor string: eta = b1*x1 + b2*x2 + ...
  coef_str <- paste(sprintf("%+.*f*%s", digits, unname(coefs), term_names), collapse = " ")
  cat("Linear predictor (eta):\n")
  cat("eta = ", coef_str, "\n", sep = "")
  cat("Hazard function:\n")
  cat("h(t | X) = h0(t) * exp(eta)\n\n")
  
  out <- list(coefficients = coefs)
  
  # If newdata is provided, compute eta and HR safely using predict() to honor factor coding
  if (!is.null(newdata)) {
    # Validate required columns
    missing_cols <- setdiff(all.vars(stats::terms(fit)), colnames(newdata))
    if (length(missing_cols) > 0) {
      warning(sprintf("newdata is missing variables used in the model: %s",
                      paste(missing_cols, collapse = ", ")))
    }
    
    # Compute linear predictor (eta) and HR relative to baseline
    lp <- tryCatch({
      stats::predict(fit, newdata = newdata, type = "lp")
    }, error = function(e) {
      cli::cli_abort("predict() failed for provided {.arg newdata}: {e$message}")
    })
    hr <- exp(lp)
    
    # Print a small summary
    cat(sprintf("Computed linear predictor and hazard ratio for %d row(s):\n", length(lp)))
    print(data.frame(eta = lp, HR = hr), row.names = FALSE)
    
    out$lp <- lp
    out$hr <- hr
  }
  
  if (isTRUE(return_values)) return(out) else invisible(NULL)
}





