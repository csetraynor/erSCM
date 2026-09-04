#' Expand interaction specification to character vector of interaction terms
#'
#' Takes an interaction specification (character vector of "A:B" terms or a list
#' of variable groups) and returns a character vector of all interaction terms
#' in R formula notation.
#'
#' @param interactions Character vector or list specifying interaction terms.
#'   - If a character vector: each element is an interaction term in R formula notation
#'     using ":" (e.g., \code{c("age:sex", "AUC:weight")}). Returned as-is after validation.
#'   - If a list of character vectors: each vector defines a set of variables from which
#'     all pairwise (or higher-order) interactions are generated.
#'   - If NULL: returns \code{character(0)}.
#' @param max_interaction_order Integer. Maximum order of interactions to generate when
#'   \code{interactions} is a list. Default 2 (pairwise only).
#'
#' @return Character vector of interaction terms in "A:B" notation.
#'
#' @examples
#' expand_interactions(c("age:sex", "AUC:weight"))
#' expand_interactions(list(c("AUC", "age", "sex")))
#' expand_interactions(list(c("A", "B", "C")), max_interaction_order = 3)
#'
#' @export
expand_interactions <- function(interactions, max_interaction_order = 2) {

  if (is.null(interactions)) return(character(0))

  if (is.character(interactions)) {
    # Already in "A:B" format -- validate format
    for (term in interactions) {
      parts <- strsplit(term, ":")[[1]]
      if (length(parts) < 2) {
        cli::cli_abort(
          "Interaction term {.val {term}} must contain at least two variables separated by {.code :}."
        )
      }
    }
    return(interactions)
  }

  if (is.list(interactions)) {
    # Generate all combinations up to max_interaction_order
    all_interactions <- character(0)
    for (group in interactions) {
      if (!is.character(group) || length(group) < 2) {
        cli::cli_abort("Each element in {.arg interactions} must be a character vector of length >= 2.")
      }
      # Only orders 2..min(max_interaction_order, length(group)) are valid
      # interactions. Guard against max_interaction_order < 2, where the naive
      # `2:min(...)` would descend to c(2, 1) and emit spurious 1-way "terms".
      upper_order <- min(max_interaction_order, length(group))
      if (upper_order >= 2) {
        for (order in seq.int(2, upper_order)) {
          combos <- utils::combn(group, order, simplify = FALSE)
          terms <- vapply(combos, function(x) paste(x, collapse = ":"), character(1))
          all_interactions <- c(all_interactions, terms)
        }
      }
    }
    return(unique(all_interactions))
  }

  cli::cli_abort("{.arg interactions} must be NULL, a character vector, or a list of character vectors.")
}

#' Get main effects required by an interaction term
#'
#' Splits an interaction term (e.g., "A:B" or "A:B:C") into its constituent
#' main effect variable names.
#'
#' @param interaction_term Character scalar, e.g. "A:B" or "A:B:C".
#'
#' @return Character vector of main effect names.
#'
#' @examples
#' \dontrun{
#' get_parent_effects("age:sex")
#' get_parent_effects("A:B:C")
#' }
#'
#' @keywords internal
#' @noRd
get_parent_effects <- function(interaction_term) {
  strsplit(interaction_term, ":")[[1]]
}

#' Check which interaction terms are eligible given current model terms
#'
#' An interaction term is eligible for inclusion only if ALL of its parent
#' main effects are already present in the current model (hierarchy/marginality
#' constraint).
#'
#' @param candidate_interactions Character vector of interaction terms to check.
#' @param current_relation Character vector of terms currently in the model.
#'
#' @return Character vector of eligible interaction terms (subset of
#'   \code{candidate_interactions}).
#'
#' @examples
#' \dontrun{
#' get_eligible_interactions(
#'   c("age:sex", "AUC:age", "AUC:weight"),
#'   current_relation = c("AUC", "age", "sex")
#' )
#' # Returns: c("age:sex", "AUC:age")
#' }
#'
#' @keywords internal
#' @noRd
get_eligible_interactions <- function(candidate_interactions, current_relation) {
  if (length(candidate_interactions) == 0) return(character(0))


  # Only main effects in current_relation count (not other interactions)
  current_main_effects <- current_relation[!grepl(":", current_relation, fixed = TRUE)]

  eligible <- vapply(candidate_interactions, function(term) {
    parents <- get_parent_effects(term)
    all(parents %in% current_main_effects)
  }, logical(1))

  candidate_interactions[eligible]
}

#' Check which main effects are protected by interaction terms in the model
#'
#' A main effect is "protected" if removing it would violate the hierarchy
#' constraint -- i.e., there exists an interaction term in the current model
#' that involves this main effect.
#'
#' @param current_relation Character vector of all terms in the model
#'   (both main effects and interaction terms).
#'
#' @return Character vector of protected main effect names.
#'
#' @examples
#' get_protected_main_effects(c("AUC", "age", "sex", "AUC:age"))
#' # Returns: c("AUC", "age")
#'
#' @export
get_protected_main_effects <- function(current_relation) {
  interaction_terms <- current_relation[grepl(":", current_relation, fixed = TRUE)]
  if (length(interaction_terms) == 0) return(character(0))

  # All main effects that appear in any interaction

  protected <- unique(unlist(lapply(interaction_terms, get_parent_effects)))
  # Only those that are actually in the model as main effects
  intersect(protected, current_relation[!grepl(":", current_relation, fixed = TRUE)])
}

#' Validate interaction specification against available covariates
#'
#' Checks that all variables referenced in interaction terms exist in the
#' provided covariate and exposure variable sets. Stops with an informative
#' error if any variable is missing.
#'
#' @param interactions Character vector of expanded interaction terms (in "A:B" notation).
#' @param covariates Character vector of candidate covariates.
#' @param exposure_vars Character vector of exposure variables.
#'
#' @return \code{TRUE} invisibly if validation passes.
#'
#' @examples
#' \dontrun{
#' validate_interactions(
#'   c("age:sex", "AUC:weight"),
#'   covariates = c("age", "sex", "weight"),
#'   exposure_vars = c("AUC", "Cmax")
#' )
#' }
#'
#' @keywords internal
#' @noRd
validate_interactions <- function(interactions, covariates, exposure_vars) {
  if (length(interactions) == 0) return(invisible(TRUE))

  all_vars <- c(covariates, exposure_vars)

  for (term in interactions) {
    parents <- get_parent_effects(term)
    missing <- parents[!parents %in% all_vars]
    if (length(missing) > 0) {
      cli::cli_abort(c(
        "Interaction term {.val {term}} references variable(s) not in {.arg covariates} or {.arg exposure_vars}:",
        "x" = "{.val {missing}}"
      ))
    }
  }
  invisible(TRUE)
}


#' Extract main-effect column names from a Relation vector
#'
#' Filters a character vector of model terms (which may contain interaction terms
#' in "A:B" notation) to only those elements that are actual data column names --
#' i.e., terms that do NOT contain a colon.
#'
#' @param relation Character vector of model terms, potentially including
#'   interaction terms like "A:B".
#'
#' @return Character vector containing only the main-effect terms (no colons).
#'
#' @examples
#' get_main_effect_columns(c("AUC", "age", "sex", "AUC:age"))
#' # Returns: c("AUC", "age", "sex")
#'
#' @export
get_main_effect_columns <- function(relation) {
  relation[!grepl(":", relation, fixed = TRUE)]
}


#' Check if a model Relation vector contains interaction terms
#'
#' Returns \code{TRUE} if any element of the character vector contains a colon,
#' indicating the presence of interaction terms.
#'
#' @param relation Character vector of model terms.
#'
#' @return Logical scalar.
#'
#' @examples
#' has_interaction_terms(c("AGE", "WT"))
#' has_interaction_terms(c("AGE", "WT", "AGE:WT"))
#'
#' @export
has_interaction_terms <- function(relation) {
  any(grepl(":", relation, fixed = TRUE))
}


#' Extract interaction terms from a Relation vector
#'
#' Returns only the elements that contain ":" (interaction terms).
#'
#' @param relation Character vector of model terms.
#'
#' @return Character vector of interaction terms only.
#'
#' @examples
#' get_interaction_terms(c("AUC", "age", "sex", "AUC:age"))
#' # Returns: c("AUC:age")
#'
#' @export
get_interaction_terms <- function(relation) {
  relation[grepl(":", relation, fixed = TRUE)]
}


#' Compute conditional effect and SE for a variable participating in interactions
#'
#' Given a main-effect variable and one or more interaction terms involving it,
#' computes the conditional coefficient and its standard error (via the delta method)
#' at specified modifier values. This is used to evaluate the main-effect forest plot
#' row at clinically meaningful modifier values (e.g., median) rather than at modifier = 0.
#'
#' @param coefs Named coefficient vector from the fitted model (from \code{stats::coef()}).
#' @param V Variance-covariance matrix from the fitted model (from \code{stats::vcov()}).
#' @param effect_var Character scalar; name of the main-effect coefficient (e.g., \code{"CMAXss"}).
#' @param int_terms Character vector; names of the interaction coefficients involving
#'   \code{effect_var} (e.g., \code{c("CMAXss:WT")}). Must match coefficient names exactly.
#' @param modifier_values Numeric vector of the same length as \code{int_terms}; the value
#'   at which each interaction modifier is evaluated (e.g., the median of the modifier
#'   variable from the data).
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{beta}: the conditional coefficient value (numeric scalar).
#'   \item \code{se}: the standard error of the conditional coefficient (numeric scalar),
#'     derived via the generalized delta method: \eqn{SE = \sqrt{w' V w}} where
#'     \eqn{w = (1, m_1, m_2, \ldots)}.
#' }
#'
#' @examples
#' data(demo_eres)
#' fit <- glm(sim_status ~ CMAXss * WT, data = demo_eres, family = binomial())
#' compute_conditional_effect(
#'   coefs = coef(fit), V = vcov(fit),
#'   effect_var = "CMAXss", int_terms = "CMAXss:WT",
#'   modifier_values = median(demo_eres$WT)
#' )
#'
#' @export
compute_conditional_effect <- function(coefs, V, effect_var, int_terms, modifier_values) {
  cond_beta <- coefs[effect_var]
  for (i in seq_along(int_terms)) {
    cond_beta <- cond_beta + coefs[int_terms[i]] * modifier_values[i]
  }

  all_terms <- c(effect_var, int_terms)
  weights <- c(1, modifier_values)
  var_cond <- as.numeric(t(weights) %*% V[all_terms, all_terms] %*% weights)

  se_cond <- sqrt(max(var_cond, 0))
  list(beta = unname(cond_beta), se = se_cond)
}
