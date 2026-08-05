# ============================================================================ #
# NOTE: This script is retained for quick interactive testing.
# For documented workflows, see vignette("interactions", package = "erSCM")
# ============================================================================ #
#
# Example: Interaction Terms in erSCM
#
# This script demonstrates how to specify, fit, and visualize models that
# include interaction terms using the erSCM package. Interaction terms allow
# you to test whether the effect of one covariate depends on the value of
# another — a common question in pharmacometric exposure-safety analyses.
#
# Prerequisites: erSCM installed
# Runtime: ~60 seconds
# ============================================================================ #

library(erSCM)
library(dplyr)

# ============================================================================ #
# ---- Section 1: Data Preparation ----
# ============================================================================ #

data("demo_eres")

# Define covariates and exposure metrics
all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
all_exposure   <- c("CMAXss", "CMINss", "AUCss")

# Variable metadata — scaling and centering are handled by center = TRUE
my_spec <- c(
  CMAXss = "Cmax,ss (ng/mL)",
  CMINss = "Cmin,ss (ng/mL)",
  AUCss  = "AUC,ss (ng*h/mL)",
  WT     = "Body Weight (kg)",
  AGE    = "Age (years)",
  EGFR   = "eGFR (mL/min/1.73m2)",
  ALB    = "Albumin (g/dL)",
  BMI    = "BMI (kg/m2)",
  SEX    = "Sex",
  AAG    = "AAG (mg/dL)",
  AST    = "AST (U/L)",
  ALT    = "ALT (U/L)"
)

# Minimal prep — no manual scaling needed
demo_eres <- demo_eres |>
  mutate(SEX = as.factor(SEX), DVID = "ALL")

# ============================================================================ #
# ---- Section 2: Baseline Model (Main Effects Only) ----
# ============================================================================ #
# First, fit a standard SCM with no interactions as our reference point.

surv_base <- scm_coxph(
  eres_data          = demo_eres,
  covariates         = all_covariates,
  exposure_vars      = all_exposure,
  dv_col             = "sim_status",
  time_col           = "sim_time",
  significance_level = c(0.05, 0.01),
  dvid               = "ALL",
  dvid_col           = "DVID",
  var_spec           = my_spec,
  center             = TRUE
)

cat("\n=== Baseline model (main effects only) ===\n")
print(surv_base)

# ============================================================================ #
# ---- Section 3: Specifying Interaction Terms ----
# ============================================================================ #
#
# The `interactions` parameter accepts two formats:
#
# FORMAT 1: Character vector of explicit terms (most control)
#   interactions = c("CMAXss:AGE", "CMAXss:SEX", "WT:SEX")
#
# FORMAT 2: List of variable groups (generates all pairwise combinations)
#   interactions = list(c("CMAXss", "AGE", "WT", "SEX"))
#   -> Generates: CMAXss:AGE, CMAXss:WT, CMAXss:SEX, AGE:WT, AGE:SEX, WT:SEX
#
# IMPORTANT: The hierarchy/marginality constraint is enforced automatically.
# An interaction A:B is only tested if BOTH A and B are already in the model.
# This means not all candidate interactions will necessarily be evaluated.
#
# Scientific rationale for our candidate interactions:
#   - CMAXss:WT  -> Does the peak exposure effect on safety differ by body size?
#   - CMAXss:SEX -> Does the exposure-safety relationship differ between sexes?
#   - CMAXss:AGE -> Does peak exposure risk change with age?
#   - WT:SEX     -> Does the protective/risk effect of weight differ by sex?

# ============================================================================ #
# ---- Section 4: Cox PH with Interactions (after_main strategy) ----
# ============================================================================ #
#
# Strategy "after_main" (default): First completes the standard main-effect
# forward selection, then tests interactions in a separate forward loop.
# This is the recommended approach — it ensures the main-effect model has
# stabilized before testing interactions.

surv_interactions <- scm_coxph(
  eres_data            = demo_eres,
  covariates           = all_covariates,
  exposure_vars        = all_exposure,
  dv_col               = "sim_status",
  time_col             = "sim_time",
  significance_level   = c(0.05, 0.01),
  dvid                 = "ALL",
  dvid_col             = "DVID",
  var_spec             = my_spec,
  center               = TRUE,

  # --- Interaction parameters ---
  interactions         = list(c("CMAXss", "AGE", "WT", "SEX")),
  interaction_strategy = "after_main",
  p_interaction        = c(1, 0.5),
  max_interaction_order = 2,

  verbose              = TRUE
)

cat("\n=== Model with interactions (after_main) ===\n")
print(surv_interactions)

# ============================================================================ #
# ---- Section 5: Understanding the Results ----
# ============================================================================ #

# The final model terms (main effects + any significant interactions)
cat("\nFinal model terms:\n")
cat(paste(surv_interactions$metadata$finalRelation, collapse = "\n"), "\n")

# Check which terms are interactions vs main effects
main_effects <- get_main_effect_columns(surv_interactions$metadata$finalRelation)
interactions_in_model <- get_interaction_terms(surv_interactions$metadata$finalRelation)

cat("\nMain effects in final model:", paste(main_effects, collapse = ", "), "\n")
cat("Interactions in final model:", paste(interactions_in_model, collapse = ", "), "\n")

# View the selection path (shows both forward addition and backward steps)
cat("\nSelection path:\n")
print(surv_interactions$path)

# View coefficient table
cat("\nCoefficient summary:\n")
print(surv_interactions$summary)

# ============================================================================ #
# ---- Section 6: Interpreting Interaction Coefficients ----
# ============================================================================ #
#
# IMPORTANT: When an interaction A:B is in the model, the interpretation of
# the main-effect coefficients CHANGES:
#
#   - Without interaction: coef(A) = average effect of A across all B values
#   - With interaction:    coef(A) = effect of A when B = 0 (or reference level)
#                          coef(B) = effect of B when A = 0 (or reference level)
#                          coef(A:B) = how much the effect of A changes per
#                                      unit increase in B (and vice versa)
#
# For a Cox model: HR = exp(beta_A * A + beta_B * B + beta_AB * A * B)
#
# The hazard ratio for a 1-unit increase in A depends on the value of B:
#   HR_A(at B=b) = exp(beta_A + beta_AB * b)
#
# COMMON PITFALL: Do not interpret main effect p-values in isolation when
# the interaction is significant. The main effect is conditional on the
# other variable being at its reference value.

# ============================================================================ #
# ---- Section 7: Visualization — Forest Plot ----
# ============================================================================ #
#
# The forest plot shows:
#   - Main effects: odds/hazard ratios at data quantiles (5th, 25th, 50th, 75th, 95th)
#   - Interactions: exponentiated coefficients (exp(beta)) with 95% CI
#
# Interaction terms appear as a separate panel labeled with "A:B" notation.

cat("\n=== Forest Plot (with interactions) ===\n")
fp_interaction <- make_forest_plot(surv_interactions, plot = FALSE)
print(fp_interaction)

# Generate the plot
make_forest_plot(surv_interactions)

ggplot2::ggsave("plot_cmax_wt.png")

# Compare with baseline (no interactions)
# make_forest_plot(surv_base)

# ============================================================================ #
# ---- Section 8: Visualization — VPC ----
# ============================================================================ #
#
# The VPC compares observed Kaplan-Meier curves against model-predicted
# survival. Interaction effects are captured through the model predictions
# (survfit with newdata) — they don't need special handling.
#
# The VPC stratifies by categorical covariates in the final model.

cat("\n=== VPC (model with interactions) ===\n")
vpc_coxph(surv_interactions)

# ============================================================================ #
# ---- Section 9: Visualization — Covariate Term Plots ----
# ============================================================================ #
#
# plot_cox_terms() shows observed vs predicted survival stratified by a
# single covariate. You can use any MAIN EFFECT term — the interaction
# effects are still captured in the model predictions.
#
# NOTE: You cannot pass an interaction term (e.g., "CMAXss:AGE") as the
# covariate_term argument. Use the parent main effects individually.

cat("\n=== Covariate Term Plots ===\n")

# Plot exposure effect (stratified into quantile bins)
plot_cox_terms(surv_interactions, covariate_term = "CMAXss", n_bins = 4)

# Plot categorical covariate effect
plot_cox_terms(surv_interactions, covariate_term = "WT")

# This would error — interaction terms cannot be used as covariate_term:
# plot_cox_terms(surv_interactions, covariate_term = "CMAXss:AGE")

# ============================================================================ #
# ---- Section 10: GLM with Interactions ----
# ============================================================================ #
#
# The same interaction support is available for logistic regression via scm_glm().
# The syntax is identical.

lr_interactions <- scm_glm(
  eres_data            = demo_eres,
  covariates           = all_covariates,
  exposure_vars        = all_exposure,
  dv_col               = "sim_status",
  significance_level   = c(0.05, 0.01),
  var_spec             = my_spec,
  center               = TRUE,

  # --- Interaction parameters ---
  interactions         = c("CMAXss:WT", "CMAXss:SEX"),
  interaction_strategy = "after_main",
  p_interaction        = c(1, 1),

  verbose              = TRUE
)

cat("\n=== GLM with interactions ===\n")
print(lr_interactions)

# Forest plot for GLM (interactions appear as OR with CI)
make_forest_plot(lr_interactions)

# ============================================================================ #
# ---- Section 11: Simultaneous Strategy ----
# ============================================================================ #
#
# Strategy "simultaneous": Interactions compete with main effects at each
# forward step. An interaction can enter the model at any point once its
# parent main effects are present. This is more aggressive and may select
# interactions earlier, but can be useful for exploratory analyses.

surv_simultaneous <- scm_coxph(
  eres_data            = demo_eres,
  covariates           = all_covariates,
  exposure_vars        = all_exposure,
  dv_col               = "sim_status",
  time_col             = "sim_time",
  significance_level   = c(0.05, 0.01),
  dvid                 = "ALL",
  dvid_col             = "DVID",
  var_spec             = my_spec,
  center               = TRUE,
  interactions         = list(c("CMAXss", "AGE", "WT", "SEX")),
  interaction_strategy = "simultaneous",
  p_interaction        = c(0.05, 0.01),
  verbose              = TRUE
)

cat("\n=== Model with simultaneous strategy ===\n")
print(surv_simultaneous)

# Compare final terms between strategies
cat("\nafter_main final terms: ",
    paste(surv_interactions$metadata$finalRelation, collapse = ", "), "\n")
cat("simultaneous final terms: ",
    paste(surv_simultaneous$metadata$finalRelation, collapse = ", "), "\n")

# ============================================================================ #
# ---- Section 12: Custom p-value Thresholds ----
# ============================================================================ #
#
# Use `p_interaction` to apply stricter thresholds for interactions than for
# main effects. This is recommended because:
#   1. Interactions inflate the multiple testing burden
#   2. Overfitting risk is higher with many interaction candidates
#   3. Pharmacometric convention often uses stricter criteria for interactions
#
# Example: main effects at 0.05/0.01, interactions at 0.01/0.001

surv_strict <- scm_coxph(
  eres_data            = demo_eres,
  covariates           = all_covariates,
  exposure_vars        = all_exposure,
  dv_col               = "sim_status",
  time_col             = "sim_time",
  significance_level   = c(0.05, 0.01),
  dvid                 = "ALL",
  dvid_col             = "DVID",
  var_spec             = my_spec,
  center               = TRUE,
  interactions         = list(c("CMAXss", "AGE", "WT", "SEX")),
  interaction_strategy = "after_main",
  p_interaction        = c(0.01, 0.001)
)

cat("\nStrict interaction threshold — final terms:\n")
cat(paste(surv_strict$metadata$finalRelation, collapse = ", "), "\n")

# ============================================================================ #
# ---- Section 13: Backward Elimination with Interactions ----
# ============================================================================ #
#
# During backward elimination, the HIERARCHY CONSTRAINT is preserved:
#   - A main effect CANNOT be removed if an interaction involving it is still
#     in the model (it is "protected")
#   - Interaction terms can be removed like any other term
#   - Once an interaction is removed, its parent main effects become eligible
#     for removal in subsequent steps
#
# This ensures the model remains interpretable at all times.

# You can verify which main effects are protected:
example_relation <- c("CMAXss", "AGE", "WT", "SEX", "CMAXss:AGE", "WT:SEX")
protected <- get_protected_main_effects(example_relation)
cat("\nExample: terms in model:", paste(example_relation, collapse = ", "), "\n")
cat("Protected main effects:", paste(protected, collapse = ", "), "\n")
cat("Removable terms:", paste(example_relation[!example_relation %in% protected], collapse = ", "), "\n")

# ============================================================================ #
# ---- End of Example ----
# ============================================================================ #
cat("\n✓ Interaction example complete.\n")
