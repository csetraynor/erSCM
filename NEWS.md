# erSCM 0.1.4

## API Changes

* Parameters that were silently ignored since the initial commit now raise errors
  via `.NotYetImplemented()` so callers get an explicit signal rather than silent
  misconfiguration:
  - `removable_terms` in `scm_backward_delete()` / `scm_backward_delete_coxph()`:
    changed from a required argument to `= NULL`; non-NULL raises an error.
  - `fixed_terms` in `scm_backward_delete()` / `scm_backward_delete_coxph()`:
    non-NULL raises an error.
  - `criterion` in `scm_backward_delete()`, `scm_backward_delete_coxph()`,
    `scm_forward_add()`, `scm_forward_add_coxph()`,
    `scm_forward_first_step()`, `scm_forward_first_step_coxph()`:
    any value other than `"AIC"` raises an error.
  - `min_terms` in `scm_backward_delete()` / `scm_backward_delete_coxph()`:
    any positive value raises an error.
  - `max_terms` in `scm_forward_add()` / `scm_forward_add_coxph()`:
    any finite value raises an error.
  - `parallel` in `scm_coxph()`, `fit_univariate_models()`,
    `fit_univariate_models_coxph()`: `TRUE` (or any non-`FALSE`) raises an error.

* `outcome` in `scm_backward_delete()` / `scm_backward_delete_coxph()`: changed
  from a required argument to `= NULL`; silently ignored. Documented as
  deprecated; will be removed in 0.2.0.

* `run_tag` in all step functions: no behaviour change; `@param` documentation
  updated to "Reserved for future audit logging; currently unused."

## New Features

* New function `vpc_coxph_sim()`: a simulation-based predictive Cox VPC for
  `scm_coxph()` fits. It draws coefficients from
  `N(beta_hat, vcov(beta_hat))`, simulates event times by inverse sampling from
  a fixed Breslow baseline hazard, applies fixed observed follow-up censoring,
  and summarizes per-stratum Kaplan-Meier curves on a common time grid (default
  unique observed event times). The output overlays observed KM curves against
  simulated pointwise prediction bands and supports seed-parity reproducibility
  between sequential and parallel execution.

* New function `vpc_logistic()`: a simulation-based visual predictive check for
  binary logistic SCM fits from `scm_glm()`, complementing `vpc_coxph()`. It
  draws coefficient vectors from the fitted sampling distribution (parametric
  bootstrap via the Cholesky factor of `vcov()`), simulates outcomes on the
  observed design, and overlays the resulting per-bin prediction band against
  the observed event proportion (with Wilson intervals) across covariate bins.
  Unlike `plot_logistic_terms()` (a covariate-effect plot at a reference
  profile), `vpc_logistic()` assesses model calibration. It is a seed-parity
  sibling of `scm_bootstrap()`: for a given `seed`, results are identical whether
  `parallel = FALSE` or `parallel = TRUE`, and the caller's RNG state is
  preserved. The per-replicate L'Ecuyer-CMRG stream generation is now shared
  between `scm_bootstrap()` and `vpc_logistic()` via an internal helper
  (behaviour of `scm_bootstrap()` is unchanged).

## Bug Fixes

* `vpc_coxph()`: added strict input validation (`eres_scm_coxph` class,
  non-NULL `final_model`, `ci` in (0, 1), integer `n_bins >= 1`), removed
  grouped-summary message noise, enforced stable ordered-facet strata display,
  and corrected documentation to describe the current analytical
  confidence-band diagnostic behavior.
* User-visible change: `vpc_coxph()` now errors explicitly when `n_bins > 1`
  and the fitted main effects include numeric terms. This path was previously
  dead/ambiguous; real numeric binning is deferred to
  `reports/vpc_coxph_v2_plan.md`.

* Bug fix (critical): `plot_cox_terms()` previously assigned predicted survival
  curves to strata by position, causing silent mislabelling when bin labels did
  not sort in numeric order (e.g., with var_spec back-transformed labels or
  n_bins >= 10). Curves are now matched by bin/level identity.
* `scm_glm()` now returns a proper `eres_scm_glm` S3 object even when no
  exposure variable reaches the entry significance threshold. Previously it
  returned a plain list, which caused downstream functions (`print`,
  `plot_nonlinear_effect`, `scm_bootstrap`) to fail or produce degraded output.
  `scm_coxph()` already had this behavior. The returned object has
  `final_model = NULL`; `print()` now describes the no-model state clearly.
* `plot_nonlinear_effect()` now produces a clear error ("Cannot plot: no model
  was fit.") when passed a fit where no exposure was selected.
* `scm_bootstrap()` now produces a clear error ("Cannot bootstrap: the original
  fit did not produce a model.") rather than a cryptic downstream failure.
* `scm_bootstrap()`: when multiple exposure variables have non-linear
  transformations enabled, transformation selection frequencies are now
  correctly tracked for all candidates. Previously only the first exposure's
  transforms were recorded; all subsequent NL-enabled exposures were silently
  dropped from `exposure_transform_freq`.
* `plot_nonlinear_effect()`: removed unused intermediate variables in the
  response-scale GLM branch. The y-axis correctly shows OR (GLM) or HR
  (CoxPH) relative to the median reference value.
* Bug fix (reproducibility): `scm_bootstrap()` now produces identical resamples
  for the same seed regardless of the `parallel` argument. Previously,
  sequential and parallel paths seeded RNG independently and diverged. The
  internal implementation now derives per-replicate L'Ecuyer-CMRG streams from
  `seed` and applies them uniformly to both paths.
* Note: results produced by `scm_bootstrap(seed = ...)` in prior versions will
  not match this version. To reproduce a prior analysis, re-run under this
  version and treat its output as the new reference.
* Fix (message noise): `plot_cox_terms()` no longer emits dplyr join and
  `summarise()` grouping messages. All internal joins now specify keys
  explicitly and all summaries specify `.groups`.
* Fix (validation): `.fp_term()` and `.get_back_transform_fn()` now abort with
  an informative message on invalid inputs (previously they returned silently
  invalid formulas or fell back to identity).

## Improvements

* `plot_cox_terms()` now renders facets in numeric/factor-level order of the
  covariate bins rather than lexical order of the label strings.
* `run_scm()` now accepts `nonlinear`, `linearity_alpha`, `fp_powers`, and
  `spline_df` as explicit named parameters with validation (previously forwarded
  silently via `...`). Invalid values now produce informative errors via
  `cli::cli_abort()`.
* `plot_nonlinear_effect()` and `plot_all_nonlinear()` are now correctly exported
  and accessible after `library(erSCM)` (previously had `@export` tags but were
  absent from `NAMESPACE`).
* `scm_glm()` and `scm_coxph()` now emit a `cli::cli_warn()` when both `nonlinear`
  and `interactions` are specified, documenting the known limitation that
  interactions whose parent variable enters in a non-linear form are silently
  skipped.
* New vignette: *Non-Linear Transformations in SCM*
  (`vignette("nonlinear-transformations")`).
* `scm_bootstrap()`'s `exposure_transform_freq` return field is now a named list
  (one entry per NL-enabled exposure variable) for consistent handling across
  single- and multi-exposure cases.
* Expanded test coverage for non-linear transformations combined with
  `adjust_for` and `start_covariates` (backward direction).
* Packaging: raised the minimum ggplot2 version to 3.4.0 to reflect actual code
  usage (linewidth aesthetic).

## Internal

* Cleaned up `utils::globalVariables()`; function names that were previously
  masked (`median`, `desc`, `order`, `all_of`) are now qualified at their call
  sites.
* Interaction-forward LRT p-value extraction is now robust to `anova`
  column-name variants (`Pr(>Chi)` and `Pr(>|Chi|)`).

## Testing

* Bootstrap test suite optimised: `n_boot` values reduced to 10–20 across all
  tests. Long-running stability assessments (n_boot ≥ 200) are demonstrated in
  the vignette rather than the test suite.
* All bootstrap tests are guarded with `skip_on_cran()`. The non-CRAN portion of
  the suite completes in under 3 minutes on a typical CPU.

## Known Limitations

* `vpc_coxph_sim()` is categorical-first in this release. Numeric
  `covariate_term` with `n_bins > 1` aborts explicitly and points to follow-up
  phase S-6 for numeric binning support.
* `vpc_coxph_sim()` aborts on time-dependent Cox model forms (`tt()`,
  `strata()`, `cluster()`, and start-stop `Surv()` responses). Support for
  these forms is deferred to a v3 design pass.

* The backward simplification for non-linear terms currently tests
  FP2 → linear directly, skipping the FP1 intermediate step from the full
  Royston & Altman MFP hierarchy. This is a deliberate simplification adequate
  for typical pharmacometric datasets.
* Interactions between a non-linearly-entered variable and other covariates are
  silently skipped; users should not specify interactions involving non-linearly
  entered variables. A warning is emitted in this case.

# erSCM 0.1.3

## New Features

* `scm_bootstrap()` now includes a two-layer safeguard for **categorical
  covariates with rare factor levels**, which can cause complete separation or
  convergence failures during bootstrap resampling.
  - New parameters `min_level_n = 5L` and `min_events_per_level = 1L`: any
    factor covariate with a level below these thresholds is excluded from that
    replicate's candidate pool (counted as "not selected" — conservative).
  - Pre-bootstrap diagnostic (`diagnostics` return field): scans the full
    dataset for problematic factor levels and emits a `cli_warn` (high risk) or
    `cli_inform` (moderate risk) before the loop begins.
  - New return fields: `convergence_issues` (tibble of failed replicates),
    `variable_exclusion_rate` (named numeric vector of per-variable exclusion
    proportions), and `diagnostics` (pre-bootstrap suitability report).
  - `print.scm_bootstrap()` now shows a "Factor Level Warnings" section when
    any variable has a >10% exclusion rate.
  - Set both `min_level_n = 0L` and `min_events_per_level = 0L` to disable all
    new safeguards and restore pre-Phase-6 behavior.

# erSCM 0.1.2

## New Features

* Added **Exposure Selection Stability** tracking to `scm_bootstrap()`. When
  multiple candidate exposure variables are provided, the bootstrap now records
  which exposure is selected in each replicate. New return fields:
  - `exposure_selection_prob`: proportion of converged replicates each exposure
    was selected (includes a `none_significant` entry when applicable).
  - `exposure_selected_per_boot`: character vector (length `n_converged`) of the
    exposure selected in each converged replicate.
  - `inclusion_prob_by_exposure`: named list of covariate inclusion probabilities
    conditional on each exposure being selected; `NULL` for single-exposure fits.
  - `plot(stability, type = "exposure")` for a dedicated selection bar chart.
  - `stable_covariates(stability, include_exposure = TRUE)` to include exposures
    in the stable-term list.
  - `summary()` now returns an `exposure_stability` tibble with columns
    `exposure_var`, `selection_prob`, and `is_original`.
  - New vignette section: *Exposure Selection Stability* (section 7).
* `modal_model` now returns a `list(exposure, covariates)` instead of a plain
  character vector, capturing the full modal model identity.
* `original_selected` now returns a `list(exposure, covariates)` instead of a
  plain character vector.
* Fixed a metadata bug in `scm_coxph()`: `metadata$candidate_exposure_vars`
  now correctly stores the **full list** of candidate exposures (previously only
  stored the selected winner). Both `scm_glm()` and `scm_coxph()` now include
  `metadata$selected_exposure` recording the winning exposure.

# erSCM 0.1.1

## New Features

* Added `scm_bootstrap()` for Bootstrap Model Selection Stability (BMS) assessment.
  Accepts a fitted `scm_glm()` or `scm_coxph()` object and resamples the dataset
  `n_boot` times, running the full SCM procedure on each replicate to quantify
  covariate selection stability.
  - `print()`, `summary()`, and `plot()` S3 methods for `"scm_bootstrap"` objects.
  - `stable_covariates()` utility to extract covariates exceeding a stability
    threshold.
  - Patient-level resampling via `subject_col` for repeated-measures or
    multi-row-per-subject data.
  - Optional parallel execution via `future.apply` (set `parallel = TRUE` with
    an active `future::plan()`).
  - Co-inclusion matrix to identify competing or complementary covariate pairs.
  - Model variability index (`1 - modal_model_frequency`) for overall stability
    assessment.
  - New vignette: *Bootstrap Model Selection Stability*
    (`vignette("bootstrap-stability")`).

* Added `run_scm()` unified dispatcher — routes to `scm_glm()` or `scm_coxph()`
  based on the `type` argument, providing a single entry point for programmatic
  use such as bootstrap resampling engines.

# erSCM 0.1.0

* Initial CRAN submission.
* Stepwise Covariate Modeling for GLM (`scm_glm()`) and Cox PH (`scm_coxph()`).
* S3 classes `eres_scm_glm` and `eres_scm_coxph` with `print()` and `plot()` methods.
* `adjust_for` parameter in `scm_glm()` and `scm_coxph()` for forced covariate
  adjustment: terms are included in the base model before forward selection and
  protected from backward elimination.
* Visual distinction between forced-adjustment and SCM-selected covariates in
  forest plots (`make_forest_plot()`) and summary tables (`get_sum_tbl()`).
* GLM covariate effect plots (`plot_glm_terms()`, `plot_glm_continuous()`,
  `plot_glm_categorical()`).
* Variable specification framework via `var_spec()` for display labels, units,
  and back-transformation of scaled covariates.
* Interaction term support in forward selection (`expand_interactions()`).
* Forest plot visualisation (`make_forest_plot()`, `plot_grouped_ratios()`).
* Visual Predictive Check for Cox models (`vpc_coxph()`, `plot_cox_terms()`).
* Summary table utilities (`get_sum_tbl()`, `get_sum_tbl_coxph()`).
* Univariate screening (`fit_univariate_models()`, `fit_univariate_models_coxph()`).
* Internal formula construction now uses `stats::as.formula()` throughout;
  all `eval(parse(...))` patterns have been replaced.
* Demo dataset `demo_eres` included.
