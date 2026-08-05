# erSCM Phase B Planning Report

**Date:** 2026-08-03  
**Branch:** `feat-advance`  
**Scope:** Phase A verification + Phase B (Emax models) design  
**Status:** Planning only — no code changes

---

## 1. Executive Summary

Phase A (Non-Linear Transformations Within SCM) is **substantially complete**. All
seven specified source files exist, are correctly structured, and have passed the QC
refactoring pass. The transformation engine (`nonlinear_spec`, `nonlinear_transforms`,
`nonlinear_selection`), forward step, backward step, visualisation, and bootstrap
integration are all present, exported or correctly kept internal, and supported by
dedicated test files. The one structural gap is that `run_scm()` forwards nonlinear
parameters only through `...` rather than explicit named parameters — functionally
correct but fragile and undiscoverable by users. The two documentation gaps
(missing vignette coverage, missing export of `plot_nonlinear_effect`) are the only
issues that need addressing before Phase A can be considered shipped.

Phase B (Emax and related mechanistic models for the exposure–response relationship)
is a significantly more complex undertaking than Phase A. The core challenge is that
Emax models are *non-linear in parameters*, requiring non-linear optimisation rather
than column transformations, and have pharmacometric interpretability constraints
(EC50 and Emax must remain on the original concentration scale). The recommended
approach is to build Phase B as a *companion function family* (`fit_emax()`) rather
than extending `scm_glm()`/`scm_coxph()`, keeping the SCM covariate-selection
machinery unchanged and adding Emax fitting on top of it. Bootstrap integration and
covariate effects on Emax parameters should be deferred to later sub-phases.

---

## 2. Part 1: Phase A Implementation Status

### Pre-Phase — Transformation Engine

#### `R/nonlinear_spec.R`
- **EXISTS: YES**
- **Functions:** `parse_nonlinear_spec()`, `validate_nonlinear_spec()`,
  `.validate_transform_names()`, `.ALLOWED_TRANSFORMS`, `.DEFAULT_TRANSFORMS`
- **QC pass applied:** YES — full roxygen, `cli::cli_abort()`/`cli::cli_warn()` error
  handling, `\dontrun{}` examples, `@keywords internal` on internals, class attribute
  set on return value.

#### `R/nonlinear_transforms.R`
- **EXISTS: YES**
- **Functions:** `.FP_POWERS`, `.fp_term()`, `.fp2_terms()`, `.fit_model()`,
  `.model_deviance()`, `.make_simple_transform()`, `.make_spline_transform()`,
  `.find_best_fp()`, `.model_df()`, `.rcs_df()`
- **QC pass applied:** YES — complete roxygen, guard clauses, tryCatch around
  fitting, all helpers marked `@keywords internal`.

#### `R/nonlinear_selection.R`
- **EXISTS: YES**
- **Functions:** `.lrt_pvalue()`, `.select_best_transform()`, `.fit_one_transform()`
- **QC pass applied:** YES — full documentation, closed-testing logic described in
  `@details`, safe NA handling throughout.

---

### Phase 1 — Forward Step with Transformations

#### `R/scm_forward_nonlinear.R`
- **EXISTS: YES**
- **Functions:** `scm_forward_first_step_nonlinear()`, `scm_forward_add_nonlinear()`,
  plus internal helpers `.response_lhs()`, `.term_to_varname()`, `.vars_in_model()`,
  `.model_nobs()`, `.build_nonlinear_step_table()`
- **QC pass applied:** YES — well-documented, error-safe, consistent return structure.

#### `scm_glm()` accepts nonlinear parameters
- **DONE** — `nonlinear`, `linearity_alpha`, `fp_powers`, `spline_df` are all explicit
  named parameters at lines 137–140.

#### `scm_coxph()` accepts nonlinear parameters
- **DONE** — same four parameters explicit at lines 200–203.

#### Both functions branch to nonlinear-aware forward step
- **DONE** — `use_nonlinear` flag set, `scm_forward_first_step_nonlinear()` called
  when TRUE (scm_glm.R:329–348, scm_coxph.R:395–413); subsequent covariate addition
  loops call `scm_forward_add_nonlinear()` (scm_glm.R:425–444).

#### `transform_registry` stored in fitted object metadata
- **DONE** — stored at `scm_glm.R:658`, `scm_coxph.R:721`. Also stores
  `nonlinear_spec`, `linearity_alpha`, `fp_powers`, `spline_df` in metadata.

---

### Phase 2 — Backward Step and Output Adaptation

#### `R/scm_backward_nonlinear.R`
- **EXISTS: YES**
- **Functions:** `scm_backward_delete_nonlinear()`, `.enhance_summary_nonlinear()`

#### Closed-testing simplification logic (FP2 → linear → remove)
- **DONE** — simplification phase at lines 95–140. Tests whether the linear form is
  adequate before full removal. Only simplifies when `simplify = TRUE` (default) and
  `df > 1`.

#### Multi-df terms removed as units
- **DONE** — the backward step removes a variable by its `formula_rhs` from
  `transform_registry`, not term-by-term.

#### `transform_registry` updated after backward step
- **DONE** — registry entry is set to NULL for the removed variable (line 298),
  simplified to `"linear"` entry when simplification occurs (lines 120–130).

#### Path table shows transformation labels
- **DONE** — `Transform` column added to `outtab`/`toto_back1` tables; labels come
  from `transform_registry[[v]]$label`.

---

### Phase 3 — Visualisation

#### `R/plot_nonlinear.R`
- **EXISTS: YES**
- **Functions:** `plot_nonlinear_effect()` (exported via `@export`),
  `plot_all_nonlinear()` (exported), `.build_prediction_grid()`,
  `.predict_partial_effect()`, `.transform_subtitle()`,
  `.compute_quantile_effects()`

#### **CRITICAL GAP — NAMESPACE NOT UPDATED**
- `plot_nonlinear_effect` and `plot_all_nonlinear` both have `@export` roxygen tags
  in the source, but **neither appears in `NAMESPACE`**. The NAMESPACE file ends at
  `vpc_coxph` without listing these two functions. This means they are not accessible
  to users after `library(erSCM)` — they exist only as unexported symbols.
- **Fix required:** Run `devtools::document()` or manually add
  `export(plot_nonlinear_effect)` and `export(plot_all_nonlinear)` to `NAMESPACE`.

#### Function handles GLM and CoxPH
- **DONE** — `is_surv <- inherits(fit, "eres_scm_coxph")` check, coxph uses `type =
  "lp"`, GLM uses `type = "link"`.

#### Function handles log, sqrt, FP1, FP2, and RCS
- **DONE** — the function works at the level of `fit$final_model` predictions, not
  by transform type; all forms are handled uniformly through `predict()`.

#### `.compute_quantile_effects()` helper
- **DONE** — present at line 546; `@keywords internal`.

#### Link and response scales, CI ribbons, rug plots
- **DONE** — `scale` argument, CI ribbon via `geom_ribbon`, rug via `geom_rug`.

---

### Phase 4 — Bootstrap Integration

#### `.extract_scm_settings()` extracts nonlinear settings
- **DONE** — `nonlinear_spec`, `linearity_alpha`, `fp_powers`, `spline_df` extracted
  at lines 117–121 of `scm_bootstrap_helpers.R`.

#### `.run_single_replicate()` passes nonlinear parameters
- **DONE** — bootstrap worker (`.run_bootstrap_replicates()` inner function) appends
  `nonlinear`, `linearity_alpha`, `fp_powers`, `spline_df` to `scm_args` when
  `nonlinear_spec` is non-NULL (lines 541–554).
- **Note:** These are passed via `do.call(run_scm, scm_args)`, and `run_scm()` passes
  them through `...` to `scm_glm()`/`scm_coxph()`. This works but is implicit.

#### Per-replicate results record transformation selected
- **DONE** — `.extract_boot_terms()` extracts `exposure_transform` and
  `covariate_transforms` from `transform_registry` (lines 199–230).

#### `scm_bootstrap()` returns `transform_selection_freq` and `exposure_transform_freq`
- **DONE** — both fields present in output list (lines 478–479 of `scm_bootstrap.R`).

#### `print.scm_bootstrap()` displays transform selection frequencies
- **DONE** — "Non-Linear Form Selection" section at lines 171–179 of
  `scm_bootstrap_methods.R`, calls `.print_transform_freq_block()`.

#### `plot.scm_bootstrap()` supports `type = "transforms"`
- **DONE** — implemented at line 592, calls `.plot_transform_stability()`.

---

## 3. Part 2: `run_scm()` Status and Recommendation

### 1. Does `run_scm()` exist?
**Yes.** Defined in `R/run_scm.R`, exported in NAMESPACE. Signature:

```r
run_scm(eres_data, exposure_vars, covariates, type = c("glm", "coxph"),
        significance_level = c(0.05, 0.001), dvid = NULL, dvid_col = NULL,
        dv_col = "DV", time_col = "TIME", family = "binomial",
        var_spec = NULL, scale_to_median = FALSE, center = TRUE,
        adjust_for = NULL, start_covariates = NULL, direction = "both",
        interactions = NULL, interaction_strategy = "after_main",
        p_interaction = NULL, max_interaction_order = 2L,
        verbose = FALSE, return_models = FALSE, ...)
```

### 2. Does `run_scm()` accept and forward nonlinear parameters?
**Partially.** `run_scm()` does **not** declare `nonlinear`, `linearity_alpha`,
`fp_powers`, or `spline_df` as explicit named parameters. They are forwarded
implicitly through `...` to `scm_glm()`/`scm_coxph()`. The bootstrap
`worker_fn` adds these to `scm_args` and calls `do.call(run_scm, scm_args)`,
which works because the extra names reach `scm_glm()`/`scm_coxph()` via `...`.

**Problem:** Users who call `run_scm()` directly will get no autocomplete hint,
no formal documentation, and no validation for these four parameters. If they
misspell `nonlinear` it will silently be ignored.

### 3. Current recommended entry point
For users who want SCM with non-linear terms, the current recommended path is:

```r
scm_glm(eres_data = ..., nonlinear = TRUE, ...)
scm_coxph(eres_data = ..., nonlinear = TRUE, ...)
```

`run_scm()` is currently positioned as the *programmatic* interface for tools
(like `scm_bootstrap()`), not as the primary user-facing function.

### 4. Recommendation
**Add the four nonlinear parameters as explicit named parameters in `run_scm()`.**
Do not redesign the function — simply promote `nonlinear`, `linearity_alpha`,
`fp_powers`, and `spline_df` from `...` to formal parameters with the same
defaults as in `scm_glm()`/`scm_coxph()`. This:

- Fixes discoverability (tab-completion, `?run_scm` shows the parameters)
- Adds formal validation (type/range checks)
- Does not break the bootstrap worker (named args in `scm_args` list still work)
- Does not require redesigning `scm_bootstrap()` which continues to use `run_scm()`

Keeping `run_scm()` as the single dispatcher is the right architecture.
`scm_bootstrap()` uses it internally, and it provides a stable API surface.
Adding nonlinear params explicitly tightens it without changing anything else.

---

## 4. Part 3: Phase A Gaps and Issues

### Critical

| # | Severity | File(s) | Issue | Suggested Fix |
|---|----------|---------|-------|---------------|
| 1 | **Critical** | `NAMESPACE` | `plot_nonlinear_effect` and `plot_all_nonlinear` are `@export`ed in source but absent from `NAMESPACE`. Users cannot call them after `library(erSCM)`. | Run `devtools::document()` or add the two `export()` lines manually to NAMESPACE. |

### Important

| # | Severity | File(s) | Issue | Suggested Fix |
|---|----------|---------|-------|---------------|
| 2 | **Important** | `R/run_scm.R` | Nonlinear params forwarded through `...` — undiscoverable, no validation, silent mistyping. | Promote `nonlinear`, `linearity_alpha`, `fp_powers`, `spline_df` to explicit named parameters with defaults `NULL`, `0.05`, `.FP_POWERS`, `3L`. |
| 3 | **Important** | `vignettes/` | No vignette covers the `nonlinear` feature. Users reading `erSCM-introduction.Rmd` or `bootstrap-stability.Rmd` get no guidance on Phase A. | Add a `vignettes/nonlinear-transformations.Rmd` vignette demonstrating `nonlinear = TRUE`, `nonlinear = list(...)`, `plot_nonlinear_effect()`, and the bootstrap transform stability plot. |
| 4 | **Important** | `R/scm_forward_nonlinear.R` | When `nonlinear` is active and interactions are also requested, the forward interaction step still uses the *linear* forward add internally. Interaction terms involving an FP2-transformed variable would be constructed as `log(CMAX):AGE` etc., which the standard `scm_forward_add_interaction()` may not handle correctly. | Document this limitation explicitly; add a test that confirms what happens and either fix or document as unsupported. |
| 5 | **Important** | `R/scm_backward_nonlinear.R` | The simplification path tests `FP2 → linear` but does NOT test `FP2 → FP1 → linear` (the full Royston & Altman MFP closed-testing hierarchy). | For completeness, the simplification should also test FP1 before linear. This is a design gap vs the MFP specification. |

### Minor

| # | Severity | File(s) | Issue | Suggested Fix |
|---|----------|---------|-------|---------------|
| 6 | **Minor** | `tests/testthat/` | No test covers `nonlinear` + `adjust_for` together — a combination where `adjust_for` terms appear in the model formula alongside transformed covariates. | Add a test in `test-nonlinear_forward.R` with `adjust_for = "RACE"` and `nonlinear = list(CMAX = c("linear","log"))`. |
| 7 | **Minor** | `tests/testthat/` | No test covers `nonlinear` + `start_covariates` (backward-only direction with pre-selected non-linear start terms). | Add a test in `test-nonlinear_backward.R` for `direction = "backward"` with a start set containing a non-linear variable. |
| 8 | **Minor** | `R/plot_nonlinear.R` | For GLM response scale, the code applies `exp()` to the centred log-OR predictions (line 329) rather than converting back to probability via the inverse logit. Comment says "OR relative to reference" but the y-axis label says "OR (relative to median)". The interpretation is correct (HR/OR relative to median), but the code path for `inv_link` (line 321) is computed but never used. | Remove the dead `ref_lp`/`inv_link` computation or document the deliberate choice to show ORs (not probabilities) on the response scale. |
| 9 | **Minor** | `R/scm_bootstrap_helpers.R` | `.compute_transform_freq()` only tracks the first exposure variable's transform when multiple NL exposure candidates exist (`ev <- exp_nl_vars[[1L]]`). For multi-exposure NL runs this silently drops exposure 2+. | Extend to track all NL exposure candidates, keyed by name. |

---

## 5. Part 4: Phase B Detailed Plan — Emax Models for Exposure Variables

### 4.1 API Design

#### Main Function Name Options

| Option | Tradeoffs |
|--------|-----------|
| `fit_emax()` | Clean, descriptive, not tied to "SCM". Makes it clear this is a separate fitting procedure, not an SCM extension. Easy to discover. Precedent from the `DoseFinding` package (`fitMod()`). **Recommended.** |
| `scm_emax()` | Implies it's an extension of the SCM framework, which is misleading — Emax fitting is mechanistic, not stepwise. Would confuse users expecting an SCM path table. |
| `er_nonlinear()` | Too generic — "non-linear" could mean Phase A column transforms. Doesn't signal mechanistic/pharmacometric intent. |
| `eres_emax()` | Package-namespaced, but `eres` already means "exposure–response", so redundant. |

**Recommendation: `fit_emax()`** with a companion `compare_er_shapes()` for the
model comparison workflow.

#### Proposed Signature

```r
fit_emax(
  data,
  exposure,               # column name of the PK metric (string)
  response,               # column name of the response (string); replaces dv_col
  time     = NULL,        # for TTE/CoxPH: column name of time variable
  type     = c("logistic", "coxph", "linear_gaussian"),
  shapes   = c("emax", "linear", "log_linear", "sigmoid"),
                          # candidate shapes to fit
  covariates      = NULL, # character vector; enter linearly on the intercept
  adjust_for      = NULL, # forced linear covariates (same semantics as scm_glm)
  emax_covariates = NULL, # named list: list(Emax = c("WT"), EC50 = c("WT"))
  family   = "binomial",  # for type = "logistic"
  start    = NULL,        # named list of user-provided starting values
  var_spec = NULL,
  ci_method = c("delta", "profile", "bootstrap"),
  ci_level  = 0.95,
  n_boot    = 200L,       # used only when ci_method = "bootstrap"
  subject_col = NULL,
  dvid     = NULL,
  dvid_col = NULL,
  verbose  = FALSE
)
```

#### Return Object Class: `"eres_emax"`

Proposed return structure:
```
$fits          # Named list: one fitted model per shape
$best_shape    # Character: name of best-fitting shape (by AIC)
$comparison    # tibble: shape, n_params, log_lik, AIC, BIC, delta_AIC
$coef          # Named list: parameter estimates per shape
$vcov          # Named list: variance-covariance matrices per shape
$ci            # Named list: CI tables per shape
$data          # The modelling dataset
$exposure      # Exposure column name
$response      # Response column name
$type          # "logistic", "coxph", etc.
$covariates    # Linear covariates
$emax_covariates  # Covariate effects on Emax parameters
$call
```

S3 methods: `print.eres_emax`, `summary.eres_emax`, `plot.eres_emax`,
`coef.eres_emax`, `confint.eres_emax`, `AIC.eres_emax`.

#### Separate companion function

```r
compare_er_shapes(data, exposure, response, type, shapes, ...)
```

Returns only the `$comparison` table, without the full fit objects — useful for
quick screening before committing to `fit_emax()`.

---

### 4.2 Algorithm Design

#### Model Shapes and Fitting Strategy

For all shapes, `link(P)` refers to the link function appropriate for the model
family (logit for logistic regression, log for log-hazard in Cox PH, identity for
Gaussian).

| Shape | Model | Fitting |
|-------|-------|---------|
| **Linear** | `link(P) = β₀ + β₁ · C` | `stats::glm()` / `survival::coxph()` — no change from current |
| **Log-linear** | `link(P) = β₀ + β₁ · log(C)` | `stats::glm()` — already in Phase A via `nonlinear = "log"` |
| **Emax** | `link(P) = β₀ + Emax · C / (EC50 + C)` | Non-linear: `stats::nls()` with custom predictor; or custom NLL optimizer |
| **Sigmoid Emax (Hill)** | `link(P) = β₀ + Emax · C^γ / (EC50^γ + C^γ)` | Same as Emax but with additional Hill coefficient γ |

**Recommended fitting approach: `stats::nls()` with `stats::glm()` linearisation
for starting values, wrapped in tryCatch with grid-search fallback.**

Rationale for `stats::nls()` over alternatives:
- `bbmle::mle2()`: more flexible but adds a suggested dependency; better suited for
  custom likelihoods that don't fit the response ~ f(pred) paradigm.
- `DoseFinding::fitMod()`: purpose-built, but only handles normal/logistic outcomes
  and uses its own model parameterisation that doesn't align with adding linear
  covariates on the intercept.
- Custom NLL optimizer: maximum flexibility (handles CoxPH naturally) but requires
  the most implementation effort and testing.

**Recommended path:**
- For **logistic** outcomes: use `stats::nls()` on the link-transformed response,
  with offset for covariates. Alternatively, maximise the binomial log-likelihood
  directly using `stats::optim()` with analytical gradient.
- For **CoxPH**: cannot use `stats::nls()`. Must build the partial likelihood
  explicitly and optimise with `stats::optim()`. This is the most complex sub-case;
  consider deferring Emax-CoxPH to Phase B.2 or later.
- For **Gaussian** (if ever added): `stats::nls()` directly.

#### Starting Value Strategy
1. Fit the log-linear GLM to get a rough relationship.
2. Back-calculate Emax as the maximum predicted `link(P)` minus the minimum,
   EC50 as the concentration at which the predicted value is mid-range.
3. Set γ = 1 (Emax start).
4. Run a coarse grid over (Emax, EC50) values (5×5 grid on observed data range) and
   pick the pair with lowest residual sum of squares.
5. Use the grid winner as `start` if user did not provide starting values.

#### Standard Errors
- **Delta method** (default): compute gradient of parameter transformations using
  numerical differentiation; requires `stats::vcov()` from the fitted nls object.
- **Profile likelihood** (optional): `stats::confint.nls()` — more accurate for
  asymmetric parameters like EC50, but slow.
- **Bootstrap** (optional): subject-level resampling, re-fit each replicate. Most
  robust but computationally expensive.

#### Model Comparison
Because Emax and linear are **not nested**, AIC/BIC are the primary comparison tools.
For nested pairs (sigmoid Emax vs Emax: one extra parameter γ), a likelihood ratio
test is valid and should be reported alongside AIC.

The comparison table should include:
- Shape name, number of parameters, log-likelihood, AIC, BIC, ΔAIC (vs best), ΔBIC
- For nested pairs: LRT chi-square and p-value
- Convergence flag

#### Convergence Failures
- If `nls()` throws an error (singular gradient, step factor reduced, etc.): record
  the failure, return `NULL` for that shape's fit, and proceed with others.
- Provide `verbose = TRUE` output that reports which shapes converged.
- If *all* non-linear shapes fail to converge, return the linear fit with a warning.

---

### 4.3 Covariate Effects on Emax Parameters

#### Standard Pharmacometric Parameterisation
The industry standard (following Holford, Sheiner, and PK/PD modelling conventions
used by NONMEM/Monolix) is the **multiplicative/exponential parameterisation** on
Emax and EC50:

```
Emax_i  = Emax_pop × exp(β_Emax_WT  × (WT_i  - WT_ref))
EC50_i  = EC50_pop × exp(β_EC50_WT  × (WT_i  - WT_ref))
γ_i     = γ_pop   (Hill coefficient usually not covariate-modelled)
```

Reference values (`WT_ref`) should be the median of the covariate in the dataset.

Alternative: additive parameterisation (`Emax_i = Emax_pop + β × (WT_i - WT_ref)`)
is simpler to estimate but risks `Emax_i < 0` for large covariate effects — avoid.

#### User Specification
Via `emax_covariates` argument (named list):
```r
emax_covariates = list(
  Emax = c("WT"),         # covariates on Emax
  EC50 = c("WT", "AGE"),  # covariates on EC50
  Hill = NULL             # Hill coefficient (no covariate effects)
)
```

#### Interaction with Linear Covariates on Intercept
The `covariates` argument controls linear additive effects on `β₀` (the baseline
response), *not* on the Emax curve parameters. These are entirely separable in the
parameterisation:

```
link(P_i) = (β₀ + β_AGE·AGE_i) + Emax_i · C_i / (EC50_i + C_i)
```

This is the recommended approach: keep intercept-level covariates (demographic
confounders) linear, and allow mechanistic PK-relevant covariates (WT, organ function)
to modify the Emax parameters.

---

### 4.4 Integration with Existing Package Features

#### Bootstrap Stability
**Recommend deferring to Phase B.5.** The uncertainty in Emax parameters is already
addressed by the CI method (`delta`, `profile`, or `bootstrap` within `fit_emax()`).
A separate SCM-style bootstrap of the *model selection* step (which shapes to include)
is less informative for Emax — there are only 4 candidate shapes, so bootstrap
instability is a weaker concept than for covariate selection.

If implemented, it should be a separate function (`bootstrap_emax()`) not an extension
of `scm_bootstrap()`, because the underlying procedure is fundamentally different.

#### Forest Plots
`fit_emax()` should NOT try to reuse `make_forest_plot()`. The forest plot for an
Emax model shows different quantities:
- Emax (maximum effect), EC50 (potency), γ (steepness) and their CIs
- Covariate effects on these parameters (if `emax_covariates` is specified)

A dedicated `plot_emax_params()` function is appropriate.

#### `subject_col` for TTE data
**Include from the start.** Pass `subject_col` through to the bootstrap resampler
(if implemented) and to the data preparation step.

#### Compatibility with `plot_nonlinear_effect()`
**Do not reuse.** `plot_nonlinear_effect()` operates on `eres_scm` objects and assumes
the model is already fitted by SCM with a `transform_registry`. An Emax fit is a
fundamentally different object. Use a dedicated `plot.eres_emax()` method.

#### `adjust_for`, `start_covariates`, `direction`
- `adjust_for`: **Include** — forced linear covariates on the intercept, same semantics.
- `start_covariates`, `direction`: **Exclude** — these apply to the SCM stepwise
  algorithm; Emax fitting is not stepwise. If the user wants covariate selection
  among the linear covariates, they should run `scm_glm()` first (see section 4.1),
  then pass the selected `adjust_for` to `fit_emax()`.

---

### 4.5 Visualisation

#### Primary Plot: Exposure–Response Curve
The single most important output. Should show:
- Predicted curve (on response scale: probability or HR) vs concentration
- Pointwise CI ribbon (from delta/profile/bootstrap SE)
- Observed data as a rug or binned proportions (with CIs) overlaid
- Reference line at the median exposure
- Parameter annotations (Emax, EC50, γ)

#### Secondary Plot: Shape Comparison Overlay
A single panel showing all converged candidate shapes overlaid (different
colours/linetypes), with AIC values in the legend. Allows visual inspection of
whether shapes differ meaningfully in the observed data range.

#### Parameter Forest Plot
A forest-plot-style display of Emax, EC50, γ (and covariate effects on these
parameters), with confidence intervals. Similar to `make_forest_plot()` in style
but with parameter-specific axis formatting (EC50 on original concentration scale,
Emax on the link scale or effect scale).

#### Implementation
Deliver as S3 `plot.eres_emax(type = c("curve", "comparison", "params"))` — a
single method with a `type` switch, following the pattern of `plot.scm_bootstrap()`.

---

### 4.6 Statistical Considerations

#### DVID Subsetting
For datasets with multiple AEs (multiple DVIDs), fit Emax models **per DVID** as
separate `fit_emax()` calls, not pooled. A pooled model with DVID as a linear
covariate is rarely meaningful for mechanistic Emax interpretation because the
Emax and EC50 will likely differ by AE type.

If the user wants a pooled model with DVID effects, they can pass `dvid` as one of
the `covariates` (linear intercept effect). Document this clearly.

#### Scaling and Centering
**Scaling the exposure variable is STRONGLY inadvisable for Emax models.** EC50 has
a physical meaning on the original concentration scale (the concentration producing
50% of maximum effect). If the user scales CMAX to median = 1, the estimated EC50
will be in units of "multiples of median CMAX" — hard to interpret clinically.

Recommendation: **Always fit Emax models on the original concentration scale.** The
function should:
1. Default `scale_to_median = FALSE` and `center = FALSE` for Emax fits.
2. Warn if the user passes a pre-scaled dataset (detectable if max(C) ≈ 1).
3. Document the scaling decision prominently.

#### Minimum Sample Size
Emax models require:
- At least 3 distinct concentration levels below EC50, at EC50, and above EC50 to
  estimate all 3 parameters reliably.
- A minimum of ~60–100 observations for reliable EC50 estimation with logistic
  response; fewer may suffice for Gaussian.
- The sigmoid Emax (4 parameters) typically requires ≥150 observations.

Emit informational warnings when `nrow(data)` is below these thresholds.

---

### 4.7 Suggested Implementation Phases for Phase B

#### Phase B Pre-Phase: Fitting Infrastructure
**What:** Non-linear model fitting utilities isolated from the user-facing function.

- `R/emax_shapes.R`: Define model shape functions (linear predictor given
  concentration and parameters), gradient functions for delta method, and starting
  value generators.
- `R/emax_fit_engine.R`: Low-level fitters for each shape (wrapping `stats::nls()`
  or `stats::optim()`), with convergence logging. Returns a canonical `emax_fit`
  list (estimates, vcov, log-likelihood, AIC, convergence status).
- Unit tests: convergence on simulated data for each shape; confirm that EC50 is
  recovered within ±20% under moderate sample sizes.

#### Phase B.1: Core Fitting — Single Exposure, Linear Covariates on Intercept
**What:** The main `fit_emax()` function.

- Fit all shapes in `shapes` argument.
- Linear covariates on `β₀` only.
- Delta-method CIs for all parameters.
- `comparison` table (AIC/BIC table).
- S3 methods: `print`, `summary`, `coef`, `confint`, `AIC`.
- Basic `plot.eres_emax(type = "curve")` for the winning shape.

#### Phase B.2: Shape Comparison Visualisation
**What:** `compare_er_shapes()` companion and `plot.eres_emax(type = "comparison")`.

- All candidate curves overlaid on observed data.
- AIC/BIC table display in `print.eres_emax`.
- `plot.eres_emax(type = "comparison")`.

#### Phase B.3: Covariate Effects on Emax Parameters
**What:** `emax_covariates` argument support.

- Extend `emax_fit_engine.R` to accept individual-level Emax/EC50 computed from
  covariate effects.
- Add `plot.eres_emax(type = "params")` parameter forest plot.
- Tests for covariate effects on Emax vs EC50.

#### Phase B.4: Profile Likelihood and Bootstrap CIs
**What:** `ci_method = "profile"` and `ci_method = "bootstrap"`.

- Profile likelihood via `MASS::confint.nls()` or custom profiling.
- Bootstrap resampling (subject-level) for CI.
- Comparison of CI methods in documentation.

#### Phase B.5: Bootstrap Stability for Shape Selection
**What:** `bootstrap_emax()` assessing which shape is selected most often.

- N bootstrap samples → fit all shapes → record which has lowest AIC.
- Return shape selection frequency table.
- Plot stacked bar analogous to `plot.scm_bootstrap(type = "transforms")`.

---

### 4.8 Risks and Open Questions

#### Technical Risks

1. **Convergence failures with small samples.** Emax models are notoriously hard
   to fit when the observed concentration range does not span EC50 (plateau not
   reached). Risk: the grid-search starting value strategy may fail if data are
   all below EC50. Mitigation: document the requirement; warn when max(C) < 2 ×
   estimated EC50 from the log-linear fit.

2. **Identifiability.** When the data do not cover the plateau, Emax is poorly
   identified and its CI becomes very wide or unbounded. The delta method will then
   produce nonsensical CIs. Profile likelihood or bootstrap is essential in these
   cases. Mitigation: check identifiability heuristically (is the upper CI of Emax
   > 10 × the estimate?) and emit a warning.

3. **CoxPH Emax fitting.** The partial likelihood has no simple closed-form
   expression for the gradient w.r.t. Emax parameters. `stats::optim()` with
   numerical gradient will be slow and numerically fragile. Mitigation: defer
   Emax-CoxPH to Phase B.3 or later; start with logistic only.

4. **Scaling conflict with Phase A.** The SCM data prep pipeline in `scm_glm()`
   applies `scale_to_median` and `center`. If a user runs Phase A and then passes
   the returned `model_data` (which is scaled) to `fit_emax()`, the EC50 will be
   on the scaled scale. Mitigation: document clearly; `fit_emax()` should default
   to `scale_to_median = FALSE` and always check whether the input data looks
   pre-scaled.

#### Design Risks

5. **Over-engineering the covariate-on-Emax-parameters feature.** The most valuable
   use case for Emax models in exposure-response analysis is the simple 3-parameter
   Emax without covariate effects. Spending significant effort on the multiplicative
   covariate parameterisation (Phase B.3) may deliver little user value compared to
   getting the core curve and comparison working robustly.

6. **Compatibility creep.** Attempting to integrate `fit_emax()` too tightly with
   `scm_bootstrap()` or reusing the existing `eres_scm` infrastructure will create
   maintenance burden. The two are better kept as parallel tools that share only the
   dataset preparation utilities.

#### Open Questions Requiring User Input

1. **Primary use case:** Should `fit_emax()` be primarily a *confirmatory* tool
   (the user already chose CMAX as exposure via SCM, now fits an Emax curve for
   reporting), or a *model-selection* tool (the user doesn't know if the relationship
   is linear or Emax, and wants automatic selection)?

2. **CoxPH support in Phase B.1:** Is Emax for TTE endpoints (CoxPH) needed
   immediately, or can Phase B.1 target logistic only?

3. **Desired CI method for primary deliverable:** Delta method is fastest and most
   familiar, profile likelihood is more accurate. Which is the default?

4. **Covariate selection for Emax:** Should `fit_emax()` integrate its own stepwise
   covariate selection for the linear intercept covariates, or should users always
   run `scm_glm()` first and pass selected covariates as `adjust_for`? The latter
   is cleaner but requires a two-step workflow.

5. **Integration with `DoseFinding` package:** The `DoseFinding` package has mature
   infrastructure for dose-response modelling (including Emax, sigmoid Emax,
   quadratic, etc.) with `fitMod()`, `MCPMod()`, etc. Should `fit_emax()` be built
   on top of `DoseFinding`, or built from scratch to retain full control over
   parameterisation, link functions, and output format? Key trade-off: building on
   `DoseFinding` is faster but constrains the API and adds a dependency that may not
   be on CRAN for all target environments.

#### Relevant Literature and Packages

- Sheiner LB, Stanski DR (1982). Simultaneous modeling of pharmacokinetics and
  pharmacodynamics. *Clin Pharmacol Ther*.
- Holford NHG, Sheiner LB (1982). Kinetics of pharmacologic response. *Pharmacol Ther*.
- Mager DE, Jusko WJ (2008). Pharmacodynamic modeling of time-dependent
  transduction systems. *Clin Pharmacol Ther*.
- **DoseFinding** package (Bretz et al.): `fitMod()`, `MCPMod()` — dose-response
  focused, normal/binary outcomes only.
- **nlme** package: `gnls()` — handles non-linear mixed effects; overkill for
  population-level Emax but relevant if individual-level PK data are available.
- Monolix/NONMEM conventions: multiplicative covariate parameterisation on CL and V;
  directly applicable to EC50 parameterisation in Phase B.3.

---

## 6. Recommended Next Steps

### Immediate (before proceeding to Phase B)

1. **Fix the NAMESPACE gap** (Gap #1, Critical): run `devtools::document()` to
   regenerate NAMESPACE from roxygen tags. Confirm `plot_nonlinear_effect` and
   `plot_all_nonlinear` are exported.

2. **Promote nonlinear params in `run_scm()`** (Gap #2, Important): add four
   explicit named parameters. Low risk, 15-minute change.

3. **Write the nonlinear vignette** (Gap #3, Important): document the full Phase A
   feature for users before Phase B work begins.

4. **Clarify FP2→FP1→linear simplification** (Gap #5, Important): decide whether
   to implement the full MFP hierarchy or document the current 2-step simplification
   as a deliberate simplification.

### Phase B Sequencing

5. Answer the **open questions** in section 4.8 before starting implementation,
   particularly: (a) primary use case (confirmatory vs selection), (b) CoxPH scope
   in Phase B.1, and (c) DoseFinding dependency decision.

6. Start with **Phase B Pre-Phase** (emax shapes infrastructure), which can be built
   and unit-tested independently before the user-facing API is finalised.

7. Keep **Phase B.3 (covariate effects on Emax)** and **Phase B.5 (bootstrap)**
   as deferred phases — deliver Phase B.1 + B.2 first and get feedback before
   investing in the more complex features.
