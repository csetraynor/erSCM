## R CMD check results

0 errors | 0 warnings | 1 note

* This is a re-submission (version 0.1.4).

## Notes

1. **"unable to verify current time"** — the build server cannot reach an
   external time server. This is an infrastructure issue unrelated to the
   package.

## Test suite

Bootstrap tests are guarded with `skip_on_cran()` and have been optimized to
use `n_boot` in the range 10–25. The full test suite (excluding
`skip_on_cran()` tests) completes in under 60 seconds. Long-running stability
assessments (n_boot >= 200) are demonstrated in the vignette rather than the
test suite.

## Changes since last submission (0.1.3 → 0.1.4)

* Nonlinear (FP / RCS / log) transformation support added to `scm_glm()` and
  `scm_coxph()`: new `nonlinear`, `linearity_alpha`, `fp_powers`, and
  `spline_df` arguments.
* `plot_nonlinear_effect()` and `plot_all_nonlinear()` exported.
* Bootstrap transformation-frequency tracking fixed for multi-exposure NL fits.
* S3 class always returned from `scm_glm()` / `scm_coxph()` even when no
  exposure enters the model.
* New vignette: *Non-Linear Transformations in SCM*.

## Changes since last submission (0.1.1 → 0.1.2)

* Added **Exposure Selection Stability** tracking to `scm_bootstrap()`:
  new fields `exposure_selection_prob`, `exposure_selected_per_boot`, and
  `inclusion_prob_by_exposure`.
* `modal_model` and `original_selected` now return `list(exposure, covariates)`
  instead of a plain character vector.
* Fixed a metadata bug in `scm_coxph()`: `metadata$candidate_exposure_vars`
  now stores the full candidate list rather than only the selected winner.
* New `plot(type = "exposure")` and `stable_covariates(include_exposure)` options.
* UTF-8 fallback in `.block_bar()` for non-UTF-8 terminals.
* All new variables registered in `utils::globalVariables()` to avoid NOTEs.

## Test environments

* Linux x86_64, R 4.3.2 (local build server)
* win-builder (devel and release) — pending
* R-hub (linux, windows, macos) — pending

## Downstream dependencies

There are no reverse dependencies.
