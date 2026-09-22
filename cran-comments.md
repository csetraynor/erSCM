# cran-comments

## Submission

This is the first CRAN submission of erSCM.

## Test environments

- Local: Linux (RHEL 8, kernel 4.18.0), R version 4.5.1 (2025-06-13)
- Additional platforms (R-hub / win-builder): to be verified on compute node before final submission

## R CMD check results

0 errors | 0 warnings | 0 notes

(Full `R CMD check` to be re-verified on compute node before submission.)

## Reverse dependencies

None (first submission).

## Notes for the reviewer

- All exported step-function `@examples` that fit models are wrapped in `\donttest{}`
  to keep example run time within CRAN's 5-second-per-example guideline.
- Heavy integration tests (model-fitting, bootstrap, VPC) are guarded with
  `skip_on_cran()` in six test files; the remaining test suite runs entirely
  within 3 minutes on a single CPU.
- `inst/WORDLIST` was created to cover domain-specific terms (pharmacokinetics,
  pharmacometrics, statistician names, R package names, abbreviations) that are
  not in the standard en-GB dictionary.
