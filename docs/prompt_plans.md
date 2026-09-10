# Prompt for Claude Code — Phase 1c

---

Phase 1b is accepted and committed. Proceed to **Phase 1c only**. Do not start subsequent phases until I confirm.

## Reminder of standing rules (still in force)

- Do **not** run `git add`, `git commit`, `git push`, or `gh pr create`. Read-only git commands are fine.
- Apply the **semantic-equivalence rule**, the **verbose-gating test-suite rule**, and the **print/cat detection rule** established in Phases 1a and 1b.
- Findings that mention a specific symbol apply **only to that symbol** — grep for callers before deleting anything.
- Keep R usage HPC-light: single-threaded, no parallel test runs, no `R CMD check` on the login node.
- Save all test output to `reports/test_phase1c.txt` (and `reports/test_phase1c_error_message.txt` if needed). Overwrite on each run.

## Phase 1c — CRAN blockers

Apply the following in a single phase:

### 1. Sync `DESCRIPTION` version with `NEWS.md` (checklist #6)

- The audit noted `DESCRIPTION: 0.1.2` vs `NEWS.md: 0.1.4 (development)`.
- Confirm the intended release version. Default assumption: bump `DESCRIPTION` to `0.1.4` and mark the top `NEWS.md` entry as the released version (drop the `(development)` tag or replace with today's date, `2026-09-10`).
- If you are unsure whether `0.1.4` or `0.1.3` is correct given the intervening Phase 1a/1b commits, **stop and ask** before editing.

### 2. Add `skip_on_cran()` to heavy test files (checklist #13)

Add `testthat::skip_on_cran()` at the top of each `test_that()` block in the 6 heavy test files identified by the audit:

- `test-scm_bootstrap_helpers.R`
- `test-vpc_coxph_sim-shape.R`
- `test-vpc_coxph_sim-validate.R`
- `test-vpc_logistic-e2e.R`
- `test-vpc_logistic-binning.R`
- `test-vpc_logistic-validate.R`

Do **not** skip on our HPC test runs — `skip_on_cran()` only fires under `NOT_CRAN != "true"`. Confirm the tests still pass locally after the additions.

### 3. Wrap heavy `@examples` in `\donttest{}` (checklist #9)

- Grep the exported low-level step functions (`scm_backward_delete`, `scm_backward_delete_coxph`, `scm_forward_add`, `scm_forward_add_coxph`, `scm_forward_first_step`, `scm_forward_first_step_coxph`, `fit_univariate_models`, `fit_univariate_models_coxph`, and any similar) for `@examples` blocks.
- If the example fits a model or otherwise takes non-trivial time, wrap the model-fitting portion in `\donttest{ ... }`. Keep any lightweight setup (data prep, `data(...)`) outside the wrapper so users still see runnable context.
- Do **not** use `\dontrun{}` — that hides examples entirely from CRAN checks. `\donttest{}` is the correct choice for "runs but is slow".
- Regenerate `.Rd` files (`devtools::document()` or `roxygen2::roxygenise()`) after edits. Confirm the resulting `man/*.Rd` files render without warnings.

### 4. Verify `LICENSE` file state (checklist #24)

- Inspect the current state: is there a bare `LICENSE`, a `LICENSE.md`, both, or neither?
- Read the `License:` field in `DESCRIPTION`.
- Apply CRAN's rules:
  - `License: GPL-3` (or `GPL-2`, `MIT`, etc. — standard templates) → **no** `LICENSE` file needed; delete it if present. A `LICENSE.md` may remain if it just points to the standard text, but `.Rbuildignore` must exclude it, or CRAN will NOTE.
  - `License: MIT + file LICENSE` → `LICENSE` file **required** with the specific two-line template.
  - `License: file LICENSE` → full text required in `LICENSE`.
- State clearly in the phase report which case applies to `erSCM` and what action was taken.

## Deliverables

1. Files changed (list, one-line summary each).
2. `git diff --stat` output.
3. Test suite run summary (`FAIL 0 | WARN 0 | SKIP N` — note that SKIP count will now be higher due to `skip_on_cran()`, which is expected and correct).
4. Confirmation that no `.Rd` warnings were introduced.
5. Proposed commit message (Conventional Commits, one commit for the whole phase is fine here, or split into `chore(desc): sync version`, `test: skip heavy tests on CRAN`, `docs: wrap heavy examples in donttest`, `chore: fix LICENSE state` if cleaner).
6. Proposed PR title and description referencing checklist items #6, #9, #13, #24.
7. Verification steps I should run locally before committing.

Then stop and wait for my review.

---

# Prompt for Claude Code — Phase 1d (queued, do not start yet)

---

**Do not start until Phase 1c is committed and I say "proceed to Phase 1d".**

## Phase 1d — Spelling + cran-comments + Rbuildignore

### 1. Spell check

- Run `spelling::spell_check_package()` on the package.
- Review the output. For each flagged word:
  - **Real typo** → fix in the source file (`.R`, `.Rd`, `.md`, vignette).
  - **Domain term / proper noun / package name** (e.g. `AstraZeneca`, `SCM`, `CoxPH`, `covariate`, `pharmacometrics` terms) → add to `inst/WORDLIST` (create the file if it does not exist).
- Save the pre-fix and post-fix spelling output to `reports/spelling_phase1d.txt`.

### 2. `cran-comments.md`

- Create `cran-comments.md` at the package root with this template, filled in appropriately:

```markdown
## Submission

This is the first CRAN submission of erSCM.

## Test environments

- Local: <OS>, R <version>
- (Additional platforms to be added after R-hub / win-builder runs, deferred to compute node)

## R CMD check results

0 errors | 0 warnings | 0 notes

(To be re-verified on compute node before submission.)

## Reverse dependencies

None (first submission).

## Notes for the reviewer

- <any known NOTEs and their justification>
```

- Fill in `<OS>` and `<R version>` from `R.version.string` on the current machine; leave the compute-node verification note in place.

### 3. `.Rbuildignore`

- Add `^cran-comments\.md$` to `.Rbuildignore` (CRAN requires this — the file must not be in the built tarball).
- Confirm `reports/`, `.claude`, `.Rproj.user`, IDE artefacts, `examples/` (if dev-only) are already excluded per the audit checklist #21.

## Deliverables

Same shape as Phase 1c (files, diff-stat, test summary, commit message, PR body). Stop after Phase 1d.

---

# Prompt for Claude Code — Phase 2 (queued, do not start yet)

---

**Do not start until Phase 1d is committed and I say "proceed to Phase 2".**

## Phase 2 — API hygiene for unused parameters (audit findings #6, #7, #8)

The affected parameters are:

- `scm_coxph`, `fit_univariate_models`, `fit_univariate_models_coxph`: `parallel`
- `scm_backward_delete`, `scm_backward_delete_coxph`: `removable_terms`, `outcome`, `fixed_terms`, `criterion`, `min_terms`, `run_tag`
- `scm_forward_add`, `scm_forward_add_coxph`: `criterion`, `max_terms`, `run_tag`

### Step 1 — Decision table (do this first, do not edit code yet)

For each parameter, produce a table with:

| Function | Parameter | Status (planned feature / abandoned / partially implemented) | Recommended action | Rationale |
|---|---|---|---|---|

Base the "Status" column on:

- `git log -p` for the function to see whether the parameter was recently added (planned feature) or has been dead since inception (abandoned).
- `NEWS.md` mentions.
- Grep the rest of the codebase for any partial implementation.

Recommended actions to choose from:

- **Keep + `.NotYetImplemented()`**: for planned features — parameter stays, non-default values trigger a clear error.
- **Keep + `lifecycle::deprecate_warn()`**: for abandoned parameters we plan to remove in the next minor version. Requires adding `lifecycle` to `Imports`.
- **Silent keep**: only if the parameter has *some* documented effect I have overlooked (state which).
- **Remove**: only if you can prove no external caller and no downstream test uses it. Even then, prefer deprecation for a first CRAN submission.

**Stop and present the table for my approval before editing any function signature or body.**

### Step 2 — Implementation (after I approve the decisions)

Apply the agreed action per parameter. Update roxygen docs to reflect the new status (e.g. add `@param parallel Not yet implemented; setting to \code{TRUE} will error.`). Update `NEWS.md` under a `## Deprecated` or `## Not yet implemented` heading.

## Deliverables

- Phase 2 report as usual.
- The decision table produced in Step 1 must appear verbatim in the report so I can audit the reasoning.

---

# Prompt for Claude Code — Phase 3 (queued, do not start yet)

---

**Do not start until Phase 2 is committed and I say "proceed to Phase 3". This is the highest-risk phase.**

## Phase 3 — GLM/CoxPH duplication collapse (audit finding #1, ~600 LOC)

### Step 1 — Design proposal (do this first, do not edit code yet)

Produce a written design covering:

1. **Shared internal API**
   - Function signature of the shared table-builder (name, arguments, return shape).
   - Function signature of the `._fit_fn` dispatch mechanism (closure? S3? simple function argument?).
   - How the `nobs` vs `n` column rename is handled (rename at source, or in the builder).
2. **File layout**
   - Where the new shared code lives (`R/scm_shared.R`? `R/scm_stepwise_core.R`?).
   - Which functions in `R/scm_glm.R` and `R/scm_coxph.R` become thin wrappers.
3. **Regression testing plan**
   - Fixed small dataset (name it) usable for both GLM and CoxPH paths.
   - Which top-level entry points will be snapshotted (`scm_glm`, `scm_coxph` end-to-end, plus at least one intermediate step function each).
   - Tolerance strategy: bit-identical where possible; for floating-point outputs, define numerical tolerance and justify.
   - Snapshot storage: `testthat::expect_snapshot_value()` with serialization or `expect_equal()` against a saved `.rds`.
4. **Rollback plan** if the collapse turns out to introduce subtle numerical differences.

**Stop and present the design for my approval before writing any refactor code or tests.**

### Step 2 — Regression tests first (after design approval)

Write the regression tests against the **current (pre-refactor)** codebase. Confirm they pass. Commit boundary marker: these tests must exist and pass before any refactor edit.

### Step 3 — Refactor

Implement the shared internals. Update `scm_glm` and `scm_coxph` code paths to use them. Delete the duplicated code only after the regression tests still pass.

### Step 4 — Verify

Full test suite green (`FAIL 0 | WARN 0`; SKIP fine). Regression snapshots unchanged. Report as usual.

### HPC compute-node marker

Add a `TODO(compute-node)` note in the phase report listing the final validation steps that must run on a compute node before the release tag:

- `R CMD check --as-cran --no-build-vignettes --no-manual --timings`
- `lintr::lint_package()`
- `styler::style_pkg(dry = "on")`
- Example timing verification (5 s CRAN limit per example)
- Full test suite without `skip_on_cran()` gates
- `spelling::spell_check_package()` final run

## Deliverables

- Design document (Step 1).
- New regression tests (Step 2) as a separate proposed commit.
- Refactor diff (Step 3) as a second proposed commit.
- Final phase report with the compute-node runbook appended (which will become `PRE-SUBMISSION.md` in a small trailing phase).

---

Start with **Phase 1c**. Stop and report.