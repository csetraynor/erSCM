# Skill: R Package Code Review, Refactoring & QC

## Purpose
Systematically review, refactor, document, and quality-control R package source code to ensure it is production-ready, readable, maintainable, and CRAN-compliant.

## When to Use
Invoke this skill when asked to "review", "refactor", "QC", "clean up", or "improve code quality" for an R package.

## Execution Protocol

### Phase 1: Inventory & Architecture Scan

1. **List all R source files** in `R/` directory, sorted by size (largest first)
2. **Map the dependency graph**: which files source/call functions from which other files
3. **Identify the public API**: all exported functions (from NAMESPACE or @export tags)
4. **Identify internal helpers**: all `.`-prefixed or non-exported functions
5. **Check file organisation**: do files group logically? (one concept per file, or sprawl?)

Report a brief architecture summary before proceeding.

### Phase 2: Per-File Deep Review

For each `.R` file in `R/`, perform these checks IN ORDER:

#### A. Structure & Organisation
- [ ] File has a clear single responsibility (or is clearly a collection of related helpers)
- [ ] Functions are ordered logically (public first, then private helpers)
- [ ] No function exceeds 150 lines (flag for refactoring if so)
- [ ] No excessive nesting (max 4 levels of indentation)
- [ ] Related utility functions are grouped, not scattered across files

#### B. Function Design
- [ ] Each function does ONE thing (single responsibility)
- [ ] Function names are descriptive and consistent (snake_case for exported, .dot_prefix for internal)
- [ ] Parameter names are consistent across similar functions
- [ ] Default values are sensible and documented
- [ ] No magic numbers — use named constants or parameters
- [ ] Early returns for validation (guard clauses) before main logic
- [ ] No deeply nested if/else — refactor to early returns or switch/helper functions

#### C. Documentation (roxygen2)
- [ ] Every exported function has complete roxygen2 block
- [ ] `@param` for every parameter with type and description
- [ ] `@return` with complete description of return value structure
- [ ] `@examples` present (in \dontrun{} if slow)
- [ ] `@seealso` links to related functions
- [ ] `@family` tag for function grouping
- [ ] `@keywords internal` for non-exported helpers
- [ ] Internal functions have at minimum: brief title, @param, @return, @keywords internal
- [ ] No outdated/stale documentation (params that don't exist, etc.)

#### D. Code Style & Readability
- [ ] Consistent spacing: `x <- 1` not `x<-1` or `x = 1`
- [ ] Consistent indentation (2 spaces per level)
- [ ] Line length ≤ 100 characters (prefer ≤ 80 for readability)
- [ ] Blank lines separate logical blocks
- [ ] Complex operations have inline comments explaining WHY (not WHAT)
- [ ] No commented-out dead code
- [ ] No `TODO` or `FIXME` without an associated plan
- [ ] Pipe usage is consistent (either `|>` or `%>%`, not mixed)

#### E. Inline Comments
Apply this commenting philosophy:
- **DO comment**: WHY something is done (business logic, algorithm choice, edge case handling)
- **DO comment**: non-obvious consequences of a code choice
- **DO comment**: references to algorithms, papers, or external specs
- **DO NOT comment**: WHAT the code does (if the code is clear, the comment is noise)
- **DO NOT over-comment**: every line doesn't need a comment
- **Section headers**: Use `# --- Section Name ---` or `# ── Section Name ──` for major logical blocks within long functions

Comment template for complex functions:
```r
#' @description
#' [Brief description]
#'
#' Algorithm:
#' 1. [Step 1]
#' 2. [Step 2]
#' ...
Expand code
F. Error Handling & Defensive Programming
•	Input validation at function entry (type, length, range, existence)
•	cli::cli_abort() for user-facing errors (not stop())
•	cli::cli_warn() for recoverable issues
•	cli::cli_inform() for informational messages
•	tryCatch() around operations that can fail (model fitting, file I/O)
•	Informative error messages: state what was expected, what was received, and what to do
•	No naked try() — always handle the error explicitly
G. Performance & Efficiency
•	No unnecessary copies of large objects
•	vapply() preferred over sapply() (type-stable)
•	Avoid growing vectors/lists in loops (pre-allocate)
•	No repeated computation that could be cached
•	seq_len(n) instead of 1:n (handles n=0 correctly)
•	seq_along(x) instead of 1:length(x)
H. Namespace & Dependencies
•	All external functions called with pkg::fun() syntax
•	No library() or require() calls
•	Suggests packages guarded with rlang::check_installed() or rlang::is_installed()
•	@importFrom used sparingly and only where :: would be excessively verbose
•	No reliance on non-standard evaluation without proper .data$ pronouns
I. Testability
•	Functions have clear inputs and outputs (testable without side effects)
•	Complex internal logic is factored into small testable helpers
•	No hard-coded paths, system calls, or environment dependencies
•	Random operations use seeds passed as parameters (not internal set.seed())
Phase 3: Cross-File Consistency
After reviewing all files individually:
•	Function naming conventions are consistent across the package
•	Error message style is consistent (cli formatting)
•	Parameter names are consistent when the same concept appears in multiple functions (e.g., significance_level vs alpha vs sig_level — pick one)
•	Return value structures are consistent for similar functions
•	Internal data structures (model state, registries) have consistent format
•	Similar operations are factored into shared helpers (DRY principle)
Phase 4: CRAN Compliance Check
•	globalVariables() declares all NSE variables (ggplot aes, dplyr, etc.)
•	No ::: calls to other packages (not allowed on CRAN)
•	All examples run in < 5 seconds or are in \dontrun{}
•	No hardcoded paths, no system-specific code
•	Encoding: UTF-8 in DESCRIPTION
•	All Suggests are conditionally loaded
•	NAMESPACE is auto-generated by roxygen2 (no manual edits)
•	Package passes R CMD check --as-cran with 0 errors, 0 warnings
Phase 5: Test Coverage Assessment
•	Every exported function has at least one test
•	Every internal helper called from multiple places has a test
•	Edge cases are tested (NULL inputs, empty vectors, single elements, NA values)
•	Error paths are tested (expect_error())
•	Warning paths are tested (expect_warning())
•	Integration tests exist for end-to-end workflows
Refactoring Patterns for R
Pattern: Long Function → Extract Helper
# BEFORE: 200-line function with embedded logic
my_function <- function(x) {
  # ... 50 lines of validation ...
  # ... 80 lines of computation ...
  # ... 70 lines of output formatting ...
}

# AFTER: Composed of focused helpers
my_function <- function(x) {
  .validate_my_input(x)
  result <- .compute_my_result(x)
  .format_my_output(result)
}
Expand code
Pattern: Repeated Conditional → Dispatch Table
# BEFORE:
if (type == "glm") {
  result <- do_glm_thing(...)
} else if (type == "coxph") {
  result <- do_coxph_thing(...)
} else if (type == "poisson") {
  result <- do_poisson_thing(...)
}

# AFTER:
.dispatch <- list(
  glm = do_glm_thing,
  coxph = do_coxph_thing,
  poisson = do_poisson_thing
)
result <- .dispatch[[type]](...)
Expand code
Pattern: Deep Nesting → Early Return
# BEFORE:
f <- function(x) {
  if (!is.null(x)) {
    if (length(x) > 0) {
      if (is.numeric(x)) {
        # actual logic here
      }
    }
  }
}

# AFTER:
f <- function(x) {
  if (is.null(x)) return(NULL)
  if (length(x) == 0) return(numeric(0))
  if (!is.numeric(x)) cli::cli_abort("x must be numeric")
  # actual logic here — at top indentation level
}
Expand code
Pattern: Repeated tryCatch → Wrapper
# BEFORE: tryCatch repeated 10 times with same error handling
result1 <- tryCatch(fit_model(data1), error = function(e) NULL)
result2 <- tryCatch(fit_model(data2), error = function(e) NULL)

# AFTER:
.safe_fit <- function(data) {
  tryCatch(fit_model(data), error = function(e) {
    cli::cli_inform("Model fitting failed: {conditionMessage(e)}")
    NULL
  })
}
result1 <- .safe_fit(data1)
result2 <- .safe_fit(data2)
Expand code
Pattern: Magic Numbers → Named Constants
# BEFORE:
if (n_converged / n_boot < 0.8) warn(...)
bar <- strrep("#", round(prob * 20))

# AFTER:
.CONVERGENCE_WARN_THRESHOLD <- 0.8
.PROGRESS_BAR_WIDTH <- 20L

if (n_converged / n_boot < .CONVERGENCE_WARN_THRESHOLD) warn(...)
bar <- strrep("#", round(prob * .PROGRESS_BAR_WIDTH))
Expand code
Output Format
For each file reviewed, produce:
## File: R/filename.R

### Summary
- Lines: N | Functions: N exported, N internal
- Complexity: [Low/Medium/High]
- Quality: [Good/Needs Work/Major Issues]

### Issues Found (by severity)

#### Critical (must fix)
1. [Issue description + location + fix]

#### Important (should fix)
1. [Issue description + location + fix]

#### Minor (nice to have)
1. [Issue description + location + fix]

### Refactoring Applied
1. [What was changed and why]

### Comments Added
1. [Summary of commenting improvements]
Expand code
Completion Criteria
The review is complete when:
1.	Every R/ file has been reviewed
2.	All Critical issues are fixed
3.	All Important issues are fixed or documented with justification for deferral
4.	devtools::document() runs cleanly
5.	devtools::test() passes (no regressions from refactoring)
6.	devtools::check() passes with 0 errors, 0 warnings, ≤ 1 NOTE
