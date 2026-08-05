# Contributing to erSCM

## TODO and Future Work Policy

### Inline code comments

`# TODO:` comments are permitted in internal implementation code for notes
invisible to package users. Every `# TODO:` must reference either a GitHub
issue number or a concrete decision criterion:

```r
# TODO(#42): Implement FP1 intermediate step when a user dataset surfaces the need.
# TODO: revisit if rcs() df selection ever needs to be auto-tuned.
```

### Roxygen documentation

**`@section Future Work:` is NOT permitted in exported function documentation.**
Future-work statements in user-facing docs create maintenance debt and mislead
users about current capabilities.

**`@section Limitations:` is also discouraged.** Use `@section Notes:` for
permanent current-behavior caveats that users need to know about. If the
content describes something that will be fixed, it belongs in a GitHub issue,
not in roxygen.

Design limitations that users need to know go in:
- The package vignettes (especially `vignette("nonlinear-transformations")`).
- `NEWS.md` under a "Known Limitations" subsection of the current release.
- `@section Notes:` in exported function docs (for permanent behavioral constraints).

### Summary

| Location | `# TODO:` | `@section Future Work:` | `@section Notes:` | Vignette/NEWS |
|----------|:---------:|:-----------------------:|:-----------------:|:-------------:|
| Internal functions | Yes | No | N/A (not user-facing) | — |
| Exported functions | No | **No** | Yes (permanent constraints only) | Yes |
| Vignette / NEWS.md | — | — | — | Yes (all limitations) |
