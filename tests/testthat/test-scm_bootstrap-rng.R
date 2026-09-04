## test-scm_bootstrap-rng.R
## Finding #5 — bootstrap RNG parity.
## scm_bootstrap() must return identical resamples for the same seed regardless
## of parallel = FALSE vs parallel = TRUE, must be deterministic under a fixed
## seed, and must not leak its internal L'Ecuyer-CMRG switch into the caller's
## RNG state.

# ---- Tiny, fast fixture ------------------------------------------------------
# Small synthetic GLM SCM fit: <= 50 rows, <= 3 covariates. Kept local to this
# file so the RNG tests are self-contained and fast.
.is_rstudio_session <- function() {
  identical(Sys.getenv("RSTUDIO"), "1")
}

.fit_rng_test_scm_glm <- function() {
  set.seed(101)
  n <- 50L
  CMAX <- rnorm(n, 100, 20)
  AGE  <- rnorm(n, 55, 12)
  WT   <- rnorm(n, 75, 15)
  # DV is driven by CMAX so the exposure reliably enters the model, giving the
  # bootstrap something to resample (a fit with final_model != NULL).
  lp   <- -6 + 0.06 * CMAX
  DV   <- rbinom(n, 1L, plogis(lp))
  test_data <- data.frame(
    DV   = DV,
    CMAX = CMAX,
    AGE  = AGE,
    WT   = WT,
    SEX  = factor(sample(c("M", "F"), n, replace = TRUE)),
    stringsAsFactors = FALSE
  )
  suppressMessages(
    scm_glm(
      eres_data          = test_data,
      exposure_vars      = "CMAX",
      covariates         = c("AGE", "WT", "SEX"),
      significance_level = c(0.5, 0.3),
      dv_col             = "DV",
      family             = "binomial",
      center             = FALSE,
      scale_to_median    = FALSE,
      verbose            = FALSE
    )
  )
}

# ---- Seed-parity (core acceptance test) --------------------------------------

test_that("scm_bootstrap() resamples are identical for parallel = FALSE vs TRUE", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

  fit <- .fit_rng_test_scm_glm()

  seq_run <- suppressMessages(
    scm_bootstrap(fit, n_boot = 10L, seed = 20260828L, parallel = FALSE)
  )

  # Choose a real parallel backend usable in the current context:
  # * multisession (as specified) spawns fresh worker R sessions that must load
  #   the INSTALLED erSCM package. That is the case under `devtools::check()`.
  # * Under `devtools::test()` / pkgload::load_all() the package is not installed,
  #   so multisession workers cannot find it; there we fall back to multicore,
  #   which forks the current process and therefore shares the loaded namespace
  #   (Linux/macOS only). This keeps the parity assertion executable in both
  #   gates. Skip only if no backend can reach erSCM.
  erscm_installed <- any(vapply(
    .libPaths(),
    function(l) dir.exists(file.path(l, "erSCM")),
    logical(1L)
  ))
  if (erscm_installed) {
    old_plan <- suppressWarnings(
      future::plan(future::multisession, workers = 2L)
    )
  } else if (.Platform$OS.type == "unix" && !.is_rstudio_session()) {
    # Forking shares the already-loaded dev namespace, so it works under
    # pkgload::load_all() where the package is not installed. Fork support is
    # unavailable on Windows and unreliable inside RStudio.
    old_plan <- suppressWarnings(
      future::plan(future::multicore, workers = 2L)
    )
  } else {
    skip("No parallel backend can access erSCM in this context.")
  }
  withr::defer(future::plan(old_plan))

  par_run <- suppressMessages(
    scm_bootstrap(fit, n_boot = 10L, seed = 20260828L, parallel = TRUE)
  )

  # boot_models is the per-replicate selected model (exposure + covariates); it
  # is a deterministic function of the resample, so identical boot_models across
  # the two paths proves the resamples themselves are identical.
  expect_identical(seq_run$boot_models,    par_run$boot_models)
  expect_identical(seq_run$boot_converged, par_run$boot_converged)

  # Aggregated outputs must therefore match too.
  expect_identical(seq_run$inclusion_prob, par_run$inclusion_prob)
})

# ---- Deterministic-seed sanity ------------------------------------------------

test_that("scm_bootstrap() is deterministic under a fixed seed (sequential)", {
  skip_on_cran()

  fit <- .fit_rng_test_scm_glm()

  r1 <- suppressMessages(scm_bootstrap(fit, n_boot = 10L, seed = 4242L))
  r2 <- suppressMessages(scm_bootstrap(fit, n_boot = 10L, seed = 4242L))
  expect_identical(r1$boot_models, r2$boot_models)

  # A different seed should generally produce a different set of resamples.
  r3 <- suppressMessages(scm_bootstrap(fit, n_boot = 10L, seed = 9999L))
  expect_false(identical(r1$boot_models, r3$boot_models))
})

# ---- RNG-state hygiene --------------------------------------------------------

test_that("scm_bootstrap() restores the caller's RNG state", {
  skip_on_cran()

  fit <- .fit_rng_test_scm_glm()

  # Establish a known caller RNG kind and state.
  set.seed(2026L)
  before_kind <- RNGkind()
  before_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  invisible(suppressMessages(scm_bootstrap(fit, n_boot = 10L, seed = 55L)))

  after_kind <- RNGkind()
  after_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  # The internal RNGkind("L'Ecuyer-CMRG") switch must not leak out.
  expect_identical(after_kind, before_kind)
  expect_identical(after_seed, before_seed)
})
