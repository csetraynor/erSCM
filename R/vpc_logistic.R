#' Visual Predictive Check (VPC) for Logistic SCM Models
#'
#' Performs a simulation-based visual predictive check (VPC) for a binary
#' logistic model fit by \code{\link{scm_glm}}. Replicate outcome datasets are
#' simulated under the fitted model (parametric bootstrap of the coefficient
#' vector), a per-bin event proportion is computed for each replicate, and the
#' resulting simulated prediction band is compared with the observed event
#' proportion within covariate bins. This assesses model \emph{calibration}
#' across the range of a covariate.
#'
#' @details
#' \code{vpc_logistic()} is complementary to, and distinct from,
#' \code{\link{plot_logistic_terms}}. \code{plot_logistic_terms()} is a
#' \emph{covariate-effect} plot: it shows the model-implied probability curve at
#' a single reference profile (medians / modes of the other covariates) and is
#' \strong{not} a simulation-based predictive check. \code{vpc_logistic()}, by
#' contrast, simulates whole datasets under the fitted model and contrasts the
#' distribution of a summary statistic (the binned event proportion) with the
#' observed value of that same statistic. Use \code{plot_logistic_terms()} to
#' read off an adjusted covariate effect; use \code{vpc_logistic()} to judge
#' whether the model reproduces the observed event rates.
#'
#' Only binary logistic models (\code{family = "binomial"}) fit by
#' \code{scm_glm()} are supported. Multinomial and ordinal logistic models are
#' out of scope for this version.
#'
#' @section Simulation mechanism:
#' For each of \code{n_sim} replicates a coefficient vector is drawn from
#' \eqn{\mathcal{N}(\hat{\beta}, \widehat{\mathrm{Var}}(\hat{\beta}))} using the
#' Cholesky factor of \code{stats::vcov(fit$final_model)} together with
#' \code{stats::rnorm()} (no external dependency). Fitted probabilities are
#' evaluated on the \emph{observed} design matrix, binary outcomes are drawn with
#' \code{stats::rbinom()}, and the per-bin event proportion is recorded. The band
#' shown is the pointwise \code{probs} quantile envelope of those simulated
#' proportions. Drawing conditional on the observed design isolates parameter
#' uncertainty (plus binomial sampling noise), which is the quantity a
#' calibration VPC targets; a non-parametric row bootstrap is deferred to a
#' future version.
#'
#' @section Reproducibility:
#' Randomness is controlled by \code{seed}. Internally the function derives one
#' independent \code{"L'Ecuyer-CMRG"} stream per replicate (via the shared
#' helper used by \code{\link{scm_bootstrap}}), so results are identical for a
#' given \code{seed} whether \code{parallel = FALSE} or \code{parallel = TRUE}.
#' The caller's RNG kind and \code{.Random.seed} are saved on entry and restored
#' on exit, so the internal switch to \code{"L'Ecuyer-CMRG"} does not leak into
#' the calling session. \code{vpc_logistic()} and \code{scm_bootstrap()} are
#' seed-parity siblings.
#'
#' @param fit A list of class \code{eres_scm_glm} from \code{\link{scm_glm}},
#'   containing \code{fit$final_model} (a binomial \code{glm}),
#'   \code{fit$model_data}, \code{fit$metadata}, and \code{fit$var_spec}.
#' @param covariate_term Character; name of the main-effect covariate to
#'   stratify on. Interaction terms (containing \code{":"}) are not supported.
#' @param n_bins Integer; number of quantile bins for a continuous covariate.
#'   Ignored for factor/character covariates (their levels are the bins).
#'   Default 4.
#' @param n_sim Integer; number of simulated replicates. Default 1000.
#' @param probs Numeric vector of length 2 or 3 giving the band quantiles. A
#'   length-3 vector \code{c(lo, mid, hi)} draws a median line; a length-2 vector
#'   \code{c(lo, hi)} omits it. Default \code{c(0.05, 0.5, 0.95)}.
#' @param seed Integer; master RNG seed. Default 42.
#' @param parallel Logical; if \code{TRUE}, replicates are dispatched with
#'   \code{future.apply::future_lapply()} using the same per-replicate streams as
#'   the sequential path. Default \code{FALSE}.
#' @param plot Logical; if \code{TRUE} (default) return a \code{ggplot} object,
#'   if \code{FALSE} return the per-bin data frame described under
#'   \emph{Value}.
#'
#' @return If \code{plot = TRUE}, a \code{ggplot2} object: observed per-bin event
#'   proportions with Wilson confidence intervals overlaid on the simulated
#'   prediction band (ribbon for continuous covariates, interval bars for
#'   categorical). If \code{plot = FALSE}, a per-bin \code{data.frame} with one
#'   row per bin and columns \code{bin}, \code{level}, \code{mean_x} (bin mean of
#'   the covariate on the model scale; \code{NA} for categorical covariates),
#'   \code{n}, \code{events}, \code{obs_rate}, \code{obs_lo}, \code{obs_hi},
#'   \code{sim_lo}, \code{sim_mid} (\code{NA} when \code{probs} has length 2), and
#'   \code{sim_hi}.
#'
#'   Note: this per-bin return shape deliberately differs from
#'   \code{\link{vpc_coxph}}, which returns a bare long data frame of curves;
#'   the binned-proportion summary is the natural object for a logistic VPC.
#'
#' @seealso \code{\link{vpc_coxph}}, \code{\link{plot_logistic_terms}},
#'   \code{\link{scm_glm}}
#'
#' @examples
#' \donttest{
#' data(demo_eres)
#' demo_eres$SEX <- as.factor(demo_eres$SEX)
#' fit <- scm_glm(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", covariates = c("AGE", "WT"),
#'   dv_col = "sim_status", significance_level = c(0.2, 0.2)
#' )
#' vpc_logistic(fit, "WT", n_sim = 200)
#' vpc_logistic(fit, "WT", n_sim = 200, plot = FALSE)
#' }
#'
#' @importFrom stats vcov coef rnorm rbinom plogis qnorm quantile model.frame
#' @importFrom stats model.matrix model.response
#' @export
vpc_logistic <- function(fit, covariate_term, n_bins = 4, n_sim = 1000,
                         probs = c(0.05, 0.5, 0.95), seed = 42L,
                         parallel = FALSE, plot = TRUE) {

  # ---- Validation ------------------------------------------------------------
  if (!inherits(fit, "eres_scm_glm")) {
    cli::cli_abort(c(
      "{.arg fit} must be an {.cls eres_scm_glm} object from {.fun scm_glm}.",
      "x" = "Got an object of class {.cls {class(fit)}}."
    ))
  }
  if (is.null(fit$final_model)) {
    cli::cli_abort(c(
      "{.arg fit} has no {.field final_model}.",
      "i" = "No covariate was selected, so there is nothing to check."
    ))
  }
  glm_fit <- fit$final_model
  fam <- tryCatch(stats::family(glm_fit)$family, error = function(e) NA_character_)
  if (!identical(fam, "binomial")) {
    cli::cli_abort(c(
      "{.fun vpc_logistic} supports binary logistic models only.",
      "x" = "{.field final_model} has family {.val {fam}}, not {.val binomial}."
    ))
  }
  if (!is.character(covariate_term) || length(covariate_term) != 1L) {
    cli::cli_abort("{.arg covariate_term} must be a single variable name.")
  }
  if (grepl(":", covariate_term, fixed = TRUE)) {
    main_effects <- get_main_effect_columns(fit$metadata$finalRelation)
    cli::cli_abort(c(
      "{.fun vpc_logistic} does not support interaction terms as {.arg covariate_term}.",
      "i" = "Use a main-effect variable name, e.g. one of: {.val {main_effects}}."
    ))
  }
  if (!is.numeric(n_bins) || length(n_bins) != 1L || n_bins < 2) {
    cli::cli_abort("{.arg n_bins} must be a single integer >= 2.")
  }
  if (!is.numeric(n_sim) || length(n_sim) != 1L || n_sim < 1 ||
      n_sim != as.integer(n_sim)) {
    cli::cli_abort("{.arg n_sim} must be a single positive integer.")
  }
  n_sim  <- as.integer(n_sim)
  n_bins <- as.integer(n_bins)
  if (!is.numeric(probs) || !length(probs) %in% c(2L, 3L) ||
      any(probs <= 0) || any(probs >= 1)) {
    cli::cli_abort("{.arg probs} must be a numeric vector of length 2 or 3 in (0, 1).")
  }
  probs      <- sort(probs)
  has_median <- length(probs) == 3L

  # ---- Extract observed design, response, covariate --------------------------
  mf <- stats::model.frame(glm_fit)
  X  <- stats::model.matrix(glm_fit)
  y  <- as.numeric(stats::model.response(mf))

  if (!covariate_term %in% colnames(mf)) {
    cli::cli_abort(c(
      "{.arg covariate_term} = {.val {covariate_term}} is not a predictor in {.field final_model}.",
      "i" = "Available predictors: {.val {setdiff(colnames(mf), colnames(mf)[1])}}."
    ))
  }
  cov_values <- mf[[covariate_term]]
  is_numeric_cov <- is.numeric(cov_values)

  # ---- Drop aliased coefficients so X, beta, vcov are conformable ------------
  beta_hat <- stats::coef(glm_fit)
  keep     <- !is.na(beta_hat)
  beta_hat <- beta_hat[keep]
  X        <- X[, keep, drop = FALSE]
  V        <- stats::vcov(glm_fit)
  # vcov already excludes aliased terms; align defensively by name.
  V        <- V[names(beta_hat), names(beta_hat), drop = FALSE]
  L        <- base::chol(V)          # V = t(L) %*% L, L upper-triangular
  p_dim    <- length(beta_hat)

  # ---- Bin assignment (done ONCE; shared by observed and simulated stats) -----
  # Identity guarantee: both the observed proportion and the simulated band for a
  # given bin are aggregated over the SAME rows via this single factor.
  if (is_numeric_cov) {
    bin_factor <- ggplot2::cut_number(cov_values, n = n_bins)
  } else {
    lev <- if (is.factor(cov_values)) {
      levels(droplevels(cov_values))
    } else {
      sort(unique(as.character(cov_values)))
    }
    bin_factor <- factor(as.character(cov_values), levels = lev, ordered = TRUE)
  }
  bin_factor <- droplevels(bin_factor)
  bin_levels <- levels(bin_factor)

  # ---- Simulation replicates -------------------------------------------------
  worker_fn <- function(i) {
    z      <- stats::rnorm(p_dim)
    beta_b <- beta_hat + as.vector(crossprod(L, z))   # t(L) %*% z ~ N(0, V)
    p_b    <- stats::plogis(as.vector(X %*% beta_b))
    y_b    <- stats::rbinom(length(p_b), size = 1L, prob = p_b)
    as.numeric(tapply(y_b, bin_factor, mean))
  }

  rng     <- .make_rng_streams(n_sim, seed)
  on.exit(rng$restore(), add = TRUE)
  streams <- rng$streams

  if (parallel) {
    sim_list <- future.apply::future_lapply(
      seq_len(n_sim), worker_fn, future.seed = streams
    )
  } else {
    sim_list <- vector("list", n_sim)
    for (i in seq_len(n_sim)) {
      assign(".Random.seed", streams[[i]], envir = .GlobalEnv)
      sim_list[[i]] <- worker_fn(i)
    }
  }
  sim_mat <- do.call(rbind, sim_list)               # n_sim x n_bins
  colnames(sim_mat) <- bin_levels

  # ---- Simulated band --------------------------------------------------------
  q       <- apply(sim_mat, 2L, stats::quantile, probs = probs,
                   names = FALSE, na.rm = TRUE)      # length(probs) x n_bins
  sim_lo  <- q[1L, ]
  sim_hi  <- q[nrow(q), ]
  sim_mid <- if (has_median) q[2L, ] else NA_real_

  sim_df <- data.frame(
    bin     = bin_levels,
    sim_lo  = as.numeric(sim_lo),
    sim_mid = as.numeric(sim_mid),
    sim_hi  = as.numeric(sim_hi),
    stringsAsFactors = FALSE
  )

  # ---- Observed statistic (Wilson CI at the band's coverage) -----------------
  conf <- probs[length(probs)] - probs[1L]
  zc   <- stats::qnorm(1 - (1 - conf) / 2)

  obs_input <- data.frame(
    bin = bin_factor,
    y   = y,
    x   = if (is_numeric_cov) as.numeric(cov_values) else NA_real_,
    stringsAsFactors = FALSE
  )
  obs_df <- obs_input %>%
    dplyr::group_by(.data$bin) %>%
    dplyr::summarise(
      n      = dplyr::n(),
      events = sum(.data$y, na.rm = TRUE),
      mean_x = mean(.data$x, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      obs_rate = .data$events / .data$n,
      obs_lo   = (.data$obs_rate + zc^2 / (2 * .data$n) -
                    zc * sqrt((.data$obs_rate * (1 - .data$obs_rate) +
                                 zc^2 / (4 * .data$n)) / .data$n)) /
        (1 + zc^2 / .data$n),
      obs_hi   = (.data$obs_rate + zc^2 / (2 * .data$n) +
                    zc * sqrt((.data$obs_rate * (1 - .data$obs_rate) +
                                 zc^2 / (4 * .data$n)) / .data$n)) /
        (1 + zc^2 / .data$n),
      bin = as.character(.data$bin)
    )

  # ---- Assemble per-bin frame (explicit join by bin identity) ----------------
  vpc_df <- dplyr::inner_join(obs_df, sim_df, by = "bin")
  vpc_df$level <- vpc_df$bin
  # Preserve the numeric / factor bin order (never lexical) -- Phase B#3.
  vpc_df$bin <- factor(vpc_df$bin, levels = bin_levels, ordered = TRUE)
  vpc_df <- vpc_df[order(vpc_df$bin), , drop = FALSE]
  rownames(vpc_df) <- NULL

  vpc_df <- vpc_df[, c("bin", "level", "mean_x", "n", "events",
                       "obs_rate", "obs_lo", "obs_hi",
                       "sim_lo", "sim_mid", "sim_hi")]

  if (!plot) {
    return(vpc_df)
  }

  # ---- Plot ------------------------------------------------------------------
  var_spec_obj  <- fit$var_spec
  display_label <- get_var_label(var_spec_obj, covariate_term)
  display_units <- get_var_units(var_spec_obj, covariate_term)
  x_label <- if (!is.na(display_units) && !is.null(display_units)) {
    paste0(display_label, " (", display_units, ")")
  } else {
    display_label
  }
  dv_label <- if (!is.null(fit$metadata$dvid_name)) {
    fit$metadata$dvid_name
  } else if (!is.null(fit$metadata$dvid)) {
    fit$metadata$dvid
  } else {
    "Event"
  }
  band_label <- paste0(round(conf * 100), "% prediction band")

  if (is_numeric_cov) {
    x_range  <- range(vpc_df$mean_x, na.rm = TRUE)
    err_w    <- diff(x_range) / (nrow(vpc_df) * 4)
    back_fn  <- .get_back_transform_fn(var_spec_obj, covariate_term)
    x_scale  <- if (identical(back_fn, identity)) {
      ggplot2::scale_x_continuous()
    } else {
      ggplot2::scale_x_continuous(labels = function(x) round(back_fn(x), 2))
    }

    p <- ggplot2::ggplot(vpc_df, ggplot2::aes(x = .data$mean_x)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$sim_lo, ymax = .data$sim_hi),
        alpha = 0.2, fill = "#2166AC"
      )
    if (has_median) {
      p <- p + ggplot2::geom_line(
        ggplot2::aes(y = .data$sim_mid),
        color = "#2166AC", linewidth = 0.8
      )
    }
    p <- p +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$obs_lo, ymax = .data$obs_hi),
        color = "grey40", width = err_w
      ) +
      ggplot2::geom_point(
        ggplot2::aes(y = .data$obs_rate),
        color = "grey40", size = 3
      ) +
      x_scale
  } else {
    p <- ggplot2::ggplot(vpc_df, ggplot2::aes(x = .data$bin)) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = .data$sim_lo, ymax = .data$sim_hi),
        color = "#2166AC", linewidth = 6, alpha = 0.3
      )
    if (has_median) {
      p <- p + ggplot2::geom_point(
        ggplot2::aes(y = .data$sim_mid),
        color = "#2166AC", shape = 95, size = 8
      )
    }
    p <- p +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$obs_lo, ymax = .data$obs_hi),
        color = "grey40", width = 0.2
      ) +
      ggplot2::geom_point(
        ggplot2::aes(y = .data$obs_rate),
        color = "grey40", size = 3
      )
  }

  p +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::labs(
      x = x_label,
      y = paste0("P(", dv_label, ")"),
      subtitle = paste0("Observed vs. simulated event rate  (", band_label, ")")
    ) +
    ggplot2::theme_bw() +
    ggpubr::labs_pubr()
}
