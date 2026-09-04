#' Extract the LRT p-value from an \code{anova()} comparison table
#'
#' Locates the p-value column by name so the extraction is robust to the
#' differing \code{anova()} layouts of \code{glm} (column \code{"Pr(>Chi)"}) and
#' \code{coxph} (column \code{"Pr(>|Chi|)"}) objects, rather than relying on a
#' hard-coded positional index that silently breaks if the table layout changes.
#'
#' @param anova_tab An \code{anova} comparison table, e.g. the result of
#'   \code{anova(null_fit, full_fit, test = "Chisq")}.
#' @return Numeric scalar: the p-value from row 2 (the "vs. null/reduced" row,
#'   matching the interaction-LRT usage).
#' @keywords internal
.anova_pvalue <- function(anova_tab) {
  cols <- colnames(anova_tab)
  # Locate the chi-square p-value column by name rather than a hard-coded
  # positional index. The exact spelling differs by model class: `glm` anova()
  # uses "Pr(>Chi)", while `survival::coxph` anova() uses "P(>|Chi|)" (some
  # survival methods spell it "Pr(>|Chi|)"). The regex matches all variants.
  pcol <- grep("^Pr?\\(>.*Chi", cols)
  if (length(pcol) == 0L) {
    cli::cli_abort(c(
      "Could not locate the LRT p-value column in the {.fun anova} table.",
      "x" = "No chi-square p-value column (e.g. {.val Pr(>Chi)}, {.val P(>|Chi|)}) is present.",
      "i" = "Available columns: {.val {cols}}."
    ))
  }
  anova_tab[2L, pcol[1L]]
}

#' Forward Addition of Interaction Terms (GLM)
#'
#' Executes a single forward-addition iteration for interaction terms in the SCM
#' workflow. Tests each eligible interaction term (whose parent main effects are
#' already in the model) and selects the best qualifying term based on LRT.
#'
#' @param current_model List. The current model state containing \code{Relation},
#'   \code{finalforward}, \code{toto_tot}, and \code{flagforst}.
#' @param eres_data data.frame. Dataset used for fitting.
#' @param candidate_interactions Character vector. Eligible interaction terms to test
#'   (already filtered by hierarchy constraint).
#' @param family Character or family object. GLM family. Default \code{"binomial"}.
#' @param significance_level Numeric. Entry threshold for interaction terms.
#' @param verbose Logical. Print progress messages. Default \code{FALSE}.
#' @param return_models Logical. Return fitted model objects. Default \code{FALSE}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{finalforward}: Run number of selected model.
#'   \item \code{toto_int}: tibble of this step's results.
#'   \item \code{toto_tot}: cumulative tibble of all steps.
#'
#'   \item \code{Relation}: updated character vector of model terms.
#'   \item \code{flagforst}: logical, whether to continue.
#'   \item \code{models} (optional): fitted candidate models.
#' }
#'
#' @examples
#' data(demo_eres)
#' demo_eres$DV <- demo_eres$sim_status
#'
#' # Build a current_model with two main-effect terms already selected (see
#' # scm_forward_add() for details; this bookkeeping is normally done for you
#' # inside scm_glm()).
#' best1 <- scm_forward_first_step(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", dvid = NULL,
#'   family = "binomial", significance_level = 0.2, run_tag = "step1"
#' )
#' x <- best1$step1_sum
#' current_model <- list(
#'   toto_tot = x,
#'   Relation = x$Model[x$Model == x$SELECT][!is.na(x$Model[x$Model == x$SELECT])],
#'   finalforward = x$Run[x$Model == x$SELECT][1],
#'   flagforst = any(!is.na(x$Run[x$Model == x$SELECT]))
#' )
#' step_add <- scm_forward_add(
#'   current_model = current_model, eres_data = demo_eres,
#'   covariates = c("AGE", "WT"), dvid = NULL, family = "binomial",
#'   significance_level = 0.2, run_tag = "fwd_add_1"
#' )
#'
#' step_int <- scm_forward_add_interaction(
#'   current_model = step_add, eres_data = demo_eres,
#'   candidate_interactions = "CMAXss:WT",
#'   family = "binomial", significance_level = 0.2
#' )
#' step_int$Relation
#'
#' @export
scm_forward_add_interaction <- function(
    current_model,
    eres_data,
    candidate_interactions,
    family = "binomial",
    significance_level = 0.01,
    verbose = FALSE,
    return_models = FALSE
) {
  modds <- eres_data

  ## Current model state

  tmp_ref <- current_model$finalforward
  Relation <- current_model$Relation
  toto <- current_model$toto_tot

  # Filter out interactions already in the model
  candidate_interactions <- candidate_interactions[!candidate_interactions %in% Relation]

  if (length(candidate_interactions) == 0) {
    return(list(
      finalforward = current_model$finalforward,
      toto_int = NULL,
      toto_tot = current_model$toto_tot,
      Relation = current_model$Relation,
      flagforst = FALSE
    ))
  }

  # Build base formula from current Relation
  tmp_cov_ref <- unique(Relation)
  tmp_base_cov <- paste(tmp_cov_ref, collapse = " + ")

  # Build candidate formulas: base + each interaction

  tmp_cov_combo <- paste(tmp_base_cov, "+", candidate_interactions)

  # Fit base (null for this step) model

  mod_nul <- stats::glm(build_model_formula("DV", tmp_base_cov), data = modds, family = family)
  mod_nul_sum <- broom::glance(mod_nul)
  df_ref <- mod_nul$df.null - mod_nul$df.residual
  glance_ref <- broom::glance(mod_nul)

  # Initialize lists

  mod_obj <- vector(mode = "list", length = length(candidate_interactions))
  mod_sum <- vector(mode = "list", length = length(candidate_interactions))
  mod_coef <- vector(mode = "list", length = length(candidate_interactions))
  mod_metric <- vector(mode = "list", length = length(candidate_interactions))
  summ.list <- vector(mode = "list", length = length(candidate_interactions))

  for (i in seq_along(candidate_interactions)) {
    mod_obj[[i]] <- stats::glm(build_model_formula("DV", tmp_cov_combo[i]), data = modds, family = family)
    mod_sum[[i]] <- broom::glance(mod_obj[[i]])
    mod_coef[[i]] <- broom::tidy(mod_obj[[i]], conf.int = TRUE)
    mod_metric[i] <- candidate_interactions[i]
    summ.list[[i]] <- summary(mod_obj[[i]])
  }

  metrics <- unlist(mod_metric)
  nsub.tab <- do.call(rbind, mod_sum) %>%
    dplyr::select(logLik, AIC, nobs) %>%
    dplyr::mutate(Metric = metrics)

  if (verbose) {
    message(paste(utils::capture.output(print(nsub.tab)), collapse = "\n"))
  }

  outtab <- NULL
  for (i in seq_along(metrics)) {
    dftot <- mod_obj[[i]]$df.null - mod_obj[[i]]$df.residual
    tabi <- summ.list[[i]]
    tabi.conf <- as.data.frame(quiet_confint(mod_obj[[i]], level = 0.95))
    tabi.coef <- data.frame(tabi$coefficients)
    out.tab <- dplyr::bind_cols(tabi.coef, tabi.conf)
    dftest <- dftot - df_ref
    loglik <- mod_sum[[i]]$logLik

    out.tab <- out.tab %>%
      dplyr::mutate(Metric = metrics[i],
                    Effect = row.names(out.tab),
                    RSE = 100 * (Std..Error / abs(Estimate)),
                    estRSE = sprintf("%s (%s)", signif(Estimate, 3), signif(RSE, 3)),
                    CI = sprintf("(%s, %s)", signif(`2.5 %`, 3), signif(`97.5 %`, 3)),
                    pvalue = signif(Pr...z.., 3),
                    dftest = dftest, dftot = dftot, logLik = loglik) %>%
      dplyr::select(Effect, Metric, CI, dftest, dftot)

    nul <- mod_nul
    fit <- mod_obj[[i]]
    pval.fitx <- round(.anova_pvalue(anova(nul, fit, test = "Chisq")), 10)
    out.tab$pvalue <- pval.fitx
    outtab <- rbind(outtab, out.tab)
  }

  # Only take covariates (remove intercept)
  outtab <- outtab[outtab$Effect != "(Intercept)", ]
  outtab <- outtab %>%
    dplyr::filter(!is.na(pvalue)) %>%
    dplyr::left_join(nsub.tab, by = "Metric") %>%
    dplyr::mutate(min2LL = -2 * logLik) %>%
    dplyr::mutate(diff2LL = -2 * glance_ref$logLik - min2LL) %>%
    dplyr::select(Metric, min2LL, diff2LL, AIC, pvalue, nobs, dftest, dftot)
  outtab <- outtab[!duplicated(outtab[, c("Metric")]), ]

  if (verbose) {
    message(paste(utils::capture.output(print(outtab)), collapse = "\n"))
  }

  toto_int <- dplyr::bind_rows(
    tibble::tibble(Model = "Base",
                   pvalue = NA_real_,
                   Signif = NA_real_,
                   AIC = mod_nul$aic,
                   min2LL = -2 * mod_nul_sum$logLik,
                   diff2LL = 0,
                   nobs = mod_nul_sum$nobs,
                   dftest = 0,
                   dftot = df_ref),
    outtab %>%
      dplyr::rename(Model = Metric) %>%
      dplyr::mutate(Model = as.character(Model),
                    SELECT = ifelse(min(pvalue, na.rm = TRUE) < significance_level[1],
                                    Model[AIC == min(AIC, na.rm = TRUE)],
                                    NA_character_),
                    Signif = ifelse(pvalue < significance_level[1], 1, 0),
                    Reference = tmp_ref)) %>%
    dplyr::mutate(Run = -1 + (max(toto$Run) + 1):(max(toto$Run) + dplyr::n())) %>%
    dplyr::relocate(Run) %>%
    tidyr::fill(SELECT, .direction = "downup") %>%
    dplyr::mutate(RUNTEST = ifelse(Run != 1, paste0("Run", Reference, " + ", Model), NA_character_))

  toto_int <- toto_int[-1, ]
  toto_tot <- dplyr::bind_rows(current_model$toto_tot, toto_int)

  if (!all(is.na(toto_int$Run[toto_int$Model == toto_int$SELECT]))) {
    finalforward <- toto_int$Run[toto_int$Model == toto_int$SELECT]
    Relation <- c(Relation, toto_int$Model[toto_int$Model == toto_int$SELECT])
  } else {
    finalforward <- current_model$finalforward
    Relation <- current_model$Relation
  }

  flagforst <- FALSE
  if (any(!is.na(toto_int$Run[toto_int$Model == toto_int$SELECT]))) {
    flagforst <- TRUE
  }

  toto_int <- toto_int %>%
    dplyr::select(Run, min2LL, diff2LL, AIC,
                  Reference, pvalue, Signif,
                  nobs, dftest, dftot, RUNTEST,
                  SELECT, Model)

  if (return_models) {
    names(mod_obj) <- metrics
    res <- list(finalforward = finalforward,
                toto_int = toto_int,
                toto_tot = toto_tot,
                Relation = Relation,
                flagforst = flagforst,
                models = mod_obj)
  } else {
    res <- list(finalforward = finalforward,
                toto_int = toto_int,
                toto_tot = toto_tot,
                Relation = Relation,
                flagforst = flagforst)
  }

  return(res)
}


#' Forward Addition of Interaction Terms (Cox PH)
#'
#' Executes a single forward-addition iteration for interaction terms in the SCM
#' workflow for Cox proportional hazards models. Tests each eligible interaction term
#' (whose parent main effects are already in the model) and selects the best qualifying
#' term based on LRT.
#'
#' @param current_model List. The current model state containing \code{Relation},
#'   \code{finalforward}, \code{toto_tot}, and \code{flagforst}.
#' @param eres_data data.frame. Dataset used for fitting (must contain TIME and DV columns).
#' @param candidate_interactions Character vector. Eligible interaction terms to test
#'   (already filtered by hierarchy constraint).
#' @param significance_level Numeric. Entry threshold for interaction terms.
#' @param verbose Logical. Print progress messages. Default \code{FALSE}.
#' @param return_models Logical. Return fitted model objects. Default \code{FALSE}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{finalforward}: Run number of selected model.
#'   \item \code{toto_int}: tibble of this step's results.
#'   \item \code{toto_tot}: cumulative tibble of all steps.
#'   \item \code{Relation}: updated character vector of model terms.
#'   \item \code{flagforst}: logical, whether to continue.
#'   \item \code{models} (optional): fitted candidate models.
#' }
#'
#' @examples
#' data(demo_eres)
#' demo_eres$DV <- demo_eres$sim_status
#' demo_eres$TIME <- demo_eres$sim_time
#'
#' # Build a current_model with two main-effect terms already selected (see
#' # scm_forward_add_coxph() for details; this bookkeeping is normally done
#' # for you inside scm_coxph()).
#' best1 <- scm_forward_first_step_coxph(
#'   eres_data = demo_eres, exposure_vars = "CMAXss", dvid = NULL,
#'   significance_level = 0.2, run_tag = "step1"
#' )
#' x <- best1$step1_sum
#' current_model <- list(
#'   toto_tot = x,
#'   Relation = x$Model[x$Model == x$SELECT][!is.na(x$Model[x$Model == x$SELECT])],
#'   finalforward = x$Run[x$Model == x$SELECT][1],
#'   flagforst = any(!is.na(x$Run[x$Model == x$SELECT]))
#' )
#' step_add <- scm_forward_add_coxph(
#'   current_model = current_model, eres_data = demo_eres,
#'   covariates = c("AGE", "WT"), dvid = NULL,
#'   significance_level = 0.2, run_tag = "fwd_add_1"
#' )
#'
#' step_int <- scm_forward_add_interaction_coxph(
#'   current_model = step_add, eres_data = demo_eres,
#'   candidate_interactions = "CMAXss:WT",
#'   significance_level = 0.2
#' )
#' step_int$Relation
#'
#' @export
scm_forward_add_interaction_coxph <- function(
    current_model,
    eres_data,
    candidate_interactions,
    significance_level = 0.01,
    verbose = FALSE,
    return_models = FALSE
) {
  modds <- eres_data

  ## Current model state
  tmp_ref <- current_model$finalforward
  Relation <- current_model$Relation
  toto <- current_model$toto_tot

  # Filter out interactions already in the model
  candidate_interactions <- candidate_interactions[!candidate_interactions %in% Relation]

  if (length(candidate_interactions) == 0) {
    return(list(
      finalforward = current_model$finalforward,
      toto_int = NULL,
      toto_tot = current_model$toto_tot,
      Relation = current_model$Relation,
      flagforst = FALSE
    ))
  }

  # Build base formula from current Relation
  tmp_cov_ref <- unique(Relation)
  tmp_base_cov <- paste(tmp_cov_ref, collapse = " + ")

  # Build candidate formulas: base + each interaction
  tmp_cov_combo <- paste(tmp_base_cov, "+", candidate_interactions)

  # Fit base model
  mod_nul <- survival::coxph(build_model_formula("survival::Surv(TIME, DV)", tmp_base_cov), data = modds)
  mod_nul_tab <- summary(mod_nul)
  mod_nul_sum <- broom::glance(mod_nul)
  df_ref <- mod_nul_tab$logtest[["df"]]
  glance_ref <- broom::glance(mod_nul)

  # Initialize lists
  mod_obj <- vector(mode = "list", length = length(candidate_interactions))
  mod_sum <- vector(mode = "list", length = length(candidate_interactions))
  mod_coef <- vector(mode = "list", length = length(candidate_interactions))
  mod_metric <- vector(mode = "list", length = length(candidate_interactions))
  summ.list <- vector(mode = "list", length = length(candidate_interactions))

  for (i in seq_along(candidate_interactions)) {
    mod_obj[[i]] <- survival::coxph(build_model_formula("survival::Surv(TIME, DV)", tmp_cov_combo[i]), data = modds)
    mod_sum[[i]] <- broom::glance(mod_obj[[i]])
    mod_coef[[i]] <- broom::tidy(mod_obj[[i]], conf.int = TRUE)
    mod_metric[i] <- candidate_interactions[i]
    summ.list[[i]] <- summary(mod_obj[[i]])
  }

  metrics <- unlist(mod_metric)
  nsub.tab <- do.call(rbind, mod_sum) %>%
    dplyr::select(logLik, AIC, nobs = n) %>%
    dplyr::mutate(Metric = metrics)

  if (verbose) {
    message(paste(utils::capture.output(print(nsub.tab)), collapse = "\n"))
  }

  outtab <- NULL
  for (i in seq_along(metrics)) {
    tabi <- summ.list[[i]]
    dftot <- tabi$logtest[["df"]]
    tabi.conf <- as.data.frame(confint(mod_obj[[i]], level = 0.95))
    tabi.coef <- data.frame(tabi$coefficients)
    out.tab <- dplyr::bind_cols(tabi.coef, tabi.conf)
    dftest <- dftot - df_ref
    loglik <- mod_sum[[i]]$logLik

    out.tab <- out.tab %>%
      dplyr::mutate(Metric = metrics[i],
                    Effect = row.names(out.tab),
                    RSE = 100 * (se.coef. / abs(coef)),
                    estRSE = sprintf("%s (%s)", signif(coef, 3), signif(RSE, 3)),
                    CI = sprintf("(%s, %s)", signif(`2.5 %`, 3), signif(`97.5 %`, 3)),
                    pvalue = signif(Pr...z.., 3),
                    dftest = dftest, dftot = dftot, logLik = loglik) %>%
      dplyr::select(Effect, Metric, CI, dftest, dftot)

    nul <- mod_nul
    fit <- mod_obj[[i]]
    pval.fitx <- round(.anova_pvalue(anova(nul, fit, test = "Chisq")), 10)
    out.tab$pvalue <- pval.fitx
    outtab <- rbind(outtab, out.tab)
  }

  # Only take covariates (remove intercept)
  outtab <- outtab[outtab$Effect != "(Intercept)", ]
  outtab <- outtab %>%
    dplyr::filter(!is.na(pvalue)) %>%
    dplyr::left_join(nsub.tab, by = "Metric") %>%
    dplyr::mutate(min2LL = -2 * logLik) %>%
    dplyr::mutate(diff2LL = -2 * glance_ref$logLik - min2LL) %>%
    dplyr::select(Metric, min2LL, diff2LL, AIC, pvalue, nobs, dftest, dftot)
  outtab <- outtab[!duplicated(outtab[, c("Metric")]), ]

  if (verbose) {
    message(paste(utils::capture.output(print(outtab)), collapse = "\n"))
  }

  toto_int <- dplyr::bind_rows(
    tibble::tibble(Model = "Base",
                   pvalue = NA_real_,
                   Signif = NA_real_,
                   AIC = mod_nul$aic,
                   min2LL = -2 * mod_nul_sum$logLik,
                   diff2LL = 0,
                   nobs = mod_nul_sum$n,
                   dftest = 0,
                   dftot = df_ref),
    outtab %>%
      dplyr::rename(Model = Metric) %>%
      dplyr::mutate(Model = as.character(Model),
                    SELECT = ifelse(min(pvalue, na.rm = TRUE) < significance_level[1],
                                    Model[AIC == min(AIC, na.rm = TRUE)],
                                    NA_character_),
                    Signif = ifelse(pvalue < significance_level[1], 1, 0),
                    Reference = tmp_ref)) %>%
    dplyr::mutate(Run = -1 + (max(toto$Run) + 1):(max(toto$Run) + dplyr::n())) %>%
    dplyr::relocate(Run) %>%
    tidyr::fill(SELECT, .direction = "downup") %>%
    dplyr::mutate(RUNTEST = ifelse(Run != 1, paste0("Run", Reference, " + ", Model), NA_character_))

  toto_int <- toto_int[-1, ]
  toto_tot <- dplyr::bind_rows(current_model$toto_tot, toto_int)

  if (!all(is.na(toto_int$Run[toto_int$Model == toto_int$SELECT]))) {
    finalforward <- toto_int$Run[toto_int$Model == toto_int$SELECT]
    Relation <- c(Relation, toto_int$Model[toto_int$Model == toto_int$SELECT])
  } else {
    finalforward <- current_model$finalforward
    Relation <- current_model$Relation
  }

  flagforst <- FALSE
  if (any(!is.na(toto_int$Run[toto_int$Model == toto_int$SELECT]))) {
    flagforst <- TRUE
  }

  toto_int <- toto_int %>%
    dplyr::select(Run, min2LL, diff2LL, AIC,
                  Reference, pvalue, Signif,
                  nobs, dftest, dftot, RUNTEST,
                  SELECT, Model)

  if (return_models) {
    names(mod_obj) <- metrics
    res <- list(finalforward = finalforward,
                toto_int = toto_int,
                toto_tot = toto_tot,
                Relation = Relation,
                flagforst = flagforst,
                models = mod_obj)
  } else {
    res <- list(finalforward = finalforward,
                toto_int = toto_int,
                toto_tot = toto_tot,
                Relation = Relation,
                flagforst = flagforst)
  }

  return(res)
}
