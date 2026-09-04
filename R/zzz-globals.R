# Suppress R CMD check NOTEs for non-standard evaluation variables
utils::globalVariables(c(
# plot_logistic_terms
"bin", "ci_low", "ci_high", "events", "level", "mean_x", "pred_lo", "pred_hi",
"pred_prob", "rate", "source",
# forest_plot
"term_type",
# plot_glm
"bin_midpoint", "ci_lower", "ci_upper", "n_obs", "n_obs.x", "observed_prob",
"predicted_prob", "term_value",

# scm_bootstrap_methods
"stability_class",
"row_var", "col_var", "co_inclusion_value",
"in_original", "label_text", "inclusion_prob",
"term_fct", "fill_val", "tile_col",
"row_idx", "col_idx", "is_diagonal",
"exposure_var", "selection_prob", "is_original", "rank",

# Column names used in dplyr/tidyr NSE
".", ".data", "AE_name", "AIC", "CI", "CI95hig", "CI95low", "CI_VAL", "CIFIN",
"COV", "DV", "DVID", "DVID_name", "Effect", "Estimate", "ID", "Metric",
"Model", "Pr...z..", "Q", "RSE", "RUNTEST", "Reference", "Relation",
"Run", "SELECT", "Signif", "Significant", "Std..Error", "TIME", "VAL", "VAL2",
"ci_label", "coef", "conf.high", "conf.low", "conf.level",
"col_name", "covname", "diff2LL", "display_covname", "dftot", "dftest", "estimate", "estRSE",
"finalforward", "flagforst", "group", "group_level", "group_level_f",
"hi", "id", "label", "lo", "logLik", "lower", "lower2", "med",
"medcov", "median_point", "mid", "min2LL", "n", "n.event", "name",
"nobs", "origin", "paramname", "pvalue", "qhigh", "qlow",
"quartile_cat", "quartile_mean", "row", "se.coef.", "std.error",
"strata", "term", "time", "toto_tot", "units", "upper", "upper2",
"val", "value", "x", "2.5 %", "97.5 %", "p.value", "statistic",
":=",

# plot_nonlinear
"x_val", "predicted", "ci_lower_nl", "ci_upper_nl", "observed",

# scm_bootstrap_methods (Phase 4 -- transform stability plot)
"prop", "transform", "transform_fct", "variable_fct",

# scm_bootstrap_helpers (Phase 6 -- categorical diagnostics)
"variable", "n_events", "issue",

# scm_bootstrap (Phase 6 -- convergence issues tibble)
"replicate", "reason",

# scm_backward_nonlinear / scm_forward_nonlinear NSE column names
"Variable", "Transform", "Label"
))

#' @importFrom stats anova as.formula confint glm median qnorm quantile
#' @importFrom splines ns
NULL
