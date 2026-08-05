# NOTE: This script is retained for quick interactive testing.
# For documented workflows, see vignette("introduction", package = "erSCM")

# library(erSCM)
devtools::document()
library(dplyr)

data("demo_eres")

all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
all_exposure   <- c("CMAXss", "CMINss", "AUCss")

# -------------------------------------------------------- #
# Variable metadata (replaces manual normalization + naming_ds)
# -------------------------------------------------------- #
my_spec <- c(
  CMAXss = "Cmax,ss (ng/mL)",
  CMINss = "Cmin,ss (ng/mL)",
  AUCss  = "AUC,ss (ng*h/mL)",
  WT     = "Body Weight (kg)",
  AGE    = "Age (years)",
  EGFR   = "eGFR (mL/min/1.73m2)",
  ALB    = "Albumin (g/dL)",
  BMI    = "BMI (kg/m2)",
  SEX    = "Sex",
  AAG    = "AAG (mg/dL)",
  AST    = "AST (U/L)",
  ALT    = "ALT (U/L)"
)

# -------------------------------------------------------- #
# Minimal data prep (no manual scaling needed)
# -------------------------------------------------------- #
demo_eres <- demo_eres |>
  mutate(SEX = as.factor(SEX), DVID = "ALL")

# -------------------------------------------------------- #
# Cox PH SCM
# -------------------------------------------------------- #
surv_search <- run_scm(
  type               = "coxph",
  eres_data          = demo_eres,
  covariates         = all_covariates,
  exposure_vars      = all_exposure,
  dv_col             = "sim_status",
  time_col           = "sim_time",
  significance_level = c(0.01, 0.05),
  dvid               = "ALL",
  dvid_col           = "DVID",
  var_spec           = my_spec,
  adjust_for         =  "EGFR",      
  center             = TRUE,
  verbose            = TRUE
)




surv_search <- scm_coxph(
  eres_data          = demo_eres,
  covariates         = all_covariates,
  exposure_vars      = all_exposure,
  dv_col             = "sim_status",
  time_col           = "sim_time",
  significance_level = c(0.01, 0.05),
  dvid               = "ALL",
  dvid_col           = "DVID",
  var_spec           = my_spec,
  center             = TRUE,
  verbose            = TRUE
)

surv_search
plot(surv_search)
vpc_coxph(surv_search)
plot_cox_terms(surv_search, "CMAXss")

# -------------------------------------------------------- #
# Logistic regression SCM
# -------------------------------------------------------- #
lr_search <- scm_glm(
  eres_data          = demo_eres,
  covariates         = all_covariates,
  exposure_vars      = all_exposure,
  dv_col             = "sim_status",
  significance_level = c(0.05, 0.01),
  var_spec           = my_spec,
  center             = TRUE,
  verbose            = TRUE
)

lr_search
plot(lr_search)

# -------------------------------------------------------- #
# Legacy approach (still works, but var_spec is preferred)
# -------------------------------------------------------- #
# naming_ds <- data.frame(
#   covname      = c("CMAXss", "WT"),
#   units        = c("ug*day/mL", "kg"),
#   median_point = c(median(demo_eres$CMAXss), median(demo_eres$WT))
# )
# make_forest_plot(surv_search, namings = naming_ds)
