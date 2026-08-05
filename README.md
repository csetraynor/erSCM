# **erSCM**

**erSCM** is an R package for performing **Stepwise Covariate Modeling (SCM)** in the **ERES** setting. It provides convenient wrappers to run SCM for common endpoint types, including **time-to-event** and **binary** outcomes.

---

## **Installation**

### **From GitHub**

```r
# install.packages("devtools")
devtools::install_github("your-org/erSCM")
```

### **Local development**

```r
devtools::install("path/to/erSCM")
library(erSCM)
```


## **Quick start**

This example uses the built-in dataset **`demo_eres`** and runs SCM for both: **Survival (Cox proportional hazards)** via **`scm_coxph()`** and **Binary (logistic GLM)** via **`scm_glm()`**.

```r

data("demo_eres")

# Optional: inspect available columns
colnames(demo_eres) |> paste(collapse = "','")

all_covariates <- c("WT", "AGE", "EGFR", "ALB", "BMI", "SEX", "AAG", "AST", "ALT")
con_covariates <- c("WT", "AGE", "HT", "EGFR", "ALB", "BMI", "AAG", "SCR", "AST", "ALT")
all_exposure   <- c("CMAXss", "CMINss", "AUCss")

# Normalize continuous covariates / exposure metrics by their medians (recommended)
demo_eres <- demo_eres |>
  mutate(across(c(con_covariates, all_exposure), ~ . / median(., na.rm = TRUE)))

# Add a DVID column (Optional)
demo_eres <- demo_eres |> mutate(DVID = "ALL")
```

## **SCM for time-to-event endpoints (Cox PH)**

Run SCM using a Cox model with a time and event column.

```r
surv_search <- scm_coxph(
  eres_data = demo_eres,
  covariates = all_covariates,
  exposure_vars = all_exposure,
  dv_col = "sim_status",
  time_col = "sim_time",
  significance_level = c(0.05, 0.01),
  dvid = "ALL",
  dvid_col = "DVID"
)

surv_search
```

**Notes**

- **`dv_col`** is the event indicator (commonly **0/1**).
- **`time_col`** is the time-to-event column.
- **`significance_level`** can be a vector (e.g., forward/backward thresholds) depending on the SCM implementation.

## **SCM for binary endpoints (GLM)**

Run SCM using a GLM (typically logistic regression when `dv_col` is binary).

```r
lr_search <- scm_glm(
  eres_data = demo_eres,
  covariates = all_covariates,
  exposure_vars = all_exposure,
  dv_col = "sim_status",
  significance_level = c(0.05, 0.01)
)

lr_search
```

## **Data requirements**

Your input `eres_data` should include:

- An **outcome column** (e.g., `sim_status`)
- For Cox models: a **time column** (e.g., `sim_time`)
- Candidate **covariate columns**
- Candidate **exposure metric columns**
- Optional **subgroup column** if using `dvid_col` (e.g., `DVID`)

## **Interaction terms**

Both `scm_coxph()` and `scm_glm()` support testing **interaction terms** while maintaining hierarchy constraints (an interaction A:B is only tested when both A and B are in the model).

```r
surv_int <- scm_coxph(
  eres_data = demo_eres,
  covariates = all_covariates,
  exposure_vars = all_exposure,
  dv_col = "sim_status",
  time_col = "sim_time",
  significance_level = c(0.05, 0.01),
  dvid = "ALL",
  dvid_col = "DVID",
  interactions = c("CMAXss:AGE", "CMAXss:WT", "CMAXss:SEX"),
  interaction_strategy = "after_main",
  p_interaction = c(0.05, 0.01)
)
```

Key parameters:
- **`interactions`** — character vector of explicit terms or a list of variable groups
- **`interaction_strategy`** — `"after_main"` (default, recommended) or `"simultaneous"`
- **`p_interaction`** — forward/backward p-value thresholds for interactions
- **`max_interaction_order`** — max interaction order (default: 2, pairwise)

See `examples/example_interactions.R` for a comprehensive walkthrough.

## Bootstrap Model Selection Stability

Assess how robust your covariate selection is to sampling variability:

```r
# Fit the SCM model
fit <- scm_glm(
  eres_data          = my_data,
  exposure_vars      = "CMAX",
  covariates         = c("AGE", "WT", "SEX", "RACE"),
  significance_level = c(0.05, 0.01)
)

# Bootstrap stability assessment (500 replicates)
stability <- scm_bootstrap(fit, n_boot = 500, seed = 42)
print(stability)
plot(stability)                         # inclusion probability bar chart
plot(stability, type = "coinclusion")  # co-inclusion heatmap
plot(stability, type = "both")         # combined (requires patchwork)

# Extract reliably selected covariates
stable_covariates(stability, threshold = 0.7)

# Multiple exposure candidates — tracks which exposure wins each replicate
fit_multi <- scm_glm(
  eres_data          = my_data,
  exposure_vars      = c("CMAX", "AUC"),
  covariates         = c("AGE", "WT", "SEX"),
  significance_level = c(0.05, 0.01)
)
stability_multi <- scm_bootstrap(fit_multi, n_boot = 500, seed = 42)
plot(stability_multi, type = "exposure")  # exposure selection bar chart
summary(stability_multi)$exposure_stability  # selection probabilities per exposure
```

The bootstrap procedure replays the full SCM algorithm on resampled datasets,
reporting how frequently each covariate is selected. Key outputs:

- **Inclusion probabilities** — proportion of replicates selecting each covariate
- **Co-inclusion matrix** — identifies competing or complementary covariate pairs
- **Modal model** — the most frequently selected covariate combination (exposure + covariates)
- **Model variability index** — overall stability summary (0 = stable, 1 = variable)
- **Exposure selection probabilities** — which exposure metric wins most often (multi-exposure fits)

For patient-level resampling (multi-row data) use `subject_col = "SUBJID"`.
For parallel execution set `parallel = TRUE` with an active `future::plan()`.

See `vignette("bootstrap-stability")` for detailed guidance and interpretation.

## **Typical workflow**

1. Define candidate **covariates** and **exposure metrics**
2. (Optional) Normalize continuous predictors (e.g., divide by median)
3. Run SCM for the endpoint type:
   - `scm_coxph()` for time-to-event
   - `scm_glm()` for binary
4. (Optional) Specify interaction candidates
5. Review returned results (selected terms, step history, model fit metrics)

