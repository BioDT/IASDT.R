# Species Distribution Modelling Workflow for Single-Species Models

This comprehensive workflow implements single-species species
distribution models (sSDMs) for invasive alien plant species in Europe
at the habitat level. It orchestrates the entire process from data
preparation to model fitting, evaluation, and prediction across current
and future climate scenarios. The workflow employs the `sdm` R package
for model fitting and handles cross-validation, parallel processing, and
various environmental predictors.

## Usage

``` r
fit_sdm_models(
  sdm_method = NULL,
  model_settings = NULL,
  model_dir = NULL,
  hab_abb = NULL,
  cv_type = "cv_dist",
  n_cores = 8L,
  n_cores_check = 8L,
  n_cores_summary = 8L,
  strategy = "multisession",
  future_max_size = 2000L,
  selected_species = NULL,
  excluded_species = NULL,
  env_file = ".env",
  clamp_pred = TRUE,
  fix_efforts = "q90",
  fix_rivers = "q90",
  climate_models = "all",
  climate_scenarios = "all",
  climate_periods = "all",
  copy_maxent_html = TRUE
)
```

## Arguments

- sdm_method:

  Character. A single SDM algorithm to use for fitting models. Valid
  values: "glm", "glmpoly", "gam", "glmnet", "mars", "mars2", "gbm",
  "rf", "ranger", "cart", "rpart", "maxent", "mlp", "rbf", "svm",
  "svm2", "mda", and "fda". These correspond to selected methods
  supported by the `sdm` package. For details and supported options, see
  [`sdm::getmethodNames()`](https://rdrr.io/pkg/sdm/man/add.html).

- model_settings:

  List or NULL. List of model-specific settings. If `NULL`, defaults to
  custom settings defined within the workflow.

- model_dir:

  Character. Path to the directory containing model data and where
  outputs and results will be saved. Model data are prepared using the
  [`mod_prepare_hpc()`](https://biodt.github.io/IASDT.R/reference/Mod_inputs.md)
  and
  [`mod_prepare_data()`](https://biodt.github.io/IASDT.R/reference/Mod_inputs.md)
  functions.

- hab_abb:

  Character. Abbreviation for a single SynHab habitat type. Valid
  values: "0", "1", "2", "3", "4a", "4b", "10", "12a", "12b". See
  [`mod_prepare_hpc()`](https://biodt.github.io/IASDT.R/reference/Mod_inputs.md)
  for more details.

- cv_type:

  Character. Cross-validation type. One of `cv_dist` (default) or
  `cv_large`. See
  [`mod_cv_fit()`](https://biodt.github.io/IASDT.R/reference/Mod_CV_Fit.md)
  for more details.

- n_cores, n_cores_check, n_cores_summary:

  Integer. Number of CPU cores for parallel processing of model fitting,
  model checking, and summarising model outputs. Default is 8.

- strategy:

  Character. The parallel processing strategy to use. Valid options are
  "sequential", "multisession" (default), "multicore", and "cluster".
  See
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  and
  [`ecokit::set_parallel()`](https://elgabbas.github.io/ecokit/reference/set_parallel.html)
  for details.

- future_max_size:

  Numeric. Maximum allowed total size (in megabytes) of global variables
  identified. See
  [`ecokit::set_parallel()`](https://elgabbas.github.io/ecokit/reference/set_parallel.html)
  and `future.globals.maxSize` argument of
  [`future::future.options()`](https://future.futureverse.org/reference/zzz-future.options.html)
  for more details.

- selected_species, excluded_species:

  Character vector or NULL. Names of species to include or exclude for
  modelling.

- env_file:

  Character. Path to a file with environment variable definitions for
  spatial datasets. Default is `".env"`.

- clamp_pred:

  Logical. Should clamping be applied to sampling efforts and river
  length predictors for prediction? Default is `TRUE`.

- fix_efforts, fix_rivers:

  Character or numeric (length 1). Method or fixed value for sampling
  effort and river length (both at log-scale) when clamping is enabled
  (`clamp_pred = TRUE`). Valid methods: "identity" (use observed, with
  no clamping), summary statistics for the sampling efforts layer
  ("median", "mean", "max", or "q90" (default; 90th percentile)), or a
  single numeric value within observed range.

- climate_models:

  Character vector or "all". Which climate change models to use for
  future projections. Valid values (case-sensitive): "GFDL-ESM4",
  "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL", or "all"
  (default, meaning all available models). If a subset, must be a subset
  of the listed valid models.

- climate_scenarios:

  Character vector or "all". Which climate change scenarios to use for
  future projections. Valid values: "ssp126", "ssp370", "ssp585", or
  "all" (default, meaning all available scenarios). If a subset, must be
  a subset of the listed valid scenarios.

- climate_periods:

  Character vector or "all". Time periods for prediction. Valid values
  are "2011-2040", "2041-2070", "2071-2100", or "all" (default), or
  subset of supported periods.

- copy_maxent_html:

  Logical. Whether to copy the directory containing HTML results from
  Maxent to the modelling directory. Default is `TRUE`.

## Value

A tibble summarizing model results for each species, including:

- Evaluation metrics for training and testing data (AUC, TSS, Kappa,
  etc.)

- Variable importance scores

- Response curves for each environmental variable

- Prediction summaries for current and future climate scenarios

- Paths to generated model files and prediction rasters

Additionally, the function saves various outputs to disk for future use:

- Fitted model objects (as .RData files)

- Extracted model information (evaluation metrics, variable importance,
  etc.)

- Prediction rasters for each species, cross-validation fold, and
  climate scenario

- Summary statistics across CV folds (mean, weighted mean, SD, and
  coefficient of variation)

- Species richness maps for each climate scenario

## Details

The `fit_sdm_models` function orchestrates a comprehensive workflow that
handles all aspects of single-species distribution modelling for
invasive alien plant species in Europe. The workflow integrates several
internal components that manage different stages of the modelling
process:

**Overall workflow:**

- *Input validation*: Checks all parameters for validity

- *Data preparation*: Loads and processes model data

- *Parallel processing setup*: Configures computational resources

- *Model fitting and prediction*: For each species and CV fold

- *Results summarization*: Compiles metrics, variable importance, and
  predictions

- *Species richness calculation*: Across all modelled species

**Core capabilities:**

- *Data preparation*: The workflow validates and prepares necessary
  input data including modelling data, environmental predictors, and
  prediction datasets. It handles species selection, data loading, and
  preprocessing of spatial predictors (including clamping of sampling
  efforts and river length when required).

- *Model parameterization*: The function provides carefully selected
  default settings for various SDM algorithms, ensuring consistent
  parameterization across models.

- *Model information extraction*: After fitting, the workflow
  automatically extracts key information from fitted SDM objects,
  including evaluation metrics, variable importance, and response
  curves.

- *Model optimization*: Technical improvements like optimizing SDM model
  object size by setting formula environments to the base environment
  address known issues in the sdm package.

- *Parallel prediction*: The workflow efficiently generates predictions
  for each species and cross-validation fold, handling model fitting,
  information extraction, prediction, and file saving in parallel.

- *Statistical summarization*: Summary statistics are calculated across
  cross-validation folds, including mean, weighted mean (by test AUC),
  standard deviation, and coefficient of variation of predictions.

## Author

Ahmed El-Gabbas
