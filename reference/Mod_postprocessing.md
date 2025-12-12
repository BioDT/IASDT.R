# Model pipeline for post-processing fitted Hmsc models

These functions post-process fitted Hmsc models on both CPU and GPU. The
main functions in the pipeline includes `mod_postprocess_1_cpu`,
`mod_prepare_tf`, and `mod_postprocess_2_cpu` for full models without
cross-validation, as well as `mod_postprocess_cv_1_cpu` and
`mod_postprocess_cv_2_cpu` for cross-validated models. See details for
more information.

## Usage

``` r
mod_postprocess_1_cpu(
  model_dir = NULL,
  hab_abb = NULL,
  strategy = "multisession",
  future_max_size = 1500L,
  n_cores = 8L,
  n_cores_pred = n_cores,
  n_cores_lf = n_cores,
  n_cores_vp = n_cores,
  env_file = ".env",
  path_hmsc = NULL,
  memory_per_cpu = "64G",
  job_runtime = "01:00:00",
  from_json = FALSE,
  gpp_dist = NULL,
  use_trees = "tree",
  mcmc_n_samples = 1000L,
  mcmc_thin = NULL,
  n_omega = 1000L,
  cv_name = c("cv_dist", "cv_large"),
  n_grid = 50L,
  use_tf = TRUE,
  tf_use_single = FALSE,
  lf_temp_cleanup = TRUE,
  lf_check = FALSE,
  temp_cleanup = TRUE,
  tf_environ = NULL,
  pred_new_sites = TRUE,
  width_omega = 26,
  height_omega = 22.5,
  width_beta = 25,
  height_beta = 35,
  spatial_model = TRUE,
  tar_predictions = TRUE,
  plot_predictions = TRUE,
  is_cv_model = FALSE,
  clamp_pred = TRUE,
  fix_efforts = "q90",
  fix_rivers = "q90",
  climate_models = c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0",
    "UKESM1-0-LL"),
  climate_scenario = c("ssp126", "ssp370", "ssp585")
)

mod_prepare_tf(
  process_vp = TRUE,
  process_lf = TRUE,
  n_batch_files = 210L,
  env_file = ".env",
  working_directory = NULL,
  partition_name = "small-g",
  lf_runtime = "01:00:00",
  model_prefix = NULL,
  vp_runtime = "02:00:00"
)

mod_postprocess_2_cpu(
  model_dir = NULL,
  hab_abb = NULL,
  strategy = "multisession",
  future_max_size = 1500L,
  n_cores = 8L,
  n_cores_pred = n_cores,
  n_cores_lf = n_cores,
  n_cores_rc = n_cores,
  n_cores_vp = n_cores,
  env_file = ".env",
  gpp_dist = NULL,
  use_trees = "tree",
  mcmc_n_samples = 1000L,
  mcmc_thin = NULL,
  use_tf = TRUE,
  tf_environ = NULL,
  tf_use_single = FALSE,
  lf_check = FALSE,
  lf_temp_cleanup = TRUE,
  temp_cleanup = TRUE,
  n_grid = 50L,
  climate_models = c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0",
    "UKESM1-0-LL"),
  climate_scenario = c("ssp126", "ssp370", "ssp585"),
  clamp_pred = TRUE,
  fix_efforts = "q90",
  fix_rivers = "q90",
  pred_new_sites = TRUE,
  tar_predictions = TRUE,
  rc_prepare = TRUE,
  rc_plot = TRUE,
  vp_prepare = TRUE,
  vp_plot = TRUE,
  predict_suitability = TRUE,
  plot_predictions = TRUE,
  plot_lf = TRUE,
  plot_internal_evaluation = TRUE,
  spatial_model = TRUE,
  is_cv_model = FALSE
)

mod_postprocess_cv_1_cpu(
  model_dir = NULL,
  cv_names = NULL,
  n_cores = 8L,
  strategy = "multisession",
  env_file = ".env",
  from_json = FALSE,
  use_tf = TRUE,
  tf_use_single = FALSE,
  tf_environ = NULL,
  n_cores_lf = n_cores,
  lf_only = TRUE,
  lf_temp_cleanup = TRUE,
  lf_check = FALSE,
  lf_runtime = "01:00:00",
  temp_cleanup = TRUE,
  n_batch_files = 210L,
  working_directory = NULL,
  partition_name = "small-g"
)

mod_postprocess_cv_2_cpu(
  model_dir = NULL,
  cv_names = NULL,
  n_cores = 8L,
  strategy = "multisession",
  env_file = ".env",
  use_tf = TRUE,
  tf_use_single = FALSE,
  temp_cleanup = TRUE,
  lf_temp_cleanup = TRUE,
  tf_environ = NULL,
  n_cores_lf = n_cores,
  lf_check = FALSE
)
```

## Arguments

- model_dir:

  Character. Path to the root directory of the fitted model.

- hab_abb:

  Character. Habitat abbreviation indicating the specific
  [SynHab](https://www.preslia.cz/article/pdf?id=11548) habitat type.
  Valid values: `0`, `1`, `2`, `3`, `4a`, `4b`, `10`, `12a`, `12b`. See
  [Pysek et al.](https://doi.org/10.23855/preslia.2022.447) for details.

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
  identified. See `future.globals.maxSize` argument of
  [future::future.options](https://future.futureverse.org/reference/zzz-future.options.html)
  for more details.

- n_cores, n_cores_pred, n_cores_lf, n_cores_vp, n_cores_rc:

  Integer. Number of cores to use for parallel processing. They are used
  for different processing steps: `n_cores` for merging chains and
  plotting convergence convergence diagnostics; `n_cores_pred` for
  predicting species' habitat suitability; `n_cores_lf` for predicting
  latent factors; `n_cores_vp` for processing variance partitioning; and
  `n_cores_rc` for response curve prediction. All default to `8L`. If
  `strategy = "sequential"`, all of these arguments are set to `1L`.

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

- path_hmsc:

  Character. Path to the Hmsc-HPC installation.

- memory_per_cpu:

  Character. Memory allocation per CPU core. Example: "32G" for 32
  gigabytes. Defaults to "64G".

- job_runtime:

  Character. Maximum allowed runtime for jobs for refitting the models
  (if needed) and cross validating models. Defaults to "01:00:00" for
  one hour. If not provided, the function throws an error.

- from_json:

  Logical. Whether to convert loaded models from JSON format before
  reading. Defaults to `FALSE`.

- gpp_dist:

  Integer. Distance in *kilometres* between knots for the selected
  model.

- use_trees:

  Character. Whether a phylogenetic tree was used in the selected model.
  Accepts "tree" (default) or "no_tree".

- mcmc_thin, mcmc_n_samples:

  Integer. Thinning value and the number of MCMC samples of the selected
  model.

- n_omega:

  Integer. The number of species to be sampled for the `Omega` parameter
  transformation. Defaults to 100.

- cv_name:

  `NULL` or character vector. Column name(s) in the model input data to
  be used to cross-validate the models (see
  [mod_prepare_data](https://biodt.github.io/IASDT.R/reference/Mod_inputs.md)
  and
  [mod_cv_prepare](https://biodt.github.io/IASDT.R/reference/mod_CV_prepare.md)).
  If `cv_name = NULL`, no cross-validation data preparation is done. See
  [mod_cv_fit](https://biodt.github.io/IASDT.R/reference/Mod_CV_Fit.md)
  for valid options.

- n_grid:

  Integer. Number of points along the gradient for continuous focal
  variables. Higher values result in smoother curves. Default: 50. See
  [Hmsc::constructGradient](https://rdrr.io/pkg/Hmsc/man/constructGradient.html)
  for details.

- use_tf:

  Logical. Whether to use `TensorFlow` for calculations. Defaults to
  `TRUE`.

- tf_use_single:

  Logical. Whether to use single precision for the `TensorFlow`
  calculations. Defaults to `FALSE`.

- lf_check:

  Logical. If `TRUE`, the function checks if the output files are
  already created and valid. If `FALSE`, the function will only check if
  the files exist without checking their integrity. Default is `FALSE`.

- temp_cleanup, lf_temp_cleanup:

  Logical. Whether to delete temporary files after finishing predicting
  latent factor or species distribution. Default: `TRUE`.

- tf_environ:

  Character. Path to the Python environment. This argument is required
  if `use_tf` is `TRUE` under Windows. Defaults to `NULL`.

- pred_new_sites:

  Logical. Whether to predict suitability at new sites. Default: `TRUE`.

- width_omega, height_omega, width_beta, height_beta:

  Integer. The width and height of the generated heatmaps of the Omega
  and Beta parameters in centimetres.

- spatial_model, plot_lf:

  Logical. Whether the model is spatial (`TRUE`) or not (`FALSE`) and
  whether to plot latent factors for spatial models as JPEG files (using
  [plot_latent_factor](https://biodt.github.io/IASDT.R/reference/plot_latent_factor.md)).
  Defaults to `TRUE`.

- is_cv_model:

  Logical. Whether the model is a cross-validated model (`TRUE`) or
  fitted with the full dataset (`FALSE`; default). If `TRUE`, the
  explanatory and predictive power of the model will be computed.

- clamp_pred:

  Logical indicating whether to clamp the sampling efforts at a single
  value. If `TRUE` (default), the `fix_efforts` argument must be
  provided.

- fix_efforts:

  Numeric or character. When `clamp_pred = TRUE`, fixes the sampling
  efforts predictor at this value during predictions. If numeric, uses
  the value directly (on log₁₀ scale). If character, must be one of
  `identity` (i.e., do not fix), `median`, `mean`, `max`, or `q90` (90%
  quantile). Using `max` may reflect extreme sampling efforts from
  highly sampled locations, while `q90` captures high sampling areas
  without extremes. Required if `clamp_pred = TRUE`.

- fix_rivers:

  Numeric, character, or `NULL`. Similar to `fix_efforts`, but for the
  river length predictor. If `NULL`, the river length is not fixed.
  Default: `q90`.

- climate_models:

  Character vector. Climate models for future predictions. Available
  options are
  `c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")`
  (default).

- climate_scenario:

  Character vector. Climate scenarios for future predictions. Available
  options are: `c("ssp126", "ssp370", "ssp585")` (default).

- process_vp, process_lf:

  Logical. Whether to prepares batch scripts for variance partitioning
  computations and latent factor predictions on GPUs. Defaults to
  `TRUE`.

- n_batch_files:

  Integer. Number of output batch files to create. Must be less than or
  equal to the maximum job limit of the HPC environment.

- working_directory:

  Character. Optionally sets the working directory in batch scripts to
  this path. If `NULL`, the directory remains unchanged.

- partition_name:

  Character. Name of the partition to submit the SLURM jobs to. Default
  is `small-g`.

- lf_runtime, vp_runtime:

  Character. Time limit for latent factor prediction and variance
  partitioning processing jobs, respectively. Defaults are `01:00:00`
  and `02:00:00` respectively.

- model_prefix:

  Character. Prefix for the model name. A directory named
  `model_prefix_TF` is created in the `model_dir` to store the
  `TensorFlow` running commands. Defaults to `NULL`. This can not be
  `NULL`.

- rc_prepare, rc_plot:

  Logical. Whether to prepare the data for response curve prediction
  (using
  [rc_prepare_data](https://biodt.github.io/IASDT.R/reference/Response_curves.md))
  and plot the response curves as JPEG files. (using
  [rc_plot_sr](https://biodt.github.io/IASDT.R/reference/Response_curves.md),
  [rc_plot_species](https://biodt.github.io/IASDT.R/reference/Response_curves.md),
  and
  [rc_plot_species_all](https://biodt.github.io/IASDT.R/reference/Response_curves.md)).
  Defaults to `TRUE`.

- vp_prepare, vp_plot:

  Logical. Whether to prepare the data for variance partitioning (using
  [variance_partitioning_compute](https://biodt.github.io/IASDT.R/reference/Variance_partitioning.md))
  and plot its results (using
  [variance_partitioning_plot](https://biodt.github.io/IASDT.R/reference/Variance_partitioning.md)).
  Defaults to `TRUE`.

- predict_suitability, tar_predictions, plot_predictions:

  Logical. Whether to predict habitat suitability across different
  climate options (using
  [predict_maps](https://biodt.github.io/IASDT.R/reference/Predict_Maps.md)),
  compress the resulted files into a single `*.tar` file (without
  compression), or to plot species and species richness predictions as
  JPEG files (using
  [plot_prediction](https://biodt.github.io/IASDT.R/reference/plot_prediction.md)).
  Defaults to `TRUE`.

- plot_internal_evaluation:

  Logical. Whether to compute and visualise model internal evaluation
  (explanatory power) using
  [plot_evaluation](https://biodt.github.io/IASDT.R/reference/plot_evaluation.md).
  Defaults to `TRUE`.

- cv_names:

  Character vector. Names of cross-validation strategies to merge,
  matching those used during model setup. Defaults to
  `c("cv_dist", "cv_large")`. The names should be one of `cv_dist`,
  `cv_large`, or `cv_sac`. Applies only to `mod_merge_chains_cv`.

- lf_only:

  Logical. Whether to predict only the latent factor. This is useful for
  distributing processing load between GPU and CPU. When
  `lf_only = TRUE`, latent factor prediction needs to be computed
  separately on GPU. When computations are finished on GPU, the function
  can later be rerun with `lf_only = FALSE` (default) to predict habitat
  suitability using the already-computed latent factor predictions.

## Details

**mod_postprocess_1_cpu**

This function performs the initial post-processing step for
habitat-specific fitted models, automating the following tasks:

- check unsuccessful models:
  [mod_slurm_refit](https://biodt.github.io/IASDT.R/reference/Mod_SLURM.md)

- merge chains and save R objects (fitted model object and coda object)
  to `qs2` or `RData` files:
  [mod_merge_chains](https://biodt.github.io/IASDT.R/reference/Mod_Merge_Chains.md)

- visualise the convergence of all model variants fitted
  [convergence_plot_all](https://biodt.github.io/IASDT.R/reference/Convergence_Plot_All.md)

- visualise the convergence of selected model, including plotting
  Gelman-Rubin-Brooks
  [plot_gelman](https://biodt.github.io/IASDT.R/reference/plot_gelman.md)
  and
  [convergence_plot](https://biodt.github.io/IASDT.R/reference/Convergence_plots.md)
  for model convergence diagnostics of the `rho`, `alpha`, `omega`, and
  `beta` parameters.

- extract and save model summary:
  [mod_summary](https://biodt.github.io/IASDT.R/reference/Mod_Summary.md)

- plotting model parameters:
  [mod_heatmap_omega](https://biodt.github.io/IASDT.R/reference/Parameter_Heatmap.md),
  [mod_heatmap_beta](https://biodt.github.io/IASDT.R/reference/Parameter_Heatmap.md)

- prepare data for cross-validation and fit initial cross-validated
  models:
  [mod_cv_fit](https://biodt.github.io/IASDT.R/reference/Mod_CV_Fit.md)

- Prepare scripts for GPU processing, including:

  - predicting latent factors of the response curves:
    [rc_prepare_data](https://biodt.github.io/IASDT.R/reference/Response_curves.md)

  - predicting latent factors for new sampling units:
    [predict_maps](https://biodt.github.io/IASDT.R/reference/Predict_Maps.md)

  - computing variance partitioning:
    [variance_partitioning_compute](https://biodt.github.io/IASDT.R/reference/Variance_partitioning.md)

  

**mod_prepare_tf**

After running `mod_postprocess_1_cpu` for all habitat types, this
function prepares batch scripts for GPU computations of all habitat
types:

- for *variance partitioning*, the function matches all files with the
  pattern ` "vp_.+command.txt"` (created by
  [variance_partitioning_compute](https://biodt.github.io/IASDT.R/reference/Variance_partitioning.md)
  and merges their contents into a single file
  (`model_prefix_TF/vp_commands.txt`). Then, it prepares a SLURM script
  for variance partitioning computations
  (`model_prefix_TF/vp_slurm.slurm`).

- for *latent factor predictions*, the function matches all files with
  the pattern `"^lf_new_sites_commands_.+.txt|^lf_rc_commands_.+txt"`
  and split their contents into multiple scripts at the
  `model_prefix_TF` directory for processing as a batch job. The
  function prepares a SLURM script for latent factor predictions
  (`lf_slurm.slurm`).

This function is tailored for the LUMI HPC environment and assumes that
the `tensorflow` module is installed and correctly configured with all
required Python packages. On other HPC systems, users may need to modify
the function to load a Python virtual environment or install the
required dependencies for `TensorFlow` and related packages.

  
  

**mod_postprocess_2_cpu**

This function continues running the analysis pipeline for
post-processing Hmsc by automating the following steps:

- process and visualise response curves:
  [response_curves](https://biodt.github.io/IASDT.R/reference/Response_curves.md)

- predict habitat suitability across different climate options:
  [predict_maps](https://biodt.github.io/IASDT.R/reference/Predict_Maps.md)

- plot species & SR predictions as JPEG:
  [plot_prediction](https://biodt.github.io/IASDT.R/reference/plot_prediction.md)

- plot latent factors as JPEG:
  [plot_latent_factor](https://biodt.github.io/IASDT.R/reference/plot_latent_factor.md)

- process and visualise variance partitioning:
  [variance_partitioning_compute](https://biodt.github.io/IASDT.R/reference/Variance_partitioning.md)
  and
  [variance_partitioning_plot](https://biodt.github.io/IASDT.R/reference/Variance_partitioning.md)

- compute and visualizing model internal evaluation (explanatory power):
  [plot_evaluation](https://biodt.github.io/IASDT.R/reference/plot_evaluation.md)

- initiate post-processing of fitted cross-validated models: prepare
  commands for latent factor predictions on GPU — **Ongoing**

This function should be run after:

- completing `mod_postprocess_1_cpu` and `mod_prepare_tf` on CPU,

- running `vp_slurm.slurm` and `lf_slurm.slurm` on GPU to process
  response curves and latent factor predictions (both scripts are
  generated by `mod_prepare_tf`).

- submitting SLURM jobs for cross-validated model fitting.

  

**mod_postprocess_cv_1_cpu**

This function is similar to `mod_postprocess_1_cpu`, but it is
specifically designed for cross-validated models. It automates merging
fitted cross-validated model chains into `Hmsc` model objects and
prepare scripts for latent factor prediction on `TensorFlow` using
[predict_maps_cv](https://biodt.github.io/IASDT.R/reference/Predict_Maps.md).

  
  

**mod_postprocess_cv_2_cpu**

The function 1) processes `*.feather` files resulted from Latent Factor
predictions (using `TensorFlow`) and saves LF predication to disk; 2)
predicts species-specific mean habitat suitability at testing
cross-validation folds and calculates testing evaluation metrics; 3)
generates plots of the evaluation metrics.

## Author

Ahmed El-Gabbas
