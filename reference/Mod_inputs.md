# Prepare initial models for model fitting with Hmsc-HPC

The **`mod_prepare_hpc`** function prepares input data and initialises
models for fitting with
[Hmsc-HPC](https://doi.org/10.1371/journal.pcbi.1011914). It performs
multiple tasks, including data preparation, defining spatial block
cross-validation folds, generating Gaussian Predictive Process (GPP)
knots ([Tikhonov et al.](https://doi.org/10.1002/ecy.2929)),
initialising models, and creating HPC execution commands. The function
supports parallel processing and offers the option to include or exclude
phylogenetic tree data.  
  
The `mod_prepare_data` function is used to prepare habitat-specific data
for Hmsc models. This function processes environmental and species
presence data, reads environment variables from a file, verifies paths,
loads and filters species data based on habitat type and minimum
presence grid cells per species, and merges various environmental layers
(e.g., CHELSA Bioclimatic variables, habitat coverage, road and railway
intensity, sampling efforts) into a single dataset. Processed data is
saved to disk as an `*.RData` file.

## Usage

``` r
mod_prepare_data(
  hab_abb = NULL,
  directory_name = NULL,
  min_efforts_n_species = 100L,
  exclude_cultivated = TRUE,
  exclude_0_habitat = TRUE,
  n_pres_per_species = 80L,
  env_file = ".env",
  verbose_progress = TRUE
)

mod_prepare_hpc(
  hab_abb = NULL,
  directory_name = NULL,
  min_efforts_n_species = 100L,
  n_pres_per_species = 80L,
  env_file = ".env",
  gpp = TRUE,
  gpp_dists = NULL,
  min_lf = NULL,
  max_lf = NULL,
  alphapw = list(Prior = NULL, Min = 10, Max = 1500, Samples = 101),
  efforts_as_predictor = TRUE,
  road_rail_as_predictor = TRUE,
  habitat_as_predictor = TRUE,
  river_as_predictor = FALSE,
  soil_as_predictor = TRUE,
  wetness_as_predictor = TRUE,
  bio_variables = c("bio3", "bio4", "bio11", "bio18", "bio19", "npp"),
  quadratic_variables = c("bio4", "bio11"),
  n_species_per_grid = 0L,
  exclude_cultivated = TRUE,
  exclude_0_habitat = TRUE,
  cv_n_folds = 4L,
  cv_n_grids = 20L,
  cv_n_rows = 2L,
  cv_n_columns = 2L,
  cv_sac = FALSE,
  cv_fit = list(cv_type = NULL, cv_fold = NULL, inherit_dir = NULL),
  use_phylo_tree = TRUE,
  no_phylo_tree = FALSE,
  overwrite_rds = TRUE,
  n_cores = 8L,
  strategy = "multisession",
  future_max_size = 1000L,
  mcmc_n_chains = 4L,
  mcmc_thin = NULL,
  mcmc_samples = 1000L,
  mcmc_transient_factor = 500L,
  mcmc_verbose = 200L,
  skip_fitted = TRUE,
  n_array_jobs = 210L,
  model_country = NULL,
  verbose_progress = TRUE,
  slurm_prepare = TRUE,
  memory_per_cpu = "64G",
  job_runtime = NULL,
  job_name = NULL,
  path_hmsc = NULL,
  check_python = FALSE,
  to_json = FALSE,
  precision = 64L,
  ...
)
```

## Arguments

- hab_abb:

  Character. Abbreviation for the habitat type (based on
  [SynHab](https://www.preslia.cz/article/pdf?id=11548)) for which to
  prepare data. Valid values are `0`, `1`, `2`, `3`, `4a`, `4b`, `10`,
  `12a`, `12b`. If `hab_abb` = `0`, data is prepared irrespective of the
  habitat type. For more details, see [Pysek et
  al.](https://doi.org/10.23855/preslia.2022.447).

- directory_name:

  Character. Directory name, without its parents, where the models will
  be saved. This directory will be created.

- min_efforts_n_species:

  Integer. Minimum number of vascular plant species per grid cell (from
  GBIF data) required for inclusion in the models. This is to exclude
  grid cells with very little sampling efforts. Defaults to `100`.

- exclude_cultivated:

  Logical. Whether to exclude countries with cultivated or casual
  observations per species. Defaults to `TRUE`.

- exclude_0_habitat:

  Logical. Whether to exclude grid cells with zero percentage habitat
  coverage. Defaults to `TRUE`.

- n_pres_per_species:

  Integer. The minimum number of presence grid cells for a species to be
  included in the analysis. The number of presence grid cells per
  species is calculated after discarding grid cells with low sampling
  efforts (`min_efforts_n_species`) and zero percentage habitat coverage
  `exclude_0_habitat`. Defaults to `80`.

- env_file:

  Character. Path to the environment file containing paths to data
  sources. Defaults to `.env`.

- verbose_progress:

  Logical. Whether to print a message upon successful saving of files.
  Defaults to `FALSE`.

- gpp:

  Logical. Whether to fit spatial random effect using Gaussian
  Predictive Process. Defaults to `TRUE`. If `FALSE`, non-spatial models
  will be fitted.

- gpp_dists:

  Integer. Spacing (in kilometres) between GPP knots, as well as the
  minimum allowable distance between a knot and the nearest sampling
  point. The knots are generated using the
  [prepare_knots](https://biodt.github.io/IASDT.R/reference/prepare_knots.md)
  function, and this value is used for both `knotDist` and `minKnotDist`
  in
  [Hmsc::constructKnots](https://rdrr.io/pkg/Hmsc/man/constructKnots.html).

- min_lf, max_lf:

  Integer. Minimum and maximum number of latent factors to be used. Both
  default to `NULL` which means that the number of latent factors will
  be estimated from the data. If either is provided, the respective
  values will be used as arguments to
  [Hmsc::setPriors](https://rdrr.io/pkg/Hmsc/man/setPriors.html).

- alphapw:

  Prior for the alpha parameter. Defaults to a list with `Prior = NULL`,
  `Min = 10`, `Max = 1500`, and `Samples = 101`. If `alphapw` is `NULL`
  or a list with all `NULL` list items, the default prior will be used.
  If `Prior` is a matrix, it will be used as the prior. If
  `Prior = NULL`, the prior will be generated using the minimum and
  maximum values of the alpha parameter (`min` and `max`, respectively;
  in kilometre) and the number of samples (`Samples`). Defaults to a
  prior with 101 samples ranging from 10 to 1500 km, with the first
  value in the second column set to 0.5.

- efforts_as_predictor:

  Logical. Whether to include the (log₁₀) sampling efforts as predictor
  to the model. Default: `TRUE`.

- road_rail_as_predictor:

  Logical. Whether to include the (log₁₀) sum of road and railway
  intensity as predictor to the model. Default: `TRUE`.

- habitat_as_predictor:

  Logical. Whether to include the (log₁₀) percentage coverage of
  respective habitat type per grid cell as predictor to the model.
  Default: `TRUE`. Only valid if `hab_abb` not equals to `0`.

- river_as_predictor:

  Logical. Whether to include the (log₁₀) total length of rivers per
  grid cell as predictor to the model. Default: `FALSE`. See
  [river_length](https://biodt.github.io/IASDT.R/reference/River_Length.md)
  for more details.

- soil_as_predictor:

  Logical. Whether to include soil bulk density at depth of 15-30 cm as
  predictor to the model. Default: `TRUE`. See
  [soil_density_process](https://biodt.github.io/IASDT.R/reference/soil_density_process.md).

- wetness_as_predictor:

  Logical. Whether to include topographic wetness index as predictor to
  the model. Default: `TRUE`. See
  [wetness_index_process](https://biodt.github.io/IASDT.R/reference/wetness_index_process.md).

- bio_variables:

  Character vector. Variables from CHELSA (bioclimatic variables
  (bio1-bio19) and additional predictors (e.g., Net Primary
  Productivity, npp)) to be used in the model. By default, six
  ecologically relevant and minimally correlated variables are selected:
  `c("bio3", "bio4", "bio11", "bio18", "bio19", "npp")`.

- quadratic_variables:

  Character vector for variables for which quadratic terms are used.
  Defaults to `c("bio4", "bio11")`. If `quadratic_variables` is `NULL`,
  no quadratic terms will be used.

- n_species_per_grid:

  Integer. Minimum number of species required for a grid cell to be
  included in the analysis. This filtering occurs after applying
  `min_efforts_n_species` (sampling effort thresholds),
  `n_pres_per_species` (minimum species presence thresholds), and
  `exclude_0_habitat` (exclude 0% habitat coverage). Default (0):
  Includes all grid cells. Positive value (\>0): Includes only grid
  cells where at least `n_species_per_grid` species are present.

- cv_n_folds:

  Integer. Number of cross-validation folds. Default: 4L.

- cv_n_grids:

  Integer. For `cv_dist` cross-validation strategy (see
  [mod_cv_prepare](https://biodt.github.io/IASDT.R/reference/mod_CV_prepare.md)),
  this argument determines the size of the blocks (how many grid cells
  in both directions).

- cv_n_rows, cv_n_columns:

  Integer. Number of rows and columns used in the `cv_large`
  cross-validation strategy (see
  [mod_cv_prepare](https://biodt.github.io/IASDT.R/reference/mod_CV_prepare.md)),
  in which the study area is divided into large blocks given the
  provided `cv_n_rows` and `cv_n_columns` values. Both default to 2
  which means to split the study area into four large blocks at the
  median latitude and longitude.

- cv_sac:

  Logical. Whether to use the spatial autocorrelation to determine the
  block size. Defaults to `FALSE`,

- cv_fit:

  A list with three elements determining if the current model is for a
  specific cross-validation fold or for full dataset.

  - `cv_type` (character): the type of cross-validation to use. Valid
    options are "cv_dist", "cv_large", and "cv_sac". Default: `NULL`,
    which means fit models on the full dataset. This can not be `NULL`
    if `cv_fold` is provided.

  - `cv_fold` (integer): the id of the cross-validation fold to fit. For
    example, `cv_fold = 4` means use the fourth fold for testing.
    Default: `NULL` which means no cross-validation is performed. This
    can not be `NULL` if `cv_type` is provided.

  - `inherit_dir` (character): the name of the directory (without its
    parents) to inherit (copy) species and cross-validation data from.
    Defaults to `NULL`, which means that a data on species and
    cross-validation will be calculated.

- use_phylo_tree, no_phylo_tree:

  Logical. Whether to fit models with (use_phylo_tree) or without
  (no_phylo_tree) phylogenetic trees. Defaults are
  `use_phylo_tree = TRUE` and `no_phylo_tree = FALSE`, meaning only
  models with phylogenetic trees are fitted by default. At least one of
  `use_phylo_tree` and `no_phylo_tree` should be `TRUE`.

- overwrite_rds:

  Logical. Whether to overwrite previously exported RDS files for
  initial models. Default: `TRUE`.

- n_cores:

  Integer. Number of CPU cores to use for parallel processing. Default:
  8.

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

- mcmc_n_chains:

  Integer. Number of model chains. Default: 4L.

- mcmc_thin:

  Integer vector. Thinning value(s) in MCMC sampling. If more than one
  value is provided, a separate model will be fitted at each value of
  thinning.

- mcmc_samples:

  Integer vector. Value(s) for the number of MCMC samples. If more than
  one value is provided, a separate model will be fitted at each value
  of number of samples. Defaults to 1000.

- mcmc_transient_factor:

  Integer. Transient multiplication factor. The value of `transient`
  will equal the multiplication of `mcmc_transient_factor` and
  `mcmc_thin`. Default: 500.

- mcmc_verbose:

  Integer. Interval at which MCMC sampling progress is reported.
  Default: `200`.

- skip_fitted:

  Logical. Whether to skip already fitted models. Default: `TRUE`.

- n_array_jobs:

  Integer. Number of jobs per SLURM script file. In LUMI HPC, there is a
  limit of 210 submitted jobs per user for the `small-g` partition. This
  argument is used to split the jobs into multiple SLURM scripts if
  needed. Default: 210. See [LUMI
  documentation](https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/partitions)
  for more details.

- model_country:

  Character. Country or countries to filter observations by. Default:
  `NULL`, which means prepare data for the whole Europe.

- slurm_prepare:

  Logical. Whether to prepare SLURM command files. If `TRUE` (default),
  the SLURM commands will be saved to disk using the
  [mod_slurm](https://biodt.github.io/IASDT.R/reference/Mod_SLURM.md)
  function.

- memory_per_cpu:

  Character. Memory per CPU for the SLURM job. This value will be
  assigned to the `#SBATCH --mem-per-cpu=` SLURM argument. Example:
  "32G" to request 32 gigabyte. Only effective if
  `slurm_prepare = TRUE`. Defaults to "64G".

- job_runtime:

  Character. Requested time for each job in the SLURM bash arrays.
  Example: "01:00:00" to request an hour. Only effective if
  `slurm_prepare = TRUE`.

- job_name:

  Character. Name of the submitted job(s) for SLURM. If `NULL`
  (Default), the job name will be prepared based on the folder path and
  the `hab_abb` value. Only effective if `slurm_prepare = TRUE`.

- path_hmsc:

  Character. Directory path to `Hmsc-HPC` extension installation. This
  will be provided as the `path_hmsc` argument of the
  [mod_slurm](https://biodt.github.io/IASDT.R/reference/Mod_SLURM.md)
  function.

- check_python:

  Logical. Whether to check if the Python executable exists.

- to_json:

  Logical. Whether to convert unfitted models to JSON before saving to
  RDS file. Default: `FALSE`.

- precision:

  Integer. Must be either 32 or 64 (default). Defines the floating-point
  precision mode for `Hmsc-HPC` sampling (–fp 32 or –fp 64).

- ...:

  Additional parameters provided to the
  [mod_slurm](https://biodt.github.io/IASDT.R/reference/Mod_SLURM.md)
  function.

## Author

Ahmed El-Gabbas
