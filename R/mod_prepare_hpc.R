## |------------------------------------------------------------------------| #
# mod_prepare_hpc ----
## |------------------------------------------------------------------------| #

#' Prepare initial models for model fitting with Hmsc-HPC
#'
#' The **`mod_prepare_hpc`** function prepares input data and initialises models
#' for fitting with [Hmsc-HPC](https://doi.org/10.1371/journal.pcbi.1011914). It
#' performs multiple tasks, including data preparation, defining spatial block
#' cross-validation folds, generating Gaussian Predictive Process (GPP) knots
#' ([Tikhonov et al.](https://doi.org/10.1002/ecy.2929)), initialising models,
#' and creating HPC execution commands. The function supports parallel
#' processing and offers the option to include or exclude phylogenetic tree
#' data.<br/><br/> The `mod_prepare_data` function is used to prepare
#' habitat-specific data for Hmsc models. This function processes environmental
#' and species presence data, reads environment variables from a file, verifies
#' paths, loads and filters species data based on habitat type and minimum
#' presence grid cells per species, and merges various environmental layers
#' (e.g., CHELSA Bioclimatic variables, habitat coverage, road and railway
#' intensity, sampling efforts) into a single dataset. Processed data is saved
#' to disk as an `*.RData` file.
#' @param directory_name Character. Directory name, without its parents, where
#'   the models will be saved. This directory will be created.
#' @param gpp Logical. Whether to fit spatial random effect using Gaussian
#'   Predictive Process. Defaults to `TRUE`. If `FALSE`, non-spatial models will
#'   be fitted.
#' @param gpp_dists Integer. Spacing (in kilometres) between GPP knots, as well
#'   as the minimum allowable distance between a knot and the nearest sampling
#'   point. The knots are generated using the [prepare_knots] function, and this
#'   value is used for both `knotDist` and `minKnotDist` in
#'   [Hmsc::constructKnots].
#' @param efforts_as_predictor Logical. Whether to include the
#'   (log<sub>10</sub>) sampling efforts as predictor to the model. Default:
#'   `TRUE`.
#' @param road_rail_as_predictor Logical. Whether to include the
#'   (log<sub>10</sub>) sum of road and railway intensity as predictor to the
#'   model. Default: `TRUE`.
#' @param habitat_as_predictor Logical. Whether to include the
#'   (log<sub>10</sub>) percentage coverage of respective habitat type per grid
#'   cell as predictor to the model. Default: `TRUE`. Only valid if `hab_abb`
#'   not equals to `0`.
#' @param river_as_predictor Logical. Whether to include the (log<sub>10</sub>)
#'   total length of rivers per grid cell as predictor to the model. Default:
#'   `FALSE`. See [river_length] for more details.
#' @param soil_as_predictor Logical. Whether to include soil bulk density at
#'   depth of 15-30 cm as predictor to the model. Default: `TRUE`. See
#'   [soil_density_process].
#' @param wetness_as_predictor Logical. Whether to include topographic wetness
#'   index as predictor to the model. Default: `TRUE`. See
#'   [wetness_index_process].
#' @param bio_variables Character vector. Variables from CHELSA (bioclimatic
#'   variables (bio1-bio19) and additional predictors (e.g., Net Primary
#'   Productivity, npp)) to be used in the model. By default, six ecologically
#'   relevant and minimally correlated variables are selected: `c("bio3",
#'   "bio4", "bio11", "bio18", "bio19", "npp")`.
#' @param quadratic_variables Character vector for variables for which quadratic
#'   terms are used. Defaults to `c("bio4", "bio11")`. If `quadratic_variables`
#'   is `NULL`, no quadratic terms will be used.
#' @param n_species_per_grid Integer. Minimum number of species required for a
#'   grid cell to be included in the analysis. This filtering occurs after
#'   applying `min_efforts_n_species` (sampling effort thresholds),
#'   `n_pres_per_species` (minimum species presence thresholds), and
#'   `exclude_0_habitat` (exclude 0% habitat coverage). Default (0): Includes
#'   all grid cells. Positive value (>0): Includes only grid cells where at
#'   least `n_species_per_grid` species are present.
#' @param cv_fit A list with three elements determining if the current model is
#'   for a specific cross-validation fold or for full dataset.
#'   - `cv_type` (character): the type of cross-validation to use. Valid
#'   options are "cv_dist", "cv_large", and "cv_sac". Default: `NULL`, which
#'   means fit models on the full dataset. This can not be `NULL` if `cv_fold`
#'   is provided.
#'  - `cv_fold` (integer): the id of the cross-validation fold to fit.
#'   For example, `cv_fold = 4` means use the fourth fold for testing. Default:
#'   `NULL` which means no cross-validation is performed. This can not be `NULL`
#'   if `cv_type` is provided.
#'   - `inherit_dir` (character): the name of the directory (without its
#'   parents) to inherit (copy) species and cross-validation data from. Defaults
#'   to `NULL`, which means that a data on species and cross-validation will be
#'   calculated.
#' @param use_phylo_tree,no_phylo_tree Logical. Whether to fit models with
#'   (use_phylo_tree) or without (no_phylo_tree) phylogenetic trees. Defaults
#'   are `use_phylo_tree = TRUE` and `no_phylo_tree = FALSE`, meaning only
#'   models with phylogenetic trees are fitted by default. At least one of
#'   `use_phylo_tree` and `no_phylo_tree` should be `TRUE`.
#' @param overwrite_rds Logical. Whether to overwrite previously exported RDS
#'   files for initial models. Default: `TRUE`.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 8.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "sequential", "multisession" (default), "multicore", and
#'   "cluster". See [future::plan()] and [ecokit::set_parallel()] for details.
#' @param mcmc_n_chains Integer. Number of model chains. Default: 4L.
#' @param mcmc_thin Integer vector. Thinning value(s) in MCMC sampling. If more
#'   than one value is provided, a separate model will be fitted at each value
#'   of thinning.
#' @param mcmc_samples Integer vector. Value(s) for the number of MCMC samples.
#'   If more than one value is provided, a separate model will be fitted at each
#'   value of number of samples. Defaults to 1000.
#' @param mcmc_transient_factor Integer. Transient multiplication factor. The
#'   value of `transient` will equal the multiplication of
#'   `mcmc_transient_factor` and `mcmc_thin`. Default: 500.
#' @param mcmc_verbose Integer. Interval at which MCMC sampling progress is
#'   reported. Default: `200`.
#' @param skip_fitted Logical. Whether to skip already fitted models. Default:
#'   `TRUE`.
#' @param n_array_jobs Integer. Number of jobs per SLURM script file. In LUMI
#'   HPC, there is a limit of 210 submitted jobs per user for the `small-g`
#'   partition. This argument is used to split the jobs into multiple SLURM
#'   scripts if needed. Default: 210. See [LUMI
#'   documentation](https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/partitions)
#'   for more details.
#' @param model_country Character. Country or countries to filter observations
#'   by. Default: `NULL`, which means prepare data for the whole Europe.
#' @param slurm_prepare Logical. Whether to prepare SLURM command files. If
#'   `TRUE` (default), the SLURM commands will be saved to disk using the
#'   [mod_slurm] function.
#' @param memory_per_cpu Character. Memory per CPU for the SLURM job. This value
#'   will be assigned to the `#SBATCH --mem-per-cpu=` SLURM argument. Example:
#'   "32G" to request 32 gigabyte. Only effective if `slurm_prepare = TRUE`.
#'   Defaults to "64G".
#' @param job_runtime Character. Requested time for each job in the SLURM bash
#'   arrays. Example: "01:00:00" to request an hour. Only effective if
#'   `slurm_prepare = TRUE`.
#' @param job_name Character. Name of the submitted job(s) for SLURM. If `NULL`
#'   (Default), the job name will be prepared based on the folder path and the
#'   `hab_abb` value. Only effective if `slurm_prepare = TRUE`.
#' @param path_hmsc  Character. Directory path to `Hmsc-HPC` extension
#'   installation. This will be provided as the `path_hmsc` argument of the
#'   [mod_slurm] function.
#' @param to_json Logical. Whether to convert unfitted models to JSON before
#'   saving to RDS file. Default: `FALSE`.
#' @param check_python Logical. Whether to check if the Python executable
#'   exists.
#' @param precision Integer. Must be either 32 or 64 (default). Defines the
#'   floating-point precision mode for `Hmsc-HPC` sampling (--fp 32 or --fp 64).
#' @param hab_abb Character. Abbreviation for the habitat type (based on
#'   [SynHab](https://www.preslia.cz/article/pdf?id=11548)) for which to prepare
#'   data. Valid values are `0`, `1`, `2`, `3`, `4a`, `4b`, `10`, `12a`, `12b`.
#'   If `hab_abb` = `0`, data is prepared irrespective of the habitat type. For
#'   more details, see [Pysek et
#'   al.](https://doi.org/10.23855/preslia.2022.447).
#' @param min_efforts_n_species Integer. Minimum number of vascular plant
#'   species per grid cell (from GBIF data) required for inclusion in the
#'   models. This is to exclude grid cells with very little sampling efforts.
#'   Defaults to `100`.
#' @param exclude_cultivated Logical. Whether to exclude countries with
#'   cultivated or casual observations per species. Defaults to `TRUE`.
#' @param exclude_0_habitat Logical. Whether to exclude grid cells with zero
#'   percentage habitat coverage. Defaults to `TRUE`.
#' @param n_pres_per_species Integer. The minimum number of presence grid cells
#'   for a species to be included in the analysis. The number of presence grid
#'   cells per species is calculated after discarding grid cells with low
#'   sampling efforts (`min_efforts_n_species`) and zero percentage habitat
#'   coverage `exclude_0_habitat`. Defaults to `80`.
#' @param cv_n_grids Integer. For `cv_dist` cross-validation strategy (see
#'   [mod_cv_prepare]), this argument determines the size of the blocks (how
#'   many grid cells in both directions).
#' @param cv_n_rows,cv_n_columns Integer. Number of rows and columns used in the
#'   `cv_large` cross-validation strategy  (see [mod_cv_prepare]), in which the
#'   study area is divided into large blocks given the provided `cv_n_rows` and
#'   `cv_n_columns` values. Both default to 2 which means to split the study
#'   area into four large blocks at the median latitude and longitude.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param verbose_progress Logical. Whether to print a message upon successful
#'   saving of files. Defaults to `FALSE`.
#' @param ... Additional parameters provided to the [mod_slurm] function.
#' @export
#' @inheritParams prepare_knots
#' @inheritParams mod_cv_prepare
#' @name mod_inputs
#' @rdname mod_inputs
#' @order 1
#' @importFrom rlang .data
#' @author Ahmed El-Gabbas

mod_prepare_hpc <- function(
    hab_abb = NULL, directory_name = NULL,
    min_efforts_n_species = 100L, n_pres_per_species = 80L, env_file = ".env",
    gpp = TRUE, gpp_dists = NULL, min_lf = NULL, max_lf = NULL,
    alphapw = list(Prior = NULL, Min = 10, Max = 1500, Samples = 101),
    efforts_as_predictor = TRUE, road_rail_as_predictor = TRUE,
    habitat_as_predictor = TRUE, river_as_predictor = FALSE,
    soil_as_predictor = TRUE, wetness_as_predictor = TRUE,
    bio_variables = c("bio3", "bio4", "bio11", "bio18", "bio19", "npp"),
    quadratic_variables = c("bio4", "bio11"),
    n_species_per_grid = 0L, exclude_cultivated = TRUE,
    exclude_0_habitat = TRUE, cv_n_folds = 4L, cv_n_grids = 20L,
    cv_n_rows = 2L, cv_n_columns = 2L, cv_sac = FALSE,
    cv_fit = list(cv_type = NULL, cv_fold = NULL, inherit_dir = NULL),
    use_phylo_tree = TRUE, no_phylo_tree = FALSE, overwrite_rds = TRUE,
    n_cores = 8L, strategy = "multisession", mcmc_n_chains = 4L,
    mcmc_thin = NULL, mcmc_samples = 1000L, mcmc_transient_factor = 500L,
    mcmc_verbose = 200L, skip_fitted = TRUE, n_array_jobs = 210L,
    model_country = NULL, verbose_progress = TRUE, slurm_prepare = TRUE,
    memory_per_cpu = "64G", job_runtime = NULL, job_name = NULL,
    path_hmsc = NULL, check_python = FALSE, to_json = FALSE, precision = 64L,
    ...) {

  .start_time <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Initial checking -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Checking input arguments", verbose = verbose_progress)

  ecokit::check_args(
    args_type = "character", args_to_check = c("directory_name", "path_hmsc"))

  logical_args <- c(
    "use_phylo_tree", "no_phylo_tree",
    "exclude_0_habitat", "skip_fitted", "verbose_progress", "to_json",
    "cv_sac", "check_python", "overwrite_rds", "slurm_prepare",
    "exclude_cultivated", "gpp", "efforts_as_predictor",
    "road_rail_as_predictor", "habitat_as_predictor", "river_as_predictor",
    "soil_as_predictor", "wetness_as_predictor")
  ecokit::check_args(args_to_check = logical_args, args_type = "logical")

  numeric_args <- c(
    "mcmc_n_chains", "mcmc_verbose", "n_array_jobs", "n_pres_per_species",
    "min_efforts_n_species", "mcmc_transient_factor", "cv_n_folds",
    "cv_n_grids", "cv_n_rows", "cv_n_columns", "precision")
  ecokit::check_args(args_to_check = numeric_args, args_type = "numeric")

  if (slurm_prepare) {
    # Validate memory_per_cpu and job_runtime
    memory_per_cpu <- .validate_slurm_ram(memory_per_cpu)
    job_runtime <- .validate_slurm_runtime(job_runtime)
  }

  # Phylogenetic tree options
  if (isFALSE(use_phylo_tree) && isFALSE(no_phylo_tree)) {
    ecokit::stop_ctx(
      "At least one of `use_phylo_tree` or `no_phylo_tree` has to be true",
      use_phylo_tree = use_phylo_tree, no_phylo_tree = no_phylo_tree,
      include_backtrace = TRUE)
  }

  rm(logical_args, numeric_args, envir = environment())

  if (!(precision %in% c(32, 64))) {
    ecokit::stop_ctx(
      "`precision` should be either of 32 or 64", precision = precision,
      include_backtrace = TRUE)
  }

  hab_abb <- .validate_hab_abb(as.character(hab_abb))
  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)

  if (!all(is.numeric(mcmc_samples)) || any(mcmc_samples <= 0)) {
    ecokit::stop_ctx(
      "`mcmc_samples` should be numeric and greater than zero",
      mcmc_samples = mcmc_samples, include_backtrace = TRUE)
  }

  if (!all(is.numeric(mcmc_thin)) || any(mcmc_thin <= 0)) {
    ecokit::stop_ctx(
      "`mcmc_thin` should be numeric and greater than zero",
      mcmc_thin = mcmc_thin, include_backtrace = TRUE)
  }

  if (min_efforts_n_species <= 0) {
    ecokit::stop_ctx(
      "`min_efforts_n_species` should be numeric and greater than zero",
      min_efforts_n_species = min_efforts_n_species, include_backtrace = TRUE)
  }

  if (n_species_per_grid < 0) {
    ecokit::stop_ctx(
      "`n_species_per_grid` has to be integer >= 0",
      n_species_per_grid = n_species_per_grid, include_backtrace = TRUE)
  }

  if (gpp) {
    if (is.null(gpp_dists)) {
      ecokit::stop_ctx(
        "`gpp_dists` can not be empty", gpp_dists = gpp_dists,
        include_backtrace = TRUE)
    }
    if (!all(is.numeric(gpp_dists)) || any(gpp_dists <= 0)) {
      ecokit::stop_ctx(
        "`gpp_dists` should be numeric and greater than zero",
        gpp_dists = gpp_dists, include_backtrace = TRUE)
    }
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  n_cells <- sp <- ias_id <- x <- y <- country <- m_thin <- rl <-
    m_name_init <- rl2 <- m_samples <- path_m_for_hpc <- m_transient <-
    m_name_fit <- chain <- post_missing <- command_hpc <- command_ws <-
    path_post <- path_mod_progress <- taxa_info_file <- path_grid <-
    eu_boundaries <- path_pa <- NAME_ENGL <- n_sp <- species_file <- file <-
    ias_id <- pa_file <- pa_model_file <- path_model <- NULL

  path_python <- fs::path(path_hmsc, "Scripts", "python.exe")

  ecokit::info_chunk(
    paste0("Preparing data for Hmsc-HPC models - ", directory_name),
    line_char = "=", verbose = verbose_progress, cat_bold = TRUE,
    line_char_rep = 70, cat_red = TRUE)

  # # |||||||||||||||||||||||||||||||||||
  # Load/check environment variables -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time(
    "Load and check environment variables", verbose = verbose_progress)

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "taxa_info_file", "DP_R_taxa_info", FALSE, TRUE,
    "eu_boundaries", "DP_R_country_boundaries", FALSE, TRUE,
    "path_grid", "DP_R_grid_processed", TRUE, FALSE,
    "path_model", "DP_R_model_root_path", FALSE, FALSE,
    "path_pa", "DP_R_pa", TRUE, FALSE)

  # Check if Python executable exists
  if (check_python && !file.exists(path_python) && Sys.info()[1] == "Windows") {
    ecokit::stop_ctx(
      "Python executable does not exist", path_python = path_python,
      include_backtrace = TRUE)
  }

  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())


  # Validate structure of cv_fit list argument
  if (!inherits(cv_fit, "list") || length(cv_fit) != 3 ||
      !setequal(names(cv_fit), c("cv_type", "cv_fold", "inherit_dir"))) {
    ecokit::stop_ctx(
      paste0(
        "`cv_fit` must be a list with exactly three elements: ",
        "`cv_type`, `cv_fold`, and `inherit_dir` (no more, no less)"),
      cv_fit = cv_fit, include_backtrace = TRUE)
  }

  cv_fit_type <- cv_fit$cv_type
  cv_fit_fold <- cv_fit$cv_fold
  cv_fit_inherit <- cv_fit$inherit_dir

  if (!is.null(cv_fit_inherit)) {
    # Validate inherited cross-validation directory name
    if (length(cv_fit_inherit) != 1 || !is.character(cv_fit_inherit) ||
        !nzchar(cv_fit_inherit)) {
      ecokit::stop_ctx(
        "Directory to inherit data from must be a non-empty character string",
        cv_fit = cv_fit, include_backtrace = TRUE)
    }

    # Validate inherited cross-validation directory
    cv_fit_inherit <- fs::path(path_model, cv_fit_inherit)
    if (!fs::dir_exists(cv_fit_inherit)) {
      ecokit::stop_ctx(
        "Directory to inherit data from does not exist",
        cv_fit = cv_fit, cv_fit_inherit = cv_fit_inherit,
        include_backtrace = TRUE)
    }

    # Validate files existence in the inherited directory
    files_to_check <- fs::path(
      cv_fit_inherit, c("model_data.RData", "cv_data.RData"))
    purrr::walk(
      files_to_check,
      ~ {
        if (!ecokit::check_data(.x, warning = FALSE)) {
          list_rdata <- list.files(cv_fit_inherit, pattern = "RData$")
          ecokit::stop_ctx(
            paste0("File does not exist in the inherited directory: ", .x),
            cv_fit = cv_fit, rdata_files = list_rdata, include_backtrace = TRUE)
        }
      })

    if (is.null(cv_fit_type) || is.null(cv_fit_fold)) {
      ecokit::stop_ctx(
        paste0(
          "`cv_fit$cv_type` and `cv_fit$cv_fold` can not be NULL",
          "if `cv_fit$inherit_dir` is provided"),
        cv_fit = cv_fit, include_backtrace = TRUE)
    }
  }

  # cv_fit_type and cv_fit_fold must be consistent
  if (
    (is.null(cv_fit_type) && !is.null(cv_fit_fold)) ||
    (!is.null(cv_fit_type) && is.null(cv_fit_fold))) {
    missing_val <- dplyr::if_else(
      is.null(cv_fit_type), "`cv_fit$cv_type`", "`cv_fit$cv_fold`")
    present_val <- dplyr::if_else(
      !is.null(cv_fit_type), "`cv_fit$cv_type`", "`cv_fit$cv_fold`")
    ecokit::stop_ctx(
      paste0(
        "Both `cv_fit$cv_type` and `cv_fit$cv_fold` must be either NULL ",
        "(both) or both non-NULL. Currently, ", missing_val,
        " is missing and ", present_val, " is provided."),
      cv_fit = cv_fit, include_backtrace = TRUE)
  }

  # Validate cv_fit_type and cv_fit_fold
  if (!is.null(cv_fit_type) && !is.null(cv_fit_fold)) {

    cv_fit_type <- .validate_cv_name(cv_fit_type)

    if (is.null(cv_fit_fold) || !is.numeric(cv_fit_fold) ||
        length(cv_fit_fold) != 1L || cv_fit_fold < 1L) {
      ecokit::stop_ctx(
        "`cv_fit_fold` must be a single positive integer",
        cv_fit = cv_fit, include_backtrace = TRUE)
    }

    if (cv_fit_fold > cv_n_folds) {
      ecokit::stop_ctx(
        "`cv_fit_fold` must be <= cv_n_folds",
        cv_fit = cv_fit, cv_n_folds = cv_n_folds, include_backtrace = TRUE)
    }

    # Flag for filtering data based on this cross-validation
    filter_cv_data <- TRUE
  } else {
    filter_cv_data <- FALSE
  }

  path_model <- fs::path(path_model, directory_name)
  if (fs::dir_exists(path_model) && length(fs::dir_ls(path_model)) > 0) {
    ecokit::stop_ctx(
      "Model directory already exists", path_model = path_model,
      include_backtrace = TRUE)
  }
  fs::dir_create(path_model)

  path_grid_r <- fs::path(path_grid, "grid_10_land_crop.RData")
  if (!file.exists(path_grid_r)) {
    ecokit::stop_ctx(
      "Path for the Europe boundaries does not exist",
      path_grid_r = path_grid_r, include_backtrace = TRUE)
  }

  # Loading country boundary data
  eu_boundaries <- ecokit::load_as(eu_boundaries) %>%
    magrittr::extract2("Bound_sf_Eur") %>%
    magrittr::extract2("L_03") %>%
    suppressWarnings()

  ecokit::record_arguments(
    out_path = fs::path(path_model, "args_mod_prepare_hpc.RData"))

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c(
      "dplyr", "sf", "Hmsc", "jsonify", "magrittr", "IASDT.R", "ecokit"),
    strategy = strategy)

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # File paths - Creating missing paths ----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time(
    "File paths - Creating missing paths", verbose = verbose_progress)
  fs::dir_create(fs::path(path_model, "init_mod_for_hpc"))
  fs::dir_create(fs::path(path_model, "model_fitting_hpc"))
  # Also create directory for SLURM outputs
  fs::dir_create(fs::path(path_model, "model_fitting_hpc", "jobs_log"))
  # Directory to save species distribution tiffs
  path_distribution <- fs::path(path_model, "species_distribution")
  fs::dir_create(path_distribution)

  path_model_data <- fs::path(path_model, "model_info.RData")

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # # Prepare list of predictors -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Prepare list of predictors", verbose = verbose_progress)

  # Check bio_variables values
  if (isFALSE(all(bio_variables %in% IASDT.R::chelsa_variables$variable))) {
    wrong_bio <- bio_variables[
      which(!(bio_variables %in% IASDT.R::chelsa_variables$variable))]
    ecokit::stop_ctx(
      "Invalid Bioclimatic variables", wrong_bio = wrong_bio,
      include_backtrace = TRUE)
  }

  x_vars <- bio_variables

  if (efforts_as_predictor) x_vars <- c(x_vars, "efforts_log")
  if (road_rail_as_predictor) x_vars <- c(x_vars, "road_rail_log")
  if (river_as_predictor) x_vars <- c(x_vars, "rivers_log")
  if (soil_as_predictor) x_vars <- c(x_vars, "soil")
  if (wetness_as_predictor) x_vars <- c(x_vars, "wetness")
  if (hab_abb != "0" && habitat_as_predictor) x_vars <- c(x_vars, "habitat_log")

  # Check quadratic variables
  if (!is.null(quadratic_variables) && !all(quadratic_variables %in% x_vars)) {
    ecokit::stop_ctx(
      "Some quadratic variables are not in the predictor variables",
      missing_vars = quadratic_variables[!quadratic_variables %in% x_vars],
      quadratic_variables = quadratic_variables, x_vars = x_vars,
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # # Preparing input data -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Preparing input data", verbose = verbose_progress)

  hab_val <- c(
    "0_all", "1_forests", "2_open_forests", "3_scrub",
    "4a_natural_grasslands", "4b_human_maintained_grasslands",
    "10_wetland", "12a_ruderal_habitats", "12b_agricultural_habitats") %>%
    stringr::str_subset(paste0("^", hab_abb, "_"))


  if (is.null(cv_fit_inherit)) {

    ecokit::info_chunk(
      "Preparing input data using IASDT.R::mod_prepare_data",
      verbose = verbose_progress, line_char_rep = 65)

    data_all <- IASDT.R::mod_prepare_data(
      hab_abb = hab_abb, directory_name = directory_name,
      min_efforts_n_species = min_efforts_n_species,
      exclude_cultivated = exclude_cultivated,
      exclude_0_habitat = exclude_0_habitat,
      n_pres_per_species = n_pres_per_species,
      env_file = env_file, verbose_progress = verbose_progress)

    ecokit::cat_sep(
      n_separators = 1L, sep_lines_before = 1L, sep_lines_after = 2L,
      verbose = verbose_progress)

  } else {

    ecokit::cat_time("Loading input data", verbose = verbose_progress)
    data_all <- ecokit::load_as(fs::path(cv_fit_inherit, "model_data.RData"))

    # copy model_data.RData
    fs::file_copy(
      fs::path(cv_fit_inherit, "model_data.RData"),
      fs::path(path_model, "model_data.RData"), overwrite = TRUE)

  }
  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Exclude grid cells with low number of presences -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time(
    "Exclude grid cells with low number of presences",
    verbose = verbose_progress)

  if (n_species_per_grid == 0) {
    ecokit::cat_time(
      paste0(
        "All grid cells, irrespective of number of species ",
        "presence, will be considered"),
      level = 1L, cat_timestamp = FALSE, verbose = verbose_progress)
  } else {
    empty_grids_id <- dplyr::select(
      data_all, tidyselect::starts_with("sp_")) %>%
      rowSums() %>%
      magrittr::is_less_than(as.integer(n_species_per_grid)) %>%
      which() %>%
      magrittr::multiply_by(-1)

    if (length(empty_grids_id) > 0) {
      ecokit::cat_time(
        paste0(
          "Excluding grid cells with < ", n_species_per_grid,
          " species presence"),
        level = 1L, verbose = verbose_progress)
      data_all <- dplyr::slice(data_all, empty_grids_id)
    }
  }

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # # Subsetting study area -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Subsetting study area", verbose = verbose_progress)

  if (!is.null(model_country)) {

    valid_countries <- model_country %in% unique(data_all$country)

    if (!all(valid_countries)) {
      ecokit::stop_ctx(
        "Invalid country names",
        invalid_countries = model_country[!valid_countries],
        include_backtrace = TRUE)
    }

    ecokit::cat_time(
      paste0(
        "Subsetting data to: ", paste(sort(model_country), collapse = " & ")),
      level = 1L, verbose = verbose_progress)

    sample_excl_sp <- dplyr::filter(data_all, country %in% model_country) %>%
      dplyr::summarise(
        dplyr::across(tidyselect::starts_with("sp_"), sum)) %>%
      tidyr::pivot_longer(
        cols = tidyselect::everything(),
        names_to = "sp", values_to = "n_cells") %>%
      dplyr::filter(n_cells < n_pres_per_species) %>%
      dplyr::pull(sp)

    ecokit::cat_time(
      paste0(length(sample_excl_sp), " species are excluded"), level = 1L,
      verbose = verbose_progress)
    data_all <- dplyr::filter(data_all, country %in% model_country) %>%
      dplyr::select(-tidyselect::all_of(sample_excl_sp))

    # # |||||||||||||||||||||||||||||||||||
    ## Plotting subsetted data -----
    # # |||||||||||||||||||||||||||||||||||

    ecokit::cat_time(
      "Plotting subsetted data", level = 1L, verbose = verbose_progress)

    n_sp_subset <- data_all %>%
      dplyr::mutate(
        n_sp = rowSums(
          dplyr::select(., tidyselect::starts_with("sp_")),
          na.rm = TRUE),
        n_sp = as.integer(n_sp)) %>%
      dplyr::select(tidyselect::all_of(c("x", "y", "n_sp"))) %>%
      sf::st_as_sf(coords = c("x", "y"), crs = 3035) %>%
      terra::rasterize(
        y = ecokit::load_as(path_grid_r, unwrap_r = TRUE), field = "n_sp") %>%
      terra::trim()

    if (hab_abb == "0") {
      habitat_column <- NULL
    } else {
      habitat_column <- c(
        "hab_1_forests", "hab_2_open_forests", "hab_3_scrub", "hab_4_grasslands",
        "hab_4a_natural_grasslands", "hab_4b_human_maintained_grasslands",
        "hab_5_sandy", "hab_6_rocky", "hab_7_dryland", "hab_8_saline",
        "hab_9_riparian", "hab_10_wetland", "hab_11_aquatic", "hab_12_man_made",
        "hab_12a_ruderal_habitats", "hab_12b_agricultural_habitats") %>%
        stringr::str_subset(paste0("_", hab_abb, "_")) %>%
        stringr::str_remove("hab_")
    }

    eu_boundaries_sub <- dplyr::filter(
      eu_boundaries, NAME_ENGL %in% model_country) %>%
      dplyr::select(tidyselect::all_of("CNTR_NAME"))

    plot_limits <- as.vector(terra::ext(n_sp_subset))
    # Relative JPEG height
    plot_dim_x <- plot_limits[2] - plot_limits[1]
    plot_dim_y <- plot_limits[4] - plot_limits[3]
    plot_height <- (plot_dim_y * 25) / (plot_dim_x * 0.95)
    n_grids <- format(nrow(data_all), big.mark = ",")
    n_sp <- length(stringr::str_subset(names(data_all), "sp_"))

    plot_caption <- paste(sort(model_country), collapse = "; ") %>%
      stringr::str_wrap(width = 110) %>%
      paste0("<strong>Selected countries</strong>: <br/>", .) %>%
      stringr::str_replace_all("\n", "<br/>")
    plot_subtitle <- paste0(
      n_sp, " IAS within \u2265", n_pres_per_species,
      " presence grid cells in the selected countries (",
      n_grids, " grid cells)")

    n_sp_per_grid_sub <- ggplot2::ggplot(environment = emptyenv()) +
      tidyterra::geom_spatraster(data = n_sp_subset) +
      tidyterra::scale_fill_whitebox_c(
        na.value = "transparent", palette = "bl_yl_rd", name = NULL) +
      ggplot2::labs(
        title = paste0(
          '<span style="color:blue; font-size:24px;"><b>',
          "Number of IAS per grid cell to be used in the models</b></span>",
          '<span style="color:black; font-size:18px;"> (',
          habitat_column, ")</span>"),
        subtitle = plot_subtitle, caption = plot_caption) +
      ggplot2::geom_sf(
        data = eu_boundaries_sub, fill = "transparent", colour = "black") +
      ggplot2::scale_y_continuous(
        expand = c(0, 0), limits = plot_limits[c(3, 4)]) +
      ggplot2::scale_x_continuous(
        expand = c(0, 0), limits = plot_limits[c(1, 2)]) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(0.1, 0, 0.1, 0, "cm"),
        plot.title = ggtext::element_markdown(
          size = 16, hjust = 0, margin = ggplot2::margin(0, 0, 0.1, 0, "cm")),
        plot.subtitle = ggplot2::element_text(
          size = 14, color = "darkgrey"),
        plot.caption = ggtext::element_markdown(
          size = 14, colour = "grey40", hjust = 0),
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.95, 0.9),
        legend.key.size = grid::unit(0.8, "cm"))

    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    ragg::agg_jpeg(
      filename = fs::path(path_model, "n_sp_per_grid_sub.jpeg"),
      width = 25, height = plot_height, res = 600, quality = 100, units = "cm")

    print(n_sp_per_grid_sub)
    grDevices::dev.off()

    rm(plot_limits, n_sp_per_grid_sub, eu_boundaries_sub, envir = environment())

  } else {
    ecokit::cat_time(
      "No data subsetting was implemented", level = 1L,
      verbose = verbose_progress)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Cross-validation ----
  # # |||||||||||||||||||||||||||||||||||

  if (is.null(cv_fit_inherit)) {

    ecokit::cat_time(
      "Prepare cross-validation folds", verbose = verbose_progress)

    data_all <- IASDT.R::mod_cv_prepare(
      input_data = data_all, env_file = env_file, x_vars = x_vars,
      cv_n_folds = cv_n_folds, cv_n_grids = cv_n_grids, cv_n_rows = cv_n_rows,
      cv_n_columns = cv_n_columns, out_path = path_model, cv_sac = cv_sac)

  } else {

    ecokit::cat_time(
      "Loading cross-validation data", verbose = verbose_progress)
    cv_data <- ecokit::load_as(fs::path(cv_fit_inherit, "cv_data.RData"))

    # copy cv_blocks.pdf and cv_data.RData files
    fs::file_copy(
      fs::path(cv_fit_inherit, "cv_blocks.pdf"),
      fs::path(path_model, "cv_blocks.pdf"), overwrite = TRUE)
    fs::file_copy(
      fs::path(cv_fit_inherit, "cv_data.RData"),
      fs::path(path_model, "cv_data.RData"), overwrite = TRUE)

    # Add cross-validation information
    if (!is.null(cv_data$cv_sac)) {
      data_all$cv_sac <- magrittr::extract2(cv_data$cv_sac, "folds_ids")
    }

    data_all$cv_dist <- magrittr::extract2(cv_data$cv_dist, "folds_ids")
    data_all$cv_large <- magrittr::extract2(cv_data$cv_large, "folds_ids")

  }

  # # |||||||||||||||||||||||||||||||||||
  # Filter data if CV info is provided ------
  # # |||||||||||||||||||||||||||||||||||

  if (filter_cv_data) {
    ecokit::cat_time(
      paste0(
        "Filtering data for testing cross-validation fold: ",
        cv_fit_type, "==", cv_fit_fold),
      verbose = verbose_progress)

    data_testing <- dplyr::filter(
      data_all, !!rlang::sym(cv_fit_type) == cv_fit_fold)
    data_all <- dplyr::filter(data_all, !!rlang::sym(cv_fit_type) != cv_fit_fold)

    # Save cross-validation training and testing datasets for model fitting
    ecokit::save_as(
      object = data_all, object_name = "model_data_training",
      out_path = fs::path(path_model, "model_data_training.RData"))
    ecokit::save_as(
      object = data_testing, object_name = "model_data_testing",
      out_path = fs::path(path_model, "model_data_testing.RData"))

    # plot training and testing data ------
    calculate_n_species_raster <- function(data) {
      species_names <- stringr::str_subset(names(data), "sp_")
      data %>%
        dplyr::select(x, y, tidyselect::starts_with("sp_")) %>%
        sf::st_as_sf(coords = c("x", "y"), crs = 3035) %>%
        terra::rasterize(
          ecokit::load_as(path_grid_r, unwrap_r = TRUE),
          field = species_names) %>%
        terra::app(sum)
    }

    n_species_all <- c(
      stats::setNames(
        calculate_n_species_raster(data_all),
        paste0(
          "**Training data** --- folds ",
          toString(setdiff(seq_len(cv_n_folds), cv_fit_fold)))),
      stats::setNames(
        calculate_n_species_raster(data_testing),
        paste0("**Testing data** --- fold ", cv_fit_fold)))

    # cross-validation folds as sf
    cv_blocks <- fs::path(path_model, "cv_data.RData")
    if (!fs::file_exists(cv_blocks)) {
      ecokit::stop_ctx(
        "CV_blocks file does not exist", cv_blocks = cv_blocks,
        include_backtrace = TRUE)
    }
    cv_blocks <- ecokit::load_as(cv_blocks) %>%
      magrittr::extract2(cv_fit_type) %>%
      magrittr::extract2("blocks")

    n_species_plot <- ggplot2::ggplot(environment = emptyenv()) +
      ggplot2::geom_sf(
        data = eu_boundaries, fill = "gray98", colour = "black",
        linewidth = 0.15) +
      tidyterra::geom_spatraster(data = n_species_all) +
      ggplot2::facet_wrap(~lyr) +
      ggplot2::geom_sf(
        data = cv_blocks, fill = "transparent",
        colour = "gray85", linewidth = 0.175) +
      tidyterra::scale_fill_whitebox_c(
        na.value = "transparent", palette = "bl_yl_rd", name = NULL) +
      ggplot2::scale_y_continuous(
        expand = c(0, 0), limits = c(1450000, 5420000)) +
      ggplot2::scale_x_continuous(
        expand = c(0, 0), limits = c(2600000, 6550000)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.border = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
        legend.position = "inside",
        legend.position.inside = c(0.97, 0.85),
        legend.key.size = grid::unit(0.5, "cm"),
        strip.text = ggtext::element_markdown(size = 12),
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank())

    ragg::agg_jpeg(
      filename = fs::path(path_model, "training_testing_data.jpeg"),
      width = 25.5, height = 13.5, res = 600, quality = 100, units = "cm")
    print(n_species_plot)
    grDevices::dev.off()

    rm(
      data_testing, n_species_all, cv_blocks, n_species_plot,
      envir = environment())

  } else {

    ecokit::cat_time(
      "No cross-validation filtering was made", verbose = verbose_progress)

  }

  invisible(gc())

  # cross-validation data to be saved
  data_cv <- data_all %>%
    dplyr::select(
      tidyselect::all_of(c("cell_num", "CellCode", "country")),
      tidyselect::starts_with("cv"))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Response - Y matrix ----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Response - Y matrix", verbose = verbose_progress)
  data_y <- dplyr::select(data_all, tidyselect::starts_with("sp_")) %>%
    as.data.frame()
  ecokit::cat_time(
    paste0(ncol(data_y), " species"), level = 1L, cat_timestamp = FALSE,
    verbose = verbose_progress)

  # Prepare species info summary only if no inherited directory is used

  if (is.null(cv_fit_inherit)) {
    ecokit::cat_time(
      "Save species summary", level = 1L, cat_timestamp = FALSE,
      verbose = verbose_progress)
    sp_summary <- fs::path(path_pa, "sp_pa_summary_df.RData")
    if (!file.exists(sp_summary)) {
      ecokit::stop_ctx(
        "sp_summary file does not exist", sp_summary = sp_summary,
        include_backtrace = TRUE)
    }
    sp_summary <- ecokit::load_as(sp_summary) %>%
      dplyr::arrange(ias_id) %>%
      dplyr::mutate(
        ias_id = stringr::str_pad(ias_id, pad = "0", width = 4),
        ias_id = paste0("sp_", ias_id)) %>%
      dplyr::filter(ias_id %in% names(data_y))

    save(sp_summary, file = fs::path(path_model, "sp_summary.RData"))

  } else {

    fs::file_copy(
      fs::path(cv_fit_inherit, "sp_summary.RData"),
      fs::path(path_model, "sp_summary.RData"), overwrite = TRUE)

  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Xformula -----
  # # |||||||||||||||||||||||||||||||||||

  # The formula object becomes too large (up to > 2GB!) if created within a
  # function. Setting the environment of the formula as an empty environment
  # release this unnecessary size. https://stackoverflow.com/questions/66241212

  ecokit::cat_time("Xformula", verbose = verbose_progress)

  if (is.null(quadratic_variables)) {
    form_vars <- x_vars
    ecokit::cat_time(
      paste0(
        "Models will be fitted using ", length(x_vars), " predictors: ",
        paste(x_vars, collapse = " + ")), level = 1L,
      cat_timestamp = FALSE, verbose = verbose_progress)
  } else {
    OnlyLinear <- setdiff(x_vars, quadratic_variables)
    form_vars <- c(
      OnlyLinear,
      paste0("stats::poly(", quadratic_variables, ", degree = 2, raw = TRUE)"))

    ecokit::cat_time(
      "Models will be fitted using:", level = 1L, cat_timestamp = FALSE,
      verbose = verbose_progress)

    ecokit::cat_time(
      paste0(length(OnlyLinear), " linear effect: "),
      level = 2L, cat_timestamp = FALSE, verbose = verbose_progress)
    ecokit::cat_time(
      paste(OnlyLinear, collapse = " + "), level = 3L, cat_timestamp = FALSE,
      verbose = verbose_progress)

    ecokit::cat_time(
      paste0(length(quadratic_variables), " linear and quadratic effects: "),
      level = 2L, cat_timestamp = FALSE, verbose = verbose_progress)
    ecokit::cat_time(
      paste(quadratic_variables, collapse = " + "),
      level = 3L, cat_timestamp = FALSE, verbose = verbose_progress)
  }

  form_x <- stringr::str_c(form_vars, collapse = " + ") %>%
    stringr::str_c("~ ", .) %>%
    stats::as.formula(env = baseenv())

  data_x <- dplyr::select(data_all, tidyselect::all_of(x_vars)) %>%
    as.data.frame()

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Phylogenetic tree data -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Phylogenetic tree data", verbose = verbose_progress)

  if (use_phylo_tree) {
    # Taxonomy as a proxy for phylogeny
    plant_tree <- readr::read_tsv(
      file = taxa_info_file, show_col_types = FALSE, progress = FALSE) %>%
      dplyr::mutate(
        ias_id = stringr::str_pad(ias_id, pad = "0", width = 4),
        ias_id = paste0("sp_", ias_id),
        taxon_name = NULL, species_name = NULL, species_name2 = NULL,
        species_file = NULL, species = NULL) %>%
      dplyr::filter(ias_id %in% names(data_y)) %>%
      dplyr::mutate(dplyr::across(tidyselect::everything(), factor)) %>%
      ape::as.phylo(
        ~class / order / family / genus / ias_id, data = ., collapse = FALSE)

    plant_tree$edge.length <- rep(1, length(plant_tree$edge))
  } else {
    plant_tree <- NULL
  }

  tree <- c("tree", "no_tree")[c(use_phylo_tree, no_phylo_tree)]

  ecokit::cat_time(
    paste0("Models will be fitted using ", paste(tree, collapse = " & ")),
    level = 1L, cat_timestamp = FALSE, verbose = verbose_progress)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Spatial info / random effect ------
  # # |||||||||||||||||||||||||||||||||||

  if (gpp) {

    ecokit::cat_time(
      "Spatial info and random effect", verbose = verbose_progress)
    study_design <- data.frame(sample = as.factor(seq_len(nrow(data_x))))

    data_xy <- as.matrix(dplyr::select(data_all, tidyselect::all_of(c("x", "y"))))
    rownames(data_xy) <- study_design$sample

    # Prepare gpp knots
    ecokit::cat_time(
      "Preparing GPP knots", level = 1L, verbose = verbose_progress)

    n_cores_gpp <- length(gpp_dists)

    if (n_cores_gpp > 1) {

      ecokit::set_parallel(
        n_cores = n_cores, level = 2L, future_max_size = 800L,
        strategy = strategy)
      withr::defer(future::plan("sequential", gc = TRUE))

      ecokit::cat_time(
        "Prepare GPP knots", level = 2L, verbose = verbose_progress)
      gpp_knots <- future.apply::future_lapply(
        X = gpp_dists * 1000,
        FUN = function(x) {
          IASDT.R::prepare_knots(
            coordinates = data_xy, min_distance = x, min_lf = min_lf,
            max_lf = max_lf, alphapw = alphapw)
        },
        future.scheduling = Inf, future.seed = TRUE,
        future.globals = c("data_xy", "gpp_dists", "max_lf", "min_lf", "alphapw"),
        future.packages = pkg_to_export) %>%
        stats::setNames(paste0("gpp_", gpp_dists))

      # Stopping cluster
      ecokit::set_parallel(stop_cluster = TRUE, level = 2L)
      future::plan("sequential", gc = TRUE)

    } else {

      ecokit::cat_time(
        "Working sequentially", cat_timestamp = FALSE, level = 2L,
        verbose = verbose_progress)

      gpp_knots <- purrr::map(
        .x = gpp_dists * 1000,
        .f = ~ IASDT.R::prepare_knots(
          coordinates = data_xy, min_distance = .x, min_lf = min_lf,
          max_lf = max_lf, alphapw = alphapw)) %>%
        stats::setNames(paste0("gpp_", gpp_dists))
    }

    ## Plotting knot location ----
    ecokit::cat_time(
      "Plotting GPP knots", level = 1L, verbose = verbose_progress)

    grid_r <- ecokit::load_as(path_grid_r, unwrap_r = TRUE)

    grid_r <- sf::st_as_sf(
      x = data.frame(data_xy), coords = c("x", "y"), crs = 3035) %>%
      terra::rasterize(grid_r) %>%
      terra::as.factor() %>%
      stats::setNames("grid_r")

    knots_plots <- purrr::map(
      .x = gpp_dists,
      .f = ~{

        Knot_sf <- gpp_knots[[paste0("gpp_", .x)]]$sKnot %>%
          sf::st_as_sf(coords = c("var_1", "var_2"), crs = 3035)

        n_knots <- nrow(gpp_knots[[paste0("gpp_", .x)]]$sKnot) %>%
          formatC(format = "d", big.mark = ",")

        plot <- ggplot2::ggplot(environment = emptyenv()) +
          ggplot2::geom_sf(
            data = eu_boundaries, fill = "gray95", colour = "darkgrey",
            linewidth = 0.4) +
          tidyterra::geom_spatraster(data = grid_r) +
          ggplot2::scale_fill_manual(
            values = "grey50", na.value = "transparent") +
          ggplot2::geom_sf(
            data = Knot_sf, colour = "blue", shape = 4,
            size = 2.5, stroke = 2) +
          ggplot2::scale_x_continuous(expand = c(0, 0)) +
          ggplot2::scale_y_continuous(expand = c(0, 0)) +
          ggplot2::coord_sf(
            xlim = c(2600000, 6570000), ylim = c(1360000, 5480000)) +
          ggplot2::labs(
            title = "GPP knots",
            subtitle = paste0(
              "Minimum distance between knots and between knots and grid ",
              "cells is ", .x, " km (", n_knots, " knots)")) +
          ggplot2::theme_void() +
          ggplot2::theme(
            plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
            plot.title = ggplot2::element_text(
              size = 20, face = "bold", color = "blue",
              margin = ggplot2::margin(0.1, 0, 0.125, 0.5, "cm")),
            plot.subtitle = ggplot2::element_text(
              size = 15, color = "darkgrey",
              margin = ggplot2::margin(0.125, 0, 0.125, 0.5, "cm")),
            legend.position = "none")
        return(plot)
      })

    grDevices::cairo_pdf(
      filename = fs::path(path_model, "knot_locations.pdf"),
      width = 12, height = 13.25, onefile = TRUE)
    invisible(purrr::map(knots_plots, print))
    grDevices::dev.off()

    ecokit::cat_time(
      "Saving GPP knots data", level = 1L, verbose = verbose_progress)
    save(gpp_knots, file = fs::path(path_model, "gpp_knots.RData"))

  } else {

    ecokit::cat_time(
      "Models will be fitted without spatial random effect",
      verbose = verbose_progress)
    gpp_knots <- study_design <- data_xy <- NULL
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Define the initial models -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Define the initial models", verbose = verbose_progress)

  if (gpp) {

    model_variants <- tidyr::expand_grid(
      m_thin = mcmc_thin, m_samples = mcmc_samples) %>%
      dplyr::mutate(m_transient = m_thin * mcmc_transient_factor)

    model_info <- tibble::tibble(rl = gpp_dists, rl2 = gpp_knots) %>%
      # Combinations of rl and tree
      tidyr::expand_grid(tree = tree) %>%
      dplyr::mutate(
        # Model name
        m_name_init = paste0("gpp", rl, "_", tree),
        # Save initial models
        path_m_init = purrr::map2_chr(
          .x = m_name_init, .y = rl2,
          .f = ~{

            path_out <- fs::path(
              path_model, paste0("init_mod_", .x, ".RData"))

            if (!file.exists(path_out)) {
              if (endsWith(.x, "_tree")) {
                tree <- plant_tree
              } else {
                tree <- NULL
              }

              init_model <- Hmsc::Hmsc(
                Y = data_y, XFormula = form_x, XData = data_x, distr = "probit",
                studyDesign = study_design, ranLevels = list(sample = .y),
                phyloTree = tree)

              ecokit::save_as(
                object = init_model, object_name = paste0("init_mod_", .x),
                out_path = path_out)
            }

            return(path_out)
          }),
        rl2 = NULL) %>%
      # add all combinations of thinning, number of samples and transient
      tidyr::expand_grid(model_variants) %>%
      dplyr::mutate(
        m_hpc = purrr::pmap(
          .l = list(m_name_init, m_thin, m_samples),
          .f = function(m_name_init, m_thin, m_samples) {

            m_name_fit <- paste0(m_name_init, "_samp", m_samples, "_th", m_thin)

            path_m_for_hpc <- fs::path(
              path_model, "init_mod_for_hpc",
              paste0("init_mod_", m_name_fit, ".rds"))

            list(m_name_fit = m_name_fit, path_m_for_hpc = path_m_for_hpc)

          })) %>%
      tidyr::unnest_wider("m_hpc")

  } else {

    # Non-spatial model
    model_info <- tidyr::expand_grid(
      m_thin = mcmc_thin, m_samples = mcmc_samples) %>%
      dplyr::mutate(m_transient = m_thin * mcmc_transient_factor) %>%
      tidyr::expand_grid(tree = tree) %>%
      dplyr::mutate(
        # Model name
        m_name_init = paste0("nonspatial_", tree),
        # Save initial models
        path_m_init = purrr::map_chr(
          .x = m_name_init,
          .f = ~{
            path_out <- fs::path(
              path_model, paste0("init_mod_", .x, ".RData"))

            init_model_exists <- ecokit::check_data(path_out, warning = FALSE)

            if (isFALSE(init_model_exists)) {
              if (endsWith(.x, "_tree")) {
                tree <- plant_tree
              } else {
                tree <- NULL
              }

              n_samples <- nrow(data_y)
              study_design <- data.frame(sample = as.factor(seq_len(n_samples)))
              rl <- Hmsc::HmscRandomLevel(units = study_design$sample)

              if (is.null(min_lf) && !is.null(max_lf)) {
                rl <- Hmsc::setPriors(rl, nfMax = max_lf)
              }
              if (!is.null(min_lf) && is.null(max_lf)) {
                rl <- Hmsc::setPriors(rl, nfMin = min_lf)
              }
              if (!is.null(min_lf) && !is.null(max_lf)) {
                rl <- Hmsc::setPriors(rl, nfMin = min_lf, nfMax = max_lf)
              }

              init_model <- Hmsc::Hmsc(
                Y = data_y, XFormula = form_x, XData = data_x,
                distr = "probit", phyloTree = tree,
                studyDesign = study_design, ranLevels = list(sample = rl))

              ecokit::save_as(
                object = init_model, object_name = paste0("init_mod_", .x),
                out_path = path_out)
            }
            return(path_out)
          }),

        m_hpc = purrr::pmap(
          .l = list(m_name_init, m_thin, m_samples),
          .f = function(m_name_init, m_thin, m_samples) {

            m_name_fit <- paste0(m_name_init, "_samp", m_samples, "_th", m_thin)

            path_m_for_hpc <- fs::path(
              path_model, "init_mod_for_hpc",
              paste0("init_mod_", m_name_fit, ".rds"))

            list(m_name_fit = m_name_fit, path_m_for_hpc = path_m_for_hpc)

          })) %>%
      tidyr::unnest_wider("m_hpc")
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Prepare and save unfitted models -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Save unfitted models", verbose = verbose_progress)

  if (overwrite_rds) {
    ecokit::cat_time(
      "Processing all model variants", level = 1L, verbose = verbose_progress)
  } else {

    n_mod_to_export <- sum(!file.exists(model_info$path_m_for_hpc))

    if (n_mod_to_export == 0) {
      ecokit::cat_time(
        "All model variants were already available as RDS files", level = 1L,
        verbose = verbose_progress)
    } else {
      ecokit::cat_time(
        paste0(
          n_mod_to_export, " model variants need to be exported as RDS files"),
        level = 1L, verbose = verbose_progress)
    }
  }

  # `init_fit_fun` - Function to start sampling
  init_fit_fun <- function(id) {
    curr_path <- model_info$path_m_for_hpc[id]

    if (overwrite_rds && file.exists(curr_path)) {
      fs::file_delete(curr_path)
    }

    if (isFALSE(overwrite_rds) && file.exists(curr_path) &&
        ecokit::check_data(curr_path, warning = FALSE)) {
      return(invisible(NULL))
    }

    try_n <- 0
    while (try_n < 6) {
      try_n <- try_n + 1
      model_obj <- Hmsc::sampleMcmc(
        hM = ecokit::load_as(model_info$path_m_init[id]),
        samples = model_info$m_samples[id],
        thin = model_info$m_thin[id],
        transient = model_info$m_transient[id],
        nChains = mcmc_n_chains, verbose = mcmc_verbose, engine = "HPC")

      if (to_json) {
        model_obj <- jsonify::to_json(model_obj)
      }

      saveRDS(model_obj, file = curr_path)

      if (file.exists(curr_path) &&
          ecokit::check_data(curr_path, warning = FALSE)) {
        break
      }
    }
    invisible(gc())
    invisible(NULL)
  }

  # Implement `init_fit_fun` function: start sampling and save output files
  if (n_cores > 1) {

    ecokit::set_parallel(
      n_cores = min(n_cores, nrow(model_info)), level = 1L,
      future_max_size = 800L, strategy = strategy)
    withr::defer(future::plan("sequential", gc = TRUE))

    model_process <- future.apply::future_lapply(
      X = seq_len(nrow(model_info)),
      FUN = init_fit_fun, future.scheduling = Inf, future.seed = TRUE,
      future.globals = c(
        "init_fit_fun", "model_info", "overwrite_rds",
        "mcmc_verbose", "mcmc_n_chains", "to_json"),
      future.packages = pkg_to_export)

    ecokit::set_parallel(stop_cluster = TRUE, level = 1L)
    future::plan("sequential", gc = TRUE)

  } else {
    model_process <- purrr::map(
      .x = seq_len(nrow(model_info)), .f = init_fit_fun)
  }
  rm(model_process, envir = environment())

  # Which models failed to be exported as RDS files after 5 trials
  failed_to_export <- dplyr::filter(model_info, !file.exists(path_m_for_hpc))
  if (nrow(failed_to_export) == 0) {
    ecokit::cat_time(
      "All model variants were exported as RDS files", level = 1L,
      verbose = verbose_progress)
  } else {
    ecokit::cat_time(
      paste0(
        nrow(failed_to_export),
        " model variants failed to be exported to rds files after 5 tries."),
      level = 1L, verbose = verbose_progress)
    save(
      failed_to_export, file = fs::path(path_model, "failed_to_export.RData"))
    readr::write_tsv(
      x = failed_to_export, file = fs::path(path_model, "failed_to_export.txt"))
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Prepare Hmsc-HPC fitting commands -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time(
    "Prepare Hmsc-HPC fitting commands", verbose = verbose_progress)

  model_info <- model_info %>%
    dplyr::mutate(chain = list(seq_len(mcmc_n_chains))) %>%
    tidyr::unnest_longer("chain") %>%
    dplyr::arrange(m_name_fit) %>%
    dplyr::mutate(
      m_chain = purrr::pmap(
        .l = list(m_name_fit, path_m_for_hpc, chain,
                  m_transient, m_samples, m_thin),
        .f = function(m_name_fit, path_m_for_hpc, chain,
                      m_transient, m_samples, m_thin) {

          # Turn off scientific notation
          withr::local_options(list(scipen = 999))

          # Input model
          path_model_2 <- fs::path(
            path_model, "init_mod_for_hpc", basename(path_m_for_hpc))

          # Path for posterior sampling
          path_post <- fs::path(
            path_model, "model_fitting_hpc",
            paste0(m_name_fit, "_chain", chain, "_post.rds"))

          # Path for progress
          path_progress <- fs::path(
            path_model, "model_fitting_hpc",
            paste0(m_name_fit, "_chain", chain, "_progress.txt"))

          post_missing <- isFALSE(file.exists(path_post))

          # File path for the python script
          path_model_2_for_cmd <- ecokit::normalize_path(path_model_2)
          path_post_for_cmd <- ecokit::normalize_path(path_post)
          path_progress_for_cmd <- ecokit::normalize_path(path_progress)
          path_python <- ecokit::normalize_path(path_python)

          # `TF_ENABLE_ONEDNN_OPTS=0` is used to disable the following warning:
          #
          # I tensorflow/core/util/port.cc:113] oneDNN custom operations are on.
          # You may see slightly different numerical results due to
          # floating-point round-off errors from different computation orders.
          # To turn them off, set the environment variable
          # `TF_ENABLE_ONEDNN_OPTS=0`.
          #
          # `export TF_CPP_MIN_LOG_LEVEL=3` is used to reduce debug output from
          # `Tensorflow`

          command_hpc <- paste0(
            # Not needed now as this now added to the `setup-env.sh` file
            # "export TF_CPP_MIN_LOG_LEVEL=3; export TF_ENABLE_ONEDNN_OPTS=0; ",

            "/usr/bin/time -v ", # nolint: absolute_paths_linter

            # Not needed as the python path is exported - check `setup-env.sh`
            # path_python,

            "python3 -m hmsc.run_gibbs_sampler",
            " --input ", path_model_2_for_cmd,
            " --output ", path_post_for_cmd,
            " --samples ", m_samples,
            " --transient ", m_transient,
            " --thin ", m_thin,
            " --verbose ", mcmc_verbose,
            " --chain ", (chain - 1),
            " --fp ", precision,
            " >& ", path_progress_for_cmd)

          command_ws <- paste0(
            "set TF_CPP_MIN_LOG_LEVEL=3 && set TF_ENABLE_ONEDNN_OPTS=0 && ",
            path_python,
            " -m hmsc.run_gibbs_sampler",
            " --input ", path_model_2_for_cmd,
            " --output ", path_post_for_cmd,
            " --samples ", m_samples,
            " --transient ", m_transient,
            " --thin ", m_thin,
            " --verbose ", mcmc_verbose,
            " --chain ", (chain - 1),
            " --fp ", precision,
            " > ", path_progress_for_cmd,
            " 2>&1")

          list(
            m4hpc_path_lumi = path_model_2, path_post = path_post,
            post_missing = post_missing, path_mod_progress = path_progress,
            command_hpc = command_hpc, command_ws = command_ws)

        })) %>%
    tidyr::unnest_wider("m_chain")

  # # |||||||||||||||||||||||||||||||||||
  # # Skip fitted models -----
  # # |||||||||||||||||||||||||||||||||||

  if (skip_fitted) {
    ecokit::cat_time("Skip fitted models", verbose = verbose_progress)
    models_to_fit_hpc <- dplyr::filter(model_info, post_missing) %>%
      dplyr::pull(command_hpc) %>%
      unlist()
    models_to_fit_ws <- dplyr::filter(model_info, post_missing) %>%
      dplyr::pull(command_ws) %>%
      unlist()
  } else {
    models_to_fit_hpc <- unlist(dplyr::pull(model_info, command_hpc))
    models_to_fit_ws <- unlist(dplyr::pull(model_info, command_ws))
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Save commands in a text file -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time(
    "Save model fitting commands to text file(s)", verbose = verbose_progress)

  ecokit::cat_time(
    "Save fitting commands for windows PC", level = 1L,
    verbose = verbose_progress)
  f <- file(
    description = fs::path(path_model, "commands_all_windows.txt"),
    open = "wb")
  on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
  cat(models_to_fit_ws, sep = "\n", append = FALSE, file = f)
  close(f)


  ecokit::cat_time(
    "Save fitting commands for HPC", level = 1L, verbose = verbose_progress)
  n_jobs <- length(models_to_fit_hpc)

  if (n_jobs > n_array_jobs) {
    n_splits <- ceiling((n_jobs / n_array_jobs))
    ids <- ecokit::split_vector(vector = seq_len(n_jobs), n_splits = n_splits)
  } else {
    n_splits <- 1
    ids <- list(seq_len(n_jobs))
  }

  # Save all fitting commands to single file
  ecokit::cat_time(
    "Save all fitting commands to single file",
    level = 1L, verbose = verbose_progress)
  f <- file(
    description = fs::path(path_model, "commands_all.txt"),
    open = "wb")
  on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
  cat(models_to_fit_hpc, sep = "\n", append = FALSE, file = f)
  close(f)

  # Save model fitting commands for batch SLURM jobs
  ecokit::cat_time(
    "Save model fitting commands for batch SLURM jobs", level = 1L,
    verbose = verbose_progress)
  ecokit::cat_time(
    paste0("Models will be fitted in ", n_splits, " SLURM job(s)"),
    level = 2L, cat_timestamp = FALSE, verbose = verbose_progress)

  purrr::walk(
    .x = seq_len(n_splits),
    .f = function(x) {

      if (n_splits > 1) {
        command_file <- fs::path(
          path_model, paste0("commands_to_fit_", x, ".txt"))
      } else {
        command_file <- fs::path(path_model, "commands_to_fit.txt")
      }

      # create connection to SLURM file. This is better than using sink to have
      # a platform independent file (here, to maintain a linux-like new line
      # ending)
      f <- file(command_file, open = "wb")
      on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
      cat(models_to_fit_hpc[ids[[x]]], sep = "\n", append = FALSE, file = f)
      on.exit(close(f))
    })

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Save data to disk -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Save data to disk", verbose = verbose_progress)

  set_chain_name <- function(object, chain) {
    if (is.null(object) || is.null(chain)) {
      ecokit::stop_ctx(
        "object and chain cannot be empty", object = object, chain = chain,
        include_backtrace = TRUE)
    }
    as.vector(unlist(object)) %>%
      stats::setNames(paste0("chain", unlist(chain)))
  }

  model_info <- model_info %>%
    tidyr::nest(
      path_post = path_post, path_mod_progress = path_mod_progress,
      chain = chain, command_hpc = command_hpc, command_ws = command_ws,
      post_missing = post_missing) %>%
    dplyr::mutate(
      path_post = purrr::map2(path_post, chain, set_chain_name),
      chain = purrr::map2(chain, chain, set_chain_name),
      command_hpc = purrr::map2(command_hpc, chain, set_chain_name),
      command_ws = purrr::map2(command_ws, chain, set_chain_name),
      path_mod_progress = purrr::map2(path_mod_progress, chain, set_chain_name),
      post_missing = purrr::map2(post_missing, chain, set_chain_name),
      post_aligned = NA)

  save(model_info, file = path_model_data)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Prepare SLURM file ------
  # # |||||||||||||||||||||||||||||||||||

  if (slurm_prepare) {
    ecokit::cat_time("Preparing SLURM file", verbose = verbose_progress)
    if (is.null(job_name)) {
      job_name <- stringr::str_remove_all(
        basename(path_model), paste0("_", hab_val))
    }

    IASDT.R::mod_slurm(
      model_dir = path_model, job_name = job_name,
      memory_per_cpu = memory_per_cpu, job_runtime = job_runtime,
      env_file = env_file, path_hmsc = path_hmsc, ...)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Save PA maps to disk ------

  # Save all data without cross-validation

  if (is.null(cv_fit_inherit) && (is.null(cv_fit_fold) || cv_fit_fold == 1)) {

    ecokit::cat_time("Save PA maps to disk", verbose = verbose_progress)

    # coordinates for all grid cells used in any models (for cross-validated
    # models, this includes coordinates for all training and testing sites)
    xy_all <- fs::path(path_model, "model_data.RData") %>%
      ecokit::load_as() %>%
      dplyr::select(tidyselect::all_of(c("x", "y")))

    # Mask layer for grid cells used in the models
    model_mask <- ecokit::load_as(path_grid_r, unwrap_r = TRUE) %>%
      terra::rasterize(xy_all, .)
    # Whether to use masked (exclude_cultivated) or full PA
    pa_layer <- dplyr::if_else(exclude_cultivated, "pa_masked", "pa")

    ecokit::cat_time(
      "Processing and exporting maps as tif files", level = 1L,
      verbose = verbose_progress)

    mod_pa <- sp_summary %>%
      # Select only species name and id
      dplyr::select(ias_id = ias_id, species_file) %>%
      dplyr::mutate(
        # Path storing PA maps as raster files
        file = fs::path(
          path_pa, "pa_raster",
          paste0(species_file, "_pa.RData")),
        map_sp = purrr::map2(
          .x = file, .y = ias_id,
          .f = ~ {
            # Path for storing PA map - full EU extent
            pa_file <- fs::path(path_distribution, paste0(.y, "_full.tif"))
            # Load PA map
            pa <- ecokit::load_as(.x, unwrap_r = TRUE) %>%
              magrittr::extract2(pa_layer)
            # Save PA map
            terra::writeRaster(pa, pa_file, overwrite = TRUE)

            # Path for storing PA map - only in modelling grid cells
            pa_model_file <- fs::path(path_distribution, paste0(.y, "_model.tif"))
            # mask PA map to modelling grid cells
            pa_model <- terra::mask(pa, model_mask)
            # Save masked PA map
            terra::writeRaster(pa_model, pa_model_file, overwrite = TRUE)

            tibble::tibble(
              pa = list(terra::wrap(pa)),
              pa_file = pa_file,
              pa_model = list(terra::wrap(pa_model)),
              pa_model_file = pa_model_file)
          })) %>%
      tidyr::unnest(cols = "map_sp") %>%
      dplyr::select(-tidyselect::all_of(c("species_file", "file")))

    ecokit::cat_time(
      "Calculate species richness - full", level = 1L,
      verbose = verbose_progress)
    sr_file <- fs::path(path_distribution, "sr_full.tif")
    sr <- purrr::map(mod_pa$pa, terra::unwrap) %>%
      terra::rast() %>%
      sum(na.rm = TRUE)
    terra::writeRaster(sr, sr_file, overwrite = TRUE)

    ecokit::cat_time(
      "Calculate species richness - modelling", level = 1L,
      verbose = verbose_progress)
    sr_model_file <- fs::path(path_distribution, "sr_model.tif")
    sr_model <- purrr::map(mod_pa$pa_model, terra::unwrap) %>%
      terra::rast() %>%
      sum(na.rm = TRUE)
    terra::writeRaster(sr_model, sr_model_file, overwrite = TRUE)

    # Bind sr and PA in the same tibble
    mod_pa <- tibble::tribble(
      ~ias_id, ~pa, ~pa_file, ~pa_model, ~pa_model_file,
      "sr", terra::unwrap(sr), sr_file,
      terra::unwrap(sr_model), sr_model_file) %>%
      dplyr::bind_rows(mod_pa)

    ecokit::cat_time(
      "Save maps as RData", level = 1L, verbose = verbose_progress)
    ecokit::save_as(
      object = mod_pa, object_name = "pa_with_maps",
      out_path = fs::path(path_distribution, "pa_with_maps.RData"))

    ecokit::cat_time(
      "Save maps without maps as RData", level = 1L, verbose = verbose_progress)
    mod_pa <- dplyr::select(mod_pa, -tidyselect::all_of(c("pa", "pa_model")))
    ecokit::save_as(
      object = mod_pa, object_name = "pa",
      out_path = fs::path(path_distribution, "pa.RData"))

    ecokit::cat_time(
      "Save only paths text file", level = 1L, verbose = verbose_progress)
    mod_pa %>%
      dplyr::mutate(
        pa_file = basename(pa_file),
        pa_model_file = basename(pa_model_file)) %>%
      utils::write.table(
        sep = "\t", row.names = FALSE, col.names = TRUE,
        file = fs::path(path_distribution, "pa.txt"), quote = FALSE,
        fileEncoding = "UTF-8")

    # Save maps as tar file
    ecokit::cat_time(
      "Save maps as single tar file", level = 1L, verbose = verbose_progress)
    # Path of the tar file
    tar_file <- fs::path(path_distribution, "pa_maps.tar")
    # list of files to tar
    files_to_tar <- c(
      mod_pa$pa_file, mod_pa$pa_model_file,
      fs::path(path_distribution, "pa.txt")) %>%
      basename() %>%
      paste(collapse = " ")
    tar_command <- stringr::str_glue(
      'tar -cf "{tar_file}" -C "{path_distribution}" {files_to_tar}')
    invisible(system(tar_command))

    # Change the permission of the tar file
    Sys.chmod(tar_file, "755", use_umask = FALSE)

    rm(files_to_tar, model_mask, pa_layer, envir = environment())

  }

  # # |||||||||||||||||||||||||||||||||||
  # # Save small datasets prepared in the function ------
  # # |||||||||||||||||||||||||||||||||||

  data_split <- list(
    data_all = data_all, data_y = data_y, form_x = form_x, data_x = data_x, x_vars = x_vars,
    data_cv = data_cv, use_phylo_tree = use_phylo_tree, plant_tree = plant_tree,
    tree = tree, study_design = study_design, data_xy = data_xy,
    gpp_knots = gpp_knots)

  ecokit::save_as(
    object = data_split, object_name = "model_data_subset",
    out_path = fs::path(path_model, "model_data_subset.RData"))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_diff(
    init_time = .start_time, prefix = "Processing modelling data took ",
    verbose = verbose_progress)

  return(invisible(NULL))
}
