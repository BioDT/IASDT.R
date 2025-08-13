## |------------------------------------------------------------------------| #
# mod_prepare_HPC ----
## |------------------------------------------------------------------------| #

#' Prepare initial models for model fitting with Hmsc-HPC
#'
#' The **`mod_prepare_HPC`** function prepares input data and initialises models
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
#' @param GPP Logical. Whether to fit spatial random effect using Gaussian
#'   Predictive Process. Defaults to `TRUE`. If `FALSE`, non-spatial models will
#'   be fitted.
#' @param GPP_dists Integer. Spacing (in kilometres) between GPP knots, as well
#'   as the minimum allowable distance between a knot and the nearest sampling
#'   point. The knots are generated using the [prepare_knots] function, and this
#'   value is used for both `knotDist` and `minKnotDist` in
#'   [Hmsc::constructKnots].
#' @param bio_variables Character vector. Variables from CHELSA (bioclimatic
#'   variables (bio1-bio19) and additional predictors (e.g., Net Primary
#'   Productivity, npp)) to be used in the model. By default, six ecologically
#'   relevant and minimally correlated variables are selected: `c("bio3",
#'   "bio4", "bio11", "bio18", "bio19", "npp")`.
#' @param quadratic_variables Character vector for variables for which quadratic
#'   terms are used. Defaults to all variables of the `bio_variables` in
#'   addition to soil bulk density and topographic wetness index (if used). If
#'   `quadratic_variables` is `NULL`, no quadratic terms will be used.
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
#' @param n_species_per_grid Integer. Minimum number of species required for a
#'   grid cell to be included in the analysis. This filtering occurs after
#'   applying `min_efforts_n_species` (sampling effort thresholds),
#'   `n_pres_per_species` (minimum species presence thresholds), and
#'   `exclude_0_habitat` (exclude 0% habitat coverage). Default (0): Includes
#'   all grid cells. Positive value (>0): Includes only grid cells where at
#'   least `n_species_per_grid` species are present.
#' @param CV_fit A list with three elements determining if the current model is
#'   for a specific cross-validation fold or for full dataset.
#'   - `CV_type` (character): the type of cross-validation to use. Valid
#'   options are "CV_Dist", "CV_Large", and "CV_SAC". Default: `NULL`, which
#'   means fit models on the full dataset. This can not be `NULL` if `CV_fold`
#'   is provided.
#'  - `CV_fold` (integer): the ID of the cross-validation fold to fit.
#'   For example, `CV_Fold = 4` means use the fourth fold for testing. Default:
#'   `NULL` which means no cross-validation is performed. This can not be `NULL`
#'   if `CV_type` is provided.
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
#' @param MCMC_n_chains Integer. Number of model chains. Default: 4.
#' @param MCMC_thin Integer vector. Thinning value(s) in MCMC sampling. If more
#'   than one value is provided, a separate model will be fitted at each value
#'   of thinning.
#' @param MCMC_samples Integer vector. Value(s) for the number of MCMC samples.
#'   If more than one value is provided, a separate model will be fitted at each
#'   value of number of samples. Defaults to 1000.
#' @param MCMC_transient_factor Integer. Transient multiplication factor. The
#'   value of `transient` will equal the multiplication of
#'   `MCMC_transient_factor` and `MCMC_thin`. Default: 500.
#' @param MCMC_verbose Integer. Interval at which MCMC sampling progress is
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
#' @param SLURM_prepare Logical. Whether to prepare SLURM command files. If
#'   `TRUE` (default), the SLURM commands will be saved to disk using the
#'   [mod_SLURM] function.
#' @param memory_per_cpu Character. Memory per CPU for the SLURM job. This value
#'   will be assigned to the `#SBATCH --mem-per-cpu=` SLURM argument. Example:
#'   "32G" to request 32 gigabyte. Only effective if `SLURM_prepare = TRUE`.
#'   Defaults to "64G".
#' @param job_runtime Character. Requested time for each job in the SLURM bash
#'   arrays. Example: "01:00:00" to request an hour. Only effective if
#'   `SLURM_prepare = TRUE`.
#' @param job_name Character. Name of the submitted job(s) for SLURM. If `NULL`
#'   (Default), the job name will be prepared based on the folder path and the
#'   `hab_abb` value. Only effective if `SLURM_prepare = TRUE`.
#' @param path_Hmsc  Character. Directory path to `Hmsc-HPC` extension
#'   installation. This will be provided as the `path_Hmsc` argument of the
#'   [mod_SLURM] function.
#' @param to_JSON Logical. Whether to convert unfitted models to JSON before
#'   saving to RDS file. Default: `FALSE`.
#' @param check_python Logical. Whether to check if the Python executable
#'   exists.
#' @param precision Integer. Must be either 32 or 64 (default). Defines the
#'  floating-point precision mode for `Hmsc-HPC` sampling (--fp 32 or --fp 64).
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
#' @param CV_n_grids Integer. For `CV_Dist` cross-validation strategy (see
#'   [mod_CV_prepare]), this argument determines the size of the blocks (how
#'   many grid cells in both directions).
#' @param CV_n_rows,CV_n_columns Integer. Number of rows and columns used in the
#'   `CV_Large` cross-validation strategy  (see [mod_CV_prepare]), in which the
#'   study area is divided into large blocks given the provided `CV_n_rows` and
#'   `CV_n_columns` values. Both default to 2 which means to split the study
#'   area into four large blocks at the median latitude and longitude.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param verbose_progress Logical. Whether to print a message upon successful
#'   saving of files. Defaults to `FALSE`.
#' @param ... Additional parameters provided to the [mod_SLURM] function.
#' @export
#' @inheritParams prepare_knots
#' @inheritParams mod_CV_prepare
#' @name mod_inputs
#' @rdname mod_inputs
#' @order 1
#' @importFrom rlang .data
#' @author Ahmed El-Gabbas

mod_prepare_HPC <- function(
    hab_abb = NULL, directory_name = NULL,
    min_efforts_n_species = 100L, n_pres_per_species = 80L, env_file = ".env",
    GPP = TRUE, GPP_dists = NULL, min_LF = NULL, max_LF = NULL,
    alphapw = list(Prior = NULL, Min = 20, Max = 1300, Samples = 150),
    bio_variables = c("bio3", "bio4", "bio11", "bio18", "bio19", "npp"),
    quadratic_variables = c(
      bio_variables,
      ifelse(soil_as_predictor, "soil", NULL),
      ifelse(wetness_as_predictor, "wetness", NULL)),
    efforts_as_predictor = TRUE, road_rail_as_predictor = TRUE,
    habitat_as_predictor = TRUE, river_as_predictor = FALSE,
    soil_as_predictor = TRUE, wetness_as_predictor = TRUE,
    n_species_per_grid = 0L, exclude_cultivated = TRUE,
    exclude_0_habitat = TRUE, CV_n_folds = 4L, CV_n_grids = 20L,
    CV_n_rows = 2L, CV_n_columns = 2L, CV_plot = TRUE,
    CV_SAC = FALSE,
    CV_fit = list(CV_type = NULL, CV_fold = NULL, inherit_dir = NULL),
    use_phylo_tree = TRUE, no_phylo_tree = FALSE, overwrite_rds = TRUE,
    n_cores = 8L, strategy = "multisession", MCMC_n_chains = 4L,
    MCMC_thin = NULL, MCMC_samples = 1000L, MCMC_transient_factor = 500L,
    MCMC_verbose = 200L, skip_fitted = TRUE, n_array_jobs = 210L,
    model_country = NULL, verbose_progress = TRUE, SLURM_prepare = TRUE,
    memory_per_cpu = "64G", job_runtime = NULL, job_name = NULL,
    path_Hmsc = NULL, check_python = FALSE, to_JSON = FALSE, precision = 64L,
    ...) {

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Initial checking -----
  # # |||||||||||||||||||||||||||||||||||

  .start_time <- lubridate::now(tzone = "CET")


  CheckNULL <- c(
    "directory_name", "n_pres_per_species", "MCMC_thin", "MCMC_samples",
    "memory_per_cpu", "path_Hmsc", "hab_abb")
  IsNull <- purrr::map_lgl(CheckNULL, ~ is.null(get(.x)))

  if (any(IsNull)) {
    ecokit::stop_ctx(
      paste0(
        paste0("`", CheckNULL[which(IsNull)], "`", collapse = ", "),
        " can not be empty"),
      sum_IsNull = sum(IsNull), directory_name = directory_name,
      n_pres_per_species = n_pres_per_species, MCMC_thin = MCMC_thin,
      MCMC_samples = MCMC_samples, memory_per_cpu = memory_per_cpu,
      path_Hmsc = path_Hmsc, hab_abb = hab_abb, include_backtrace = TRUE)
  }

  if (!(precision %in% c(32, 64))) {
    ecokit::stop_ctx(
      "`precision` should be either of 32 or 64", precision = precision,
      include_backtrace = TRUE)
  }

  hab_abb <- .validate_hab_abb(as.character(hab_abb))
  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)

  if (!all(is.numeric(MCMC_samples)) || any(MCMC_samples <= 0)) {
    ecokit::stop_ctx(
      "`MCMC_samples` should be numeric and greater than zero",
      MCMC_samples = MCMC_samples, include_backtrace = TRUE)
  }

  if (!all(is.numeric(MCMC_thin)) || any(MCMC_thin <= 0)) {
    ecokit::stop_ctx(
      "`MCMC_thin` should be numeric and greater than zero",
      MCMC_thin = MCMC_thin, include_backtrace = TRUE)
  }

  if (!all(is.numeric(n_pres_per_species)) || n_pres_per_species <= 0) {
    ecokit::stop_ctx(
      "`n_pres_per_species` should be numeric and greater than zero",
      n_pres_per_species = n_pres_per_species, include_backtrace = TRUE)
  }

  if (!all(is.numeric(min_efforts_n_species)) || min_efforts_n_species <= 0) {
    ecokit::stop_ctx(
      "`min_efforts_n_species` should be numeric and greater than zero",
      min_efforts_n_species = min_efforts_n_species, include_backtrace = TRUE)
  }

  if (!is.numeric(n_species_per_grid) || n_species_per_grid < 0) {
    ecokit::stop_ctx(
      "`n_species_per_grid` has to be integer >= 0",
      n_species_per_grid = n_species_per_grid, include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  NCells <- Sp <- IAS_ID <- x <- y <- Country <- M_thin <- rL <-
    M_Name_init <- rL2 <- M_samples <- M4HPC_Path <- M_transient <-
    M_Name_Fit <- Chain <- Post_Missing <- Command_HPC <- Command_WS <-
    Post_Path <- Path_ModProg <- TaxaInfoFile <- Path_Grid <- EU_Bound <-
    Path_PA <- NAME_ENGL <- NSp <- Species_File <- File <- ias_id <-
    PA <- PA_model <- PA_file <- PA_model_file <- path_model <- NULL

  Path_Python <- fs::path(path_Hmsc, "Scripts", "python.exe")

  ecokit::info_chunk(
    "Preparing data for Hmsc-HPC models", line_char = "=",
    verbose = verbose_progress)

  # # |||||||||||||||||||||||||||||||||||
  # Load/check environment variables -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time(
    "Load and check environment variables", verbose = verbose_progress)

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "TaxaInfoFile", "DP_R_Taxa_info", FALSE, TRUE,
    "EU_Bound", "DP_R_EUBound", FALSE, TRUE,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE,
    "path_model", "DP_R_Model_path", FALSE, FALSE,
    "Path_PA", "DP_R_PA", TRUE, FALSE)

  # Check if Python executable exists
  if (check_python && !file.exists(Path_Python) && Sys.info()[1] == "Windows") {
    ecokit::stop_ctx(
      "Python executable does not exist", Path_Python = Path_Python,
      include_backtrace = TRUE)
  }

  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())


  # Validate structure of CV_fit list argument
  if (!inherits(CV_fit, "list") || length(CV_fit) != 3 ||
      !setequal(names(CV_fit), c("CV_type", "CV_fold", "inherit_dir"))) {
    ecokit::stop_ctx(
      paste0(
        "`CV_fit` must be a list with exactly three elements: ",
        "`CV_type`, `CV_fold`, and `inherit_dir` (no more, no less)"),
      CV_fit = CV_fit, include_backtrace = TRUE)
  }

  CV_fit_type <- CV_fit$CV_type
  CV_fit_fold <- CV_fit$CV_fold
  CV_fit_inherit <- CV_fit$inherit_dir

  if (!is.null(CV_fit_inherit)) {
    # Validate inherited cross-validation directory name
    if (length(CV_fit_inherit) != 1 || !is.character(CV_fit_inherit) ||
        !nzchar(CV_fit_inherit)) {
      ecokit::stop_ctx(
        "Directory to inherit data from must be a non-empty character string",
        CV_fit = CV_fit, include_backtrace = TRUE)
    }

    # Validate inherited cross-validation directory
    CV_fit_inherit <- fs::path(path_model, CV_fit_inherit)
    if (!fs::dir_exists(CV_fit_inherit)) {
      ecokit::stop_ctx(
        "Directory to inherit data from does not exist",
        CV_fit = CV_fit, include_backtrace = TRUE)
    }

    # Validate files existence in the inherited directory
    files_to_check <- fs::path(
      CV_fit_inherit, c("ModDT.RData", "CV_data.RData"))
    purrr::walk(
      files_to_check,
      ~ {
        if (!ecokit::check_data(.x, warning = FALSE)) {
          list_rdata <- list.files(CV_fit_inherit, pattern = "RData$")
          ecokit::stop_ctx(
            paste0("File does not exist in the inherited directory: ", .x),
            CV_fit = CV_fit, rdata_files = list_rdata, include_backtrace = TRUE)
        }
      })
  }

  # CV_fit_type and CV_fit_fold must be consistent
  if (
    (is.null(CV_fit_type) && !is.null(CV_fit_fold)) ||
    (!is.null(CV_fit_type) && is.null(CV_fit_fold))) {
    missing_val <- dplyr::if_else(
      is.null(CV_fit_type), "`CV_fit$CV_type`", "`CV_fit$CV_fold`")
    present_val <- dplyr::if_else(
      !is.null(CV_fit_type), "`CV_fit$CV_type`", "`CV_fit$CV_fold`")
    ecokit::stop_ctx(
      paste0(
        "Both `CV_fit$CV_type` and `CV_fit$CV_fold` must be either NULL ",
        "(both) or both non-NULL. Currently, ", missing_val,
        " is missing and ", present_val, " is provided."),
      CV_fit = CV_fit, include_backtrace = TRUE)
  }

  # Validate CV_fit_type and CV_fit_fold
  if (!is.null(CV_fit_type) && !is.null(CV_fit_fold)) {
    valid_cv_types <- c("CV_Dist", "CV_Large", "CV_SAC")
    if (length(CV_fit_type) != 1L || !is.character(CV_fit_type) ||
        !nzchar(CV_fit_type)) {
      ecokit::stop_ctx(
        "`CV_fit_type` must be a single character string",
        CV_fit = CV_fit, include_backtrace = TRUE)
    }

    if (!CV_fit_type %in% valid_cv_types) {
      ecokit::stop_ctx(
        "Invalid `CV_fit_type` value",
        CV_fit = CV_fit, valid_cv_types = valid_cv_types,
        include_backtrace = TRUE)
    }

    if (is.null(CV_fit_fold) || !is.numeric(CV_fit_fold) ||
        length(CV_fit_fold) != 1L || CV_fit_fold < 1L) {
      ecokit::stop_ctx(
        "`CV_fit_fold` must be a single positive integer",
        CV_fit = CV_fit, include_backtrace = TRUE)
    }

    # Flag for filtering data based on this cross-validation
    filter_CV_data <- TRUE
  } else {
    filter_CV_data <- FALSE
  }

  path_model <- fs::path(path_model, directory_name)
  if (fs::dir_exists(path_model) && length(fs::dir_ls(path_model)) > 0) {
    ecokit::stop_ctx(
      "Model directory already exists", path_model = path_model,
      include_backtrace = TRUE)
  }
  fs::dir_create(path_model)

  Path_GridR <- fs::path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(Path_GridR)) {
    ecokit::stop_ctx(
      "Path for the Europe boundaries does not exist", Path_GridR = Path_GridR,
      include_backtrace = TRUE)
  }

  ecokit::record_arguments(
    out_path = fs::path(path_model, "Args_mod_prepare_HPC.RData"))

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c(
      "dplyr", "sf", "Hmsc", "jsonify", "magrittr", "IASDT.R", "ecokit"),
    strategy = strategy)

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Check input arguments ----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Checking input arguments", verbose = verbose_progress)

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  CharArgs <- c(
    "hab_abb", "directory_name", "path_Hmsc", "env_file", "bio_variables")
  ecokit::check_args(
    args_all = AllArgs, args_to_check = CharArgs, args_type = "character")

  LogicArgs <- c(
    "use_phylo_tree", "no_phylo_tree",
    "exclude_0_habitat", "skip_fitted", "verbose_progress", "to_JSON",
    "CV_SAC", "check_python", "overwrite_rds", "SLURM_prepare",
    "exclude_cultivated", "GPP", "efforts_as_predictor",
    "road_rail_as_predictor", "habitat_as_predictor", "river_as_predictor",
    "soil_as_predictor", "wetness_as_predictor")
  ecokit::check_args(
    args_all = AllArgs, args_to_check = LogicArgs, args_type = "logical")

  NumericArgs <- c(
    "n_cores", "MCMC_n_chains", "MCMC_thin", "MCMC_samples", "MCMC_verbose",
    "n_array_jobs", "n_pres_per_species", "min_efforts_n_species",
    "MCMC_transient_factor", "CV_n_folds", "CV_n_grids", "CV_n_rows",
    "CV_n_columns", "precision")
  if (GPP) {
    NumericArgs <- c(NumericArgs, "GPP_dists")
  }
  ecokit::check_args(
    args_all = AllArgs, args_to_check = NumericArgs, args_type = "numeric")


  if (SLURM_prepare) {
    ecokit::check_args(
      args_all = AllArgs, args_type = "character",
      args_to_check = c("memory_per_cpu", "job_runtime"))
  }

  # Phylogenetic tree options
  if (isFALSE(use_phylo_tree) && isFALSE(no_phylo_tree)) {
    ecokit::stop_ctx(
      "At least one of `use_phylo_tree` or `no_phylo_tree` has to be true",
      use_phylo_tree = use_phylo_tree, no_phylo_tree = no_phylo_tree,
      include_backtrace = TRUE)
  }

  NumArgsInvalid <- purrr::map_lgl(.x = NumericArgs, .f = ~all(get(.x) < 1))
  if (any(NumArgsInvalid)) {
    ecokit::stop_ctx(
      "Some parameter(s) can not be < 1",
      NumArgsInvalid = NumericArgs[NumArgsInvalid], include_backtrace = TRUE)
  }

  rm(AllArgs, CharArgs, LogicArgs, NumericArgs, envir = environment())

  # # ..................................................................... ###

  if (GPP) {
    if (is.null(GPP_dists)) {
      ecokit::stop_ctx(
        "`GPP_dists` can not be empty", GPP_dists = GPP_dists,
        include_backtrace = TRUE)
    }
    if (!all(is.numeric(GPP_dists)) || any(GPP_dists <= 0)) {
      ecokit::stop_ctx(
        "`GPP_dists` should be numeric and greater than zero",
        GPP_dists = GPP_dists, include_backtrace = TRUE)
    }
  }

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # File paths - Creating missing paths ----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time(
    "File paths - Creating missing paths", verbose = verbose_progress)
  fs::dir_create(fs::path(path_model, "InitMod4HPC"))
  fs::dir_create(fs::path(path_model, "Model_Fitting_HPC"))
  # Also create directory for SLURM outputs
  fs::dir_create(fs::path(path_model, "Model_Fitting_HPC", "JobsLog"))
  # Directory to save species distribution tiffs
  Path_Dist <- fs::path(path_model, "Species_distribution")
  fs::dir_create(Path_Dist)

  Path_ModelDT <- fs::path(path_model, "Model_Info.RData")

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # # Prepare list of predictors -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Prepare list of predictors", verbose = verbose_progress)

  # Check bio_variables values
  if (isFALSE(all(bio_variables %in% IASDT.R::CHELSA_variables$Variable))) {
    WrongBio <- bio_variables[
      which(!(bio_variables %in% IASDT.R::CHELSA_variables$Variable))]
    ecokit::stop_ctx(
      "Invalid Bioclimatic variables", WrongBio = WrongBio,
      include_backtrace = TRUE)
  }

  XVars <- bio_variables

  if (efforts_as_predictor) XVars <- c(XVars, "EffortsLog")
  if (road_rail_as_predictor) XVars <- c(XVars, "RoadRailLog")
  if (river_as_predictor) XVars <- c(XVars, "RiversLog")
  if (soil_as_predictor) XVars <- c(XVars, "soil")
  if (wetness_as_predictor) XVars <- c(XVars, "wetness")

  if (hab_abb != "0" && habitat_as_predictor) XVars <- c(XVars, "HabLog")

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # # Preparing input data -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Preparing input data", verbose = verbose_progress)

  HabVal <- c(
    "0_All", "1_Forests", "2_Open_forests", "3_Scrub",
    "4a_Natural_grasslands", "4b_Human_maintained_grasslands",
    "10_Wetland", "12a_Ruderal_habitats", "12b_Agricultural_habitats") %>%
    stringr::str_subset(paste0("^", as.character(hab_abb), "_"))


  if (is.null(CV_fit_inherit)) {

    ecokit::info_chunk(
      "Preparing input data using IASDT.R::mod_prepare_data",
      verbose = verbose_progress, line_char_rep = 65)

    DT_All <- IASDT.R::mod_prepare_data(
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
    DT_All <- ecokit::load_as(fs::path(CV_fit_inherit, "ModDT.RData"))

    # copy ModDT.RData
    fs::file_copy(
      fs::path(CV_fit_inherit, "ModDT.RData"),
      fs::path(path_model, "ModDT.RData"), overwrite = TRUE)

  }

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # # Subsetting study area -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Subsetting study area", verbose = verbose_progress)

  if (!is.null(model_country)) {

    ValidCountries <- model_country %in% unique(DT_All$Country)

    if (!all(ValidCountries)) {
      ecokit::stop_ctx(
        "Invalid country names",
        invalid_countries = model_country[!ValidCountries],
        include_backtrace = TRUE)
    }

    ecokit::cat_time(
      paste0(
        "Subsetting data to: ", paste(sort(model_country), collapse = " & ")),
      level = 1L, verbose = verbose_progress)

    Sample_ExclSp <- dplyr::filter(DT_All, Country %in% model_country) %>%
      dplyr::summarise(
        dplyr::across(tidyselect::starts_with("Sp_"), sum)) %>%
      tidyr::pivot_longer(
        cols = tidyselect::everything(),
        names_to = "Sp", values_to = "NCells") %>%
      dplyr::filter(NCells < n_pres_per_species) %>%
      dplyr::pull(Sp)

    ecokit::cat_time(
      paste0(length(Sample_ExclSp), " species are excluded"), level = 1L,
      verbose = verbose_progress)
    DT_All <- dplyr::filter(DT_All, Country %in% model_country) %>%
      dplyr::select(-tidyselect::all_of(Sample_ExclSp))

    # # |||||||||||||||||||||||||||||||||||
    ## Plotting subsetted data -----
    # # |||||||||||||||||||||||||||||||||||

    ecokit::cat_time(
      "Plotting subsetted data", level = 1L, verbose = verbose_progress)

    NSpSubset <- DT_All %>%
      dplyr::mutate(
        NSp = rowSums(
          dplyr::select(., tidyselect::starts_with("Sp_")),
          na.rm = TRUE),
        NSp = as.integer(NSp)) %>%
      dplyr::select("x", "y", "NSp") %>%
      sf::st_as_sf(coords = c("x", "y"), crs = 3035) %>%
      terra::rasterize(
        y = ecokit::load_as(Path_GridR, unwrap_r = TRUE), field = "NSp") %>%
      terra::trim()

    if (hab_abb == "0") {
      Hab_column <- NULL
    } else {
      Hab_column <- c(
        "Hab_1_Forests", "Hab_2_Open_forests", "Hab_3_Scrub",
        "Hab_4a_Natural_grasslands", "Hab_4b_Human_maintained_grasslands",
        "Hab_10_Wetland", "Hab_12a_Ruderal_habitats",
        "Hab_12b_Agricultural_habitats") %>%
        stringr::str_subset(paste0("_", as.character(hab_abb), "_")) %>%
        stringr::str_remove("Hab_")
    }

    EU_Bound_sub <- ecokit::load_as(EU_Bound) %>%
      magrittr::extract2("Bound_sf_Eur") %>%
      magrittr::extract2("L_03") %>%
      dplyr::filter(NAME_ENGL %in% model_country)

    Limits <- as.vector(terra::ext(NSpSubset))
    # Relative JPEG height
    DimX <- Limits[2] - Limits[1]
    DimY <- Limits[4] - Limits[3]
    plot_height <- (DimY * 25) / (DimX * 0.95)
    NGrids <- format(nrow(DT_All), big.mark = ",")
    NSp <- length(stringr::str_subset(names(DT_All), "Sp_"))

    Caption <- paste(sort(model_country), collapse = "; ") %>%
      stringr::str_wrap(width = 110) %>%
      paste0("<strong>Selected countries</strong>: <br/>", .) %>%
      stringr::str_replace_all("\n", "<br/>")
    Subtitle <- paste0(
      NSp, " IAS within \u2265", n_pres_per_species,
      " presence grid cells in the selected countries (",
      NGrids, " grid cells)")

    NSpPerGrid_Sub <- ggplot2::ggplot() +
      tidyterra::geom_spatraster(data = NSpSubset) +
      tidyterra::scale_fill_whitebox_c(
        na.value = "transparent", palette = "bl_yl_rd", name = NULL) +
      ggplot2::labs(
        title = paste0(
          '<span style="color:blue; font-size:24px;"><b>',
          "Number of IAS per grid cell to be used in the models</b></span>",
          '<span style="color:black; font-size:18px;"> (',
          Hab_column, ")</span>"),
        subtitle = Subtitle, caption = Caption) +
      ggplot2::geom_sf(
        data = EU_Bound_sub, fill = "transparent", colour = "black") +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = Limits[c(3, 4)]) +
      ggplot2::scale_x_continuous(expand = c(0, 0), limits = Limits[c(1, 2)]) +
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
      filename = fs::path(path_model, "NSpPerGrid_Sub.jpeg"),
      width = 25, height = plot_height, res = 600, quality = 100, units = "cm")

    print(NSpPerGrid_Sub)
    grDevices::dev.off()

    rm(Limits, NSpPerGrid_Sub, EU_Bound_sub, envir = environment())

  } else {
    ecokit::cat_time(
      "No data subsetting was implemented", level = 1L,
      verbose = verbose_progress)
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
    EmptyGridsID <- dplyr::select(DT_All, tidyselect::starts_with("Sp_")) %>%
      rowSums() %>%
      magrittr::is_less_than(as.integer(n_species_per_grid)) %>%
      which() %>%
      magrittr::multiply_by(-1)

    if (length(EmptyGridsID) > 0) {
      ecokit::cat_time(
        paste0(
          "Excluding grid cells with < ", n_species_per_grid,
          " species presence"),
        level = 1L, verbose = verbose_progress)
      DT_All <- dplyr::slice(DT_All, EmptyGridsID)
    }
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Cross-validation ----
  # # |||||||||||||||||||||||||||||||||||

  if (is.null(CV_fit_inherit)) {

    ecokit::cat_time(
      "Prepare cross-validation folds", verbose = verbose_progress)

    DT_All <- IASDT.R::mod_CV_prepare(
      input_data = DT_All, env_file = env_file, x_vars = XVars,
      CV_n_folds = CV_n_folds, CV_n_grids = CV_n_grids, CV_n_rows = CV_n_rows,
      CV_n_columns = CV_n_columns, out_path = path_model, CV_plot = CV_plot,
      CV_SAC = CV_SAC)

  } else {

    ecokit::cat_time(
      "Loading cross-validation data", verbose = verbose_progress)
    CV_data <- ecokit::load_as(fs::path(CV_fit_inherit, "CV_data.RData"))

    # copy CV_Blocks.pdf and CV_data.RData files
    fs::file_copy(
      fs::path(CV_fit_inherit, "CV_Blocks.pdf"),
      fs::path(path_model, "CV_Blocks.pdf"), overwrite = TRUE)
    fs::file_copy(
      fs::path(CV_fit_inherit, "CV_data.RData"),
      fs::path(path_model, "CV_data.RData"), overwrite = TRUE)

    # Add cross-validation information
    if (!is.null(CV_data$CV_SAC)) {
      DT_All$CV_SAC <- magrittr::extract2(CV_data$CV_SAC, "folds_ids")
    }

    DT_All$CV_Dist <- magrittr::extract2(CV_data$CV_Dist, "folds_ids")
    DT_All$CV_Large <- magrittr::extract2(CV_data$CV_Large, "folds_ids")

  }

  # Filter data if CV info is provided
  if (filter_CV_data) {
    ecokit::cat_time(
      paste0(
        "Filtering data for testing cross-validation fold: ",
        CV_fit_type, " -- ", CV_fit_fold),
      verbose = verbose_progress)

    DT_All <- dplyr::filter(DT_All, !!rlang::sym(CV_fit_type) != CV_fit_fold)

  } else {

    ecokit::cat_time(
      "No cross-validation filtering was made", verbose = verbose_progress)

  }

  # cross-validation data to be saved
  DT_CV <- DT_All %>%
    dplyr::select(
      "CellNum", "CellCode", "Country", tidyselect::starts_with("CV"))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Response - Y matrix ----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Response - Y matrix", verbose = verbose_progress)
  DT_y <- dplyr::select(DT_All, tidyselect::starts_with("Sp_")) %>%
    as.data.frame()
  ecokit::cat_time(
    paste0(ncol(DT_y), " species"), level = 1L, cat_timestamp = FALSE,
    verbose = verbose_progress)

  # Prepare species info summary only if no inherited directory is used

  if (is.null(CV_fit_inherit)) {
    ecokit::cat_time(
      "Save species summary", level = 1L, cat_timestamp = FALSE,
      verbose = verbose_progress)
    SpSummary <- fs::path(Path_PA, "Sp_PA_Summary_DF.RData")
    if (!file.exists(SpSummary)) {
      ecokit::stop_ctx(
        "SpSummary file does not exist", SpSummary = SpSummary,
        include_backtrace = TRUE)
    }
    SpSummary <- ecokit::load_as(SpSummary) %>%
      dplyr::arrange(IAS_ID) %>%
      dplyr::mutate(
        IAS_ID = stringr::str_pad(IAS_ID, pad = "0", width = 4),
        IAS_ID = paste0("Sp_", IAS_ID)) %>%
      dplyr::filter(IAS_ID %in% names(DT_y))

    save(SpSummary, file = fs::path(path_model, "SpSummary.RData"))

  } else {

    fs::file_copy(
      fs::path(CV_fit_inherit, "SpSummary.RData"),
      fs::path(path_model, "SpSummary.RData"), overwrite = TRUE)

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
    FormVars <- XVars
    ecokit::cat_time(
      paste0(
        "Models will be fitted using ", length(XVars), " predictors: ",
        paste(XVars, collapse = " + ")), level = 1L,
      cat_timestamp = FALSE, verbose = verbose_progress)
  } else {
    OnlyLinear <- setdiff(XVars, quadratic_variables)
    FormVars <- c(
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

  Form_x <- stringr::str_c(FormVars, collapse = " + ") %>%
    stringr::str_c("~ ", .) %>%
    stats::as.formula(env = baseenv())

  DT_x <- dplyr::select(DT_All, tidyselect::all_of(XVars)) %>%
    as.data.frame()

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Phylogenetic tree data -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Phylogenetic tree data", verbose = verbose_progress)

  if (use_phylo_tree) {
    # Taxonomy as a proxy for phylogeny
    plant.tree <- readr::read_tsv(
      file = TaxaInfoFile, show_col_types = FALSE, progress = FALSE) %>%
      dplyr::mutate(
        IAS_ID = stringr::str_pad(IAS_ID, pad = "0", width = 4),
        IAS_ID = paste0("Sp_", IAS_ID),
        taxon_name = NULL, Species_name = NULL, Species_name2 = NULL,
        Species_File = NULL, Species = NULL) %>%
      dplyr::filter(IAS_ID %in% names(DT_y)) %>%
      dplyr::mutate(dplyr::across(tidyselect::everything(), factor)) %>%
      ape::as.phylo(
        ~Class / Order / Family / Genus / IAS_ID, data = ., collapse = FALSE)

    plant.tree$edge.length <- rep(1, length(plant.tree$edge))
  } else {
    plant.tree <- NULL
  }

  Tree <- c("Tree", "NoTree")[c(use_phylo_tree, no_phylo_tree)]

  ecokit::cat_time(
    paste0("Models will be fitted using ", paste(Tree, collapse = " & ")),
    level = 1L, cat_timestamp = FALSE, verbose = verbose_progress)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Spatial info / random effect ------
  # # |||||||||||||||||||||||||||||||||||

  if (GPP) {

    ecokit::cat_time(
      "Spatial info and random effect", verbose = verbose_progress)
    studyDesign <- data.frame(sample = as.factor(seq_len(nrow(DT_x))))

    DT_xy <- as.matrix(dplyr::select(DT_All, x, y))
    rownames(DT_xy) <- studyDesign$sample

    # Prepare GPP knots
    ecokit::cat_time(
      "Preparing GPP knots", level = 1L, verbose = verbose_progress)

    NCores_GPP <- length(GPP_dists)

    if (NCores_GPP > 1) {

      ecokit::set_parallel(
        n_cores = n_cores, level = 2L, future_max_size = 800L,
        strategy = strategy)
      withr::defer(future::plan("sequential", gc = TRUE))

      ecokit::cat_time(
        "Prepare GPP knots", level = 2L, verbose = verbose_progress)
      GPP_Knots <- future.apply::future_lapply(
        X = GPP_dists * 1000,
        FUN = function(x) {
          IASDT.R::prepare_knots(
            coordinates = DT_xy, min_distance = x, min_LF = min_LF,
            max_LF = max_LF, alphapw = alphapw)
        },
        future.scheduling = Inf, future.seed = TRUE,
        future.globals = c("DT_xy", "GPP_dists", "max_LF", "min_LF", "alphapw"),
        future.packages = pkg_to_export) %>%
        stats::setNames(paste0("GPP_", GPP_dists))

      # Stopping cluster
      ecokit::set_parallel(stop_cluster = TRUE, level = 2L)
      future::plan("sequential", gc = TRUE)

    } else {

      ecokit::cat_time(
        "Working sequentially", cat_timestamp = FALSE, level = 2L,
        verbose = verbose_progress)

      GPP_Knots <- purrr::map(
        .x = GPP_dists * 1000,
        .f = ~ IASDT.R::prepare_knots(
          coordinates = DT_xy, min_distance = .x, min_LF = min_LF,
          max_LF = max_LF, alphapw = alphapw)) %>%
        stats::setNames(paste0("GPP_", GPP_dists))
    }

    ## Plotting knot location ----
    ecokit::cat_time(
      "Plotting GPP knots", level = 1L, verbose = verbose_progress)

    GridR <- ecokit::load_as(Path_GridR, unwrap_r = TRUE)

    GridR <- sf::st_as_sf(
      x = data.frame(DT_xy), coords = c("x", "y"), crs = 3035) %>%
      terra::rasterize(GridR) %>%
      terra::as.factor() %>%
      stats::setNames("GridR")

    EU_Bound <- ecokit::load_as(EU_Bound) %>%
      magrittr::extract2("Bound_sf_Eur_s") %>%
      magrittr::extract2("L_03") %>%
      suppressWarnings()

    Knots_Plots <- purrr::map(
      .x = GPP_dists,
      .f = ~{

        Knot_sf <- GPP_Knots[[paste0("GPP_", .x)]]$sKnot %>%
          sf::st_as_sf(coords = c("Var1", "Var2"), crs = 3035)

        NKnots <- nrow(GPP_Knots[[paste0("GPP_", .x)]]$sKnot) %>%
          formatC(format = "d", big.mark = ",")

        Plot <- ggplot2::ggplot() +
          ggplot2::geom_sf(
            data = EU_Bound, fill = "gray95", colour = "darkgrey",
            linewidth = 0.4) +
          tidyterra::geom_spatraster(data = GridR) +
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
              "cells is ", .x, " km (", NKnots, " knots)")) +
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

        return(Plot)
      })

    grDevices::cairo_pdf(
      filename = fs::path(path_model, "knot_Locations.pdf"),
      width = 12, height = 13.25, onefile = TRUE)
    invisible(purrr::map(Knots_Plots, print))
    grDevices::dev.off()

    ecokit::cat_time(
      "Saving GPP knots data", level = 1L, verbose = verbose_progress)
    save(GPP_Knots, file = fs::path(path_model, "GPP_Knots.RData"))

  } else {

    ecokit::cat_time(
      "Models will be fitted without spatial random effect",
      verbose = verbose_progress)
    GPP_Knots <- studyDesign <- DT_xy <- NULL
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Define the initial models -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Define the initial models", verbose = verbose_progress)

  if (GPP) {

    ModelVariants <- tidyr::expand_grid(
      M_thin = MCMC_thin, M_samples = MCMC_samples) %>%
      dplyr::mutate(M_transient = M_thin * MCMC_transient_factor)

    Model_Info <- tibble::tibble(rL = GPP_dists, rL2 = GPP_Knots) %>%
      # Combinations of rL and Tree
      tidyr::expand_grid(Tree = Tree) %>%
      dplyr::mutate(
        # Model name
        M_Name_init = paste0("GPP", rL, "_", Tree),
        # Save initial models
        M_Init_Path = purrr::map2_chr(
          .x = M_Name_init, .y = rL2,
          .f = ~{

            PathOut <- fs::path(
              path_model, paste0("InitMod_", .x, ".RData"))

            if (!file.exists(PathOut)) {
              if (endsWith(.x, "_Tree")) {
                Tree <- plant.tree
              } else {
                Tree <- NULL
              }

              InitModel <- Hmsc::Hmsc(
                Y = DT_y, XFormula = Form_x, XData = DT_x, distr = "probit",
                studyDesign = studyDesign, ranLevels = list(sample = .y),
                phyloTree = Tree)

              ecokit::save_as(
                object = InitModel, object_name = paste0("InitMod_", .x),
                out_path = PathOut)
            }

            return(PathOut)
          }),
        rL2 = NULL) %>%
      # add all combinations of thinning, number of samples and transient
      tidyr::expand_grid(ModelVariants) %>%
      dplyr::mutate(
        M_HPC = purrr::pmap(
          .l = list(M_Name_init, M_thin, M_samples),
          .f = function(M_Name_init, M_thin, M_samples) {

            M_Name_Fit <- paste0(M_Name_init, "_samp", M_samples, "_th", M_thin)

            M4HPC_Path <- fs::path(
              path_model, "InitMod4HPC", paste0("InitMod_", M_Name_Fit, ".rds"))

            return(list(M_Name_Fit = M_Name_Fit, M4HPC_Path = M4HPC_Path))

          })) %>%
      tidyr::unnest_wider("M_HPC")

  } else {

    # Non-spatial model
    Model_Info <- tidyr::expand_grid(
      M_thin = MCMC_thin, M_samples = MCMC_samples) %>%
      dplyr::mutate(M_transient = M_thin * MCMC_transient_factor) %>%
      tidyr::expand_grid(Tree = Tree) %>%
      dplyr::mutate(
        # Model name
        M_Name_init = paste0("NonSpatial_", Tree),
        # Save initial models
        M_Init_Path = purrr::map_chr(
          .x = M_Name_init,
          .f = ~{
            PathOut <- fs::path(
              path_model, paste0("InitMod_", .x, ".RData"))

            InitModExists <- ecokit::check_data(PathOut)

            if (isFALSE(InitModExists)) {
              if (endsWith(.x, "_Tree")) {
                Tree <- plant.tree
              } else {
                Tree <- NULL
              }

              InitModel <- Hmsc::Hmsc(
                Y = DT_y, XFormula = Form_x, XData = DT_x,
                distr = "probit", phyloTree = Tree)

              ecokit::save_as(
                object = InitModel, object_name = paste0("InitMod_", .x),
                out_path = PathOut)
            }
            return(PathOut)
          }),

        M_HPC = purrr::pmap(
          .l = list(M_Name_init, M_thin, M_samples),
          .f = function(M_Name_init, M_thin, M_samples) {

            M_Name_Fit <- paste0(M_Name_init, "_samp", M_samples, "_th", M_thin)

            M4HPC_Path <- fs::path(
              path_model, "InitMod4HPC", paste0("InitMod_", M_Name_Fit, ".rds"))

            return(list(M_Name_Fit = M_Name_Fit, M4HPC_Path = M4HPC_Path))

          })) %>%
      tidyr::unnest_wider("M_HPC")
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Prepare and save unfitted models -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Save unfitted models", verbose = verbose_progress)

  if (overwrite_rds) {
    ecokit::cat_time(
      "Processing all model variants", level = 1L, verbose = verbose_progress)
  } else {

    NMod2Export <- sum(!file.exists(Model_Info$M4HPC_Path))

    if (NMod2Export == 0) {
      ecokit::cat_time(
        "All model variants were already available as RDS files", level = 1L,
        verbose = verbose_progress)
    } else {
      ecokit::cat_time(
        paste0(
          NMod2Export, " model variants need to be exported as RDS files"),
        level = 1L, verbose = verbose_progress)
    }
  }

  # `InitFitFun` - Function to start sampling
  InitFitFun <- function(ID) {
    CurrPath <- Model_Info$M4HPC_Path[ID]

    if (overwrite_rds && file.exists(CurrPath)) {
      fs::file_delete(CurrPath)
    }

    if (isFALSE(overwrite_rds) && file.exists(CurrPath) &&
        ecokit::check_data(CurrPath)) {
      return(invisible(NULL))
    }

    Try <- 0
    while (Try < 6) {
      Try <- Try + 1
      Model <- Hmsc::sampleMcmc(
        hM = ecokit::load_as(Model_Info$M_Init_Path[ID]),
        samples = Model_Info$M_samples[ID],
        thin = Model_Info$M_thin[ID],
        transient = Model_Info$M_transient[ID],
        nChains = MCMC_n_chains, verbose = MCMC_verbose, engine = "HPC")

      if (to_JSON) {
        Model <- jsonify::to_json(Model)
      }

      saveRDS(Model, file = CurrPath)

      if (file.exists(CurrPath) && ecokit::check_rds(CurrPath)) {
        break
      }
    }
    invisible(gc())
    return(invisible(NULL))
  }

  # Implement `InitFitFun` function: start sampling and save output files
  if (n_cores > 1) {

    ecokit::set_parallel(
      n_cores = min(n_cores, nrow(Model_Info)), level = 1L,
      future_max_size = 800L, strategy = strategy)
    withr::defer(future::plan("sequential", gc = TRUE))

    Model_Process <- future.apply::future_lapply(
      X = seq_len(nrow(Model_Info)),
      FUN = InitFitFun, future.scheduling = Inf, future.seed = TRUE,
      future.globals = c(
        "InitFitFun", "Model_Info", "overwrite_rds",
        "MCMC_verbose", "MCMC_n_chains", "to_JSON"),
      future.packages = pkg_to_export)

    ecokit::set_parallel(stop_cluster = TRUE, level = 1L)
    future::plan("sequential", gc = TRUE)

  } else {
    Model_Process <- purrr::map(
      .x = seq_len(nrow(Model_Info)), .f = InitFitFun)
  }
  rm(Model_Process, envir = environment())

  # Which models failed to be exported as RDS files after 5 trials
  Failed2Export <- dplyr::filter(Model_Info, !file.exists(M4HPC_Path))
  if (nrow(Failed2Export) == 0) {
    ecokit::cat_time(
      "All model variants were exported as RDS files", level = 1L,
      verbose = verbose_progress)
  } else {
    ecokit::cat_time(
      paste0(
        nrow(Failed2Export),
        " model variants failed to be exported to rds files after 5 tries."),
      level = 1L, verbose = verbose_progress)
    save(Failed2Export, file = fs::path(path_model, "Failed2Export.RData"))
    readr::write_tsv(
      x = Failed2Export, file = fs::path(path_model, "Failed2Export.txt"))
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Prepare Hmsc-HPC fitting commands -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time(
    "Prepare Hmsc-HPC fitting commands", verbose = verbose_progress)

  Model_Info <- Model_Info %>%
    dplyr::mutate(Chain = list(seq_len(MCMC_n_chains))) %>%
    tidyr::unnest_longer("Chain") %>%
    dplyr::arrange(M_Name_Fit) %>%
    dplyr::mutate(
      M_Chain = purrr::pmap(
        .l = list(M_Name_Fit, M4HPC_Path, Chain,
                  M_transient, M_samples, M_thin),
        .f = function(M_Name_Fit, M4HPC_Path, Chain,
                      M_transient, M_samples, M_thin) {

          # Turn off scientific notation
          withr::local_options(list(scipen = 999))

          # Input model
          Path_Model2 <- fs::path(
            path_model, "InitMod4HPC", basename(M4HPC_Path))

          # Path for posterior sampling
          Path_Post <- fs::path(
            path_model, "Model_Fitting_HPC",
            paste0(M_Name_Fit, "_Chain", Chain, "_post.rds"))

          # Path for progress
          Path_Prog <- fs::path(
            path_model, "Model_Fitting_HPC",
            paste0(M_Name_Fit, "_Chain", Chain, "_Progress.txt"))

          Post_Missing <- isFALSE(file.exists(Path_Post))

          # File path for the python script
          Path_Model2_4cmd <- ecokit::normalize_path(Path_Model2)
          Path_Post_4cmd <- ecokit::normalize_path(Path_Post)
          Path_Prog_4cmd <- ecokit::normalize_path(Path_Prog)
          Path_Python <- ecokit::normalize_path(Path_Python)

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

          Command_HPC <- paste0(
            # Not needed now as this now added to the `setup-env.sh` file
            # "export TF_CPP_MIN_LOG_LEVEL=3; export TF_ENABLE_ONEDNN_OPTS=0; ",

            "/usr/bin/time -v ", # nolint: absolute_paths_linter

            # Not needed as the python path is exported - check `setup-env.sh`
            # Path_Python,

            "python3 -m hmsc.run_gibbs_sampler",
            " --input ", Path_Model2_4cmd,
            " --output ", Path_Post_4cmd,
            " --samples ", M_samples,
            " --transient ", M_transient,
            " --thin ", M_thin,
            " --verbose ", MCMC_verbose,
            " --chain ", (Chain - 1),
            " --fp ", precision,
            " >& ", Path_Prog_4cmd)

          Command_WS <- paste0(
            "set TF_CPP_MIN_LOG_LEVEL=3 && set TF_ENABLE_ONEDNN_OPTS=0 && ",
            Path_Python,
            " -m hmsc.run_gibbs_sampler",
            " --input ", Path_Model2_4cmd,
            " --output ", Path_Post_4cmd,
            " --samples ", M_samples,
            " --transient ", M_transient,
            " --thin ", M_thin,
            " --verbose ", MCMC_verbose,
            " --chain ", (Chain - 1),
            " --fp ", precision,
            " > ", Path_Prog_4cmd,
            " 2>&1")

          return(
            list(
              M4HPC_Path_LUMI = Path_Model2, Post_Path = Path_Post,
              Post_Missing = Post_Missing, Path_ModProg = Path_Prog,
              Command_HPC = Command_HPC, Command_WS = Command_WS))

        })) %>%
    tidyr::unnest_wider("M_Chain")

  # # |||||||||||||||||||||||||||||||||||
  # # Skip fitted models -----
  # # |||||||||||||||||||||||||||||||||||

  if (skip_fitted) {
    ecokit::cat_time("Skip fitted models", verbose = verbose_progress)
    Models2Fit_HPC <- dplyr::filter(Model_Info, Post_Missing) %>%
      dplyr::pull(Command_HPC) %>%
      unlist()
    Models2Fit_WS <- dplyr::filter(Model_Info, Post_Missing) %>%
      dplyr::pull(Command_WS) %>%
      unlist()
  } else {
    Models2Fit_HPC <- unlist(dplyr::pull(Model_Info, Command_HPC))
    Models2Fit_WS <- unlist(dplyr::pull(Model_Info, Command_WS))
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
    description = fs::path(path_model, "Commands_All_Windows.txt"),
    open = "wb")
  on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
  cat(Models2Fit_WS, sep = "\n", append = FALSE, file = f)
  close(f)


  ecokit::cat_time(
    "Save fitting commands for HPC", level = 1L, verbose = verbose_progress)
  NJobs <- length(Models2Fit_HPC)

  if (NJobs > n_array_jobs) {
    NSplits <- ceiling((NJobs / n_array_jobs))
    IDs <- ecokit::split_vector(vector = seq_len(NJobs), n_splits = NSplits)
  } else {
    NSplits <- 1
    IDs <- list(seq_len(NJobs))
  }

  # Save all fitting commands to single file
  ecokit::cat_time(
    "Save all fitting commands to single file",
    level = 1L, verbose = verbose_progress)
  f <- file(
    description = fs::path(path_model, "Commands_All.txt"),
    open = "wb")
  on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
  cat(Models2Fit_HPC, sep = "\n", append = FALSE, file = f)
  close(f)

  # Save model fitting commands for batch SLURM jobs
  ecokit::cat_time(
    "Save model fitting commands for batch SLURM jobs", level = 1L,
    verbose = verbose_progress)
  ecokit::cat_time(
    paste0("Models will be fitted in ", NSplits, " SLURM job(s)"),
    level = 2L, cat_timestamp = FALSE, verbose = verbose_progress)

  purrr::walk(
    .x = seq_len(NSplits),
    .f = function(x) {

      if (NSplits > 1) {
        CommandFile <- fs::path(
          path_model, paste0("Commands2Fit_", x, ".txt"))
      } else {
        CommandFile <- fs::path(path_model, "Commands2Fit.txt")
      }

      # create connection to SLURM file. This is better than using sink to have
      # a platform independent file (here, to maintain a linux-like new line
      # ending)
      f <- file(CommandFile, open = "wb")
      on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
      cat(Models2Fit_HPC[IDs[[x]]], sep = "\n", append = FALSE, file = f)
      on.exit(close(f))
    })

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Save data to disk -----
  # # |||||||||||||||||||||||||||||||||||

  ecokit::cat_time("Save data to disk", verbose = verbose_progress)

  SetChainName <- function(Obj, Chain) {
    if (is.null(Obj) || is.null(Chain)) {
      ecokit::stop_ctx(
        "Obj and Chain cannot be empty", Obj = Obj, Chain = Chain,
        include_backtrace = TRUE)
    }
    Obj %>%
      unlist() %>%
      as.vector() %>%
      stats::setNames(paste0("Chain", unlist(Chain)))
  }

  Model_Info <- Model_Info %>%
    tidyr::nest(
      Post_Path = Post_Path, Path_ModProg = Path_ModProg,
      Chain = Chain, Command_HPC = Command_HPC, Command_WS = Command_WS,
      Post_Missing = Post_Missing) %>%
    dplyr::mutate(
      Post_Path = purrr::map2(Post_Path, Chain, SetChainName),
      Chain = purrr::map2(Chain, Chain, SetChainName),
      Command_HPC = purrr::map2(Command_HPC, Chain, SetChainName),
      Command_WS = purrr::map2(Command_WS, Chain, SetChainName),
      Path_ModProg = purrr::map2(Path_ModProg, Chain, SetChainName),
      Post_Missing = purrr::map2(Post_Missing, Chain, SetChainName),
      Post_Aligned = NA)

  save(Model_Info, file = Path_ModelDT)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Prepare SLURM file ------
  # # |||||||||||||||||||||||||||||||||||

  if (SLURM_prepare) {
    ecokit::cat_time("Preparing SLURM file", verbose = verbose_progress)
    if (is.null(job_name)) {
      job_name <- stringr::str_remove_all(
        basename(path_model), paste0("_", HabVal))
    }

    IASDT.R::mod_SLURM(
      model_dir = path_model, job_name = job_name,
      memory_per_cpu = memory_per_cpu, job_runtime = job_runtime,
      env_file = env_file, path_Hmsc = path_Hmsc, ...)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Save PA maps to disk ------

  if (is.null(CV_fit_inherit) && (is.null(CV_fit_fold) || CV_fit_fold == 1)) {

    ecokit::cat_time("Save PA maps to disk", verbose = verbose_progress)

    # Mask layer for grid cells used in the models
    Model_Mask <- ecokit::load_as(Path_GridR, unwrap_r = TRUE) %>% # nolint: object_name_linter
      terra::rasterize(DT_xy, .)
    # Whether to use masked (exclude_cultivated) or full PA
    PA_Layer <- dplyr::if_else(exclude_cultivated, "PA_Masked", "PA")      # nolint: object_name_linter

    ecokit::cat_time(
      "Processing and exporting maps as tif files", level = 1L,
      verbose = verbose_progress)

    Mod_PA <- SpSummary %>%
      # Select only species name and ID
      dplyr::select(ias_id = IAS_ID, Species_File) %>%
      dplyr::mutate(
        # Path storing PA maps as raster files
        File = fs::path(Path_PA, "PA_RData", paste0(Species_File, "_PA.RData")),
        Map_Sp = purrr::map2(
          .x = File, .y = ias_id,
          .f = ~ {
            # Path for storing PA map - full EU extent
            PA_file <- fs::path(Path_Dist, paste0(.y, "_full.tif"))
            # Load PA map
            PA <- ecokit::load_as(.x, unwrap_r = TRUE) %>%
              magrittr::extract2(PA_Layer)
            # Save PA map
            terra::writeRaster(PA, PA_file, overwrite = TRUE)

            # Path for storing PA map - only in modelling grid cells
            PA_model_file <- fs::path(Path_Dist, paste0(.y, "_model.tif"))
            # mask PA map to modelling grid cells
            PA_model <- terra::mask(PA, Model_Mask)
            # Save masked PA map
            terra::writeRaster(PA_model, PA_model_file, overwrite = TRUE)

            tibble::tibble(
              PA = list(terra::wrap(PA)),
              PA_file = PA_file,
              PA_model = list(terra::wrap(PA_model)),
              PA_model_file = PA_model_file)
          })) %>%
      tidyr::unnest(cols = "Map_Sp") %>%
      dplyr::select(-Species_File, -File)

    ecokit::cat_time(
      "Calculate species richness - full", level = 1L, verbose = verbose_progress)
    SR_file <- fs::path(Path_Dist, "SR_full.tif")
    SR <- purrr::map(Mod_PA$PA, terra::unwrap) %>%
      terra::rast() %>%
      sum(na.rm = TRUE)
    terra::writeRaster(SR, SR_file, overwrite = TRUE)

    ecokit::cat_time(
      "Calculate species richness - modelling", level = 1L,
      verbose = verbose_progress)
    SR_model_file <- fs::path(Path_Dist, "SR_model.tif")
    SR_model <- purrr::map(Mod_PA$PA_model, terra::unwrap) %>%
      terra::rast() %>%
      sum(na.rm = TRUE)
    terra::writeRaster(SR_model, SR_model_file, overwrite = TRUE)

    # Bind SR and PA in the same tibble
    Mod_PA <- tibble::tribble(
      ~ias_id, ~PA, ~PA_file, ~PA_model, ~PA_model_file,
      "SR", terra::unwrap(SR), SR_file,
      terra::unwrap(SR_model), SR_model_file) %>%
      dplyr::bind_rows(Mod_PA)

    ecokit::cat_time("Save maps as RData", level = 1L, verbose = verbose_progress)
    ecokit::save_as(
      object = Mod_PA, object_name = "PA_with_maps",
      out_path = fs::path(Path_Dist, "PA_with_maps.RData"))

    ecokit::cat_time(
      "Save maps without maps as RData", level = 1L, verbose = verbose_progress)
    Mod_PA <- dplyr::select(Mod_PA, -PA, -PA_model)
    ecokit::save_as(
      object = Mod_PA, object_name = "PA",
      out_path = fs::path(Path_Dist, "PA.RData"))

    ecokit::cat_time(
      "Save only paths text file", level = 1L, verbose = verbose_progress)
    Mod_PA %>%
      dplyr::mutate(
        PA_file = basename(PA_file),
        PA_model_file = basename(PA_model_file)) %>%
      utils::write.table(
        sep = "\t", row.names = FALSE, col.names = TRUE,
        file = fs::path(Path_Dist, "PA.txt"), quote = FALSE,
        fileEncoding = "UTF-8")

    # Save maps as tar file
    ecokit::cat_time(
      "Save maps as single tar file", level = 1L, verbose = verbose_progress)
    # Path of the tar file
    tar_file <- fs::path(Path_Dist, "PA_maps.tar")
    # list of files to tar
    Files2Tar <- c(      # nolint: object_name_linter
      Mod_PA$PA_file, Mod_PA$PA_model_file,
      fs::path(Path_Dist, "PA.txt")) %>%
      basename() %>%
      paste(collapse = " ")
    tar_command <- stringr::str_glue(
      'tar -cf "{tar_file}" -C "{Path_Dist}" {Files2Tar}')
    invisible(system(tar_command))

    # Change the permission of the tar file
    Sys.chmod(tar_file, "755", use_umask = FALSE)

  }

  # # |||||||||||||||||||||||||||||||||||
  # # Save small datasets prepared in the function ------
  # # |||||||||||||||||||||||||||||||||||

  DT_Split <- list(
    DT_All = DT_All, DT_y = DT_y, Form_x = Form_x, DT_x = DT_x, XVars = XVars,
    DT_CV = DT_CV, use_phylo_tree = use_phylo_tree, plant.tree = plant.tree,
    Tree = Tree, studyDesign = studyDesign, DT_xy = DT_xy,
    GPP_Knots = GPP_Knots)

  ecokit::save_as(
    object = DT_Split, object_name = "ModDT_subset",
    out_path = fs::path(path_model, "ModDT_subset.RData"))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_diff(
    init_time = .start_time, prefix = "Processing modelling data took ",
    verbose = verbose_progress)

  return(invisible(NULL))
}
