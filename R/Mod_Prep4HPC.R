# # |------------------------------------------------------------------------| #
# Mod_Prep4HPC ----
## |------------------------------------------------------------------------| #

#' Prepare initial models in R for model fitting by Hmsc-HPC
#'
#' This function prepares initial models in R for model fitting by Hmsc-HPC. It
#' involves data preparation, define spatial block cross-validation folds, model
#' initialization, and generating commands for running models on HPC. It
#' supports parallel processing, options to include/not include phylogenetic
#' tree data. The models are fitted with Gaussian Predictive Process (GPP; see
#' [Tikhonov et al.](https://doi.org/10.1002/ecy.2929)) for more details) using
#' the [Hmsc-HPC](https://doi.org/10.1101/2024.02.13.580046) extension.
#'
#' @param Hab_Abb Character. Habitat abbreviation indicating the specific
#'   `SynHab` habitat type to prepare data for. Valid values include "0", "1",
#'   "2", "3", "4a", "4b", "10", "12a", "12b". If `Hab_Abb` = "0", data is
#'   prepared irrespective of the habitat type.
#' @param Path_Data String specifying the path where modeling data is read from.
#'   Only effective if `PrepareData` = `FALSE`; see below.
#' @param Path_Model String (without trailing slash) specifying the path where
#'   all output, including models to be fitted, will be saved.
#' @param MinPresGrids Integer. Indicating the minimum number of presence grid
#'   cells per species for a species to be used in the model. Default: 50.
#' @param EnvFile String specifying the path to read environment variables from,
#'   with a default value of `.env`.
#' @param PrepareData Logical indicating whether to prepare input data or load
#'   it from disk. Defaults to `TRUE` which means the input data will be
#'   prepared using the [IASDT.R::Mod_PrepData] function.
#' @param GPP_Dists Integer specifying the distance in kilometers for both the
#'   distance between knots and the minimum distance of a knot to the nearest
#'   data point. The GPP knots are prepared by the [IASDT.R::PrepKnots]
#'   function. The same value will be used for the `knotDist` and
#'   `minKnotDist`	arguments of the [Hmsc::constructKnots] function.
#' @param GPP_Save Logical indicating whether to save the resulted knots as
#'   `RData` Default: `TRUE`.
#' @param GPP_Plot Logical indicating whether to plot the coordinates of the
#'   sampling units and the knots in a pdf file. Default: `TRUE`.
#' @param rLMinLF,rLMaxLF Integer. Minimum and maximum number of latent factors
#'   used for the spatial random effect. See [PrepKnots] and [Hmsc::setPriors]
#'   for more details.
#' @param BioVars Vector of strings specifying bioclimatic variables to be used
#'   in the model. Default value: `NULL`, which means to use the following
#'   variables: `bio4`, `bio6`, `bio8`, `bio12`, `bio15`, and `bio18`.
#' @param EffortsAsPredictor Logical indicating whether to include the (log10)
#'   sampling efforts as predictor to the model. Default: `TRUE`.
#' @param RoadRailAsPredictor Logical indicating whether to include the (log10)
#'   sum of road and railway intensity as predictor to the model. Default:
#'   `TRUE`.
#' @param HabAsPredictor Logical indicating whether to include the (log10)
#'   percentage coverage of respective habitat type per grid cell as predictor
#'   to the model. Default: `TRUE`. Only valid if `Hab_Abb` not equals to "0".
#' @param NspPerGrid Integer. Indicating the minimum number of species per grid
#'   cell for a grid cell to be include in the analysis. Default to 0 resulting
#'   in exclusion of any grid cell with no species presence. This parameter can
#'   be set to `NULL` to include all grid cells irrespective of the number of
#'   species.
#' @param NGrids For `CV_Dist` cross-validation strategy, how many grid cells in
#'   both directions to be used in cross-validation. See [IASDT.R::GetCV] for
#'   more details.
#' @param NFolds Number of cross-validation folds. Default: 4. See
#'   [IASDT.R::GetCV] for more details.
#' @param NR,NC Integer, the number of rows and columns to divide the spatial
#'   area into. Defaults to 2 row and 2 columns. See [IASDT.R::GetCV] for more
#'   details.
#' @param PlotCV Logical. Indicating whether to plot the block cross-validation
#'   folds.
#' @param PhyloTree,NoPhyloTree Logical indicating whether to fit model variants
#'   with or without phylogenetic trees, respectively. The default of both
#'   arguments is `TRUE`, which means to fit a model variant with the respective
#'   option. If both `PhyloTree` and `NoPhyloTree` are `TRUE` (Default), models
#'   for both options will be fitted. At least one of `PhyloTree` and
#'   `NoPhyloTree` should be `TRUE`.
#' @param OverwriteInitMod Logical. Indicating whether to overwrite previously
#'   exported RDS files for initial models. Default: `TRUE`.
#' @param NParallel Integer specifying the number of parallel cores for
#'   parallelization. Default: 8 cores.
#' @param nChains Integer specifying the number of model chains. Default: 4.
#' @param thin Integer specifying the value(s) for thinning in MCMC sampling. If
#'   more than one value is provided, a separate model will be fitted at each
#'   value of thinning.
#' @param samples Integer specifying the value(s) for the number of MCMC
#'   samples. If more than one value is provided, a separate model will be
#'   fitted at each value of number of samples.
#' @param transientFactor Integer specifying the transient multiplication
#'   factor. The value of `transient` will equal  the multiplication of
#'   `transientFactor` and `thin`. Default: 300.
#' @param verbose Integer specifying how often the results of the MCMC sampling
#'   should be reported. Default: `200`.
#' @param SkipFitted Logical indicating whether to skip already fitted models.
#'   Default: `TRUE`.
#' @param MaxJobCounts Integer specifying the maximum allowed number of array
#'   jobs per SLURM file. Default: 210. See
#'   [here](https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/partitions)
#'   for more details.
#' @param ModelCountry String or vector of strings specifying the country or
#'   countries to filter observations by. Default: `NULL`, which means prepare
#'   data for the whole Europe.
#' @param MinPresPerCountry Integer specifying the minimum number of grid cells
#'   for the selected country/countries for species to be considered in the
#'   models. Effective only if a valid `ModelCountry` is provided. Default: 50.
#' @param VerboseProgress Logical indicating whether to show messages for the
#'   progress of creating files. Default: `FALSE`.
#' @param FromHPC Logical indicating whether the work is being done from HPC, to
#'   adjust file paths accordingly. Default: `TRUE`.
#' @param PrepSLURM Logical indicating whether to prepare SLURM command files.
#'   If `TRUE` (default), the SLURM commands will be saved to disk using the
#'   [IASDT.R::Mod_SLURM] function.
#' @param MemPerCpu String specifying the memory per CPU for the SLURM job. This
#'   value will be assigned to the `#SBATCH --mem-per-cpu=` SLURM argument.
#'   Example: "32G" to request 32 gigabyte. Only effective if `PrepSLURM =
#'   TRUE`.
#' @param Time String specifying the requested time for each job in the SLURM
#'   bash arrays. Example: "01:00:00" to request an hour. Only effective if
#'   `PrepSLURM = TRUE`.
#' @param JobName String specifying the name of the submitted job(s) for SLURM.
#'   If `NULL` (Default), the job name will be prepared based on the folder path
#'   and the `Hab_Abb` value. Only effective if `PrepSLURM = TRUE`.
#' @param Path_Hmsc String specifying the path for the Hmsc-HPC. This will be
#'   provided as the `Path_Hmsc` argument of the [IASDT.R::Mod_SLURM] function.
#' @param ToJSON Logical indicating whether to convert unfitted models to JSON
#'   before saving to RDS file. Default: `FALSE`.
#' @param ... Additional parameters provided to the [IASDT.R::Mod_SLURM]
#'   function.
#' @name Mod_Prep4HPC
#' @author Ahmed El-Gabbas
#' @return The function is used for its side effects of preparing data and
#'   models for HPC and does not return any value.
#' @details The function provides options for:
#'
#' - for which habitat types the models will be fitted.
#' - selection of species based on minimum number of presence-grid cells
#'   (`MinPresGrids`)
#' - optionally model fitting on specified list of countries: (`ModelCountry`
#'   and `MinPresPerCountry`)
#' - whether to exclude grid cells with few species (`NspPerGrid`)
#' - number of cross-validation folds
#' - options for whether or not to include phylogenetic information to the model
#' - different values for knot distance for GPP (`GPP_Dists`)
#' - which bioclimatic variables to be uses in the models (`BioVars`)
#' - whether to include sampling efforts `EffortsAsPredictor`, percentage of
#'   respective habitat type per grid cell `HabAsPredictor`, and railway and
#'   road intensity per grid cell `RoadRailAsPredictor`
#' - Hmsc options (`nChains`, `thin`, `samples`, `transientFactor`, and `verbose`)
#' - prepare SLURM commands (`PrepSLURM`) and some specifications (e.g.
#'   `MaxJobCounts`, `MemPerCpu`, `Time`, `JobName`)
#'
#' The function reads the following environment variables:
#'   - **`DP_R_Grid`** (if `FromHPC = TRUE`) or
#'    **`DP_R_Grid_Local`** (if `FromHPC = FALSE`). The function reads
#' the content of the `Grid_10_Land_Crop.RData` file from this path.
#'   - **`DP_R_Path_Python`**: Python path on LUMI
#'   - **`DP_R_TaxaInfo`** or **`DP_R_TaxaInfo_Local`** for the location of the
#' `Species_List_ID.txt` file representing species information.
#'   - **`DP_R_EUBound_sf`** or **`DP_R_EUBound_sf_Local`** for the path of the
#' `RData` file containing the country boundaries (`sf` object)
#'   - **`DP_R_PA`** or **`DP_R_PA_Local`**: The function reads the contents of
#'  the `Sp_PA_Summary_DF.RData` file from this path

#' @export

Mod_Prep4HPC <- function(
    Hab_Abb = NULL, Path_Data = NULL, Path_Model = NULL, MinPresGrids = 50,
    EnvFile = ".env", PrepareData = TRUE, GPP_Dists = NULL, GPP_Save = TRUE,
    GPP_Plot = TRUE, rLMinLF = NULL, rLMaxLF = NULL,
    BioVars = c("bio4", "bio6", "bio8", "bio12", "bio15", "bio18"),
    EffortsAsPredictor = TRUE, RoadRailAsPredictor = TRUE,
    HabAsPredictor = TRUE, NspPerGrid = 0, NFolds = 4,
    NGrids = 20, NR = 2, NC = 2, PlotCV = TRUE, PhyloTree = TRUE,
    NoPhyloTree = TRUE, OverwriteInitMod = TRUE, NParallel = 8, nChains = 4,
    thin = NULL, samples = NULL, transientFactor = 300, verbose = 200,
    SkipFitted = TRUE, MaxJobCounts = 210, ModelCountry = NULL,
    MinPresPerCountry = 50, VerboseProgress = FALSE, FromHPC = TRUE,
    PrepSLURM = TRUE, MemPerCpu = NULL, Time = NULL, JobName = NULL,
    Path_Hmsc = NULL, ToJSON = FALSE, ...) {

  .StartTime <- lubridate::now(tzone = "CET")

  # # |||||||||||||||||||||||||||||||||||
  # # Initial checking -----
  # # |||||||||||||||||||||||||||||||||||

  if (is.null(Path_Model) || is.null(MinPresGrids)) {
    stop("Path_Model and MinPresGrids cannot be empty")
  }

  if (!PrepareData && is.null(Path_Data)) {
    stop("Path_Data can not be empty if PrepareData = FALSE")
  }

  if (is.null(thin) || is.null(samples) || is.null(GPP_Dists)) {
    stop("thin, samples and GPP_Dists cannot be empty")
  }

  if (is.null(MemPerCpu) || is.null(Time) || is.null(Path_Hmsc)) {
    stop("MemPerCpu, Time and Path_Hmsc cannot be empty")
  }

  if (!all(is.numeric(GPP_Dists)) || any(GPP_Dists <= 0)) {
    stop("GPP_Dists should be numeric and greater than zero")
  }

  if (!all(is.numeric(samples)) || any(samples <= 0)) {
    stop("samples should be numeric and greater than zero")
  }

  if (!all(is.numeric(thin)) || any(thin <= 0)) {
    stop("thin should be numeric and greater than zero")
  }

  if (!all(is.numeric(MinPresGrids)) || any(MinPresGrids <= 0)) {
    stop("MinPresGrids should be numeric and greater than zero")
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  NCells <- Sp <- IAS_ID <- x <- y <- dplyr <- sf <- Country <-
    Hmsc <- jsonify <- magrittr <- M_thin <- rL <- M_Name_init <- rL2 <-
    M_samples <- M4HPC_Path <- M_transient <- M_Init_Path <- M_Name_Fit <-
    Chain <- Post_Missing <- Command_HPC <- Command_WS <- Post_Path <-
    Path_ModProg <- TaxaInfoFile <- Path_Python <- Path_Grid <- EU_Bound <-
    Path_PA <- NULL

  if (magrittr::not(VerboseProgress)) {
    sink(file = nullfile())
    on.exit(sink(), add = TRUE)
  }

  IASDT.R::CatSep(Rep = 1, Extra1 = 1, Extra2 = 0, Char = "=")
  IASDT.R::CatTime("Preparing data for Hmsc-HPC models")
  IASDT.R::CatSep(Rep = 1, Extra1 = 0, Extra2 = 1, Char = "=")

  # # |||||||||||||||||||||||||||||||||||
  # # Load/check environment variables -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Load/check environment variables")

  if (magrittr::not(file.exists(EnvFile))) {
    stop(paste0(
      "Path for environment variables: ", EnvFile, " was not found"))
  }

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Python", "DP_R_Path_Python", FALSE, FALSE,
      "TaxaInfoFile", "DP_R_TaxaInfo", FALSE, TRUE,
      "EU_Bound", "DP_R_EUBound_sf", FALSE, TRUE,
      "Path_Grid", "DP_R_Grid", TRUE, FALSE,
      "Path_PA", "DP_R_PA", TRUE, FALSE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Python", "DP_R_Path_Python", FALSE, FALSE,
      "TaxaInfoFile", "DP_R_TaxaInfo_Local", FALSE, TRUE,
      "EU_Bound", "DP_R_EUBound_sf_Local", FALSE, TRUE,
      "Path_Grid", "DP_R_Grid_Local", TRUE, FALSE,
      "Path_PA", "DP_R_PA_Local", TRUE, FALSE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  if (GPP_Plot) {
    Path_GridR <- file.path(Path_Grid, "Grid_10_Land_Crop.RData")
    if (magrittr::not(file.exists(Path_GridR))) {
      stop(paste0("Path for the Europe boundaries does not exist: ", Path_GridR))
    }
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Check input arguments ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Checking input arguments")

  AllArgs <- ls()
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  CharArgs <- c(
    "Hab_Abb", "Path_Model", "Path_Hmsc")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = CharArgs, Type = "character")

  LogicArgs <- c("GPP_Save", "GPP_Plot", "PhyloTree",
                 "NoPhyloTree", "VerboseProgress", "ToJSON")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = LogicArgs, Type = "logical")

  NumericArgs <- c(
    "GPP_Dists", "NParallel", "nChains", "thin", "samples", "verbose",
    "MinPresGrids", "transientFactor", "NFolds", "NGrids", "NR", "NC")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = NumericArgs, Type = "numeric")

  if (PrepSLURM) {
    IASDT.R::CheckArgs(
      AllArgs = AllArgs, Args = c("MemPerCpu", "Time"), Type = "character")
  }

  # Phylogenetic tree options
  if (PhyloTree == FALSE && NoPhyloTree == FALSE) {
    stop("At least one of PhyloTree or NoPhyloTree has to be true")
  }

  NumArgsInvalid <- purrr::map_lgl(.x = NumericArgs, .f = ~all(get(.x) < 1))
  if (any(NumArgsInvalid)) {
    paste0(
      "The following parameter(s) can not be < 1\n  >>  ",
      paste0(NumericArgs[NumArgsInvalid], collapse = " | ")) %>%
      stop(call. = FALSE)
  }

  rm(AllArgs, CharArgs, LogicArgs, NumericArgs)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # File paths - Creating missing paths ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("File paths - Creating missing paths")
  fs::dir_create(file.path(Path_Model, "InitMod_HPC"))
  fs::dir_create(file.path(Path_Model, "Model_Fitting_HPC"))
  # Also create directory for SLURM outputs
  fs::dir_create(file.path(Path_Model, "Model_Fitting_HPC", "SLURM_Results"))
  Path_ModelDT <- file.path(Path_Model, "Model_Info.RData")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Prepare list of predictors -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Prepare list of predictors")
  if (is.null(BioVars)) {
    BioVars <- c(
      "bio4", "bio6", "bio8", "bio12", "bio15", "bio18")
  }

  XVars <- BioVars

  if (EffortsAsPredictor) {
    XVars <- c(XVars, "RoadRailLog")
  }
  if (RoadRailAsPredictor) {
    XVars <- c(XVars, "BiasLog")
  }
  if (Hab_Abb != "0" && HabAsPredictor) {
    XVars <- c(XVars, "HabLog")
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Loading/preparing input data -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Loading/preparing input data")
  ValidHabAbbs <- c(0:3, "4a", "4b", 10, "12a", "12b")
  if (magrittr::not(as.character(Hab_Abb) %in% ValidHabAbbs)) {
    stop(paste0("Hab_Abb has to be one of the following:\n >> ",
                paste0(ValidHabAbbs, collapse = " | ")))
  }

  HabVal <- c(
    "0_All", "1_Forests", "2_Open_forests", "3_Scrub",
    "4a_Natural_grasslands", "4b_Human_maintained_grasslands",
    "10_Wetland", "12a_Ruderal_habitats", "12b_Agricultural_habitats") %>%
    stringr::str_subset(paste0("^", as.character(Hab_Abb), "_"))

  if (PrepareData) {
    IASDT.R::CatSep(Rep = 1, Extra1 = 1, Extra2 = 0)
    IASDT.R::CatTime("Preparing input data using IASDT.R::Mod_PrepData")
    IASDT.R::CatSep(Rep = 1, Extra1 = 0, Extra2 = 0)

    DT_All <- IASDT.R::Mod_PrepData(
      Hab_Abb = Hab_Abb, MinPresGrids = MinPresGrids, EnvFile = EnvFile,
      BioVars = BioVars, ReturnData = TRUE, OutputPath = Path_Model,
      VerboseProgress = VerboseProgress, FromHPC = FromHPC)

    IASDT.R::CatSep(Rep = 1, Extra1 = 0, Extra2 = 1)
  } else {
    IASDT.R::CatTime("   >>>   Loading input data")
    Path_Data2 <- file.path(
      Path_Data,
      paste0(
        "ModDT_", MinPresGrids, "Grids_",
        stringr::str_remove(HabVal, "Hab_"), ".RData"))

    if (!file.exists(Path_Data2)) {
      stop(paste0("Data file ", Path_Data2, " does not exist"))
    }

    DT_All <- IASDT.R::LoadAs(Path_Data2)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Subsetting study area -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Subsetting study area")

  if (magrittr::not(is.null(ModelCountry))) {

    ValidCountries <- ModelCountry %in% unique(DT_All$Country)

    if (magrittr::not(all(ValidCountries))) {
      stop(paste0(
        "The following are invalid country names: ",
        paste0(ModelCountry[!ValidCountries], collapse = " & ")))
    }

    IASDT.R::CatTime(
      paste0("   >>>   Subsetting data to ",
             paste0(ModelCountry, collapse = " & ")))

    Sample_ExclSp <- dplyr::filter(DT_All, Country %in% ModelCountry) %>%
      dplyr::summarise(
        dplyr::across(tidyselect::starts_with("Sp_"), sum)) %>%
      tidyr::pivot_longer(
        cols = tidyselect::everything(),
        names_to = "Sp", values_to = "NCells") %>%
      dplyr::filter(NCells < MinPresPerCountry) %>%
      dplyr::pull(Sp)

    DT_All <- dplyr::filter(DT_All, Country %in% ModelCountry) %>%
      dplyr::select(
        -tidyselect::all_of(Sample_ExclSp), -"CellCode", -"Country", -"Country_Nearest")

  } else {
    IASDT.R::CatTime("   >>>   No data subsetting was implemented")
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Exclude grid cells with no / low number of presences -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Exclude grid cells with no / low number of presences")

  if (is.null(NspPerGrid)) {
    IASDT.R::CatTime("   >>>   All grid cells will be included in the analysis")
    SubTitle <- "All grid cells will be used in the models"

  } else {

    if (magrittr::not(is.numeric(NspPerGrid)) || NspPerGrid < 0) {
      stop(
        "`NspPerGrid` has to be either `NULL` or integer >= 0", call. = FALSE)
    }

    SubTitle <- paste0(
      "Data from grid cells with &#8805; ", NspPerGrid,
      " species will be used in the models")

    EmptyGridsID <- dplyr::select(DT_All, tidyselect::starts_with("Sp_")) %>%
      rowSums() %>%
      magrittr::is_weakly_less_than(as.integer(NspPerGrid)) %>%
      which() %>%
      magrittr::multiply_by(-1)

    if (length(EmptyGridsID) > 0) {
      IASDT.R::CatTime(
        paste0("   >>>   Excluding grid cells with < ", NspPerGrid))
      DT_All <- dplyr::slice(DT_All, EmptyGridsID)
    }
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime(
    "   >>>   Plotting grid cells used / excluded in the models")

  GridR <- terra::unwrap(IASDT.R::LoadAs(Path_GridR))
  GridsIn <- terra::rasterize(
    x = as.matrix(DT_All[, c("x", "y")]), y = GridR)
  GridR <- GridR * 2
  GridsOut <- terra::classify(terra::cover(GridsIn, GridR), cbind(1, NA))
  GridsInOut <- as.factor(sum(GridsIn, GridsOut, na.rm = TRUE))

  Limits <- terra::trim(GridR) %>%
    terra::ext() %>%
    as.vector()

  GridsInOutPlot <- ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = GridsInOut) +
    ggplot2::geom_sf() +
    ggplot2::scale_fill_manual(
      values = c("blue", "orange"), na.value = "transparent") +
    ggplot2::labs(
      title = paste0(
        "Grid cells to be used (blue) or excluded (orange) ",
        "in the model fitting"),
      subtitle = SubTitle) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = Limits[3:4]) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = Limits[1:2]) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(0.05, 0, 0, 0, "cm"),
      plot.title = ggplot2::element_text(
        size = 16, color = "blue", face = "bold", hjust = 0,
        margin = ggplot2::margin(0, 0, 0.1, 0, "cm")),
      plot.subtitle = ggtext::element_markdown(size = 14, hjust = 0),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank())

  ggplot2::ggsave(
    plot = GridsInOutPlot, width = 24, height = 25, units = "cm", dpi = 600,
    filename = file.path(Path_Model, "GridsInOut.jpeg"))

  rm(GridR, GridsIn, GridsOut, GridsInOut, Limits, SubTitle, GridsInOutPlot)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Cross-validation ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Prepare cross-validation folds")

  DT_All <- IASDT.R::GetCV(
    DT = DT_All, EnvFile = EnvFile, XVars = XVars,
    NFolds = NFolds, FromHPC = FromHPC, NGrids = NGrids, NR = NR,
    NC = NC, OutPath = Path_Model, PlotCV = PlotCV)


  if (Hab_Abb == "0") {
    OutObjName <- "ModDT_0_All_subset"
  } else {
    OutObjName <- c(
      "1_Forests", "2_Open_forests", "3_Scrub",
      "4a_Natural_grasslands", "4b_Human_maintained_grasslands",
      "10_Wetland", "12a_Ruderal_habitats",
      "12b_Agricultural_habitats") %>%
      stringr::str_subset(paste0("^", as.character(Hab_Abb), "_")) %>%
      paste0("ModDT_", ., "_subset")
  }

  IASDT.R::SaveAs(
    InObj = DT_All, OutObj = OutObjName,
    OutPath = file.path(Path_Model, paste0(OutObjName, ".RData")))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Response - Y matrix ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Response - Y matrix")
  DT_y <- dplyr::select(DT_All, tidyselect::starts_with("Sp_")) %>%
    as.data.frame()
  IASDT.R::CatTime(paste0("   >>>   ", ncol(DT_y), " species"))

  IASDT.R::CatTime("   >>>   Save species summary")
  SpSummary <- file.path(Path_PA, "Sp_PA_Summary_DF.RData")
  if (magrittr::not(file.exists(SpSummary))) {
    stop(paste0(SpSummary, " file does not exist"))
  }
  SpSummary <- IASDT.R::LoadAs(SpSummary) %>%
    dplyr::mutate(
      IAS_ID = stringr::str_pad(IAS_ID, pad = "0", width = 4),
      IAS_ID = paste0("Sp_", IAS_ID)) %>%
    dplyr::filter(IAS_ID %in% names(DT_y))

  save(SpSummary, file = file.path(Path_Model, "SpSummary.RData"))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Xformula -----
  # # |||||||||||||||||||||||||||||||||||

  # The formula object becomes too large (up to > 2GB!) if created within a
  # function. Setting the environment of the formula as an empty environment
  # release this unnecessary size. https://stackoverflow.com/questions/66241212
  IASDT.R::CatTime("Xformula")
  Form_x <- stringr::str_c(XVars, collapse = " + ") %>%
    stringr::str_c("~ ", .) %>%
    stats::as.formula(env = baseenv())
  DT_x <- dplyr::select(DT_All, tidyselect::all_of(XVars)) %>%
    as.data.frame()

  IASDT.R::CatTime(
    paste0(
      "   >>>   Models will be fitted using ", length(XVars),
      " predictors: \n   >>>   ", paste0(XVars, collapse = " + ")))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Save cross-validation data -----
  # # |||||||||||||||||||||||||||||||||||

  DT_CV <- DT_All %>%
    dplyr::select(
      "CellNum", "CellCode", "Country", tidyselect::starts_with("CV"))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Phylogenetic tree data -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Phylogenetic tree data")

  if (PhyloTree) {
    # Taxonomy as a proxy for phylogeny
    plant.tree <- readr::read_tsv(TaxaInfoFile, show_col_types = FALSE) %>%
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

  Tree <- c("Tree", "NoTree")[c(PhyloTree, NoPhyloTree)]

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Preparing working on parallel -----
  # # |||||||||||||||||||||||||||||||||||

  if (NParallel > 1) {
    IASDT.R::CatTime(paste0(
      "Preparing working on parallel (", NParallel, " cores)"))
    c1 <- snow::makeSOCKcluster(NParallel)
    on.exit(invisible(try(snow::stopCluster(c1), silent = TRUE)), add = TRUE)
    future::plan(future::cluster, workers = c1, gc = TRUE)
  } else {
    IASDT.R::CatTime("Working sequentially")
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Spatial info / random effect ------
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Spatial info / random effect")
  studyDesign <- data.frame(sample = as.factor(seq_len(nrow(DT_x))))

  DT_xy <- as.matrix(dplyr::select(DT_All, x, y))
  rownames(DT_xy) <- studyDesign$sample
  rm(DT_All)

  # Prepare GPP knots
  IASDT.R::CatTime("   >>>   Prepare GPP knots")

  if (NParallel > 1) {
    invisible(snow::clusterEvalQ(
      cl = c1, IASDT.R::LoadPackages(dplyr, sf, Hmsc, jsonify, magrittr)))
    snow::clusterExport(
      cl = c1, list = c("DT_xy", "GPP_Dists"), envir = environment())

    GPP_Knots <- snow::parLapply(
      cl = c1, x = GPP_Dists * 1000,
      fun = function(x) {
        IASDT.R::PrepKnots(
          Coords = DT_xy, MinDist = x, MinLF = rLMinLF, MaxLF = rLMaxLF)
      }) %>%
      stats::setNames(paste0("GPP_", GPP_Dists))
  } else {
    GPP_Knots <- purrr::map(
      .x = GPP_Dists * 1000,
      .f = ~IASDT.R::PrepKnots(
        Coords = DT_xy, MinDist = .x, MinLF = rLMinLF, MaxLF = rLMaxLF)) %>%
      stats::setNames(paste0("GPP_", GPP_Dists))
  }

  ## Plotting knot location ----
  if (GPP_Plot) {
    IASDT.R::CatTime("   >>>   Plotting GPP knots")

    KnotExt <- rbind(
      as.matrix(purrr::map_dfr(GPP_Knots, ~.x$sKnot)), as.matrix(DT_xy)) %>%
      as.data.frame() %>%
      stats::setNames(c("x", "y")) %>%
      sf::st_as_sf(coords = c("x", "y"), crs = 3035) %>%
      sf::st_buffer(10000) %>%
      sf::st_bbox()
    AspectRatio <- (KnotExt[3] - KnotExt[1]) / (KnotExt[4] - KnotExt[2])

    GridR <- IASDT.R::LoadAs(Path_GridR) %>%
      terra::unwrap()
    GridR <- sf::st_as_sf(
      x = data.frame(DT_xy), coords = c("x", "y"), crs = 3035) %>%
      terra::rasterize(GridR) %>%
      terra::trim() %>%
      stats::setNames("GridR")

    EU_Bound <- IASDT.R::LoadAs(EU_Bound) %>%
      magrittr::extract2("Bound_sf_Eur_s") %>%
      magrittr::extract2("L_01") %>%
      sf::st_crop(KnotExt) %>%
      suppressWarnings()

    Knots_Plots <- purrr::map(
      .x = GPP_Dists,
      .f = ~{

        Knot_sf <- GPP_Knots[[paste0("GPP_", .x)]]$sKnot %>%
          sf::st_as_sf(coords = c("Var1", "Var2"), crs = 3035)
        NKnots <- nrow(GPP_Knots[[paste0("GPP_", .x)]]$sKnot) %>%
          formatC(format = "d", big.mark = ",")

        Plot <- ggplot2::ggplot() +
          tidyterra::geom_spatraster(data = GridR) +
          ggplot2::geom_sf(
            data = EU_Bound, fill = "transparent", colour = "darkgrey",
            linewidth = 1.5) +
          ggplot2::geom_sf(
            data = Knot_sf, colour = "black", shape = 19,
            stroke = 1.75, size = dplyr::if_else(.x < 50, 1.5, 3)) +
          ggplot2::scale_x_continuous(
            limits = sf::st_bbox(EU_Bound)[c(1, 3)], expand = c(0, 0)) +
          ggplot2::scale_y_continuous(
            limits = sf::st_bbox(EU_Bound)[c(2, 4)], expand = c(0, 0)) +
          ggplot2::ggtitle(
            paste0(
              "<span style='font-size: 35pt;'>GPP knots</span><br>",
              "<span style='font-size: 30pt;'>  Minimum distance between ",
              "knots and between knots and grid ",
              " cells is ", .x, " km  &mdash; ", NKnots, " knots</span>")) +
          ggplot2::scale_fill_continuous(na.value = "transparent") +
          ggplot2::theme_void() +
          ggplot2::theme(
            plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
            plot.title = ggtext::element_markdown(),
            legend.position = "none")

        return(Plot)
      })

    ggplot2::ggsave(
      filename = file.path(Path_Model, "knot_Locations.pdf"),
      plot = gridExtra::marrangeGrob(
        Knots_Plots, nrow = 1, ncol = 1, top = NULL),
      width = 25 * AspectRatio, height = 25, unit = "in")
  }

  if (GPP_Save) {
    IASDT.R::CatTime("   >>>   Saving GPP knots data")
    save(GPP_Knots, file = file.path(Path_Model, "GPP_Knots.RData"))
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Define the models -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Define the models")

  ModelVariants <- tidyr::expand_grid(M_thin = thin, M_samples = samples) %>%
    dplyr::mutate(M_transient = M_thin * transientFactor)

  Model_Info <- tibble::tibble(rL = GPP_Dists, rL2 = GPP_Knots) %>%
    tidyr::expand_grid(Tree = Tree) %>%
    dplyr::mutate(
      M_Name_init = paste0("GPP", rL, "_", Tree),
      M_Init_Path = purrr::map2_chr(
        .x = M_Name_init, .y = rL2,
        .f = ~{

          PathOut <- file.path(Path_Model, paste0("InitMod_", .x, ".RData"))

          if (file.exists(PathOut)) {
            InitModel <- IASDT.R::LoadAs(PathOut)
          } else {
            if (stringr::str_detect(.x, "_Tree$")) {
              Tree <- plant.tree
            } else {
              Tree <- NULL
            }

            InitModel <- Hmsc::Hmsc(
              Y = DT_y, XFormula = Form_x, XData = DT_x, distr = "probit",
              studyDesign = studyDesign, ranLevels = list("sample" = .y),
              phyloTree = Tree)

            IASDT.R::SaveAs(
              InObj = InitModel, OutObj = paste0("InitMod_", .x),
              OutPath = PathOut)
          }
          return(PathOut)
        })) %>%
    tidyr::expand_grid(ModelVariants) %>%
    dplyr::mutate(
      rL2 = NULL,
      M_HPC = purrr::pmap(
        .l = list(M_Name_init, M_thin, M_samples),
        .f = function(M_Name_init, M_thin, M_samples) {

          M_Name_Fit <- paste0(M_Name_init, "_samp", M_samples, "_th", M_thin)

          M4HPC_Path <- file.path(
            Path_Model, "InitMod_HPC",
            paste0("InitMod_HPC_", M_Name_Fit, ".rds"))

          return(list(M_Name_Fit = M_Name_Fit, M4HPC_Path = M4HPC_Path))

        })) %>%
    tidyr::unnest_wider("M_HPC")

  # # |||||||||||||||||||||||||||||||||||
  # # Prepare and save unfitted models -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Save unfitted models")

  if (OverwriteInitMod) {
    IASDT.R::CatTime(paste0("   >>>   Processing all model variants"))
  } else {
    NMod2Export <- Model_Info %>%
      dplyr::filter(magrittr::not(file.exists(M4HPC_Path))) %>%
      nrow()
    if (NMod2Export == 0) {
      IASDT.R::CatTime(
        "   >>>   All model variants were already available as RDS files")
    } else {
      IASDT.R::CatTime(
        paste0("   >>>   ", NMod2Export,
               " model variants need to be exported as RDS files"))
    }
  }

  InitFitFun <- function(ID) {
    CurrPath <- Model_Info$M4HPC_Path[ID]
    if (OverwriteInitMod && file.exists(CurrPath)) {
      file.remove(CurrPath)
    }
    Try <- 0
    while (Try < 6 && magrittr::not(file.exists(CurrPath))) {
      Try <- Try + 1
      Model <- Hmsc::sampleMcmc(
        hM = IASDT.R::LoadAs(Model_Info$M_Init_Path[ID]),
        samples = Model_Info$M_samples[ID],
        thin = Model_Info$M_thin[ID],
        transient = Model_Info$M_transient[ID],
        nChains = nChains,
        verbose = verbose, engine = "HPC")
      if (ToJSON) {
        Model <- jsonify::to_json(Model)
      }
      saveRDS(Model, file = CurrPath)
      if (file.exists(CurrPath)) {
        break
      }
    }
    return(invisible(NULL))
  }


  if (NParallel > 1) {
    invisible(snow::clusterEvalQ(
      cl = c1, IASDT.R::LoadPackages(dplyr, Hmsc, jsonify, magrittr)))
    snow::clusterExport(
      cl = c1, envir = environment(),
      list = c("InitFitFun", "Model_Info", "OverwriteInitMod",
               "verbose", "nChains", "ToJSON"))

    Model_Process <- snow::parLapply(
      cl = c1, x = seq_len(nrow(Model_Info)), fun = InitFitFun)

  } else {
    Model_Process <- purrr::map(
      .x = seq_len(nrow(Model_Info)), .f = InitFitFun)
  }

  Failed2Export <- Model_Info %>%
    dplyr::filter(magrittr::not(file.exists(M4HPC_Path)))

  if (nrow(Failed2Export) == 0) {
    IASDT.R::CatTime(
      "   >>>   All model variants were exported as RDS files")
  } else {
    IASDT.R::CatTime(
      paste0(
        "   >>>   ", nrow(Failed2Export),
        " model variants failed to be exported to rds files after 5 tries."))
    save(Failed2Export, file = file.path(Path_Model, "Failed2Export.RData"))
    readr::write_tsv(
      x = Failed2Export, file = file.path(Path_Model, "Failed2Export.txt"))
  }

  # # # |||||||||||||||||||||||||||||||||||
  # # # Stopping cluster -----
  # # # |||||||||||||||||||||||||||||||||||

  if (NParallel > 1) {
    IASDT.R::CatTime("Stopping cluster")
    snow::stopCluster(c1)
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Prepare Hmsc-HPC fitting commands -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Prepare Hmsc-HPC fitting commands")

  Model_Info <- dplyr::mutate(Model_Info, Chain = list(seq_len(nChains))) %>%
    tidyr::unnest_longer("Chain") %>%
    dplyr::arrange(M_Name_Fit) %>%

    dplyr::mutate(
      M_Chain = purrr::pmap(
        .l = list(M_Name_Fit, M4HPC_Path, Chain,
                  M_transient, M_samples, M_thin),
        .f = function(M_Name_Fit, M4HPC_Path, Chain,
                      M_transient, M_samples, M_thin) {

          # Input model
          M4HPC_Path2 <- file.path(
            Path_Model, "InitMod_HPC", basename(M4HPC_Path))

          # Path for posterior sampling
          Post_Path <- file.path(
            Path_Model, "Model_Fitting_HPC",
            paste0(M_Name_Fit, "_Chain", Chain, "_post.rds"))

          # Path for progress
          Path_ModProg <- file.path(
            Path_Model, "Model_Fitting_HPC",
            paste0(M_Name_Fit, "_Chain", Chain, "_Progress.txt"))

          Post_Missing <- magrittr::not(file.exists(Post_Path))

          Command_HPC <- paste0(
            "export TF_CPP_MIN_LOG_LEVEL=3; /usr/bin/time -v ", Path_Python,
            " -m hmsc.run_gibbs_sampler",
            " --input ", shQuote(M4HPC_Path2),
            " --output ", shQuote(Post_Path),
            " --samples ", M_samples,
            " --transient ", M_transient,
            " --thin ", M_thin,
            " --verbose ", verbose,
            " --chain ", (Chain - 1),
            " >& ", shQuote(Path_ModProg))

          Command_WS <- paste0(
            "-m hmsc.run_gibbs_sampler",
            " --input ", shQuote(M4HPC_Path2),
            " --output ", shQuote(Post_Path),
            " --samples ", M_samples,
            " --transient ", M_transient,
            " --thin ", M_thin,
            " --verbose ", verbose,
            " --chain ", (Chain - 1))

          list(
            M4HPC_Path_LUMI = M4HPC_Path2,
            Post_Path = Post_Path, Post_Missing = Post_Missing,
            Path_ModProg = Path_ModProg, Command_HPC = Command_HPC,
            Command_WS = Command_WS) %>%
            return()
        })) %>%
    tidyr::unnest_wider("M_Chain")

  # # |||||||||||||||||||||||||||||||||||
  # # Skip fitted models -----
  # # |||||||||||||||||||||||||||||||||||

  if (SkipFitted) {
    IASDT.R::CatTime("Skip fitted models")
    Models2Fit_HPC <- Model_Info %>%
      dplyr::filter(Post_Missing) %>%
      dplyr::pull(Command_HPC) %>%
      unlist()
    Models2Fit_WS <- Model_Info %>%
      dplyr::filter(Post_Missing) %>%
      dplyr::pull(Command_WS) %>%
      unlist()
  } else {
    Models2Fit_HPC <- Model_Info %>%
      dplyr::pull(Command_HPC) %>%
      unlist()
    Models2Fit_WS <- Model_Info %>%
      dplyr::pull(Command_WS) %>%
      unlist()
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Save commands in a text file -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Save commands in a text file")

  NJobs <- length(Models2Fit_HPC)

  if (NJobs > MaxJobCounts) {
    NSplits <- ceiling((NJobs / MaxJobCounts))
    IDs <- IASDT.R::SplitVector(Vector = seq_len(NJobs), NSplit = NSplits)
  } else {
    NSplits <- 1
    IDs <- list(seq_len(NJobs))
  }

  # Save all fitting commands to single file
  f <- file(
    description = file.path(Path_Model, "Commands_All.txt"),
    open = "wb")
  on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
  cat(Models2Fit_HPC, sep = "\n", append = FALSE, file = f)
  close(f)

  # Save fitting commands into separate files, when necessary
  purrr::walk(
    .x = seq_len(NSplits),
    .f = function(x) {
      CurrIDs <- IDs[[x]]
      if (NSplits > 1) {
        CommandFile <- file.path(Path_Model, paste0("Commands2Fit_", x, ".txt"))
      } else {
        CommandFile <- file.path(Path_Model, "Commands2Fit.txt")
      }

      # create connection to SLURM file This is better than using sink to have a
      # platform independent file (here, to maintain a linux-like new line
      # ending)
      f <- file(CommandFile, open = "wb")
      on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
      cat(Models2Fit_HPC[CurrIDs], sep = "\n", append = FALSE, file = f)
      close(f)
    })

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Save data to disk -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Save data to disk")

  Model_Info <- Model_Info %>%
    tidyr::nest(
      Post_Path = Post_Path, Path_ModProg = Path_ModProg,
      Chain = Chain, Command_HPC = Command_HPC, Command_WS = Command_WS,
      Post_Missing = Post_Missing) %>%
    dplyr::mutate(
      Post_Path = purrr::map2(Post_Path, Chain, IASDT.R::SetChainName),
      Chain = purrr::map2(Chain, Chain, IASDT.R::SetChainName),
      Command_HPC = purrr::map2(Command_HPC, Chain, IASDT.R::SetChainName),
      Command_WS = purrr::map2(Command_WS, Chain, IASDT.R::SetChainName),
      Path_ModProg = purrr::map2(Path_ModProg, Chain, IASDT.R::SetChainName),
      Post_Missing = purrr::map2(Post_Missing, Chain, IASDT.R::SetChainName),
      Post_Aligned = NA)

  save(Model_Info, file = Path_ModelDT)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Prepare SLURM file ------
  # # |||||||||||||||||||||||||||||||||||

  if (PrepSLURM) {
    IASDT.R::CatTime("Preparing SLURM file")
    if (is.null(JobName)) {
      JobName <- stringr::str_remove_all(
        basename(Path_Model), paste0("_", HabVal))
    }
    IASDT.R::Mod_SLURM(
      Path_Model = Path_Model, JobName = JobName, MemPerCpu = MemPerCpu,
      Time = Time, EnvFile = EnvFile, FromHPC = FromHPC,
      Path_Hmsc = Path_Hmsc, ...)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Save small datasets prepared in the function ------
  # # |||||||||||||||||||||||||||||||||||

  DT_Split <- list(
    DT_y = DT_y, Form_x = Form_x, DT_x = DT_x, XVars = XVars, DT_CV = DT_CV,
    PhyloTree = PhyloTree, plant.tree = plant.tree, Tree = Tree,
    studyDesign = studyDesign, DT_xy = DT_xy, GPP_Knots = GPP_Knots)

  OutObjName <- stringr::str_replace(OutObjName, "_subset", "_split")

  IASDT.R::SaveAs(
    InObj = DT_Split, OutObj = OutObjName,
    OutPath = file.path(Path_Model, paste0(OutObjName, ".RData")))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatDiff(
    InitTime = .StartTime, ChunkText = "Function summary", CatInfo = TRUE)

  return(invisible(NULL))
}
