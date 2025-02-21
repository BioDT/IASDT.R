## |------------------------------------------------------------------------| #
# Mod_Prep4HPC ----
## |------------------------------------------------------------------------| #

#' Prepare initial models for model fitting with Hmsc-HPC
#'
#' The **`Mod_Prep4HPC`** function prepares input data and initializes models
#' for fitting with [Hmsc-HPC](https://doi.org/10.1371/journal.pcbi.1011914). It
#' performs multiple tasks, including data preparation, defining spatial block
#' cross-validation folds, generating Gaussian Predictive Process (GPP) knots
#' ([Tikhonov et al.](https://doi.org/10.1002/ecy.2929)), initializing models,
#' and creating HPC execution commands. The function supports parallel
#' processing and offers the option to include or exclude phylogenetic tree
#' data.<br/><br/> The internal **`Mod_PrepData`** function is used to prepare
#' habitat-specific data for Hmsc models. This function processes environmental
#' and species presence data, reads environment variables from a file, verifies
#' paths, loads and filters species data based on habitat type and minimum
#' presence grid cells per species, and merges various environmental layers
#' (e.g., CHELSA Bioclimatic variables, habitat coverage, road and railway
#' intensity, sampling efforts) into a single dataset. Processed data is saved
#' to disk as an `*.RData` file.
#' @param Path_Model Character. Directory path where all output files, including
#'   the models to be fitted, will be saved.
#' @param GPP Logical. Whether to fit spatial random effect using Gaussian
#'   Predictive Process. Defaults to `TRUE`. If `FALSE`, non-spatial models will
#'   be fitted.
#' @param GPP_Dists Integer. Spacing (in kilometers) between GPP knots, as well
#'   as the minimum allowable distance between a knot and the nearest sampling
#'   point. The knots are generated using the [Mod_PrepKnots] function, and this
#'   value is used for both `knotDist` and `minKnotDist` in
#'   [Hmsc::constructKnots].
#' @param GPP_Save Logical. Whether `RData` file. Default: `TRUE`.
#' @param GPP_Plot Logical. Whether to plot the coordinates of the sampling
#'   units and the knots in a pdf file. Default: `TRUE`.
#' @param BioVars Character vector. Variables from CHELSA (bioclimatic variables
#'   (bio1–bio19) and additional predictors (e.g., Net Primary Productivity,
#'   npp)) to be used in the model. By default, six ecologically relevant and
#'   minimally correlated variables are selected: c("bio3", "bio4", "bio11",
#'   "bio18", "bio19", "npp").
#' @param QuadraticVars Character vector for variables for which quadratic terms
#'   are used. Defaults to all variables of the `BioVars`. If `QuadraticVars` is
#'   `NULL`, no quadratic terms will be used.
#' @param EffortsAsPredictor Logical. Whether to include the (log<sub>10</sub>)
#'   sampling efforts as predictor to the model. Default: `TRUE`.
#' @param RoadRailAsPredictor Logical. Whether to include the (log<sub>10</sub>)
#'   sum of road and railway intensity as predictor to the model. Default:
#'   `TRUE`.
#' @param HabAsPredictor Logical. Whether to include the (log<sub>10</sub>)
#'   percentage coverage of respective habitat type per grid cell as predictor
#'   to the model. Default: `TRUE`. Only valid if `Hab_Abb` not equals to `0`.
#' @param RiversAsPredictor Logical. Whether to include the (log<sub>10</sub>)
#'   total length of rivers per grid cell as predictor to the model. Default:
#'   `TRUE`. See [River_Length] for more details.
#' @param NspPerGrid Integer. Minimum number of species required for a grid cell
#'   to be included in the analysis. This filtering occurs after applying
#'   `MinEffortsSp` (sampling effort thresholds), `PresPerSpecies` (minimum
#'   species presence thresholds), and `ExcludeZeroHabitat` (exclude 0% habitat
#'   coverage). Default (0): Includes all grid cells. Positive value (>0):
#'   Includes only grid cells where at least `NspPerGrid` species are present.
#' @param PhyloTree,NoPhyloTree Logical. Whether to fit models with (PhyloTree)
#'   or without (NoPhyloTree) phylogenetic trees. Defaults are `PhyloTree =
#'   TRUE` and `NoPhyloTree = FALSE`, meaning only models with phylogenetic
#'   trees are fitted by default. At least one of `PhyloTree` and `NoPhyloTree`
#'   should be `TRUE`.
#' @param OverwriteRDS Logical. Whether to overwrite previously exported RDS
#'   files for initial models. Default: `TRUE`.
#' @param NCores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 8.
#' @param NChains Integer. Number of model chains. Default: 4.
#' @param thin Integer vector. Thinning value(s) in MCMC sampling. If more than one
#'   value is provided, a separate model will be fitted at each value of
#'   thinning.
#' @param samples Integer vector. Value(s) for the number of MCMC samples. If more than
#'   one value is provided, a separate model will be fitted at each value of
#'   number of samples. Defaults to 1000.
#' @param transientFactor Integer. Transient multiplication factor. The value of
#'   `transient` will equal the multiplication of `transientFactor` and `thin`.
#'   Default: 500.
#' @param verbose Integer. Interval at which MCMC sampling progress is reported.
#'   Default: `200`.
#' @param SkipFitted Logical. Whether to skip already fitted models. Default:
#'   `TRUE`.
#' @param NumArrayJobs Integer. Number of jobs per SLURM script file. In LUMI
#'   HPC, there is a limit of 210 submitted jobs per user for the `small-g`
#'   partition. This argument is used to split the jobs into multiple SLURM
#'   scripts if needed. Default: 210. See [LUMI
#'   documentation](https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/partitions)
#'   for more details.
#' @param ModelCountry Character. Country or countries to filter observations
#'   by. Default: `NULL`, which means prepare data for the whole Europe.
#' @param PrepSLURM Logical. Whether to prepare SLURM command files. If `TRUE`
#'   (default), the SLURM commands will be saved to disk using the
#'   [Mod_SLURM] function.
#' @param MemPerCpu Character. Memory per CPU for the SLURM job. This value will
#'   be assigned to the `#SBATCH --mem-per-cpu=` SLURM argument. Example: "32G"
#'   to request 32 gigabyte. Only effective if `PrepSLURM = TRUE`.
#' @param Time Character. Requested time for each job in the SLURM bash arrays.
#'   Example: "01:00:00" to request an hour. Only effective if `PrepSLURM =
#'   TRUE`.
#' @param JobName Character. Name of the submitted job(s) for SLURM. If `NULL`
#'   (Default), the job name will be prepared based on the folder path and the
#'   `Hab_Abb` value. Only effective if `PrepSLURM = TRUE`.
#' @param Path_Hmsc  Character. Directory path to `Hmsc-HPC` extension
#'   installation. This will be provided as the `Path_Hmsc` argument of the
#'   [Mod_SLURM] function.
#' @param ToJSON Logical. Whether to convert unfitted models to JSON before
#'   saving to RDS file. Default: `FALSE`.
#' @param CheckPython Logical. Whether to check if the Python executable exists.
#'   Only valid if FromHPC = `FALSE`.
#' @param Precision Integer (either 32 or 64). Defines the floating-point
#'   precision mode for `Hmsc-HPC` sampling (--fp 32 or --fp 64). The default is
#'   64, which is the default precision in `Hmsc-HPC`.
#' @param Hab_Abb Character. Abbreviation for the habitat type (based on
#'   [SynHab](https://www.preslia.cz/article/pdf?id=11548)) for which to prepare
#'   data. Valid values are `0`, `1`, `2`, `3`, `4a`, `4b`, `10`, `12a`, `12b`.
#'   If `Hab_Abb` = `0`, data is prepared irrespective of the habitat type. For
#'   more details, see [Pysek et
#'   al.](https://doi.org/10.23855/preslia.2022.447).
#' @param MinEffortsSp Integer. Minimum number of vascular plant species per
#'   grid cell (from GBIF data) required for inclusion in the models. This is to
#'   exclude grid cells with very little sampling efforts. Defaults to `100`.
#' @param ExcludeCult Logical. Whether to exclude countries with cultivated or
#'   casual observations per species. Defaults to `TRUE`.
#' @param ExcludeZeroHabitat Logical. Whether to exclude grid cells with zero
#'   percentage habitat coverage. Defaults to `TRUE`.
#' @param PresPerSpecies Integer. The minimum number of presence grid cells for
#'   a species to be included in the analysis. The number of presence grid cells
#'   per species is calculated after discarding grid cells with low sampling
#'   efforts (`MinEffortsSp`) and zero percentage habitat coverage
#'   `ExcludeZeroHabitat`. Defaults to `80`.
#' @param CV_NGrids Integer. For `CV_Dist` cross-validation strategy (see
#'   [Mod_GetCV]), this argument determines the size of the blocks (how many
#'   grid cells in both directions).
#' @param CV_NR,CV_NC Integer. Number of rows and columns used in the `CV_Large`
#'   cross-validation strategy  (see [Mod_GetCV]), in which the study area is
#'   divided into large blocks given the provided `CV_NR` and `CV_NC` values.
#'   Both default to 2 which means to split the study area into four large
#'   blocks at the median latitude and longitude.
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param VerboseProgress Logical. Whether to print a message upon successful
#'   saving of files. Defaults to `FALSE`.
#' @param FromHPC Logical. Whether the processing is being done on an
#'   High-Performance Computing (HPC) environment, to adjust file paths
#'   accordingly. Default: `TRUE`.
#' @param ... Additional parameters provided to the [Mod_SLURM] function.
#' @export
#' @inheritParams Mod_PrepKnots
#' @inheritParams Mod_GetCV
#' @name Mod_inputs
#' @rdname Mod_inputs
#' @order 1
#' @importFrom rlang .data
#' @author Ahmed El-Gabbas
#' @details The current models are fitted for 8 habitat types see [Pysek et
#'   al.](https://doi.org/10.23855/preslia.2022.447):
#' - **1. Forests** -- closed vegetation dominated by deciduous or evergreen
#'   trees
#' - **2. Open forests** -- woodlands with canopy openings created by
#'   environmental stress or disturbance, including forest edges
#' - **3. Scrub** -- shrublands maintained by environmental stress (aridity) or
#'   disturbance
#' - **4a. Natural grasslands** -- grasslands maintained by climate (aridity,
#'   unevenly distributed precipitation), herbivores or environmental stress
#'   (aridity, instability or toxicity of substrate)
#' - **4b. Human-maintained grasslands** -- grasslands dependent on regular
#'   human-induced management (mowing, grazing by livestock, artificial burning)
#' - **10. Wetland** -- sites with the permanent or seasonal influence of
#'   moisture, ranging from oligotrophic to eutrophic
#' - **12a. Ruderal habitats** -- anthropogenically disturbed or eutrophicated
#'   sites, where the anthropogenic disturbance or fertilization is typically a
#'   side-product and not the aim of the management
#' - **12b. Agricultural habitats** -- synanthropic habitats directly
#'   associated with growing of agricultural products, thus dependent on
#'   specific type of management (ploughing, fertilization)
#'
#'   <br/>The following habitat types are excluded from the analysis:
#' - **5. Sandy** -- dunes and other habitats on unstable sandy substrate,
#'   stressed by low nutrients, drought and disturbed by sand movement
#' - **6. Rocky** -- cliffs and rock outcrops with very shallow or no soil
#' - **7. Dryland** -- habitats in which drought stress limits vegetation
#'   development
#' - **8. Saline** -- habitats stressed by high soil salinity
#' - **9. Riparian** -- a mosaic of wetlands, grasslands, tall-forb stands,
#'   scrub and open forests in stream corridors
#' - **11. Aquatic** -- water bodies and streams with submerged and floating
#'   plant species

Mod_Prep4HPC <- function(
    Hab_Abb = NULL, Path_Model = NULL,
    MinEffortsSp = 100L, PresPerSpecies = 80L, EnvFile = ".env",
    GPP = TRUE, GPP_Dists = NULL, GPP_Save = TRUE,
    GPP_Plot = TRUE, MinLF = NULL, MaxLF = NULL,
    Alphapw = list(Prior = NULL, Min = 20, Max = 1200, Samples = 200),
    BioVars = c("bio3", "bio4", "bio11", "bio18", "bio19", "npp"),
    QuadraticVars = BioVars,
    EffortsAsPredictor = TRUE, RoadRailAsPredictor = TRUE,
    HabAsPredictor = TRUE, RiversAsPredictor = TRUE,
    NspPerGrid = 0L, ExcludeCult = TRUE, ExcludeZeroHabitat = TRUE,
    CV_NFolds = 4L, CV_NGrids = 20L, CV_NR = 2L, CV_NC = 2L, CV_Plot = TRUE,
    CV_SAC = FALSE,
    PhyloTree = TRUE, NoPhyloTree = FALSE,
    OverwriteRDS = TRUE, NCores = 8L, NChains = 4L,
    thin = NULL, samples = 1000L, transientFactor = 500L, verbose = 200L,
    SkipFitted = TRUE, NumArrayJobs = 210L, ModelCountry = NULL,
    VerboseProgress = TRUE, FromHPC = TRUE, PrepSLURM = TRUE, MemPerCpu = NULL,
    Time = NULL, JobName = NULL, Path_Hmsc = NULL,
    CheckPython = FALSE, ToJSON = FALSE, Precision = 64L, ...) {

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Initial checking -----
  # # |||||||||||||||||||||||||||||||||||

  .StartTime <- lubridate::now(tzone = "CET")

  CheckNULL <- c(
    "Path_Model", "PresPerSpecies", "thin", "samples",
    "MemPerCpu", "Path_Hmsc", "Hab_Abb")
  IsNull <- purrr::map_lgl(CheckNULL, ~ is.null(get(.x)))

  if (any(IsNull)) {
    stop(
      paste0(
        paste0("`", CheckNULL[which(IsNull)], "`", collapse = ", "),
        " can not be empty"),
      call. = FALSE)
  }

  if (!(Precision %in% c(32, 64))) {
    stop(
      paste0(
        "Precision should be either of 32 or 64, not ", Precision),
      call. = FALSE)
  }

  Hab_Abb <- as.character(Hab_Abb)


  if (!all(is.numeric(samples)) || any(samples <= 0)) {
    stop("samples should be numeric and greater than zero", call. = FALSE)
  }

  if (!all(is.numeric(thin)) || any(thin <= 0)) {
    stop("thin should be numeric and greater than zero", call. = FALSE)
  }

  if (!all(is.numeric(PresPerSpecies)) || PresPerSpecies <= 0) {
    stop(
      "`PresPerSpecies` should be numeric and greater than zero",
      call. = FALSE)
  }

  if (!all(is.numeric(MinEffortsSp)) || MinEffortsSp <= 0) {
    stop(
      "`MinEffortsSp` should be numeric and greater than zero",
      call. = FALSE)
  }

  if (!is.numeric(NspPerGrid) || NspPerGrid < 0) {
    stop(
      "`NspPerGrid` has to be integer >= 0", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  NCells <- Sp <- IAS_ID <- x <- y <- Country <- M_thin <- rL <-
    M_Name_init <- rL2 <- M_samples <- M4HPC_Path <- M_transient <-
    M_Name_Fit <- Chain <- Post_Missing <- Command_HPC <-
    Command_WS <- Post_Path <- Path_ModProg <- TaxaInfoFile <-
    Path_Grid <- EU_Bound <- Path_PA <- NAME_ENGL <- NSp <- NULL

  if (isFALSE(VerboseProgress)) {
    sink(file = nullfile())
    on.exit(sink(), add = TRUE)
  }

  Path_Python <- IASDT.R::Path(Path_Hmsc, "Scripts", "python.exe")

  IASDT.R::InfoChunk("Preparing data for Hmsc-HPC models", Char = "=")

  # # |||||||||||||||||||||||||||||||||||
  # Load/check environment variables -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Load/check environment variables")

  if (!file.exists(EnvFile)) {
    stop(paste0(
      "Path to environment variables: ", EnvFile, " was not found"),
      call. = FALSE)
  }

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "TaxaInfoFile", "DP_R_TaxaInfo", FALSE, TRUE,
      "EU_Bound", "DP_R_EUBound_sf", FALSE, TRUE,
      "Path_Grid", "DP_R_Grid", TRUE, FALSE,
      "Path_PA", "DP_R_PA", TRUE, FALSE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "TaxaInfoFile", "DP_R_TaxaInfo_Local", FALSE, TRUE,
      "EU_Bound", "DP_R_EUBound_sf_Local", FALSE, TRUE,
      "Path_Grid", "DP_R_Grid_Local", TRUE, FALSE,
      "Path_PA", "DP_R_PA_Local", TRUE, FALSE)

    # Check if Python executable exists
    if (CheckPython && !file.exists(Path_Python)) {
      stop(
        paste0("Python executable does not exist: ", Path_Python),
        call. = FALSE)
    }
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Check input arguments ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Checking input arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  CharArgs <- c("Hab_Abb", "Path_Model", "Path_Hmsc", "EnvFile", "BioVars")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = CharArgs, Type = "character")

  LogicArgs <- c(
    "GPP_Save", "GPP_Plot", "PhyloTree", "NoPhyloTree", "ExcludeZeroHabitat",
    "SkipFitted", "VerboseProgress", "ToJSON", "CV_SAC", "CheckPython",
    "OverwriteRDS", "PrepSLURM", "ExcludeCult", "GPP",
    "EffortsAsPredictor", "RoadRailAsPredictor", "HabAsPredictor",
    "RiversAsPredictor")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = LogicArgs, Type = "logical")

  NumericArgs <- c(
    "NCores", "NChains", "thin", "samples", "verbose", "NumArrayJobs",
    "PresPerSpecies", "MinEffortsSp", "transientFactor", "CV_NFolds",
    "CV_NGrids", "CV_NR", "CV_NC", "Precision")
  if (GPP) {
    NumericArgs <- c(NumericArgs, "GPP_Dists")
  }
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = NumericArgs, Type = "numeric")


  if (PrepSLURM) {
    IASDT.R::CheckArgs(
      AllArgs = AllArgs, Args = c("MemPerCpu", "Time"), Type = "character")
  }

  # Phylogenetic tree options
  if (isFALSE(PhyloTree) && isFALSE(NoPhyloTree)) {
    stop(
      "At least one of PhyloTree or NoPhyloTree has to be true", call. = FALSE)
  }

  NumArgsInvalid <- purrr::map_lgl(.x = NumericArgs, .f = ~all(get(.x) < 1))
  if (any(NumArgsInvalid)) {
    paste0(
      "The following parameter(s) can not be < 1\n  >>  ",
      paste0(NumericArgs[NumArgsInvalid], collapse = " | ")) %>%
      stop(call. = FALSE)
  }

  rm(AllArgs, CharArgs, LogicArgs, NumericArgs, envir = environment())

  # # ..................................................................... ###

  if (GPP) {
    if (is.null(GPP_Dists)) {
      stop("`GPP_Dists` can not be empty", call. = FALSE)
    }
    if (!all(is.numeric(GPP_Dists)) || any(GPP_Dists <= 0)) {
      stop("`GPP_Dists` should be numeric and greater than zero", call. = FALSE)
    }

    if (GPP_Plot) {
      Path_GridR <- IASDT.R::Path(Path_Grid, "Grid_10_Land_Crop.RData")
      if (!file.exists(Path_GridR)) {
        stop(
          paste0("Path for the Europe boundaries does not exist: ", Path_GridR),
          call. = FALSE)
      }
    }
  }

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # File paths - Creating missing paths ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("File paths - Creating missing paths")
  fs::dir_create(IASDT.R::Path(Path_Model, "InitMod4HPC"))
  fs::dir_create(IASDT.R::Path(Path_Model, "Model_Fitting_HPC"))
  # Also create directory for SLURM outputs
  fs::dir_create(IASDT.R::Path(Path_Model, "Model_Fitting_HPC", "JobsLog"))

  Path_ModelDT <- IASDT.R::Path(Path_Model, "Model_Info.RData")

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # # Prepare list of predictors -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Prepare list of predictors")

  # Check BioVars values
  if (isFALSE(all(BioVars %in% IASDT.R::CHELSA_Vars$Variable))) {
    WrongBio <- BioVars[which(!(BioVars %in% IASDT.R::CHELSA_Vars$Variable))]
    stop(
      paste0(
        "The following are invalid Bioclimatic variables: ",
        paste0(WrongBio, collapse = " | ")),
      call. = FALSE)
  }

  XVars <- BioVars

  if (EffortsAsPredictor) {
    XVars <- c(XVars, "EffortsLog")
  }
  if (RoadRailAsPredictor) {
    XVars <- c(XVars, "RoadRailLog")
  }
  if (RiversAsPredictor) {
    XVars <- c(XVars, "RiversLog")
  }

  if (Hab_Abb != "0" && HabAsPredictor) {
    XVars <- c(XVars, "HabLog")
  }

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # # Preparing input data -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Preparing input data")
  ValidHabAbbs <- c(0:3, "4a", "4b", 10, "12a", "12b")
  if (isFALSE(as.character(Hab_Abb) %in% ValidHabAbbs)) {
    stop(
      paste0(
        "Hab_Abb has to be one of the following:\n >> ",
        paste0(ValidHabAbbs, collapse = " | ")),
      call. = FALSE)
  }

  HabVal <- c(
    "0_All", "1_Forests", "2_Open_forests", "3_Scrub",
    "4a_Natural_grasslands", "4b_Human_maintained_grasslands",
    "10_Wetland", "12a_Ruderal_habitats", "12b_Agricultural_habitats") %>%
    stringr::str_subset(paste0("^", as.character(Hab_Abb), "_"))

  IASDT.R::InfoChunk("Preparing input data using IASDT.R::Mod_PrepData")

  DT_All <- IASDT.R::Mod_PrepData(
    Hab_Abb = Hab_Abb, MinEffortsSp = MinEffortsSp, ExcludeCult = ExcludeCult,
    ExcludeZeroHabitat = ExcludeZeroHabitat, PresPerSpecies = PresPerSpecies,
    EnvFile = EnvFile, Path_Model = Path_Model,
    VerboseProgress = VerboseProgress, FromHPC = FromHPC)

  IASDT.R::CatSep(Rep = 1, Extra1 = 1, Extra2 = 2)

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # # Subsetting study area -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Subsetting study area")

  if (!is.null(ModelCountry)) {

    ValidCountries <- ModelCountry %in% unique(DT_All$Country)

    if (!all(ValidCountries)) {
      stop(paste0(
        "The following are invalid country names: ",
        paste0(ModelCountry[!ValidCountries], collapse = " & ")),
        call. = FALSE)
    }

    IASDT.R::CatTime(
      paste0(
        "Subsetting data to: ", paste0(sort(ModelCountry), collapse = " & ")),
      Level = 1)

    Sample_ExclSp <- dplyr::filter(DT_All, Country %in% ModelCountry) %>%
      dplyr::summarise(
        dplyr::across(tidyselect::starts_with("Sp_"), sum)) %>%
      tidyr::pivot_longer(
        cols = tidyselect::everything(),
        names_to = "Sp", values_to = "NCells") %>%
      dplyr::filter(NCells < PresPerSpecies) %>%
      dplyr::pull(Sp)

    IASDT.R::CatTime(
      paste0(length(Sample_ExclSp), " species are excluded"), Level = 1)
    DT_All <- dplyr::filter(DT_All, Country %in% ModelCountry) %>%
      dplyr::select(-tidyselect::all_of(Sample_ExclSp))

    # # |||||||||||||||||||||||||||||||||||
    ## Plotting subsetted data -----
    # # |||||||||||||||||||||||||||||||||||

    IASDT.R::CatTime("Plotting subsetted data", Level = 1)

    NSpSubset <- DT_All %>%
      dplyr::mutate(
        NSp = rowSums(
          dplyr::select(., tidyselect::starts_with("Sp_")),
          na.rm = TRUE),
        NSp = as.integer(NSp)) %>%
      dplyr::select("x", "y", "NSp") %>%
      sf::st_as_sf(coords = c("x", "y"), crs = 3035) %>%
      terra::rasterize(
        y = terra::unwrap(IASDT.R::LoadAs(Path_GridR)),
        field = "NSp") %>%
      terra::trim()

    if (Hab_Abb == "0") {
      Hab_column <- NULL
    } else {
      Hab_column <- c(
        "Hab_1_Forests", "Hab_2_Open_forests", "Hab_3_Scrub",
        "Hab_4a_Natural_grasslands", "Hab_4b_Human_maintained_grasslands",
        "Hab_10_Wetland", "Hab_12a_Ruderal_habitats",
        "Hab_12b_Agricultural_habitats") %>%
        stringr::str_subset(paste0("_", as.character(Hab_Abb), "_")) %>%
        stringr::str_remove("Hab_")
    }

    EU_Bound_sub <- IASDT.R::LoadAs(EU_Bound) %>%
      magrittr::extract2("Bound_sf_Eur") %>%
      magrittr::extract2("L_03") %>%
      dplyr::filter(NAME_ENGL %in% ModelCountry)

    Limits <- as.vector(terra::ext(NSpSubset))
    # Relative JPEG height
    DimX <- Limits[2] - Limits[1]
    DimY <- Limits[4] - Limits[3]
    PlotHeight <- (DimY * 25) / (DimX * 0.95)
    NGrids <- format(nrow(DT_All), big.mark = ",")
    NSp <- length(stringr::str_subset(names(DT_All), "Sp_"))

    Caption <- paste0(sort(ModelCountry), collapse = "; ") %>%
      stringr::str_wrap(width = 110) %>%
      paste0("<strong>Selected countries</strong>: <br/>", .) %>%
      stringr::str_replace_all("\n", "<br/>")
    Subtitle <- paste0(
      NSp, " IAS within \u2265", PresPerSpecies,
      " presence grid cells in the selected country/countries (",
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
    grDevices::jpeg(
      filename = IASDT.R::Path(Path_Model, "NSpPerGrid_Sub.jpeg"),
      width = 25, height = PlotHeight, units = "cm", quality = 100, res = 600)
    print(NSpPerGrid_Sub)
    grDevices::dev.off()

    rm(Limits, NSpPerGrid_Sub, EU_Bound_sub, envir = environment())

  } else {
    IASDT.R::CatTime("No data subsetting was implemented", Level = 1)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Exclude grid cells with low number of presences -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Exclude grid cells with low number of presences")

  if (NspPerGrid == 0) {
    IASDT.R::CatTime(
      paste0(
        "All grid cells, irrespective of number of species ",
        "presence, will be considered"),
      Level = 1, Time = FALSE)
  } else {
    EmptyGridsID <- dplyr::select(DT_All, tidyselect::starts_with("Sp_")) %>%
      rowSums() %>%
      magrittr::is_less_than(as.integer(NspPerGrid)) %>%
      which() %>%
      magrittr::multiply_by(-1)

    if (length(EmptyGridsID) > 0) {
      IASDT.R::CatTime(
        paste0("Excluding grid cells with < ", NspPerGrid, " species presence"),
        Level = 1)
      DT_All <- dplyr::slice(DT_All, EmptyGridsID)
    }
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Cross-validation ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Prepare cross-validation folds")

  DT_All <- IASDT.R::Mod_GetCV(
    Data = DT_All, EnvFile = EnvFile, XVars = XVars,
    CV_NFolds = CV_NFolds, FromHPC = FromHPC, CV_NGrids = CV_NGrids,
    CV_NR = CV_NR, CV_NC = CV_NC, OutPath = Path_Model, CV_Plot = CV_Plot,
    CV_SAC = CV_SAC)

  # cross-validation data to be saved
  DT_CV <- DT_All %>%
    dplyr::select(
      "CellNum", "CellCode", "Country", tidyselect::starts_with("CV"))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Response - Y matrix ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Response - Y matrix")
  DT_y <- dplyr::select(DT_All, tidyselect::starts_with("Sp_")) %>%
    as.data.frame()
  IASDT.R::CatTime(paste0(ncol(DT_y), " species"), Level = 1)

  IASDT.R::CatTime("Save species summary", Level = 1)
  SpSummary <- IASDT.R::Path(Path_PA, "Sp_PA_Summary_DF.RData")
  if (!file.exists(SpSummary)) {
    stop(paste0(SpSummary, " file does not exist"), call. = FALSE)
  }
  SpSummary <- IASDT.R::LoadAs(SpSummary) %>%
    dplyr::arrange(IAS_ID) %>%
    dplyr::mutate(
      IAS_ID = stringr::str_pad(IAS_ID, pad = "0", width = 4),
      IAS_ID = paste0("Sp_", IAS_ID)) %>%
    dplyr::filter(IAS_ID %in% names(DT_y))

  save(SpSummary, file = IASDT.R::Path(Path_Model, "SpSummary.RData"))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Xformula -----
  # # |||||||||||||||||||||||||||||||||||

  # The formula object becomes too large (up to > 2GB!) if created within a
  # function. Setting the environment of the formula as an empty environment
  # release this unnecessary size. https://stackoverflow.com/questions/66241212

  IASDT.R::CatTime("Xformula")

  if (is.null(QuadraticVars)) {
    FormVars <- XVars
    IASDT.R::CatTime(
      paste0(
        "Models will be fitted using ", length(XVars), " predictors: ",
        paste0(XVars, collapse = " + ")), Level = 1)
  } else {
    OnlyLinear <- setdiff(XVars, QuadraticVars)
    FormVars <- c(
      OnlyLinear,
      paste0("stats::poly(", QuadraticVars, ", degree = 2, raw = TRUE)"))

    IASDT.R::CatTime("Models will be fitted using:", Level = 1)

    IASDT.R::CatTime(paste0(length(OnlyLinear), " linear effect: "), Level = 2)
    IASDT.R::CatTime(paste0(OnlyLinear, collapse = " + "), Level = 3)

    IASDT.R::CatTime(
      paste0(length(QuadraticVars), " linear and quadratic effects: "),
      Level = 2)
    IASDT.R::CatTime(paste0(QuadraticVars, collapse = " + "), Level = 3)

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

  IASDT.R::CatTime(
    paste0("Models will be fitted using ", paste0(Tree, collapse = " & ")),
    Level = 1)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Spatial info / random effect ------
  # # |||||||||||||||||||||||||||||||||||

  if (GPP) {

    IASDT.R::CatTime("Spatial info / random effect")
    studyDesign <- data.frame(sample = as.factor(seq_len(nrow(DT_x))))

    DT_xy <- as.matrix(dplyr::select(DT_All, x, y))
    rownames(DT_xy) <- studyDesign$sample

    # Prepare GPP knots
    IASDT.R::CatTime("Preparing GPP knots", Level = 1)

    NCores_GPP <- length(GPP_Dists)

    if (NCores_GPP > 1) {

      IASDT.R::CatTime(
        paste0("Prepare working on parallel using `", NCores_GPP, "` cores."),
        Level = 1)

      withr::local_options(
        future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
        future.seed = TRUE)
      c1 <- snow::makeSOCKcluster(NCores_GPP)
      on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
      future::plan("future::cluster", workers = c1, gc = TRUE)
      on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)

      IASDT.R::CatTime("Prepare GPP knots", Level = 2)
      GPP_Knots <- future.apply::future_lapply(
        X = GPP_Dists * 1000,
        FUN = function(x) {
          IASDT.R::Mod_PrepKnots(
            Coords = DT_xy, MinDist = x, MinLF = MinLF, MaxLF = MaxLF,
            Alphapw = Alphapw)
        },
        future.scheduling = Inf, future.seed = TRUE,
        future.globals = c("DT_xy", "GPP_Dists", "MaxLF", "MinLF", "Alphapw"),
        future.packages = c("dplyr", "sf", "Hmsc", "jsonify", "magrittr")) %>%
        stats::setNames(paste0("GPP_", GPP_Dists))

      # Stopping cluster
      snow::stopCluster(c1)
      future::plan("future::sequential", gc = TRUE)

    } else {

      IASDT.R::CatTime("Working sequentially")

      GPP_Knots <- purrr::map(
        .x = GPP_Dists * 1000,
        .f = ~IASDT.R::Mod_PrepKnots(
          Coords = DT_xy, MinDist = .x, MinLF = MinLF, MaxLF = MaxLF,
          Alphapw = Alphapw)) %>%
        stats::setNames(paste0("GPP_", GPP_Dists))
    }


    ## Plotting knot location ----
    if (GPP_Plot) {
      IASDT.R::CatTime("Plotting GPP knots", Level = 1)

      GridR <- terra::unwrap(IASDT.R::LoadAs(Path_GridR))

      GridR <- sf::st_as_sf(
        x = data.frame(DT_xy), coords = c("x", "y"), crs = 3035) %>%
        terra::rasterize(GridR) %>%
        terra::as.factor() %>%
        stats::setNames("GridR")

      EU_Bound <- IASDT.R::LoadAs(EU_Bound) %>%
        magrittr::extract2("Bound_sf_Eur_s") %>%
        magrittr::extract2("L_03") %>%
        suppressWarnings()

      Knots_Plots <- purrr::map(
        .x = GPP_Dists,
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
        filename = IASDT.R::Path(Path_Model, "knot_Locations.pdf"),
        width = 12, height = 13.25, onefile = TRUE)
      invisible(purrr::map(Knots_Plots, print))
      grDevices::dev.off()
    }

    if (GPP_Save) {
      IASDT.R::CatTime("Saving GPP knots data", Level = 1)
      save(GPP_Knots, file = IASDT.R::Path(Path_Model, "GPP_Knots.RData"))
    }
  } else {

    IASDT.R::CatTime("Models will be fitted without spatial random effect")
    GPP_Knots <- studyDesign <- DT_xy <- NULL
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Define the initial models -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Define the initial models")

  if (GPP) {

    ModelVariants <- tidyr::expand_grid(M_thin = thin, M_samples = samples) %>%
      dplyr::mutate(M_transient = M_thin * transientFactor)

    Model_Info <- tibble::tibble(rL = GPP_Dists, rL2 = GPP_Knots) %>%
      # Combinations of rL and Tree
      tidyr::expand_grid(Tree = Tree) %>%
      dplyr::mutate(
        # Model name
        M_Name_init = paste0("GPP", rL, "_", Tree),
        # Save initial models
        M_Init_Path = purrr::map2_chr(
          .x = M_Name_init, .y = rL2,
          .f = ~{

            PathOut <- IASDT.R::Path(
              Path_Model, paste0("InitMod_", .x, ".RData"))

            if (!file.exists(PathOut)) {
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
          }),
        rL2 = NULL) %>%
      # add all combinations of thinning, number of samples and transient
      tidyr::expand_grid(ModelVariants) %>%
      dplyr::mutate(
        M_HPC = purrr::pmap(
          .l = list(M_Name_init, M_thin, M_samples),
          .f = function(M_Name_init, M_thin, M_samples) {

            M_Name_Fit <- paste0(M_Name_init, "_samp", M_samples, "_th", M_thin)

            M4HPC_Path <- IASDT.R::Path(
              Path_Model, "InitMod4HPC", paste0("InitMod_", M_Name_Fit, ".rds"))

            return(list(M_Name_Fit = M_Name_Fit, M4HPC_Path = M4HPC_Path))

          })) %>%
      tidyr::unnest_wider("M_HPC")

  } else {

    # Non-spatial model
    Model_Info <- tidyr::expand_grid(M_thin = thin, M_samples = samples) %>%
      dplyr::mutate(M_transient = M_thin * transientFactor) %>%
      tidyr::expand_grid(Tree = Tree) %>%
      dplyr::mutate(
        # Model name
        M_Name_init = paste0("NonSpatial_", Tree),
        # Save initial models
        M_Init_Path = purrr::map_chr(
          .x = M_Name_init,
          .f = ~{
            PathOut <- IASDT.R::Path(
              Path_Model, paste0("InitMod_", .x, ".RData"))

            InitModExists <- IASDT.R::CheckData(PathOut)

            if (isFALSE(InitModExists)) {
              if (stringr::str_detect(.x, "_Tree$")) {
                Tree <- plant.tree
              } else {
                Tree <- NULL
              }

              InitModel <- Hmsc::Hmsc(
                Y = DT_y, XFormula = Form_x, XData = DT_x,
                distr = "probit", phyloTree = Tree)

              IASDT.R::SaveAs(
                InObj = InitModel, OutObj = paste0("InitMod_", .x),
                OutPath = PathOut)
            }
            return(PathOut)
          }),

        M_HPC = purrr::pmap(
          .l = list(M_Name_init, M_thin, M_samples),
          .f = function(M_Name_init, M_thin, M_samples) {

            M_Name_Fit <- paste0(M_Name_init, "_samp", M_samples, "_th", M_thin)

            M4HPC_Path <- IASDT.R::Path(
              Path_Model, "InitMod4HPC", paste0("InitMod_", M_Name_Fit, ".rds"))

            return(list(M_Name_Fit = M_Name_Fit, M4HPC_Path = M4HPC_Path))

          })) %>%
      tidyr::unnest_wider("M_HPC")
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Prepare and save unfitted models -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Save unfitted models")

  if (OverwriteRDS) {
    IASDT.R::CatTime("Processing all model variants", Level = 1)
  } else {
    NMod2Export <- Model_Info %>%
      dplyr::filter(!file.exists(M4HPC_Path)) %>%
      nrow()
    if (NMod2Export == 0) {
      IASDT.R::CatTime(
        "All model variants were already available as RDS files", Level = 1)
    } else {
      IASDT.R::CatTime(
        paste0(
          NMod2Export, " model variants need to be exported as RDS files"),
        Level = 1)
    }
  }

  # `InitFitFun` - Function to start sampling
  InitFitFun <- function(ID) {
    CurrPath <- Model_Info$M4HPC_Path[ID]

    if (OverwriteRDS && file.exists(CurrPath)) {
      fs::file_delete(CurrPath)
    }

    if (isFALSE(OverwriteRDS) && file.exists(CurrPath) &&
        IASDT.R::CheckData(CurrPath)) {
      return(invisible(NULL))
    }

    Try <- 0
    while (Try < 6) {
      Try <- Try + 1
      Model <- Hmsc::sampleMcmc(
        hM = IASDT.R::LoadAs(Model_Info$M_Init_Path[ID]),
        samples = Model_Info$M_samples[ID],
        thin = Model_Info$M_thin[ID],
        transient = Model_Info$M_transient[ID],
        nChains = NChains, verbose = verbose, engine = "HPC")

      if (ToJSON) {
        Model <- jsonify::to_json(Model)
      }

      saveRDS(Model, file = CurrPath)

      if (file.exists(CurrPath) && IASDT.R::CheckRDS(CurrPath)) {
        break
      }
    }
    invisible(gc())
    return(invisible(NULL))
  }

  # Implement `InitFitFun` function: start sampling and save output files
  if (NCores > 1) {

    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
      future.seed = TRUE)
    c1 <- snow::makeSOCKcluster(min(NCores, nrow(Model_Info)))
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("future::cluster", workers = c1, gc = TRUE)
    on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)

    Model_Process <- future.apply::future_lapply(
      X = seq_len(nrow(Model_Info)),
      FUN = InitFitFun, future.scheduling = Inf, future.seed = TRUE,
      future.globals = c(
        "InitFitFun", "Model_Info", "OverwriteRDS",
        "verbose", "NChains", "ToJSON"),
      future.packages = c("Hmsc", "jsonify", "IASDT.R"))

    snow::stopCluster(c1)
    future::plan("future::sequential", gc = TRUE)

  } else {
    Model_Process <- purrr::map(
      .x = seq_len(nrow(Model_Info)), .f = InitFitFun)
  }
  rm(Model_Process, envir = environment())

  # Which models failed to be exported as RDS files after 5 trials
  Failed2Export <- dplyr::filter(Model_Info, !file.exists(M4HPC_Path))
  if (nrow(Failed2Export) == 0) {
    IASDT.R::CatTime("All model variants were exported as RDS files", Level = 1)
  } else {
    IASDT.R::CatTime(
      paste0(
        nrow(Failed2Export),
        " model variants failed to be exported to rds files after 5 tries."),
      Level = 1)
    save(Failed2Export, file = IASDT.R::Path(Path_Model, "Failed2Export.RData"))
    readr::write_tsv(
      x = Failed2Export, file = IASDT.R::Path(Path_Model, "Failed2Export.txt"))
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Prepare Hmsc-HPC fitting commands -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Prepare Hmsc-HPC fitting commands")

  Model_Info <- dplyr::mutate(Model_Info, Chain = list(seq_len(NChains))) %>%
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
          Path_Model2 <- IASDT.R::Path(
            Path_Model, "InitMod4HPC", basename(M4HPC_Path))

          # Path for posterior sampling
          Path_Post <- IASDT.R::Path(
            Path_Model, "Model_Fitting_HPC",
            paste0(M_Name_Fit, "_Chain", Chain, "_post.rds"))

          # Path for progress
          Path_Prog <- IASDT.R::Path(
            Path_Model, "Model_Fitting_HPC",
            paste0(M_Name_Fit, "_Chain", Chain, "_Progress.txt"))

          Post_Missing <- isFALSE(file.exists(Path_Post))

          # File path for the python script
          Path_Model2_4cmd <- IASDT.R::NormalizePath(Path_Model2)
          Path_Post_4cmd <- IASDT.R::NormalizePath(Path_Post)
          Path_Prog_4cmd <- IASDT.R::NormalizePath(Path_Prog)
          Path_Python <- IASDT.R::NormalizePath(Path_Python)

          # `TF_ENABLE_ONEDNN_OPTS=0` is used to disable the following warning:
          #
          # I tensorflow/core/util/port.cc:113] oneDNN custom operations are on.
          # You may see slightly different numerical results due to
          # floating-point round-off errors from different computation orders.
          # To turn them off, set the environment variable
          # `TF_ENABLE_ONEDNN_OPTS=0`.
          #
          # `export TF_CPP_MIN_LOG_LEVEL=3` is used to reduce debug output from
          # Tensorflow

          Command_HPC <- paste0(
            # Not needed now as this now added to the `setup-env.sh` file
            # "export TF_CPP_MIN_LOG_LEVEL=3; export TF_ENABLE_ONEDNN_OPTS=0; ",

            "/usr/bin/time -v ",

            # Not needed as the python path is exported - check `setup-env.sh`
            # Path_Python,

            "python3 -m hmsc.run_gibbs_sampler",
            " --input ", Path_Model2_4cmd,
            " --output ", Path_Post_4cmd,
            " --samples ", M_samples,
            " --transient ", M_transient,
            " --thin ", M_thin,
            " --verbose ", verbose,
            " --chain ", (Chain - 1),
            " --fp ", Precision,
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
            " --verbose ", verbose,
            " --chain ", (Chain - 1),
            " --fp ", Precision,
            " > ", Path_Prog_4cmd,
            " 2>&1")

          list(
            M4HPC_Path_LUMI = Path_Model2, Post_Path = Path_Post,
            Post_Missing = Post_Missing, Path_ModProg = Path_Prog,
            Command_HPC = Command_HPC, Command_WS = Command_WS) %>%
            return()

        })) %>%
    tidyr::unnest_wider("M_Chain")

  # # |||||||||||||||||||||||||||||||||||
  # # Skip fitted models -----
  # # |||||||||||||||||||||||||||||||||||

  if (SkipFitted) {
    IASDT.R::CatTime("Skip fitted models")
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

  IASDT.R::CatTime("Save model fitting commands to text file(s)")

  IASDT.R::CatTime("Save fitting commands for windows PC", Level = 1)
  f <- file(
    description = IASDT.R::Path(Path_Model, "Commands_All_Windows.txt"),
    open = "wb")
  on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
  cat(Models2Fit_WS, sep = "\n", append = FALSE, file = f)
  close(f)


  IASDT.R::CatTime("Save fitting commands for HPC", Level = 1)
  NJobs <- length(Models2Fit_HPC)

  if (NJobs > NumArrayJobs) {
    NSplits <- ceiling((NJobs / NumArrayJobs))
    IDs <- IASDT.R::SplitVector(Vector = seq_len(NJobs), NSplit = NSplits)
  } else {
    NSplits <- 1
    IDs <- list(seq_len(NJobs))
  }

  # Save all fitting commands to single file
  IASDT.R::CatTime("Save all fitting commands to single file", Level = 1)
  f <- file(
    description = IASDT.R::Path(Path_Model, "Commands_All.txt"),
    open = "wb")
  on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
  cat(Models2Fit_HPC, sep = "\n", append = FALSE, file = f)
  close(f)

  # Save model fitting commands for batch SLURM jobs
  IASDT.R::CatTime(
    "Save model fitting commands for batch SLURM jobs", Level = 1)
  IASDT.R::CatTime(
    paste0("Models will be fitted in ", NSplits, " SLURM job(s)"), Level = 2)

  purrr::walk(
    .x = seq_len(NSplits),
    .f = function(x) {

      if (NSplits > 1) {
        CommandFile <- IASDT.R::Path(
          Path_Model, paste0("Commands2Fit_", x, ".txt"))
      } else {
        CommandFile <- IASDT.R::Path(Path_Model, "Commands2Fit.txt")
      }

      # create connection to SLURM file. This is better than using sink to have
      # a platform independent file (here, to maintain a linux-like new line
      # ending)
      f <- file(CommandFile, open = "wb")
      on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
      cat(Models2Fit_HPC[IDs[[x]]], sep = "\n", append = FALSE, file = f)
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
      ModelDir = Path_Model, JobName = JobName, MemPerCpu = MemPerCpu,
      Time = Time, EnvFile = EnvFile, FromHPC = FromHPC,
      Path_Hmsc = Path_Hmsc, ...)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Save small datasets prepared in the function ------
  # # |||||||||||||||||||||||||||||||||||

  DT_Split <- list(
    DT_All = DT_All, DT_y = DT_y, Form_x = Form_x, DT_x = DT_x, XVars = XVars,
    DT_CV = DT_CV, PhyloTree = PhyloTree, plant.tree = plant.tree, Tree = Tree,
    studyDesign = studyDesign, DT_xy = DT_xy, GPP_Knots = GPP_Knots)

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
    InObj = DT_Split, OutObj = OutObjName,
    OutPath = IASDT.R::Path(Path_Model, paste0(OutObjName, ".RData")))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatDiff(InitTime = .StartTime)

  return(invisible(NULL))
}
