## |------------------------------------------------------------------------| #
# Mod_Prep4HPC ----
## |------------------------------------------------------------------------| #

#' Prepare initial models in R for model fitting by Hmsc-HPC
#'
#' This function prepares initial models in R for model fitting by Hmsc-HPC. It
#' involves data preparation, define spatial block cross-validation folds, model
#' initialization, GPP knots, and generating commands for running models on HPC.
#' It supports parallel processing, options to include/not include phylogenetic
#' tree data. The models will be fitted using Gaussian Predictive Process (GPP;
#' see [Tikhonov et al.](https://doi.org/10.1002/ecy.2929)) for more details)
#' via the [Hmsc-HPC](https://doi.org/10.1101/2024.02.13.580046) extension.
#'
#' @param Path_Model String (without trailing slash) specifying the path where
#'   all output, including models to be fitted, will be saved.
#' @param GPP_Dists Integer specifying the distance in kilometers for both the
#'   distance between knots and the minimum distance of a knot to the nearest
#'   data point. The GPP knots are prepared by the [IASDT.R::PrepKnots]
#'   function. The same value will be used for the `knotDist` and
#'   `minKnotDist`	arguments of the [Hmsc::constructKnots] function.
#' @param GPP_Save Logical indicating whether to save the resulted knots as
#'   `RData` file. Default: `TRUE`.
#' @param GPP_Plot Logical indicating whether to plot the coordinates of the
#'   sampling units and the knots in a pdf file. Default: `TRUE`.
#' @param EffortsAsPredictor Logical indicating whether to include the
#'   (log<sub>10</sub>) sampling efforts as predictor to the model. Default:
#'   `FALSE`.
#' @param RoadRailAsPredictor Logical indicating whether to include the
#'   (log<sub>10</sub>) sum of road and railway intensity as predictor to the
#'   model. Default: `TRUE`.
#' @param HabAsPredictor Logical indicating whether to include the
#'   (log<sub>10</sub>) percentage coverage of respective habitat type per grid
#'   cell as predictor to the model. Default: `TRUE`. Only valid if `Hab_Abb`
#'   not equals to `0`.
#' @param NspPerGrid Integer. Indicating the minimum number of species per grid
#'   cell for a grid cell to be include in the analysis. This is calculated
#'   after filtering grid cells by sampling efforts (`MinEffortsSp`) and
#'   filtering species by the number of presence grid cells per predictor
#'   (`PresPerVar`). If `NspPerGrid` = `NULL` or 1 (default), grid cells with
#'   at least one species presence will be considered in the models.
#' @param PhyloTree,NoPhyloTree Logical indicating whether to fit model variants
#'   with or without phylogenetic trees, respectively. The default of both
#'   arguments is `TRUE`, which means to fit a model variant with the respective
#'   option. If both `PhyloTree` and `NoPhyloTree` are `TRUE` (Default), models
#'   for both options will be fitted. At least one of `PhyloTree` and
#'   `NoPhyloTree` should be `TRUE`.
#' @param OverwriteRDS Logical. Indicating whether to overwrite previously
#'   exported RDS files for initial models. Default: `TRUE`.
#' @param NCores Integer specifying the number of parallel cores for
#'   parallelization. Default: 8 cores.
#' @param NChains Integer specifying the number of model chains. Default: 4.
#' @param thin Integer specifying the value(s) for thinning in MCMC sampling. If
#'   more than one value is provided, a separate model will be fitted at each
#'   value of thinning.
#' @param samples Integer specifying the value(s) for the number of MCMC
#'   samples. If more than one value is provided, a separate model will be
#'   fitted at each value of number of samples. Defaults to 1000.
#' @param transientFactor Integer specifying the transient multiplication
#'   factor. The value of `transient` will equal  the multiplication of
#'   `transientFactor` and `thin`. Default: 300.
#' @param verbose Integer specifying how often the results of the MCMC sampling
#'   should be reported. Default: `200`.
#' @param SkipFitted Logical indicating whether to skip already fitted models.
#'   Default: `TRUE`.
#' @param NArrayJobs Integer specifying the maximum allowed number of array
#'   jobs per SLURM file. Default: 210. See
#'   [LUMI documentation](https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/partitions)
#'   for more details.
#' @param ModelCountry String or vector of strings specifying the country or
#'   countries to filter observations by. Default: `NULL`, which means prepare
#'   data for the whole Europe.
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
#' @param Path_Python String specifying the path for Python.
#' @param ToJSON Logical indicating whether to convert unfitted models to JSON
#'   before saving to RDS file. Default: `FALSE`.
#' @param ... Additional parameters provided to the [IASDT.R::Mod_SLURM]
#'   function.
#' @name Mod_Prep4HPC
#' @inheritParams Mod_PrepData
#' @inheritParams PrepKnots
#' @inheritParams GetCV
#' @author Ahmed El-Gabbas
#' @return The function is used for its side effects of preparing data and
#'   models for HPC and does not return any value.
#' @details The function provides options for:
#'
#' - for which habitat types the models will be fitted
#' - excluding grid cells with very low sampling efforts (`MinEffortsSp`)
#' - selection of species based on minimum number of presence-grid cells
#'   (`PresPerVar` * number of predictors)
#' - optionally model fitting on specified list of countries: (`ModelCountry`)
#' - whether to exclude grid cells with few species (`NspPerGrid`)
#' - number of cross-validation folds
#' - options for whether or not to include phylogenetic information to the model
#' - different values for knot distance for GPP (`GPP_Dists`)
#' - which Bioclimatic variables to be uses in the models (`BioVars`)
#' - whether to include sampling efforts `EffortsAsPredictor`, percentage of
#'   respective habitat type per grid cell `HabAsPredictor`, and railway and
#'   road intensity per grid cell `RoadRailAsPredictor`
#' - Hmsc options (`NChains`, `thin`, `samples`, `transientFactor`, and `verbose`)
#' - prepare SLURM commands (`PrepSLURM`) and some specifications (e.g.
#'   `NArrayJobs`, `MemPerCpu`, `Time`, `JobName`)
#'
#'   The function reads the following environment variables:
#'   - **`DP_R_Grid`** (if `FromHPC = TRUE`) or
#'    **`DP_R_Grid_Local`** (if `FromHPC = FALSE`). The function reads
#'   the content of the `Grid_10_Land_Crop.RData` file from this path.
#'   - **`DP_R_TaxaInfo`** or **`DP_R_TaxaInfo_Local`** for the location of the
#'   `Species_List_ID.txt` file representing species information.
#'   - **`DP_R_EUBound_sf`** or **`DP_R_EUBound_sf_Local`** for the path of the
#'   `RData` file containing the country boundaries (`sf` object)
#'   - **`DP_R_PA`** or **`DP_R_PA_Local`**: The function reads the contents of
#'   the `Sp_PA_Summary_DF.RData` file from this path
#' @export

Mod_Prep4HPC <- function(
    Hab_Abb = NULL, Path_Model = NULL,
    MinEffortsSp = 100L, PresPerVar = 10L,
    EnvFile = ".env", GPP_Dists = NULL, GPP_Save = TRUE,
    GPP_Plot = TRUE, MinLF = NULL, MaxLF = NULL,
    BioVars = c("bio4", "bio6", "bio8", "bio12", "bio15", "bio18"),
    EffortsAsPredictor = FALSE, RoadRailAsPredictor = TRUE,
    HabAsPredictor = TRUE, NspPerGrid = 1L, ExcludeCult = TRUE,
    CV_NFolds = 4L, CV_NGrids = 20L, CV_NR = 2L, CV_NC = 2L, CV_Plot = TRUE,
    PhyloTree = TRUE, SaveData = TRUE,
    NoPhyloTree = TRUE, OverwriteRDS = TRUE, NCores = 8L, NChains = 4L,
    thin = NULL, samples = 1000L, transientFactor = 300L, verbose = 200L,
    SkipFitted = TRUE, NArrayJobs = 210L, ModelCountry = NULL,
    VerboseProgress = TRUE, FromHPC = TRUE, PrepSLURM = TRUE, MemPerCpu = NULL,
    Time = NULL, JobName = NULL, Path_Hmsc = NULL, Path_Python = NULL,
    ToJSON = FALSE, ...) {

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # Initial checking -----
  # # |||||||||||||||||||||||||||||||||||

  .StartTime <- lubridate::now(tzone = "CET")

  CheckNULL <- c(
    "Path_Model", "PresPerVar", "thin", "samples", "GPP_Dists",
    "MemPerCpu", "Path_Hmsc", "Path_Python", "Hab_Abb")
  IsNull <- purrr::map_lgl(CheckNULL, ~is.null(get(.x)))

  if (any(IsNull)) {
    stop(
      paste0(
        paste0("`", CheckNULL[which(IsNull)], "`", collapse = ", "),
        " can not be empty"),
      call. = FALSE)
  }

  Hab_Abb <- as.character(Hab_Abb)

  if (!all(is.numeric(GPP_Dists)) || any(GPP_Dists <= 0)) {
    stop("GPP_Dists should be numeric and greater than zero", call. = FALSE)
  }

  if (!all(is.numeric(samples)) || any(samples <= 0)) {
    stop("samples should be numeric and greater than zero", call. = FALSE)
  }

  if (!all(is.numeric(thin)) || any(thin <= 0)) {
    stop("thin should be numeric and greater than zero", call. = FALSE)
  }

  if (!all(is.numeric(PresPerVar)) || PresPerVar <= 0) {
    stop("`PresPerVar` should be numeric and greater than zero", call. = FALSE)
  }

  if (!all(is.numeric(MinEffortsSp)) || MinEffortsSp <= 0) {
    stop(
      "`MinEffortsSp` should be numeric and greater than zero",
      call. = FALSE)
  }

  if (!is.null(NspPerGrid) && (!is.numeric(NspPerGrid) || NspPerGrid < 1)) {
    stop(
      "`NspPerGrid` has to be either `NULL` or positive integer",
      call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  NCells <- Sp <- IAS_ID <- x <- y <- Country <- M_thin <- rL <-
    M_Name_init <- rL2 <- M_samples <- M4HPC_Path <- M_transient <-
    M_Init_Path <- M_Name_Fit <- Chain <- Post_Missing <- Command_HPC <-
    Command_WS <- Post_Path <- Path_ModProg <- TaxaInfoFile <-
    Path_Grid <- EU_Bound <- Path_PA <- SpeciesID <- Species_name <- PA <-
    Species_File <- NAME_ENGL <- NULL

  if (isFALSE(VerboseProgress)) {
    sink(file = nullfile())
    on.exit(sink(), add = TRUE)
  }

  IASDT.R::CatSep(Rep = 1, Extra1 = 1, Extra2 = 0, Char = "=")
  IASDT.R::CatTime("Preparing data for Hmsc-HPC models")
  IASDT.R::CatSep(Rep = 1, Extra1 = 0, Extra2 = 1, Char = "=")

  # # |||||||||||||||||||||||||||||||||||
  # Load/check environment variables -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Load/check environment variables")

  if (!file.exists(EnvFile)) {
    stop(paste0(
      "Path for environment variables: ", EnvFile, " was not found"),
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
    if (!file.exists(Path_Python)) {
      stop(
        paste0("Python executable does not exist: ", Path_Python),
        call. = FALSE)
    }
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  if (GPP_Plot) {
    Path_GridR <- file.path(Path_Grid, "Grid_10_Land_Crop.RData")
    if (!file.exists(Path_GridR)) {
      stop(
        paste0("Path for the Europe boundaries does not exist: ", Path_GridR),
        call. = FALSE)
    }
  }

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

  CharArgs <- c(
    "Hab_Abb", "Path_Model", "Path_Hmsc", "Path_Python")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = CharArgs, Type = "character")

  LogicArgs <- c(
    "GPP_Save", "GPP_Plot", "PhyloTree", "NoPhyloTree",
    "VerboseProgress", "ToJSON")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = LogicArgs, Type = "logical")

  NumericArgs <- c(
    "GPP_Dists", "NCores", "NChains", "thin", "samples", "verbose",
    "PresPerVar", "MinEffortsSp", "transientFactor", "CV_NFolds", "CV_NGrids",
    "CV_NR", "CV_NC")
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

  rm(AllArgs, CharArgs, LogicArgs, NumericArgs)

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # File paths - Creating missing paths ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("File paths - Creating missing paths")
  fs::dir_create(file.path(Path_Model, "InitMod4HPC"))
  fs::dir_create(file.path(Path_Model, "Model_Fitting_HPC"))
  # Also create directory for SLURM outputs
  fs::dir_create(file.path(Path_Model, "Model_Fitting_HPC", "JobsLog"))
  Path_ModelDT <- file.path(Path_Model, "Model_Info.RData")

  # # ..................................................................... ###

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
    XVars <- c(XVars, "BiasLog")
  }
  if (RoadRailAsPredictor) {
    XVars <- c(XVars, "RoadRailLog")
  }
  if (Hab_Abb != "0" && HabAsPredictor) {
    XVars <- c(XVars, "HabLog")
  }

  # minimum number of presence grids per species
  MinPresGrids <- PresPerVar * length(XVars)

  # # ..................................................................... ###

  # # |||||||||||||||||||||||||||||||||||
  # # Preparing input data -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Preparing input data")
  ValidHabAbbs <- c(0:3, "4a", "4b", 10, "12a", "12b")
  if (!(as.character(Hab_Abb) %in% ValidHabAbbs)) {
    stop(paste0("Hab_Abb has to be one of the following:\n >> ",
                paste0(ValidHabAbbs, collapse = " | ")), call. = FALSE)
  }

  HabVal <- c(
    "0_All", "1_Forests", "2_Open_forests", "3_Scrub",
    "4a_Natural_grasslands", "4b_Human_maintained_grasslands",
    "10_Wetland", "12a_Ruderal_habitats", "12b_Agricultural_habitats") %>%
    stringr::str_subset(paste0("^", as.character(Hab_Abb), "_"))

  IASDT.R::CatSep(Rep = 1, Extra1 = 1, Extra2 = 0)
  IASDT.R::CatTime("Preparing input data using IASDT.R::Mod_PrepData")
  IASDT.R::CatSep(Rep = 1, Extra1 = 0, Extra2 = 0)

  DT_All <- IASDT.R::Mod_PrepData(
    Hab_Abb = Hab_Abb, MinEffortsSp = MinEffortsSp, PresPerVar = PresPerVar,
    NVars = length(XVars), EnvFile = EnvFile, BioVars = BioVars,
    Path_Model = Path_Model, VerboseProgress = VerboseProgress,
    FromHPC = FromHPC, SaveData = SaveData, ExcludeCult = ExcludeCult)

  IASDT.R::CatSep(Rep = 1, Extra1 = 0, Extra2 = 1)

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
      dplyr::filter(NCells < MinPresGrids) %>%
      dplyr::pull(Sp)

    IASDT.R::CatTime(
      paste0(length(Sample_ExclSp), " species are excluded"), Level = 1)
    DT_All <- dplyr::filter(DT_All, Country %in% ModelCountry) %>%
      dplyr::select(-tidyselect::all_of(Sample_ExclSp))

    # # |||||||||||||||||||||||||||||||||||
    ## Plotting subsetted data -----
    # # |||||||||||||||||||||||||||||||||||

    IASDT.R::CatTime("Plotting subsetted data", Level = 1)
    GridSubset <- terra::rasterize(
      x = as.matrix(DT_All[, c("x", "y")]),
      y = terra::unwrap(IASDT.R::LoadAs(Path_GridR)))

    DT_Sp <- file.path(Path_PA, "Sp_PA_Summary_DF.RData")
    if (!file.exists(DT_Sp)) {
      stop(paste0(DT_Sp, " file does not exist"), call. = FALSE)
    }
    DT_Sp <- IASDT.R::LoadAs(DT_Sp)

    if (Hab_Abb == "0") {
      Hab_column <- NULL
    } else {
      Hab_column <- c(
        "Hab_1_Forests", "Hab_2_Open_forests", "Hab_3_Scrub",
        "Hab_4a_Natural_grasslands", "Hab_4b_Human_maintained_grasslands",
        "Hab_10_Wetland", "Hab_12a_Ruderal_habitats",
        "Hab_12b_Agricultural_habitats") %>%
        stringr::str_subset(paste0("_", as.character(Hab_Abb), "_"))

      DT_Sp <- dplyr::filter(DT_Sp, !!as.symbol(Hab_column))
    }

    R_Sp_sum <- DT_Sp %>%
      dplyr::mutate(SpeciesID = paste0("Sp_", SpeciesID)) %>%
      dplyr::select(SpeciesID, Species_name, Species_File) %>%
      dplyr::filter(
        SpeciesID %in% stringr::str_subset(names(DT_All), "^Sp_")) %>%
      dplyr::mutate(
        PA = purrr::map2(
          .x = Species_File, .y = SpeciesID,
          .f = ~{
            R <- file.path(
              Path_PA, "RData", stringr::str_c(.x, "_PA.RData")) %>%
              IASDT.R::LoadAs() %>%
              terra::unwrap() %>%
              magrittr::extract2("PA") %>%
              stats::setNames(paste0("Sp_", .y))
          })) %>%
      dplyr::pull(PA) %>%
      terra::rast() %>%
      sum(na.rm = TRUE) %>%
      terra::crop(GridSubset) %>%
      terra::mask(GridSubset)

    EU_Bound_sub <- IASDT.R::LoadAs(EU_Bound) %>%
      magrittr::extract2("Bound_sf_Eur") %>%
      magrittr::extract2("L_03") %>%
      dplyr::filter(NAME_ENGL %in% ModelCountry)

    R_Sp_sumP <- terra::classify(R_Sp_sum, cbind(0, NA))
    Limits <- terra::trim(R_Sp_sumP) %>%
      terra::ext() %>%
      as.vector()
    NSpPerGrid_Sub <- ggplot2::ggplot() +
      tidyterra::geom_spatraster(data = R_Sp_sumP) +
      tidyterra::scale_fill_whitebox_c(
        na.value = "transparent", palette = "bl_yl_rd", name = NULL) +
      ggplot2::labs(
        # title = "Number of presence species per grid cell",
        title = paste0(
          '<span style="color:blue; font-size:20px;"><b>',
          "Number of IAS per grid cell to be used in the models</b></span>",
          '<span style="color:black; font-size:16px;"> (',
          stringr::str_remove(Hab_column, "Hab_"), ")</span>"),
        subtitle = paste0(
          "Only species within &#8805;", MinPresGrids,
          " presence grid cells in the selected country/countries are shown")) +
      ggplot2::geom_sf(
        data = EU_Bound_sub, fill = "transparent", colour = "black") +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = Limits[c(3, 4)]) +
      ggplot2::scale_x_continuous(expand = c(0, 0), limits = Limits[c(1, 2)]) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(0.05, 0, 0, 0, "cm"),
        plot.title = ggtext::element_markdown(
          size = 16, hjust = 0,
          margin = ggplot2::margin(0, 0, 0.1, 0, "cm")),
        plot.subtitle = ggtext::element_markdown(size = 14, hjust = 0),
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.95, 0.9),
        legend.key.size = grid::unit(0.8, "cm"))

    # Relative JPEG height
    DimX <- Limits[2] - Limits[1]
    DimY <- Limits[4] - Limits[3]
    PlotHeight <- (DimY * 25) / DimX

    ggplot2::ggsave(
      plot = NSpPerGrid_Sub, width = 25, height = PlotHeight,
      units = "cm", dpi = 600,
      filename = file.path(Path_Model, "NSpPerGrid_Sub.jpeg"))

    rm(
      Limits, NSpPerGrid_Sub, R_Sp_sum, R_Sp_sumP,
      EU_Bound_sub, DT_Sp, GridSubset)

  } else {
    IASDT.R::CatTime("No data subsetting was implemented", Level = 1)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Exclude grid cells with low number of presences -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Exclude grid cells with low number of presences")

  if (is.null(NspPerGrid) || NspPerGrid == 1) {
    IASDT.R::CatTime(
      "All grid cells with at least single species presence will be considered",
      Level = 1)
  } else {
    if (NspPerGrid > 1) {
      EmptyGridsID <- dplyr::select(DT_All, tidyselect::starts_with("Sp_")) %>%
        rowSums() %>%
        magrittr::is_less_than(as.integer(NspPerGrid)) %>%
        which() %>%
        magrittr::multiply_by(-1)

      if (length(EmptyGridsID) > 0) {
        IASDT.R::CatTime(
          paste0("Excluding grid cells with < ", NspPerGrid), Level = 1)
        DT_All <- dplyr::slice(DT_All, EmptyGridsID)
      }
    }
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Cross-validation ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Prepare cross-validation folds")

  DT_All <- IASDT.R::GetCV(
    DT = DT_All, EnvFile = EnvFile, XVars = XVars,
    CV_NFolds = CV_NFolds, FromHPC = FromHPC, CV_NGrids = CV_NGrids,
    CV_NR = CV_NR, CV_NC = CV_NC, OutPath = Path_Model, CV_Plot = CV_Plot)

  # Save cross-validation data
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
  SpSummary <- file.path(Path_PA, "Sp_PA_Summary_DF.RData")
  if (!file.exists(SpSummary)) {
    stop(paste0(SpSummary, " file does not exist"), call. = FALSE)
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
    paste0("Models will be fitted using ", length(XVars), " predictors"),
    Level = 1)
  IASDT.R::CatTime(paste0(XVars, collapse = " + "), Level = 1)

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
    paste0(
      "Models will be fitted using ", paste0(Tree, collapse = " & ")),
    Level = 1)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Spatial info / random effect ------
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Spatial info / random effect")
  studyDesign <- data.frame(sample = as.factor(seq_len(nrow(DT_x))))

  DT_xy <- as.matrix(dplyr::select(DT_All, x, y))
  rownames(DT_xy) <- studyDesign$sample

  # Prepare GPP knots
  IASDT.R::CatTime("Preparing GPP knots", Level = 1)

  if (NCores > 1) {

    IASDT.R::CatTime(
      paste0("Prepare working on parallel using `", NCores, "` cores."),
      Level = 1)

    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE)

    future::plan("multisession", workers = NCores, gc = TRUE)
    on.exit(future::plan("sequential"), add = TRUE)

    GPP_Knots <- future.apply::future_lapply(
      X = GPP_Dists * 1000,
      FUN = function(x) {
        IASDT.R::PrepKnots(
          Coords = DT_xy, MinDist = x, MinLF = MinLF, MaxLF = MaxLF)
      },
      future.scheduling = Inf, future.seed = TRUE,
      future.globals = c("DT_xy", "GPP_Dists", "MaxLF", "MinLF"),
      future.packages = c("dplyr", "sf", "Hmsc", "jsonify", "magrittr")) %>%
      stats::setNames(paste0("GPP_", GPP_Dists))
    future::plan("sequential")

  } else {
    IASDT.R::CatTime("Working sequentially")

    GPP_Knots <- purrr::map(
      .x = GPP_Dists * 1000,
      .f = ~IASDT.R::PrepKnots(
        Coords = DT_xy, MinDist = .x, MinLF = MinLF, MaxLF = MaxLF)) %>%
      stats::setNames(paste0("GPP_", GPP_Dists))
  }

  ## Plotting knot location ----
  if (GPP_Plot) {
    IASDT.R::CatTime("Plotting GPP knots", Level = 1)

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
      magrittr::extract2("L_03") %>%
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
              " cells is ", .x, " km  &#8212; ", NKnots, " knots</span>")) +
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
    IASDT.R::CatTime("Saving GPP knots data", Level = 1)
    save(GPP_Knots, file = file.path(Path_Model, "GPP_Knots.RData"))
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Define the initial models -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Define the initial models")

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

          PathOut <- file.path(Path_Model, paste0("InitMod_", .x, ".RData"))

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

          M4HPC_Path <- file.path(
            Path_Model, "InitMod4HPC", paste0("InitMod_", M_Name_Fit, ".rds"))

          return(list(M_Name_Fit = M_Name_Fit, M4HPC_Path = M4HPC_Path))

        })) %>%
    tidyr::unnest_wider("M_HPC")

  # # |||||||||||||||||||||||||||||||||||
  # # Prepare and save unfitted models -----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Save unfitted models")

  if (OverwriteRDS) {
    IASDT.R::CatTime(paste0("Processing all model variants"), Level = 1)
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
        IASDT.R::CheckRDS(CurrPath)) {
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

    future::plan("multisession", workers = NCores, gc = TRUE)
    on.exit(future::plan("sequential"), add = TRUE)

    Model_Process <- future.apply::future_lapply(
      X = seq_len(nrow(Model_Info)),
      FUN = InitFitFun,
      future.scheduling = Inf, future.seed = TRUE,
      future.globals = c(
        "InitFitFun", "Model_Info", "OverwriteRDS",
        "verbose", "NChains", "ToJSON"),
      future.packages = c("Hmsc", "jsonify", "IASDT.R"))

    future::plan("sequential")

  } else {
    Model_Process <- purrr::map(
      .x = seq_len(nrow(Model_Info)), .f = InitFitFun)
  }


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
    save(Failed2Export, file = file.path(Path_Model, "Failed2Export.RData"))
    readr::write_tsv(
      x = Failed2Export, file = file.path(Path_Model, "Failed2Export.txt"))
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
          M4HPC_Path2 <- file.path(
            Path_Model, "InitMod4HPC", basename(M4HPC_Path))
          M4HPC_Path2_Win <- stringr::str_replace_all(M4HPC_Path2, "/", "\\\\")

          # Path for posterior sampling
          Post_Path <- file.path(
            Path_Model, "Model_Fitting_HPC",
            paste0(M_Name_Fit, "_Chain", Chain, "_post.rds"))
          Post_Path_Win <- stringr::str_replace_all(Post_Path, "/", "\\\\")

          # Path for progress
          Path_ModProg <- file.path(
            Path_Model, "Model_Fitting_HPC",
            paste0(M_Name_Fit, "_Chain", Chain, "_Progress.txt"))

          Post_Missing <- !file.exists(Post_Path)

          Exports <- paste0(
            "export TF_CPP_MIN_LOG_LEVEL=3; ",
            "export TF_ENABLE_ONEDNN_OPTS=0; ")

          # `TF_ENABLE_ONEDNN_OPTS=0` is used to disable the following warning:
          #
          # I tensorflow/core/util/port.cc:113] oneDNN custom operations are on.
          # You may see slightly different numerical results due to
          # floating-point round-off errors from different computation orders.
          # To turn them off, set the environment variable
          # `TF_ENABLE_ONEDNN_OPTS=0`.
          #
          # `export TF_CPP_MIN_LOG_LEVEL=3` is used to reduce debug output from
          # tensorflow

          Command_HPC <- paste0(
            Exports,
            "/usr/bin/time -v ",
            Path_Python,
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
            Exports,
            Path_Python,
            " -m hmsc.run_gibbs_sampler",
            " --input ", M4HPC_Path2_Win,
            " --output ", Post_Path_Win,
            " --samples ", M_samples,
            " --transient ", M_transient,
            " --thin ", M_thin,
            " --verbose ", verbose,
            " --chain ", (Chain - 1),
            " >& ", shQuote(Path_ModProg))

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
    description = file.path(Path_Model, "Commands_All_Windows.txt"),
    open = "wb")
  on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)
  cat(Models2Fit_WS, sep = "\n", append = FALSE, file = f)
  close(f)


  IASDT.R::CatTime("Save fitting commands for HPC", Level = 1)
  NJobs <- length(Models2Fit_HPC)

  if (NJobs > NArrayJobs) {
    NSplits <- ceiling((NJobs / NArrayJobs))
    IDs <- IASDT.R::SplitVector(Vector = seq_len(NJobs), NSplit = NSplits)
  } else {
    NSplits <- 1
    IDs <- list(seq_len(NJobs))
  }

  # Save all fitting commands to single file
  IASDT.R::CatTime("Save all fitting commands to single file", Level = 1)
  f <- file(
    description = file.path(Path_Model, "Commands_All.txt"),
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
        CommandFile <- file.path(Path_Model, paste0("Commands2Fit_", x, ".txt"))
      } else {
        CommandFile <- file.path(Path_Model, "Commands2Fit.txt")
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
      Path_Model = Path_Model, JobName = JobName, MemPerCpu = MemPerCpu,
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
    OutPath = file.path(Path_Model, paste0(OutObjName, ".RData")))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatDiff(
    InitTime = .StartTime, ChunkText = "Function summary", CatInfo = TRUE)

  return(invisible(NULL))
}
