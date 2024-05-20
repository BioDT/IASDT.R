# |---------------------------------------------------| #
# PrepMod4HPC ----
# |---------------------------------------------------| #

#' Prepare initial models in R to be used by the Hmsc-HPC
#'
#' Prepare initial models in R to be used by the Hmsc-HPC
#'
#' @param Hab_Abb String. Habitat type. This has to be one of the following: c("0", "1", "2", "3", "4a", "4b", "5", "6", "8", "10", "12a", "12b"). "0" means prepare data irrespective of the habitat type
#' @param Path_Data String. Path to read modelling data (without trailing slash)
#' @param Path_Model String. Path to save all the output, including the to be fitted models (without trailing slash)
#' @param MinPresGrids Integer. Minimum number of presence grid cells per species. Only species with ≥ this number will be considered
#' @param Path_EnvFile String. Path to read the environment variables. Default value: `.env`
#' @param GPP_Dists Integer. Distance in kilometer for the distance between the knots and the minimum distance of a knot to the nearest data point. The same value will be used for the `knotDist`	and `minKnotDist`	arguments of the `Hmsc::constructKnots` function.
#' @param GPP_Save Logical. Save the resulted knots as RData. Default: `TRUE`
#' @param GPP_Plot Logical. Plot the coordinates of the sampling units and the knots. Default: `TRUE`
#' @param XVars String. Variables to be used in the model. Default value: `NULL`, which means to use the following variables: bio4, bio6, bio8, bio12, bio15, bio18, RoadRailLog, BiasLog
#' @param PhyloTree Logical. Fit model variants with phylogenetic trees. Default: `TRUE`. At least one of `PhyloTree` and `NoPhyloTree` can not be `FALSE`.
#' @param NoPhyloTree Logical. Fit model variants without phylogenetic trees. Default: `TRUE`. At least one of `PhyloTree` and `NoPhyloTree` can not be `FALSE`.
#' @param NParallel Number of parallel cores for palatalization
#' @param nChains Integer. Number of model chains
#' @param thin Integer. Value(s) for thinning
#' @param samples Integer. Value(s) for the number of MCMC samples
#' @param transientFactor Integer. Transient multiplication factor. The value of transient = `transientFactor * thin`.
#' @param verbose Integer. How often the results of the MCMC sampling be reported
#' @param ModelCountry String. Filter observations to the following country list. Default: `NULL`, which means prepare data for the whole Europe
#' @param MinPresPerCountry Integer. Minimum of grid cells for the selected country/countries for species to be considered in the models. Effective only if a valid `ModelCountry` is provided.
#' @name PrepMod4HPC
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

PrepMod4HPC <- function(
    Hab_Abb = NULL, Path_Data = NULL, Path_Model = NULL, MinPresGrids = NULL,
    Path_EnvFile = ".env",
    GPP_Dists = c(20, 30, 40, 50, 60), GPP_Save = TRUE, GPP_Plot = TRUE,
    XVars = NULL, PhyloTree = TRUE, NoPhyloTree = TRUE, NParallel = 8,
    nChains = 4, thin = c(5, 10, 20), samples = c(1000, 2000, 3000),
    transientFactor = 300, verbose = 1000,
    ModelCountry = NULL, MinPresPerCountry = 50) {

  # https://github.com/hmsc-r/HMSC/issues/180
  # https://github.com/hmsc-r/HMSC/issues/139

  .StartTime <- lubridate::now(tzone = "CET")

  Hab_Abb <- as.character(Hab_Abb)

  # # |||||||||||||||||||||||||||||||||||
  # # Load environment variables
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("\nLoad environment variables")
  if (file.exists(Path_EnvFile)) {
    readRenviron(Path_EnvFile)
  } else {
    MSG <- paste0(
      "Path for environment variables: ", Path_EnvFile, " was not found")
    stop(MSG)
  }

  Path_Hmsc <- Sys.getenv("DP_R_Mod_Path_Hmsc")
  Path_Python <- Sys.getenv("DP_R_Mod_Path_Python")
  Path_TaxaList <- Sys.getenv("DP_R_Mod_Path_TaxaList")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # CHECK input arguments
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Checking input arguments")

  AllArgs <- ls()
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  CharArgs <- c(
    "Path_Data", "Hab_Abb", "Path_Model", "Path_TaxaList",
    "Path_Hmsc", "Path_Python")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = CharArgs, Type = "character")

  LogicArgs <- c("GPP_Save", "GPP_Plot", "PhyloTree", "NoPhyloTree")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = LogicArgs, Type = "logical")

  NumericArgs <- c(
    "GPP_Dists", "NParallel", "nChains", "thin", "samples", "verbose",
    "MinPresGrids", "transientFactor")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = NumericArgs, Type = "numeric")

  rm(AllArgs, CharArgs, LogicArgs, NumericArgs)

  # # character arguments
  # MissingArgs <- list(
  #   Path_Data = Path_Data, Hab_Abb = Hab_Abb, Path_Model = Path_Model,
  #   Path_TaxaList = Path_TaxaList,
  #   Path_Hmsc = Path_Hmsc, Path_Python = Path_Python) %>%
  #   purrr::map(~inherits(.x, "character") && nchar(.x) > 0) %>%
  #   purrr::discard(.p = isTRUE) %>%
  #   names() %>%
  #   sort()
  # if (length(MissingArgs) > 0) {
  #   MSG <- paste0("The following argument(s) must be provided\n  >>  ",
  #                 paste0(MissingArgs, collapse = " | "))
  #   stop(MSG)
  # }
  #
  # # logical arguments
  # ArgsLgl <- list(
  #   GPP_Save = GPP_Save, GPP_Plot = GPP_Plot,
  #   PhyloTree = PhyloTree, NoPhyloTree = NoPhyloTree) %>%
  #   purrr::map(~inherits(.x, "logical")) %>%
  #   purrr::discard(.p = isTRUE) %>%
  #   names() %>%
  #   sort()
  # if (length(ArgsLgl) > 0) {
  #   MSG <- paste0("The following argument(s) must be logical\n  >>  ", ArgsLgl)
  #   stop(MSG)
  # }
  #
  # # numeric arguments
  # ArgsInt <- list(
  #   GPP_Dists = GPP_Dists, NParallel = NParallel, nChains = nChains,
  #   thin = thin, samples = samples, verbose = verbose,
  #   MinPresGrids = MinPresGrids, transientFactor = transientFactor) %>%
  #   purrr::map(~inherits(.x, "numeric")) %>%
  #   purrr::discard(.p = isTRUE) %>%
  #   names() %>%
  #   sort()
  # if (length(ArgsInt) > 0) {
  #   MSG <- paste0(
  #     "The following argument(s) must be numeric (integer)\n  >>  ", ArgsInt)
  #   stop(MSG)
  # }

  # Phylogenetic tree options
  if (PhyloTree == FALSE && NoPhyloTree == FALSE) {
    stop("At least one of PhyloTree or NoPhyloTree has to be true")
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # File paths
  # # |||||||||||||||||||||||||||||||||||

  fs::dir_create(file.path(Path_Model, "InitMod_HPC"))
  Path_ModelDT <- file.path(Path_Model, "Model_Info.RData")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Loading input data
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Loading input data")
  ValidHabAbbs <- c(0:3, "4a", "4b", 5, 6, 8, 10, "12a", "12b")
  if (magrittr::not(as.character(Hab_Abb) %in% ValidHabAbbs)) {
    MSG <- paste0("Hab_Abb has to be one of the following:\n >> ",
                  paste0(ValidHabAbbs, collapse = " | "))
    stop(MSG)
  }
  HabVal <- c(
    "0_All", "1_Forests", "2_Open_forests", "3_Scrub",
    "4a_Natural_grasslands", "4b_Human_maintained_grasslands",
    "5_Sandy", "6_Rocky", "8_Saline", "10_Wetland",
    "12a_Ruderal_habitats", "12b_Agricultural_habitats") %>%
    stringr::str_subset(paste0("^", as.character(Hab_Abb), "_"))

  Path_Data2 <- file.path(
    Path_Data,
    paste0("ModelDT_", MinPresGrids, "Grids_",
           stringr::str_remove(HabVal, "Hab_"), ".RData"))

  DT_All <- IASDT.R::LoadAs(Path_Data2)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Loading input data
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Subsetting data")

  if (magrittr::not(is.null(ModelCountry))) {

    ValidCountries <- ModelCountry %in% unique(DT_All$Country)
    if (all(ValidCountries)) {
      IASDT.R::CatTime(paste0("  >>  Subsetting data to ", paste0(ModelCountry, collapse = " & ")))
      Sample_ExclSp <- DT_All %>%
        dplyr::filter(Country %in% ModelCountry) %>%
        dplyr::summarise(dplyr::across(tidyselect::starts_with("Sp_"), sum)) %>%
        tidyr::pivot_longer(
          cols = tidyselect::everything(),
          names_to = "Sp", values_to = "NCells") %>%
        dplyr::filter(NCells < MinPresPerCountry) %>%
        dplyr::pull(Sp)
    } else {
      MSG <- paste0(
        "The following are invalid country names: ",
        paste0(ModelCountry[!ValidCountries], collapse = " & "))
      stop(MSG)
    }

    DT_All <- DT_All %>%
      dplyr::filter(Country %in% ModelCountry) %>%
      dplyr::select(-tidyselect::all_of(Sample_ExclSp), -CellCode, -Country)
  } else {
    IASDT.R::CatTime("  >>  No data subsetting was implemented")
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Cross-validation
  # IASDT.R::GetCV(NR = 4, NC = 3)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Response - Y matrix
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Response - Y matrix")
  DT_y <- DT_All %>%
    dplyr::select(tidyselect::starts_with("Sp_")) %>%
    as.data.frame()

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Predictors - X matrix
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Predictors - X matrix")
  if (is.null(XVars)) {
    XVars <- IASDT.R::cc(
      bio4, bio6, bio8, bio12, bio15, bio18, RoadRailLog, BiasLog)
  }

  # The formula object becomes too large (up to > 2GB!) if created within a function. Setting the environment of the formula as an empty environment release this unnecessary size.
  # https://stackoverflow.com/questions/66241212
  Form_x <- XVars %>%
    stringr::str_c(collapse = " + ") %>%
    stringr::str_c("~ ", .) %>%
    stats::as.formula(env = baseenv())

  DT_x <- DT_All %>%
    dplyr::select(tidyselect::all_of(XVars)) %>%
    as.data.frame()

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Phylogenetic tree data
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Phylogenetic tree data")

  if (PhyloTree) {
    # Taxonomy as a proxy for phylogeny
    plant.tree <- file.path(Path_TaxaList, "Species_List_ID.txt") %>%
      readr::read_tsv(show_col_types = FALSE) %>%
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
  # # Preparing working on parallel
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime(paste0(
    "Preparing working on parallel (", NParallel, " cores)"))
  c1 <- snow::makeSOCKcluster(NParallel)
  future::plan(future::cluster, workers = c1, gc = TRUE)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # spatial info / random effect
  # # |||||||||||||||||||||||||||||||||||

  # Hmsc-hpc does not help much to speed up NNGP models
  # I will only proceed with GPP models

  IASDT.R::CatTime("Spatial info / random effect")
  studyDesign <- data.frame(sample = as.factor(seq_len(nrow(DT_x))))

  DT_xy <- DT_All %>%
    dplyr::select(x, y) %>%
    as.matrix()
  rownames(DT_xy) <- studyDesign$sample
  rm(DT_All)

  IASDT.R::CatTime("  >>  Prepare GPP knots")

  invisible(snow::clusterEvalQ(
    cl = c1, IASDT.R::LoadPackages(dplyr, sf, Hmsc, jsonify, magrittr)))
  snow::clusterExport(
    cl = c1, list = c("DT_xy", "GPP_Dists"), envir = environment())

  GPP_Knots <- snow::parLapply(
    cl = c1, x = GPP_Dists * 1000,
    fun = (\(x) IASDT.R::PrepKnots(Coords = DT_xy, MinDist = x))) %>%
    stats::setNames(paste0("GPP_", GPP_Dists))

  if (GPP_Plot) {
    IASDT.R::CatTime("  >>  Plotting GPP knots")
    grDevices::pdf(file.path(Path_Model, "knot_Locations.pdf"), width = 16, height = 20)
    graphics::par(oma = c(0.5, 0.5, 0.5, 0.5), mar = c(0.5, 0.5, 2, 0.5))
    purrr::map(
      .x = GPP_Dists,
      .f = ~{
        plot(DT_xy, pch = 19, cex = 0.3, axes = FALSE,
             main = paste0(.x, " km"), cex.main = 2)
        points(GPP_Knots[[paste0("GPP_", .x)]]$sKnot,
               pch = 4, col = "blue", cex = 1, font = 1)
      })
    grDevices::dev.off()
  }

  if (GPP_Save) {
    IASDT.R::CatTime("  >>  Saving GPP knots data")
    save(GPP_Knots, file = file.path(Path_Model, "GPP_Knots.RData"))
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Define the models
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Define the models")

  ModelVariants <- tidyr::expand_grid(M_thin = thin, M_samples = samples) %>%
    dplyr::mutate(M_transient = M_thin * transientFactor)

  Model_Info <- tibble::tibble(rL = GPP_Dists, rL2 = GPP_Knots) %>%
    tidyr::expand_grid(Tree = c("Tree", "NoTree")) %>%
    dplyr::mutate(
      M_Name_init = paste0("GPP", rL, "_", Tree),
      M_Init_Path = purrr::map2_chr(
        .x = M_Name_init, .y = rL2,
        .f = ~{

          PathOut <- file.path(Path_Model, paste0("InitMod_", .x, ".RData"))

          if (file.exists(PathOut)) {
            InitModel <- LoadAs(PathOut)
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
        }, .progress = FALSE)) %>%
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

  rm(GPP_Knots)
  # invisible(snow::clusterEvalQ(c1, gc()))

  # # |||||||||||||||||||||||||||||||||||
  # # Save unfitted models
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Save unfitted models")
  DT2Export <- Model_Info %>%
    dplyr::filter(magrittr::not(file.exists(M4HPC_Path))) %>%
    dplyr::select(M4HPC_Path, M_samples, M_thin, M_transient, M_Init_Path)

  if (nrow(DT2Export) > 0) {

    ExportModel <- function(x) {
      M4HPC_Path <- DT2Export$M4HPC_Path[x]
      M_samples <- DT2Export$M_samples[x]
      M_thin <- DT2Export$M_thin[x]
      M_transient <- DT2Export$M_transient[x]
      M_Init <- IASDT.R::LoadAs(DT2Export$M_Init_Path[x])

      Model <- Hmsc::sampleMcmc(
        hM = M_Init, samples = M_samples, thin = M_thin,
        transient = M_transient, nChains = nChains,
        verbose = verbose, engine = "HPC") %>%
        jsonify::to_json()

      saveRDS(Model, file = M4HPC_Path)
      return(pryr::object_size(Model))
    }

    Export_Results <- seq_len(nrow(DT2Export)) %>%
      furrr::future_map(
        .f = purrr::safely(.f = ExportModel, otherwise = NULL, quiet = FALSE),
        .options = furrr::furrr_options(seed = TRUE, scheduling = Inf),
        .progress = FALSE)
    save(Export_Results,
         file = file.path(Path_Model, "Export_Results.RData"))

    FailedExports <- DT2Export %>%
      dplyr::filter(magrittr::not(file.exists(M4HPC_Path)))

    if (nrow(FailedExports) > 0) {
      paste0("  >>  ", nrow(FailedExports),
             "model variants failed to be exported to rds files.") %>%
        IASDT.R::CatTime()
      save(FailedExports, file = file.path(Path_Model, "FailedExports.RData"))
    }
  }

  # # # |||||||||||||||||||||||||||||||||||
  # # # Stopping cluster
  # # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Stopping cluster")
  snow::stopCluster(c1)
  invisible(gc())

  # # |||||||||||||||||||||||||||||||||||
  # # Prepare commands for running on HPC
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Prepare commands for running on HPC")

  Model_Info <- Model_Info %>%
    dplyr::mutate(Chain = list(seq_len(nChains))) %>%
    tidyr::unnest_longer("Chain") %>%
    dplyr::arrange(M_Name_Fit) %>%

    dplyr::mutate(
      M_Chain = purrr::pmap(
        .l = list(M_Name_Fit, M4HPC_Path, Chain,
                  M_transient, M_samples, M_thin),
        .f = function(M_Name_Fit, M4HPC_Path, Chain,
                      M_transient, M_samples, M_thin) {

          M4HPC_Path2 <- file.path(
            Path_Model, "InitMod_HPC", basename(M4HPC_Path))
          Post_Path <- file.path(
            Path_Model, paste0(M_Name_Fit, "_Chain", Chain, "_post.rds"))
          Path_ModPorg <- file.path(
            Path_Model, paste0(M_Name_Fit, "_Chain", Chain, "_Progress.txt"))
          Post_Missing <- magrittr::not(file.exists(Post_Path))

          Command <- paste0(
            "export TF_CPP_MIN_LOG_LEVEL=3; export PYTHONPATH=", Path_Hmsc,
            ":$PYTHONPATH; /usr/bin/time -v ", Path_Python,
            " -m hmsc.run_gibbs_sampler",
            " --input ", shQuote(M4HPC_Path2),
            " --output ", shQuote(Post_Path),
            " --samples ", M_samples,
            " --transient ", M_transient,
            " --thin ", M_thin,
            " --verbose ", verbose,
            " --chain ", Chain - 1,
            " >& ", shQuote(Path_ModPorg))

          list(
            M4HPC_Path_LUMI = M4HPC_Path2,
            Post_Path = Post_Path, Post_Missing = Post_Missing,
            Path_ModPorg = Path_ModPorg, Command = Command) %>%
            return()
        })) %>%
    tidyr::unnest_wider("M_Chain")

  IASDT.R::CatTime("Save commands in a text file")

  cat(Model_Info$Command, sep = "\n", append = FALSE,
      file = file.path(Path_Model, "Commands_All.txt"))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # # |||||||||||||||||||||||||||||||||||
  # # Save data to disk
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Save data to disk")
  Model_Info <- Model_Info %>%
    tidyr::nest(
      Post_Path = Post_Path, Path_ModPorg = Path_ModPorg,
      Chain = Chain, Command = Command) %>%
    dplyr::mutate(
      Post_Path = purrr::map2(Post_Path, Chain, IASDT.R::SetChainName),
      Chain = purrr::map2(Chain, Chain, IASDT.R::SetChainName),
      Command = purrr::map2(Command, Chain, IASDT.R::SetChainName),
      Path_ModPorg = purrr::map2(Path_ModPorg, Chain, IASDT.R::SetChainName))

  save(Model_Info, file = Path_ModelDT)

  IASDT.R::CatDiff(.StartTime, CatInfo = FALSE)
  return(invisible(NULL))
}
