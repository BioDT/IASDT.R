# # |------------------------------------------------------------------------| #
# IAS_Process ----
## |------------------------------------------------------------------------| #

#' Process IAS data
#'
#' This function processes Invasive Alien Species (IAS) data. The function
#' merges pre-processed distribution data from 3 data sources: GBIF
#' ([GBIF_Process]), EASIN ([EASIN_Process]), eLTER ([elTER_Process]).
#' The function prepares final species outputs in the form of 1) species
#' distribution as SpatRaster (`.RData`) and `.tif` using [IAS_Distribution]; 2)
#' summary table on the distribution of the species; and 3) JPEG files for the
#' species distribution showing the the sources of presence data [IAS_Plot]. It
#' further generates summary results in the form of maps and figures.
#' @param FromHPC Logical indicating whether the work is being done from HPC, to
#'   adjust file paths accordingly. Default: `TRUE`.
#' @param EnvFile Character. The path to the environment file containing
#'   variables required by the function. Default is ".env".
#' @param NCores Numeric. Number of cores to use for parallel processing.
#' @param Overwrite Logical. If `TRUE`, existing JPEG files will be overwritten
#'   during processing. See [IAS_Plot]
#' @name IAS_Process
#' @export
#' @author Ahmed El-Gabbas
#' @return The function returns `NULL` invisibly.
#' @note
#'   - The function should be used only after data from the three data sources
#' were processed: GBIF ([GBIF_Process]), EASIN ([EASIN_Process]), eLTER
#' ([elTER_Process]).
#'   - The function depends on the following functions: [IAS_Distribution] for
#' prepare species-specific final maps and [IAS_Plot] for plotting species
#' distribution as JPEG.

IAS_Process <- function(
    FromHPC = TRUE, EnvFile = ".env", NCores = 6, Overwrite = TRUE) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # Checking arguments ----
  IASDT.R::CatTime("Checking arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, ~ get(.x, envir = environment())) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "character", Args = "EnvFile")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "logical", Args = c("FromHPC", "Overwrite"))
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "numeric", Args = "NCores")

  rm(AllArgs)

  withr::local_options(future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE)

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_PA <- Path_TaxaInfo_RData <- taxon_name <- Path_TaxaStand <-
    Path_HabAff <- Species <- PA_Map <- NCells <- Species_name <- NCells_All <-
    synhab_name <- IAS_ID <- Count <- Hab <- EU_Bound <- Threshold <- NSp <-
    taxon_name <- PA_Masked_Map <- NCells_Naturalized <- NULL

  # # ..................................................................... ###

  # Environment variables ----
  IASDT.R::CatTime("Environment variables")

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_GBIF", "DP_R_GBIF", TRUE, FALSE,
      "Path_EASIN", "DP_R_EASIN", TRUE, FALSE,
      "Path_eLTER", "DP_R_eLTER_Out", FALSE, TRUE,
      "EU_Bound", "DP_R_EUBound_sf", FALSE, TRUE,
      "Path_PA", "DP_R_PA", FALSE, FALSE,
      "Path_HabAff", "DP_R_HabAff", FALSE, TRUE,
      "Path_TaxaStand", "DP_R_TaxaStand", FALSE, TRUE,
      "Path_TaxaInfo_RData", "DP_R_TaxaInfo_RData", FALSE, TRUE,

      # The following are needed for other called functions
      "Path_Grid", "DP_R_Grid", TRUE, FALSE,
      "Path_Grid_Ref", "DP_R_Grid_Ref", TRUE, FALSE,
      "Path_TaxaCNT", "DP_R_Taxa_Country", FALSE, TRUE,
      "Path_TaxaInfo", "DP_R_TaxaInfo", FALSE, TRUE,
      "Path_BioReg", "DP_R_BioReg", FALSE, TRUE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_GBIF", "DP_R_GBIF_Local", TRUE, FALSE,
      "Path_EASIN", "DP_R_EASIN_Local", TRUE, FALSE,
      "Path_eLTER", "DP_R_eLTER_Out_Local", FALSE, TRUE,
      "EU_Bound", "DP_R_EUBound_sf_Local", FALSE, TRUE,
      "Path_PA", "DP_R_PA_Local", FALSE, FALSE,
      "Path_HabAff", "DP_R_HabAff_Local", FALSE, TRUE,
      "Path_TaxaStand", "DP_R_TaxaStand_Local", FALSE, TRUE,
      "Path_TaxaInfo_RData", "DP_R_TaxaInfo_RData_Local", FALSE, TRUE,

      # The following are needed for other called functions
      "Path_Grid", "DP_R_Grid_Local", TRUE, FALSE,
      "Path_Grid_Ref", "DP_R_Grid_Ref_Local", TRUE, FALSE,
      "Path_TaxaCNT", "DP_R_Taxa_Country_Local", FALSE, TRUE,
      "Path_TaxaInfo", "DP_R_TaxaInfo_Local", FALSE, TRUE,
      "Path_BioReg", "DP_R_BioReg_Local", FALSE, TRUE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  rm(EnvVars2Read)

  # # ..................................................................... ###

  # Reading input data and check/create directories ------

  IASDT.R::CatTime("Reading input data and check/create directories")

  Path_PA_JPEG <- file.path(Path_PA, "JPEG_Maps")
  fs::dir_create(c(Path_PA_JPEG))

  # last update info
  LastUpdate <- stringr::str_glue(
    'Last update: {format(Sys.Date(), "%d %B %Y")}')

  ## Standardized taxonomy -----
  IASDT.R::CatTime("Standardized taxonomy", Level = 1)
  # list of original taxonomy (including dummy ID column `taxon_id_`)
  TaxaList_Original <- readRDS(Path_TaxaStand) %>%
    dplyr::select(tidyselect::all_of(c("taxon_id_", "taxon_name")))

  # # .................................... ###

  ## TaxaInfo ----
  IASDT.R::CatTime("TaxaInfo", Level = 1)
  TaxaList <- IASDT.R::LoadAs(Path_TaxaInfo_RData)
  TaxaList_Distinct <- dplyr::distinct(TaxaList, taxon_name, IAS_ID)

  # # .................................... ###

  ## Species habitat affinity -----
  IASDT.R::CatTime("Species habitat affinity", Level = 1)
  SynHab_List <- tibble::tribble(
    ~SynHab_code, ~SynHab_name,
    "1", "Forests",
    "2", "Open_forests",
    "3", "Scrub",
    "4", "Grasslands",
    "4a", "Natural_grasslands",
    "4b", "Human_maintained_grasslands",
    "5", "Sandy",
    "6", "Rocky",
    "7", "Dryland",
    "8", "Saline",
    "9", "Riparian",
    "10", "Wetland",
    "11", "Aquatic",
    "12", "Man_made",
    "12a", "Ruderal_habitats",
    "12b", "Agricultural_habitats")

  Sp_SynHab <- readRDS(Path_HabAff) %>%
    dplyr::rename(SynHab_name = synhab_name) %>%
    dplyr::left_join(SynHab_List, by = "SynHab_name") %>%
    dplyr::mutate(Count = TRUE) %>%
    tidyr::pivot_wider(
      names_from = c("SynHab_code", "SynHab_name"),
      names_prefix = "Hab_", values_from = Count) %>%
    # add taxon_name column
    dplyr::left_join(TaxaList_Original, by = "taxon_id_") %>%
    # add IAS_ID column
    dplyr::left_join(TaxaList_Distinct, by = "taxon_name") %>%
    dplyr::select(-tidyselect::all_of(c("taxon_id_", "taxon_name"))) %>%
    dplyr::select(IAS_ID, gtools::mixedsort(tidyselect::peek_vars())) %>%
    dplyr::distinct() %>%
    dplyr::arrange(IAS_ID)

  rm(TaxaList_Original, TaxaList_Distinct)

  # # ..................................................................... ###

  # Species-specific data ------

  IASDT.R::CatTime("Species-specific data")
  .StartTimeDist <- lubridate::now(tzone = "CET")

  ## Prepare working on parallel -----
  IASDT.R::CatTime(
    paste0("Prepare working on parallel using `", NCores, "` cores."),
    Level = 1)
  future::plan("multisession", workers = NCores, gc = TRUE)
  on.exit(future::plan("sequential"), add = TRUE)

  # # .................................... ###

  ## Species-specific data on parallel ----
  IASDT.R::CatTime("Species-specific data on parallel", Level = 1)

  Sp_PA_Data <- future.apply::future_lapply(
    X = sort(unique(TaxaList$Species_name)),
    FUN = function(x) {
      IASDT.R::IAS_Distribution(
        Species = x, FromHPC = FromHPC, EnvFile = EnvFile, Verbose = FALSE,
        Overwrite = Overwrite)
    },
    future.scheduling = Inf,
    future.seed = TRUE,
    future.packages =   c(
      "dplyr", "lubridate", "IASDT.R", "purrr", "stringr", "readr", "fs",
      "sf", "terra", "readxl", "tidyr", "tidyselect", "ggplot2", "ggtext",
      "grid", "tidyterra", "cowplot", "scales"),
    future.globals = c("FromHPC", "EnvFile", "Overwrite")) %>%
    dplyr::bind_rows()

  # # .................................... ###

  ## Stopping cluster ----
  IASDT.R::CatTime("Stopping cluster", Level = 1)
  future::plan("sequential")

  IASDT.R::CatDiff(
    InitTime = .StartTimeDist,
    Prefix = "Processing Species-specific data took ", NLines = 1, Level = 2)

  # # ..................................................................... ###

  # Merge Species-specific summary info -----
  IASDT.R::CatTime("Merge Species-specific summary info")

  Sp_PA_Data <- Sp_PA_Data %>%
    # Split number of grid cell per country / biogeographical region as separate
    # column
    tidyr::unnest_wider(c("BioRegs_DT", "BioRegsMask_DT", "SpCountry")) %>%
    dplyr::mutate(
      # replace NA for biogeographical regions and country analysis with 0
      dplyr::across(
        tidyselect::matches("^CNT_|^BioReg_"),
        ~ dplyr::if_else(is.na(.x), 0L, .x))) %>%
    # change column order
    dplyr::select(
      Species:NCells_Naturalized, PA_Map, PA_Masked_Map,
      tidyselect::matches("^GBIF"),
      tidyselect::matches("^EASIN"),
      tidyselect::matches("^eLTER"),
      tidyselect::matches("^iNat"),
      tidyselect::matches("^BioReg_"),
      tidyselect::matches("^BioRegSumm_"),
      tidyselect::matches("^BioReg_Masked_"),
      tidyselect::matches("^BioRegMaskSumm_"),
      tidyselect::matches("^CNT_"),
      tidyselect::everything()) %>%
    dplyr::rename(Species_name = Species) %>%
    dplyr::full_join(
      y = dplyr::distinct(
        dplyr::select(TaxaList, -tidyselect::all_of("speciesKey"))),
      by = "Species_name") %>%
    dplyr::select(
      tidyselect::all_of(
        c(
          "IAS_ID", "taxon_name", "Species_name",
          "Species_name2", "Species_File")),
      tidyselect::everything()
    ) %>%
    dplyr::left_join(Sp_SynHab, by = "IAS_ID")

  rm(TaxaList)

  ## Save summary results -----
  IASDT.R::CatTime("Save summary results", Level = 1)

  IASDT.R::CatTime("`RData`", Level = 2)
  save(Sp_PA_Data, file = file.path(Path_PA, "Sp_PA_Data.RData"))

  IASDT.R::CatTime("Summary data without maps", Level = 1)
  IASDT.R::CatTime("csv format", Level = 2)
  Sp_PA_Summary_DF <- Sp_PA_Data %>%
    dplyr::select(
      -tidyselect::all_of(
        c("GBIF_R", "EASIN_R", "eLTER_R", "PA_Map", "PA_Masked_Map")))

  readr::write_excel_csv(
    Sp_PA_Summary_DF, file = file.path(Path_PA, "Sp_PA_Summary_DF.csv"),
    progress = FALSE)

  IASDT.R::CatTime("RData format", Level = 2)
  save(Sp_PA_Summary_DF, file = file.path(Path_PA, "Sp_PA_Summary_DF.RData"))

  # # # ..................................................................... ###
  #
  # # Plotting species-specific maps ------
  #
  # IASDT.R::CatTime("Plotting species-specific maps")
  # .StartTimeDist <- lubridate::now(tzone = "CET")
  #
  # ## Prepare working on parallel -----
  # IASDT.R::CatTime(
  #   paste0("Prepare working on parallel using `", NCores, "` cores."),
  #   Level = 1)
  # future::plan("multisession", workers = NCores, gc = TRUE)
  # on.exit(future::plan("sequential"), add = TRUE)
  #
  # # # .................................... ###
  #
  # ## Plotting species-specific maps on parallel ----
  # IASDT.R::CatTime("Plotting species-specific maps on parallel", Level = 1)
  #
  # Sp_PA_Data <- future.apply::future_lapply(
  #   X = sort(unique(TaxaList$Species_name)),
  #   FUN = function(x) {
  #     IASDT.R::IAS_Plot(
  #       Species = x, FromHPC = FromHPC, EnvFile = EnvFile,
  #       Overwrite = Overwrite)
  #   },
  #   future.scheduling = Inf,
  #   future.seed = TRUE,
  #   future.packages = c(
  #     "dplyr", "lubridate", "IASDT.R", "purrr", "stringr", "readr", "fs",
  #     "sf", "terra", "readxl", "tidyr", "tidyselect", "ggplot2", "ggtext",
  #     "grid", "tidyterra", "cowplot", "scales"),
  #   future.globals = c("FromHPC", "EnvFile", "Overwrite")) %>%
  #   dplyr::bind_rows()
  #
  # # # .................................... ###
  #
  # ## Stopping cluster ----
  # IASDT.R::CatTime("Stopping cluster", Level = 1)
  # future::plan("sequential")
  #
  # IASDT.R::CatDiff(
  #   InitTime = .StartTimeDist,
  #   Prefix = "Processing Species-specific data took ", NLines = 1, Level = 2)

  # # ..................................................................... ###

  # Summary of merged data ------
  IASDT.R::CatTime("Summary of merged data")

  ## Number of IAS per grid cell -----
  IASDT.R::CatTime("# IAS per grid cell", Level = 1)

  IAS_NumSp <- dplyr::filter(Sp_PA_Data, NCells_All > 0) %>%
    dplyr::pull(PA_Map) %>%
    purrr::map(terra::unwrap) %>%
    terra::rast() %>%
    sum(na.rm = TRUE) %>%
    stats::setNames("IAS_NumSp") %>%
    # Ensure that values are read from memory
    IASDT.R::setRastVals()

  # save as RData
  IASDT.R::SaveAs(
    InObj = terra::wrap(IAS_NumSp), OutObj = "IAS_NumSp",
    OutPath = file.path(Path_PA, "IAS_NumSp.RData"))

  # save as tif
  raster::writeRaster(
    x = IAS_NumSp, overwrite = TRUE,
    filename = file.path(Path_PA, "IAS_NumSp.tif"))

  # # .................................... ###

  ## Number of IAS per grid cell - masked data -----

  IAS_NumSp_Masked <- dplyr::filter(Sp_PA_Data, NCells_Naturalized > 0) %>%
    dplyr::pull(PA_Masked_Map) %>%
    purrr::map(terra::unwrap) %>%
    terra::rast() %>%
    sum(na.rm = TRUE) %>%
    stats::setNames("IAS_NumSp_Masked") %>%
    # Ensure that values are read from memory
    IASDT.R::setRastVals()

  # save as RData
  IASDT.R::SaveAs(
    InObj = terra::wrap(IAS_NumSp_Masked), OutObj = "IAS_NumSp_Masked",
    OutPath = file.path(Path_PA, "IAS_NumSp_Masked.RData"))

  # save as tif
  raster::writeRaster(
    x = IAS_NumSp_Masked, overwrite = TRUE,
    filename = file.path(Path_PA, "IAS_NumSp_Masked.tif"))

  # # .................................... ###

  ## Plotting summary of IAS data -----
  IASDT.R::CatTime("Plotting summary of IAS data", Level = 1)

  EUBound <- IASDT.R::LoadAs(EU_Bound) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_03")
  MapLimX <- c(2600000, 6550000)
  MapLimY <- c(1450000, 5410000)

  # # +++++++++++++++++++++++++++++++++ ###

  PlottingTheme <- ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0.1, "cm"),
      plot.title = ggplot2::element_text(
        size = 10, color = "blue", face = "bold", hjust = 0,
        margin = ggplot2::margin(2, 0, 2, 0)),
      plot.subtitle = ggplot2::element_text(
        size = 8, color = "darkgrey", face = "italic", hjust = 0,
        margin = ggplot2::margin(1, 0, 1, 0)),
      strip.text = ggplot2::element_text(size = 6, face = "bold"),
      legend.key.size = grid::unit(0.6, "cm"),
      legend.key.width = grid::unit(0.5, "cm"),
      legend.background = ggplot2::element_rect(fill = "transparent"),
      legend.text = ggplot2::element_text(size = 6),
      legend.box.spacing = grid::unit(0, "pt"),
      legend.title = ggplot2::element_text(
        color = "blue", size = 6, face = "bold", hjust = 0),
      legend.position = "inside",
      legend.position.inside = c(0.92, 0.825),
      axis.text.x = ggplot2::element_text(size = 7),
      axis.text.y = ggplot2::element_text(size = 7, hjust = 0.5, angle = 90),
      axis.ticks = ggplot2::element_line(colour = "blue", linewidth = 0.25),
      axis.ticks.length = grid::unit(0.04, "cm"),
      panel.spacing = grid::unit(0.3, "lines"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        linewidth = 0.1, colour = "grey40", linetype = 2),
      panel.border = ggplot2::element_blank(),
      plot.tag.position = c(0.94, 0.011),
      plot.tag = ggtext::element_markdown(colour = "grey", size = 4),
      panel.ontop = TRUE, panel.background = ggplot2::element_rect(fill = NA))

  # # +++++++++++++++++++++++++++++++++ ###

  Plot_Nsp <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      EUBound, mapping = ggplot2::aes(), color = "grey30",
      linewidth = 0.1, fill = "grey95", inherit.aes = TRUE) +
    tidyterra::geom_spatraster(
      data = terra::classify(IAS_NumSp, cbind(0, NA)), maxcell = Inf) +
    paletteer::scale_fill_paletteer_c(
      na.value = "transparent", "viridis::plasma",
      breaks = IASDT.R::integer_breaks()) +
    ggplot2::geom_sf(
      EUBound, mapping = ggplot2::aes(), color = "grey40",
      linewidth = 0.075, fill = "transparent", inherit.aes = TRUE) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)), limits = MapLimX) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)), limits = MapLimY) +
    ggplot2::labs(
      title = "Number of IAS per 10\u00D710 km grid",
      subtitle = "All data (including cultivated or casual observations)",
      fill = "# IAS") +
    PlottingTheme

  Plot_Nsp_log <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      EUBound, mapping = ggplot2::aes(), color = "grey30",
      linewidth = 0.1, fill = "grey95", inherit.aes = TRUE) +
    tidyterra::geom_spatraster(
      data = log10(terra::classify(IAS_NumSp, cbind(0, NA))), maxcell = Inf) +
    paletteer::scale_fill_paletteer_c(
      na.value = "transparent", "viridis::plasma",
      breaks = IASDT.R::integer_breaks()) +
    ggplot2::geom_sf(
      EUBound, mapping = ggplot2::aes(), color = "grey40",
      linewidth = 0.075, fill = "transparent", inherit.aes = TRUE) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)), limits = MapLimX) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)), limits = MapLimY) +
    ggplot2::labs(
      title = "Number of IAS per 10\u00D710 km grid (log10 scale)",
      subtitle = "All data (including cultivated or casual observations)",
      fill = "# IAS\n(log10)") +
    PlottingTheme

  (cowplot::plot_grid(Plot_Nsp, Plot_Nsp_log, ncol = 2, nrow = 1) +
      cowplot::draw_label(
        label = LastUpdate, color = "grey", x = 0.99, y = 0.98,
        size = 7, hjust = 1)) %>%
    ggplot2::ggsave(
      filename = file.path(Path_PA, "IAS_NumSpecies.jpeg"),
      width = 30, height = 15.5, units = "cm", dpi = 600)

  # # +++++++++++++++++++++++++++++++++ ###

  Plot_Nsp_Masked <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      EUBound, mapping = ggplot2::aes(), color = "grey30",
      linewidth = 0.1, fill = "grey95", inherit.aes = TRUE) +
    tidyterra::geom_spatraster(
      data = terra::classify(IAS_NumSp_Masked, cbind(0, NA)), maxcell = Inf) +
    paletteer::scale_fill_paletteer_c(
      na.value = "transparent", "viridis::plasma",
      breaks = IASDT.R::integer_breaks()) +
    ggplot2::geom_sf(
      EUBound, mapping = ggplot2::aes(), color = "grey40",
      linewidth = 0.075, fill = "transparent", inherit.aes = TRUE) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)), limits = MapLimX) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)), limits = MapLimY) +
    ggplot2::labs(
      title = "Number of IAS per 10\u00D710 km grid",
      subtitle = "Excluding cultivated or casual observations",
      fill = "# IAS") +
    PlottingTheme

  Plot_Nsp_Masked_log <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      EUBound, mapping = ggplot2::aes(), color = "grey30",
      linewidth = 0.1, fill = "grey95", inherit.aes = TRUE) +
    tidyterra::geom_spatraster(
      data = log10(terra::classify(IAS_NumSp_Masked, cbind(0, NA))),
      maxcell = Inf) +
    paletteer::scale_fill_paletteer_c(
      na.value = "transparent", "viridis::plasma",
      breaks = IASDT.R::integer_breaks()) +
    ggplot2::geom_sf(
      EUBound, mapping = ggplot2::aes(), color = "grey40",
      linewidth = 0.075, fill = "transparent", inherit.aes = TRUE) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)), limits = MapLimX) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)), limits = MapLimY) +
    ggplot2::labs(
      title = "Number of IAS per 10\u00D710 km grid (log10 scale)",
      subtitle = "Excluding cultivated or casual observations",
      fill = "# IAS\n(log10)") +
    PlottingTheme

  (cowplot::plot_grid(
    Plot_Nsp_Masked, Plot_Nsp_Masked_log, ncol = 2, nrow = 1) +
      cowplot::draw_label(
        label = LastUpdate, color = "grey", x = 0.99, y = 0.98,
        size = 7, hjust = 1)) %>%
    ggplot2::ggsave(
      filename = file.path(Path_PA, "IAS_NumSpecies_Masked.jpeg"),
      width = 30, height = 15.5, units = "cm", dpi = 600)

  # # +++++++++++++++++++++++++++++++++ ###

  rm(
    IAS_NumSp, Plot_Nsp, Plot_Nsp_log, PlottingTheme, EUBound, MapLimX,
    MapLimY, Plot_Nsp_Masked_log, Plot_Nsp_Masked)

  # # .................................... ###

  ## # Number of IAS to be used in the model - per habitat type ------

  Threshold_Theme <- ggplot2::theme(
    plot.margin = ggplot2::margin(0.05, 0.15, 0.05, 0, "cm"),
    plot.title = ggplot2::element_text(
      size = 11, color = "black", hjust = 0.05, face = "italic",
      margin = ggplot2::margin(0.05, 0, 0.05, 0, "cm")),
    legend.position = "bottom",
    legend.text = ggplot2::element_text(
      size = 6, face = "bold", margin = ggplot2::margin(0, 0, 0, 0, "cm")),
    legend.box.spacing = grid::unit(-0.05, "cm"),
    legend.key.spacing.x = grid::unit(0.75, "cm"),
    legend.key = ggplot2::element_rect(fill = scales::alpha("white", 0)),
    legend.key.height = grid::unit(0.4, "cm"),
    legend.key.width = grid::unit(0.1, "cm"),
    axis.title.y.left = ggplot2::element_text(size = 9, vjust = -1),
    axis.title.x = ggplot2::element_text(size = 9, vjust = 1),
    axis.text.y.left = ggplot2::element_text(size = 7),
    axis.text.x.bottom = ggplot2::element_text(size = 7),
    axis.text.y.right = ggplot2::element_blank(),
    axis.text.x.top = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(
      linewidth = 0.1, colour = "grey40", linetype = 2),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank())

  # # +++++++++++++++++++++++++++++++++ ###

  NSp_DT <- dplyr::select(Sp_PA_Data, Species_name, NCells_All)
  NSp_DT <- purrr::map_dfr(
    .x = sort(unique(NSp_DT$NCells_All)),
    .f = ~ tibble::tibble(
      Threshold = .x, NSp = sum(NSp_DT$NCells_All >= .x, na.rm = TRUE))) %>%
    dplyr::filter(Threshold < 450)

  NSp_DT_Masked <- dplyr::select(Sp_PA_Data, Species_name, NCells_Naturalized)
  NSp_DT_Masked <- purrr::map_dfr(
    .x = sort(unique(NSp_DT_Masked$NCells_Naturalized)),
    .f = ~ tibble::tibble(
      Threshold = .x,
      NSp = sum(NSp_DT_Masked$NCells_Naturalized >= .x, na.rm = TRUE))) %>%
    dplyr::filter(Threshold < 450)

  # # +++++++++++++++++++++++++++++++++ ###

  NSp_Hab_DT <- purrr::map_dfr(
    .x = c(
      "Hab_1_Forests", "Hab_2_Open_forests", "Hab_3_Scrub",
      "Hab_4a_Natural_grasslands", "Hab_4b_Human_maintained_grasslands",
      "Hab_10_Wetland", "Hab_12a_Ruderal_habitats",
      "Hab_12b_Agricultural_habitats"),
    .f = function(Hab) {
      CurrDT <- dplyr::filter(Sp_PA_Data, !!as.symbol(Hab))
      NCells <- sort(unique(Sp_PA_Data$NCells_All))
      purrr::map_dfr(
        sort(unique(Sp_PA_Data$NCells_All)),
        ~ tibble::tibble(
          Hab = Hab, Threshold = .x,
          NSp = sum(CurrDT$NCells_All >= .x, na.rm = TRUE)))
    }
  ) %>%
    dplyr::mutate(
      Hab = stringr::str_remove(Hab, "^Hab_"),
      Hab = stringr::str_replace(Hab, "_", " - "),
      Hab = stringr::str_replace_all(Hab, "_", " "),
      Hab = forcats::fct_inorder(Hab)) %>%
    dplyr::slice(gtools::mixedorder(Hab)) %>%
    dplyr::filter(Threshold < 450) %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(
      ggplot2::aes(x = Threshold, y = NSp, group = Hab, colour = Hab),
      linewidth = 0.75) +
    ggplot2::geom_line(
      data = NSp_DT, inherit.aes = FALSE, colour = "black",
      ggplot2::aes(x = Threshold, y = NSp), linewidth = 0.75) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0, 0, 0), limits = c(0, NA), oob = scales::oob_keep) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0, 0, 0), limits = c(-5, NA), oob = scales::oob_keep) +
    ggplot2::xlab(
      "Threshold (# of 10\u00D710 km presence grid cells per species)") +
    ggplot2::ylab("Number of species") +
    ggplot2::scale_color_discrete(name = NULL) +
    ggplot2::labs(
      title = "All data (including cultivated or casual observations)") +
    Threshold_Theme +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        nrow = 1, label.position = "bottom",
        override.aes = list(linewidth = 1.5)))

  NSp_Hab_Masked_DT <- purrr::map_dfr(
    .x = c(
      "Hab_1_Forests", "Hab_2_Open_forests", "Hab_3_Scrub",
      "Hab_4a_Natural_grasslands", "Hab_4b_Human_maintained_grasslands",
      "Hab_10_Wetland", "Hab_12a_Ruderal_habitats",
      "Hab_12b_Agricultural_habitats"),
    .f = function(Hab) {
      CurrDT <- dplyr::filter(Sp_PA_Data, !!as.symbol(Hab))
      NCells <- sort(unique(Sp_PA_Data$NCells_Naturalized))
      purrr::map_dfr(
        sort(unique(Sp_PA_Data$NCells_Naturalized)),
        ~ tibble::tibble(
          Hab = Hab, Threshold = .x,
          NSp = sum(CurrDT$NCells_Naturalized >= .x, na.rm = TRUE)))
    }
  ) %>%
    dplyr::mutate(
      Hab = stringr::str_remove(Hab, "^Hab_"),
      Hab = stringr::str_replace(Hab, "_", " - "),
      Hab = stringr::str_replace_all(Hab, "_", " "),
      Hab = forcats::fct_inorder(Hab)) %>%
    dplyr::slice(gtools::mixedorder(Hab)) %>%
    dplyr::filter(Threshold < 450) %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(
      ggplot2::aes(x = Threshold, y = NSp, group = Hab, colour = Hab),
      linewidth = 0.75) +
    ggplot2::geom_line(
      data = NSp_DT_Masked, inherit.aes = FALSE, colour = "black",
      ggplot2::aes(x = Threshold, y = NSp), linewidth = 0.75) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0, 0, 0), limits = c(0, NA), oob = scales::oob_keep) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0, 0, 0), limits = c(-5, NA), oob = scales::oob_keep) +
    ggplot2::xlab(
      "Threshold (# of 10\u00D710 km presence grid cells per species)") +
    ggplot2::ylab("Number of species") +
    ggplot2::scale_color_discrete(name = NULL) +
    ggplot2::labs(title = "Excluding cultivated or casual observations") +
    Threshold_Theme +
    ggplot2::theme(legend.position = "none")

  Legend <- ggpubr::as_ggplot(ggpubr::get_legend(NSp_Hab_DT))
  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste0(
        "Number of IAS to be used in the models based on the arbitrary ",
        "selection of # of presence grid cells per species"),
      fontface = "bold", colour = "blue")

  cowplot::plot_grid(
    title,
    cowplot::plot_grid(
      (NSp_Hab_DT + ggplot2::theme(legend.position = "none")),
      NSp_Hab_Masked_DT,
      ncol = 2, nrow = 1),
    Legend,
    ncol = 1, rel_heights = c(0.05, 1, 0.05)
  ) %>%
    ggplot2::ggsave(
      filename = file.path(Path_PA, "IAS_NSp_threshold_Hab.jpeg"),
      width = 30, height = 17, units = "cm", dpi = 600)

  # # ..................................................................... ###

  # Function Summary ----
  IASDT.R::CatDiff(
    InitTime = .StartTime,
    Prefix = "\nProcessing species data was finished in ", ... = "\n")

  return(invisible(NULL))
}
