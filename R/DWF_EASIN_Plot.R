## |------------------------------------------------------------------------| #
# EASIN_Plot ----
## |------------------------------------------------------------------------| #

#' EASIN_Plot
#'
#' Generates a series of plots for cleaned EASIN data, including the number of
#' observations, number of species, and the distribution of these metrics per
#' data partner. The function requires environment variables for paths to data
#' and expects specific files to be present at these paths.
#'
#' @param EnvFile Character. The path to the `.env` file containing environment
#'   variables. Defaults to ".env".
#' @param FromHPC Logical indicating whether the work is being done from HPC, to
#'   adjust file paths accordingly. Default: `TRUE`.
#' @return The function returns NULL invisibly. Plots are saved as JPEG files in
#'   the specified directory.
#' @details The function first loads environment variables either from a
#'   high-performance computing (HPC) environment or a local environment. Then,
#'   it loads required data files, checks for missing files, and proceeds to
#'   generate and save plots as JPEG files. The plots include the number of
#'   observations and species and their distribution per data partner.
#' @note This function is not intended to be used directly by the user or in the
#'   IAS-pDT, but only used inside the [EASIN_Process] function.
#' @name EASIN_Plot
#' @author Ahmed El-Gabbas
#' @export

EASIN_Plot <- function(EnvFile = ".env", FromHPC = TRUE) {

  .PlotStartTime <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Checking arguments ----
  IASDT.R::CatTime("Checking arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, ~get(.x, envir = environment())) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "character", Args = "EnvFile")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "logical", Args = "FromHPC")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_EASIN_Summary <- EU_Bound <- NULL

  # # |||||||||||||||||||||||||||||||||||
  # # Environment variables ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Environment variables", Level = 1)
  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "EU_Bound", "DP_R_EUBound_sf", FALSE, TRUE,
      "Path_EASIN", "DP_R_EASIN", TRUE, FALSE,
      "Path_EASIN_Interim", "DP_R_EASIN_Interim", TRUE, FALSE,
      "Path_EASIN_Summary", "DP_R_EASIN_Summary", TRUE, FALSE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "EU_Bound", "DP_R_EUBound_sf_Local", FALSE, TRUE,
      "Path_EASIN", "DP_R_EASIN_Local", TRUE, FALSE,
      "Path_EASIN_Interim", "DP_R_EASIN_Interim_Local", TRUE, FALSE,
      "Path_EASIN_Summary", "DP_R_EASIN_Summary_Local", TRUE, FALSE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  # # |||||||||||||||||||||||||||||||||||
  # # Input maps ----
  # # |||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Loading input maps", Level = 1)

  ## Country boundaries ----
  IASDT.R::CatTime("Country boundaries", Level = 2)
  EuroBound <- IASDT.R::LoadAs(EU_Bound) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_03")

  ## Check input summary maps -----
  IASDT.R::CatTime("Check input summary maps", Level = 2)
  Path_NSp <- file.path(Path_EASIN_Summary, "EASIN_NSp.RData")
  Path_NSp_PerPartner <- file.path(
    Path_EASIN_Summary, "EASIN_NSp_PerPartner.RData")

  Path_NObs <- file.path(Path_EASIN_Summary, "EASIN_NObs.RData")
  Path_NObs_PerPartner <- file.path(
    Path_EASIN_Summary, "EASIN_NObs_PerPartner.RData")

  PathSummaryMaps <- c(
    Path_NSp, Path_NSp_PerPartner, Path_NObs, Path_NObs_PerPartner)
  SummaryMapsMissing <- !file.exists(PathSummaryMaps)

  if (any(SummaryMapsMissing)) {
    stop(
      paste0(
        "The following input files are missing: \n",
        paste0(
          " >> ", PathSummaryMaps[which(SummaryMapsMissing)],collapse = "\n")),
      call. = FALSE)
  }

  # # |||||||||||||||||||||||||||||||||||
  # # Plotting ----
  # # |||||||||||||||||||||||||||||||||||

  ## NObs + NSp ----

  IASDT.R::CatTime("Number of observations and species", Level = 1)

  Plot_EASIN_All <- function(
    MapPath, Title, EuroBound, addTag = FALSE, Legend = FALSE) {

    LastUpdate <- paste0(
      "<b>Last update:</b> <i>", format(Sys.Date(), "%d %B %Y"), "</i>")

    Map <- IASDT.R::LoadAs(MapPath) %>%
      terra::unwrap() %>%
      log10()

    NCells <- terra::ncell(Map)

    PlottingTheme <- ggplot2::theme_bw() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
        plot.title = ggtext::element_markdown(
          size = 8, color = "blue", face = "bold", hjust = 0.5,
          margin = ggplot2::margin(0, 0, 0, 0)),
        legend.key.size = grid::unit(0.3, "cm"),
        legend.key.width = grid::unit(0.3, "cm"),
        legend.background = ggplot2::element_rect(fill = "transparent"),
        legend.text = ggplot2::element_text(size = 4),
        legend.box.spacing = grid::unit(0, "pt"),
        legend.title = ggtext::element_markdown(
          color = "blue", size = 6, face = "bold"),
        legend.position	= "inside",
        legend.position.inside = c(0.92, 0.85),
        axis.text.x = ggplot2::element_text(size = 4),
        axis.text.y = ggplot2::element_text(size = 4, hjust = 0.5, angle = 90),
        axis.ticks = ggplot2::element_line(colour = "blue", linewidth = 0.25),
        axis.ticks.length = grid::unit(0.04, "cm"),
        panel.spacing = grid::unit(0.3, "lines"),
        panel.grid.minor = ggplot2::element_line(linewidth = 0.125),
        panel.grid.major = ggplot2::element_line(linewidth = 0.25),
        panel.border = ggplot2::element_blank(),
        plot.tag.position = c(0.88, 0.999),
        plot.tag = ggtext::element_markdown(colour = "grey", size = 5))

    Plot <- ggplot2::ggplot() +
      tidyterra::geom_spatraster(data = Map, maxcell = NCells) +
      paletteer::scale_fill_paletteer_c(
        na.value = "transparent", "viridis::plasma") +
      ggplot2::geom_sf(
        EuroBound, mapping = ggplot2::aes(), color = "grey30",
        linewidth = 0.04, fill = scales::alpha("grey80", 0.2),
        inherit.aes = TRUE) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0, 0)),
        limits = c(2600000, 6700000)) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0)),
        limits = c(1450000, 5420000)) +
      PlottingTheme

    if (addTag) {
      Plot <- Plot +
        ggplot2::labs(
          title = Title, fill = "log<sub>10</sub>", tag = LastUpdate)
    } else {
      Plot <- Plot +
        ggplot2::labs(title = Title, fill = "log<sub>10</sub>")
    }
    return(Plot)
  }

  ### Number of observations ----
  IASDT.R::CatTime("Number of observations", Level = 2)
  Plot_NObs <- Plot_EASIN_All(
    MapPath = Path_NObs, Title = "Number of observations",
    EuroBound = EuroBound, addTag = FALSE, Legend = FALSE)

  ### Number of species ----
  IASDT.R::CatTime("Number of species", Level = 2)
  Plot_NSp <- Plot_EASIN_All(
    MapPath = Path_NSp, Title = "Number of species",
    EuroBound = EuroBound, addTag = TRUE, Legend = TRUE)

  ### Combine maps ----
  IASDT.R::CatTime(
    Text = "Merge maps side by side and save as JPEG", Level = 2)
  (ggpubr::ggarrange(
    Plot_NObs, (ggplot2::ggplot() + ggplot2::theme_void()), Plot_NSp,
    widths = c(1, 0, 1), nrow = 1) +
      patchwork::plot_annotation(
        title = "EASIN data",
        theme = ggplot2::theme(
          plot.margin = ggplot2::margin(0.1, 0, 0, 0, "cm"),
          plot.title = ggtext::element_markdown(
            size = 9, face = "bold", hjust = 0.5,
            margin = ggplot2::margin(0, 0, 0, 0))))) %>%
    ggplot2::ggsave(
      filename = file.path(Path_EASIN_Summary, "EASIN_Data.jpeg"),
      width = 20, height = 10.3, units = "cm", dpi = 600)

  rm(Plot_NSp, Plot_NObs)

  ## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ## Number of species/observations per partner ----

  IASDT.R::CatTime("Number of species/observations per partner", Level = 1)

  Plot_EASIN_Partner <- function(MapPath, File_prefix, Title) {

    LastUpdate <- paste0(
      "<b>Last update:</b> <i>", format(Sys.Date(), "%d %B %Y"), "</i>")

    PlottingTheme2 <-  ggplot2::theme_bw() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(0.125, 0, 0, 0, "cm"),
        plot.title = ggtext::element_markdown(
          size = 14, color = "blue", face = "bold", hjust = 0.5,
          margin = ggplot2::margin(2, 0, 4, 0)),
        panel.spacing = grid::unit(0.15, "lines"),
        panel.grid.minor = ggplot2::element_line(linewidth = 0.125),
        panel.grid.major = ggplot2::element_line(linewidth = 0.25),
        panel.border = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(size = 10, face = "bold"),
        axis.text.x = ggplot2::element_text(size = 4),
        axis.text.y = ggplot2::element_text(size = 4, hjust = 0.5, angle = 90),
        axis.ticks = ggplot2::element_line(colour = "blue", linewidth = 0.25),
        axis.ticks.length = grid::unit(0.04, "cm"),
        legend.text = ggplot2::element_text(size = 6),
        legend.title = ggtext::element_markdown(
          color = "blue", size = 9, face = "bold"),
        legend.position	= "inside",
        legend.position.inside = c(0.97, 0.90),
        legend.key.size = grid::unit(0.35, "cm"),
        legend.key.width = grid::unit(0.4, "cm"),
        legend.background = ggplot2::element_rect(fill = "transparent"),
        legend.box.spacing = grid::unit(0, "pt"),
        plot.tag.position = c(0.92, 0.975),
        plot.tag = ggtext::element_markdown(colour = "grey", size = 9))

    Map <- IASDT.R::LoadAs(MapPath) %>%
      terra::unwrap() %>%
      log10()
    NCells <- terra::ncell(Map)

    LegLimit <- c(
      min(terra::global(Map, min, na.rm = TRUE), na.rm = TRUE),
      max(terra::global(Map, max, na.rm = TRUE), na.rm = TRUE))

    for (i in seq_len(2)) {
      Start_Lyr <- (i - 1) * 8 + 1
      End_Lyr <- min(i * 8, terra::nlyr(Map))

      plot <- ggplot2::ggplot() +
        tidyterra::geom_spatraster(
          data = Map[[Start_Lyr:End_Lyr]], maxcell = NCells) +
        ggplot2::facet_wrap(~lyr, nrow = 2, ncol = 4) +
        paletteer::scale_fill_paletteer_c(
          na.value = "transparent", "viridis::plasma", limits = LegLimit) +
        ggplot2::geom_sf(
          EuroBound, mapping = ggplot2::aes(), color = "grey30",
          linewidth = 0.04, fill = scales::alpha("grey80", 0.2),
          inherit.aes = TRUE) +
        ggplot2::scale_x_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(2600000, 6700000)) +
        ggplot2::scale_y_continuous(
          expand = ggplot2::expansion(mult = c(0, 0)),
          limits = c(1450000, 5420000)) +
        ggplot2::labs(
          title = paste(Title, "[p", i, "]", sep = ""),
          fill = "log<sub>10</sub>", tag = LastUpdate) +
        PlottingTheme2

      ggplot2::ggsave(
        filename = file.path(
          Path_EASIN_Summary, paste0(File_prefix, "_p", i, ".jpeg")),
        width = 30, height = 16.5, units = "cm", dpi = 600)
    }

    return(invisible(NULL))
  }


  ### Number of observations per partner ----
  IASDT.R::CatTime("Number of observations per partner", Level = 2)
  Plot_EASIN_Partner(
    MapPath = Path_NObs_PerPartner,
    File_prefix = "EASIN_NObs_per_partner",
    Title = "EASIN data - Number of observations per data partner ")


  ### Number of species per partner ----
  IASDT.R::CatTime("Number of species per partner", Level = 2)
  Plot_EASIN_Partner(
    MapPath = Path_NSp_PerPartner,
    File_prefix = "EASIN_NSp_per_partner",
    Title = "EASIN data - Number of species per data partner ")

  ## |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatDiff(
    InitTime = .PlotStartTime, 
    Prefix = "Plotting EASIN data was finished in ", Level = 1)

  return(invisible(NULL))
}
