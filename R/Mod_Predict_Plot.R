## |------------------------------------------------------------------------| #
# Mod_Predict_Plot ----
## |------------------------------------------------------------------------| #

#' Plot species and level of invasion predictions as JPEG files using `ggplot2`
#'
#' Generate predictions for species and habitat models and saves the output as
#' JPEG files.
#'
#' @param ModelDir Path to the model directory containing predictions.
#' @param EnvFile Path to the environment file (`.env`) for setting paths.
#' @param FromHPC Boolean indicating whether the environment is an HPC system.
#' @return Saves prediction plots as JPEG files in the specified output
#'   directory.
#' @name Predict_Hmsc
#' @author Ahmed El-Gabbas
#' @export

Mod_Predict_Plot <- function(
  ModelDir, EnvFile = ".env", FromHPC = TRUE, NCores = 8) {

  .StartTime <- lubridate::now(tzone = "CET")

  # Set null device for ragg. This is to properly render the plots using
  # ggtext::geom_richtext
  cowplot::set_null_device("cairo")

  tif_path_mean <- tif_path_sd <- tif_path_cov <- Path_CLC <- SpeciesID <-
    IAS_ID <- Species_File <- NULL

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Assign environment variables ----

  if (!file.exists(EnvFile)) {
    stop(paste0(
      "Path to environment variables: ", EnvFile, " was not found"),
      call. = FALSE)
  }

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_CLC", "DP_R_CLC", TRUE, FALSE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_CLC", "DP_R_CLC_Local", TRUE, FALSE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  rm(EnvVars2Read, envir = environment())

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Plotting function ----

  # helper function for generating ggplot objects
  # - Map: A SpatRaster object for plotting
  # - Title: Title of the plot
  # - Observed: Whether the plot is for observed data (default: FALSE)
  # - LegendTitle: Title of the legend (default: NULL)
  # - ShowLegend: Whether to display the legend (default: FALSE)
  # - breaks, limits: Optional parameters for legend customization (default:
  # - NULL)

  PrepPlots <- function(
    Map, Title = NULL, Observed = FALSE, LegendTitle = NULL,
    ShowLegend = FALSE, breaks = NULL, limits = NULL) {

    Xlim <- c(2600000, 6500000)
    Ylim <- c(1450000, 5420000)

    Plot <- ggplot2::ggplot() +
      tidyterra::geom_spatraster(
        data = Map, maxcell = Inf, show.legend = ShowLegend)

    if (is.null(breaks)) {
      LegendBreaks <- Legendlabels <- ggplot2::waiver()
    } else {
      LegendBreaks <- Legendlabels <- breaks
    }

    if (is.null(limits)) {
      PlotLimits <- unlist(terra::global(Map, range, na.rm = TRUE))
    } else {
      PlotLimits <- limits
    }

    if (Observed) {
      Plot <- Plot +
        ggplot2::scale_fill_manual(
          breaks = c(0, 1, 3), values = c("grey90", "blue", "red"),
          labels = c("Not observed", "Present", "Excluded"),
          na.value = "transparent", name = NULL)
    } else {
      Plot <- Plot +
        paletteer::scale_fill_paletteer_c(
          na.value = "transparent", palette = "viridis::plasma",
          limits = PlotLimits, breaks = LegendBreaks, name = LegendTitle,
          labels = Legendlabels)
    }

    Plot <- Plot +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0, 0)), limits = Xlim,
        oob = scales::oob_keep) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0)), limits = Ylim) +
      ggplot2::labs(title = Title) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.margin = ggplot2::margin(0, 0.1, 0, 0.1, "cm"),
        plot.title = ggplot2::element_text(
          size = 9, color = "grey60", face = "bold", hjust = 0.5,
          margin = ggplot2::margin(0, 0, 0, 0)),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_blank(),
        legend.key.size = grid::unit(0.4, "cm"),
        legend.key.width = grid::unit(0.35, "cm"),
        legend.position = dplyr::if_else(ShowLegend, "inside", "none"),
        legend.position.inside = c(0.875, 0.775),
        legend.background = ggplot2::element_rect(fill = "transparent"),
        legend.text = ggplot2::element_text(size = 5),
        legend.box.spacing = grid::unit(0, "pt"),
        legend.title = ggplot2::element_text(
          color = "blue", size = 6, face = "bold", hjust = 0, vjust = 0),
        axis.text.x = ggplot2::element_text(size = 5),
        axis.text.y = ggplot2::element_text(size = 5, hjust = 0.5, angle = 90),
        axis.ticks = ggplot2::element_line(colour = "blue", linewidth = 0.25),
        axis.ticks.length = grid::unit(0.04, "cm"),
        axis.title = ggplot2::element_blank(),
        panel.spacing = grid::unit(0.3, "lines"),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(
          linewidth = 0.05, colour = "grey40", linetype = 2),
        panel.border = ggplot2::element_blank(),
        panel.ontop = TRUE,
        panel.background = ggplot2::element_rect(fill = NA))

    return(Plot)
  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Load summary of prediction maps ----

  IASDT.R::CatTime("Load summary of prediction maps")

  # Without clamping
  Map_summary_NoClamp <- file.path(
    ModelDir,
    "Model_Prediction/NoClamp/Prediction_Current_Summary.RData") %>%
    IASDT.R::LoadAs() %>%
    dplyr::rename(
      tif_path_mean_no_clamp = tif_path_mean,
      tif_path_sd_no_clamp = tif_path_sd,
      tif_path_cov_no_clamp = tif_path_cov)

  # With clamping
  Map_summary_Clamp <- file.path(
    ModelDir,
    "Model_Prediction/Clamp/Prediction_Current_Summary.RData") %>%
    IASDT.R::LoadAs() %>%
    dplyr::rename(
      tif_path_mean_clamp = tif_path_mean,
      tif_path_sd_clamp = tif_path_sd,
      tif_path_cov_clamp = tif_path_cov)

  # combined data
  Map_summary <- dplyr::full_join(
    Map_summary_NoClamp, Map_summary_Clamp,
    by = c(
      "hab_abb", "hab_name", "time_period", "climate_model",
      "climate_scenario", "ias_id", "taxon_name", "species_name", "class",
      "order", "family")) %>%
    dplyr::select(-c(
      "hab_abb", "hab_name", "time_period", "climate_model",
      "climate_scenario", "taxon_name"))

  HabAbb <- Map_summary_Clamp$hab_abb[[1]]
  Hab_Name <- Map_summary_Clamp$hab_name[[1]]
  AllSpID <- stringr::str_subset(Map_summary$ias_id, "^Sp_") %>%
    stringr::str_remove("^Sp_")
  rm(Map_summary_NoClamp, Map_summary_Clamp, envir = environment())

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Calculate observed species richness ----

  IASDT.R::CatTime("Calculate observed species richness")

  R_SR <- "datasets/processed/IAS_PA/Sp_PA_Data.RData" %>%
    IASDT.R::LoadAs() %>%
    dplyr::filter(SpeciesID %in% AllSpID) %>%
    dplyr::pull("PA_Masked_Map") %>%
    purrr::map(terra::unwrap) %>%
    terra::rast() %>%
    terra::app(sum, na.rm = TRUE) %>%
    terra::wrap()

  invisible(gc())

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Load habitat map ----

  IASDT.R::CatTime("Load habitat map")

  Path_Hab <- file.path(Path_CLC, "Summary_RData", "PercCov_SynHab_Crop.RData")
  if (!file.exists(Path_Hab)) {
    stop(
      paste0("Path_Hab file: ", Path_Hab, " does not exist"), call. = FALSE)
  }
  R_habitat <- IASDT.R::LoadAs(Path_Hab) %>%
    terra::unwrap() %>%
    terra::subset(paste0("SynHab_", HabAbb)) %>%
    terra::wrap()

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  Path_Plots <- file.path(ModelDir, "Model_Prediction", "Plots")
  fs::dir_create(Path_Plots)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Plotting ----

  IASDT.R::CatTime("Plotting")

  IASDT.R::CatTime(
    paste0("Preparing working on parallel using ", NCores, " cores"),
    Level = 1)
  c1 <- parallel::makePSOCKcluster(NCores)
  on.exit(try(parallel::stopCluster(c1), silent = TRUE), add = TRUE)

  IASDT.R::CatTime("Exporting variables to parallel cores", Level = 1)
  parallel::clusterExport(
    cl = c1,
    varlist = c("Map_summary", "PrepPlots", "R_SR", "R_habitat"),
    envir = environment())

  IASDT.R::CatTime("Loading packages at parallel cores", Level = 1)
  invisible(parallel::clusterEvalQ(
    cl = c1,
    expr = {
      sapply(
        c("dplyr", "terra", "ggplot2", "stringr", "cowplot", "tidyterra",
          "purrr", "ggtext", "ragg", "paletteer", "grid"),
        library, character.only = TRUE)
    }))

  IASDT.R::CatTime("Prepare and save plots on parallel", Level = 1)
  Plots <- parallel::parLapplyLB(
    cl = c1,
    X = seq_len(nrow(Map_summary)),
    fun = function(ID) {

      SpID <- Map_summary$ias_id[[ID]]
      SpName <- Map_summary$species_name[[ID]]
      ClassName <- Map_summary$class[[ID]]
      OrderName <- Map_summary$order[[ID]]
      FamilyName <- Map_summary$family[[ID]]
      Rank <- dplyr::if_else(SpID == "SR", "SR", "Species")
      Species <- stringr::str_detect(SpID, "^Sp")

      R_mean_NoClamp <- terra::rast(Map_summary$tif_path_mean_no_clamp[[ID]])
      R_mean_Clamp <- terra::rast(Map_summary$tif_path_mean_clamp[[ID]])
      R_sd_NoClamp <- terra::rast(Map_summary$tif_path_sd_no_clamp[[ID]])
      R_sd_Clamp <- terra::rast(Map_summary$tif_path_sd_clamp[[ID]])
      R_cov_NoClamp <- Map_summary$tif_path_cov_no_clamp[[ID]] %>%
        terra::rast() %>%
        "+"(0.001) %>%
        log10()
      R_cov_Clamp <- Map_summary$tif_path_cov_clamp[[ID]] %>%
        terra::rast() %>%
        "+"(0.001) %>%
        log10()

      if (Species) {

        Range_Mean <- c(R_mean_Clamp, R_mean_NoClamp) %>%
          terra::global(max, na.rm = TRUE) %>%
          max() %>%
          c(0, .)
        Breaks_Mean <- NULL

        SpID2 <- stringr::str_remove(SpID, "^Sp_")
        Path_JPEG <- file.path(
          Path_Plots,
          paste0("Pred_Current_Sp", SpID2, "_", SpName, ".jpeg"))
        SpID2 <- as.integer(SpID2)

      } else {

        Path_JPEG <- file.path(Path_Plots, "Pred_Current_SR.jpeg")

        Range_Mean <- c(terra::unwrap(R_SR), R_mean_NoClamp) %>%
          terra::global(max, na.rm = TRUE) %>%
          max() %>%
          c(0, .)
        Breaks_Mean <- NULL
      }

      Range_sd <- c(R_sd_NoClamp, R_sd_Clamp) %>%
        terra::global(range, na.rm = TRUE) %>%
        range()
      Range_cov <- c(R_cov_NoClamp, R_cov_Clamp) %>%
        terra::global(range, na.rm = TRUE) %>%
        range()

      # mean_no_clamp
      Plot_mean_no_clamp <- PrepPlots(
        Map = R_mean_NoClamp, Title = "Mean", ShowLegend = TRUE,
        breaks = Breaks_Mean, limits = Range_Mean)
      # mean_clamp
      Plot_mean_clamp <- PrepPlots(
        Map = R_mean_Clamp, breaks = Breaks_Mean, limits = Range_Mean)

      # sd_no_clamp
      Plot_sd_no_clamp <- PrepPlots(
        Map = R_sd_NoClamp, Title = "Standard deviation", ShowLegend = TRUE,
        breaks = NULL, limits = Range_sd)
      # sd_clamp
      Plot_sd_clamp <- PrepPlots(
        Map = R_sd_Clamp, breaks = NULL, limits = Range_sd)

      # cov_no_clamp
      Plot_cov_no_clamp <- PrepPlots(
        Map = R_cov_NoClamp, Title = "Coefficient of variation",
        ShowLegend = TRUE, breaks = NULL, limits = Range_cov,
        LegendTitle = "log10")
      # cov_clamp
      Plot_cov_clamp <- PrepPlots(
        Map = R_cov_Clamp, breaks = NULL, limits = Range_cov,
        LegendTitle = "log10")

      # Observed data
      if (Species) {

        # Observed species presence

        # Files containing observed data maps
        Path_observed <- "datasets/processed/IAS_PA/Sp_PA_Summary_DF.RData" %>%
          IASDT.R::LoadAs() %>%
          dplyr::filter(IAS_ID == SpID2) %>%
          dplyr::pull(Species_File) %>%
          paste0(., c("_Masked.tif", "_All.tif")) %>%
          file.path("datasets/processed/IAS_PA/tif", .)

        # Check if observed data files exist
        if (!all(file.exists(Path_observed))) {
          stop(
            paste0("Observed data for species: ", SpName, " not found"),
            call. = FALSE)
        }

        Plot_observed <- terra::rast(Path_observed) %>%
          terra::app(fun = function(vals) {
            ifelse(vals[1] == 0 & vals[2] == 1, 3, vals[1])
          }) %>%
          as.factor() %>%
          PrepPlots(
            Title = "Species observations", Observed = TRUE,
            ShowLegend = TRUE) +
          ggplot2::theme(
            panel.border = ggplot2::element_rect(
              color = "black", linewidth = 0.25, fill = NA,
              linetype = "dashed"))

      } else {

        # Observed species richness
        Plot_observed <- R_SR %>%
          terra::unwrap() %>%
          PrepPlots(
            Title = "Observed species richness", Observed = FALSE,
            ShowLegend = TRUE, breaks = Breaks_Mean, limits = Range_Mean) +
          ggplot2::theme(
            panel.border = ggplot2::element_rect(
              color = "black", linewidth = 0.25, fill = NA,
              linetype = "dashed"))
      }

      Plot_habitat <- PrepPlots(
        Map = terra::unwrap(R_habitat), breaks = seq(0, 100, 20),
        limits = c(0, 100), ShowLegend = TRUE, Title = "% Habitat cover",
        LegendTitle = "%") +
        ggplot2::theme(
          panel.border = ggplot2::element_rect(
            color = "black", linewidth = 0.25, fill = NA, linetype = "dashed"))

      plot_grid_main <- cowplot::plot_grid(
        Plot_mean_no_clamp, Plot_mean_clamp, Plot_sd_no_clamp, Plot_sd_clamp,
        Plot_cov_no_clamp, Plot_cov_clamp, Plot_observed, Plot_habitat,
        ncol = 4, nrow = 2, byrow = FALSE)

      YLab <- cowplot::ggdraw() +
        ggtext::geom_richtext(
          ggplot2::aes(
            x = 0.5, y = 0.725,
            label = stringr::str_glue(
              '<SPAN STYLE="font-size:11pt; color: darkgrey"><b>Without \\
              clamping</b></SPAN>')),
          fill = NA, label.color = NA, hjust = 0.5, vjust = 0.25,
          angle = 90, color = "black") +
        ggtext::geom_richtext(
          ggplot2::aes(
            x = 0.5, y = 0.275,
            label = stringr::str_glue(
              '<SPAN STYLE="font-size:11pt; color: darkgrey"><b>With \\
              clamping</b></SPAN>')),
          fill = NA, label.color = NA, hjust = 0.5, vjust = 0.25,
          angle = 90, color = "black") +
        ggtext::geom_richtext(
          ggplot2::aes(
            x = 1.1, y = 0.275,
            label = stringr::str_glue(
              '<span style="color: grey; font-size:8pt"> \\
              (efforts predictor fixed at 90% quantile)</span>')
          ),
          fill = NA, label.color = NA, hjust = 0.5, vjust = 0,
          angle = 90, color = "black") +
        ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))

      if (Species) {
        PlotTitle1 <- stringr::str_glue(
          '<SPAN STYLE="font-size:16pt">Predicted habitat suitability of \\
          <b><i>{SpName}</i></b></SPAN>')

        PlotTitle2 <- stringr::str_glue(
          '<SPAN STYLE="font-size:8pt; color: darkgrey"><b>Class:</b> \\
          {ClassName}; <b>Order:</b> {OrderName}; <b>Family:</b> \\
          {FamilyName}; <b>IAS_ID:</b> {SpID2}</SPAN>')
      } else {
        PlotTitle1 <- "Predicted level of invasion (mean species richness)"
        PlotTitle2 <- ""
      }

      Date <- format(Sys.Date(), "%d %B %Y")
      MainTitle <- cowplot::ggdraw() +
        ggtext::geom_richtext(
          ggplot2::aes(x = 0.01, y = 0.6, label = PlotTitle1),
          fill = NA, label.color = NA, hjust = 0, vjust = 0.5, size = 5,
          color = "black") +
        ggtext::geom_richtext(
          ggplot2::aes(x = 0.01, y = 0.2, label = PlotTitle2),
          fill = NA, label.color = NA, hjust = 0, vjust = 0.5, size = 5,
          color = "black") +
        ggtext::geom_richtext(
          ggplot2::aes(
            x = 1, y = 0.55,
            label = stringr::str_glue(
              '<SPAN STYLE="font-size:12.5pt; color: red"><b>{Hab_Name} \\
              habitat</b></SPAN> &#8212; <SPAN STYLE="font-size:12.5pt; \\
              color: blue"><b>Current climate</b></SPAN>')),
          fill = NA, label.color = NA, hjust = 1, vjust = 0.5) +
        ggtext::geom_richtext(
          ggplot2::aes(
            x = 1, y = 0.25,
            label = stringr::str_glue(
              '<span style="color: grey; font-size:6pt">Last update: {Date} \\
              </span>')),
          fill = NA, label.color = NA, hjust = 1, vjust = 0.5) +
        ggplot2::theme_void() +
        ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))

      final_plot <- cowplot::plot_grid(
        MainTitle,
        cowplot::plot_grid(
          YLab, plot_grid_main, ncol = 2, rel_widths = c(0.03, 1), align = "h"),
        ncol = 1, rel_heights = c(0.07, 1))

      ragg::agg_jpeg(
        filename = Path_JPEG, width = 30, height = 15.5, res = 600,
        quality = 100, units = "cm")
      print(final_plot)
      grDevices::dev.off()

      return(invisible(NULL))
    })

  rm(Plots, envir = environment())

  IASDT.R::CatDiff(InitTime = .StartTime)

  return(invisible(NULL))
}
