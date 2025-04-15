## |------------------------------------------------------------------------| #
# plot_prediction ----
## |------------------------------------------------------------------------| #

#' Plot species and level of invasion predictions as JPEG files using `ggplot2`
#'
#' Generate predictions for species and habitat models and saves the output as
#' JPEG files.
#'
#' @param model_dir Character. Path to the model directory containing
#'   predictions.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 8.
#' @return Saves prediction plots as JPEG files in the specified output
#'   directory.
#' @name plot_prediction
#' @author Ahmed El-Gabbas
#' @export

plot_prediction <- function(model_dir = NULL, env_file = ".env", n_cores = 8L) {

  # # ..................................................................... ###
  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  tif_path_mean <- tif_path_sd <- tif_path_cov <- Path_CLC <- Path_Grid <-
    IAS_ID <- Species_File <- Observed <- Clamp <- NoClamp <- Path_PA <- NULL

  # # ..................................................................... ###
  # # ..................................................................... ###

  if (is.null(model_dir) || !is.character(model_dir) || !nzchar(model_dir)) {
    stop("`model_dir` has to be a character with length > 0", call. = FALSE)
  }
  if (!fs::dir_exists(model_dir)) {
    stop("`model_dir` is not a valid directory", call. = FALSE)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Assign environment variables ----

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_CLC", "DP_R_CLC_processed", TRUE, FALSE,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE,
    "Path_PA", "DP_R_PA", TRUE, FALSE)
  # Assign environment variables and check file and paths
  IASDT.R::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())
  invisible(gc())

  # Reference grid
  Gird10 <- IASDT.R::path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(Gird10)) {
    stop(
      "Path for the Europe boundaries does not exist: ", Gird10, call. = FALSE)
  }
  Gird10 <- IASDT.R::load_as(Gird10) %>%
    terra::unwrap()

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Load summary of prediction maps ----

  IASDT.R::cat_time("Load summary of prediction maps")

  # Without clamping
  Map_summary_NoClamp <- IASDT.R::path(
    model_dir, "Model_Prediction", "NoClamp",
    "Prediction_Current_Summary.RData")

  if (!file.exists(Map_summary_NoClamp)) {
    stop(
      "`Map_summary_NoClamp` file: ", Map_summary_NoClamp,
      " does not exist", call. = FALSE)
  }
  Map_summary_NoClamp <- IASDT.R::load_as(Map_summary_NoClamp) %>%
    dplyr::rename(
      tif_path_mean_no_clamp = tif_path_mean,
      tif_path_sd_no_clamp = tif_path_sd,
      tif_path_cov_no_clamp = tif_path_cov)

  # With clamping
  Map_summary_Clamp <- IASDT.R::path(
    model_dir, "Model_Prediction", "Clamp", "Prediction_Current_Summary.RData")
  if (!file.exists(Map_summary_Clamp)) {
    stop(
      "`Map_summary_Clamp` file: ", Map_summary_Clamp,
      " does not exist", call. = FALSE)
  }
  Map_summary_Clamp <- IASDT.R::load_as(Map_summary_Clamp) %>%
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
  # nolint start
  Hab_Name <- paste0(HabAbb, ". ", Map_summary_Clamp$hab_name[[1]])
  # nolint end

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Calculate observed species richness ----

  IASDT.R::cat_time("Calculate observed species richness")

  # Modelling data
  Model_Data <- list.files(
    model_dir, pattern = "^ModDT_.*subset.RData$", full.names = TRUE)

  if (length(Model_Data) != 1) {
    stop("Model data does not exist", call. = FALSE)
  }
  Model_Data <- IASDT.R::load_as(Model_Data)

  # Observed species richness
  R_SR <- tibble::tibble(
    as.data.frame(Model_Data$DT_xy), SR = rowSums(Model_Data$DT_y)) %>%
    sf::st_as_sf(coords = c("x", "y"), crs = 3035) %>%
    terra::rasterize(Gird10, field = "SR") %>%
    terra::wrap()

  rm(Model_Data, Map_summary_NoClamp, Map_summary_Clamp, envir = environment())
  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Load habitat map ----

  IASDT.R::cat_time("Load habitat map")

  Path_Hab <- IASDT.R::path(
    Path_CLC, "Summary_RData", "PercCov_SynHab_Crop.RData")
  if (!file.exists(Path_Hab)) {
    stop("Path_Hab file: ", Path_Hab, " does not exist", call. = FALSE)
  }
  R_habitat <- IASDT.R::load_as(Path_Hab) %>%
    terra::unwrap() %>%
    terra::classify(c(0, NA)) %>%
    terra::subset(paste0("SynHab_", HabAbb)) %>%
    terra::wrap()

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  Path_Plots <- IASDT.R::path(model_dir, "Model_Prediction", "Plots_Current")
  fs::dir_create(Path_Plots)

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Helper functions ----

  ## PrepPlots -----

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

    # Convert to SpatRaster if character
    if (inherits(Map, "character")) {
      Map <- terra::rast(Map)
    }

    if (is.null(breaks)) {
      LegendBreaks <- Legendlabels <- ggplot2::waiver()
    } else {
      LegendBreaks <- Legendlabels <- breaks
    }

    if (is.null(limits)) {
      PlotLimits <- unlist(terra::global(Map, "range", na.rm = TRUE))
    } else {
      PlotLimits <- limits
    }

    Plot <- ggplot2::ggplot() +
      tidyterra::geom_spatraster(
        data = Map, maxcell = Inf, show.legend = ShowLegend)

    rm(Map, envir = environment())

    if (Observed) {
      Plot <- Plot +
        ggplot2::scale_fill_manual(
          breaks = c(0, 1, 3), values = c("grey80", "blue", "red"),
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
        plot.margin = ggplot2::margin(0, 0.05, 0, 0.05, "cm"),
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
        panel.spacing = grid::unit(0.2, "lines"),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(
          linewidth = 0.05, colour = "grey40", linetype = 2),
        panel.border = ggplot2::element_blank(),
        panel.ontop = TRUE,
        panel.background = ggplot2::element_rect(fill = NA))

    return(Plot)
  }

  # # ..................................................................... ###

  ## PlotMaps ----

  # helper function for plotting and saving the maps

  PlotMaps <- function(ID) {

    # nolint start
    SpID <- Map_summary$ias_id[[ID]]
    SpName <- Map_summary$species_name[[ID]]
    ClassName <- Map_summary$class[[ID]]
    OrderName <- Map_summary$order[[ID]]
    FamilyName <- Map_summary$family[[ID]]
    Rank <- dplyr::if_else(SpID == "SR", "SR", "Species")
    Species <- stringr::str_detect(SpID, "^Sp")
    Date <- format(Sys.Date(), "%d %B %Y")
    # nolint end

    # prediction map - mean
    R_mean_NoClamp <- terra::rast(Map_summary$tif_path_mean_no_clamp[[ID]])
    R_mean_Clamp <- terra::rast(Map_summary$tif_path_mean_clamp[[ID]])

    # prediction map - sd
    R_sd_NoClamp <- terra::rast(Map_summary$tif_path_sd_no_clamp[[ID]])
    R_sd_Clamp <- terra::rast(Map_summary$tif_path_sd_clamp[[ID]])

    # prediction map - cov
    R_cov_NoClamp <- Map_summary$tif_path_cov_no_clamp[[ID]] %>%
      terra::rast() %>%
      "+"(0.001) %>%
      log10()
    R_cov_Clamp <- Map_summary$tif_path_cov_clamp[[ID]] %>%
      terra::rast() %>%
      "+"(0.001) %>%
      log10()

    # Plotting range and breaks
    if (Species) {
      Range_Mean <- c(R_mean_Clamp, R_mean_NoClamp) %>%
        terra::global("max", na.rm = TRUE) %>%
        dplyr::pull("max") %>%
        max() %>%
        c(0, .)
      Breaks_Mean <- NULL
      SpID2 <- stringr::str_remove(SpID, "^Sp_")
      path_JPEG <- IASDT.R::path(
        Path_Plots, paste0("Pred_Current_Sp", SpID2, "_", SpName, ".jpeg"))
      SpID2 <- as.integer(SpID2)

    } else {
      path_JPEG <- IASDT.R::path(Path_Plots, "Pred_Current_SR.jpeg")
      Range_Mean <- c(terra::unwrap(R_SR), R_mean_NoClamp) %>%
        terra::global("max", na.rm = TRUE) %>%
        dplyr::pull("max") %>%
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

    invisible(gc())

    # ggplot objects

    ## mean_no_clamp
    Plot_mean_no_clamp <- PrepPlots(
      Map = R_mean_NoClamp, Title = "Mean", ShowLegend = TRUE,
      breaks = Breaks_Mean, limits = Range_Mean)
    ## mean_clamp
    Plot_mean_clamp <- PrepPlots(
      Map = R_mean_Clamp, breaks = Breaks_Mean, limits = Range_Mean)

    ## sd_no_clamp
    Plot_sd_no_clamp <- PrepPlots(
      Map = R_sd_NoClamp, Title = "Standard deviation", ShowLegend = TRUE,
      breaks = NULL, limits = Range_sd)
    ## sd_clamp
    Plot_sd_clamp <- PrepPlots(
      Map = R_sd_Clamp, breaks = NULL, limits = Range_sd)
    rm(R_sd_NoClamp, R_sd_Clamp, envir = environment())

    ## cov_no_clamp
    Plot_cov_no_clamp <- PrepPlots(
      Map = R_cov_NoClamp, Title = "Coefficient of variation",
      ShowLegend = TRUE, breaks = NULL, limits = Range_cov,
      LegendTitle = "log10")
    ## cov_clamp
    Plot_cov_clamp <- PrepPlots(
      Map = R_cov_Clamp, breaks = NULL, limits = Range_cov,
      LegendTitle = "log10")
    rm(R_cov_NoClamp, R_cov_Clamp, envir = environment())

    invisible(gc())

    # Observed data
    if (Species) {

      # Observed species presence

      # Files containing observed data maps
      Path_observed <- IASDT.R::path(Path_PA, "Sp_PA_Summary_DF.RData")
      if (!file.exists(Path_observed)) {
        stop(
          "Path_observed file: ", Path_observed, " does not exist",
          call. = FALSE)
      }
      Path_observed <- IASDT.R::load_as(Path_observed) %>%
        dplyr::filter(IAS_ID == SpID2) %>%
        dplyr::pull(Species_File) %>%
        paste0(c("_Masked.tif", "_All.tif")) %>%
        IASDT.R::path(Path_PA, "tif", .)

      # Check if observed data files exist
      if (!all(file.exists(Path_observed))) {
        stop(
          "Observed data for species: ", SpName, " not found", call. = FALSE)
      }

      Plot_observed <- terra::rast(Path_observed)
      Plot_observed <- terra::ifel(
        Plot_observed[[1]] == 0 & Plot_observed[[2]] == 1,
        3, Plot_observed[[1]]) %>%
        terra::mask(terra::unwrap(R_SR)) %>%
        terra::as.factor() %>%
        PrepPlots(
          Title = "Species observations", Observed = TRUE, ShowLegend = TRUE) +
        ggplot2::theme(
          panel.border = ggplot2::element_rect(
            color = "black", linewidth = 0.25, fill = NA,
            linetype = "dashed"))

      # Percentage habitat coverage
      Plot_Final <- PrepPlots(
        Map = terra::unwrap(R_habitat), breaks = seq(0, 100, 20),
        limits = c(0, 100), ShowLegend = TRUE, Title = "% Habitat cover",
        LegendTitle = "%") +
        ggplot2::theme(
          panel.border = ggplot2::element_rect(
            color = "black", linewidth = 0.25, fill = NA, linetype = "dashed"))

    } else {

      # Observed species richness
      Plot_observed <- terra::unwrap(R_SR) %>%
        PrepPlots(
          Title = "Observed species richness", Observed = FALSE,
          ShowLegend = TRUE, breaks = Breaks_Mean, limits = Range_Mean) +
        ggplot2::theme(
          panel.border = ggplot2::element_rect(
            color = "black", linewidth = 0.25, fill = NA,
            linetype = "dashed"))

      Range_q1 <- as.vector(stats::quantile(Range_Mean, 0.025))
      Range_q2 <- as.vector(stats::quantile(Range_Mean, 0.7))
      Plot_limit <- Range_Mean + c(-2, 2)

      # Observed vs predicted SR
      Plot_Final <- c(terra::unwrap(R_SR), R_mean_Clamp, R_mean_NoClamp) %>%
        terra::as.data.frame(na.rm = TRUE) %>%
        stats::setNames(c("Observed", "Clamp", "NoClamp")) %>%
        ggplot2::ggplot(mapping = ggplot2::aes(x = Observed)) +
        ggplot2::geom_point(
          ggplot2::aes(y = Clamp, colour = "with clamping"),
          shape = 17, size = 0.03, alpha = 0.2) +
        ggplot2::geom_point(
          ggplot2::aes(y = NoClamp, colour = "without clamping"),
          shape = 16, size = 0.03, alpha = 0.2) +
        ggplot2::scale_colour_manual(
          name = NULL, drop = FALSE,
          values = c("with clamping" = "red", "without clamping" = "blue")) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2) +
        ggplot2::coord_equal(
          xlim = Plot_limit, ylim = Plot_limit, expand = FALSE, clip = "off") +
        ggplot2::annotate(
          "text", x = Range_q2, y = Range_q1,
          angle = 0, size = 3, color = "darkgrey",
          label = "Observed species richness", hjust = 0.5, vjust = 0.5) +
        ggplot2::annotate(
          "text", x = Range_q1, y = Range_q2,
          angle = 90, size = 3, color = "darkgrey",
          label = "Predicted species richness", hjust = 0.5, vjust = 0.5) +
        ggplot2::labs(title =  "Observed vs. predicted species richness") +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.margin = ggplot2::margin(0, 0.05, 0, 0.05, "cm"),
          plot.title = ggplot2::element_text(
            size = 8, color = "grey60", face = "bold", hjust = 0.5,
            margin = ggplot2::margin(0, 0, 0, 0)),
          legend.position = "inside",
          legend.position.inside = c(0.4, 0.96),
          legend.direction = "horizontal",
          legend.background =  ggplot2::element_rect(
            fill = "transparent", colour = "transparent"),
          legend.margin = ggplot2::margin(0, 4, 0, 0),
          legend.text = ggplot2::element_text(size = 6, vjust = 0.5),
          legend.box.spacing = grid::unit(0, "pt"),
          legend.title = ggplot2::element_text(
            color = "blue", size = 6, face = "bold", hjust = 0, vjust = 0),
          legend.spacing.x = grid::unit(0.05, "cm"),
          legend.key.width = grid::unit(0.15, "cm"),
          axis.text.x = ggplot2::element_text(size = 5),
          axis.text.y = ggplot2::element_text(
            size = 5, hjust = 0.5, angle = 90),
          axis.ticks = ggplot2::element_blank(),
          axis.title = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_line(
            linewidth = 0.025, colour = "grey40", linetype = 2),
          panel.border = ggplot2::element_rect(
            color = "black", linewidth = 0.25, fill = NA,
            linetype = "dashed"),
          panel.background = ggplot2::element_rect(fill = NA)) +
        ggplot2::guides(
          colour = ggplot2::guide_legend(
            override.aes = list(size = 1.5, shape = c(17, 16), alpha = 1)))
    }

    plot_grid_main <- cowplot::plot_grid(
      Plot_mean_no_clamp, Plot_mean_clamp, Plot_sd_no_clamp, Plot_sd_clamp,
      Plot_cov_no_clamp, Plot_cov_clamp, Plot_observed, Plot_Final,
      ncol = 4, nrow = 2, byrow = FALSE) +
      ggplot2::theme(
        plot.margin = grid::unit(c(0.05, 0.05, 0.05, 0.05), "cm"),
        panel.spacing = grid::unit(0.1, "lines"))

    rm(
      Plot_mean_no_clamp, Plot_mean_clamp, Plot_sd_no_clamp, Plot_sd_clamp,
      Plot_cov_no_clamp, Plot_cov_clamp, Plot_observed, Plot_Final,
      envir = environment())

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
              (efforts/rivers predictors fixed at 90% quantile)</span>')
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

    Hab_Name0 <- stringr::str_remove(Hab_Name, " habitats")
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
            '<SPAN STYLE="font-size:12.5pt; color: red"><b>{Hab_Name0} \\
              habitats </b></SPAN>&#8212;<SPAN STYLE="font-size:12.5pt; \\
              color: blue"><b> current climate</b></SPAN>')),
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
      filename = path_JPEG, width = 30, height = 15.5, res = 600,
      quality = 100, units = "cm")
    print(final_plot)
    grDevices::dev.off()

    return(invisible(NULL))
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Plotting ----

  IASDT.R::cat_time("Plotting")

  IASDT.R::cat_time(
    paste0("Preparing working in parallel using ", n_cores, " cores"),
    level = 1)
  c1 <- parallel::makePSOCKcluster(n_cores)
  on.exit(try(parallel::stopCluster(c1), silent = TRUE), add = TRUE)

  IASDT.R::cat_time("Exporting variables to parallel cores", level = 1)
  parallel::clusterExport(
    cl = c1,
    varlist = c(
      "Map_summary", "PrepPlots", "R_SR", "Path_PA",
      "R_habitat", "Path_Plots", "Hab_Name", "PlotMaps"),
    envir = environment())

  IASDT.R::cat_time("Loading packages at parallel cores", level = 1)
  invisible(parallel::clusterEvalQ(
    cl = c1,
    expr = {
      sapply(
        c("dplyr", "terra", "ggplot2", "stringr", "cowplot", "tidyterra",
          "purrr", "ggtext", "ragg", "paletteer", "grid", "scales"),
        library, character.only = TRUE)

      # Set null device for `cairo`. This is to properly render the plots using
      # ggtext - https://github.com/wilkelab/cowplot/issues/73
      cowplot::set_null_device("cairo")
    }))

  IASDT.R::cat_time("Prepare and save plots in parallel", level = 1)
  Plots <- parallel::parLapply(
    cl = c1, X = seq_len(nrow(Map_summary)), fun = PlotMaps)

  rm(Plots, envir = environment())

  IASDT.R::cat_diff(init_time = .StartTime)

  return(invisible(NULL))
}
