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
#' @param is_cv_model Logical. Whether the model is a cross-validated model
#'   (`TRUE`) or fitted with the full dataset (`FALSE`; default).
#' @return Saves prediction plots as JPEG files in the specified output
#'   directory.
#' @name plot_prediction
#' @author Ahmed El-Gabbas
#' @export

plot_prediction <- function(
    model_dir = NULL, env_file = ".env", n_cores = 8L, is_cv_model = FALSE) {

  .start_time <- lubridate::now(tzone = "CET")

  tif_path_mean <- tif_path_sd <- tif_path_cov <- path_clc <- path_grid <- x <-
    ias_id <- species_file <- observed <- clamp <- no_clamp <- path_pa <- NULL

  # # ..................................................................... ###
  # # ..................................................................... ###

  n_cores <- .validate_n_cores(n_cores)

  if (is.null(model_dir) || !is.character(model_dir) || !nzchar(model_dir)) {
    ecokit::stop_ctx(
      "`model_dir` has to be a character with length > 0",
      model_dir = model_dir, include_backtrace = TRUE)
  }
  if (!fs::dir_exists(model_dir)) {
    ecokit::stop_ctx(
      "`model_dir` is not a valid directory", model_dir = model_dir,
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Assign environment variables ----

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_clc", "DP_R_clc_processed", TRUE, FALSE,
    "path_grid", "DP_R_grid_processed", TRUE, FALSE,
    "path_pa", "DP_R_pa", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())
  invisible(gc())

  # Reference grid
  gird_10 <- fs::path(path_grid, "grid_10_land_crop.RData")
  if (!file.exists(gird_10)) {
    ecokit::stop_ctx(
      "Path for the Europe boundaries does not exist", gird_10 = gird_10,
      include_backtrace = TRUE)
  }
  gird_10 <- ecokit::load_as(gird_10, unwrap_r = TRUE)

  # # ..................................................................... ###

  # Load summary of prediction maps ----

  ecokit::cat_time("Load summary of prediction maps")

  # Without clamping
  map_summary_no_clamp <- fs::path(
    model_dir, "model_prediction", "no_clamp",
    "prediction_current_summary.RData")

  if (!file.exists(map_summary_no_clamp)) {
    ecokit::stop_ctx(
      "`map_summary_no_clamp` file does not exist",
      map_summary_no_clamp = map_summary_no_clamp, include_backtrace = TRUE)
  }
  map_summary_no_clamp <- ecokit::load_as(map_summary_no_clamp) %>%
    dplyr::rename(
      tif_path_mean_no_clamp = tif_path_mean,
      tif_path_sd_no_clamp = tif_path_sd,
      tif_path_cov_no_clamp = tif_path_cov)

  # With clamping
  map_summary_clamp <- fs::path(
    model_dir, "model_prediction", "clamp", "prediction_current_summary.RData")
  if (!file.exists(map_summary_clamp)) {
    ecokit::stop_ctx(
      "`map_summary_clamp` file does not exist",
      map_summary_clamp = map_summary_clamp, include_backtrace = TRUE)
  }
  map_summary_clamp <- ecokit::load_as(map_summary_clamp) %>%
    dplyr::rename(
      tif_path_mean_clamp = tif_path_mean,
      tif_path_sd_clamp = tif_path_sd,
      tif_path_cov_clamp = tif_path_cov)

  # combined data
  map_summary <- dplyr::full_join(
    map_summary_no_clamp, map_summary_clamp,
    by = c(
      "hab_abb", "hab_name", "time_period", "climate_model",
      "climate_scenario", "ias_id", "taxon_name", "species_name", "class",
      "order", "family")) %>%
    dplyr::select(-c(
      "hab_abb", "hab_name", "time_period", "climate_model",
      "climate_scenario", "taxon_name"))

  hab_abb <- map_summary_clamp$hab_abb[[1]]
  # nolint start
  hab_name <- paste0(hab_abb, ". ", map_summary_clamp$hab_name[[1]])
  # nolint end

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Calculate observed species richness ----

  ecokit::cat_time("Calculate observed species richness")

  if (is_cv_model) {
    model_data <- fs::path(model_dir, "model_data_training.RData")
    if (!ecokit::check_data(model_data)) {
      ecokit::stop_ctx(
        "Model data file not found",
        model_data = model_data, include_backtrace = TRUE)
    }
    r_sr <- ecokit::load_as(model_data) %>%
      dplyr::mutate(
        sr = rowSums(
          dplyr::select(., tidyselect::starts_with("sp")), na.rm = TRUE)) %>%
      dplyr::select(tidyselect::all_of(c("x", "y", "sr"))) %>%
      sf::st_as_sf(coords = c("x", "y"), crs = 3035) %>%
      terra::rasterize(gird_10, field = "sr") %>%
      terra::wrap()

  } else {

    model_data <- fs::path(model_dir, "model_data_subset.RData")
    if (!ecokit::check_data(model_data)) {
      ecokit::stop_ctx(
        "Model data file not found",
        model_data = model_data, include_backtrace = TRUE)
    }
    model_data <- ecokit::load_as(model_data)
    r_sr <- tibble::tibble(
      as.data.frame(model_data$data_xy), sr = rowSums(model_data$data_y)) %>%
      sf::st_as_sf(coords = c("x", "y"), crs = 3035) %>%
      terra::rasterize(gird_10, field = "sr") %>%
      terra::wrap()

  }

  rm(model_data, map_summary_no_clamp, map_summary_clamp, envir = environment())
  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Load habitat map ----

  ecokit::cat_time("Load habitat map")

  path_hab <- fs::path(
    path_clc, "summary_rdata", "perc_cover_synhab_crop.RData")
  if (!file.exists(path_hab)) {
    ecokit::stop_ctx(
      "path_hab file does not exist", path_hab = path_hab,
      include_backtrace = TRUE)
  }
  r_habitat <- ecokit::load_as(path_hab, unwrap_r = TRUE) %>%
    terra::classify(rcl = cbind(0, NA)) %>%
    terra::subset(paste0("SynHab_", hab_abb)) %>%
    terra::wrap()

  invisible(gc())

  # # ..................................................................... ###
  # # ..................................................................... ###

  path_plots <- fs::path(model_dir, "model_prediction", "plots_current")
  fs::dir_create(path_plots)

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Helper functions ----

  ## prepare_plots -----

  # helper function for generating ggplot objects
  # - map: A SpatRaster object for plotting
  # - plot_title: plot_title of the plot
  # - observed: Whether the plot is for observed data (default: FALSE)
  # - legend_title: plot_title of the legend (default: NULL)
  # - show_legend: Whether to display the legend (default: FALSE)
  # - breaks, limits: Optional parameters for legend customization (default:
  # - NULL)

  prepare_plots <- function(
    map, plot_title = NULL, observed = FALSE, legend_title = NULL,
    show_legend = FALSE, breaks = NULL, limits = NULL) {

    x_lim <- c(2600000, 6500000)
    y_lim <- c(1450000, 5420000)

    # Convert to SpatRaster if character
    if (inherits(map, "character")) {
      map <- terra::rast(map)
    }

    if (is.null(breaks)) {
      legend_breaks <- legend_labels <- ggplot2::waiver()
    } else {
      legend_breaks <- legend_labels <- breaks
    }

    if (is.null(limits)) {
      plot_limits <- unlist(terra::global(map, "range", na.rm = TRUE))
    } else {
      plot_limits <- limits
    }

    plot <- ggplot2::ggplot(environment = emptyenv()) +
      tidyterra::geom_spatraster(
        data = map, maxcell = Inf, show.legend = show_legend)
    rm(map, envir = environment())

    if (observed) {
      plot <- plot +
        ggplot2::scale_fill_manual(
          breaks = c(0, 1, 3), values = c("grey80", "blue", "red"),
          labels = c("Not observed", "Present", "Excluded"),
          na.value = "transparent", name = NULL)
    } else {
      plot <- plot +
        paletteer::scale_fill_paletteer_c(
          na.value = "transparent", palette = "viridis::plasma",
          limits = plot_limits, breaks = legend_breaks, name = legend_title,
          labels = legend_labels)
    }

    plot <- plot +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0, 0)), limits = x_lim,
        oob = scales::oob_keep) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0)), limits = y_lim) +
      ggplot2::labs(title = plot_title) +
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
        legend.position = dplyr::if_else(show_legend, "inside", "none"),
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

    return(plot)
  }

  # # ..................................................................... ###

  ## plot_maps ----

  # helper function for plotting and saving the maps

  plot_maps <- function(id) {

    # Set null device for `cairo`. This is to properly render the plots using
    # ggtext - https://github.com/wilkelab/cowplot/issues/73
    cowplot::set_null_device("cairo")

    # nolint start
    sp_id <- map_summary$ias_id[[id]]
    sp_name <- map_summary$species_name[[id]]
    class_name <- map_summary$class[[id]]
    order_name <- map_summary$order[[id]]
    family_name <- map_summary$family[[id]]
    rank_name <- dplyr::if_else(sp_id == "sr", "sr", "species")
    species <- stringr::str_detect(sp_id, "^sp")
    current_date <- format(Sys.Date(), "%d %B %Y")
    # nolint end

    # prediction map - mean
    r_mean_no_clamp <- terra::rast(map_summary$tif_path_mean_no_clamp[[id]])
    r_mean_clamp <- terra::rast(map_summary$tif_path_mean_clamp[[id]])

    # prediction map - sd
    r_sd_no_clamp <- terra::rast(map_summary$tif_path_sd_no_clamp[[id]])
    r_sd_clamp <- terra::rast(map_summary$tif_path_sd_clamp[[id]])

    # prediction map - cov
    r_cov_no_clamp <- map_summary$tif_path_cov_no_clamp[[id]] %>%
      terra::rast() %>%
      "+"(0.001) %>%
      log10()
    r_cov_clamp <- map_summary$tif_path_cov_clamp[[id]] %>%
      terra::rast() %>%
      "+"(0.001) %>%
      log10()

    # Plotting range and breaks
    if (species) {
      range_mean <- c(r_mean_clamp, r_mean_no_clamp) %>%
        terra::global("max", na.rm = TRUE) %>%
        dplyr::pull("max") %>%
        max() %>%
        c(0, .)
      breaks_mean <- NULL
      sp_id_2 <- stringr::str_remove(sp_id, "^sp_")
      path_jpeg <- fs::path(
        path_plots, paste0("pred_current_sp", sp_id_2, "_", sp_name, ".jpeg"))
      sp_id_2 <- as.integer(sp_id_2)

    } else {
      path_jpeg <- fs::path(path_plots, "pred_current_sr.jpeg")
      range_mean <- c(terra::unwrap(r_sr), r_mean_no_clamp) %>%
        terra::global("max", na.rm = TRUE) %>%
        dplyr::pull("max") %>%
        max() %>%
        c(0, .)
      breaks_mean <- NULL
    }

    range_sd <- c(r_sd_no_clamp, r_sd_clamp) %>%
      terra::global(range, na.rm = TRUE) %>%
      range()
    range_cov <- c(r_cov_no_clamp, r_cov_clamp) %>%
      terra::global(range, na.rm = TRUE) %>%
      range()

    invisible(gc())

    # ggplot objects

    ## mean_no_clamp
    plot_mean_no_clamp <- prepare_plots(
      map = r_mean_no_clamp, plot_title = "Mean", show_legend = TRUE,
      breaks = breaks_mean, limits = range_mean)
    ## mean_clamp
    plot_mean_clamp <- prepare_plots(
      map = r_mean_clamp, breaks = breaks_mean, limits = range_mean)

    ## sd_no_clamp
    plot_sd_no_clamp <- prepare_plots(
      map = r_sd_no_clamp, plot_title = "Standard deviation",
      show_legend = TRUE, breaks = NULL, limits = range_sd)
    ## sd_clamp
    plot_sd_clamp <- prepare_plots(
      map = r_sd_clamp, breaks = NULL, limits = range_sd)
    rm(r_sd_no_clamp, r_sd_clamp, envir = environment())

    ## cov_no_clamp
    plot_cov_no_clamp <- prepare_plots(
      map = r_cov_no_clamp, plot_title = "Coefficient of variation",
      show_legend = TRUE, breaks = NULL, limits = range_cov,
      legend_title = "log10")
    ## cov_clamp
    plot_cov_clamp <- prepare_plots(
      map = r_cov_clamp, breaks = NULL, limits = range_cov,
      legend_title = "log10")
    rm(r_cov_no_clamp, r_cov_clamp, envir = environment())

    invisible(gc())

    # observed data

    if (species) {

      # observed species presence

      # Files containing observed data maps
      path_observed <- fs::path(path_pa, "sp_pa_summary_df.RData")
      if (!file.exists(path_observed)) {
        ecokit::stop_ctx(
          "path_observed file does not exist", path_observed = path_observed,
          include_backtrace = TRUE)
      }
      path_observed <- ecokit::load_as(path_observed) %>%
        dplyr::filter(ias_id == sp_id_2) %>%
        dplyr::pull(species_file) %>%
        paste0(c("_masked.tif", "_all.tif")) %>%
        fs::path(path_pa, "pa_tif", .)

      # Check if observed data files exist
      if (!all(file.exists(path_observed))) {
        ecokit::stop_ctx(
          paste0("Observed data for species: ", sp_name, " not found"),
          path_observed = path_observed, sp_name = sp_name,
          include_backtrace = TRUE)
      }

      plot_observed <- terra::rast(path_observed)
      plot_observed <- terra::ifel(
        plot_observed[[1]] == 0 & plot_observed[[2]] == 1,
        3, plot_observed[[1]]) %>%
        terra::mask(terra::unwrap(r_sr)) %>%
        terra::as.factor() %>%
        prepare_plots(
          plot_title = "Species observations", observed = TRUE,
          show_legend = TRUE) +
        ggplot2::theme(
          panel.border = ggplot2::element_rect(
            color = "black", linewidth = 0.25, fill = NA,
            linetype = "dashed"))

      # Percentage habitat coverage
      plot_final <- prepare_plots(
        map = terra::unwrap(r_habitat), breaks = seq(0, 100, 20),
        limits = c(0, 100), show_legend = TRUE, plot_title = "% Habitat cover",
        legend_title = "%") +
        ggplot2::theme(
          panel.border = ggplot2::element_rect(
            color = "black", linewidth = 0.25, fill = NA, linetype = "dashed"))

    } else {

      # observed species richness
      plot_observed <- terra::unwrap(r_sr) %>%
        prepare_plots(
          plot_title = "Observed species richness", observed = FALSE,
          show_legend = TRUE, breaks = breaks_mean, limits = range_mean) +
        ggplot2::theme(
          panel.border = ggplot2::element_rect(
            color = "black", linewidth = 0.25, fill = NA,
            linetype = "dashed"))

      range_q1 <- as.vector(stats::quantile(range_mean, 0.025))
      range_q2 <- as.vector(stats::quantile(range_mean, 0.7))
      plot_limit <- range_mean + c(-2, 2)

      # observed vs predicted sr
      plot_final <- c(terra::unwrap(r_sr), r_mean_clamp, r_mean_no_clamp) %>%
        terra::as.data.frame(na.rm = TRUE) %>%
        stats::setNames(c("observed", "clamp", "no_clamp")) %>%
        ggplot2::ggplot(
          mapping = ggplot2::aes(x = observed), environment = emptyenv()) +
        ggplot2::geom_point(
          ggplot2::aes(y = clamp, colour = "with clamping"),
          shape = 17, size = 0.03, alpha = 0.2) +
        ggplot2::geom_point(
          ggplot2::aes(y = no_clamp, colour = "without clamping"),
          shape = 16, size = 0.03, alpha = 0.2) +
        ggplot2::scale_colour_manual(
          name = NULL, drop = FALSE,
          values = c("with clamping" = "red", "without clamping" = "blue")) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2) +
        ggplot2::coord_equal(
          xlim = plot_limit, ylim = plot_limit, expand = FALSE, clip = "off") +
        ggplot2::annotate(
          "text", x = range_q2, y = range_q1,
          angle = 0, size = 3, color = "darkgrey",
          label = "Observed species richness", hjust = 0.5, vjust = 0.5) +
        ggplot2::annotate(
          "text", x = range_q1, y = range_q2,
          angle = 90, size = 3, color = "darkgrey",
          label = "Predicted species richness", hjust = 0.5, vjust = 0.5) +
        ggplot2::labs(title = "Observed vs. predicted species richness") +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.margin = ggplot2::margin(0, 0.05, 0, 0.05, "cm"),
          plot.title = ggplot2::element_text(
            size = 8, color = "grey60", face = "bold", hjust = 0.5,
            margin = ggplot2::margin(0, 0, 0, 0)),
          legend.position = "inside",
          legend.position.inside = c(0.4, 0.96),
          legend.direction = "horizontal",
          legend.background = ggplot2::element_rect(
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
      plot_mean_no_clamp, plot_mean_clamp, plot_sd_no_clamp, plot_sd_clamp,
      plot_cov_no_clamp, plot_cov_clamp, plot_observed, plot_final,
      ncol = 4, nrow = 2, byrow = FALSE) +
      ggplot2::theme(
        plot.margin = grid::unit(c(0.05, 0.05, 0.05, 0.05), "cm"),
        panel.spacing = grid::unit(0.1, "lines"))

    rm(
      plot_mean_no_clamp, plot_mean_clamp, plot_sd_no_clamp, plot_sd_clamp,
      plot_cov_no_clamp, plot_cov_clamp, plot_observed, plot_final,
      envir = environment())

    y_lab <- cowplot::ggdraw() +
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

    if (species) {
      plot_title_1 <- stringr::str_glue(
        '<SPAN STYLE="font-size:16pt">Predicted habitat suitability of \\
          <b><i>{sp_name}</i></b></SPAN>')
      plot_title_2 <- stringr::str_glue(
        '<SPAN STYLE="font-size:8pt; color: darkgrey"><b>class:</b> \\
          {class_name}; <b>Order:</b> {order_name}; <b>Family:</b> \\
          {family_name}; <b>ias_id:</b> {sp_id_2}</SPAN>')
    } else {
      plot_title_1 <- "Predicted level of invasion (mean species richness)"
      plot_title_2 <- ""
    }

    hab_name_0 <- stringr::str_remove(hab_name, " habitats")      # nolint: object_name_linter
    main_title <- cowplot::ggdraw() +
      ggtext::geom_richtext(
        ggplot2::aes(x = 0.01, y = 0.6, label = plot_title_1),
        fill = NA, label.color = NA, hjust = 0, vjust = 0.5, size = 5,
        color = "black") +
      ggtext::geom_richtext(
        ggplot2::aes(x = 0.01, y = 0.2, label = plot_title_2),
        fill = NA, label.color = NA, hjust = 0, vjust = 0.5, size = 5,
        color = "black") +
      ggtext::geom_richtext(
        ggplot2::aes(
          x = 1, y = 0.55,
          label = stringr::str_glue(
            '<SPAN STYLE="font-size:12.5pt; color: red"><b>{hab_name_0} \\
              habitats </b></SPAN>&#8212;<SPAN STYLE="font-size:12.5pt; \\
              color: blue"><b> current climate</b></SPAN>')),
        fill = NA, label.color = NA, hjust = 1, vjust = 0.5) +
      ggtext::geom_richtext(
        ggplot2::aes(
          x = 1, y = 0.25,
          label = stringr::str_glue(
            paste0(
              '<span style="color: grey; font-size:6pt">Last',
              "update: {current_date}</span>"))),
        fill = NA, label.color = NA, hjust = 1, vjust = 0.5) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))

    final_plot <- cowplot::plot_grid(
      main_title,
      cowplot::plot_grid(
        y_lab, plot_grid_main, ncol = 2, rel_widths = c(0.03, 1), align = "h"),
      ncol = 1, rel_heights = c(0.07, 1))

    ragg::agg_jpeg(
      filename = path_jpeg, width = 30, height = 15.5, res = 600,
      quality = 100, units = "cm")
    print(final_plot)
    grDevices::dev.off()

    invisible(NULL)
  }

  # # ..................................................................... ###
  # # ..................................................................... ###

  # Plotting ----

  ecokit::cat_time("Plotting")

  # set up parallel processing
  doParallel::registerDoParallel(cores = n_cores)
  ecokit::load_packages(package_list = "foreach")
  withr::defer(doParallel::stopImplicitCluster())

  plots <- foreach::foreach(
    x = seq_len(nrow(map_summary)),
    .export = c(
      "map_summary", "prepare_plots", "r_sr", "path_pa",
      "r_habitat", "path_plots", "hab_name", "plot_maps"),
    .packages = c(
      "dplyr", "terra", "ggplot2", "stringr", "cowplot", "tidyterra",
      "purrr", "ggtext", "ragg", "paletteer", "grid", "scales", "ecokit",
      "magrittr")) %dopar% {
        plot_maps(x)
      }

  rm(plots, envir = environment())
  doParallel::stopImplicitCluster()

  ecokit::cat_diff(init_time = .start_time)

  return(invisible(NULL))
}
