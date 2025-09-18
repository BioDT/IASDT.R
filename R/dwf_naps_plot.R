# # |------------------------------------------------------------------------| #
# naps_plot ----
## |------------------------------------------------------------------------| #

#' @export
#' @author Ahmed El-Gabbas
#' @name naps_data
#' @rdname naps_data
#' @order 3

naps_plot <- function(species = NULL, env_file = ".env") {

  # # ..................................................................... ###

  # Checking arguments ----
  ecokit::check_args(args_to_check = "species", args_type = "character")

  species2 <- ecokit::replace_space(species)

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  path_grid_ref <- path_grid <- path_taxa_info_rdata <- path_pa <- NAME_ENGL <-
    path_taxa_info <- eu_boundaries <- species_name2 <- CellCode <- x <-
    y <- NULL

  # # ..................................................................... ###

  # Environment variables ----

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_grid", "DP_R_grid_processed", TRUE, FALSE,
    "path_grid_ref", "DP_R_grid_raw", TRUE, FALSE,
    "eu_boundaries", "DP_R_country_boundaries", FALSE, TRUE,
    "path_pa", "DP_R_pa", FALSE, FALSE,
    "path_taxa_info_rdata", "DP_R_taxa_info_rdata", FALSE, TRUE,
    "path_taxa_info", "DP_R_taxa_info", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)

  rm(env_vars_to_read, envir = environment())
  invisible(gc())

  # Summary information for the current species
  species_info <- ecokit::load_as(path_taxa_info_rdata)
  if (nrow(species_info) == 0) {
    ecokit::stop_ctx(
      "Species information could not be found",
      path_taxa_info_rdata = path_taxa_info_rdata,
      species_info = species_info, nrow_species_info = nrow(species_info),
      include_backtrace = TRUE)
  }
  if (!species %in% species_info$species_name) {
    ecokit::stop_ctx(
      "Species name not found in the list of species",
      species = species, include_backtrace = TRUE)
  }

  species_info <- dplyr::filter(species_info, species_name2 == species2) %>%
    dplyr::select(-"speciesKey") %>%
    dplyr::distinct()

  # # ..................................................................... ###

  # Check / create directories
  path_summary <- fs::path(path_pa, "pa_summary")
  path_jpeg <- fs::path(path_pa, "distribution_jpeg")
  if (!fs::dir_exists(path_jpeg)) {
    fs::dir_create(path_jpeg)
  }
  out_path <- fs::path(path_jpeg, paste0(species_info$species_file[1], ".jpeg"))

  # Loading species data
  selected_columns <- c(
    "n_cells_all", "gbif_grid100", "easin_grid100", "elter_grid100",
    "species_id", "pa_map", "countries_to_exclude",
    "gbif", "gbif_unique", "easin", "easin_unique", "elter", "elter_unique",
    "n_cells_naturalized", "gbif_masked", "gbif_masked_unique",
    "easin_masked", "easin_masked_unique", "elter_masked",
    "elter_masked_unique")

  species_data <- fs::path(
    path_summary, paste0(species_info$species_file[1], ".qs2"))
  if (!file.exists(species_data)) {
    return(invisible(NULL))
  }
  species_data <- ecokit::load_as(species_data) %>%
    dplyr::select(tidyselect::all_of(selected_columns))
  invisible(gc())

  if (species_data$n_cells_all == 0) {
    return(invisible(NULL))
  }

  country_boundaries <- ecokit::load_as(eu_boundaries) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_10")
  last_update <- paste0("Last update: ", format(Sys.Date(), "%d %B %Y"))


  grid_100_sf <- fs::path(path_grid_ref, "grid_100_sf.RData") %>%
    ecokit::load_as() %>%
    magrittr::extract2("grid_100_sf_s")

  # Location of legend
  legend_gbif <- dplyr::filter(grid_100_sf, CellCode == "100kmE27N45")
  legend_easin <- dplyr::filter(grid_100_sf, CellCode == "100kmE27N44")
  legend_elter <- dplyr::filter(grid_100_sf, CellCode == "100kmE27N43")

  rm(species_info, grid_100_sf, envir = environment())
  invisible(gc())

  # the study area as simple feature object for plotting
  grid10_sf <- fs::path(path_grid, "grid_10_land_sf.RData") %>%
    ecokit::load_as()

  gbif_grid100 <- species_data$gbif_grid100[[1]]
  easin_grid100 <- species_data$easin_grid100[[1]]
  elter_grid100 <- species_data$elter_grid100[[1]]

  species_data <- dplyr::select(
    species_data,
    -tidyselect::all_of(c("gbif_grid100", "easin_grid100", "elter_grid100")))
  invisible(gc())

  # # |||||||||||||||||||||||||||||||||||||||||||||||||

  # Figure title ----

  ias_id <- unique(species_data$species_id)
  map_title <- readr::read_tsv(
    file = path_taxa_info, show_col_types = FALSE, progress = FALSE) %>%
    dplyr::filter(species_name2 == species2) %>%
    dplyr::select(tidyselect::all_of(c("class", "order", "family"))) %>%
    unlist() %>%
    stringr::str_c(collapse = " / ") %>%
    paste0(
      "   <span style='font-size: 14pt; color:blue;'><b><i>", ias_id,
      " \u2014 ", species, "</i></b></span>   (", ., ")")

  rm(path_taxa_info, envir = environment())
  invisible(gc())

  # # |||||||||||||||||||||||||||||||||||||||||||||||||

  # plotting theme ----

  plot_theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
      plot.title = ggtext::element_markdown(
        hjust = 0.04, margin = ggplot2::margin(4, 0, 2, 0)),
      plot.subtitle = ggtext::element_markdown(
        hjust = 0.05, margin = ggplot2::margin(2, 0, 0, 0)),
      strip.text = ggplot2::element_text(size = 6, face = "bold"),
      legend.position = "none",
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.spacing = grid::unit(0.3, "lines"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      plot.tag.position = c(0.94, 0.011),
      plot.tag = ggtext::element_markdown(colour = "grey", size = 4),
      panel.ontop = TRUE,
      panel.background = ggplot2::element_rect(fill = NA))

  # # |||||||||||||||||||||||||||||||||||||||||||||||||

  # map to be plotted ----
  sp_pa <- terra::unwrap(species_data$pa_map[[1]]) %>%
    terra::classify(cbind(0, NA)) %>%
    terra::as.factor() %>%
    as.data.frame(xy = TRUE) %>%
    stats::setNames(c("x", "y", "species"))
  species_data$pa_map <- NULL

  exclude_boundary <- ecokit::load_as(eu_boundaries) %>%
    magrittr::extract2("Bound_sf_Eur") %>%
    magrittr::extract2("L_10") %>%
    dplyr::filter(NAME_ENGL %in% species_data$countries_to_exclude[[1]]) %>%
    dplyr::select("NAME_ENGL")

  rm(eu_boundaries, envir = environment())
  invisible(gc())

  # # |||||||||||||||||||||||||||||||||||||||||||||||||

  # Figure subtitle ----

  n_grids_all <- paste0(
    "<span style='font-size: 12pt; color:red;'><b>All data:</b></span> ",
    scales::label_comma()(species_data$n_cells_all),
    " presence grid cells</b> \u2014 <b>GBIF</b> (",
    scales::label_comma()(species_data$gbif), " / ",
    scales::label_comma()(species_data$gbif_unique), ") \u2014 <b>EASIN</b> (",
    scales::label_comma()(species_data$easin), " / ",
    scales::label_comma()(species_data$easin_unique), ") \u2014 <b>eLTER</b> (",
    scales::label_comma()(species_data$elter), " / ",
    scales::label_comma()(species_data$elter_unique), ")<br>",
    "<span style='font-size: 12pt; color:red;'><b>Final data: </span>",
    scales::label_comma()(species_data$n_cells_Naturalized),
    " presence grid cells</b> \u2014 <b>GBIF</b> (",
    scales::label_comma()(species_data$gbif_masked), " / ",
    scales::label_comma()(species_data$gbif_masked_unique),
    ") \u2014 <b>EASIN</b> (",
    scales::label_comma()(species_data$easin_masked), " / ",
    scales::label_comma()(species_data$easin_masked_unique),
    ") \u2014 <b>eLTER</b> (",
    scales::label_comma()(species_data$elter_masked), " / ",
    scales::label_comma()(species_data$elter_masked_unique), ")\n")

  # # |||||||||||||||||||||||||||||||||||||||||||||||||

  # ggplot object -----

  plot <- ggplot2::ggplot() +
    # country boundaries
    ggplot2::geom_sf(
      country_boundaries, mapping = ggplot2::aes(), color = "grey",
      linewidth = 0.2, fill = scales::alpha("grey", 0.2)) +
    # the study area
    ggplot2::geom_sf(
      grid10_sf, mapping = ggplot2::aes(), color = "lightgrey",
      fill = "lightgrey", linewidth = 0.15)

  if (nrow(exclude_boundary) > 0) {
    plot <- plot +
      # Countries to exclude
      ggplot2::geom_sf(
        exclude_boundary, mapping = ggplot2::aes(), color = "grey90",
        linewidth = 0.1, fill = scales::alpha("red", 0.2))
  }

  # presence grids at 100 km resolution
  if (nrow(gbif_grid100) > 0) {
    plot <- plot +
      ggplot2::geom_sf(
        gbif_grid100, mapping = ggplot2::aes(), color = "transparent",
        fill = scales::alpha("blue", 0.2), linewidth = 0.15)
  }

  # PA grid at 10 km resolution
  plot <- plot +
    ggplot2::geom_tile(
      data = sp_pa, mapping = ggplot2::aes(x = x, y = y, fill = species)) +
    ggplot2::scale_fill_manual(
      values = c("blue", "transparent"), na.value = "transparent")

  if (nrow(easin_grid100) > 0) {
    plot <- plot +
      ggplot2::geom_sf(
        easin_grid100, mapping = ggplot2::aes(), color = "red",
        fill = "transparent", linewidth = 0.65)
  }
  if (nrow(elter_grid100) > 0) {
    plot <- plot +
      ggplot2::geom_sf(
        elter_grid100, mapping = ggplot2::aes(), color = "darkgreen",
        fill = "transparent", linewidth = 0.65, linetype = "dotdash")
  }

  plot <- plot +
    # country boundaries
    ggplot2::geom_sf(
      country_boundaries, mapping = ggplot2::aes(), color = "black",
      linewidth = 0.2, fill = "transparent") +
    # legends for data type
    ggplot2::geom_sf(
      legend_gbif, mapping = ggplot2::aes(), color = "transparent",
      fill = scales::alpha("blue", 0.2), linewidth = 0.15) +
    ggplot2::geom_sf(
      legend_easin, mapping = ggplot2::aes(), color = "red",
      fill = "transparent", linewidth = 0.65) +
    ggplot2::geom_sf(
      legend_elter, mapping = ggplot2::aes(), color = "darkgreen",
      fill = "transparent", linewidth = 0.65, linetype = "dotdash") +
    ggplot2::geom_text(
      mapping = ggplot2::aes(
        x = 2650000, y = 4650000, label = "Species presence per 100 km grid"),
      hjust = 0, fontface = "bold") +
    ggplot2::geom_text(
      mapping = ggplot2::aes(x = 2850000, y = 4650000 - 1e5, label = "GBIF"),
      hjust = 0) +
    ggplot2::geom_text(
      mapping = ggplot2::aes(x = 2850000, y = 4650000 - 2e5, label = "EASIN"),
      hjust = 0) +
    ggplot2::geom_text(
      mapping = ggplot2::aes(x = 2850000, y = 4650000 - 3e5, label = "eLTER"),
      hjust = 0) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)),
      limits = c(2600000, 6550000)) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)),
      limits = c(1450000, 5410000)) +
    ggplot2::labs(title = map_title, subtitle = n_grids_all, fill = NULL) +
    plot_theme

  rm(
    species_data, grid10_sf, gbif_grid100, easin_grid100, elter_grid100,
    exclude_boundary, legend_gbif, legend_easin, legend_elter,
    envir = environment())
  invisible(gc())

  # # |||||||||||||||||||||||||||||||||||||||||||||||||

  # Save the plot as JPEG file ----

  # Use `grid::grid.text` instead of `cowplot::draw_label` to avoid these
  # warnings
  #
  # 1: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), ... :
  # conversion failure on '—' in 'mbcsToSbcs': dot substituted for <e2>
  # 2: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), ... :
  # conversion failure on '—' in 'mbcsToSbcs': dot substituted for <80>
  # 3: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), ... :
  # conversion failure on '—' in 'mbcsToSbcs': dot substituted for <94>
  #
  # Using ggplot2::ggsave directly does not show non-ascii characters correctly

  ragg::agg_jpeg(
    filename = out_path, width = 25, height = 26.5,
    res = 600, quality = 100, units = "cm")
  print(plot)
  grid::grid.text(
    label = last_update, x = 0.98, y = 0.975, hjust = 1, vjust = 1,
    gp = grid::gpar(col = "grey65", fontsize = 12))
  grDevices::dev.off()

  return(out_path)
}
