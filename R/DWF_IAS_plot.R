# # |------------------------------------------------------------------------| #
# IAS_plot ----
## |------------------------------------------------------------------------| #

#' @export
#' @author Ahmed El-Gabbas
#' @name IAS_data
#' @rdname IAS_data
#' @order 3

IAS_plot <- function(species = NULL, env_file = ".env") {

  # # ..................................................................... ###

  Species2 <- ecokit::replace_space(species)

  # Checking arguments ----

  if (is.null(species)) {
    ecokit::stop_ctx(
      "species cannot be empty", species = species, include_backtrace = TRUE)
  }

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("species", "env_file"))

  rm(AllArgs, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Grid_Ref <- Path_Grid <- Path_TaxaInfo_RData <- Path_PA <- NAME_ENGL <-
    Path_TaxaInfo <- EU_Bound <- Species_name2 <- CellCode <- x <- y <- NULL

  # # ..................................................................... ###

  # Environment variables ----

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE,
    "Path_Grid_Ref", "DP_R_Grid_raw", TRUE, FALSE,
    "EU_Bound", "DP_R_EUBound", FALSE, TRUE,
    "Path_PA", "DP_R_PA", FALSE, FALSE,
    "Path_TaxaInfo_RData", "DP_R_Taxa_info_rdata", FALSE, TRUE,
    "Path_TaxaInfo", "DP_R_Taxa_info", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)

  rm(EnvVars2Read, envir = environment())
  invisible(gc())

  # Summary information for the current species
  SpInfo <- ecokit::load_as(Path_TaxaInfo_RData)
  if (nrow(SpInfo) == 0) {
    ecokit::stop_ctx(
      "Species information could not be found",
      Path_TaxaInfo_RData = Path_TaxaInfo_RData,
      SpInfo = SpInfo, nrow_SpInfo = nrow(SpInfo), include_backtrace = TRUE)
  }
  if (!species %in% SpInfo$Species_name) {
    ecokit::stop_ctx(
      "Species name not found in the list of species",
      species = species, include_backtrace = TRUE)
  }

  SpInfo <- dplyr::filter(SpInfo, Species_name2 == Species2) %>%
    dplyr::select(-"speciesKey") %>%
    dplyr::distinct()

  # # ..................................................................... ###

  # Check / create directories
  Path_Summary <- fs::path(Path_PA, "PA_summary")
  path_JPEG <- fs::path(Path_PA, "Distribution_JPEG")
  if (!fs::dir_exists(path_JPEG)) {
    fs::dir_create(path_JPEG)
  }
  out_path <- fs::path(path_JPEG, paste0(SpInfo$Species_File[1], ".jpeg"))

  # Loading species data
  selected_columns <- c(
    "NCells_All", "GBIF_Gr100", "EASIN_Gr100", "eLTER_Gr100",
    "BioRegsSumm_N", "BioRegsSumm_Min", "BioRegsSumm_Max", "BioRegsSumm_Mean",
    "species_ID", "PA_Map", "Countries2Exclude",
    "GBIF", "GBIF_Unique", "EASIN", "EASIN_Unique", "eLTER", "eLTER_Unique",
    "NCells_Naturalized", "GBIF_Masked", "GBIF_Masked_Unique",
    "EASIN_Masked", "EASIN_Masked_Unique", "eLTER_Masked",
    "eLTER_Masked_Unique")

  SpData <- fs::path(Path_Summary, paste0(SpInfo$Species_File[1], ".qs2"))
  if (!file.exists(SpData)) {
    return(invisible(NULL))
  }
  SpData <- ecokit::load_as(SpData) %>%
    dplyr::select(tidyselect::all_of(selected_columns))
  invisible(gc())

  if (SpData$NCells_All == 0) {
    return(invisible(NULL))
  }

  CountryBound <- ecokit::load_as(EU_Bound) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_10")
  LastUpdate <- paste0("Last update: ", format(Sys.Date(), "%d %B %Y"))


  Grid_100_sf <- fs::path(Path_Grid_Ref, "Grid_100_sf.RData") %>%
    ecokit::load_as() %>%
    magrittr::extract2("Grid_100_sf_s")

  # Location of legend
  Legend_GBIF <- dplyr::filter(Grid_100_sf, CellCode == "100kmE27N45")
  Legend_EASIN <- dplyr::filter(Grid_100_sf, CellCode == "100kmE27N44")
  Legend_eLTER <- dplyr::filter(Grid_100_sf, CellCode == "100kmE27N43")

  rm(SpInfo, Grid_100_sf, envir = environment())
  invisible(gc())

  # the study area as simple feature object for plotting
  Grid10_Sf <- fs::path(Path_Grid, "Grid_10_Land_Crop_sf.RData") %>%
    ecokit::load_as()

  GBIF_Gr100 <- SpData$GBIF_Gr100[[1]]
  EASIN_Gr100 <- SpData$EASIN_Gr100[[1]]
  eLTER_Gr100 <- SpData$eLTER_Gr100[[1]]

  SpData <- dplyr::select(
    SpData, -tidyselect::all_of(c("GBIF_Gr100", "EASIN_Gr100", "eLTER_Gr100")))
  invisible(gc())

  # # |||||||||||||||||||||||||||||||||||||||||||||||||

  # Biogeographical regions ----
  BioRegAnnotation <- paste0(
    "Observed from ", SpData$BioRegsSumm_N,
    " biogeographical regions\n# presence grid cells per region ranges from ",
    SpData$BioRegsSumm_Min, " to ", SpData$BioRegsSumm_Max, "; mean: ",
    SpData$BioRegsSumm_Mean)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||

  # Figure title ----

  IAS_ID <- unique(SpData$species_ID)
  MapTitle <- readr::read_tsv(
    file = Path_TaxaInfo, show_col_types = FALSE, progress = FALSE) %>%
    dplyr::filter(Species_name2 == Species2) %>%
    dplyr::select(tidyselect::all_of(c("Class", "Order", "Family"))) %>%
    unlist() %>%
    stringr::str_c(collapse = " / ") %>%
    paste0(
      "   <span style='font-size: 14pt; color:blue;'><b><i>", IAS_ID,
      " \u2014 ", species, "</i></b></span>   (", ., ")")

  rm(Path_TaxaInfo, envir = environment())
  invisible(gc())

  # # |||||||||||||||||||||||||||||||||||||||||||||||||

  # Plotting theme ----

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
  sp_PA <- terra::unwrap(SpData$PA_Map[[1]]) %>%
    terra::classify(cbind(0, NA)) %>%
    terra::as.factor() %>%
    as.data.frame(xy = TRUE) %>%
    stats::setNames(c("x", "y", "species"))
  SpData$PA_Map <- NULL

  BoundExclude <- ecokit::load_as(EU_Bound) %>%
    magrittr::extract2("Bound_sf_Eur") %>%
    magrittr::extract2("L_10") %>%
    dplyr::filter(NAME_ENGL %in% SpData$Countries2Exclude[[1]]) %>%
    dplyr::select("NAME_ENGL")

  rm(EU_Bound, envir = environment())
  invisible(gc())

  # # |||||||||||||||||||||||||||||||||||||||||||||||||

  # Figure subtitle ----

  NGrids_All <- paste0(
    "<span style='font-size: 12pt; color:red;'><b>All data:</b></span> ",
    scales::label_comma()(SpData$NCells_All),
    " presence grid cells</b> \u2014 <b>GBIF</b> (",
    scales::label_comma()(SpData$GBIF), " / ",
    scales::label_comma()(SpData$GBIF_Unique), ") \u2014 <b>EASIN</b> (",
    scales::label_comma()(SpData$EASIN), " / ",
    scales::label_comma()(SpData$EASIN_Unique), ") \u2014 <b>eLTER</b> (",
    scales::label_comma()(SpData$eLTER), " / ",
    scales::label_comma()(SpData$eLTER_Unique), ")<br>",
    "<span style='font-size: 12pt; color:red;'><b>Final data: </span>",
    scales::label_comma()(SpData$NCells_Naturalized),
    " presence grid cells</b> \u2014 <b>GBIF</b> (",
    scales::label_comma()(SpData$GBIF_Masked), " / ",
    scales::label_comma()(SpData$GBIF_Masked_Unique),
    ") \u2014 <b>EASIN</b> (",
    scales::label_comma()(SpData$EASIN_Masked), " / ",
    scales::label_comma()(SpData$EASIN_Masked_Unique),
    ") \u2014 <b>eLTER</b> (",
    scales::label_comma()(SpData$eLTER_Masked), " / ",
    scales::label_comma()(SpData$eLTER_Masked_Unique), ")\n")

  # # |||||||||||||||||||||||||||||||||||||||||||||||||

  # ggplot object -----

  Plot <- ggplot2::ggplot() +
    # country boundaries
    ggplot2::geom_sf(
      CountryBound, mapping = ggplot2::aes(), color = "grey",
      linewidth = 0.2, fill = scales::alpha("grey", 0.2)) +
    # the study area
    ggplot2::geom_sf(
      Grid10_Sf, mapping = ggplot2::aes(), color = "lightgrey",
      fill = "lightgrey", linewidth = 0.15)

  if (nrow(BoundExclude) > 0) {
    Plot <- Plot +
      # Countries to exclude
      ggplot2::geom_sf(
        BoundExclude, mapping = ggplot2::aes(), color = "grey90",
        linewidth = 0.1, fill = scales::alpha("red", 0.2))
  }

  # presence grids at 100 km resolution
  if (nrow(GBIF_Gr100) > 0) {
    Plot <- Plot +
      ggplot2::geom_sf(
      GBIF_Gr100, mapping = ggplot2::aes(), color = "transparent",
      fill = scales::alpha("blue", 0.2), linewidth = 0.15)
  }

  # PA grid at 10 km resolution
  Plot <- Plot +
    ggplot2::geom_tile(
    data = sp_PA, mapping = ggplot2::aes(x = x, y = y, fill = species)) +
    ggplot2::scale_fill_manual(
      values = c("blue", "transparent"), na.value = "transparent")

  if (nrow(EASIN_Gr100) > 0) {
    Plot <- Plot +
      ggplot2::geom_sf(
        EASIN_Gr100, mapping = ggplot2::aes(), color = "red",
        fill = "transparent", linewidth = 0.65)
  }
  if (nrow(eLTER_Gr100) > 0) {
    Plot <- Plot +
      ggplot2::geom_sf(
        eLTER_Gr100, mapping = ggplot2::aes(), color = "darkgreen",
        fill = "transparent", linewidth = 0.65, linetype = "dotdash")
  }

  Plot <- Plot +
    # country boundaries
    ggplot2::geom_sf(
      CountryBound, mapping = ggplot2::aes(), color = "black",
      linewidth = 0.2, fill = "transparent") +
    # legends for data type
    ggplot2::geom_sf(
      Legend_GBIF, mapping = ggplot2::aes(), color = "transparent",
      fill = scales::alpha("blue", 0.2), linewidth = 0.15) +
    ggplot2::geom_sf(
      Legend_EASIN, mapping = ggplot2::aes(), color = "red",
      fill = "transparent", linewidth = 0.65) +
    ggplot2::geom_sf(
      Legend_eLTER, mapping = ggplot2::aes(), color = "darkgreen",
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
    ggplot2::labs(title = MapTitle, subtitle = NGrids_All, fill = NULL) +
    plot_theme

  rm(
    SpData, Grid10_Sf, GBIF_Gr100, EASIN_Gr100, eLTER_Gr100,
    BoundExclude, Legend_GBIF, Legend_EASIN, Legend_eLTER,
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
  print(Plot)
  grid::grid.text(
    label = LastUpdate, x = 0.98, y = 0.975, hjust = 1, vjust = 1,
    gp = grid::gpar(col = "grey65", fontsize = 12))
  grid::grid.text(
    label = BioRegAnnotation, x = 0.02, y = 0.9, hjust = 0, vjust = 0,
    gp = grid::gpar(col = "grey65", fontsize = 10, lineheight = 0.7))

  grDevices::dev.off()

  return(out_path)
}
