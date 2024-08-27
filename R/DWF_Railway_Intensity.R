# # |------------------------------------------------------------------------| #
# Railway_Intensity ----
## |------------------------------------------------------------------------| #

#' Calculate Railway Intensity Based on `OpenStreetMap` Data
#'
#' This function downloads, processes, and analyzes railway data extracted from
#' [OpenRailwayMap](https://www.openrailwaymap.org) available from
#' [OpenStreetMap Data Extracts](https://download.geofabrik.de/). It supports
#' parallel processing for faster execution and can calculate the total length
#' of railways and distance to the nearest railway for each grid cell in Europe.
#' @param FromHPC Logical indicating whether the work is being done from HPC, to
#'   adjust file paths accordingly. Default: `TRUE`.
#' @param EnvFile Character. The path to the environment file containing
#'   variables required by the function. Default is ".env".
#' @param Download Logical. If `TRUE` (default), downloads railway data from the
#'   specified source.
#' @param Extract Logical. If `TRUE` (default), extracts railway data.
#' @param GetDownLinks Logical. If `TRUE`, fetches download links for railway
#'   data. Default is `TRUE`.
#' @param CheckZip Logical. If `TRUE` (default), validates the integrity of
#'   downloaded ZIP files.
#' @param NCores Numeric. Number of CPU cores to use for parallel processing.
#'   Default is 6.
#' @return `NULL`. Outputs processed files to the directories specified in the
#'   environment file.
#' @name Railway_Intensity
#' @export
#' @author Ahmed El-Gabbas

Railway_Intensity <- function(
    FromHPC = TRUE, EnvFile = ".env", Download = TRUE, Extract = TRUE,
    GetDownLinks = TRUE, CheckZip = TRUE, NCores = 6) {

  .StartTime <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Checking arguments ----
  IASDT.R::CatTime("Checking arguments")

  AllArgs <- ls()
  AllArgs <- purrr::map(AllArgs, ~get(.x, envir = environment())) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "character", Args = "EnvFile")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "logical",
    Args = c("FromHPC", "Download", "GetDownLinks", "CheckZip"))
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "numeric", Args = "NCores")

  rm(AllArgs)

  # # ..................................................................... ###

  IASDT.R::CatTime("Check system commands")
  IASDT.R::CheckCommands("unzip")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Railways <- Path_Railways_Raw <- Path_Railways_Interim <- RefGrid <-
    Country <- URL <- URL2 <- Path <- Shp <- bridge <- tunnel <- EU_Bound <-
    fclass <- Path_Grid <- CellCode <- NULL

  # # ..................................................................... ###

  # Environment variables ----
  IASDT.R::CatTime("Environment variables")

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Railways", "DP_R_Railways", FALSE, FALSE,
      "Path_Railways_Raw", "DP_R_Railways_Raw", FALSE, FALSE,
      "Path_Railways_Interim", "DP_R_Railways_Interim", FALSE, FALSE,
      "Path_Grid", "DP_R_Grid", TRUE, FALSE,
      "EU_Bound", "DP_R_EUBound_sf", FALSE, TRUE,
      "Railways_URL", "DP_R_Railways_URL", FALSE, FALSE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_Railways", "DP_R_Railways_Local", FALSE, FALSE,
      "Path_Railways_Raw", "DP_R_Railways_Raw_Local", FALSE, FALSE,
      "Path_Railways_Interim", "DP_R_Railways_Interim_Local", FALSE, FALSE,
      "Path_Grid", "DP_R_Grid_Local", TRUE, FALSE,
      "EU_Bound", "DP_R_EUBound_sf_Local", FALSE, TRUE,
      "Railways_URL", "DP_R_Railways_URL", FALSE, FALSE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  fs::dir_create(
    c(Path_Railways, Path_Railways_Raw, Path_Railways_Interim))

  RefGrid <- file.path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(RefGrid)) {
    stop(
      paste0("The reference grid file does not exist: ", RefGrid),
      call. = FALSE)
  }

  withr::local_options(
    list(timeout = 1200, future.globals.maxSize = 8000 * 1024^2))

  # # ..................................................................... ###

  # Prepare working on parallel ----
  IASDT.R::CatTime("Prepare working on parallel")
  c1 <- snow::makeSOCKcluster(NCores)
  on.exit({
    invisible(try(snow::stopCluster(c1), silent = TRUE))
    future::plan(future::sequential, gc = TRUE)
  }, add = TRUE)
  future::plan(future::cluster, workers = c1, gc = TRUE)

  # # ..................................................................... ###

  # Prepare download links -----

  if (GetDownLinks) {

    IASDT.R::CatTime("Prepare download links")
    .StartTimeDown <- lubridate::now(tzone = "CET")

    # OpenStreetMap Data Extracts
    BaseURL <- "https://download.geofabrik.de/"

    # We download European data at country level. For most countries, the data
    # are available in single file, while for others the data are divided into
    # sub-regions. Data on 3 federal states in Germany are not available in
    # single link but at one level below state

    German_L3 <- dplyr::tribble(
      ~URL, ~Country,
      "europe/germany/baden-wuerttemberg.html", "Germany",
      "europe/germany/bayern.html", "Germany",
      "europe/germany/nordrhein-westfalen.html", "Germany") %>%
      dplyr::mutate(URL = paste0(BaseURL, URL))

    Railways_Links <- paste0(BaseURL, "europe.html") %>%
      rvest::session() %>%
      rvest::html_elements(css = "table") %>%
      magrittr::extract(2) %>%
      rvest::html_elements(css = "a") %>%
      rvest::html_attr("href") %>%
      stringr::str_subset(".html$") %>%
      dplyr::tibble(URL = .) %>%
      dplyr::mutate(
        Country = purrr::map_chr(
          .x = URL,
          .f = ~{
            stringr::str_remove_all(.x, "europe/|.html") %>%
              stringr::str_to_title()
          }),
        URL = paste0(BaseURL, URL)) %>%
      dplyr::bind_rows(German_L3) %>%
      dplyr::arrange(Country, URL) %>%
      dplyr::mutate(
        URL2 = purrr::map(
          .x = URL,
          .f = ~{
            BaseURL2 <- stringr::str_extract(.x, "^.+/")
            rvest::session(.x) %>%
              rvest::html_elements(css = "a") %>%
              rvest::html_attr(name = "href") %>%
              stringr::str_subset("-latest-free.shp.zip$") %>%
              paste0(BaseURL2, .) %>%
              unique()
          })) %>%
      tidyr::unnest(cols = "URL2") %>%
      dplyr::mutate(
        Path = purrr::map2(
          .x = URL2, .y = Country,
          .f = ~{
            stringr::str_remove_all(.x, "^.+/|-latest-free.shp") %>%
              paste0(.y, "_", .) %>%
              file.path(Path_Railways_Raw, .)
          }),
        ModDate = purrr::map(
          .x = URL2,
          .f = ~{
            system(paste0("curl -sI ", .x), intern = TRUE) %>%
              stringr::str_subset("Last-Modified") %>%
              stringr::str_remove_all("Last-Modified: ") %>%
              readr::parse_datetime(format = "%a, %d %b %Y %H:%M:%S %Z") %>%
              lubridate::as_date()
          })) %>%
      tidyr::unnest(c = "ModDate")

    IASDT.R::CatTime(
      paste0("There are ", nrow(Railways_Links), " files to be downloaded"),
      Level = 1)

    save(
      Railways_Links, file = file.path(Path_Railways, "Railways_Links.RData"))


    IASDT.R::CatDiff(
      InitTime = .StartTimeDown, CatInfo = FALSE,
      Prefix = "Preparing railways download links took ", NLines = 1, Level = 1)

  } else {
    IASDT.R::CatTime("Loading download links")
    Railways_Links <- IASDT.R::LoadAs(
      file.path(Path_Railways, "Railways_Links.RData"))
  }

  # # ..................................................................... ###

  # Download railway data ------
  IASDT.R::CatTime("Download railway data")

  if (Download) {

    IASDT.R::CatTime("Downloading files", Level = 1)
    .StartTimeDown <- lubridate::now(tzone = "CET")

    invisible(snow::clusterEvalQ(
      cl = c1, IASDT.R::LoadPackages(List = c("dplyr"))))
    snow::clusterExport(
      cl = c1, list = c("Railways_Links", "CheckZip"), envir = environment())

    DownRail <- future.apply::future_lapply(
      X = seq_len(nrow(Railways_Links)),
      FUN = function(ID) {

        URL <- Railways_Links$URL2[[ID]]
        Path <- Railways_Links$Path[[ID]]

        if (file.exists(Path)) {

          # Check if zip file is a valid file
          if (CheckZip) {
            ZipStatus <- system2(
              "unzip", args = c("-t", Path), stdout = TRUE, stderr = TRUE) %>%
              stringr::str_detect("No errors detected in compressed data") %>%
              any()

            if (ZipStatus) {
              Success <- TRUE
            } else {
              Success <- FALSE
              fs::file_delete(Path)
            }

          } else {
            Success <- TRUE
          }

        } else {
          Success <- FALSE
        }

        # Try downloading data for a max of 3 attempts, each with 20 mins time
        # out
        withr::local_options(list(timeout = 1200))

        Attempt <- 1
        Attempts <- 3

        while (!Success && (Attempt <= Attempts)) {
          tryCatch({
            utils::download.file(
              url = URL, destfile = Path, mode = "wb", quiet = TRUE)

            ZipStatus <- system2(
              "unzip", args = c("-t", Path), stdout = TRUE, stderr = TRUE) %>%
              stringr::str_detect("No errors detected in compressed data") %>%
              any()
            Success <- all(ZipStatus, file.exists(Path))
          },
          error = function(e) {
            if (Attempt < Attempts) {
              Attempt <- Attempt + 1
            } else {
              stop(
                paste0(
                  "Failed to download data from ", URL, " after ", Attempts,
                  " attempts: ", conditionMessage(e)),
                call. = FALSE)
            }
          })
        }
        return(invisible(NULL))
      },
      future.scheduling = Inf, future.seed = TRUE)

    rm(DownRail)

    IASDT.R::CatDiff(
      InitTime = .StartTimeDown, CatInfo = FALSE,
      Prefix = "Downloading railway data took ", NLines = 1, Level = 2)
  }

  # # ..................................................................... ###

  # Extract shape files for railways ----
  IASDT.R::CatTime("Extract shape files for railways")
  RefGrid <- terra::unwrap(IASDT.R::LoadAs(RefGrid))

  if (Extract) {
    .StartTimeExtract <- lubridate::now(tzone = "CET")

    ## Extracting / project railways files ----
    IASDT.R::CatTime("Extracting / project railways files", Level = 1)

    Railways_3035 <- dplyr::mutate(
      Railways_Links,
      Shp = furrr::future_map2_chr(
        .x = Path, .y = Country,
        .f = ~{

          Prefix <- stringr::str_remove_all(.x, "^.+/|.zip")

          # Filter only railways files
          InFileN <- dplyr::tibble(unzip(.x, list = TRUE)) %>%
            dplyr::filter(stringr::str_detect(Name, "railways")) %>%
            dplyr::pull(Name) %>%
            unique()

          # Extract selected files
          unzip(zipfile = .x, files = InFileN, exdir = Path_Railways_Interim)

          Path_Shp <- dplyr::tibble(
            OldName = file.path(Path_Railways_Interim, InFileN),
            NewName = file.path(
              Path_Railways_Interim,
              paste0(Prefix, ".", tools::file_ext(InFileN)))) %>%
            dplyr::mutate(Ren = purrr::map2(OldName, NewName, file.rename)) %>%
            dplyr::pull(NewName) %>%
            stringr::str_subset(".shp$")

          return(Path_Shp)
        }, .options = furrr::furrr_options(seed = TRUE))) %>%
      dplyr::mutate(
        Rail = purrr::map(
          .x = Shp,
          .f = ~sf::st_transform(
            sf::st_read(.x, quiet = TRUE), crs = 3035))) %>%
      tidyr::unnest(cols = "Rail") %>%
      dplyr::select(
        -tidyselect::all_of(c("URL", "Path", "Shp", "ModDate", "layer"))) %>%
      sf::st_as_sf() %>%
      dplyr::mutate(bridge = as.logical(bridge), tunnel = as.logical(tunnel))

    IASDT.R::CatDiff(
      InitTime = .StartTimeExtract, CatInfo = FALSE,
      Prefix = "Extracting railway data took ", NLines = 1, Level = 2)

    # # .................................... ###

    ## Saving - RData -----
    IASDT.R::CatTime("Saving - RData", Level = 1)
    save(Railways_3035, file = file.path(Path_Railways, "Railways_3035.RData"))

    # # .................................... ###

    ## Railways_3035_2plot ----
    RefGridSF <- IASDT.R::LoadAs(file.path(Path_Grid, "Grid_10_Land_sf.RData"))
    Railways_3035_2plot <- dplyr::filter(Railways_3035, fclass == "rail") %>%
      sf::st_join(RefGridSF) %>%
      dplyr::filter(!is.na(CellCode)) %>%
      dplyr::select("geometry")
    save(
      Railways_3035_2plot,
      file = file.path(Path_Railways, "Railways_3035_2plot.RData"))

    # # .................................... ###

    ## Each railway class to separate file ----
    IASDT.R::CatTime("Each railway class to separate file", Level = 1)

    sf::st_drop_geometry(Railways_3035) %>%
      dplyr::distinct(fclass) %>%
      dplyr::pull(fclass) %>%
      purrr::walk(
        .f = ~{
          IASDT.R::CatTime(.x, Level = 2)
          dplyr::filter(Railways_3035, fclass == .x) %>%
            IASDT.R::SaveAs(
              OutObj = paste0("Railways_sf_", .x),
              OutPath = file.path(
                Path_Railways, paste0("Railways_sf_", .x, ".RData")))
        })

    rm(Railways_3035, RefGridSF)

  } else {

    Railways_3035_2plot <- IASDT.R::LoadAs(file.path(Path_Railways, "Railways_3035_2plot.RData"))

  }

  snow::stopCluster(c1)
  future::plan(future::sequential, gc = TRUE)
  invisible(gc())

  # # ..................................................................... ###

  # Calculate length of railways -----
  IASDT.R::CatTime("Calculate length of railways")

  Railways_Length <- Path_Railways %>%
    list.files(pattern = "^Railways_sf_", full.names = TRUE) %>%
    purrr::map(
      .f = ~{
        Name <- stringr::str_remove_all(basename(.x), "Railways_sf_|.RData")

        IASDT.R::CatTime(Name, Level = 2)
        IASDT.R::LoadAs(.x) %>%
          terra::vect() %>%
          terra::rasterizeGeom(y = RefGrid, fun = "length", unit = "km") %>%
          terra::mask(mask = RefGrid) %>%
          stats::setNames(Name) %>%
          # Ensure that values are read from memory
          IASDT.R::setRastVals()
      }, .progress = FALSE) %>%
    terra::rast()

  # Sum of railways length of any type
  Railways_Length$Sum <- sum(Railways_Length)

  ## Saving - RData -----
  IASDT.R::CatTime("Saving - RData", Level = 1)
  IASDT.R::SaveAs(
    InObj = terra::wrap(Railways_Length), OutObj = "Railways_Length",
    OutPath = file.path(Path_Railways, "Railways_Length.RData"))

  ## Saving - tif ------
  IASDT.R::CatTime("Saving - tif", Level = 1)
  terra::writeRaster(
    x = Railways_Length, overwrite = TRUE,
    filename = file.path(
      Path_Railways, paste0("Railways_Length_", names(Railways_Length), ".tif")))

  # # ..................................................................... ###

  # Calculate distance to rail ------
  IASDT.R::CatTime("Calculate distance to rail")

  # This calculates the distance from each grid cell to the nearest grid cell
  # overlapping with a railway. This can be different than calculating the
  # actual distance to nearest railway line; which is expected to take too much
  # time to calculate

  Railways_Distance <- purrr::map(
    .x = as.list(Railways_Length),
    .f = ~{

      IASDT.R::CatTime(names(.x), Level = 2)

      Railways_Points <- terra::as.points(terra::classify(.x, cbind(0, NA)))

      terra::distance(x = .x, y = Railways_Points, unit = "km") %>%
        terra::mask(RefGrid) %>%
        stats::setNames(paste0("Railways_Distance_", names(.x))) %>%
        # Ensure that values are read from memory
        IASDT.R::setRastVals()
    }) %>%
    terra::rast()

  IASDT.R::CatTime("Save distance to railways as tif files", Level = 1)
  terra::writeRaster(
    x = Railways_Distance, overwrite = TRUE,
    filename = file.path(Path_Railways, paste0(names(Railways_Distance), ".tif")))

  IASDT.R::CatTime("Save distance to railways as RData", Level = 1)
  IASDT.R::SaveAs(
    InObj = terra::wrap(Railways_Distance), OutObj = "Railways_Distance",
    OutPath = file.path(Path_Railways, "Railways_Distance.RData"))

  # # ..................................................................... ###

  # Plotting ------
  IASDT.R::CatTime("Plotting")

  EU_Bound <- IASDT.R::LoadAs(EU_Bound) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_03")

  PlottingTheme <- ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
      plot.title = ggplot2::element_text(
        size = 20, color = "blue", face = "bold", hjust = 0.5,
        margin = ggplot2::margin(2, 0, 2, 0)),
      strip.text = ggplot2::element_text(size = 5.5, face = "bold"),
      strip.background = ggplot2::element_rect(
        fill = "transparent", colour = "transparent"),
      legend.key.size = grid::unit(0.8, "cm"),
      legend.key.width = grid::unit(0.8, "cm"),
      legend.background = ggplot2::element_rect(fill = "transparent"),
      legend.text = ggplot2::element_text(size = 12),
      legend.position	= "inside",
      legend.position.inside = c(0.94, 0.9),
      legend.title = ggplot2::element_text(
        color = "black", size = 12, face = "bold"),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.tag.position = c(0.94, 0.011),
      plot.tag = ggtext::element_markdown(colour = "grey", size = 4),
      panel.ontop = TRUE,
      panel.spacing = grid::unit(0.05, "lines"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = NA))

  # # .................................. ###

  ## Plotting length of railways -----
  IASDT.R::CatTime("Plotting length of railways", Level = 1)

  Rail2Plot <- terra::subset(Railways_Length, "rail") %>%
    terra::classify(cbind(0, NA))

  RailPlot <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      EU_Bound, mapping = ggplot2::aes(), color = "grey75",
      linewidth = 0.075, fill = "grey98") +
    tidyterra::geom_spatraster(data = log10(Rail2Plot + 1), maxcell = Inf) +
    ggplot2::geom_sf(
      EU_Bound, mapping = ggplot2::aes(), color = "grey30",
      linewidth = 0.075, fill = "transparent", inherit.aes = TRUE) +
    paletteer::scale_fill_paletteer_c(
      na.value = "transparent", "viridis::plasma") +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)),
      limits = c(2600000, 6700000)) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)),
      limits = c(1450000, 5420000)) +
    ggplot2::labs(title = "Railways length", fill = "log10") +
    PlottingTheme

  ggplot2::ggsave(
    filename = file.path(Path_Railways, "Railways_Length.jpeg"),
    plot = RailPlot, width = 31, height = 30, units = "cm", dpi = 600)

  # # .................................. ###

  ## Plotting railways -----

  RailPlotShp <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      EU_Bound, mapping = ggplot2::aes(), color = "grey75",
      linewidth = 0.075, fill = "grey98") +
    ggplot2::geom_sf(
      Railways_3035_2plot, mapping = ggplot2::aes(), color = "blue",
      linewidth = 0.05) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)),
      limits = c(2600000, 6700000)) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)),
      limits = c(1450000, 5420000)) +
    ggplot2::labs(title = "Railways in Europe", fill = NA) +
    PlottingTheme

  ggplot2::ggsave(
    filename = file.path(Path_Railways, "Railways_Lines.jpeg"),
    plot = RailPlotShp, width = 31, height = 30, units = "cm", dpi = 900)

  # # ..................................................................... ###

  # Function Summary ----

  IASDT.R::CatDiff(
    InitTime = .StartTime, CatInfo = FALSE,
    Prefix = "\nProcessing railway data was finished in ", ... = "\n")

  return(invisible(NULL))
}
