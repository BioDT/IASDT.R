# # |------------------------------------------------------------------------| #
# Railway_Intensity ----
## |------------------------------------------------------------------------| #

#' Calculate railway intensity based on `OpenStreetMap` data
#'
#' This function downloads, processes, and analyzes railway data extracted from
#' [OpenRailwayMap](https://www.openrailwaymap.org) available from
#' [OpenStreetMap Data Extracts](https://download.geofabrik.de/). It supports
#' parallel processing for faster execution and can calculate the total length
#' of railways and distance to the nearest railway for each grid cell in Europe.
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param NCores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 8.
#' @param DeleteProcessed Logical indicating whether to delete the raw
#'   downloaded railways files after processing them. This helps to free large
#'   unnecessary file space (> 55 GB). Defaults to `TRUE`.
#' @return `NULL`. Outputs processed files to the directories specified in the
#'   environment file.
#' @name Railway_Intensity
#' @export
#' @author Ahmed El-Gabbas

Railway_Intensity <- function(
    EnvFile = ".env", NCores = 6L, DeleteProcessed = TRUE) {

  .StartTime <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Checking arguments ----
  IASDT.R::CatTime("Checking arguments")

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "character", Args = "EnvFile")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "logical", Args = "Download")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "numeric", Args = "NCores")

  rm(AllArgs, envir = environment())

  # # ..................................................................... ###

  # Check `unzip` system command
  if (isFALSE(IASDT.R::CheckCommands("unzip"))) {
    stop("The system command 'unzip' is not available", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message

  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_Railways <- Path_Railways_Raw <- Path_Railways_Interim <- RefGrid <-
    Country <- URL <- URL2 <- Path <- EU_Bound <- fclass <- Path_Grid <-
    CellCode <- Name <- OldName <- NewName <- DT <- Railways_URL <- NULL

  # # ..................................................................... ###

  # Environment variables ----
  IASDT.R::CatTime("Environment variables")

  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_Railways", "DP_R_Railways_processed", FALSE, FALSE,
    "Path_Railways_Raw", "DP_R_Railways_raw", FALSE, FALSE,
    "Path_Railways_Interim", "DP_R_Railways_interim", FALSE, FALSE,
    "Path_Grid", "DP_R_Grid_processed", TRUE, FALSE,
    "EU_Bound", "DP_R_EUBound", FALSE, TRUE,
    "Railways_URL", "DP_R_Railways_url", FALSE, FALSE)
  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  fs::dir_create(
    c(Path_Railways, Path_Railways_Raw, Path_Railways_Interim)
  )

  RefGrid <- IASDT.R::Path(Path_Grid, "Grid_10_Land_Crop.RData")
  if (!file.exists(RefGrid)) {
    stop("The reference grid file does not exist: ", RefGrid, call. = FALSE)
  }
  RefGrid <- terra::unwrap(IASDT.R::LoadAs(RefGrid))


  RefGridSF <- IASDT.R::Path(Path_Grid, "Grid_10_Land_Crop_sf.RData")
  if (!file.exists(RefGridSF)) {
    stop("The reference grid file does not exist: ", RefGridSF, call. = FALSE)
  }
  RefGridSF <- IASDT.R::LoadAs(RefGridSF)

  # # ..................................................................... ###

  # Scrap download links -----
  IASDT.R::CatTime("Scrap download links")
  .StartTimeDown <- lubridate::now(tzone = "CET")

  if (!IASDT.R::CheckURL(Railways_URL)) {
    stop(
      "The base URL for railways data is not valid: ", Railways_URL,
      call. = FALSE)
  }

  IASDT.R::CatTime(paste0("Base URL is: ", Railways_URL), Level = 1)

  # We download European data at country level. For most countries, the data are
  # available in single file, while for others the data are divided into
  # sub-regions. Data on 3 federal states in Germany are not available in single
  # link but at one level below state

  German_L3 <- dplyr::tribble(
    ~URL, ~Country,
    "europe/germany/baden-wuerttemberg.html", "Germany",
    "europe/germany/bayern.html", "Germany",
    "europe/germany/nordrhein-westfalen.html", "Germany") %>%
    dplyr::mutate(URL = paste0(Railways_URL, URL))

  # scrap initial railways links
  Attempt <- 1
  Attempts <- 5
  Success <- FALSE
  while (isFALSE(Success) && Attempt <= Attempts) {
    Railways_Links <- try(
      expr = {
        paste0(Railways_URL, "europe.html") %>%
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
              .f = ~ {
                stringr::str_remove_all(.x, "europe/|.html") %>%
                  stringr::str_to_title()
              }
            ),
            URL = paste0(Railways_URL, URL)
          ) %>%
          dplyr::filter(!(Country %in% c("Russia", "Turkey", "Ukraine"))) %>%
          dplyr::bind_rows(German_L3) %>%
          dplyr::arrange(Country, URL)
      }, silent = TRUE)

    if (inherits(Railways_Links, "tibble")) {
      Success <- TRUE
    } else if (inherits(Railways_Links, "try-error")) {
      Success <- FALSE
      if (Attempt == Attempts) {
        stop(
          "Initial scraping of railways links failed after ", Attempts,
          call. = FALSE)
      }
    }

    Attempt <- Attempt + 1
  }

  # Scraping final railways links
  Railways_Links <- Railways_Links %>%
    dplyr::mutate(
      URL2 = purrr::map(
        .x = URL,
        .f = ~ {
          Success <- FALSE
          Attempt <- 1
          while (isFALSE(Success) && Attempt <= Attempts) {
            ScrapedLinks <- try(
              expr = {
                BaseURL2 <- stringr::str_extract(.x, "^.+/")
                rvest::session(.x) %>%
                  rvest::html_elements(css = "a") %>%
                  rvest::html_attr(name = "href") %>%
                  stringr::str_subset("-latest-free.shp.zip$") %>%
                  paste0(BaseURL2, .) %>%
                  unique()
              }, silent = TRUE)

            if (inherits(ScrapedLinks, "character")) {
              Success <- TRUE
            } else if (inherits(ScrapedLinks, "try-error")) {
              Success <- FALSE
              if (Attempt == Attempts) {
                stop(
                  "Scraping railways links failed for: ", ScrapedLinks,
                  "after ", Attempts, call. = FALSE)
              }
            }
            Attempt <- Attempt + 1
          }

          return(ScrapedLinks)
        }
      )
    ) %>%
    tidyr::unnest(cols = "URL2") %>%
    dplyr::mutate(
      # Download path
      Path = purrr::map2(
        .x = URL2, .y = Country,
        .f = ~ {
          stringr::str_remove_all(.x, "^.+/|-latest-free.shp") %>%
            paste0(.y, "_", .) %>%
            IASDT.R::Path(Path_Railways_Raw, .)
        }
      ),
      ModDate = purrr::map(
        .x = URL2,
        .f = ~ {
          # Last modified date
          paste0("curl -sI ", .x) %>%
            system(intern = TRUE) %>%
            stringr::str_subset("Last-Modified") %>%
            stringr::str_remove_all("Last-Modified: ") %>%
            readr::parse_datetime(format = "%a, %d %b %Y %H:%M:%S %Z") %>%
            lubridate::as_date()
        }
      )
    ) %>%
    tidyr::unnest(c = "ModDate")

  IASDT.R::CatTime(
    paste0("There are ", nrow(Railways_Links), " files to be downloaded"),
    Level = 1)

  save(
    Railways_Links,
    file = IASDT.R::Path(Path_Railways, "Railways_Links.RData"))

  IASDT.R::CatDiff(
    InitTime = .StartTimeDown, NLines = 1, Level = 1,
    Prefix = "Preparing railways download links took ")

  # # ..................................................................... ###

  # Processing railway data ------
  IASDT.R::CatTime("Processing railway data")
  .StartTimeProcess <- lubridate::now(tzone = "CET")

  ## Prepare working on parallel ----
  IASDT.R::CatTime(
    paste0("Prepare working on parallel using `", NCores, "` cores."),
    Level = 1)

  if (NCores == 1) {
    future::plan("future::sequential", gc = TRUE)
  } else {
    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
      future.seed = TRUE)
    c1 <- snow::makeSOCKcluster(NCores)
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("future::cluster", workers = c1, gc = TRUE)
    on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
  }

  ## Processing railway data ----
  IASDT.R::CatTime("Processing railway data", Level = 1)

  Railways_3035 <- future.apply::future_lapply(
    X = seq_len(nrow(Railways_Links)),
    FUN = function(ID) {
      URL <- Railways_Links$URL2[[ID]]
      Path <- Railways_Links$Path[[ID]]
      Country <- Railways_Links$Country[[ID]]
      Prefix <- stringr::str_remove_all(basename(Path), ".zip$")
      Path_Temp <- IASDT.R::Path(
        Path_Railways_Interim, paste0(Prefix, ".RData"))

      withr::local_options(
        future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
        future.seed = TRUE, timeout = 1800)

      # Check if zip file is a valid file
      if (file.exists(Path)) {
        Success <- IASDT.R::CheckZip(Path)
        if (isFALSE(Success)) {
          fs::file_delete(Path)
        }
      } else {
        Success <- FALSE
      }

      # Try downloading data for a max of 5 attempts
      Attempt <- 1
      Attempts <- 5

      while (isFALSE(Success) && Attempt <= Attempts) {
        Down <- try(
          expr = {

            stringr::str_glue(
              'curl -k -L --connect-timeout 120 --max-time 1800 --retry 5 \\
            "{URL}" -o "{Path}" --silent') %>%
              system()

            Success <- IASDT.R::CheckZip(Path)
            
            if (isFALSE(Success)) {
              fs::file_delete(Path)
            }

            Success
          }, silent = TRUE)

        if (inherits(Down, "try-error")) {
          Success <- FALSE
        }
        Attempt <- Attempt + 1
      }

      # Filter only railways files
      InFileN <- dplyr::tibble(utils::unzip(Path, list = TRUE)) %>%
        dplyr::filter(stringr::str_detect(Name, "railways")) %>%
        dplyr::pull(Name) %>%
        unique()

      utils::unzip(
        zipfile = Path, files = InFileN,
        exdir = IASDT.R::Path(Path_Railways_Interim, Prefix))

      Path_Extract <- dplyr::tibble(
        OldName = IASDT.R::Path(Path_Railways_Interim, Prefix, InFileN),
        NewName = IASDT.R::Path(
          Path_Railways_Interim,
          paste0(Prefix, ".", tools::file_ext(InFileN)))) %>%
        dplyr::mutate(Ren = purrr::map2(OldName, NewName, file.rename))

      Railway <- dplyr::pull(Path_Extract, NewName) %>%
        stringr::str_subset(".shp$") %>%
        sf::st_read(quiet = TRUE) %>%
        tibble::tibble() %>%
        sf::st_as_sf() %>%
        sf::st_transform(crs = 3035) %>%
        dplyr::select(
          -tidyselect::all_of(c("bridge", "tunnel", "layer", "name"))) %>%
        sf::st_filter(y = RefGridSF, .predicate = sf::st_intersects)

      save(Railway, file = Path_Temp)

      # Clean up
      fs::dir_delete(IASDT.R::Path(Path_Railways_Interim, Prefix))
      fs::file_delete(Path_Extract$NewName)
      if (DeleteProcessed) {
        fs::file_delete(Path)
      }

      return(
        tibble::tibble(
          URL = URL, Country = Country, Area = Prefix, Path = Path_Temp))
    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = c("dplyr", "fs", "sf", "IASDT.R", "stringr"),
    future.globals = c("Railways_Links", "RefGridSF")
  ) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(DT = purrr::map(Path, IASDT.R::LoadAs)) %>%
    tidyr::unnest(DT) %>%
    sf::st_as_sf()

  if (NCores > 1) {
    snow::stopCluster(c1)
    future::plan("future::sequential", gc = TRUE)
  }

  # Delete raw directory if empty
  RawFileList <- list.files(Path_Railways_Raw, full.names = TRUE)
  if (length(RawFileList) == 0) {
    fs::dir_delete(Path_Railways_Raw)
  }

  # # .................................... ###

  ## Saving - RData -----
  IASDT.R::CatTime("Saving - RData", Level = 1)
  save(
    Railways_3035,
    file = IASDT.R::Path(Path_Railways, "Railways_3035.RData"))

  # # .................................... ###

  ## Railways_3035_2plot ----
  Railways_3035_2plot <- dplyr::filter(Railways_3035, fclass == "rail") %>%
    sf::st_join(RefGridSF) %>%
    dplyr::filter(!is.na(CellCode)) %>%
    dplyr::select("geometry")

  save(
    Railways_3035_2plot,
    file = IASDT.R::Path(Path_Railways, "Railways_3035_2plot.RData"))

  # # .................................... ###

  ## Saving each railway class to separate file ----
  IASDT.R::CatTime("Saving each railway class to separate file", Level = 1)

  sf::st_drop_geometry(Railways_3035) %>%
    dplyr::distinct(fclass) %>%
    dplyr::pull(fclass) %>%
    purrr::walk(
      .f = ~ {
        IASDT.R::CatTime(.x, Level = 2)
        dplyr::filter(Railways_3035, fclass == .x) %>%
          IASDT.R::SaveAs(
            OutObj = paste0("Railways_sf_", .x),
            OutPath = IASDT.R::Path(
              Path_Railways, paste0("Railways_sf_", .x, ".RData")))
      }
    )

  rm(Railways_3035, RefGridSF, envir = environment())
  invisible(gc())

  IASDT.R::CatDiff(
    InitTime = .StartTimeProcess,
    Prefix = "Processing railway data took ", NLines = 1, Level = 1)

  # # ..................................................................... ###

  # Calculate length of railways for each railway class per grid cell -----
  IASDT.R::CatTime(
    "Calculate length of railways for each railway class per grid cell")

  Railways_Length <- Path_Railways %>%
    list.files(pattern = "^Railways_sf_", full.names = TRUE) %>%
    purrr::map(
      .f = ~ {
        Name <- stringr::str_remove_all(basename(.x), "Railways_sf_|.RData")

        IASDT.R::CatTime(Name, Level = 2)
        IASDT.R::LoadAs(.x) %>%
          terra::vect() %>%
          terra::rasterizeGeom(y = RefGrid, fun = "length", unit = "km") %>%
          terra::mask(mask = RefGrid) %>%
          stats::setNames(Name) %>%
          IASDT.R::setRastVals()
      }, .progress = FALSE
    ) %>%
    terra::rast()

  # Sum of railways length of any type
  Railways_Length$Sum <- sum(Railways_Length)

  ## Saving - RData -----
  IASDT.R::CatTime("Saving - RData", Level = 1)
  IASDT.R::SaveAs(
    InObj = terra::wrap(Railways_Length), OutObj = "Railways_Length",
    OutPath = IASDT.R::Path(Path_Railways, "Railways_Length.RData"))

  ## Saving - tif ------
  IASDT.R::CatTime("Saving - tif", Level = 1)
  terra::writeRaster(
    x = Railways_Length, overwrite = TRUE,
    filename = IASDT.R::Path(
      Path_Railways,
      paste0("Railways_Length_", names(Railways_Length), ".tif")))

  # # ..................................................................... ###

  # Calculate distance to rail ------
  IASDT.R::CatTime("Calculate distance to rail")

  # This calculates the distance from each grid cell to the nearest grid cell
  # overlapping with a railway. This can be different than calculating the
  # actual distance to nearest railway line; which is expected to take too much
  # time to calculate

  Railways_Distance <- purrr::map(
    .x = as.list(Railways_Length),
    .f = ~ {
      IASDT.R::CatTime(names(.x), Level = 2)

      # suppress progress bar
      terra::terraOptions(progress = 0)

      Railways_Points <- terra::as.points(terra::classify(.x, cbind(0, NA)))

      terra::distance(x = .x, y = Railways_Points, unit = "km") %>%
        terra::mask(RefGrid) %>%
        stats::setNames(paste0("Railways_Distance_", names(.x))) %>%
        # Ensure that values are read from memory
        IASDT.R::setRastVals()
    }
  ) %>%
    terra::rast()

  IASDT.R::CatTime("Save distance to railways - tif", Level = 1)
  terra::writeRaster(
    x = Railways_Distance, overwrite = TRUE,
    filename = IASDT.R::Path(
      Path_Railways, paste0(names(Railways_Distance), ".tif")))

  IASDT.R::CatTime("Save distance to railways - RData", Level = 1)
  IASDT.R::SaveAs(
    InObj = terra::wrap(Railways_Distance), OutObj = "Railways_Distance",
    OutPath = IASDT.R::Path(Path_Railways, "Railways_Distance.RData"))

  rm(RefGrid, envir = environment())

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
      legend.position = "inside",
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

  # Using ggplot2::ggsave directly does not show non-ascii characters
  # correctly
  grDevices::jpeg(
    filename = IASDT.R::Path(Path_Railways, "Railways_Length.jpeg"),
    width = 31, height = 30, units = "cm", quality = 100, res = 600)
  print(RailPlot)
  grDevices::dev.off()

  rm(RailPlot, envir = environment())

  # # .................................. ###

  ## Plotting railways -----
  IASDT.R::CatTime("Plotting European railways", Level = 1)

  RailPlotShp <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      EU_Bound, mapping = ggplot2::aes(), color = "grey75",
      linewidth = 0.075, fill = "grey98") +
    ggplot2::geom_sf(
      Railways_3035_2plot,
      mapping = ggplot2::aes(), color = "blue", linewidth = 0.05) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)),
      limits = c(2600000, 6700000)) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)),
      limits = c(1450000, 5420000)) +
    ggplot2::labs(title = "Railways in Europe", fill = NA) +
    PlottingTheme

  # Using ggplot2::ggsave directly does not show non-ascii characters
  # correctly
  grDevices::jpeg(
    filename = IASDT.R::Path(Path_Railways, "Railways_Lines.jpeg"),
    width = 31, height = 30, units = "cm", quality = 100, res = 600)
  print(RailPlotShp)
  grDevices::dev.off()

  rm(RailPlotShp, envir = environment())

  # # ..................................................................... ###

  # Function Summary ----

  IASDT.R::CatDiff(
    InitTime = .StartTime,
    Prefix = "\nProcessing railway data was finished in ", ... = "\n")

  return(invisible(NULL))
}
