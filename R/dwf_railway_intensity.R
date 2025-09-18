# # |------------------------------------------------------------------------| #
# railway_intensity ----
## |------------------------------------------------------------------------| #

#' Calculate railway intensity based on `OpenStreetMap` data
#'
#' This function downloads, processes, and analyses railway data extracted from
#' [OpenRailwayMap](https://www.openrailwaymap.org) available from
#' [OpenStreetMap Data Extracts](https://download.geofabrik.de/). It supports
#' parallel processing for faster execution and can calculate the total length
#' of railways and distance to the nearest railway for each grid cell in Europe.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 8.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "sequential", "multisession" (default), "multicore", and
#'   "cluster". See [future::plan()] and [ecokit::set_parallel()] for details.
#' @param delete_processed Logical indicating whether to delete the raw
#'   downloaded railways files after processing them. This helps to free large
#'   unnecessary file space (> 55 GB). Defaults to `TRUE`.
#' @return `NULL`. Outputs processed files to the directories specified in the
#'   environment file.
#' @name railway_intensity
#' @export
#' @author Ahmed El-Gabbas

railway_intensity <- function(
    env_file = ".env", n_cores = 6L, strategy = "multisession",
    delete_processed = TRUE) {

  .start_time <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Checking arguments ----
  ecokit::cat_time("Checking arguments")

  ecokit::check_args(args_to_check = "delete_processed", args_type = "logical")

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)

  # # ..................................................................... ###

  # Check `unzip` system command
  if (isFALSE(ecokit::check_system_command("unzip"))) {
    ecokit::stop_ctx(
      "The system command 'unzip' is not available", include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message

  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  path_railway <- path_railway_raw <- path_railway_interim <- grid_ref <-
    country <- url <- url2 <- path <- eu_boundaries <- fclass <- path_grid <-
    CellCode <- old_name <- new_name <- data <- railway_url <- Name <- NULL

  # # ..................................................................... ###

  # Environment variables ----
  ecokit::cat_time("Environment variables")

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_railway", "DP_R_railway_processed", FALSE, FALSE,
    "path_railway_raw", "DP_R_railway_raw", FALSE, FALSE,
    "path_railway_interim", "DP_R_railway_interim", FALSE, FALSE,
    "path_grid", "DP_R_grid_processed", TRUE, FALSE,
    "eu_boundaries", "DP_R_country_boundaries", FALSE, TRUE,
    "railway_url", "DP_R_railway_url", FALSE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  fs::dir_create(c(path_railway, path_railway_raw, path_railway_interim))

  grid_ref <- fs::path(path_grid, "grid_10_land_crop.RData")
  if (!file.exists(grid_ref)) {
    ecokit::stop_ctx(
      "The reference grid file does not exist", grid_ref = grid_ref,
      include_backtrace = TRUE)
  }
  grid_ref <- ecokit::load_as(grid_ref, unwrap_r = TRUE)

  grid_sf <- fs::path(path_grid, "grid_10_land_sf.RData")
  if (!file.exists(grid_sf)) {
    ecokit::stop_ctx(
      "The reference grid file does not exist", grid_sf = grid_sf,
      include_backtrace = TRUE)
  }
  grid_sf <- ecokit::load_as(grid_sf)

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c(
      "dplyr", "fs", "sf", "IASDT.R", "stringr", "withr",
      "tibble", "ecokit", "magrittr"),
    strategy = strategy)

  # # ..................................................................... ###

  # Scrap download links -----
  ecokit::cat_time("Scrap download links")
  .start_time_down <- lubridate::now(tzone = "CET")

  if (!ecokit::check_url(railway_url)) {
    ecokit::stop_ctx(
      "The base URL for railways data is not valid",
      railway_url = railway_url, include_backtrace = TRUE)
  }

  ecokit::cat_time(paste0("Base url is: ", railway_url), level = 1L)

  # We download European data at country level. For most countries, the data are
  # available in single file, while for others the data are divided into
  # sub-regions. Data on 3 federal states in Germany are not available in single
  # link but at one level below state

  german_l3 <- dplyr::tribble(
    ~url, ~country,
    "europe/germany/baden-wuerttemberg.html", "Germany",
    "europe/germany/bayern.html", "Germany",
    "europe/germany/nordrhein-westfalen.html", "Germany") %>%
    dplyr::mutate(url = paste0(railway_url, url))

  # scrap initial railways links
  attempt <- 1
  attempts <- 5
  success <- FALSE
  while (isFALSE(success) && attempt <= attempts) {
    railway_links <- try(
      expr = {
        paste0(railway_url, "europe.html") %>%
          rvest::session() %>%
          rvest::html_elements(css = "table") %>%
          magrittr::extract(2) %>%
          rvest::html_elements(css = "a") %>%
          rvest::html_attr("href") %>%
          stringr::str_subset(".html$") %>%
          stringr::str_remove("^/") %>%
          dplyr::tibble(url = .) %>%
          dplyr::mutate(
            country = purrr::map_chr(
              .x = url,
              .f = ~ {
                stringr::str_remove_all(.x, "europe/|.html") %>%
                  stringr::str_to_title()
              }
            ),
            url = paste0(railway_url, url)) %>%
          dplyr::filter(!(country %in% c("Russia", "Turkey", "Ukraine"))) %>%
          dplyr::bind_rows(german_l3) %>%
          dplyr::arrange(country, url)
      }, silent = TRUE)

    if (inherits(railway_links, "tibble")) {
      success <- TRUE
    } else if (inherits(railway_links, "try-error")) {
      success <- FALSE
      if (attempt == attempts) {
        ecokit::stop_ctx(
          paste0("Initial scraping of railways links failed after ", attempts),
          include_backtrace = TRUE)
      }
    }

    attempt <- attempt + 1
  }

  # Scraping final railways links
  railway_links <- railway_links %>%
    dplyr::mutate(
      url2 = purrr::map(
        .x = url,
        .f = ~ {

          success <- FALSE
          attempt <- 1
          while (isFALSE(success) && attempt <= attempts) {

            scraped_links <- try(
              expr = {
                base_url_2 <- stringr::str_extract(.x, "^.+/")
                rvest::session(.x) %>%
                  rvest::html_elements(css = "a") %>%
                  rvest::html_attr(name = "href") %>%
                  stringr::str_subset("-latest-free.shp.zip$") %>%
                  paste0(base_url_2, .) %>%
                  unique()
              }, silent = TRUE)

            if (inherits(scraped_links, "character")) {
              success <- TRUE
            } else if (inherits(scraped_links, "try-error")) {
              success <- FALSE
              if (attempt == attempts) {
                ecokit::stop_ctx(
                  paste0(
                    "Scraping railways links failed for: ", scraped_links,
                    "after ", attempts),
                  include_backtrace = TRUE)
              }
            }
            attempt <- attempt + 1
          }

          return(scraped_links)
        }
      )) %>%
    tidyr::unnest(cols = "url2") %>%
    dplyr::mutate(
      # download path
      path = purrr::map2(
        .x = url2, .y = country,
        .f = ~ {
          stringr::str_remove_all(.x, "^.+/|-latest-free.shp") %>%
            paste0(.y, "_", .) %>%
            fs::path(path_railway_raw, .)
        }),
      modified_date = purrr::map(
        .x = url2,
        .f = ~ {
          # Last modified date
          paste0("curl -sIL ", .x) %>%
            system(intern = TRUE) %>%
            stringr::str_subset("Last-Modified") %>%
            stringr::str_remove_all("Last-Modified: ") %>%
            readr::parse_datetime(format = "%a, %d %b %Y %H:%M:%S %Z") %>%
            lubridate::as_date()
        }
      )) %>%
    tidyr::unnest(c = "modified_date")

  ecokit::cat_time(
    paste0("There are ", nrow(railway_links), " files to be downloaded"),
    level = 1L)

  save(railway_links, file = fs::path(path_railway, "railway_links.RData"))

  ecokit::cat_diff(
    init_time = .start_time_down, msg_n_lines = 1, level = 1L,
    prefix = "Preparing railways download links took ")

  # # ..................................................................... ###

  # Processing railway data ------
  ecokit::cat_time("Processing railway data")
  .start_time_process <- lubridate::now(tzone = "CET")

  ## Prepare working in parallel ----
  if (n_cores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = n_cores, level = 1L, future_max_size = 800L,
      strategy = strategy)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  ## Processing railway data ----
  ecokit::cat_time("Processing railway data", level = 1L)

  railway_3035 <- future.apply::future_lapply(
    X = seq_len(nrow(railway_links)),
    FUN = function(id) {

      url <- railway_links$url2[[id]]
      zip_path <- railway_links$path[[id]]
      country <- railway_links$country[[id]]
      path_prefix <- stringr::str_remove_all(basename(zip_path), ".zip$")
      path_temp <- fs::path(path_railway_interim, paste0(path_prefix, ".RData"))

      withr::local_options(
        future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
        future.seed = TRUE, timeout = 1800)

      # Check if zip file is a valid file
      if (file.exists(zip_path)) {
        success <- ecokit::check_zip(zip_path)
        if (isFALSE(success)) {
          fs::file_delete(zip_path)
        }
      } else {
        success <- FALSE
      }

      # Try downloading data for a max of 5 attempts
      attempt <- 1
      attempts <- 5

      while (isFALSE(success) && attempt <= attempts) {
        down <- try(
          expr = {

            stringr::str_glue(
              'curl -k -L --connect-timeout 180 --max-time 3600 --retry 5 \\
            "{url}" -o "{zip_path}" --silent') %>%
              system()

            success <- ecokit::check_zip(zip_path)

            if (isFALSE(success)) {
              fs::file_delete(zip_path)
            }

            success
          }, silent = TRUE)

        if (inherits(down, "try-error")) {
          success <- FALSE
        }
        attempt <- attempt + 1
      }

      # Filter only railways files
      in_file_n <- dplyr::tibble(utils::unzip(zip_path, list = TRUE)) %>%
        dplyr::filter(stringr::str_detect(Name, "railways")) %>%
        dplyr::pull(Name) %>%
        unique()

      utils::unzip(
        zipfile = zip_path, files = in_file_n,
        exdir = fs::path(path_railway_interim, path_prefix))

      path_extract <- dplyr::tibble(
        old_name = fs::path(path_railway_interim, path_prefix, in_file_n),
        new_name = fs::path(
          path_railway_interim,
          paste0(path_prefix, ".", tools::file_ext(in_file_n)))) %>%
        dplyr::mutate(Ren = purrr::map2(old_name, new_name, file.rename))

      railway_sf <- dplyr::pull(path_extract, new_name) %>%
        stringr::str_subset(".shp$") %>%
        sf::st_read(quiet = TRUE) %>%
        tibble::tibble() %>%
        sf::st_as_sf() %>%
        sf::st_transform(crs = 3035) %>%
        dplyr::select(
          -tidyselect::all_of(c("bridge", "tunnel", "layer", "name"))) %>%
        sf::st_filter(y = grid_sf, .predicate = sf::st_intersects)

      save(railway_sf, file = path_temp)

      # Clean up
      fs::dir_delete(fs::path(path_railway_interim, path_prefix))
      fs::file_delete(path_extract$new_name)
      if (delete_processed) {
        fs::file_delete(zip_path)
      }

      tibble::tibble(
        url = url, country = country, area = path_prefix, path = path_temp)
    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = pkg_to_export,
    future.globals = c("railway_links", "grid_sf", "path_railway_interim")) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(data = purrr::map(path, ecokit::load_as)) %>%
    tidyr::unnest(data) %>%
    sf::st_as_sf()


  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 1L)
    future::plan("sequential", gc = TRUE)
  }

  # Delete raw directory if empty
  raw_file_list <- list.files(path_railway_raw, full.names = TRUE)
  if (length(raw_file_list) == 0) {
    fs::dir_delete(path_railway_raw)
  }

  # # .................................... ###

  ## Saving - RData -----
  ecokit::cat_time("Saving - RData", level = 1L)
  save(railway_3035, file = fs::path(path_railway, "railway_3035.RData"))

  # # .................................... ###

  ## railway_3035_to_plot ----
  railway_3035_to_plot <- dplyr::filter(railway_3035, fclass == "rail") %>%
    sf::st_join(grid_sf) %>%
    dplyr::filter(!is.na(CellCode)) %>%
    dplyr::select("geometry")

  save(
    railway_3035_to_plot,
    file = fs::path(path_railway, "railway_3035_to_plot.RData"))

  # # .................................... ###

  ## Saving each railway class to separate file ----
  ecokit::cat_time("Saving each railway class to separate file", level = 1L)

  sf::st_drop_geometry(railway_3035) %>%
    dplyr::distinct(fclass) %>%
    dplyr::pull(fclass) %>%
    purrr::walk(
      .f = ~ {
        ecokit::cat_time(.x, level = 2L)
        dplyr::filter(railway_3035, fclass == .x) %>%
          ecokit::save_as(
            object_name = paste0("railway_sf_", .x),
            out_path = fs::path(
              path_railway, paste0("railway_sf_", .x, ".RData")))
      }
    )

  rm(railway_3035, grid_sf, envir = environment())
  invisible(gc())

  ecokit::cat_diff(
    init_time = .start_time_process,
    prefix = "Processing railway data took ", msg_n_lines = 1, level = 1L)

  # # ..................................................................... ###

  # Calculate length of railways for each railway class per grid cell -----
  ecokit::cat_time(
    "Calculate length of railways for each railway class per grid cell")

  railway_length <- path_railway %>%
    list.files(pattern = "^railway_sf_", full.names = TRUE) %>%
    purrr::map(
      .f = ~ {
        map_name <- stringr::str_remove_all(basename(.x), "railway_sf_|.RData")

        ecokit::cat_time(map_name, level = 2L)

        ecokit::load_as(.x) %>%
          terra::vect() %>%
          terra::rasterizeGeom(y = grid_ref, fun = "length", unit = "km") %>%
          terra::mask(mask = grid_ref) %>%
          stats::setNames(map_name) %>%
          terra::toMemory()
      }, .progress = FALSE
    ) %>%
    terra::rast()

  # Sum of railways length of any type
  railway_length$Sum <- sum(railway_length)

  ## Saving - RData -----
  ecokit::cat_time("Saving - RData", level = 1L)
  ecokit::save_as(
    object = terra::wrap(railway_length), object_name = "railway_length",
    out_path = fs::path(path_railway, "railway_length.RData"))

  ## Saving - tif ------
  ecokit::cat_time("Saving - tif", level = 1L)
  terra::writeRaster(
    x = railway_length, overwrite = TRUE,
    filename = fs::path(
      path_railway,
      paste0("railway_length_", names(railway_length), ".tif")))

  # # ..................................................................... ###

  # Calculate distance to rail ------
  ecokit::cat_time("Calculate distance to rail")

  # This calculates the distance from each grid cell to the nearest grid cell
  # overlapping with a railway. This can be different than calculating the
  # actual distance to nearest railway line; which is expected to take too much
  # time to calculate

  railway_distance <- purrr::map(
    .x = as.list(railway_length),
    .f = ~ {
      ecokit::cat_time(names(.x), level = 2L)

      # suppress progress bar
      terra::terraOptions(progress = 0)

      railway_points <- terra::as.points(terra::classify(.x, cbind(0, NA)))

      terra::distance(x = .x, y = railway_points, unit = "km") %>%
        terra::mask(grid_ref) %>%
        stats::setNames(paste0("railway_distance_", names(.x))) %>%
        # Ensure that values are read from memory
        terra::toMemory()
    }
  ) %>%
    terra::rast()

  ecokit::cat_time("Save distance to railways - tif", level = 1L)
  terra::writeRaster(
    x = railway_distance, overwrite = TRUE,
    filename = fs::path(
      path_railway, paste0(names(railway_distance), ".tif")))

  ecokit::cat_time("Save distance to railways - RData", level = 1L)
  ecokit::save_as(
    object = terra::wrap(railway_distance), object_name = "railway_distance",
    out_path = fs::path(path_railway, "railway_distance.RData"))

  rm(grid_ref, envir = environment())

  # # ..................................................................... ###

  # Plotting ------
  ecokit::cat_time("Plotting")

  eu_boundaries <- ecokit::load_as(eu_boundaries) %>%
    magrittr::extract2("Bound_sf_Eur_s") %>%
    magrittr::extract2("L_03")

  plot_theme <- ggplot2::theme_bw() +
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
  ecokit::cat_time("Plotting length of railways", level = 1L)

  rail_to_plot <- terra::subset(railway_length, "rail") %>%
    terra::classify(cbind(0, NA))

  rail_plot <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      eu_boundaries, mapping = ggplot2::aes(), color = "grey75",
      linewidth = 0.075, fill = "grey98") +
    tidyterra::geom_spatraster(data = log10(rail_to_plot + 1), maxcell = Inf) +
    ggplot2::geom_sf(
      eu_boundaries, mapping = ggplot2::aes(), color = "grey30",
      linewidth = 0.075, fill = "transparent", inherit.aes = TRUE) +
    paletteer::scale_fill_paletteer_c(
      na.value = "transparent", "viridis::plasma") +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)),
      limits = c(2600000, 6700000)) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)),
      limits = c(1450000, 5420000)) +
    ggplot2::labs(title = "Railway length", fill = "log10") +
    plot_theme

  # Using ggplot2::ggsave directly does not show non-ascii characters
  # correctly
  ragg::agg_jpeg(
    filename = fs::path(path_railway, "railway_length.jpeg"),
    width = 31, height = 30, res = 600, quality = 100, units = "cm")
  print(rail_plot)
  grDevices::dev.off()

  rm(rail_plot, envir = environment())

  # # .................................. ###

  ## Plotting railways -----
  ecokit::cat_time("Plotting European railways", level = 1L)

  rail_plot_shp <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      eu_boundaries, mapping = ggplot2::aes(), color = "grey75",
      linewidth = 0.075, fill = "grey98") +
    ggplot2::geom_sf(
      railway_3035_to_plot,
      mapping = ggplot2::aes(), color = "blue", linewidth = 0.05) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)),
      limits = c(2600000, 6700000)) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0)),
      limits = c(1450000, 5420000)) +
    ggplot2::labs(title = "Railways in Europe", fill = NA) +
    plot_theme

  # Using ggplot2::ggsave directly does not show non-ascii characters
  # correctly
  ragg::agg_jpeg(
    filename = fs::path(path_railway, "railway_lines.jpeg"),
    width = 31, height = 30, res = 600, quality = 100, units = "cm")
  print(rail_plot_shp)
  grDevices::dev.off()

  rm(rail_plot_shp, envir = environment())

  # # ..................................................................... ###

  # Function Summary ----

  ecokit::cat_diff(
    init_time = .start_time,
    prefix = "\nProcessing railway data was finished in ", ... = "\n")

  return(invisible(NULL))
}
