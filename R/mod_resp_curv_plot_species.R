
## |------------------------------------------------------------------------| #
# rc_plot_species ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname response_curves
#' @name response_curves
#' @order 2
#' @author Ahmed El-Gabbas

rc_plot_species <- function(
    model_dir = NULL, n_cores = 20, env_file = ".env", return_data = FALSE) {

  .start_time <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Check arguments

  ecokit::cat_time("Check arguments")
  ecokit::check_args(args_to_check = "model_dir", args_type = "character")
  ecokit::check_args(args_to_check = "return_data", args_type = "logical")

  n_cores <- .validate_n_cores(n_cores)

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  path_pa <- n_cells_naturalized <- nfv <- coords <- species <- prefix <-
    data <- path_rc_prob <- variable <- ias_id <- var_desc <- var_desc_2 <-
    rc_file <- NULL

  # # ..................................................................... ###

  ecokit::cat_time("Plotting species response curves")

  path_rc_data <- fs::path(
    model_dir, "model_postprocessing", "response_curves_data")
  if (!dir.exists(path_rc_data)) {
    ecokit::stop_ctx(
      "Response curve data subfolder is missing.", path_rc_data = path_rc_data,
      include_backtrace = TRUE)
  }
  path_rc_sp <- fs::path(model_dir, "model_postprocessing", "rc_sp")
  path_rc_sp_data <- fs::path(model_dir, "model_postprocessing", "rc_sp_data")
  fs::dir_create(c(path_rc_sp, path_rc_sp_data))

  # # ..................................................................... ###

  # # Load species summary
  ecokit::cat_time("Load species summary")

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_pa", "DP_R_pa", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  sp_summary <- fs::path(path_pa, "sp_pa_summary_df.RData")
  if (!file.exists(sp_summary)) {
    ecokit::stop_ctx(
      "sp_summary file does not exist", sp_summary = sp_summary,
      include_backtrace = TRUE)
  }

  sp_summary <- ecokit::load_as(sp_summary) %>%
    dplyr::select(tidyselect::all_of(c("ias_id", "n_cells_naturalized"))) %>%
    dplyr::rename(n_cells = n_cells_naturalized)

  # # ..................................................................... ###

  # Load species names
  ecokit::cat_time("Load species names")
  sp_names <- IASDT.R::get_species_name(env_file = env_file)

  # # ..................................................................... ###

  # Prepare species-specific data in parallel

  ecokit::cat_time("Prepare species-specific data in parallel")

  sp_data_all <- fs::path(path_rc_data, "rc_data.RData") %>%
    ecokit::load_as() %>%
    dplyr::select(tidyselect::all_of(c("coords", "path_rc_prob"))) %>%
    dplyr::mutate(
      data = purrr::map(.x = path_rc_prob, .f = ecokit::load_as)) %>%
    tidyr::unnest(data) %>%
    dplyr::select(-path_rc_prob) %>%
    tidyr::nest(
      data0 = tidyselect::everything(), .by = c(nfv, coords, species)) %>%
    dplyr::mutate(ias_id = as.numeric(stringr::str_remove(species, "^sp_"))) %>%
    dplyr::left_join(sp_summary, by = "ias_id") %>%
    dplyr::mutate(
      prefix = paste0(species, "_nfv_", nfv, "_coords_", coords),
      path_jpeg_fixed = fs::path(path_rc_sp, paste0(prefix, "_fixed.jpeg")),
      path_jpeg_free = fs::path(path_rc_sp, paste0(prefix, "_free.jpeg")),
      path_sp_data = fs::path(path_rc_sp_data, paste0(prefix, ".qs2")))


  ecokit::cat_time("Export species-specific data", level = 1L)
  purrr::walk(
    .x = seq_len(nrow(sp_data_all)),
    .f = function(id) {
      data0 <- dplyr::slice(sp_data_all, id)
      if (isFALSE(ecokit::check_data(data0$path_sp_data, warning = FALSE))) {
        ecokit::save_as(object = data0, out_path = data0$path_sp_data)
      }
    })

  sp_data_all <- gtools::mixedsort(sp_data_all$path_sp_data)
  invisible(gc())

  # # ..................................................................... ###

  ecokit::cat_time("Plotting species-specific data in parallel")

  # set up parallel processing
  doParallel::registerDoParallel(cores = n_cores)
  ecokit::load_packages(package_list = "foreach")
  withr::defer(doParallel::stopImplicitCluster())

  pkg_to_export <- c(
    "dplyr", "purrr", "tidyr", "gtools", "ggtext", "patchwork", "magrittr",
    "ggplot2", "tibble", "ecokit", "ragg", "stringr", "scales", "stats",
    "fs", "grDevices", "grid", "rlang", "qs2")

  plots <- foreach::foreach(
    rc_file = sp_data_all, .export = c("sp_names", "sp_data_all"),
    .packages = pkg_to_export) %dopar% {

      data0 <- ecokit::load_as(rc_file)

      coords <- data0$coords
      species <- data0$species
      n_cells <- data0$n_cells
      nfv <- data0$nfv
      path_jpeg_fixed <- data0$path_jpeg_fixed
      path_jpeg_free <- data0$path_jpeg_free

      data0 <- dplyr::select(data0, data0) %>%
        tidyr::unnest(data0) %>%
        dplyr::slice(gtools::mixedorder(variable)) %>%
        dplyr::mutate(
          var_desc = dplyr::case_when(
            startsWith(variable, "bio") ~ stringr::str_to_sentence(variable),
            variable == "road_rail_log" ~ "Road + Rail intensity",
            variable == "efforts_log" ~ "Sampling efforts",
            variable == "rivers_log" ~ "River length",
            variable == "habitat_log" ~ "% habitat coverage",
            variable == "soil" ~ "Soil density",
            variable == "wetness" ~ "Wetness index",
            .default = variable),
          var_desc = paste0(
            "<span style='font-size: 10pt;'><b>", var_desc, "</b></span>"),

          var_desc_2 = dplyr::case_when(
            variable == "bio1" ~ "annual mean temperature",
            variable == "bio2" ~ "mean diurnal range",
            variable == "bio3" ~ "isothermality (bio2/bio7) (&times;100)",
            variable == "bio4" ~ "temperature seasonality",
            variable == "bio5" ~ "max temperature of warmest month",
            variable == "bio6" ~ "temperature of the coldest month",
            variable == "bio7" ~ "temperature annual range (bio5-bio6)",
            variable == "bio8" ~ "temperatures of the wettest quarter",
            variable == "bio9" ~ "mean temperature of driest quarter",
            variable == "bio10" ~ "mean temperature of warmest quarter",
            variable == "bio11" ~ "mean temperature of coldest quarter",
            variable == "bio12" ~ "annual precipitation amount",
            variable == "bio13" ~ "precipitation of wettest month",
            variable == "bio14" ~ "precipitation of driest month",
            variable == "bio15" ~ "precipitation seasonality",
            variable == "bio16" ~ "precipitation of wettest quarter",
            variable == "bio17" ~ "precipitation of driest quarter",
            variable == "bio18" ~ "precipitation of the warmest quarter",
            variable == "bio19" ~ "precipitation of coldest quarter",
            variable == "npp" ~ "net primary productivity",
            variable == "rivers_log" ~ " (log<sub>10</sub>(x + 0.1))",
            variable == "road_rail_log" ~ " (log<sub>10</sub>(x + 1))",
            variable == "efforts_log" ~ " (log<sub>10</sub>(x + 1))",
            variable == "soil" ~ "Soil bulk density",
            variable == "wetness" ~ "Topographic wetness index",
            variable == "habitat_log" ~ " (log<sub>10</sub>(x + 0.1))",
            .default = variable),
          var_desc_2 = paste0(
            "<span style='font-size: 8pt;'>", var_desc_2, "</span>"))

      # nolint start
      species2 <- dplyr::filter(sp_names, ias_id == !!species)
      species_name <- species2$species_name
      species_id <- stringr::str_remove(species, "^sp_")
      class <- species2$class
      order <- species2$order
      family <- species2$family
      title_text <- stringr::str_glue(
        '<span style="font-size:13pt;"><b> Response curves of </b></span>\\
        <span style="color:blue; font-size:13pt;">\\
        <b><i>{species_name}</i></b></span>\\
        <span style="font-size:8pt;"> (\\
        <b>Class:</b> {class} &#8212; <b>Order:</b> {order} &#8212; \\
        <b>Family:</b> {family} &#8212; <b>ID:</b> {species_id} \\
        &#8212; <b># presence grids:</b> {n_cells})</span>')
      # nolint end

      subtitle_text <- dplyr::if_else(
        nfv == 1,
        paste0(
          "Non-focal variables are set to most likely value <i>",
          "[non.focalVariables = 1]</i>"),
        paste0(
          "Non-focal variables = most likely value given ",
          "value of focal variable <i>[non.focalVariables = 2]</i>"))
      plot_caption <- dplyr::if_else(
        coords == "c", "Mean coordinates", "No spatial dependence")
      plot_caption <- paste0(plot_caption, " --- ", subtitle_text)

      if (nrow(data0) <= 9) {
        n_rows <- n_columns <- 3
        plot_width <- 24
        plot_height <- 22
      } else {
        n_rows <- 3
        n_columns <- 4
        plot_width <- 30
        plot_height <- 22
      }

      plots <- purrr::map(
        .x = seq_len(nrow(data0)),
        .f = ~ {

          positive_trend_prob <- data0$positive_trend_prob[[.x]]
          trend <- tibble::tibble(
            trend_2 = stringr::str_glue(
              "\n     Pr[pred(Var=max)] > Pr[pred(Var=min)] = \\
              {round(positive_trend_prob, 2)}"),
            X = -Inf, Y = Inf)

          quant_data <- data0$plot_data_quant[[.x]] %>%
            tidyr::pivot_wider(
              id_cols = XVals, names_from = quantile, values_from = pred) %>%
            stats::setNames(c("XVals", "q_25", "q_50", "q_975"))

          rug_0 <- dplyr::filter(data0$observed_pa[[.x]], pred == 0)
          rug_1 <- dplyr::filter(data0$observed_pa[[.x]], pred == 1)

          plot_theme <- ggplot2::theme_bw() +
            ggplot2::theme(
              legend.position = "none",
              axis.title = ggtext::element_markdown(size = 12, face = "bold"),
              axis.text = ggplot2::element_text(size = 8),
              plot.title = ggtext::element_markdown(
                size = 24, hjust = 0, margin = ggplot2::margin(-5, 0, -5, 0)),
              plot.subtitle = ggtext::element_markdown(
                margin = ggplot2::margin(0, 0, 0, 0), hjust = 0),
              plot.caption = ggtext::element_markdown(
                size = 10, color = "grey", hjust = 0),
              panel.grid.major = ggplot2::element_line(linewidth = 0.25),
              panel.grid.minor = ggplot2::element_line(linewidth = 0.1),
              plot.margin = ggplot2::unit(c(0.1, 0.2, 0.1, 0.2), "lines"))

          # Fixed y-axis
          plot_fixed <- ggplot2::ggplot(
            data = quant_data, mapping = ggplot2::aes(x = XVals),
            environment = emptyenv()) +
            ggplot2::geom_line(
              ggplot2::aes(y = q_975), data = quant_data,
              linetype = 2, linewidth = 0.3, colour = "blue") +
            ggplot2::geom_line(
              ggplot2::aes(y = q_25), data = quant_data,
              linetype = 2, linewidth = 0.3, colour = "blue") +
            ggplot2::geom_ribbon(
              mapping = ggplot2::aes(ymin = q_25, ymax = q_975),
              data = quant_data, fill = "blue", alpha = 0.1) +
            ggplot2::geom_line(
              mapping = ggplot2::aes(y = q_50), data = quant_data,
              linetype = 1, linewidth = 0.6, colour = "blue") +
            ggplot2::geom_rug(
              sides = "t", data = rug_1, ggplot2::aes(x = XVals),
              color = "blue", linewidth = 0.025, alpha = 0.25,
              length = grid::unit(0.03, "npc")) +
            ggplot2::geom_rug(
              sides = "b", data = rug_0, ggplot2::aes(x = XVals),
              color = "red", linewidth = 0.025, alpha = 0.25,
              length = grid::unit(0.03, "npc")) +
            ggplot2::geom_text(
              data = trend, vjust = 0.5, hjust = -0.05,
              mapping = ggplot2::aes(x = X, y = Y, label = trend_2),
              colour = "grey30", size = 2.75) +
            ggplot2::scale_y_continuous(
              limits = c(0, 1), oob = scales::squish_infinite) +
            ggplot2::scale_x_continuous(expand = c(0.015, 0.015)) +
            ggplot2::labs(
              x = NULL, y = NULL, title = data0$var_desc[[.x]],
              subtitle = data0$var_desc_2[[.x]]) +
            plot_theme

          # Free y-axis
          plot_free <- ggplot2::ggplot(
            data = quant_data, mapping = ggplot2::aes(x = XVals),
            environment = emptyenv()) +
            ggplot2::geom_line(
              ggplot2::aes(y = q_975), data = quant_data,
              linetype = 2, linewidth = 0.3, colour = "blue") +
            ggplot2::geom_line(
              ggplot2::aes(y = q_25), data = quant_data,
              linetype = 2, linewidth = 0.3, colour = "blue") +
            ggplot2::geom_ribbon(
              mapping = ggplot2::aes(ymin = q_25, ymax = q_975),
              data = quant_data, fill = "blue", alpha = 0.1) +
            ggplot2::geom_line(
              mapping = ggplot2::aes(y = q_50), data = quant_data,
              linetype = 1, linewidth = 0.6, colour = "blue") +
            ggplot2::geom_rug(
              sides = "t", data = rug_1, ggplot2::aes(x = XVals),
              color = "blue", linewidth = 0.025, alpha = 0.25,
              length = grid::unit(0.03, "npc")) +
            ggplot2::geom_rug(
              sides = "b", data = rug_0, ggplot2::aes(x = XVals),
              color = "red", linewidth = 0.025, alpha = 0.25,
              length = grid::unit(0.03, "npc")) +
            ggplot2::geom_text(
              data = trend, vjust = 0.5, hjust = -0.05,
              mapping = ggplot2::aes(x = X, y = Y, label = trend_2),
              colour = "grey30", size = 2.75) +
            ggplot2::scale_y_continuous(oob = scales::squish_infinite) +
            ggplot2::scale_x_continuous(expand = c(0.015, 0.015)) +
            ggplot2::labs(
              x = NULL, y = NULL, title = data0$var_desc[[.x]],
              subtitle = data0$var_desc_2[[.x]]) +
            plot_theme +
            ggplot2::theme(
              axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5))

          return(
            tibble::tibble(
              plot_fixed = list(plot_fixed), plot_free = list(plot_free)))
        }) %>%
        dplyr::bind_rows()

      plot_theme_2 <- ggplot2::theme(
        plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
        plot.title = ggtext::element_markdown(
          hjust = 0, margin = ggplot2::margin(0, 0, -0.125, 0, "cm")),
        plot.caption = ggtext::element_markdown(
          size = 12, color = "grey", hjust = 0))

      # Fixed y-axis
      plot_fixed <- patchwork::wrap_plots(
        plots$plot_fixed, nrow = n_rows, ncol = n_columns) +
        patchwork::plot_annotation(
          title = stringr::str_glue(
            "{title_text}<span style = 'color:#ffffff;'>...........\\
            </span><span style='font-size:10pt; color:grey;'>Fixed \\
            y-axis</span>"),
          caption = plot_caption, theme = plot_theme_2) +
        patchwork::plot_layout(axes = "collect")

      plot_fixed <- patchwork::wrap_elements(panel = plot_fixed) +
        ggplot2::labs(tag = "<b>Predicted habitat suitability</b>") +
        ggplot2::theme(
          plot.tag = ggtext::element_markdown(
            size = ggplot2::rel(1), angle = 90),
          plot.tag.position = "left")

      ragg::agg_jpeg(
        filename = path_jpeg_fixed, width = plot_width, height = plot_height,
        res = 600, quality = 100, units = "cm")
      print(plot_fixed)
      grDevices::dev.off()


      # Free y-axis
      plot_free <- patchwork::wrap_plots(
        plots$plot_free, nrow = n_rows, ncol = n_columns) +
        patchwork::plot_annotation(
          title = stringr::str_glue(
            "{title_text}<span style = 'color:#ffffff;'>...........</span>\\",
            "<span style='font-size:10pt; color:grey;'>Free y-axis</span>"),
          caption = plot_caption, theme = plot_theme_2) +
        patchwork::plot_layout(axes = "collect")

      plot_free <- patchwork::wrap_elements(panel = plot_free) +
        ggplot2::labs(tag = "<b>Predicted habitat suitability</b>") +
        ggplot2::theme(
          plot.tag = ggtext::element_markdown(
            size = ggplot2::rel(1), angle = 90),
          plot.tag.position = "left")

      ragg::agg_jpeg(
        filename = path_jpeg_free, width = plot_width, height = plot_height,
        res = 600, quality = 100, units = "cm")
      print(plot_free)
      grDevices::dev.off()

      out_data <- tibble::tibble(
        coords = coords, species = species, n_cells = n_cells, nfv = nfv,
        path_jpeg_fixed = path_jpeg_fixed, path_jpeg_free = path_jpeg_free,
        plot_height = plot_height, plot_width = plot_width)

      return(out_data)
    }

  plots <- dplyr::bind_rows(plots)

  # stopping the cluster
  doParallel::stopImplicitCluster()
  invisible(gc())

  # # ..................................................................... ###

  # Save data
  ecokit::cat_time("Save data")
  save(plots, file = fs::path(path_rc_sp, "sp_data_all.RData"))

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time, prefix = "Plotting species response curves took ")

  # # ..................................................................... ###

  if (return_data) {
    return(sp_data_all)
  } else {
    return(invisible(NULL))
  }
}
