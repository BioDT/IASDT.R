## |------------------------------------------------------------------------| #
# rc_plot_species_all ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname response_curves
#' @name response_curves
#' @order 3
#' @author Ahmed El-Gabbas

rc_plot_species_all <- function(
    model_dir = NULL, n_cores = 8L, return_data = FALSE, plotting_alpha = 0.3) {

  ecokit::cat_time("Plotting species response curves")

  # # ..................................................................... ###
  # Check arguments

  ecokit::cat_time("Check arguments", level = 1L)

  ecokit::check_args(args_to_check = "model_dir", args_type = "character")
  ecokit::check_args(args_to_check = "plotting_alpha", args_type = "numeric")
  ecokit::check_args(args_to_check = "return_data", args_type = "logical")

  if (plotting_alpha < 0 || plotting_alpha > 1) {
    ecokit::stop_ctx(
      "`plotting_alpha` must be between 0 and 1",
      plotting_alpha = plotting_alpha, include_backtrace = TRUE)
  }

  n_cores <- .validate_n_cores(n_cores)
  .start_time <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  coords <- path_rc_prob <- nfv <- data <- plot_data <- variable <-
    variable_2 <- i <- var_desc <- var_desc_2 <- NULL

  # # ..................................................................... ###

  path_rc_data <- fs::path(
    model_dir, "model_postprocessing", "response_curves_data")
  path_rc_all <- fs::path(
    model_dir, "model_postprocessing", "response_curves_all")

  if (!dir.exists(path_rc_data)) {
    ecokit::stop_ctx(
      "Response curve data subfolder is missing.", path_rc_data = path_rc_data,
      include_backtrace = TRUE)
  }

  fs::dir_create(path_rc_all)

  # # ..................................................................... ###

  # Loading & processing species response curve data in parallel

  ecokit::cat_time(
    "Loading & processing species response curve data in parallel", level = 1L)

  sp_data_all <- fs::path(path_rc_data, "rc_data.RData") %>%
    ecokit::load_as() %>%
    dplyr::select(coords, path_rc_prob)

  doParallel::registerDoParallel(cores = n_cores)
  ecokit::load_packages(package_list = "foreach")
  withr::defer(doParallel::stopImplicitCluster())

  rc_data <- foreach::foreach(
    i = sp_data_all$path_rc_prob, .packages = c("ecokit", "magrittr", "dplyr")
  ) %dopar% { # nolint: object_usage_linter
    ecokit::load_as(i) %>%
      dplyr::select(
        tidyselect::all_of(c("variable", "nfv", "species", "plot_data_quant")))
  }

  sp_data_all <- sp_data_all %>%
    dplyr::mutate(data = rc_data) %>%
    tidyr::unnest(data) %>%
    dplyr::select(-path_rc_prob) %>%
    dplyr::slice(gtools::mixedorder(variable)) %>%
    tidyr::nest(data = -c(nfv, coords))

  doParallel::stopImplicitCluster()

  invisible(gc())

  # # ..................................................................... ###

  # plot all species response curves

  ecokit::cat_time("Plot all species response curves", level = 1L)

  plots <- purrr::map_dfr(
    .x = seq_len(nrow(sp_data_all)),
    .f = function(id) {

      nfv <- sp_data_all$nfv[[id]]
      coords <- sp_data_all$coords[[id]]

      file_prefix <- paste0("rc_all_nfv_", nfv, "_coords_", coords)
      path_jpeg <- fs::path(path_rc_all, paste0(file_prefix, ".jpeg"))

      plot_data <- sp_data_all$data[[id]] %>%
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
            "<span style='font-size: 8pt;'>", var_desc_2, "</span>"),

          variable_2 = factor(variable, levels = unique(variable))) %>%
        dplyr::group_split(variable_2)

      plot_subtitle <- dplyr::if_else(
        nfv == 1,
        paste0(
          "Non-focal variables are set to most likely value <i>",
          "[non.focalVariables = 1]</i>"),
        paste0(
          "Non-focal variables = most likely value given ",
          "value of focal variable <i>[non.focalVariables = 2]</i>"))
      plot_caption <- dplyr::if_else(
        coords == "c", "Mean coordinates", "No spatial dependence")
      plot_caption <- paste0(plot_caption, " --- ", plot_subtitle)

      if (length(plot_data) <= 9) {
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
        .x = plot_data,
        .f = ~ {
          plot <- dplyr::select(.x, species, plot_data_quant) %>%
            tidyr::unnest("plot_data_quant") %>%
            dplyr::filter(quantile == 0.5)
          plot <- plot %>%
            ggplot2::ggplot(
              mapping = ggplot2::aes(x = x_vals, y = pred, group = species),
              environment = emptyenv()) +
            ggplot2::geom_line(
              linetype = 1, linewidth = 0.3,
              colour = "blue", alpha = plotting_alpha) +
            ggplot2::scale_y_continuous(
              limits = c(-0.01, 1.01), oob = scales::squish,
              expand = c(0, 0)) +
            ggplot2::scale_x_continuous(expand = c(0.015, 0.015)) +
            ggplot2::labs(
              x = NULL, y = NULL,
              title = .x$var_desc[[1]], subtitle = .x$var_desc_2[[1]]) +
            ggplot2::theme_bw() +
            ggplot2::theme(
              legend.position = "none",
              axis.title = ggtext::element_markdown(
                size = 12, face = "bold"),
              axis.text = ggplot2::element_text(size = 8),
              plot.title = ggtext::element_markdown(
                size = 24, hjust = 0,
                margin = ggplot2::margin(-5, 0, -5, 0)),
              plot.subtitle = ggtext::element_markdown(
                margin = ggplot2::margin(0, 0, 0, 0), hjust = 0),
              plot.caption = ggtext::element_markdown(
                size = 12, color = "grey", hjust = 0),
              panel.grid.major = ggplot2::element_line(linewidth = 0.25),
              panel.grid.minor = ggplot2::element_line(linewidth = 0.1),
              plot.margin = ggplot2::unit(c(0.1, 0.2, 0.1, 0.2), "lines"))
          plot

        }) %>%
        patchwork::wrap_plots(nrow = n_rows, ncol = n_columns) +
        patchwork::plot_annotation(
          title = "Response curves for all species",
          caption = plot_caption,
          theme = ggplot2::theme(
            plot.margin = ggplot2::margin(0.1, 0, 0, 0, "cm"),
            plot.title = ggtext::element_markdown(face = "bold"),
            plot.caption = ggtext::element_markdown(
              size = 12, color = "grey", hjust = 0))) +
        patchwork::plot_layout(axes = "collect")

      plots <- patchwork::wrap_elements(plots) +
        ggplot2::labs(tag = "<b>Predicted habitat suitability</b>") +
        ggplot2::theme(
          plot.tag = ggtext::element_markdown(
            size = ggplot2::rel(1.25), angle = 90),
          plot.tag.position = "left")


      # Using ggplot2::ggsave directly does not show non-ascii characters
      # correctly
      ragg::agg_jpeg(
        filename = path_jpeg, width = plot_width, height = plot_height,
        res = 600, quality = 100, units = "cm")
      print(plots)
      grDevices::dev.off()

      tibble::tibble(
        path_jpeg = path_jpeg, plot_height = plot_height,
        plot_width = plot_width)
    })

  # # ..................................................................... ###

  # Save data
  ecokit::cat_time("Save data", level = 1L)

  sp_data_all <- dplyr::select(sp_data_all, -plot_data) %>%
    dplyr::bind_cols(plots = plots)

  save(sp_data_all, file = fs::path(path_rc_all, "sp_data_all.RData"))

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time,
    prefix = "Plotting all species response curves took ", level = 1L)

  # # ..................................................................... ###

  if (return_data) {
    return(sp_data_all)
  } else {
    return(invisible(NULL))
  }
}
