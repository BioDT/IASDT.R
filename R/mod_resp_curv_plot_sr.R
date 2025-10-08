## |------------------------------------------------------------------------| #
# rc_plot_sr ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname response_curves
#' @name response_curves
#' @order 4
#' @author Ahmed El-Gabbas

rc_plot_sr <- function(
    model_dir, verbose = TRUE, n_cores = 8L, future_max_size = 1000L,
    strategy = "multisession") {

  .start_time <- lubridate::now(tzone = "CET")

  ecokit::check_args(args_to_check = "model_dir", args_type = "character")
  ecokit::check_args(args_to_check = "verbose", args_type = "logical")

  strategy <- .validate_strategy(strategy)
  if (strategy == "sequential") n_cores <- 1L
  n_cores <- .validate_n_cores(n_cores)

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  trend_2 <- variable <- quant <- observed <- trend <- nfv <- coords <-
    rc_path_sr <- path_rc_orig <- path_rc_prob <- data <- x_vals <-
    pred <- q_975 <- q_25 <- q_50 <- X <- Y <- variable_2 <- var_2 <-
    variable_1 <- NULL

  # # ..................................................................... ###

  path_rc_data <- fs::path(
    model_dir, "model_postprocessing", "response_curves_data")
  path_rc_sr <- fs::path(model_dir, "model_postprocessing", "rc_sr")

  if (!dir.exists(path_rc_data)) {
    ecokit::stop_ctx(
      "Response curve data subfolder is missing.", path_rc_data = path_rc_data,
      include_backtrace = TRUE)
  }

  fs::dir_create(path_rc_sr)

  # # ..................................................................... ###

  ecokit::cat_time(
    "Create species richness response curves", verbose = verbose)

  sr_data_all <- fs::path(path_rc_data, "rc_data.RData") %>%
    ecokit::load_as() %>%
    dplyr::select(-path_rc_orig, -path_rc_prob)

  n_cores <- max(min(n_cores, nrow(sr_data_all)), 1)

  if (n_cores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = n_cores, level = 1L, future_max_size = future_max_size,
      strategy = strategy)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  ecokit::cat_time("Prepare data", level = 1L, verbose = verbose)
  pkg_to_export <- ecokit::load_packages_future(
    packages = c("tibble", "dplyr", "magrittr", "ecokit", "tidyr"),
    strategy = strategy)

  sr_data_all <- sr_data_all %>%
    dplyr::mutate(
      data = furrr::future_map(
        .x = rc_path_sr,
        .f = ~ {
          data <- ecokit::load_as(.x) %>%
            magrittr::inset("rc_data_sr", NULL)

          quant <- data$rc_data_sr_quant %>%
            dplyr::mutate(
              variable = data$variable, nfv = data$nfv, .before = 1) %>%
            tidyr::pivot_wider(
              id_cols = c("variable", "nfv", "x_vals"),
              names_from = quantile, values_from = sr) %>%
            setNames(c("variable", "nfv", "x_vals", "q_25", "q_50", "q_975"))

          observed <- dplyr::mutate(
            .data = data$observed_sr,
            variable = data$variable, nfv = data$nfv, .before = 1)

          trend <- tibble::tibble(
            variable = data$variable, nfv = data$nfv,
            trend = data$sr_positive_trend_prob)

          return(list(quant = quant, observed = observed, trend = list(trend)))
        },
        .options = furrr::furrr_options(
          seed = TRUE, chunk_size = 1, packages = pkg_to_export))) %>%
    dplyr::select(-nfv, -rc_path_sr) %>%
    tidyr::unnest_wider(data) %>%
    tidyr::nest(.by = c(variable, coords)) %>%
    dplyr::mutate(
      quant = purrr::map(.x = data, .f = ~ dplyr::bind_rows(.x$quant)),
      observed = purrr::map(.x = data, .f = ~ dplyr::bind_rows(.x$observed)),
      trend = purrr::map(.x = data, .f = ~ dplyr::bind_rows(.x$trend))) %>%
    dplyr::select(-data)

  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 1L)
    future::plan("sequential", gc = TRUE)
  }

  invisible(gc())

  # # ..................................................................... ###

  # plot species richness response curves

  ecokit::cat_time(
    "plot species richness response curves", level = 1L, verbose = verbose)

  var_label <- tibble::tribble(
    ~variable_1, ~var_2, ~variable_2,
    "bio1", "Bio1", "Annual mean temperature",
    "bio2", "Bio2", "Mean diurnal range",
    "bio3", "Bio3", "Isothermality (bio2/bio7) (&times;100)",
    "bio4", "Bio4", "Temperature seasonality [standard deviation &times;100]",
    "bio5", "Bio5", "Max temperature of warmest month",
    "bio6", "Bio6", "Temperature of the coldest month",
    "bio7", "Bio7", "Temperature annual range (bio5-bio6)",
    "bio8", "Bio8", "Temperatures of the wettest quarter",
    "bio9", "Bio9", "Mean temperature of driest quarter",
    "bio10", "Bio10", "Mean temperature of warmest quarter",
    "bio11", "Bio11", "Mean temperature of coldest quarter",
    "bio12", "Bio12", "Annual precipitation amount",
    "bio13", "Bio13", "Precipitation of wettest month",
    "bio14", "Bio14", "Precipitation of driest month",
    "bio15", "Bio15", "Precipitation seasonality [Coefficient of Variation]",
    "bio16", "Bio16", "Precipitation of wettest quarter",
    "bio17", "Bio17", "Precipitation of driest quarter",
    "bio18", "Bio18", "Monthly precipitation amount of the warmest quarter",
    "bio19", "Bio19", "Precipitation of coldest quarter",
    "soil", "soil density", "Soil bulk density",
    "wetness", "Wetness index", "Topographic wetness index",
    "npp", "NPP", "net primary productivity",
    "rivers_log", "River length", "log10(x + 0.1)",
    "road_rail_log", "Road + Rail intensity", "log10(x + 1)",
    "habitat_log", "% habitat coverage", "log10(x + 0.1)",
    "efforts_log", "Sampling efforts", "log10(x + 1)") %>%
    dplyr::mutate(
      variable_2 = paste0(
        "<span style='font-size: 10pt;'><b>", var_2, "</b></span>",
        "<span style='font-size: 8pt;'> (", variable_2, ")</span>")) %>%
    dplyr::select(-var_2)

  sr_data_all <- sr_data_all %>%
    dplyr::mutate(
      plot = purrr::pmap(
        .l = list(variable, quant, observed, trend, coords),
        .f = function(variable, quant, observed, trend, coords) {

          ecokit::cat_time(
            paste0(variable, " - coords = ", coords),
            level = 2L, verbose = verbose, cat_timestamp = FALSE)

          # Maximum value on the y-axis
          plot_max <- max(observed$pred, quant$q_975) * 1.05

          # Label on the y axis
          y_axis_label <- paste0(
            "<span style='font-size: 12pt;'><b>Predicted ",
            "species richness</span></b>")

          # trend text
          data_trend <- trend %>%
            dplyr::mutate(
              trend_2 = stringr::str_glue(
                "\n     Pr[pred(Var=max)] > ",
                "Pr[pred(Var=min)] = {round(trend, 2)}"),
              trend_2 = as.character(trend_2), X = -Inf, Y = Inf)

          # variable long name (x-axis label)
          variable_2 <- dplyr::filter(var_label, variable_1 == variable) %>%
            dplyr::pull(variable_2)

          var_name <- dplyr::case_when(
            variable == "habitat_log" ~ "% Habitat coverage",
            variable == "road_rail_log" ~ "Road + Rail intensity",
            variable == "efforts_log" ~ "Sampling efforts",
            variable == "rivers_log" ~ "River length",
            variable == "soil" ~ "Soil density",
            variable == "wetness" ~ "Wetness index",
            .default = variable)

          # faceting labels
          label_1 <- paste0(
            '<span style="font-size:8pt; color:red;"><b>non.focalVariables',
            '= 1</b></span><br><span style="font-size:6pt;">values of ',
            "non-focal variables are set to the most likely value</span>")
          label_2 <- paste0(
            '<span style="font-size:8pt; color:red;"><b>non.focalVariables',
            '= 2</b></span><br><span style="font-size:6pt;">values of ',
            "non-focal variables are set to the most likely value given the",
            " value of focal variable</span>")
          facet_label <- ggplot2::as_labeller(c(`1` = label_1, `2` = label_2))

          # plot title
          title_text <- paste0(
            '<span style="font-size:12pt; color:blue;">',
            "Predicted species richness for <b>", var_name, "</b></span>")

          plot_caption <- dplyr::if_else(
            coords == "c",
            "Predictions are made at mean coordinates",
            paste0(
              "Predictions are made for infinite coordinates ",
              "without effect of spatial dependence"))

          # Plotting theme
          plot_theme <- ggplot2::theme(
            plot.margin = ggplot2::margin(1.5, 6, -2, 3),
            strip.text = ggtext::element_markdown(
              hjust = 0,
              margin = ggplot2::margin(0.05, 0.1, 0.05, 0.1, "cm")),
            strip.background = ggplot2::element_rect(
              colour = NA, fill = "white"),
            legend.position = "none",
            plot.caption = ggplot2::element_text(
              size = 8, color = "grey", hjust = 0, vjust = 4),
            axis.title = ggtext::element_markdown(),
            axis.text = ggplot2::element_text(size = 8),
            plot.title = ggtext::element_markdown(
              margin = ggplot2::margin(1, 0, 1, 0)),
            plot.subtitle = ggtext::element_markdown(
              size = 7, colour = "darkgrey",
              margin = ggplot2::margin(4, 0, 4, 0)),
            panel.grid.major = ggplot2::element_line(linewidth = 0.25),
            panel.grid.minor = ggplot2::element_line(linewidth = 0.1),
            panel.spacing = ggplot2::unit(0.15, "lines"))

          # plot
          plot <- ggplot2::ggplot(
            data = observed, mapping = ggplot2::aes(x = x_vals, y = pred),
            environment = emptyenv()) +
            ggplot2::geom_jitter(
              shape = 16, width = 0, height = 0.02, alpha = 0.2, size = 0.75,
              colour = "darkgrey") +
            ggplot2::geom_line(
              ggplot2::aes(x = x_vals, y = q_975),
              data = quant, linetype = 2, linewidth = 0.3, colour = "blue") +
            ggplot2::geom_line(
              ggplot2::aes(x = x_vals, y = q_25),
              data = quant, linetype = 2, linewidth = 0.3, colour = "blue") +
            ggplot2::geom_ribbon(
              data = quant, ggplot2::aes(x = x_vals, ymin = q_25, ymax = q_975),
              inherit.aes = FALSE, fill = "blue", alpha = 0.1) +
            ggplot2::geom_line(
              mapping = ggplot2::aes(x = x_vals, y = q_50), data = quant,
              linetype = 1, linewidth = 0.6, colour = "blue") +
            ggplot2::geom_text(
              data = data_trend,
              mapping = ggplot2::aes(x = X, y = Y, label = trend_2),
              colour = "darkgrey", size = 2.5, vjust = 1.4, hjust = -0.05) +
            ggplot2::facet_grid(~nfv, labeller = facet_label) +
            ggplot2::scale_y_continuous(
              limits = c(-1, plot_max),
              oob = scales::squish, expand = c(0, 0)) +
            ggplot2::scale_x_continuous(expand = c(0.0125, 0.0125)) +
            ggplot2::xlab(variable_2) +
            ggplot2::ylab(y_axis_label) +
            ggplot2::labs(title = title_text, caption = plot_caption) +
            ggplot2::theme_bw() +
            plot_theme

          # Using ggplot2::ggsave directly does not show non-ascii characters
          # correctly
          ragg::agg_jpeg(
            filename = fs::path(
              path_rc_sr,
              paste0("rc_sr_", variable, "_coords_", coords, ".jpeg")),
            width = 20, height = 12.5, res = 600, quality = 100, units = "cm")
          print(plot)
          grDevices::dev.off()


          # Back-transforming variables
          if (variable %in% c(
            "efforts_log", "road_rail_log", "habitat_log", "rivers_log")) {

            ecokit::cat_time(
              paste0(variable, " - coords = ", coords, " - original scale"),
              level = 2L, verbose = verbose, cat_timestamp = FALSE)

            observed_2 <- dplyr::mutate(observed, x_vals = 10 ^ x_vals - 0.1)
            quant_2 <- dplyr::mutate(quant, x_vals = 10 ^ x_vals - 0.1)

            # Maximum value on the y-axis
            plot_max <- max(observed_2$pred, quant_2$q_975) * 1.05

            # variable long name (x-axis label)
            variable_2 <- dplyr::case_when(
              variable == "road_rail_log" ~ "Road + Rail intensity",
              variable == "habitat_log" ~ "% habitat coverage",
              variable == "rivers_log" ~ "River length",
              variable == "efforts_log" ~ "Sampling efforts",
              .default = variable) %>%
              paste0(
                "<span style='font-size: 10pt;'><b>", .,
                "</b></span><span style='font-size: 7pt;'> (original scale)",
                "</span>")

            plot_2 <- ggplot2::ggplot(
              data = observed_2, mapping = ggplot2::aes(x = x_vals, y = pred),
              environment = emptyenv()) +
              ggplot2::geom_jitter(
                shape = 16, width = 0, height = 0.02, alpha = 0.2, size = 0.75,
                colour = "darkgrey") +
              ggplot2::geom_line(
                ggplot2::aes(x = x_vals, y = q_975), data = quant_2,
                linetype = 2, linewidth = 0.3, colour = "blue") +
              ggplot2::geom_line(
                ggplot2::aes(x = x_vals, y = q_25), data = quant_2,
                linetype = 2, linewidth = 0.3, colour = "blue") +
              ggplot2::geom_ribbon(
                data = quant_2,
                ggplot2::aes(x = x_vals, ymin = q_25, ymax = q_975),
                inherit.aes = FALSE, fill = "blue", alpha = 0.1) +
              ggplot2::geom_line(
                mapping = ggplot2::aes(x = x_vals, y = q_50), data = quant_2,
                linetype = 1, linewidth = 0.6, colour = "blue") +
              ggplot2::geom_text(
                data = data_trend,
                mapping = ggplot2::aes(x = X, y = Y, label = trend_2),
                colour = "darkgrey", size = 2.5, vjust = 1.4, hjust = -0.05) +
              ggplot2::facet_grid(~nfv, labeller = facet_label) +
              ggplot2::scale_y_continuous(
                limits = c(-1, plot_max),
                oob = scales::squish, expand = c(0, 0)) +
              ggplot2::scale_x_continuous(
                expand = c(0.0125, 0.0125), labels = scales::comma) +
              ggplot2::xlab(variable_2) +
              ggplot2::ylab(y_axis_label) +
              ggplot2::labs(title = title_text, caption = plot_caption) +
              ggplot2::theme_bw() +
              plot_theme

            # Using ggplot2::ggsave directly does not show non-ascii characters
            # correctly
            ragg::agg_jpeg(
              filename = fs::path(
                path_rc_sr,
                paste0(
                  "rc_sr_", variable, "_coords_", coords,
                  "_original_scale.jpeg")),
              width = 20, height = 12.5, res = 600, quality = 100, units = "cm")
            print(plot_2)
            grDevices::dev.off()
          }

          invisible(NULL)
        }))

  save(
    sr_data_all,
    file = fs::path(path_rc_sr, "sr_data_all.RData"))

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time,
    prefix = "Plotting response curves for species richness took ",
    verbose = verbose)

  # # ..................................................................... ###

  return(invisible(NULL))
}
