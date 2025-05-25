## |------------------------------------------------------------------------| #
# resp_curv_plot_species_all ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname response_curves
#' @name response_curves
#' @order 3
#' @author Ahmed El-Gabbas

resp_curv_plot_species_all <- function(
    model_dir = NULL, n_cores = 8L, strategy = "multisession",
    return_data = FALSE, plotting_alpha = 0.3) {

  # # ..................................................................... ###

  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores <= 0) {
    ecokit::stop_ctx(
      "n_cores must be a single positive integer.", n_cores = n_cores,
      include_backtrace = TRUE)
  }

  if (!is.character(strategy)) {
    ecokit::stop_ctx(
      "`strategy` must be a character vector",
      strategy = strategy, class_strategy = class(strategy))
  }
  if (strategy == "sequential") {
    n_cores <- 1L
  }
  if (length(strategy) != 1L) {
    ecokit::stop_ctx(
      "`strategy` must be a character vector of length 1",
      strategy = strategy, length_strategy = length(strategy))
  }
  valid_strategy <- c("sequential", "multisession", "multicore", "cluster")
  if (!strategy %in% valid_strategy) {
    ecokit::stop_ctx("Invalid `strategy` value", strategy = strategy)
  }

  ecokit::cat_time("Plotting species response curves")

  .start_time <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Coords <- RC_Path_Prob <- NFV <- Data <- DT <- Variable <- Variable2 <-
    VarDesc <- VarDesc2 <- NULL

  # # ..................................................................... ###

  # Check arguments

  ecokit::cat_time("Check arguments", level = 1L)

  if (is.null(model_dir)) {
    ecokit::stop_ctx(
      "`model_dir` cannot be NULL", model_dir = model_dir,
      include_backtrace = TRUE)
  }

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = "model_dir")
  ecokit::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = c("n_cores", "plotting_alpha"))
  rm(AllArgs, envir = environment())

  if (plotting_alpha < 0 || plotting_alpha > 1) {
    ecokit::stop_ctx(
      "`plotting_alpha` must be between 0 and 1",
      plotting_alpha = plotting_alpha, include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c("ecokit", "dplyr", "magrittr"), strategy = strategy)

  # # ..................................................................... ###

  Path_RC_DT <- fs::path(model_dir, "Model_Postprocessing", "RespCurv_DT")
  Path_RC_All <- fs::path(model_dir, "Model_Postprocessing", "RespCurv_All")

  if (!dir.exists(Path_RC_DT)) {
    ecokit::stop_ctx(
      "Response curve data subfolder is missing.", Path_RC_DT = Path_RC_DT,
      include_backtrace = TRUE)
  }

  fs::dir_create(Path_RC_All)

  # # ..................................................................... ###

  # Loading & processing species response curve data in parallel

  ecokit::cat_time(
    "Loading & processing species response curve data in parallel", level = 1L)

  Sp_DT_All <- fs::path(Path_RC_DT, "ResCurvDT.RData") %>%
    ecokit::load_as() %>%
    dplyr::select(Coords, RC_Path_Prob)

  if (n_cores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = min(n_cores, nrow(Sp_DT_All)), level = 1L,
      future_max_size = 800L, strategy = strategy)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  Sp_DT_All <- Sp_DT_All %>%
    dplyr::mutate(
      Data = furrr::future_map(
        .x = RC_Path_Prob,
        .f = ~ {
          ecokit::load_as(.x) %>%
            dplyr::select(Variable, NFV, Species, PlotData_Quant)
        },
        .options = furrr::furrr_options(
          seed = TRUE, chunk_size = 1, packages = pkg_to_export))) %>%
    tidyr::unnest(Data) %>%
    dplyr::select(-RC_Path_Prob) %>%
    dplyr::slice(gtools::mixedorder(Variable)) %>%
    tidyr::nest(DT = -c(NFV, Coords))

  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 1L)
    future::plan("sequential", gc = TRUE)
  }
  invisible(gc())

  # # ..................................................................... ###

  # Plot all species response curves

  ecokit::cat_time("Plot all species response curves", level = 1L)

  Plots <- purrr::map_dfr(
    .x = seq_len(nrow(Sp_DT_All)),
    .f = function(ID) {

      NFV <- Sp_DT_All$NFV[[ID]]
      Coords <- Sp_DT_All$Coords[[ID]]

      FilePrefix <- paste0("RespCurv_All_NFV_", NFV, "_Coords_", Coords)
      path_JPEG <- fs::path(Path_RC_All, paste0(FilePrefix, ".jpeg"))

      DT <- Sp_DT_All$DT[[ID]] %>%
        dplyr::mutate(
          VarDesc = dplyr::case_when(
            startsWith(Variable, "bio") ~ stringr::str_to_sentence(Variable),
            Variable == "RoadRailLog" ~ "Road + Rail intensity",
            Variable == "EffortsLog" ~ "Sampling efforts",
            Variable == "RiversLog" ~ "River length",
            Variable == "HabLog" ~ "% habitat coverage",
            .default = Variable),
          VarDesc = paste0(
            "<span style='font-size: 10pt;'><b>", VarDesc, "</b></span>"),

          VarDesc2 = dplyr::case_when(
            Variable == "bio1" ~ "annual mean temperature",
            Variable == "bio2" ~ "mean diurnal range",
            Variable == "bio3" ~ "isothermality (bio2/bio7) (&times;100)",
            Variable == "bio4" ~ "temperature seasonality",
            Variable == "bio5" ~ "max temperature of warmest month",
            Variable == "bio6" ~ "temperature of the coldest month",
            Variable == "bio7" ~ "temperature annual range (bio5-bio6)",
            Variable == "bio8" ~ "temperatures of the wettest quarter",
            Variable == "bio9" ~ "mean temperature of driest quarter",
            Variable == "bio10" ~ "mean temperature of warmest quarter",
            Variable == "bio11" ~ "mean temperature of coldest quarter",
            Variable == "bio12" ~ "annual precipitation amount",
            Variable == "bio13" ~ "precipitation of wettest month",
            Variable == "bio14" ~ "precipitation of driest month",
            Variable == "bio15" ~ "precipitation seasonality",
            Variable == "bio16" ~ "precipitation of wettest quarter",
            Variable == "bio17" ~ "precipitation of driest quarter",
            Variable == "bio18" ~ "precipitation of the warmest quarter",
            Variable == "bio19" ~ "precipitation of coldest quarter",
            Variable == "npp" ~ "net primary productivity",
            Variable == "RiversLog" ~ " (log<sub>10</sub>(x + 0.1))",
            Variable == "RoadRailLog" ~ " (log<sub>10</sub>(x + 0.1))",
            Variable == "EffortsLog" ~ " (log<sub>10</sub>(x + 0.1))",
            Variable == "HabLog" ~ " (log<sub>10</sub>(x + 0.1))",
            .default = Variable),
          VarDesc2 = paste0(
            "<span style='font-size: 8pt;'>", VarDesc2, "</span>"),

          Variable2 = factor(Variable, levels = unique(Variable))) %>%
        dplyr::group_split(Variable2)

      SubTitleTxt <- dplyr::if_else(
        NFV == 1,
        paste0(
          "Non-focal variables are set to most likely value <i>",
          "[non.focalVariables = 1]</i>"),
        paste0(
          "Non-focal variables = most likely value given ",
          "value of focal variable <i>[non.focalVariables = 2]</i>"))
      Caption <- dplyr::if_else(
        Coords == "c", "Mean coordinates", "No spatial dependence")
      Caption <- paste0(Caption, " --- ", SubTitleTxt)

      if (length(DT) <= 9) {
        NR <- NC <- 3
        plot_width <- 24
        plot_height <- 22
      } else {
        NR <- 3
        NC <- 4
        plot_width <- 30
        plot_height <- 22
      }

      Plots <- purrr::map(
        .x = DT,
        .f = ~ {
          .x %>%
            dplyr::select(Species, PlotData_Quant) %>%
            tidyr::unnest("PlotData_Quant") %>%
            dplyr::filter(Quantile == 0.5) %>%
            ggplot2::ggplot(
              mapping = ggplot2::aes(x = XVals, y = Pred, group = Species)) +
            ggplot2::geom_line(
              linetype = 1, linewidth = 0.3,
              colour = "blue", alpha = plotting_alpha) +
            ggplot2::scale_y_continuous(
              limits = c(-0.01, 1.01), oob = scales::squish,
              expand = c(0, 0)) +
            ggplot2::scale_x_continuous(expand = c(0.015, 0.015)) +
            ggplot2::labs(
              x = NULL, y = NULL,
              title = .x$VarDesc[[1]], subtitle = .x$VarDesc2[[1]]) +
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
        }) %>%
        patchwork::wrap_plots(nrow = NR, ncol = NC) +
        patchwork::plot_annotation(
          title = "Response curves for all species",
          caption = Caption,
          theme = ggplot2::theme(
            plot.margin = ggplot2::margin(0.1, 0, 0, 0, "cm"),
            plot.title = ggtext::element_markdown(face = "bold"),
            plot.caption = ggtext::element_markdown(
              size = 12, color = "grey", hjust = 0))) +
        patchwork::plot_layout(axes = "collect")

      Plots <- patchwork::wrap_elements(Plots) +
        ggplot2::labs(tag = "<b>Predicted habitat suitability</b>") +
        ggplot2::theme(
          plot.tag = ggtext::element_markdown(
            size = ggplot2::rel(1.25), angle = 90),
          plot.tag.position = "left")


      # Using ggplot2::ggsave directly does not show non-ascii characters
      # correctly
      ragg::agg_jpeg(
        filename = path_JPEG, width = plot_width, height = plot_height,
        res = 600, quality = 100, units = "cm")
      print(Plots)
      grDevices::dev.off()

      OutDF <- tibble::tibble(
        path_JPEG = path_JPEG, plot_height = plot_height,
        plot_width = plot_width)

      return(OutDF)
    })

  # # ..................................................................... ###

  # Save data
  ecokit::cat_time("Save data", level = 1L)

  Sp_DT_All <- dplyr::select(Sp_DT_All, -DT) %>%
    dplyr::bind_cols(Plots = Plots)

  save(Sp_DT_All, file = fs::path(Path_RC_All, "Sp_DT_All.RData"))

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time,
    prefix = "Plotting all species response curves took ", level = 1L)

  # # ..................................................................... ###

  if (return_data) {
    return(Sp_DT_All)
  } else {
    return(invisible(NULL))
  }
}
