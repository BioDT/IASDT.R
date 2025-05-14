
## |------------------------------------------------------------------------| #
# resp_curv_plot_species ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname response_curves
#' @name response_curves
#' @order 2
#' @author Ahmed El-Gabbas

resp_curv_plot_species <- function(
    model_dir = NULL, n_cores = 20, strategy = "future::multicore",
    env_file = ".env", return_data = FALSE) {

  # # ..................................................................... ###

  ecokit::cat_time("Plotting species response curves")

  .start_time <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_PA <- NCells_Naturalized <- NFV <- Coords <- Species <- Prefix <-
    Data <- RC_Path_Prob <- Variable <- IAS_ID <- VarDesc <- VarDesc2 <- NULL

  # # ..................................................................... ###

  # Check arguments

  ecokit::cat_time("Check arguments")

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
    args_to_check = c("model_dir", "env_file", "strategy"))

  ecokit::check_args(
    args_all = AllArgs, args_type = "numeric", args_to_check = "n_cores")
  rm(AllArgs, envir = environment())

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
  if (strategy == "future::sequential") {
    n_cores <- 1L
  }
  if (length(strategy) != 1L) {
    ecokit::stop_ctx(
      "`strategy` must be a character vector of length 1",
      strategy = strategy, length_strategy = length(strategy))
  }
  valid_strategy <- c(
    "future::sequential", "future::multisession", "future::multicore",
    "future::cluster")
  if (!strategy %in% valid_strategy) {
    ecokit::stop_ctx("Invalid `strategy` value", strategy = strategy)
  }

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c(
      "dplyr", "purrr", "tidyr", "gtools", "ggtext", "patchwork",
      "ggplot2", "tibble", "ecokit", "ragg", "stringr", "scales"),
    strategy = strategy)

  # # ..................................................................... ###

  Path_RC_DT <- fs::path(model_dir, "Model_Postprocessing", "RespCurv_DT")
  if (!dir.exists(Path_RC_DT)) {
    ecokit::stop_ctx(
      "Response curve data subfolder is missing.", Path_RC_DT = Path_RC_DT,
      include_backtrace = TRUE)
  }
  Path_RC_Sp <- fs::path(model_dir, "Model_Postprocessing", "RespCurv_Sp")
  Path_RC_Sp_DT <- fs::path(model_dir, "Model_Postprocessing", "RespCurv_Sp_DT")
  fs::dir_create(c(Path_RC_Sp, Path_RC_Sp_DT))

  # # ..................................................................... ###

  # # Load species summary
  ecokit::cat_time("Load species summary")

  EnvVars2Read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "Path_PA", "DP_R_PA", TRUE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  SpSummary <- fs::path(Path_PA, "Sp_PA_Summary_DF.csv")
  if (!file.exists(SpSummary)) {
    ecokit::stop_ctx(
      "SpSummary file does not exist", SpSummary = SpSummary,
      include_backtrace = TRUE)
  }

  SpSummary <- readr::read_csv(
    file = SpSummary, show_col_types = FALSE, progress = FALSE) %>%
    dplyr::select(tidyselect::all_of(c("IAS_ID", "NCells_Naturalized"))) %>%
    dplyr::rename(NCells = NCells_Naturalized)

  # # ..................................................................... ###

  # Load species names
  ecokit::cat_time("Load species names")
  SpeciesNames <- IASDT.R::get_species_name(env_file = env_file)

  # # ..................................................................... ###

  # Prepare species-specific data in parallel

  ecokit::cat_time("Prepare species-specific data in parallel")

  if (n_cores == 1) {
    future::plan("future::sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = n_cores, level = 1L, future_max_size = 800L,
      strategy = strategy)
    withr::defer(future::plan("future::sequential", gc = TRUE))
  }

  ecokit::cat_time("Processing in parallel", level = 1L)

  Sp_DT_All <- fs::path(Path_RC_DT, "ResCurvDT.RData") %>%
    ecokit::load_as() %>%
    dplyr::select(tidyselect::all_of(c("Coords", "RC_Path_Prob"))) %>%
    dplyr::mutate(
      Data = furrr::future_map(
        .x = RC_Path_Prob, .f = ecokit::load_as,
        .options = furrr::furrr_options(
          globals = "Path_RC_DT", packages = pkg_to_export))) %>%
    tidyr::unnest(Data) %>%
    dplyr::select(-RC_Path_Prob) %>%
    tidyr::nest(
      DT = tidyselect::everything(), .by = c(NFV, Coords, Species)) %>%
    dplyr::mutate(IAS_ID = as.numeric(stringr::str_remove(Species, "^Sp_"))) %>%
    dplyr::left_join(SpSummary, by = "IAS_ID") %>%
    dplyr::mutate(
      Prefix = paste0(Species, "_NFV_", NFV, "_Coords_", Coords),
      path_JPEG_fixed = fs::path(Path_RC_Sp, paste0(Prefix, "_Fixed.jpeg")),
      path_JPEG_free = fs::path(Path_RC_Sp, paste0(Prefix, "_Free.jpeg")),
      Path_Sp_DT = fs::path(Path_RC_Sp_DT, paste0(Prefix, ".qs2")))

  # stopping the cluster
  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 2L)
    future::plan("future::sequential", gc = TRUE)
  }

  ecokit::cat_time("Export species-specific data", level = 1L)
  purrr::walk(
    .x = seq_len(nrow(Sp_DT_All)),
    .f = function(ID) {
      DT <- dplyr::slice(Sp_DT_All, ID)
      if (isFALSE(ecokit::check_data(DT$Path_Sp_DT, warning = FALSE))) {
        ecokit::save_as(object = DT, out_path = DT$Path_Sp_DT)
      }
    })

  Sp_DT_All <- gtools::mixedsort(Sp_DT_All$Path_Sp_DT)
  invisible(gc())

  # # ..................................................................... ###

  ecokit::cat_time("Plotting species-specific data")

  if (n_cores == 1) {
    future::plan("future::sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = n_cores, level = 1L, future_max_size = 800L,
      strategy = strategy)
    withr::defer(future::plan("future::sequential", gc = TRUE))
  }

  ecokit::cat_time("Plotting in parallel", level = 1L)

  Plots <- future.apply::future_lapply(
    X = Sp_DT_All,
    FUN = function(RC_File) {

      DT <- ecokit::load_as(RC_File)

      Coords <- DT$Coords
      Species <- DT$Species
      NCells <- DT$NCells
      NFV <- DT$NFV
      path_JPEG_fixed <- DT$path_JPEG_fixed
      path_JPEG_free <- DT$path_JPEG_free

      DT <- dplyr::select(DT, DT) %>%
        tidyr::unnest(DT) %>%
        dplyr::slice(gtools::mixedorder(Variable)) %>%
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
            "<span style='font-size: 8pt;'>", VarDesc2, "</span>"))

      # nolint start
      Species2 <- dplyr::filter(SpeciesNames, IAS_ID == !!Species)
      Species_name <- Species2$Species_name
      Species_ID <- stringr::str_remove(Species, "^Sp_")
      Class <- Species2$Class
      Order <- Species2$Order
      Family <- Species2$Family
      TitleTxt <- stringr::str_glue(
        '<span style="font-size:13pt;"><b> Response curves of </b></span>\\
        <span style="color:blue; font-size:13pt;">\\
        <b><i>{Species_name}, "</i></b></span>\\
        <span style="font-size:8pt;"> (\\
        <b>Class:</b> {Class} &#8212; <b>Order:</b> {Order} &#8212; \\
        <b>Family:</b> {Family} &#8212; <b>ID:</b> {Species_ID} \\
        &#8212; <b># presence grids:</b> {NCells})</span>')
      # nolint end

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

      if (nrow(DT) <= 9) {
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
        .x = seq_len(nrow(DT)),
        .f = ~ {

          PositiveTrendProb <- DT$PositiveTrendProb[[.x]]
          Trend <- tibble::tibble(
            Trend2 = stringr::str_glue(
              "\n     Pr[pred(Var=max)] > Pr[pred(Var=min)] = \\
              {round(PositiveTrendProb, 2)}"),
            X = -Inf, Y = Inf)

          Quant <- DT$PlotData_Quant[[.x]] %>%
            tidyr::pivot_wider(
              id_cols = XVals, names_from = Quantile, values_from = Pred) %>%
            stats::setNames(c("XVals", "Q25", "Q50", "Q975"))

          Rug_0 <- dplyr::filter(DT$Observed_PA[[.x]], Pred == 0)
          Rug_1 <- dplyr::filter(DT$Observed_PA[[.x]], Pred == 1)

          PlottingTheme <- ggplot2::theme_bw() +
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
          Plot_Fixed <- ggplot2::ggplot(
            data = Quant, mapping = ggplot2::aes(x = XVals)) +
            ggplot2::geom_line(
              ggplot2::aes(y = Q975), data = Quant,
              linetype = 2, linewidth = 0.3, colour = "blue") +
            ggplot2::geom_line(
              ggplot2::aes(y = Q25), data = Quant,
              linetype = 2, linewidth = 0.3, colour = "blue") +
            ggplot2::geom_ribbon(
              mapping = ggplot2::aes(ymin = Q25, ymax = Q975),
              data = Quant, fill = "blue", alpha = 0.1) +
            ggplot2::geom_line(
              mapping = ggplot2::aes(y = Q50), data = Quant,
              linetype = 1, linewidth = 0.6, colour = "blue") +
            ggplot2::geom_rug(
              sides = "t", data = Rug_1, ggplot2::aes(x = XVals),
              color = "blue", linewidth = 0.025, alpha = 0.25,
              length = grid::unit(0.03, "npc")) +
            ggplot2::geom_rug(
              sides = "b", data = Rug_0, ggplot2::aes(x = XVals),
              color = "red", linewidth = 0.025, alpha = 0.25,
              length = grid::unit(0.03, "npc")) +
            ggplot2::geom_text(
              data = Trend, vjust = 0.5, hjust = -0.05,
              mapping = ggplot2::aes(x = X, y = Y, label = Trend2),
              colour = "grey30", size = 2.75) +
            ggplot2::scale_y_continuous(
              limits = c(0, 1), oob = scales::squish_infinite) +
            ggplot2::scale_x_continuous(expand = c(0.015, 0.015)) +
            ggplot2::labs(
              x = NULL, y = NULL, title = DT$VarDesc[[.x]],
              subtitle = DT$VarDesc2[[.x]]) +
            PlottingTheme

          # Free y-axis
          Plot_Free <- ggplot2::ggplot(
            data = Quant, mapping = ggplot2::aes(x = XVals)) +
            ggplot2::geom_line(
              ggplot2::aes(y = Q975), data = Quant,
              linetype = 2, linewidth = 0.3, colour = "blue") +
            ggplot2::geom_line(
              ggplot2::aes(y = Q25), data = Quant,
              linetype = 2, linewidth = 0.3, colour = "blue") +
            ggplot2::geom_ribbon(
              mapping = ggplot2::aes(ymin = Q25, ymax = Q975),
              data = Quant, fill = "blue", alpha = 0.1) +
            ggplot2::geom_line(
              mapping = ggplot2::aes(y = Q50), data = Quant,
              linetype = 1, linewidth = 0.6, colour = "blue") +
            ggplot2::geom_rug(
              sides = "t", data = Rug_1, ggplot2::aes(x = XVals),
              color = "blue", linewidth = 0.025, alpha = 0.25,
              length = grid::unit(0.03, "npc")) +
            ggplot2::geom_rug(
              sides = "b", data = Rug_0, ggplot2::aes(x = XVals),
              color = "red", linewidth = 0.025, alpha = 0.25,
              length = grid::unit(0.03, "npc")) +
            ggplot2::geom_text(
              data = Trend, vjust = 0.5, hjust = -0.05,
              mapping = ggplot2::aes(x = X, y = Y, label = Trend2),
              colour = "grey30", size = 2.75) +
            ggplot2::scale_y_continuous(oob = scales::squish_infinite) +
            ggplot2::scale_x_continuous(expand = c(0.015, 0.015)) +
            ggplot2::labs(
              x = NULL, y = NULL, title = DT$VarDesc[[.x]],
              subtitle = DT$VarDesc2[[.x]]) +
            PlottingTheme +
            ggplot2::theme(
              axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5))

          return(tibble::tibble(
            Plot_Fixed = list(Plot_Fixed), Plot_Free = list(Plot_Free)))
        }) %>%
        dplyr::bind_rows()

      PlottingTheme2 <- ggplot2::theme(
        plot.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
        plot.title = ggtext::element_markdown(
          hjust = 0, margin = ggplot2::margin(0, 0, -0.125, 0, "cm")),
        plot.caption = ggtext::element_markdown(
          size = 12, color = "grey", hjust = 0))

      # Fixed y-axis
      Plot_Fixed <- patchwork::wrap_plots(
        Plots$Plot_Fixed, nrow = NR, ncol = NC) +
        patchwork::plot_annotation(
          title = stringr::str_glue(
            "{TitleTxt}<span style = 'color:#ffffff;'>.......................\\
            ..........</span><span style='font-size:10pt; color:grey;'>Fixed \\
            y-axis</span>"),
          caption = Caption, theme = PlottingTheme2) +
        patchwork::plot_layout(axes = "collect")

      Plot_Fixed <- patchwork::wrap_elements(panel = Plot_Fixed) +
        ggplot2::labs(tag = "<b>Predicted habitat suitability</b>") +
        ggplot2::theme(
          plot.tag = ggtext::element_markdown(
            size = ggplot2::rel(1), angle = 90),
          plot.tag.position = "left")

      ragg::agg_jpeg(
        filename = path_JPEG_fixed, width = plot_width, height = plot_height,
        res = 600, quality = 100, units = "cm")
      print(Plot_Fixed)
      grDevices::dev.off()


      # Free y-axis
      Plot_Free <- patchwork::wrap_plots(
        Plots$Plot_Free, nrow = NR, ncol = NC) +
        patchwork::plot_annotation(
          title = stringr::str_glue(
            "{TitleTxt}<span style = 'color:#ffffff;'>...........</span><span \\
            style='font-size:10pt; color:grey;'>Free y-axis</span>"),
          caption = Caption, theme = PlottingTheme2) +
        patchwork::plot_layout(axes = "collect")

      Plot_Free <- patchwork::wrap_elements(panel = Plot_Free) +
        ggplot2::labs(tag = "<b>Predicted habitat suitability</b>") +
        ggplot2::theme(
          plot.tag = ggtext::element_markdown(
            size = ggplot2::rel(1), angle = 90),
          plot.tag.position = "left")

      ragg::agg_jpeg(
        filename = path_JPEG_free, width = plot_width, height = plot_height,
        res = 600, quality = 100, units = "cm")
      print(Plot_Free)
      grDevices::dev.off()

      OutDF <- tibble::tibble(
        Coords = Coords, Species = Species, NCells = NCells, NFV = NFV,
        path_JPEG_fixed = path_JPEG_fixed, path_JPEG_free = path_JPEG_free,
        plot_height = plot_height, plot_width = plot_width)

      return(OutDF)

    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = pkg_to_export,
    future.globals = c("SpeciesNames", "Sp_DT_All")) %>%
    dplyr::bind_rows()

  # stopping the cluster
  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 2L)
    future::plan("future::sequential", gc = TRUE)
  }
  invisible(gc())

  # # ..................................................................... ###

  # Save data
  ecokit::cat_time("Save data")
  save(Plots, file = fs::path(Path_RC_Sp, "Sp_DT_All.RData"))

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time, prefix = "Plotting species response curves took ")

  # # ..................................................................... ###

  if (return_data) {
    return(Sp_DT_All)
  } else {
    return(invisible(NULL))
  }
}
