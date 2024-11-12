## |------------------------------------------------------------------------| #
# RespCurv_PlotSpAll ----
## |------------------------------------------------------------------------| #

#' Plot all species response curves
#'
#' Generates and saves response curve plots for all species together in a single
#' plot.
#' @param Path_Postprocessing String. The path to the directory for model
#'   post-processing. The function reads data from the `RespCurv_Sp`
#'   sub-directory, created by [RespCurv_PrepData].
#' @param NCores Integer. Number of cores to use for parallel processing.
#'   Defaults to 20.
#' @param ReturnData Logical. Indicates whether the output data be returned as
#'   an R object.
#' @param PlotAlpha Numeric. The alpha value (transparency) for the response
#'   curve lines. Defaults to 0.3.
#' @export
#' @name RespCurv_PlotSpAll

RespCurv_PlotSpAll <- function(
    Path_Postprocessing = NULL, NCores = 20,
    ReturnData = FALSE, PlotAlpha = 0.3) {

  # # ..................................................................... ###

  IASDT.R::CatTime("Plotting species response curves")

  .StartTime <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Coords <- RC_Path_Prob <- NFV <- Data <- DT <- Variable <- NULL

  # # ..................................................................... ###

  # Check arguments

  IASDT.R::CatTime("Check arguments", Level = 1)

  if (is.null(Path_Postprocessing)) {
    stop("`Path_Postprocessing` cannot be NULL", call. = FALSE)
  }

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character", Args = "Path_Postprocessing")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric", Args = c("NCores", "PlotAlpha"))
  rm(AllArgs)

  if (PlotAlpha < 0 || PlotAlpha > 1) {
    stop("`PlotAlpha` must be between 0 and 1", call. = FALSE)
  }

  # # ..................................................................... ###

  Path_RespCurv_All <- file.path(Path_Postprocessing, "RespCurv_All")
  fs::dir_create(c(Path_RespCurv_All))

  # # ..................................................................... ###

  # Loading / process species response curve data on parallel

  IASDT.R::CatTime(
    "Loading / process species response curve data on parallel", Level = 1)

  Sp_DT_All <- file.path(
    Path_Postprocessing, "RespCurv_DT", "ResCurvDT.RData") %>%
    IASDT.R::LoadAs() %>%
    dplyr::select(Coords, RC_Path_Prob)

  withr::local_options(
    future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE)

  if (NCores == 1) {
    future::plan("future::sequential", gc = TRUE)
  } else {
    c1 <- snow::makeSOCKcluster(min(NCores, nrow(Sp_DT_All)))
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("future::cluster", workers = c1, gc = TRUE)
    on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
  }

  Sp_DT_All <- Sp_DT_All %>%
    dplyr::mutate(
      Data = furrr::future_map(
        .x = RC_Path_Prob,
        .f = ~ {
          IASDT.R::LoadAs(.x) %>%
            dplyr::select(
              Variable, NFV, Species, PlotData_Quant, Observed_PA)
        },
        .options = furrr::furrr_options(seed = TRUE, chunk_size = 1))) %>%
    tidyr::unnest(Data) %>%
    dplyr::select(-RC_Path_Prob) %>%
    tidyr::nest(DT = -c(NFV, Coords))

  snow::stopCluster(c1)
  future::plan("future::sequential", gc = TRUE)
  invisible(gc())

  # # ..................................................................... ###

  # Plot all species response curves

  IASDT.R::CatTime("Plot all species response curves", Level = 1)

  PlotAllSpRC <- function(ID) {

    NFV <- Sp_DT_All$NFV[[ID]]
    Coords <- Sp_DT_All$Coords[[ID]]

    FilePrefix <- paste0("RespCurv_All_NFV_", NFV, "_Coords_", Coords)
    Path_JPEG <- file.path(Path_RespCurv_All, paste0(FilePrefix, ".jpeg"))

    DT <- Sp_DT_All$DT[[ID]] %>%
      dplyr::mutate(
        VarDesc = dplyr::case_when(
          Variable == "bio2" ~ paste0(
            "<span style='font-size: 10pt;'><b>Bio2</b></span>"),
          Variable == "bio4" ~ paste0(
            "<span style='font-size: 10pt;'><b>Bio4</b></span>"),
          Variable == "bio6" ~ paste0(
            "<span style='font-size: 10pt;'><b>Bio6</b></span>"),
          Variable == "bio8" ~ paste0(
            "<span style='font-size: 10pt;'><b>Bio8</b></span>"),
          Variable == "bio12" ~ paste0(
            "<span style='font-size: 10pt;'><b>Bio12</b></span>"),
          Variable == "bio15" ~ paste0(
            "<span style='font-size: 10pt;'><b>Bio15</b></span>"),
          Variable == "bio18" ~ paste0(
            "<span style='font-size: 10pt;'><b>Bio18</b></span>"),
          Variable == "RoadRailLog" ~ paste0(
            "<span style='font-size: 10pt;'><b>Road + Rail intensity</b>",
            "</span>"),
          Variable == "EffortsLog" ~ paste0(
            "<span style='font-size: 10pt;'><b>Sampling efforts</b></span>"),
          Variable == "HabLog" ~ paste0(
            "<span style='font-size: 10pt;'><b>% habitat coverage</b></span>"),
          .default = Variable),

        VarDesc2 = dplyr::case_when(
          Variable == "bio2" ~ paste0(
            "<span style='font-size: 8pt;'>Mean diurnal range</span>"),
          Variable == "bio4" ~ paste0(
            "<span style='font-size: 8pt;'>temperature seasonality</span>"),
          Variable == "bio6" ~ paste0(
            "<span style='font-size: 8pt;'>temperature of the coldest ",
            "month</span>"),
          Variable == "bio8" ~ paste0(
            "<span style='font-size: 8pt;'>temperatures of the wettest ",
            "quarter</span>"),
          Variable == "bio12" ~ paste0(
            "<span style='font-size: 8pt;'>annual precipitation amount",
            "</span>"),
          Variable == "bio15" ~ paste0(
            "<span style='font-size: 8pt;'>precipitation seasonality</span>"),
          Variable == "bio18" ~ paste0(
            "<span style='font-size: 8pt;'>monthly precipitation amount ",
            "of the warmest quarter</span>"),
          Variable == "RoadRailLog" ~ paste0(
            "<span style='font-size: 8pt;'> (log<sub>10",
            "</sub>(x + 0.1))</span>"),
          Variable == "EffortsLog" ~ paste0(
            "<span style='font-size: 8pt;'> (log<sub>10",
            "</sub>(x + 0.1))</span>"),
          Variable == "HabLog" ~ paste0(
            "<span style='font-size: 8pt;'> (log<sub>10",
            "</sub>(x + 0.1))</span>"))) %>%
      dplyr::group_split(Variable)

    SubTitleTxt <- dplyr::if_else(
      NFV == 1,
      paste0(
        "Non-focal variables are set to most likely value <i>",
        "[non.focalVariables = 1]</i>"),
      paste0(
        "Non-focal variables are set to most likely value given ",
        "the value of focal variable <i>[non.focalVariables = 2]</i>"))
    Caption <- dplyr::if_else(
      Coords == "c", "Predictions are made at mean coordinates",
      paste0(
        "Predictions are made for infinite coordinates ",
        "without effect of spatial dependence"))
    Caption <- paste0(Caption, " --- ", SubTitleTxt)

    if (length(DT) <= 9) {
      NR <- NC <- 3
      PlotWidth <- 24
      PlotHeight <- 22
    } else {
      NR <- 3
      NC <- 4
      PlotWidth <- 30
      PlotHeight <- 22
    }

    Plots <- purrr::map(
      .x = DT,
      .f = ~ {

        # Observed_PA <- dplyr::bind_rows(.x$Observed_PA) %>%
        #   dplyr::mutate(
        #     Col = dplyr::if_else(Pred == 1, "red", "darkgreen"),
        #     Pred = dplyr::case_when(
        #       Pred == 1 ~ 0.97, Pred == 0 ~ 0.03, .default = Pred))

        .x %>%
          dplyr::select(Species, PlotData_Quant) %>%
          tidyr::unnest("PlotData_Quant") %>%
          dplyr::filter(Quantile == 0.5) %>%
          ggplot2::ggplot(
            mapping = ggplot2::aes(x = XVals, y = Pred, group = Species)) +
          ggplot2::geom_line(
            linetype = 1, linewidth = 0.3, colour = "blue", alpha = PlotAlpha) +
          # ggplot2::geom_jitter(
          #   data = Observed_PA, mapping = ggplot2::aes(colour = Col),
          #   shape = 19, width = 0, height = 0.02,
          #   alpha = 0.4, size = 1) +
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
    grDevices::jpeg(
      filename = Path_JPEG, width = PlotWidth, height = PlotHeight,
      units = "cm", quality = 100, res = 600)
    print(Plots)
    grDevices::dev.off()

    OutDF <- tibble::tibble(
      Path_JPEG = Path_JPEG, PlotHeight = PlotHeight, PlotWidth = PlotWidth)

    return(OutDF)
  }

  # # ..................................................................... ###

  if (NCores == 1) {
    future::plan("future::sequential", gc = TRUE)
  } else {
    c1 <- snow::makeSOCKcluster(min(NCores, nrow(Sp_DT_All)))
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("future::cluster", workers = c1, gc = TRUE)
    on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
  }

  Plots <- future.apply::future_lapply(
    X = seq_len(nrow(Sp_DT_All)),
    FUN = PlotAllSpRC, future.seed = TRUE,
    future.chunk.size = 1,
    future.packages = c(
      "dplyr", "purrr", "tidyr", "gtools", "ggtext", "patchwork",
      "ggplot2", "tibble", "IASDT.R"),
    future.globals = c(
      "Sp_DT_All", "Path_RespCurv_All", "PlotAllSpRC")) %>%
    dplyr::bind_rows()

  if (NCores > 1) {
    snow::stopCluster(c1)
    future::plan("future::sequential", gc = TRUE)
  }

  invisible(gc())

  # # ..................................................................... ###

  # Save data
  IASDT.R::CatTime("Save data", Level = 1)

  Sp_DT_All <- dplyr::select(Sp_DT_All, -DT) %>%
    dplyr::bind_cols(Plots = Plots)

  save(Sp_DT_All, file = file.path(Path_RespCurv_All, "Sp_DT_All.RData"))

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTime,
    Prefix = "Plotting all species response curves took ", Level = 1)

  # # ..................................................................... ###

  if (ReturnData) {
    return(Sp_DT_All)
  } else {
    return(invisible(NULL))
  }
}
