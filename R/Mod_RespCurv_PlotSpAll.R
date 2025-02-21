## |------------------------------------------------------------------------| #
# RespCurv_PlotSpAll ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname Response_curves
#' @name Response_curves
#' @order 3
#' @author Ahmed El-Gabbas

RespCurv_PlotSpAll <- function(
    ModelDir = NULL, NCores = 8L, ReturnData = FALSE, PlottingAlpha = 0.3) {

  # # ..................................................................... ###

  IASDT.R::CatTime("Plotting species response curves")

  .StartTime <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Coords <- RC_Path_Prob <- NFV <- Data <- DT <- Variable <- Variable2 <-
    VarDesc <- VarDesc2 <- NULL

  # # ..................................................................... ###

  # Check arguments

  IASDT.R::CatTime("Check arguments", Level = 1)

  if (is.null(ModelDir)) {
    stop("`ModelDir` cannot be NULL", call. = FALSE)
  }

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "character", Args = "ModelDir")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric", Args = c("NCores", "PlottingAlpha"))
  rm(AllArgs, envir = environment())

  if (PlottingAlpha < 0 || PlottingAlpha > 1) {
    stop("`PlottingAlpha` must be between 0 and 1", call. = FALSE)
  }

  # # ..................................................................... ###


  Path_RC_DT <- IASDT.R::Path(ModelDir, "Model_Postprocessing", "RespCurv_DT")
  Path_RC_All <- IASDT.R::Path(ModelDir, "Model_Postprocessing", "RespCurv_All")

  if (!dir.exists(Path_RC_DT)) {
    stop("Response curve data subfolder is missing.", call. = FALSE)
  }

  fs::dir_create(Path_RC_All)

  # # ..................................................................... ###

  # Loading & processing species response curve data on parallel

  IASDT.R::CatTime(
    "Loading & processing species response curve data on parallel", Level = 1)

  Sp_DT_All <- IASDT.R::Path(Path_RC_DT, "ResCurvDT.RData") %>%
    IASDT.R::LoadAs() %>%
    dplyr::select(Coords, RC_Path_Prob)

  if (NCores == 1) {
    future::plan("future::sequential", gc = TRUE)
  } else {
    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
      future.seed = TRUE)
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
              Variable, NFV, Species, PlotData_Quant)
        },
        .options = furrr::furrr_options(seed = TRUE, chunk_size = 1))) %>%
    tidyr::unnest(Data) %>%
    dplyr::select(-RC_Path_Prob) %>%
    dplyr::slice(gtools::mixedorder(Variable)) %>%
    tidyr::nest(DT = -c(NFV, Coords))


  snow::stopCluster(c1)
  future::plan("future::sequential", gc = TRUE)
  invisible(gc())

  # # ..................................................................... ###

  # Plot all species response curves

  IASDT.R::CatTime("Plot all species response curves", Level = 1)

  Plots <- purrr::map_dfr(
    .x = seq_len(nrow(Sp_DT_All)),
    .f = function(ID) {

      NFV <- Sp_DT_All$NFV[[ID]]
      Coords <- Sp_DT_All$Coords[[ID]]

      FilePrefix <- paste0("RespCurv_All_NFV_", NFV, "_Coords_", Coords)
      Path_JPEG <- IASDT.R::Path(Path_RC_All, paste0(FilePrefix, ".jpeg"))

      DT <- Sp_DT_All$DT[[ID]] %>%
        dplyr::mutate(
          VarDesc = dplyr::case_when(
            stringr::str_detect(Variable, "^bio") ~
              stringr::str_to_sentence(Variable),
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
          .x %>%
            dplyr::select(Species, PlotData_Quant) %>%
            tidyr::unnest("PlotData_Quant") %>%
            dplyr::filter(Quantile == 0.5) %>%
            ggplot2::ggplot(
              mapping = ggplot2::aes(x = XVals, y = Pred, group = Species)) +
            ggplot2::geom_line(
              linetype = 1, linewidth = 0.3,
              colour = "blue", alpha = PlottingAlpha) +
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
        filename = Path_JPEG, width = PlotWidth, height = PlotHeight,
        res = 600, quality = 100, units = "cm")
      print(Plots)
      grDevices::dev.off()

      OutDF <- tibble::tibble(
        Path_JPEG = Path_JPEG, PlotHeight = PlotHeight, PlotWidth = PlotWidth)

      return(OutDF)
    })

  # # ..................................................................... ###

  # Save data
  IASDT.R::CatTime("Save data", Level = 1)

  Sp_DT_All <- dplyr::select(Sp_DT_All, -DT) %>%
    dplyr::bind_cols(Plots = Plots)

  save(Sp_DT_All, file = IASDT.R::Path(Path_RC_All, "Sp_DT_All.RData"))

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
