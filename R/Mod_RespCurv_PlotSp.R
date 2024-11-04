## |------------------------------------------------------------------------| #
# RespCurv_PlotSp ----
## |------------------------------------------------------------------------| #

#' Plot species response curves
#'
#' Generates and saves response curve plots for species based on environmental
#' variables and other factors. Plots show predicted habitat suitability across
#' different values of each variable.
#' @param Path_Postprocessing String. The path to the directory for model
#'   post-processing. The function reads data from the `RespCurv_Sp`
#'   sub-directory, created by [RespCurv_PrepData].
#' @param NCores Integer. Number of cores to use for parallel processing.
#'   Defaults to 20.
#' @param EnvFile String. Path to the environment variables file. Defaults to
#'   ".env".
#' @param FromHPC Logical. Indicates whether the function is being run on an HPC
#'   environment, affecting file path handling. Default: `TRUE`.
#' @param ReturnData Logical. Indicates whether the output data be returned as
#'   an R object.
#' @export
#' @name RespCurv_PlotSp

RespCurv_PlotSp <- function(
    Path_Postprocessing = NULL, NCores = 20, EnvFile = ".env",
    FromHPC = TRUE, ReturnData = FALSE) {

  # # ..................................................................... ###

  IASDT.R::CatTime("Plotting species response curves")

  .StartTime <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  RC_Path_Prob <- NFV <- Species <- Coords <- Data <- PlotData <-
    DT <- Variable <- VarDesc <- Trend2 <- PositiveTrendProb <-
    VarDesc2 <- Pred <- XVals <- PlotData_Quant <- Quantile <-
    Observed_PA <- Col <- Q975 <- Q25 <- Q50 <- X <- Y <- IAS_ID <- NULL

  # # ..................................................................... ###

  # Check arguments

  IASDT.R::CatTime("Check arguments")

  if (is.null(Path_Postprocessing)) {
    stop("`Path_Postprocessing` cannot be NULL", call. = FALSE)
  }

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character",
    Args = c("Path_Postprocessing", "EnvFile"))

  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "numeric", Args = "NCores")
  rm(AllArgs)

  Path_RespCurv_Sp <- file.path(Path_Postprocessing, "RespCurv_Sp")
  Path_GG <- file.path(Path_Postprocessing, "RespCurv_Sp_GG")
  fs::dir_create(c(Path_RespCurv_Sp, Path_GG))

  # # ..................................................................... ###

  # Load species names
  IASDT.R::CatTime("Load species names")
  SpeciesNames <- IASDT.R::GetSpeciesName(EnvFile = EnvFile, FromHPC = FromHPC)

  # # ..................................................................... ###

  # Loading / process species response curve data on parallel

  IASDT.R::CatTime("Loading / process species response curve data on parallel")

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
        .f = ~ dplyr::select(IASDT.R::LoadAs(.x), -PlotData),
        .options = furrr::furrr_options(seed = TRUE, chunk_size = 1))) %>%
    tidyr::unnest(Data) %>%
    dplyr::select(-RC_Path_Prob) %>%
    tidyr::nest(
      DT = tidyselect::everything(), .by = c(NFV, Coords, Species))

  snow::stopCluster(c1)
  future::plan("future::sequential", gc = TRUE)

  # # ..................................................................... ###

  # Plot species response curves

  IASDT.R::CatTime("Plot species response curves")

  PlotSpRC <- function(ID) {

    NFV <- Sp_DT_All$NFV[[ID]]
    Coords <- Sp_DT_All$Coords[[ID]]
    Species <- Sp_DT_All$Species[[ID]]
    DT <- Sp_DT_All$DT[[ID]]

    FilePrefix <- paste0("RespCurv_", Species, "_NFV_", NFV, "_Coords_", Coords)
    Path_JPEG <- file.path(Path_RespCurv_Sp, paste0(FilePrefix, ".jpeg"))

    DT <- DT %>%
      dplyr::slice(gtools::mixedorder(Variable)) %>%
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
            "</sub>(x + 0.1))</span>")))

    Species2 <- dplyr::filter(SpeciesNames, IAS_ID == !!Species)
    TitleTxt <- paste0(
      '<span style="font-size:13pt;"><b> Response curves of </b></span>',
      '<span style="color:blue; font-size:13pt;"><b><i>',
      Species2$Species_name, "</i></b></span>",
      '<span style="font-size:8pt;"> (<b>Class:</b> ', Species2$Class,
      "  &#8212; <b>Order:</b> ", Species2$Order,
      "  &#8212; <b>Family:</b> ", Species2$Family,
      "  &#8212; <b>ID:</b></span>",
      '<span style="font-size:8pt; color:blue;"> ',
      Species, '</span><span style="font-size:8pt;">)</span>')

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

    if (nrow(DT) <= 9) {
      NR <- NC <- 3
      PlotWidth <- 24
      PlotHeight <- 22
    } else {
      NR <- 3
      NC <- 4
      PlotWidth <- 30
      PlotHeight <- 22
    }

    Plots <- DT %>%
      dplyr::mutate(
        Plots = purrr::pmap(
          .l = list(PlotData_Quant, Observed_PA, PositiveTrendProb,
                    VarDesc, VarDesc2),
          .f = function(PlotData_Quant, Observed_PA, PositiveTrendProb,
                        VarDesc, VarDesc2) {

            Trend <- tibble::tibble(
              Trend2 = stringr::str_glue(
                "\n     Pr[pred(Var=max)] > Pr[pred(Var=min)] = ",
                "{round(PositiveTrendProb, 2)}"),
              X = -Inf, Y = Inf)

            Quant <- PlotData_Quant %>%
              tidyr::pivot_wider(
                id_cols = c(XVals), names_from = Quantile,
                values_from = Pred) %>%
              stats::setNames(c("XVals", "Q25", "Q50", "Q975"))

            ObsPA <- Observed_PA %>%
              dplyr::mutate(
                Col = dplyr::if_else(Pred == 1, "red", "darkgreen"),
                Pred = dplyr::case_when(
                  Pred == 1 ~ 0.97, Pred == 0 ~ 0.03, .default = Pred))

            Plot <- ggplot2::ggplot(
              data = ObsPA,
              mapping = ggplot2::aes(x = XVals, y = Pred, colour = Col)) +
              ggplot2::geom_jitter(
                shape = 16, width = 0, height = 0.02,
                alpha = 0.4, size = 1) +
              ggplot2::geom_line(
                ggplot2::aes(x = XVals, y = Q975), data = Quant,
                linetype = 2, linewidth = 0.3, colour = "blue") +
              ggplot2::geom_line(
                ggplot2::aes(x = XVals, y = Q25), data = Quant,
                linetype = 2, linewidth = 0.3, colour = "blue") +
              ggplot2::geom_ribbon(
                data = Quant, ggplot2::aes(x = XVals, ymin = Q25, ymax = Q975),
                inherit.aes = FALSE, fill = "blue", alpha = 0.1) +
              ggplot2::geom_line(
                mapping = ggplot2::aes(x = XVals, y = Q50), data = Quant,
                linetype = 1, linewidth = 0.6, colour = "blue") +
              ggplot2::geom_text(
                data = Trend, vjust = 1.4, hjust = -0.05,
                mapping = ggplot2::aes(x = X, y = Y, label = Trend2),
                colour = "grey30", size = 2.75) +
              ggplot2::scale_y_continuous(
                limits = c(-0.005, 1.075), oob = scales::squish,
                expand = c(0, 0)) +
              ggplot2::scale_x_continuous(expand = c(0.015, 0.015)) +
              ggplot2::labs(
                x = NULL, y = NULL, title = VarDesc, subtitle = VarDesc2) +
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

            return(Plot)

          })) %>%
      dplyr::pull(Plots)

    Plot <- patchwork::wrap_plots(Plots, nrow = NR, ncol = NC) +
      patchwork::plot_annotation(
        title = TitleTxt,
        caption = Caption,
        theme = ggplot2::theme(
          plot.margin = ggplot2::margin(0.1, 0, 0, 0, "cm"),
          plot.title = ggtext::element_markdown(),
          plot.caption = ggtext::element_markdown(
            size = 12, color = "grey", hjust = 0))) +
      patchwork::plot_layout(axes = "collect")

    Plot <- patchwork::wrap_elements(panel = Plot) +
      ggplot2::labs(tag = "<b>Predicted habitat suitability</b>") +
      ggplot2::theme(
        plot.tag = ggtext::element_markdown(
          size = ggplot2::rel(1), angle = 90),
        plot.tag.position = "left")


    ggplot2::ggsave(
      filename = Path_JPEG, plot = Plot,
      width = PlotWidth, height = PlotHeight,
      units = "cm", dpi = 600, quality = 100)

    Path_GG_Obj <- file.path(Path_GG, paste0(FilePrefix, ".RData"))
    IASDT.R::SaveAs(
      InObj = Plot, OutObj = FilePrefix, OutPath = Path_GG_Obj)
    OutDF <- tibble::tibble(
      Path_JPEG = Path_JPEG, Path_GG = Path_GG_Obj,
      PlotHeight = PlotHeight, PlotWidth = PlotWidth)

    return(OutDF)
  }

  # # ..................................................................... ###

  if (NCores == 1) {
    future::plan("future::sequential", gc = TRUE)
  } else {
    c1 <- snow::makeSOCKcluster(NCores)
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("future::cluster", workers = c1, gc = TRUE)
    on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
  }

  Plots <- future.apply::future_lapply(
    X = seq_len(nrow(Sp_DT_All)),
    FUN = PlotSpRC, future.seed = TRUE,
    future.chunk.size = 1,
    future.packages = c(
      "dplyr", "purrr", "tidyr", "gtools", "ggtext", "patchwork",
      "ggplot2", "tibble", "IASDT.R"),
    future.globals = c(
      "Sp_DT_All", "Path_RespCurv_Sp", "SpeciesNames",
      "Path_Postprocessing", "PlotSpRC", "Path_GG")) %>%
    dplyr::bind_rows()

  if (NCores > 1) {
    snow::stopCluster(c1)
    future::plan("future::sequential", gc = TRUE)
  }

  # # ..................................................................... ###

  # Save data
  IASDT.R::CatTime("Save data")

  Sp_DT_All <- Sp_DT_All %>%
    dplyr::select(-DT) %>%
    dplyr::bind_cols(Plots = Plots) %>%
    dplyr::select(-Path_GG)

  save(Sp_DT_All, file = file.path(Path_RespCurv_Sp, "Sp_DT_All.RData"))


  fs::dir_delete(Path_GG)

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "Plotting species response curves took ")

  # # ..................................................................... ###

  if (ReturnData) {
    return(Sp_DT_All)
  } else {
    return(invisible(NULL))
  }
}
