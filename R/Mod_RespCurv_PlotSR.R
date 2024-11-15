## |------------------------------------------------------------------------| #
# RespCurv_PrepData ----
## |------------------------------------------------------------------------| #

#' Plot species richness response curves
#'
#' Generate species richness response curves as JPEG images for each variable
#' used in the model.
#'
#' @param ModelDir String. Path to the root directory of the fitted models
#'   without the trailing slash. The function reads data from `RespCurv_DT`
#'   subdirectory created by [RespCurv_PrepData].
#' @param NCores Integer. Number of cores to use for parallel processing.
#'   Defaults to 8.
#' @return This function does not return a value but saves JPEG images of the
#'   response curves in a subdirectory within the specified path.
#' @name RespCurv_PlotSR
#' @inheritParams RespCurv_PrepData
#' @author Ahmed El-Gabbas
#' @seealso RespCurv_PlotSp RespCurv_PrepData
#' @export

RespCurv_PlotSR <- function(ModelDir, Verbose = TRUE, NCores = 8) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  if (isFALSE(Verbose)) {
    sink(file = nullfile())
    on.exit(try(sink(), silent = TRUE), add = TRUE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Trend2 <- Variable <- Quant <- Observed <- Trend <- NFV <- Coords <-
    RC_Path_SR <- RC_Path_Orig <- RC_Path_Prob <- DT <- data <- XVals <-
    Pred <- Q975 <- Q25 <- Q50 <- X <- Y <- NULL

  # # ..................................................................... ###

  IASDT.R::CatTime("Check input arguments")
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character", Args = "ModelDir")
  rm(AllArgs, envir = environment())

  # # ..................................................................... ###

  IASDT.R::CatTime("Check the existence of response curve directory")

  Path_RC_DT <- file.path(ModelDir, "Model_Postprocessing", "RespCurv_DT")
  Path_RC_SR <- file.path(ModelDir, "Model_Postprocessing", "RespCurv_SR")

  if (!dir.exists(Path_RC_DT)) {
    stop("Response curve data subfolder is missing.", call. = FALSE)
  }

  fs::dir_create(Path_RC_SR)

  # # ..................................................................... ###

  IASDT.R::CatTime("Create species richness response curves")

  SR_DT_All <- file.path(Path_RC_DT, "ResCurvDT.RData") %>%
    IASDT.R::LoadAs() %>%
    dplyr::select(-RC_Path_Orig, -RC_Path_Prob)

  NCores <- max(min(NCores, nrow(SR_DT_All)), 1)

  IASDT.R::CatTime(
    paste0("Prepare working on parallel, using ", NCores, " cores"),
    Level = 1)

  withr::local_options(
    future.globals.maxSize = 8000 * 1024^2,
    future.gc = TRUE, future.seed = TRUE)

  c1 <- snow::makeSOCKcluster(NCores)
  on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
  future::plan("future::cluster", workers = c1, gc = TRUE)
  on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)

  IASDT.R::CatTime("Prepare data", Level = 1)

  SR_DT_All <- SR_DT_All %>%
    dplyr::mutate(
      DT = furrr::future_map(
        .x = RC_Path_SR,
        .f = ~ {
          DT <- IASDT.R::LoadAs(.x) %>%
            magrittr::inset(c("RC_Data_SR"), NULL)

          Quant <- DT$RC_Data_SR_Quant %>%
            dplyr::mutate(
              Variable = DT$Variable, NFV = DT$NFV, .before = 1) %>%
            tidyr::pivot_wider(
              id_cols = c("Variable", "NFV", "XVals"),
              names_from = Quantile, values_from = SR) %>%
            setNames(c("Variable", "NFV", "XVals", "Q25", "Q50", "Q975"))

          Observed <- DT$Observed_SR %>%
            dplyr::mutate(Variable = DT$Variable, NFV = DT$NFV, .before = 1)

          Trend <- tibble::tibble(
            Variable = DT$Variable, NFV = DT$NFV,
            Trend = DT$SR_PositiveTrendProb)

          return(list(Quant = Quant, Observed = Observed, Trend = list(Trend)))
        },
        .options = furrr::furrr_options(seed = TRUE, chunk_size = 1))) %>%
    dplyr::select(-NFV, -RC_Path_SR) %>%
    tidyr::unnest_wider(DT) %>%
    tidyr::nest(.by = c(Variable, Coords)) %>%
    dplyr::mutate(
      Quant = purrr::map(.x = data, .f = ~ dplyr::bind_rows(.x$Quant)),
      Observed = purrr::map(.x = data, .f = ~ dplyr::bind_rows(.x$Observed)),
      Trend = purrr::map(.x = data, .f = ~ dplyr::bind_rows(.x$Trend))) %>%
    dplyr::select(-data)

  snow::stopCluster(c1)
  future::plan("future::sequential", gc = TRUE)
  invisible(gc())

  # # ..................................................................... ###

  # Plot species richness response curves

  IASDT.R::CatTime("Plot species richness response curves", Level = 1)

  SR_DT_All <- SR_DT_All %>%
    dplyr::mutate(
      Plot = purrr::pmap(
        .l = list(Variable, Quant, Observed, Trend, Coords),
        .f = function(Variable, Quant, Observed, Trend, Coords) {

          # Maximum value on the y-axis
          PlotMax <- max(Observed$Pred, Quant$Q975) * 1.05

          # Label on the y axis
          YAxisLabel <- paste0(
            "<span style='font-size: 12pt;'><b>Predicted ",
            "species richness</span></b>")

          # Trend text
          DT_Trend <- Trend %>%
            dplyr::mutate(
              Trend2 = stringr::str_glue(
                "\n     Pr[pred(Var=max)] > ",
                "Pr[pred(Var=min)] = {round(Trend, 2)}"),
              Trend2 = as.character(Trend2), X = -Inf, Y = Inf)

          # Variable long name (x-axis label)
          Variable2 <- dplyr::case_when(
            Variable == "bio2" ~ paste0(
              "<span style='font-size: 10pt;'><b>Bio2</b></span><span ",
              "style='font-size: 7pt;'> (Mean diurnal range [mean of monthly ",
              "(max temperature - min temperature)])</span>"),
            Variable == "bio4" ~ paste0(
              "<span style='font-size: 10pt;'><b>Bio4</b></span><span ",
              "style='font-size: 7pt;'> (temperature seasonality [standard ",
              "deviation &times;100])</span>"),
            Variable == "bio6" ~ paste0(
              "<span style='font-size: 10pt;'><b>Bio6</b></span><span ",
              "style='font-size: 7pt;'> (temperature of the coldest ",
              "month)</span>"),
            Variable == "bio8" ~ paste0(
              "<span style='font-size: 10pt;'><b>Bio8</b></span><span ",
              "style='font-size: 7pt;'> (temperatures of the wettest ",
              "quarter)</span>"),
            Variable == "bio12" ~ paste0(
              "<span style='font-size: 10pt;'><b>Bio12</b></span><span ",
              "style='font-size: 7pt;'> (annual precipitation amount)</span>"),
            Variable == "bio15" ~ paste0(
              "<span style='font-size: 10pt;'><b>Bio15</b></span><span ",
              "style='font-size: 7pt;'> (precipitation seasonality ",
              "[Coefficient of Variation])</span>"),
            Variable == "bio18" ~ paste0(
              "<span style='font-size: 10pt;'><b>Bio18</b></span><span ",
              "style='font-size: 7pt;'> (monthly precipitation amount of ",
              "the warmest quarter)</span>"),
            Variable == "RoadRailLog" ~ paste0(
              "<span style='font-size: 10pt;'><b>Road + Rail intensity</b>",
              "</span><span style='font-size: 7pt;'> (log<sub>10</sub>(x + ",
              "0.1))</span>"),
            Variable == "HabLog" ~ paste0(
              "<span style='font-size: 10pt;'><b>% habitat coverage</b></span>",
              "<span style='font-size: 8pt;'> (log<sub>10",
              "</sub>(x + 0.1))</span>"),
            Variable == "EffortsLog" ~ paste0(
              "<span style='font-size: 10pt;'><b>Sampling efforts</b>",
              "</span><span style='font-size: 7pt;'> (log<sub>10</sub>(x + ",
              "0.1))</span>"),
            .default = Variable)

          VarName <- dplyr::case_when(
            Variable == "HabLog" ~ "% Habitat coverage",
            Variable == "RoadRailLog" ~ "Road + Rail intensity",
            Variable == "EffortsLog" ~ "Sampling efforts",
            .default = Variable)

          # facetting labels
          Label1 <- paste0(
            '<span style="font-size:8pt; color:red;"><b>non.focalVariables',
            '= 1</b></span><br><span style="font-size:6pt;">values of ',
            "non-focal variables are set to the most likely value</span>")
          Label2 <- paste0(
            '<span style="font-size:8pt; color:red;"><b>non.focalVariables',
            '= 2</b></span><br><span style="font-size:6pt;">values of ',
            "non-focal variables are set to the most likely value given the",
            " value of focal variable</span>")
          FacetLabel <- ggplot2::as_labeller(c(`1` = Label1, `2` = Label2))

          # Plot title
          TitleTxt <- paste0(
            '<span style="font-size:12pt; color:blue;">',
            "Predicted species richness for <b>", VarName, "</b></span>")

          Caption <- dplyr::if_else(
            Coords == "c",
            "Predictions are made at mean coordinates",
            paste0(
              "Predictions are made for infinite coordinates ",
              "without effect of spatial dependence"))

          # Plot
          Plot <- ggplot2::ggplot(
            data = Observed, mapping = ggplot2::aes(x = XVals, y = Pred)) +
            ggplot2::geom_jitter(
              shape = 16, width = 0, height = 0.02, alpha = 0.2, size = 1,
              colour = "blue") +
            ggplot2::geom_line(
              ggplot2::aes(x = XVals, y = Q975),
              data = Quant, linetype = 2, linewidth = 0.3, colour = "blue") +
            ggplot2::geom_line(
              ggplot2::aes(x = XVals, y = Q25),
              data = Quant, linetype = 2, linewidth = 0.3, colour = "blue") +
            ggplot2::geom_ribbon(
              data = Quant, ggplot2::aes(x = XVals, ymin = Q25, ymax = Q975),
              inherit.aes = FALSE, fill = "blue", alpha = 0.1) +
            ggplot2::geom_line(
              mapping = ggplot2::aes(x = XVals, y = Q50), data = Quant,
              linetype = 1, linewidth = 0.6, colour = "blue") +
            ggplot2::geom_text(
              data = DT_Trend,
              mapping = ggplot2::aes(x = X, y = Y, label = Trend2),
              colour = "darkgrey", size = 2.5, vjust = 1.4, hjust = -0.05) +
            ggplot2::facet_grid(~NFV, labeller = FacetLabel) +
            ggplot2::scale_y_continuous(
              limits = c(-1, PlotMax), oob = scales::squish, expand = c(0, 0)) +
            ggplot2::scale_x_continuous(expand = c(0.015, 0.015)) +
            ggplot2::xlab(Variable2) +
            ggplot2::ylab(YAxisLabel) +
            ggplot2::labs(title = TitleTxt, caption = Caption) +
            ggplot2::theme_bw() +
            ggplot2::theme(
              plot.margin = ggplot2::margin(5, 20, 5, 5),
              strip.text = ggtext::element_markdown(
                hjust = 0,
                margin = ggplot2::margin(0.05, 0.1, 0.05, 0.1, "cm")),
              strip.background = ggplot2::element_rect(
                colour = NA, fill = "white"),
              legend.position = "none",
              plot.caption = ggplot2::element_text(
                size = 8, color = "grey", hjust = 0),
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

          # Using ggplot2::ggsave directly does not show non-ascii characters
          # correctly
          ragg::agg_jpeg(
            filename = file.path(
              Path_RC_SR,
              paste0("RespCurv_SR_", Variable, "_Coords_", Coords, ".jpeg")),
            width = 20, height = 12.5, res = 600, quality = 100, units = "cm")
          print(Plot)
          grDevices::dev.off()


          # Back-transforming variables
          if (Variable %in% c("EffortsLog", "RoadRailLog", "HabLog")) {

            Observed2 <- dplyr::mutate(Observed, XVals = 10 ^ XVals - 0.1)
            Quant2 <- dplyr::mutate(Quant, XVals = 10 ^ XVals - 0.1)

            # Maximum value on the y-axis
            PlotMax <- max(Observed2$Pred, Quant2$Q975) * 1.05

            # Variable long name (x-axis label)
            Variable2 <- dplyr::case_when(
              Variable == "RoadRailLog" ~ paste0(
                "<span style='font-size: 10pt;'><b>Road + Rail intensity</b>",
                "</span><span style='font-size: 7pt;'> (back-transformed)",
                "</span>"),
              Variable == "HabLog" ~ paste0(
                "<span style='font-size: 10pt;'>",
                "<b>% habitat coverage</b></span>",
                "<span style='font-size: 8pt;'> (back-transformed)</span>"),
              Variable == "EffortsLog" ~ paste0(
                "<span style='font-size: 10pt;'><b>Sampling efforts</b>",
                "</span><span style='font-size: 7pt;'> (back-transformed)",
                "</span>"),
              .default = Variable)

            Plot2 <- ggplot2::ggplot(
              data = Observed2, mapping = ggplot2::aes(x = XVals, y = Pred)) +
              ggplot2::geom_jitter(
                shape = 16, width = 0, height = 0.02, alpha = 0.2, size = 1,
                colour = "blue") +
              ggplot2::geom_line(
                ggplot2::aes(x = XVals, y = Q975),
                data = Quant2, linetype = 2, linewidth = 0.3, colour = "blue") +
              ggplot2::geom_line(
                ggplot2::aes(x = XVals, y = Q25),
                data = Quant2, linetype = 2, linewidth = 0.3, colour = "blue") +
              ggplot2::geom_ribbon(
                data = Quant2, ggplot2::aes(x = XVals, ymin = Q25, ymax = Q975),
                inherit.aes = FALSE, fill = "blue", alpha = 0.1) +
              ggplot2::geom_line(
                mapping = ggplot2::aes(x = XVals, y = Q50), data = Quant2,
                linetype = 1, linewidth = 0.6, colour = "blue") +
              ggplot2::geom_text(
                data = DT_Trend,
                mapping = ggplot2::aes(x = X, y = Y, label = Trend2),
                colour = "darkgrey", size = 2.5, vjust = 1.4, hjust = -0.05) +
              ggplot2::facet_grid(~NFV, labeller = FacetLabel) +
              ggplot2::scale_y_continuous(
                limits = c(-1, PlotMax),
                oob = scales::squish, expand = c(0, 0)) +
              ggplot2::scale_x_continuous(
                expand = c(0.015, 0.015), labels = scales::comma) +
              ggplot2::xlab(Variable2) +
              ggplot2::ylab(YAxisLabel) +
              ggplot2::labs(title = TitleTxt, caption = Caption) +
              ggplot2::theme_bw() +
              ggplot2::theme(
                strip.text = ggtext::element_markdown(
                  hjust = 0,
                  margin = ggplot2::margin(0.05, 0.1, 0.05, 0.1, "cm")),
                strip.background = ggplot2::element_rect(
                  colour = NA, fill = "white"),
                legend.position = "none",
                plot.caption = ggplot2::element_text(
                  size = 8, color = "grey", hjust = 0),
                axis.title = ggtext::element_markdown(),
                axis.text = ggplot2::element_text(size = 8),
                plot.title = ggtext::element_markdown(
                  margin = ggplot2::margin(1, 0, 1, 0)),
                plot.subtitle = ggtext::element_markdown(
                  size = 7, colour = "darkgrey",
                  margin = ggplot2::margin(4, 0, 4, 0)),
                panel.grid.major = ggplot2::element_line(linewidth = 0.25),
                panel.grid.minor = ggplot2::element_line(linewidth = 0.1),
                panel.spacing = ggplot2::unit(0.15, "lines"),
                plot.margin = ggplot2::unit(c(0.1, 0.2, 0.1, 0.2), "lines"))


            # Using ggplot2::ggsave directly does not show non-ascii characters
            # correctly
            ragg::agg_jpeg(
              filename = file.path(
                Path_RC_SR,
                paste0(
                  "RespCurv_SR_", Variable, "_Coords_", Coords,
                  "_OriginalScale.jpeg")),
              width = 20, height = 12.5, res = 600, quality = 100, units = "cm")
            print(Plot2)
            grDevices::dev.off()
          }

          return(invisible(NULL))
        }))

  save(
    SR_DT_All,
    file = file.path(Path_RC_SR, "SR_DT_All.RData"))

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTime,
    Prefix = "Plotting response curves for species richness took ")

  # # ..................................................................... ###

  return(invisible(NULL))
}
