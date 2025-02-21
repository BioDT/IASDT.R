## |------------------------------------------------------------------------| #
# RespCurv_PlotSR ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname Response_curves
#' @name Response_curves
#' @order 4
#' @author Ahmed El-Gabbas

RespCurv_PlotSR <- function(ModelDir, Verbose = TRUE, NCores = 8L) {

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
    Pred <- Q975 <- Q25 <- Q50 <- X <- Y <- Variable2 <- Var2 <-
    Variable1 <- NULL

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

  Path_RC_DT <- IASDT.R::Path(ModelDir, "Model_Postprocessing", "RespCurv_DT")
  Path_RC_SR <- IASDT.R::Path(ModelDir, "Model_Postprocessing", "RespCurv_SR")

  if (!dir.exists(Path_RC_DT)) {
    stop("Response curve data subfolder is missing.", call. = FALSE)
  }

  fs::dir_create(Path_RC_SR)

  # # ..................................................................... ###

  IASDT.R::CatTime("Create species richness response curves")

  SR_DT_All <- IASDT.R::Path(Path_RC_DT, "ResCurvDT.RData") %>%
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


  VarLabel <- tibble::tribble(
    ~Variable1, ~Var2, ~Variable2,
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
    "npp", "NPP", "net primary productivity",
    "RiversLog", "River length", "log10(x + 0.1)",
    "RoadRailLog", "Road + Rail intensity", "log10(x + 0.1)",
    "HabLog", "% habitat coverage", "log10(x + 0.1)",
    "EffortsLog", "Sampling efforts", "log10(x + 0.1)") %>%
    dplyr::mutate(
      Variable2 = paste0(
        "<span style='font-size: 10pt;'><b>", Var2, "</b></span>",
        "<span style='font-size: 8pt;'> (", Variable2, ")</span>")) %>%
    dplyr::select(-Var2)

  SR_DT_All <- SR_DT_All %>%
    dplyr::mutate(
      Plot = purrr::pmap(
        .l = list(Variable, Quant, Observed, Trend, Coords),
        .f = function(Variable, Quant, Observed, Trend, Coords) {

          IASDT.R::CatTime(
            paste0(Variable, " - coords = ", Coords), Level = 2)

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
          Variable2 <- dplyr::filter(VarLabel, Variable1 == Variable) %>%
            dplyr::pull(Variable2)

          VarName <- dplyr::case_when(
            Variable == "HabLog" ~ "% Habitat coverage",
            Variable == "RoadRailLog" ~ "Road + Rail intensity",
            Variable == "EffortsLog" ~ "Sampling efforts",
            Variable == "RiversLog" ~ "River length",
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

          # Plotting theme
          Theme <- ggplot2::theme(
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

          # Plot
          Plot <- ggplot2::ggplot(
            data = Observed, mapping = ggplot2::aes(x = XVals, y = Pred)) +
            ggplot2::geom_jitter(
              shape = 16, width = 0, height = 0.02, alpha = 0.2, size = 0.75,
              colour = "darkgrey") +
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
            Theme

          # Using ggplot2::ggsave directly does not show non-ascii characters
          # correctly
          ragg::agg_jpeg(
            filename = IASDT.R::Path(
              Path_RC_SR,
              paste0("RespCurv_SR_", Variable, "_Coords_", Coords, ".jpeg")),
            width = 20, height = 12.5, res = 600, quality = 100, units = "cm")
          print(Plot)
          grDevices::dev.off()


          # Back-transforming variables
          if (Variable %in% c(
            "EffortsLog", "RoadRailLog", "HabLog", "RiversLog")) {

            IASDT.R::CatTime(
              paste0(Variable, " - coords = ", Coords, " - original scale"),
              Level = 2)

            Observed2 <- dplyr::mutate(Observed, XVals = 10 ^ XVals - 0.1)
            Quant2 <- dplyr::mutate(Quant, XVals = 10 ^ XVals - 0.1)

            # Maximum value on the y-axis
            PlotMax <- max(Observed2$Pred, Quant2$Q975) * 1.05

            # Variable long name (x-axis label)
            Variable2 <- dplyr::case_when(
              Variable == "RoadRailLog" ~ "Road + Rail intensity",
              Variable == "HabLog" ~ "% habitat coverage",
              Variable == "RiversLog" ~ "River length",
              Variable == "EffortsLog" ~ "Sampling efforts",
              .default = Variable) %>%
              paste0(
                "<span style='font-size: 10pt;'><b>", .,
                "</b></span><span style='font-size: 7pt;'> (original scale)",
                "</span>")

            Plot2 <- ggplot2::ggplot(
              data = Observed2, mapping = ggplot2::aes(x = XVals, y = Pred)) +
              ggplot2::geom_jitter(
                shape = 16, width = 0, height = 0.02, alpha = 0.2, size = 0.75,
                colour = "darkgrey") +
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
              Theme

            # Using ggplot2::ggsave directly does not show non-ascii characters
            # correctly
            ragg::agg_jpeg(
              filename = IASDT.R::Path(
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
    file = IASDT.R::Path(Path_RC_SR, "SR_DT_All.RData"))

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTime,
    Prefix = "Plotting response curves for species richness took ")

  # # ..................................................................... ###

  return(invisible(NULL))
}
