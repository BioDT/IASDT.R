
## |------------------------------------------------------------------------| #
# RespCurv_PlotSp ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname Response_curves
#' @name Response_curves
#' @order 2
#' @author Ahmed El-Gabbas

RespCurv_PlotSp <- function(
    ModelDir = NULL, NCores = 20, EnvFile = ".env",
    FromHPC = TRUE, ReturnData = FALSE) {

  # # ..................................................................... ###

  IASDT.R::CatTime("Plotting species response curves")

  .StartTime <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Path_PA <- NCells_Naturalized <- NFV <- Coords <- Species <- Prefix <-
    Data <- RC_Path_Prob <- Variable <- IAS_ID <- VarDesc <- VarDesc2 <- NULL

  # # ..................................................................... ###

  # Check arguments

  IASDT.R::CatTime("Check arguments")

  if (is.null(ModelDir)) {
    stop("`ModelDir` cannot be NULL", call. = FALSE)
  }

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character", Args = c("ModelDir", "EnvFile"))

  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "numeric", Args = "NCores")
  rm(AllArgs, envir = environment())

  # # ..................................................................... ###

  Path_RC_DT <- IASDT.R::Path(ModelDir, "Model_Postprocessing", "RespCurv_DT")
  if (!dir.exists(Path_RC_DT)) {
    stop("Response curve data subfolder is missing.", call. = FALSE)
  }
  Path_RC_Sp <- IASDT.R::Path(ModelDir, "Model_Postprocessing", "RespCurv_Sp")
  Path_RC_Sp_DT <- IASDT.R::Path(
    ModelDir, "Model_Postprocessing", "RespCurv_Sp_DT")
  fs::dir_create(c(Path_RC_Sp, Path_RC_Sp_DT))

  # # ..................................................................... ###

  # # Load species summary
  IASDT.R::CatTime("Load species summary")

  if (FromHPC) {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_PA", "DP_R_PA", TRUE, FALSE)
  } else {
    EnvVars2Read <- tibble::tribble(
      ~VarName, ~Value, ~CheckDir, ~CheckFile,
      "Path_PA", "DP_R_PA_Local", TRUE, FALSE)
  }

  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)

  SpSummary <- IASDT.R::Path(Path_PA, "Sp_PA_Summary_DF.csv")
  if (!file.exists(SpSummary)) {
    stop(paste0(SpSummary, " file does not exist"), call. = FALSE)
  }

  SpSummary <- readr::read_csv(SpSummary, show_col_types = FALSE) %>%
    dplyr::select(tidyselect::all_of(c("IAS_ID", "NCells_Naturalized"))) %>%
    dplyr::rename(NCells = NCells_Naturalized)

  # # ..................................................................... ###

  # Load species names
  IASDT.R::CatTime("Load species names")
  SpeciesNames <- IASDT.R::GetSpeciesName(EnvFile = EnvFile, FromHPC = FromHPC)

  # # ..................................................................... ###

  # Prepare species-specific data on parallel

  IASDT.R::CatTime("Prepare species-specific data on parallel")

  IASDT.R::CatTime("Prepare working on parallel", Level = 1)
  c1 <- parallel::makePSOCKcluster(NCores)
  on.exit(try(parallel::stopCluster(c1), silent = TRUE), add = TRUE)

  IASDT.R::CatTime("Exporting objects to cores", Level = 1)
  parallel::clusterExport(
    cl = c1, varlist = "Path_RC_DT", envir = environment())

  IASDT.R::CatTime("Load packages at each core", Level = 1)
  invisible(parallel::clusterEvalQ(
    cl = c1, expr = sapply("IASDT.R", library, character.only = TRUE)))

  IASDT.R::CatTime("Loading species richness data", Level = 1)
  Sp_DT_All <- IASDT.R::Path(Path_RC_DT, "ResCurvDT.RData") %>%
    IASDT.R::LoadAs() %>%
    dplyr::select(tidyselect::all_of(c("Coords", "RC_Path_Prob"))) %>%
    dplyr::mutate(
      Data = parallel::parLapplyLB(
        cl = c1, X = RC_Path_Prob, fun = IASDT.R::LoadAs)) %>%
    tidyr::unnest(Data) %>%
    dplyr::select(-RC_Path_Prob) %>%
    tidyr::nest(
      DT = tidyselect::everything(), .by = c(NFV, Coords, Species)) %>%
    dplyr::mutate(IAS_ID = as.numeric(stringr::str_remove(Species, "^Sp_"))) %>%
    dplyr::left_join(SpSummary, by = "IAS_ID") %>%
    dplyr::mutate(
      Prefix = paste0(Species, "_NFV_", NFV, "_Coords_", Coords),
      Path_JPEG_Fixed = IASDT.R::Path(
        Path_RC_Sp, paste0(Prefix, "_Fixed.jpeg")),
      Path_JPEG_Free = IASDT.R::Path(Path_RC_Sp, paste0(Prefix, "_Free.jpeg")),
      Path_Sp_DT = IASDT.R::Path(Path_RC_Sp_DT, paste0(Prefix, ".qs2")))

  snow::stopCluster(c1)

  IASDT.R::CatTime("Export species-specific data", Level = 1)
  purrr::walk(
    .x = seq_len(nrow(Sp_DT_All)),
    .f = function(ID) {
      DT <- dplyr::slice(Sp_DT_All, ID)
      if (isFALSE(IASDT.R::CheckData(DT$Path_Sp_DT, warning = FALSE))) {
        IASDT.R::SaveAs(InObj = DT, OutPath = DT$Path_Sp_DT)
      }
    })

  Sp_DT_All <- gtools::mixedsort(Sp_DT_All$Path_Sp_DT)
  invisible(gc())

  # # ..................................................................... ###

  IASDT.R::CatTime("Plotting species-specific data")

  IASDT.R::CatTime("Prepare working on parallel", Level = 1)
  c1 <- parallel::makePSOCKcluster(NCores)
  on.exit(try(parallel::stopCluster(c1), silent = TRUE), add = TRUE)

  IASDT.R::CatTime("Exporting objects to cores", Level = 1)
  parallel::clusterExport(
    cl = c1, varlist = c("SpeciesNames", "Sp_DT_All"), envir = environment())

  IASDT.R::CatTime("Load packages at each core", Level = 1)
  invisible(parallel::clusterEvalQ(
    cl = c1,
    expr = {
      sapply(
        c(
          "dplyr", "purrr", "tidyr", "gtools", "ggtext", "patchwork",
          "ggplot2", "tibble", "IASDT.R", "ragg"),
        library, character.only = TRUE)
    }))

  IASDT.R::CatTime("Plotting on parallel", Level = 1)
  Plots <- parallel::clusterApplyLB(
    cl = c1,
    x = Sp_DT_All,
    fun = function(RC_File) {

      DT <- IASDT.R::LoadAs(RC_File)

      Coords <- DT$Coords
      Species <- DT$Species
      NCells <- DT$NCells
      NFV <- DT$NFV
      Path_JPEG_Fixed <- DT$Path_JPEG_Fixed
      Path_JPEG_Free <- DT$Path_JPEG_Free

      DT <- dplyr::select(DT, DT) %>%
        tidyr::unnest(DT) %>%
        dplyr::slice(gtools::mixedorder(Variable)) %>%
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
        PlotWidth <- 24
        PlotHeight <- 22
      } else {
        NR <- 3
        NC <- 4
        PlotWidth <- 30
        PlotHeight <- 22
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
        filename = Path_JPEG_Fixed, width = PlotWidth, height = PlotHeight,
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
        filename = Path_JPEG_Free, width = PlotWidth, height = PlotHeight,
        res = 600, quality = 100, units = "cm")
      print(Plot_Free)
      grDevices::dev.off()

      OutDF <- tibble::tibble(
        Coords = Coords, Species = Species, NCells = NCells, NFV = NFV,
        Path_JPEG_Fixed = Path_JPEG_Fixed, Path_JPEG_Free = Path_JPEG_Free,
        PlotHeight = PlotHeight, PlotWidth = PlotWidth)

      return(OutDF)
    }) %>%
    dplyr::bind_rows()

  snow::stopCluster(c1)
  invisible(gc())

  # # ..................................................................... ###

  # Save data
  IASDT.R::CatTime("Save data")
  save(Plots, file = IASDT.R::Path(Path_RC_Sp, "Sp_DT_All.RData"))

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
