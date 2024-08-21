## |------------------------------------------------------------------------| #
# RespCurv_PlotSp ----
## |------------------------------------------------------------------------| #

#' Plot species response curves
#'
#' Generates and saves response curve plots for species based on environmental
#' variables and other factors. Plots show predicted habitat suitability across
#' different values of each variable.
#' @param Path_Model String. The path to the directory containing the models.
#'   The function reads data from the `Model_Postprocessing/RespCurv_SR`
#'   sub-directory, resulted from [RespCurv_PrepData].
#' @param NCores Integer. Number of cores to use for parallel processing.
#'   Defaults to 20.
#' @param EnvFile String. Path to the environment variables file. Defaults to
#'   ".env".
#' @param SaveGG Logical. Indicates whether to save the ggplot objects as RData
#'   files. Defaults to `FALSE`.
#' @param ShowProgress Logical. Indicates whether to show a progress bar during
#'   execution. Defaults to FALSE.
#' @param FromHPC Logical. Indicates whether the function is being run on an HPC
#'   environment, affecting file path handling. Default: `TRUE`.
#' @param ReturnData Logical. Indicates whether the output data be returned as
#'   an R object.
#' @return A data frame containing paths to the JPEG images and, if requested,
#'   the ggplot objects saved as RData files.
#' @export
#' @name RespCurv_PlotSp

RespCurv_PlotSp <- function(
    Path_Model = NULL, NCores = 20, EnvFile = ".env", SaveGG = FALSE,
    ShowProgress = FALSE, FromHPC = TRUE, ReturnData = FALSE) {

  if (is.null(Path_Model)) {
    stop("Path_Model cannot be NULL", .call = FALSE)
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  RC_Path_Orig <- RC_Path_SR <- RC_Path_Prob <- Coords <- NFV <- Species <-
    DT <- Variable <- VariableDesc <- Trend2 <- PositiveTrendProb <-
    VariableDesc <- Pred <- XVals <- PlotData_Quant <- Quantile <-
    Observed_PA <- Col <- Q975 <- Q25 <- Q50 <- X <- Y <- IAS_ID <- NULL

  AllArgs <- ls()
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character", Args = c("Path_Model", "EnvFile"))

  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "numeric", Args = "NCores")

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "logical", Args = c("SaveGG", "ShowProgress"))
  rm(AllArgs)

  fs::dir_create(file.path(Path_Model, "Model_Postprocessing", "RespCurv_Sp"))

  SpeciesNames <- IASDT.R::GetSpeciesName(EnvFile = EnvFile, FromHPC = FromHPC)

  c1 <- snow::makeSOCKcluster(NCores)
  on.exit(invisible(try(snow::stopCluster(c1), silent = TRUE)), add = TRUE)
  future::plan(future::cluster, workers = c1, gc = TRUE)

  SR_DT_All <- file.path(
    Path_Model, "Model_Postprocessing", "RespCurv_DT", "ResCurvDT.RData") %>%
    IASDT.R::LoadAs() %>%
    dplyr::select(-RC_Path_Orig, -RC_Path_SR) %>%
    dplyr::pull(RC_Path_Prob) %>%
    furrr::future_map_dfr(
      .f =  ~dplyr::select(IASDT.R::LoadAs(.x), -PlotData),
      .options = furrr::furrr_options(seed = TRUE, scheduling = Inf),
      .progress = ShowProgress) %>%
    tidyr::nest(
      DT = tidyselect::everything(), .by = c(Coords, NFV, Species)) %>%
    dplyr::mutate(
      Plot = furrr::future_pmap(
        .l = list(Coords, NFV, Species, DT),
        .f = function(Coords, NFV, Species, DT) {

          FilePrefix <- paste0(
            "RespCurv_", Species, "_Coords_", Coords, "_NFV_", NFV)
          Path_JPEG <- file.path(
            Path_Model, "Model_Postprocessing", "RespCurv_Sp",
            paste0(FilePrefix, ".jpeg"))

          DT <- DT %>%
            dplyr::mutate(
              VariableDesc = dplyr::case_when(
                Variable == "bio4" ~ "<span style='font-size: 10pt;'><b>Bio4</b></span><br/><span style='font-size: 8pt;'>temperature seasonality</span>",
                Variable == "bio6" ~ "<span style='font-size: 10pt;'><b>Bio6</b></span><br/><span style='font-size: 8pt;'>temperature of the coldest month</span>",
                Variable == "bio8" ~ "<span style='font-size: 10pt;'><b>Bio8</b></span><br/><span style='font-size: 8pt;'>temperatures of the wettest quarter</span>",
                Variable == "bio12" ~ "<span style='font-size: 10pt;'><b>Bio12</b></span><br/><span style='font-size: 8pt;'>annual precipitation amount</span>",
                Variable == "bio15" ~ "<span style='font-size: 10pt;'><b>Bio15</b></span><br/><span style='font-size: 8pt;'>precipitation seasonality</span>",
                Variable == "bio18" ~ "<span style='font-size: 10pt;'><b>Bio18</b></span><br/><span style='font-size: 8pt;'>monthly precipitation amount of the warmest quarter</span>",
                Variable == "RoadRailLog" ~ "<span style='font-size: 10pt;'><b>Road + Rail intensity</b></span><br/><span style='font-size: 8pt;'> (log<sub>10</sub>(x + 0.1))</span>",
                Variable == "BiasLog" ~ "<span style='font-size: 10pt;'><b>Sampling intensity</b></span><br/><span style='font-size: 8pt;'> (log<sub>10</sub>(x + 0.1))</span>",
                Variable == "HabLog" ~ "<span style='font-size: 10pt;'><b>% habitat coverage</b></span><br/><span style='font-size: 8pt;'> (log<sub>10</sub>(x + 0.1))</span>",

                .default = Variable),
              VariableDesc = factor(VariableDesc, levels = VariableDesc))

          Trend <- DT %>%
            dplyr::select(Variable, VariableDesc, Trend = PositiveTrendProb) %>%
            dplyr::mutate(
              Trend2 = stringr::str_glue(
                "\n     Pr[pred(Var=max)] > Pr[pred(Var=min)] = {round(Trend, 2)}"),
              X = -Inf, Y = Inf, Trend2 = as.character(Trend2))

          Quant <- DT %>%
            dplyr::select(Variable, VariableDesc, Quant = PlotData_Quant) %>%
            tidyr::unnest(cols = "Quant") %>%
            tidyr::pivot_wider(
              id_cols = c(Variable, VariableDesc, XVals), names_from = Quantile,
              values_from = Pred) %>%
            stats::setNames(
              c("Variable", "VariableDesc", "XVals", "Q25", "Q50", "Q975"))

          ObsPA <- DT %>%
            dplyr::select(Variable,  VariableDesc, PA = Observed_PA) %>%
            tidyr::unnest(cols = "PA") %>%
            dplyr::mutate(
              Col = dplyr::if_else(Pred == 1, "red", "darkgreen"),
              Pred = dplyr::case_when(
                Pred == 1 ~ 0.97, Pred == 0 ~ 0.03, .default = Pred))

          Species2 <- dplyr::filter(SpeciesNames, IAS_ID == !!Species)

          NFV_Label <- dplyr::if_else(
            NFV == 1,
            "non-focal variables are set to most likely value <i>[non.focalVariables = 1]</i>",
            "non-focal variables are set to most likely value given the value of focal variable <i>[non.focalVariables = 2]</i>")

          Coords_Label <- dplyr::if_else(
            Coords == "c",
            'Predictions at mean coordinates <i>[coordinates = "c"]</i>',
            " Predictions without effect of spatial ",
            'dependence <i>[coordinates = "i"]</i>')

          TitleTxt <- paste0(
            '<span style="font-size:13pt;">',
            "<b> Response curves of </b></span>",
            '<span style="color:blue; font-size:13pt;"><b><i>',
            Species2$Species_name, "</i></b></span>",
            '<span style="font-size:8pt;"> (<b>Class:</b> ', Species2$Class,
            "  &#8212; <b>Order:</b> ", Species2$Order,
            "  &#8212; <b>Family:</b> ", Species2$Family,
            '  &#8212; <b>ID:</b></span><span style="font-size:8pt; color:blue;"> ',
            Species,  '</span><span style="font-size:8pt;">)</span>')
          SubTitleTxt <- paste0(
            Coords_Label, "  &#8212;&#8212;  ",
            NFV_Label, "</span>")

          Plot <- ggplot2::ggplot(
            data = ObsPA,
            mapping = ggplot2::aes(x = XVals, y = Pred, colour = Col)) +
            ggplot2::geom_jitter(
              shape = 16, width = 0, height = 0.02, alpha = 0.4, size = 0.8) +
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
              data = Trend, mapping = ggplot2::aes(x = X, y = Y, label = Trend2),
              colour = "grey30", size = 2.75, vjust = 1.4, hjust = -0.05) +
            ggplot2::scale_y_continuous(
              limits = c(-0.005, 1.05), oob = scales::squish, expand = c(0, 0)) +
            ggplot2::scale_x_continuous(expand = c(0.015, 0.015)) +
            ggplot2::facet_wrap(
              ggplot2::vars(VariableDesc),
              scales = "free_x", nrow = 3, ncol = 3) +
            ggplot2::xlab("Predictor value") +
            ggplot2::ylab("Predicted habitat suitability") +
            ggplot2::labs(title = TitleTxt, subtitle = SubTitleTxt) +
            ggplot2::theme_bw() +
            ggplot2::theme(
              strip.text = ggtext::element_markdown(
                hjust = 0,
                margin = ggplot2::margin(0.05, 0.1, 0.05, 0.1, "cm")),
              strip.background = ggplot2::element_rect(
                colour = NA, fill = "white"),
              legend.position = "none",
              axis.title = ggtext::element_markdown(size = 12, face = "bold"),
              axis.text = ggplot2::element_text(size = 8),
              plot.title = ggtext::element_markdown(
                margin = ggplot2::margin(1, 0, 1, 0), hjust = 0.1),
              plot.title.position = "plot",
              plot.subtitle = ggtext::element_markdown(
                size = 8, colour = "darkgrey", hjust = 0.1,
                margin = ggplot2::margin(4, 0, 4, 0)),
              panel.grid.major = ggplot2::element_line(linewidth = 0.25),
              panel.grid.minor = ggplot2::element_line(linewidth = 0.1),
              panel.spacing = ggplot2::unit(0.15, "lines"),
              plot.margin = ggplot2::unit(c(0.1, 0.2, 0.1, 0.2), "lines"))

          ggplot2::ggsave(
            filename = Path_JPEG, plot = Plot, width = 24, height = 22,
            dpi = 600, units = "cm")

          if (SaveGG) {
            PathGG <- file.path(
              Path_Model, "Model_Postprocessing", "RespCurv_Sp_GG")
            fs::dir_create(PathGG)
            Path_GG <- file.path(PathGG, paste0(FilePrefix, ".RData"))
            IASDT.R::SaveAs(
              InObj = Plot, OutObj = FilePrefix, OutPath = Path_GG)
            OutDF <- tibble::tibble(Path_JPEG = Path_JPEG, Path_GG = Path_GG)
          } else {
            OutDF <- tibble::tibble(Path_JPEG = Path_JPEG)
          }

          return(OutDF)
        },
        .options = furrr::furrr_options(seed = TRUE, scheduling = Inf),
        .progress = ShowProgress)) %>%
    dplyr::select(-DT) %>%
    tidyr::unnest(cols = "Plot")

  if (ReturnData) {
    return(SR_DT_All)
  } else {
    return(invisible(NULL))
  }
}
