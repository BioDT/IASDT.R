## |------------------------------------------------------------------------| #
# RespCurv_PrepData ----
## |------------------------------------------------------------------------| #

#' Plot species richness response curves
#'
#' This function generates and saves species richness response curves as JPEG images for each variable in the dataset. It processes data from a specified directory, creates plots, and saves them in a subdirectory.
#'
#' @param Path_RC String. The path to the directory containing the response curve data.
#' @return This function does not return a value but saves JPEG images of the response curves in a subdirectory within the specified path.
#' @name RespCurv_PlotSR
#' @author Ahmed El-Gabbas
#' @seealso RespCurv_PlotSp RespCurv_PrepData
#' @export


RespCurv_PlotSR <- function(Path_RC) {

  if (is.null(Path_RC)) {
    stop("Path_RC cannot be NULL")
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Trend2 <- Variable <- Quant <- Observed <- Trend <- Coords <- NFV <-
    RC_Path_SR <- RC_Path_Orig <- RC_Path_Prob <- DT <- data <- XVals <-
    Pred <- Q975 <- Q25 <- Q50 <- X <- Y <- NULL

  IASDT.R::CatTime("Check input arguments")
  AllArgs <- ls()
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "character", Args = "Path_RC")
  rm(AllArgs)

  IASDT.R::CatTime("Check the existence of response curve directory")
  DirsMissing <- c(Path_RC, file.path(Path_RC, "RespCurv_DT")) %>%
    dir.exists() %>%
    all() %>%
    magrittr::not()

  if (DirsMissing) {
    stop(
      "Response curve directory or data subfolder does not exist", call. = FALSE)
  }
  rm(DirsMissing)

  fs::dir_create(file.path(Path_RC, "RespCurv_SR"))

  IASDT.R::CatTime("Sequentially create species richness response curves")
  SR_DT_All <- file.path(Path_RC, "RespCurv_DT/ResCurvDT.RData") %>%
    IASDT.R::LoadAs() %>%
    dplyr::select(-RC_Path_Orig, -RC_Path_Prob) %>%
    dplyr::mutate(
      DT = purrr::map(
        .x = RC_Path_SR,
        .f = ~{
          DT <- IASDT.R::LoadAs(.x) %>%
            magrittr::inset(c("RC_Data_SR"), NULL)

          Quant <- DT$RC_Data_SR_Quant %>%
            dplyr::mutate(
              Variable = DT$Variable, Coords = DT$Coords, NFV = DT$NFV,
              .before = 1) %>%
            tidyr::pivot_wider(
              id_cols = c("Variable", "Coords", "NFV", "XVals"),
              names_from = Quantile, values_from = SR) %>%
            setNames(c("Variable", "Coords", "NFV", "XVals",
                       "Q25", "Q50", "Q975"))

          Observed <- DT$Observed_SR %>%
            dplyr::mutate(
              Variable = DT$Variable, Coords = DT$Coords, NFV = DT$NFV,
              .before = 1)
          Trend <- tibble::tibble(
            Variable = DT$Variable, Coords = DT$Coords, NFV = DT$NFV,
            Trend = DT$SR_PositiveTrendProb)

          return(list(Quant = Quant, Observed = Observed, Trend = list(Trend)))
        })) %>%
    dplyr::select(-Coords, -NFV, -RC_Path_SR) %>%
    tidyr::unnest_wider(DT) %>%
    tidyr::nest(.by = Variable) %>%
    dplyr::mutate(
      Quant = purrr::map(.x = data, .f = ~dplyr::bind_rows(.x$Quant)),
      Observed = purrr::map(.x = data, .f = ~dplyr::bind_rows(.x$Observed)),
      Trend = purrr::map(.x = data, .f = ~dplyr::bind_rows(.x$Trend))) %>%
    dplyr::select(-data) %>%
    dplyr::mutate(
      Plot = purrr::pmap(
        .l = list(Variable, Quant, Observed, Trend),
        .f = function(Variable, Quant, Observed, Trend) {

          IASDT.R::CatTime(paste0("  >>  ", Variable))
          # Maximum value on the y-axis
          PlotMax <- max(Observed$Pred, Quant$Q975 + 3)

          # Label on the y axis
          YAxisLabel <- "<span style='font-size: 12pt;'><b>Predicted species richness</span></b>"

          # Trend text
          DT_Trend <- Trend %>%
            dplyr::mutate(
              Trend2 = stringr::str_glue(
                "\n     Pr[pred(Var=max)] > Pr[pred(Var=min)] = {round(Trend, 2)}"),
              Trend2 = as.character(Trend2), X = -Inf, Y = Inf)

          # Variable long name (x-axis label)
          Variable2 <- dplyr::case_when(
            Variable == "bio4" ~ "<span style='font-size: 12pt;'><b>Bio4</b></span><span style='font-size: 9pt;'> (temperature seasonality)</span>",
            Variable == "bio6" ~ "<span style='font-size: 12pt;'><b>Bio6</b></span><span style='font-size: 9pt;'> (temperature of the coldest month)</span>",
            Variable == "bio8" ~ "<span style='font-size: 12pt;'><b>Bio8</b></span><span style='font-size: 9pt;'> (temperatures of the wettest quarter)</span>",
            Variable == "bio12" ~ "<span style='font-size: 12pt;'><b>Bio12</b></span><span style='font-size: 9pt;'> (annual precipitation amount)</span>",
            Variable == "bio15" ~ "<span style='font-size: 12pt;'><b>Bio15</b></span><span style='font-size: 9pt;'> (precipitation seasonality)</span>",
            Variable == "bio18" ~ "<span style='font-size: 12pt;'><b>Bio18</b></span><span style='font-size: 9pt;'> (monthly precipitation amount of the warmest quarter)</span>",
            Variable == "RoadRailLog" ~ "<span style='font-size: 12pt;'><b>Road + Rail intensity</b></span><span style='font-size: 9pt;'> (log<sub>10</sub>(x + 0.1))</span>",
            Variable == "BiasLog" ~ "<span style='font-size: 12pt;'><b>Sampling intensity</b></span><span style='font-size: 9pt;'> (log<sub>10</sub>(x + 0.1))</span>",
            .default = Variable)

          # facetting labels
          FacetLabel <- ggplot2::as_labeller(c(
            `1` = '<span style="font-size:12pt; color:red;"><b>non.focalVariables = 1</b></span><br><span style="font-size:7pt;">values of non-focal variables are set to the most likely value</span>',
            `2` = '<span style="font-size:12pt; color:red;"><b>non.focalVariables = 2</b></span><br><span style="font-size:7pt;">values of non-focal variables are set to the most likely value given the value of focal variable</span>',
            `c` = '<span style="font-size:10pt; color:red;"><b>Predictions at mean coordinates</b></span>',
            `i` = '<span style="font-size:10pt; color:red;"><b>Predictions without effect of spatial dependence</b></span>'))

          # Plot title
          TitleTxt <- paste0(
            '<span style="font-size:14pt; color:blue;"><b>Response curve  (predicted species richness) for ',
            Variable, "</b></span>")

          # Plot
          Plot <- ggplot2::ggplot(
            data = Observed, mapping = ggplot2::aes(x = XVals, y = Pred)) +
            ggplot2::geom_jitter(
              shape = 16, width = 0, height = 0.02, alpha = 0.2, size = 1,
              colour = "blue") +
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
              data = DT_Trend, mapping = ggplot2::aes(x = X, y = Y, label = Trend2),
              colour = "darkgrey", size = 2.75, vjust = 1.4, hjust = -0.05) +
            ggplot2::facet_grid(Coords ~ NFV, labeller = FacetLabel) +
            ggplot2::scale_y_continuous(
              limits = c(-1, PlotMax), oob = scales::squish, expand = c(0, 0)) +
            ggplot2::scale_x_continuous(expand = c(0.015, 0.015)) +
            ggplot2::xlab(Variable2) +
            ggplot2::ylab(YAxisLabel) +
            ggplot2::labs(title = TitleTxt) +
            ggplot2::theme_bw() +
            ggplot2::theme(
              strip.text = ggtext::element_markdown(
                margin = ggplot2::margin(0.05, 0.1, 0.05, 0.1, "cm")),
              strip.background = ggplot2::element_rect(colour = NA, fill = "white"),
              legend.position = "none",
              axis.title = ggtext::element_markdown(size = 10),
              plot.title = ggtext::element_markdown(
                margin = ggplot2::margin(2, 0, 6, 0)),
              plot.subtitle = ggtext::element_markdown(
                size = 7, colour = "darkgrey", margin = ggplot2::margin(4, 0, 4, 0)),
              axis.text = ggplot2::element_text(size = 8),
              panel.grid.major = ggplot2::element_line(linewidth = 0.25),
              panel.grid.minor = ggplot2::element_line(linewidth = 0.1),
              panel.spacing = ggplot2::unit(0.15, "lines"),
              plot.margin = ggplot2::unit(c(0.1, 0.2, 0.1, 0.2), "lines"))

          Path_JPEG <- file.path(
            Path_RC, "RespCurv_SR", paste0("RespCurv_SR_", Variable, ".jpeg"))
          ggplot2::ggsave(
            filename = Path_JPEG, plot = Plot, width = 22, height = 21,
            dpi = 600, units = "cm")

          return(invisible(NULL))
        }
      ))

  return(invisible(NULL))
}
