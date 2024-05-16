# |---------------------------------------------------| #
# PlotRespCurv ----
# |---------------------------------------------------| #

#' PlotRespCurv
#'
#' PlotRespCurv
#' @param Model Model
#' @param ResCurvePath ResCurvePath
#' @param ResCurvePrefix ResCurvePrefix
#' @param Verbose Verbose
#' @name PlotRespCurv
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

PlotRespCurv <- function(
    Model, ResCurvePath = getwd(), ResCurvePrefix = "ResCurv_", Verbose = TRUE) {

  IASDT.R::CatTime("Extract names of the variables")
  ModelVars <- stringr::str_split(
    as.character(Model$XFormula)[2], "\\+", simplify = TRUE) %>%
    stringr::str_trim()

  IASDT.R::CatTime("Prepare working on parallel")
  c1 <- snow::makeSOCKcluster(min(.NCores, nrow(ResCurvDT)))
  future::plan(cluster, workers = c1, gc = TRUE)

  IASDT.R::CatTime("Prepare response curve data on parallel - coordinates = i")
  ResCurvDT <- tidyr::expand_grid(Var = ModelVars, NFV = 1:2)

  ResCurvDT_i <- ResCurvDT %>%
    dplyr::mutate(
      ResCurv = furrr::future_map2(
        .x = Var, .y = NFV,
        .f = ~{
          Gradient <- Hmsc::constructGradient(
            hM = Model, focalVariable = .x, NFV = .y,
            coordinates = list(sample = "i"))
          Pred <- predict(
            object = Model, Gradient = Gradient, expected = TRUE,
            nParallel = 2, verbose = FALSE)
          return(list(Gradient = Gradient, Pred = Pred))
        },
        .options = furrr::furrr_options(seed = TRUE), .progress = FALSE))

  # IASDT.R::CatTime("Prepare response curve data on parallel - default value")
  # ResCurvDTDef <- ResCurvDT %>%
  #   dplyr::mutate(
  #     ResCurv = furrr::future_map2(
  #       .x = Var, .y = NFV,
  #       .f = ~{
  #         Gradient <- Hmsc::constructGradient(
  #           hM = Model, focalVariable = .x, NFV = .y)
  #         Pred <- predict(
  #           object = Model, Gradient = Gradient, expected = TRUE,
  #           nParallel = 2, verbose = FALSE)
  #         list(Gradient = Gradient, Pred = Pred)
  #       },
  #       .options = furrr::furrr_options(seed = TRUE), .progress = TRUE)) %>%
  #   tidyr::unnest_wider(ResCurv)

  IASDT.R::CatTime("Stopping the cluster")
  snow::stopCluster(c1)


  IASDT.R::CatTime("Prepare species response curves")
  SpNames <- ResCurvDT_i$Pred[[1]] %>%
    magrittr::extract2(1) %>%
    colnames()

  seq_along(SpNames) %>%
    purrr::walk(
      .f = ~{
        CurrSp <- SpNames[.x] %>%
          stringr::str_remove_all("Sp_") %>%
          as.numeric()

        SpName <- TaxaList %>%
          dplyr::filter(IAS_ID == CurrSp) %>%
          dplyr::pull(Species_name)

        if (Verbose) cat(paste0("  >>  ", CurrSp, " (", SpName, ")\n"))

        paste0(ResCurvePrefix,  SpNames[.x], ".pdf") %>%
          file.path(ResCurvePath, .) %>%
          grDevices::pdf(height = 9, width = 9, pointsize = 12)

        graphics::par(
          mfrow = c(2, 2),
          oma = c(0.05, 0.05, 0.05, 0.05), mar = c(2, 3, 3, 0.5))

        ResCurvDT_i %>%
          dplyr::mutate(
            Plot = purrr::pmap(
              .l = list(Var, NFV, Gradient, Pred),
              .f = function(Var, NFV, Gradient, Pred) {
                MainLabel <- paste0(Var, " - non.focalVariables = ", NFV)
                Support <- Hmsc::plotGradient(
                  hM = Model, Gradient = Gradient, predY = Pred, measure = "Y",
                  showPosteriorSupport = FALSE, index = .x,  showData = TRUE,
                  jigger = 0.0125, cicol = rgb(0, 0, 1, alpha = 0.2),
                  main = MainLabel, las = 1, ylabel = NA)
                mtext(
                  text = paste0("Posterior support = ", Support),
                  padj = 0, adj = 0, cex = 0.65, col = "darkgrey")
              })) %>%
          invisible()

        grDevices::dev.off()

      })
  return(ResCurvDT_i)
}
