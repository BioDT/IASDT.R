## |------------------------------------------------------------------------| #
# Convergence_Plot_All ----
## |------------------------------------------------------------------------| #

#' Plot model convergence of multiple modelling alternatives
#'
#' This function generates and saves a series of diagnostic plots to assess the
#' convergence of Hmsc models across multiple modelling alternatives. It checks
#' model convergence using trace plots and Gelman-Rubin diagnostics for key
#' model parameters.
#'
#' @param ModelDir Character. Path to the root directory of the fitted model.
#'   The convergence outputs will be saved to the
#'   `Model_Convergence_All` subfolder.
#' @param NCores Integer. Number of CPU cores to use for parallel processing.
#' @name Convergence_Plot_All
#' @inheritParams Convergence_plots
#' @author Ahmed El-Gabbas
#' @return The function does not return anything but saves a series of
#'   diagnostic plots in the specified path.
#' @export

Convergence_Plot_All <- function(
    ModelDir = NULL, NOmega = 1000L, NCores = NULL, MarginType = "histogram") {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  if (is.null(ModelDir) || is.null(NCores)) {
    stop("ModelDir and NCores must not be NULL", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  GPP_Thin <- M_Name_Fit <- Tree <- rL <-
    M_thin <- M_samples <- Omega_Gelman <- Omega_ESS <- Beta_Gelman <-
    Beta_ESS <- ESS2 <- Path_Trace_Rho <- Rho <- Path_Trace_Alpha <-
    Path_Trace_Rho <- NULL

  # # ..................................................................... ###

  # Check input arguments ------

  IASDT.R::CatTime("Check input arguments")
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(AllArgs, get, envir = environment()) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character", Args = c("ModelDir", "MarginType"))
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric",  Args = c("NOmega", "NCores"))
  rm(AllArgs, envir = environment())

  if (length(MarginType) != 1) {
    stop("`MarginType` must be a single value.", call. = FALSE)
  }

  if (!MarginType %in% c("histogram", "density")) {
    stop(
      "`MarginType` must be either 'histogram' or 'density'.", call. = FALSE)
  }

  # # ..................................................................... ###

  # Prepare/load convergence data ------

  IASDT.R::CatTime("Prepare/load convergence data")

  Path_Convergence_All <- IASDT.R::Path(ModelDir, "Model_Convergence_All")
  Path_ConvDT <- IASDT.R::Path(Path_Convergence_All, "DT")
  fs::dir_create(c(Path_ConvDT, Path_Convergence_All))

  Model_Info <- IASDT.R::Path(ModelDir, "Model_Info.RData")
  if (!file.exists(Model_Info)) {
    stop("Model info file `", Model_Info, "` does not exist", call. = FALSE)
  }
  Model_Info <- IASDT.R::LoadAs(Model_Info)

  # Extract number of chains
  NChains <- length(Model_Info$Chain[[1]])

  # # ..................................................................... ###

  PrepConvergence <- function(ID) {

    Path_Coda <- Model_Info$Path_Coda[[ID]]
    Path_FittedMod <- Model_Info$Path_FittedMod[[ID]]
    M_Name_Fit <- Model_Info$M_Name_Fit[[ID]]
    Tree <- Model_Info$Tree[[ID]]

    CodaModelExist <- all(file.exists(c(Path_Coda, Path_FittedMod)))

    # Prepare traceplot ----

    if (isFALSE(CodaModelExist)) {

      Path_Trace_Rho <- Path_Trace_Alpha <- Beta_Gelman <-
        Beta_ESS <- Omega_Gelman <- Omega_ESS  <- NULL

    } else {

      Obj_Rho <- paste0(M_Name_Fit, "_TraceRho")
      Path_Trace_Rho <- IASDT.R::Path(Path_ConvDT, paste0(Obj_Rho, ".RData"))

      Obj_Alpha <- paste0(M_Name_Fit, "_TraceAlpha")
      Path_Trace_Alpha <- IASDT.R::Path(
        Path_ConvDT, paste0(Obj_Alpha, ".RData"))

      Obj_Beta_Omega <- paste0(M_Name_Fit, "_Beta_Omega")
      Path_Beta_Omega <- IASDT.R::Path(
        Path_ConvDT, paste0(Obj_Beta_Omega, ".RData"))


      if ((Tree == "Tree" && !file.exists(Path_Trace_Rho)) ||
          !file.exists(Path_Trace_Alpha) || !file.exists(Path_Beta_Omega)) {
        Model_Obj <- IASDT.R::LoadAs(Path_FittedMod)
        Coda_Obj <- IASDT.R::LoadAs(Path_Coda)
      }

      # Rho -----
      if (!file.exists(Path_Trace_Rho)) {
        if (Tree == "Tree") {
          RhoTitle <- stringr::str_remove_all(
            string = basename(Path_Coda), pattern = "_Tree|_Coda|.RData$|.qs2")

          PlotObj_Rho <- IASDT.R::Convergence_Rho(
            Post = Coda_Obj, Model = Model_Obj, Title = RhoTitle,
            MarginType = MarginType)

          IASDT.R::SaveAs(
            InObj = PlotObj_Rho, OutObj = Obj_Rho, OutPath = Path_Trace_Rho)

          rm(PlotObj_Rho, envir = environment())

        } else {
          Path_Trace_Rho <- NULL
        }
      }

      # Alpha -----
      if (!file.exists(Path_Trace_Alpha)) {
        PlotObj_Alpha <- IASDT.R::Convergence_Alpha(
          Post = Coda_Obj, Model = Model_Obj,
          Title = stringr::str_remove_all(
            basename(Path_Coda), "_Tree|_Coda|.RData$|.qs2"),
          MarginType = MarginType)

        IASDT.R::SaveAs(
          InObj = PlotObj_Alpha, OutObj = Obj_Alpha,
          OutPath = Path_Trace_Alpha)

        rm(PlotObj_Alpha, Model_Obj, envir = environment())
      }

      # Beta + Omega -----
      if (file.exists(Path_Beta_Omega)) {
        Beta_Omega <- IASDT.R::LoadAs(Path_Beta_Omega)
        Beta_Gelman <- Beta_Omega$Beta_Gelman
        Beta_ESS <- Beta_Omega$Beta_ESS
        Omega_ESS <- Beta_Omega$Omega_ESS
        Omega_Gelman <- Beta_Omega$Omega_Gelman
        rm(Beta_Omega, envir = environment())

      } else {

        Beta <- magrittr::extract2(Coda_Obj, "Beta")
        Omega <- magrittr::extract2(Coda_Obj, "Omega") %>%
          magrittr::extract2(1)

        rm(Coda_Obj, envir = environment())

        # BETA - effectiveSize
        Beta_ESS <- coda::effectiveSize(Beta)

        # BETA - gelman.diag
        Beta_Gelman <- coda::gelman.diag(Beta, multivariate = FALSE) %>%
          magrittr::extract2("psrf") %>%
          as.data.frame() %>%
          dplyr::pull(1) %>%
          magrittr::set_names(NULL)

        # OMEGA - effectiveSize
        Omega_ESS <- coda::effectiveSize(Omega)

        # OMEGA - gelman.diag
        sel <- sample.int(n = dim(Omega[[1]])[2], size = NOmega)

        Omega_Gelman <- purrr::map(.x = Omega, .f = ~ .x[, sel]) %>%
          coda::gelman.diag(multivariate = FALSE) %>%
          magrittr::extract2("psrf") %>%
          as.data.frame() %>%
          dplyr::pull(1) %>%
          magrittr::set_names(NULL)

        Beta_Omega <- list(
          Beta_Gelman = Beta_Gelman, Beta_ESS = Beta_ESS,
          Omega_Gelman = Omega_Gelman, Omega_ESS = Omega_ESS)
        save(Beta_Omega, file = Path_Beta_Omega)
        rm(Beta_Omega, envir = environment())
      }
    }

    invisible(gc())

    return(
      list(
        Path_Trace_Alpha = Path_Trace_Alpha,
        Path_Trace_Rho = Path_Trace_Rho,
        Beta_Gelman = Beta_Gelman, Beta_ESS = Beta_ESS,
        Omega_Gelman = Omega_Gelman, Omega_ESS = Omega_ESS))
  }

  # # ..................................................................... ###

  # Processing convergence data -----

  Path_DT <- IASDT.R::Path(Path_Convergence_All, "Convergence_DT.RData")

  if (IASDT.R::CheckData(Path_DT, warning = FALSE)) {

    IASDT.R::CatTime("Loading convergence data", Level = 1)
    Convergence_DT <- IASDT.R::LoadAs(Path_DT)

  } else {

    IASDT.R::CatTime("Processing convergence data", Level = 1)

    IASDT.R::CatTime(
      paste0("Prepare working on parallel using ", NCores, " cores."),
      Level = 2)

    if (NCores == 1) {
      future::plan("future::sequential", gc = TRUE)
    } else {
      withr::local_options(
        future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
        future.seed = TRUE)
      c1 <- snow::makeSOCKcluster(min(NCores, nrow(Model_Info)))
      on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
      future::plan("future::cluster", workers = c1, gc = TRUE)
      on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
    }

    Convergence_DT <- Model_Info %>%
      dplyr::mutate(
        Plots = future.apply::future_lapply(
          X = seq_len(nrow(Model_Info)), FUN = PrepConvergence,
          future.scheduling = Inf, future.seed = TRUE,
          future.packages = c(
            "dplyr", "sf", "Hmsc", "coda", "magrittr", "ggplot2",
            "magrittr", "IASDT.R"),
          future.globals = c(
            "Model_Info", "Path_ConvDT", "NOmega", "PrepConvergence"))) %>%
      dplyr::select(tidyselect::all_of(c("M_Name_Fit", "Plots"))) %>%
      tidyr::unnest_wider("Plots") %>%
      # arrange data alphanumerically by model name
      dplyr::arrange(gtools::mixedorder(M_Name_Fit)) %>%
      # discard non existed data
      dplyr::filter(!is.na(Path_Trace_Alpha))

    save(
      Convergence_DT,
      file = IASDT.R::Path(Path_Convergence_All, "Convergence_DT.RData"))

    if (NCores > 1) {
      snow::stopCluster(c1)
      future::plan("future::sequential", gc = TRUE)
    }
  }

  # # ..................................................................... ###

  # Plotting theme -----

  Theme <-  ggplot2::theme(
    strip.text = ggplot2::element_text(size = 16, face = "bold"),
    axis.title = ggtext::element_markdown(
      size = 20, colour = "darkgrey", face = "bold"),
    axis.text = ggplot2::element_text(size = 16),
    title = ggtext::element_markdown(size = 20, face = "bold", color = "blue"),
    axis.text.y = ggtext::element_markdown(
      hjust = 0, margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 5)),
    panel.spacing = ggplot2::unit(0.75, "lines"))

  Label <- ggplot2::as_labeller(c(
    `2000` = "2000 samples",
    `1000` = "1000 samples",
    `3000` = "3000 samples",
    `4000` = "4000 samples",
    `5000` = "5000 samples",
    Tree = "Phylogenetic (taxonomic) tree",
    NoTree = "No phylogenetic (taxonomic) tree"))

  # # ..................................................................... ###

  # Alpha - trace plots ------
  IASDT.R::CatTime("Alpha - trace plots")

  grDevices::cairo_pdf(
    filename = IASDT.R::Path(Path_Convergence_All, "TracePlots_Alpha.pdf"),
    width = 18, height = 12, onefile = TRUE)
  purrr::walk(
    .x = Convergence_DT$Path_Trace_Alpha,
    .f = purrr::safely(~{
      gridExtra::grid.arrange(IASDT.R::LoadAs(.x)[[1]])
    }))
  grDevices::dev.off()

  # # ..................................................................... ###

  # Rho - trace plots ------
  IASDT.R::CatTime("Rho - trace plots")

  layout_matrix <- matrix(seq_len(2 * 2), nrow = 2, byrow = TRUE)

  Plot <- Convergence_DT %>%
    dplyr::filter(stringr::str_detect(M_Name_Fit, "_Tree_")) %>%
    dplyr::mutate(
      Rho = purrr::map_if(
        .x = Path_Trace_Rho,
        .p = ~is.na(.x),
        .f = ~grid::grid.rect(gp = grid::gpar(col = "white")),
        .else = ~IASDT.R::LoadAs(.x))) %>%
    dplyr::pull(Rho) %>%
    gridExtra::marrangeGrob(
      bottom = bquote(paste0("page ", g, " of ", npages)),
      top = grid::textGrob(
        label = "Convergence of the rho parameter",
        gp = grid::gpar(fontface = "bold", fontsize = 20)),
      nrow = 2, ncol = 2, layout_matrix = layout_matrix)

  # Using ggplot2::ggsave directly does not show non-ascii characters correctly
  grDevices::cairo_pdf(
    filename = IASDT.R::Path(
      Path_Convergence_All, "TracePlots_Rho_Phylogenetic.pdf"),
    width = 18, height = 15, onefile = TRUE)
  invisible(print(Plot))
  grDevices::dev.off()

  # # ..................................................................... ###

  # Omega - Gelman convergence ------
  IASDT.R::CatTime("Omega - Gelman convergence")

  Plot_Path <- IASDT.R::Path(
    Path_Convergence_All, paste0("Convergence_Omega_Gelman.pdf"))

  Plot_Title <- paste0(
    "Gelman convergence diagnostic --- Omega (", NOmega, " samples)")

  Plot <- dplyr::left_join(Convergence_DT, Model_Info, by = "M_Name_Fit") %>%
    dplyr::select(rL, Tree, M_thin, M_samples, Omega_Gelman) %>%
    dplyr::mutate(
      GPP_Thin = paste0("GPP", rL, " | Th", M_thin),
      GPP_Thin = factor(
        GPP_Thin, levels = rev(gtools::mixedsort(unique(GPP_Thin))))) %>%
    tidyr::unnest("Omega_Gelman") %>%
    ggplot2::ggplot(ggplot2::aes(GPP_Thin, Omega_Gelman)) +
    ggplot2::geom_violin() +
    ggplot2::scale_y_log10() +
    ggplot2::facet_grid(Tree ~ M_samples, labeller = Label) +
    ggplot2::labs(title = Plot_Title) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab(
      "Gelman and Rubin's convergence diagnostic (log<sub>10</sub>)") +
    ggplot2::coord_flip(expand = FALSE) +
    Theme

  # Suppress the following message: Scale for y is already present. Adding
  # another scale for y, which will replace the existing scale.
  suppressMessages(suppressWarnings({
    Plot2 <- Plot +
      ggplot2::ylab(
        paste0(
          "Gelman and Rubin's convergence diagnostic ",
          "<sub>(only values between 0.9 and 1.1)</sub>")) +
      ggplot2::ylim(c(0.9, 1.1)) +
      Theme

    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    grDevices::cairo_pdf(
      filename = Plot_Path, width = 18, height = 12, onefile = TRUE)
    print(Plot)
    print(Plot2)
    grDevices::dev.off()

  }))

  # # ..................................................................... ###

  # Omega - Effective sample size -----
  IASDT.R::CatTime("Omega - Effective sample size")

  Plot_Path <- IASDT.R::Path(
    Path_Convergence_All, paste0("Convergence_Omega_ESS.pdf"))

  Plot_Title <- paste0(
    "Effective sample size --- Omega (", NOmega, " samples)")

  Plot <- Convergence_DT %>%
    dplyr::left_join(Model_Info, by = "M_Name_Fit") %>%
    dplyr::select(rL, Tree, M_thin, M_samples, Omega_ESS) %>%
    dplyr::mutate(
      GPP_Thin = paste0("GPP", rL, " | Th", M_thin),
      GPP_Thin = factor(
        GPP_Thin, levels = rev(gtools::mixedsort(unique(GPP_Thin))))) %>%
    tidyr::unnest("Omega_ESS") %>%
    ggplot2::ggplot(ggplot2::aes(GPP_Thin, Omega_ESS)) +
    ggplot2::geom_violin() +
    ggplot2::facet_grid(Tree ~ M_samples, labeller = Label) +
    ggplot2::labs(title = Plot_Title) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab(paste0("Effective sample size (", NChains, " chains)")) +
    ggplot2::coord_flip(expand = FALSE) +
    Theme

  Plot2 <- Convergence_DT %>%
    dplyr::left_join(Model_Info, by = "M_Name_Fit") %>%
    dplyr::select(rL, Tree, M_thin, M_samples, Omega_ESS) %>%
    dplyr::mutate(
      GPP_Thin = paste0("GPP", rL, " | Th", M_thin),
      GPP_Thin = factor(
        GPP_Thin, levels = rev(gtools::mixedsort(unique(GPP_Thin))))) %>%
    tidyr::unnest("Omega_ESS") %>%
    dplyr::mutate(ESS2 = (100 * Omega_ESS / (M_samples * NChains))) %>%
    ggplot2::ggplot(ggplot2::aes(GPP_Thin, ESS2)) +
    ggplot2::geom_violin() +
    ggplot2::facet_grid(Tree ~ M_samples, labeller = Label) +
    ggplot2::labs(title = Plot_Title) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("Mean effective sample size (%)") +
    ggplot2::coord_flip(expand = FALSE) +
    Theme

  # Using ggplot2::ggsave directly does not show non-ascii characters
  # correctly
  grDevices::cairo_pdf(
    filename = Plot_Path, width = 18, height = 12, onefile = TRUE)
  print(Plot)
  print(Plot2)
  grDevices::dev.off()

  # # ..................................................................... ###

  # Beta - Gelman convergence ------
  IASDT.R::CatTime("Beta - Gelman convergence")

  Plot_Title <- paste0("Gelman convergence diagnostic --- Beta")

  Plot_Path <- IASDT.R::Path(
    Path_Convergence_All, paste0("Convergence_Beta_Gelman.pdf"))

  Plot <- Convergence_DT %>%
    dplyr::left_join(Model_Info, by = "M_Name_Fit") %>%
    dplyr::select(rL, Tree, M_thin, M_samples, Beta_Gelman) %>%
    dplyr::mutate(
      GPP_Thin = paste0("GPP", rL, " | Th", M_thin),
      GPP_Thin = factor(
        GPP_Thin, levels = rev(gtools::mixedsort(unique(GPP_Thin))))) %>%
    tidyr::unnest("Beta_Gelman") %>%
    ggplot2::ggplot(ggplot2::aes(GPP_Thin, Beta_Gelman)) +
    ggplot2::geom_violin() +
    ggplot2::scale_y_log10() +
    ggplot2::facet_grid(Tree ~ M_samples, labeller = Label) +
    ggplot2::labs(title = Plot_Title) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab(
      "Gelman and Rubin's convergence diagnostic (log<sub>10</sub>)") +
    ggplot2::coord_flip(expand = FALSE) +
    Theme

  suppressMessages(suppressWarnings({
    Plot2 <- Plot +
      ggplot2::ylab(
        paste0(
          "Gelman and Rubin's convergence diagnostic ",
          "<sub>(only values between 0.9 and 1.1)</sub>")) +
      ggplot2::ylim(c(0.9, 1.1)) +
      Theme

    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    grDevices::cairo_pdf(
      filename = Plot_Path, width = 18, height = 12, onefile = TRUE)
    print(Plot)
    print(Plot2)
    grDevices::dev.off()
  }))

  # # ..................................................................... ###

  # Beta - Effective sample size -----
  IASDT.R::CatTime("Beta - Effective sample size")

  Plot_Path <- IASDT.R::Path(
    Path_Convergence_All, paste0("Convergence_Beta_ESS.pdf"))
  Plot_Title <- "Effective sample size --- Beta"

  Plot <- Convergence_DT %>%
    dplyr::left_join(Model_Info, by = "M_Name_Fit") %>%
    dplyr::select(rL, Tree, M_thin, M_samples, Beta_ESS) %>%
    dplyr::mutate(
      GPP_Thin = paste0("GPP", rL, " | Th", M_thin),
      GPP_Thin = factor(
        GPP_Thin, levels = rev(gtools::mixedsort(unique(GPP_Thin))))) %>%
    tidyr::unnest("Beta_ESS") %>%
    ggplot2::ggplot(ggplot2::aes(GPP_Thin, Beta_ESS)) +
    ggplot2::geom_violin() +
    ggplot2::facet_grid(Tree ~ M_samples, labeller = Label) +
    ggplot2::labs(title = Plot_Title) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab(paste0("Effective sample size (", NChains, " chains)")) +
    ggplot2::coord_flip(expand = FALSE) +
    Theme

  Plot2 <- Convergence_DT %>%
    dplyr::left_join(Model_Info, by = "M_Name_Fit") %>%
    dplyr::select(rL, Tree, M_thin, M_samples, Beta_ESS) %>%
    dplyr::mutate(
      GPP_Thin = paste0("GPP", rL, " | Th", M_thin),
      GPP_Thin = factor(
        GPP_Thin, levels = rev(gtools::mixedsort(unique(GPP_Thin))))) %>%
    tidyr::unnest("Beta_ESS") %>%
    dplyr::mutate(ESS2 = (100 * Beta_ESS / (M_samples * NChains))) %>%
    ggplot2::ggplot(ggplot2::aes(GPP_Thin, ESS2)) +
    ggplot2::geom_violin() +
    ggplot2::facet_grid(Tree ~ M_samples, labeller = Label) +
    ggplot2::labs(title = Plot_Title) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("Mean effective sample size (%)") +
    ggplot2::coord_flip(expand = FALSE) +
    Theme

  # # ..................................................................... ###

  # Using ggplot2::ggsave directly does not show non-ascii characters
  # correctly
  grDevices::cairo_pdf(
    filename = Plot_Path, width = 18, height = 12, onefile = TRUE)
  print(Plot)
  print(Plot2)
  grDevices::dev.off()

  # # ..................................................................... ###

  IASDT.R::CatDiff(
    InitTime = .StartTime, Prefix = "Plotting model convergence took ")

  # # ..................................................................... ###

  return(invisible(NULL))
}
