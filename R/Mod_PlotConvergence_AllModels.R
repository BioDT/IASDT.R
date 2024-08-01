## |------------------------------------------------------------------------| #
# PlotConvergence_AllModels ----
## |------------------------------------------------------------------------| #

#' Plot model convergence of multiple modelling alternatives
#'
#' This function generates and saves a series of diagnostic plots to assess the
#' convergence of Hmsc models across multiple modelling alternatives. It checks
#' model convergence using trace plots and Gelman-Rubin diagnostics for key
#' model parameters.
#' @param Path_Model String. Path to save all the output, including the to be
#'   fitted models (without trailing slash). Must not be `NULL`.
#' @param EnvFile String. Path to read the environment variables. Default value:
#'   `.env`.
#' @param FromHPC Logical. Indicates whether the function is being run on a
#'   High-Performance Computing (HPC) environment. Adjusts file paths
#'   accordingly.
#' @param NChains Integer. Number of MCMC chains used in the model fitting
#'   process.
#' @param maxOmega Integer. Maximum number of species interactions to sample for
#'   convergence diagnostics.
#' @param NCores Integer. Number of cores to use for parallel processing.
#' @name PlotConvergence_AllModels
#' @author Ahmed El-Gabbas
#' @return The function does not return anything but saves a series of
#'   diagnostic plots in the specified path.
#' @details The function reads the following environment variable:
#'    - **`Path_LUMI_Scratch`** (only if `FromHPC` = `TRUE`) for the path of
#'    the scratch folder of the `BioDT` project on LUMI.
#' @export

PlotConvergence_AllModels <- function(
    Path_Model = NULL, EnvFile = ".env", FromHPC = TRUE, NChains = 4,
    maxOmega = 1000, NCores = NULL) {

  if (is.null(Path_Model) || is.null(NCores)) {
    stop("Path_Model and NCores must not be NULL")
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  GPP_Thin <- Path_Coda <- Path_FittedMod <- M_Name_Fit <- Tree <- Plots <-
    rL <- M_thin <- M_samples <- Omega_Gelman <- Omega_ESS <- Beta_Gelman <-
    Beta_ESS <- ESS2 <- Path_Trace_Rho <- Rho <- NULL

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Check input arguments
  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Check input arguments")
  AllArgs <- ls()
  AllArgs <- purrr::map(
    AllArgs, ~get(.x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character", Args = c("Path_Model", "EnvFile"))
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "logical", Args = "FromHPC")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric",
    Args = c("NChains", "maxOmega", "NCores"))

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Load environment variables
  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Load environment variables")

  if (FromHPC) {
    if (magrittr::not(file.exists(EnvFile))) {
      stop(paste0(
        "Path for environment variables: ", EnvFile, " was not found"))
    }

    readRenviron(EnvFile)
    Path_Scratch <- Sys.getenv("Path_LUMI_Scratch")

    if (Path_Scratch == "") {
      stop("Path_Scratch environment variable not found.")
    }
    if (magrittr::not(dir.exists(Path_Scratch))) {
      stop("The scratch folder does not exist")
    }

    InitialWD <- getwd()
    setwd(Path_Scratch)
    on.exit(setwd(InitialWD), add = TRUE)

  }

  Path_Convergence_All <- file.path(Path_Model, "Model_Convergence_All")
  Path_ConvDT <- file.path(Path_Convergence_All, "DT")
  fs::dir_create(c(Path_ConvDT, Path_Convergence_All))

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  ## Prepare convergence data ------
  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Prepare convergence data")

  c1 <- snow::makeSOCKcluster(NCores)
  on.exit(invisible(try(snow::stopCluster(c1), silent = TRUE)), add = TRUE)
  future::plan(future::cluster, workers = c1, gc = TRUE)

  Model_Info <- IASDT.R::LoadAs(file.path(Path_Model, "Model_Info.RData"))

  Convergence_DT <- Model_Info %>%
    dplyr::mutate(
      Plots = furrr::future_pmap(
        .l = list(Path_Coda, Path_FittedMod, M_Name_Fit, Tree),
        .f = function(Path_Coda, Path_FittedMod, M_Name_Fit, Tree) {

          CodaModelExist <- all(file.exists(c(Path_Coda, Path_FittedMod)))

          # prepare traceplot ----
          if (CodaModelExist) {
            Model_Obj <- IASDT.R::LoadAs(Path_FittedMod)
            Coda_Obj <- IASDT.R::LoadAs(Path_Coda)

            # Rho -----
            ObjName_Rho <- paste0(M_Name_Fit, "_TraceRho")
            Path_Trace_Rho <- file.path(
              Path_ConvDT, paste0(ObjName_Rho, ".RData"))


            if (magrittr::not(file.exists(Path_Trace_Rho))) {
              if (Tree == "Tree") {
                RhoTitle <- stringr::str_remove_all(
                  string = basename(Path_Coda), pattern = "_Tree|_Coda.RData$")

                PlotObj_Rho <- IASDT.R::PlotRho(
                  Post = Coda_Obj, Model = Model_Obj, Title = RhoTitle)

                IASDT.R::SaveAs(
                  InObj = PlotObj_Rho, OutObj = ObjName_Rho,
                  OutPath = Path_Trace_Rho)

                rm(PlotObj_Rho)
              } else {
                Path_Trace_Rho <- NULL
              }
            }

            # Alpha -----
            ObjName_Alpha <- paste0(M_Name_Fit, "_TraceAlpha")
            Path_Trace_Alpha <- file.path(
              Path_ConvDT, paste0(ObjName_Alpha, ".RData"))

            if (magrittr::not(file.exists(Path_Trace_Alpha))) {
              PlotObj_Alpha <- IASDT.R::PlotAlpha(
                Post = Coda_Obj, Model = Model_Obj, NRC = c(2, 3),
                Title = stringr::str_remove_all(
                  basename(Path_Coda), "_Tree|_Coda.RData$"))

              IASDT.R::SaveAs(
                InObj = PlotObj_Alpha, OutObj = ObjName_Alpha,
                OutPath = Path_Trace_Alpha)

              rm(PlotObj_Alpha, Model_Obj)
            }

            # Beta + Omega -----
            ObjName_Beta_Omega <- paste0(M_Name_Fit, "_Beta_Omega")
            Path_Beta_Omega <- file.path(
              Path_ConvDT, paste0(ObjName_Beta_Omega, ".RData"))

            if (file.exists(Path_Beta_Omega)) {
              Beta_Omega <- IASDT.R::LoadAs(Path_Beta_Omega)
              Beta_Gelman <- Beta_Omega$Beta_Gelman
              Beta_ESS <- Beta_Omega$Beta_ESS
              Omega_ESS <- Beta_Omega$Omega_ESS
              Omega_Gelman <- Beta_Omega$Omega_Gelman
              rm(Beta_Omega)
            } else {

              Beta <- magrittr::extract2(Coda_Obj, "Beta")
              Omega <- magrittr::extract2(Coda_Obj, "Omega") %>%
                magrittr::extract2(1)
              rm(Coda_Obj)

              # BETA - effectiveSize
              Beta_ESS <- coda::effectiveSize(Beta)

              # BETA - gelman.diag
              Beta_Gelman <- Beta %>%
                coda::gelman.diag(multivariate = FALSE) %>%
                magrittr::extract2("psrf") %>%
                as.data.frame() %>%
                dplyr::pull(1) %>%
                magrittr::set_names(NULL)

              # OMEGA - effectiveSize
              Omega_ESS <- coda::effectiveSize(Omega)

              # OMEGA - gelman.diag
              sel <- sample(seq_len(dim(Omega[[1]])[2]), size = maxOmega)
              Omega_Gelman <- Omega %>%
                purrr::map(~ .x[, sel]) %>%
                coda::gelman.diag(multivariate = FALSE) %>%
                magrittr::extract2("psrf") %>%
                as.data.frame() %>%
                dplyr::pull(1) %>%
                magrittr::set_names(NULL)

              Beta_Omega <- list(
                Beta_Gelman = Beta_Gelman, Beta_ESS = Beta_ESS,
                Omega_Gelman = Omega_Gelman, Omega_ESS = Omega_ESS)
              save(Beta_Omega, file = Path_Beta_Omega)
              rm(Beta_Omega)
            }
          } else {
            Path_Trace_Rho <- Path_Trace_Alpha <- Beta_Gelman <-
              Beta_ESS <- Omega_Gelman <- Omega_ESS  <- NULL
          }

          list(
            Path_Trace_Alpha = Path_Trace_Alpha,
            Path_Trace_Rho = Path_Trace_Rho,
            Beta_Gelman = Beta_Gelman, Beta_ESS = Beta_ESS,
            Omega_Gelman = Omega_Gelman, Omega_ESS = Omega_ESS) %>%
            return()
        },
        .progress = FALSE, .options = furrr::furrr_options(seed = TRUE))) %>%
    dplyr::select(M_Name_Fit, Plots) %>%
    tidyr::unnest_wider("Plots")

  save(Convergence_DT,
       file = file.path(Path_Convergence_All, "Convergence_DT.RData"))

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Plotting theme -----
  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  Theme <-  ggplot2::theme(
    strip.text = ggplot2::element_text(size = 16, face = "bold"),
    axis.title = ggplot2::element_text(
      size = 20, colour = "darkgrey", face = "bold"),
    axis.text = ggplot2::element_text(size = 16),
    title = ggplot2::element_text(size = 20, face = "bold", color = "blue"),
    axis.text.y = ggplot2::element_text(
      hjust = 0, margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 5)),
    panel.spacing = ggplot2::unit(0.75, "lines"))

  Label <- ggplot2::as_labeller(c(
    `2000` = "2000 samples",
    `1000` = "1000 samples",
    `3000` = "3000 samples",
    `4000` = "4000 samples",
    `5000` = "5000 samples",
    `Tree` = "Phylogenetic (taxonomic) tree",
    `NoTree` = "No phylogenetic (taxonomic) tree"))

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  ## Alpha - trace plots ------
  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Alpha - trace plots")

  invisible({
    grDevices::pdf(
      file = file.path(Path_Convergence_All, "TracePlots_Alpha.pdf"),
      width = 18, height = 12)
    Convergence_DT$Path_Trace_Alpha %>%
      purrr::map(purrr::safely(~print(IASDT.R::LoadAs(.x))))
    grDevices::dev.off()
  })

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  ## Rho - trace plots ------
  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Rho - trace plots")

  Convergence_DT %>%
    dplyr::filter(stringr::str_detect(M_Name_Fit, "_Tree_")) %>%
    dplyr::mutate(
      Rho = purrr::map_if(
        .x = Path_Trace_Rho,
        .p = ~is.na(.x),
        .f = ~grid::grid.rect(gp = grid::gpar(col = "white")),
        .else = ~IASDT.R::LoadAs(.x))
    ) %>%
    dplyr::pull(Rho) %>%
    gridExtra::marrangeGrob(
      bottom = bquote(paste0("page ", g, " of ", npages)),
      top = grid::textGrob(
        label = "Convergence of the rho parameter",
        gp = grid::gpar(fontface = "bold", fontsize = 20)),
      nrow = 2, ncol = 3) %>%
    ggplot2::ggsave(
      dpi = 600, device = "pdf", width = 18, height = 12,
      filename = file.path(
        Path_Convergence_All, "TracePlots_Rho_Phylogenetic.pdf"))

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  ## Omega - Gelman convergence ------
  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Omega - Gelman convergence")

  Plot_Path <- file.path(
    Path_Convergence_All, paste0("Convergence_Omega_Gelman.pdf"))

  Plot_Title <- paste0(
    "Gelman convergence diagnostic - Omega (", maxOmega, " samples)")

  Plot <- dplyr::left_join(Convergence_DT, Model_Info, by = "M_Name_Fit") %>%
    dplyr::select(rL, Tree, M_thin, M_samples, Omega_Gelman) %>%
    dplyr::mutate(
      GPP_Thin = paste0("GPP", rL, " | Th", M_thin),
      GPP_Thin = factor(
        GPP_Thin, levels = rev(gtools::mixedsort(unique(GPP_Thin))))) %>%
    tidyr::unnest("Omega_Gelman") %>%
    ggplot2::ggplot(ggplot2::aes(GPP_Thin, Omega_Gelman)) +
    ggplot2::geom_violin() +
    ggplot2::geom_vline(
      xintercept = c(4.5, 8.5, 12.5), linetype = "dashed", color = "blue") +
    ggplot2::scale_y_log10() +
    ggplot2::facet_grid(Tree ~ M_samples, labeller = Label) +
    ggplot2::labs(title = Plot_Title) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("Gelman and Rubin's convergence diagnostic (log10)") +
    ggplot2::coord_flip(expand = FALSE) +
    Theme

  Plot2 <- Plot +
    ggplot2::ylab("Gelman and Rubin's convergence diagnostic - (cropped)") +
    ggplot2::coord_flip(expand = FALSE, ylim = c(0.995, 1.05))

  ggplot2::ggsave(
    filename = Plot_Path,
    plot = gridExtra::marrangeGrob(
      grobs = list(Plot, Plot2), nrow = 1, ncol = 1, top = NULL),
    width = 18, height = 12)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  ## Omega - Effective sample size -----
  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Omega - Effective sample size")

  Plot_Path <- file.path(
    Path_Convergence_All, paste0("Convergence_Omega_ESS.pdf"))

  Plot_Title <- paste0("Effective sample size - Omega (", maxOmega, " samples)")

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
    ggplot2::geom_vline(
      xintercept = c(4.5, 8.5, 12.5), linetype = "dashed", color = "blue") +
    ggplot2::facet_grid(Tree ~ M_samples, labeller = Label) +
    ggplot2::labs(title = Plot_Title) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab(paste0("Effective sample size - ", NChains, " chains")) +
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
    ggplot2::geom_vline(
      xintercept = c(4.5, 8.5, 12.5), linetype = "dashed", color = "blue") +
    ggplot2::facet_grid(Tree ~ M_samples, labeller = Label) +
    ggplot2::labs(title = Plot_Title) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("Mean effective sample size (%)") +
    ggplot2::coord_flip(expand = FALSE) +
    Theme

  ggplot2::ggsave(
    filename = Plot_Path,
    plot = gridExtra::marrangeGrob(
      grobs = list(Plot, Plot2), nrow = 1, ncol = 1, top = NULL),
    width = 18, height = 12)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  ## Beta - Gelman convergence ------
  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Beta - Gelman convergence")

  Plot_Title <- paste0("Gelman convergence diagnostic - Beta")

  Plot_Path <- file.path(
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
    ggplot2::geom_vline(
      xintercept = c(4.5, 8.5, 12.5), linetype = "dashed", color = "blue") +
    ggplot2::scale_y_log10() +
    ggplot2::facet_grid(Tree ~ M_samples, labeller = Label) +
    ggplot2::labs(title = Plot_Title) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("Gelman and Rubin's convergence diagnostic (log10)") +
    ggplot2::coord_flip(expand = FALSE) +
    Theme

  Plot2 <- Plot +
    ggplot2::ylab("Gelman and Rubin's convergence diagnostic (cropped)") +
    ggplot2::coord_flip(expand = FALSE, ylim = c(0.995, 1.05))

  ggplot2::ggsave(
    filename = Plot_Path,
    plot = gridExtra::marrangeGrob(
      grobs = list(Plot, Plot2), nrow = 1, ncol = 1, top = NULL),
    width = 18, height = 12)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  ## Beta - Effective sample size -----
  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Beta - Effective sample size")

  Plot_Path <- file.path(
    Path_Convergence_All, paste0("Convergence_Beta_ESS.pdf"))
  Plot_Title <- "Effective sample size - Beta"

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
    ggplot2::geom_vline(
      xintercept = c(4.5, 8.5, 12.5), linetype = "dashed", color = "blue") +
    ggplot2::facet_grid(Tree ~ M_samples, labeller = Label) +
    ggplot2::labs(title = Plot_Title) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab(paste0("Effective sample size - ", NChains, " chains")) +
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
    ggplot2::geom_vline(
      xintercept = c(4.5, 8.5, 12.5), linetype = "dashed", color = "blue") +
    ggplot2::facet_grid(Tree ~ M_samples, labeller = Label) +
    ggplot2::labs(title = Plot_Title) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("Mean effective sample size (%)") +
    ggplot2::coord_flip(expand = FALSE) +
    Theme

  ggplot2::ggsave(
    filename = Plot_Path,
    plot = gridExtra::marrangeGrob(
      grobs = list(Plot, Plot2), nrow = 1, ncol = 1, top = NULL),
    width = 18, height = 12)

  return(invisible(NULL))
}
