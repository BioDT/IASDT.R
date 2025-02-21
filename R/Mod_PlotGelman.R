## |------------------------------------------------------------------------| #
# PlotGelman ----
## |------------------------------------------------------------------------| #

#' Plot Gelman-Rubin-Brooks
#'
#' The `PlotGelman_*()` functions generate plots visualizing the evolution of
#' the Gelman-Rubin-Brooks shrink factor for different model parameters as the
#' number of iterations increases. These plots help assess whether MCMC chains
#' have converged to a common distribution. Each plot includes: median (solid
#' line) and 97.5<sup>th</sup> percentile (dashed line) of the shrink factor and
#' a dashed horizontal line at 1.1, representing the common convergence
#' threshold. The primary function for users is `PlotGelman()`, which internally
#' calls:
#' - `PlotGelman_Alpha()`: Plots shrink factor for the **Alpha** parameter
#' - `PlotGelman_Beta()`: Plots shrink factor for the **Beta** parameters
#' - `PlotGelman_Omega()`: Plots shrink factor for the **Omega** parameter
#' - `PlotGelman_Rho()`: Plots shrink factor for the **Rho** parameter
#' @param Path_Coda Character or `mcmc.list`. Path to a file containing the coda
#'   object, or an `mcmc.list` object representing MCMC samples.
#' @param Alpha,Beta,Omega,Rho Logical. If `TRUE`, plots the Gelman-Rubin
#'   statistic for the respective model parameters (Alpha, Beta, Omega, or Rho).
#'   Default: `TRUE` for all parameters.
#' @param FromHPC Logical. Whether the processing is being done on an 
#'   High-Performance Computing (HPC) environment, to adjust file paths 
#'   accordingly. Default: `TRUE`.
#' @param NOmega Integer. Number of species sampled for the Omega parameter.
#'   Default: 1000L.
#' @param PlottingAlpha Numeric. Transparency level (alpha) for plot lines (0 =
#'   fully transparent, 1 = fully opaque). Default: 0.25.
#' @param EnvFile Character. Path to the environment file containing paths to 
#'   data sources. Defaults to `.env`.
#' @param ReturnPlots Character. Path to the folder where the output plots will
#'   be saved.
#' @param CodaObj `mcmc.list`. An MCMC sample object containing posterior
#'   distributions from an Hmsc model.
#' @rdname PlotGelman
#' @name PlotGelman
#' @order 1
#' @export
#' @author Ahmed El-Gabbas

PlotGelman <- function(
    Path_Coda = NULL, Alpha = TRUE, Beta = TRUE, Omega = TRUE, Rho = TRUE,
    NOmega = 1000L, FromHPC = TRUE, PlottingAlpha = 0.25,
    EnvFile = ".env", ReturnPlots = FALSE) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Checking arguments --------

  if (sum(Alpha, Beta, Omega, Rho) == 0) {
    stop(
      "At least one of Alpha, Beta, Omega, and Rho must be `TRUE`",
      call. = FALSE)
  }

  if (is.null(Path_Coda)) {
    stop("Path_Coda cannot be empty", call. = FALSE)
  }

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(AllArgs, c("NOmega", "PlottingAlpha"), "numeric")
  IASDT.R::CheckArgs(
    AllArgs, c("Beta", "Rho", "Omega", "Alpha", "ReturnPlots"),
    "logical")

  rm(AllArgs, envir = environment())

  # # ..................................................................... ###

  # Loading coda object ------

  if (inherits(Path_Coda, "character")) {

    IASDT.R::CatTime("Loading coda object")
    CodaObj <- IASDT.R::LoadAs(Path_Coda)

  } else {

    if (!inherits(Path_Coda, "list")) {
      stop("`Path_Coda` is neither character path or a list", call. = FALSE)
    }
    if (!inherits(Path_Coda[[1]], "mcmc.list")) {
      stop("`Path_Coda` has no mcmc.list items", call. = FALSE)
    }

    CodaObj <- Path_Coda
    rm(Path_Coda, envir = environment())
  }

  # # ..................................................................... ###

  OutPath <- IASDT.R::Path(dirname(dirname(Path_Coda)), "Model_Convergence")
  fs::dir_create(OutPath)

  # # ..................................................................... ###

  # Alpha -----

  if (Alpha) {
    IASDT.R::CatTime("Alpha")
    PlotObj_Alpha <- PlotGelman_Alpha(
      CodaObj = CodaObj$Alpha[[1]], PlottingAlpha = PlottingAlpha)
  } else {
    PlotObj_Alpha <- NULL
  }

  # # ..................................................................... ###

  # Beta -----

  if (Beta) {
    IASDT.R::CatTime("Beta")
    PlotObj_Beta <- IASDT.R::PlotGelman_Beta(
      CodaObj = CodaObj$Beta, EnvFile = EnvFile,
      PlottingAlpha = PlottingAlpha, FromHPC = FromHPC)
  } else {
    PlotObj_Beta <- NULL
  }

  # # ..................................................................... ###

  # Omega -----

  if (Omega) {
    IASDT.R::CatTime("Omega")
    PlotObj_Omega <- IASDT.R::PlotGelman_Omega(
      CodaObj = CodaObj$Omega[[1]], NOmega = NOmega,
      PlottingAlpha = PlottingAlpha)
  } else {
    PlotObj_Omega <- NULL
  }

  # # ..................................................................... ###

  # Rho -----

  if (Rho && ("Rho" %in% names(CodaObj))) {
    IASDT.R::CatTime("Rho")
    PlotObj_Rho <- IASDT.R::PlotGelman_Rho(CodaObj$Rho)
  } else {
    PlotObj_Rho <- NULL
  }

  # # ..................................................................... ###

  # Saving plots as PDF -----

  PlotList <- list(
    Alpha = PlotObj_Alpha, Beta = PlotObj_Beta,
    Omega = PlotObj_Omega, Rho = PlotObj_Rho)

  PlotList4Plot <- purrr::list_flatten(purrr::discard(PlotList, is.null))

  if (length(PlotList4Plot) > 0) {
    # Using ggplot2::ggsave directly does not show non-ascii characters
    # correctly
    grDevices::cairo_pdf(
      filename = IASDT.R::Path(OutPath, "GelmanPlots.pdf"),
      width = 13, height = 7, onefile = TRUE)
    purrr::walk(PlotList4Plot, grid::grid.draw)
    grDevices::dev.off()
  } else {
    warning("No plots to save")
  }

  # # ..................................................................... ###

  # Saving plots as qs2 -----
  IASDT.R::SaveAs(
    InObj = PlotList, OutObj = "GelmanPlots",
    OutPath = IASDT.R::Path(OutPath, "GelmanPlots.qs2"))

  # # ..................................................................... ###

  IASDT.R::CatDiff(InitTime = .StartTime)

  if (ReturnPlots) {
    return(PlotList)
  } else {
    return(invisible(NULL))
  }
}

# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

## |------------------------------------------------------------------------| #
# PlotGelman_Alpha ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname PlotGelman
#' @name PlotGelman
#' @order 2
#' @author Ahmed El-Gabbas

PlotGelman_Alpha <- function(CodaObj, PlottingAlpha = 0.25) {

  # # ..................................................................... ###

  if (!inherits(CodaObj, "mcmc.list")) {
    stop("CodaObj has to be of class mcmc.list", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Var_LV <- Type <- Iter <- Plot <- data <- Median <- NULL

  # # ..................................................................... ###

  Gelman_Alpha_DT <- magrittr::extract2(CodaObj, 1) %>%
    attr("dimnames") %>%
    magrittr::extract2(2) %>%
    sort() %>%
    purrr::map_dfr(
      .f = function(x) {

        Alpha_Preplot1 <- lapply(CodaObj, function(Y) {
          Y[, x, drop = TRUE]
        }) %>%
          coda::mcmc.list()

        Alpha_Preplot2 <- try(
          gelman.preplot(
            x = Alpha_Preplot1,
            bin.width = 10, max.bins = 50, confidence = 0.95,
            transform = FALSE, autoburnin = TRUE),
          silent = TRUE)

        if (inherits(Alpha_Preplot2, "try-error")) {
          Alpha_Preplot2 <- gelman.preplot(
            x = Alpha_Preplot1,
            bin.width = 10, max.bins = 50, confidence = 0.95,
            transform = FALSE, autoburnin = FALSE)
        }

        Alpha_Preplot2 %>%
          magrittr::extract2("shrink") %>%
          tibble::as_tibble(rownames = "Iter") %>%
          purrr::set_names(c("Iter", "Median", "Q97_5")) %>%
          dplyr::filter(!is.nan(Median)) %>%
          dplyr::mutate(Iter = as.integer(Iter)) %>%
          tidyr::pivot_longer(
            cols = -Iter, names_to = "Type", values_to = "ShrinkFactor") %>%
          dplyr::arrange(Type, Iter) %>%
          dplyr::mutate(Type = factor(Type), Var_LV = x)
      })

  Gelman_Alpha_Plot <- Gelman_Alpha_DT %>%
    dplyr::mutate(
      group = paste0(Var_LV, "_", Type),
      Var_LV = purrr::map_chr(
        .x = Var_LV, .f = ~stringr::str_remove_all(.x, "Alpha1\\[|\\]"))) %>%
    tidyr::nest(.by = "Var_LV") %>%
    dplyr::mutate(
      Plot = purrr::map2(
        .x = Var_LV, .y = data,
        .f = ~{
          .y %>%
            ggplot2::ggplot() +
            ggplot2::geom_line(
              mapping = ggplot2::aes(
                x = Iter, y = ShrinkFactor, group = group, color = Type),
              alpha = PlottingAlpha) +
            ggplot2::scale_color_manual(
              values = c("Median" = "red", "Q97_5" = "black")) +
            ggplot2::geom_hline(
              yintercept = 1.1, linetype = "dashed", col = "darkgrey",
              linewidth = 0.8) +
            ggplot2::facet_grid(~ Type, labeller = ggplot2::label_parsed) +
            ggplot2::scale_x_continuous(
              limits = range(Gelman_Alpha_DT$Iter), expand = c(0, 0)) +
            ggplot2::coord_cartesian(expand = FALSE, clip = "off") +
            ggplot2::labs(
              title = "Gelman-Rubin-Brooks plot --- Alpha",
              subtitle = .x,
              caption = paste0(
                "This plot shows the evolution of Gelman and Rubin's shrink ",
                "factor as the number of iterations increases.")) +
            ggplot2::xlab(NULL) +
            ggplot2::ylab("Shrink factor") +
            ggplot2::theme_bw() +
            ggplot2::theme(
              plot.margin = ggplot2::margin(5, 20, 5, 5),
              legend.position = "none",
              legend.background = ggplot2::element_blank(),
              legend.title = ggplot2::element_blank(),
              strip.text = ggplot2::element_text(size = 12, face = "bold"),
              axis.title = ggplot2::element_text(
                size = 12, colour = "darkgrey", face = "bold"),
              axis.text = ggplot2::element_text(size = 12),
              plot.title = ggtext::element_markdown(
                size = 18, face = "bold", color = "blue"),
              plot.subtitle = ggplot2::element_text(
                size = 12, face = "italic", color = "darkgrey"),
              panel.spacing = ggplot2::unit(0.85, "lines"),
              plot.caption = ggplot2::element_text(
                color = "darkgrey", face = "italic", size = 10))
        })
    ) %>%
    dplyr::pull(Plot)

  return(Gelman_Alpha_Plot)
}


# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

## |------------------------------------------------------------------------| #
# PlotGelman_Beta ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname PlotGelman
#' @name PlotGelman
#' @order 3
#' @author Ahmed El-Gabbas

PlotGelman_Beta <- function(
    CodaObj, EnvFile = ".env", PlottingAlpha = 0.25, FromHPC = TRUE) {

  # # ..................................................................... ###

  if (is.null(CodaObj)) {
    stop("CodaObj cannot be empty", call. = FALSE)
  }

  if (!inherits(CodaObj, "mcmc.list")) {
    stop("CodaObj has to be of class mcmc.list", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Iter <- Type <- Var_Sp <- ShrinkFactor <- group <- NULL

  # # ..................................................................... ###

  Beta_Coda <- IASDT.R::Coda_to_tibble(
    CodaObj = CodaObj, Type = "beta", EnvFile = EnvFile, FromHPC = FromHPC)

  NVars <- length(unique(Beta_Coda$Variable))
  NSp <- length(unique(Beta_Coda$Species))
  SubTitle <- paste0(NVars, " covariates - ", NSp, " species")
  rm(Beta_Coda, envir = environment())

  Gelman_Beta_Vals <- magrittr::extract2(CodaObj, 1) %>%
    attr("dimnames") %>%
    magrittr::extract2(2) %>%
    sort() %>%
    purrr::map_dfr(
      .f = function(x) {

        Beta_Preplot1 <- lapply(CodaObj, function(Y) {
          Y[, x, drop = TRUE]
        }) %>%
          coda::mcmc.list()

        Beta_Preplot2 <- try(
          gelman.preplot(
            x = Beta_Preplot1,
            bin.width = 10, max.bins = 50, confidence = 0.95,
            transform = FALSE, autoburnin = TRUE),
          silent = TRUE)

        if (inherits(Beta_Preplot2, "try-error")) {
          Beta_Preplot2 <- gelman.preplot(
            x = Beta_Preplot1,
            bin.width = 10, max.bins = 50, confidence = 0.95,
            transform = FALSE, autoburnin = FALSE)
        }

        Beta_Preplot2 %>%
          magrittr::extract2("shrink") %>%
          tibble::as_tibble(rownames = "Iter") %>%
          purrr::set_names(c("Iter", "Median", "Q97_5")) %>%
          dplyr::mutate(Iter = as.integer(Iter)) %>%
          tidyr::pivot_longer(
            cols = -Iter, names_to = "Type", values_to = "ShrinkFactor") %>%
          dplyr::arrange(Type, Iter) %>%
          dplyr::mutate(Type = factor(Type), Var_Sp = x)
      }) %>%
    dplyr::mutate(group = paste0(Var_Sp, "_", Type))

  # # ..................................................................... ###

  Gelman_Beta_Plot <- Gelman_Beta_Vals %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(
      mapping = ggplot2::aes(
        x = Iter, y = ShrinkFactor, group = group, color = Type),
      alpha = PlottingAlpha) +
    ggplot2::scale_color_manual(
      values = c("Median" = "red", "Q97_5" = "black")) +
    ggplot2::geom_hline(
      yintercept = 1.1, linetype = "dashed", col = "darkgrey",
      linewidth = 0.8) +
    ggplot2::facet_grid(
      ~Type,
      labeller = ggplot2::as_labeller(
        c(`Median` = "Median", `Q97_5` = "97.5%"))) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Gelman-Rubin-Brooks plot --- Beta", subtitle = SubTitle,
      caption = paste0(
        "This plot shows the evolution of Gelman and Rubin's shrink factor as ",
        "the number of iterations increases.")) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("Shrink factor") +
    ggplot2::theme(
      legend.position = "none",
      legend.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      axis.title = ggplot2::element_text(
        size = 12, colour = "darkgrey", face = "bold"),
      axis.text = ggplot2::element_text(size = 12),
      plot.title = ggtext::element_markdown(
        size = 18, face = "bold", color = "blue"),
      plot.subtitle = ggplot2::element_text(
        size = 12, face = "italic", color = "darkgrey"),
      panel.spacing = ggplot2::unit(0.85, "lines"),
      plot.caption = ggplot2::element_text(
        color = "darkgrey", face = "italic", size = 10))

  # # ..................................................................... ###

  return(Gelman_Beta_Plot)
}

# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

## |------------------------------------------------------------------------| #
# PlotGelman_Omega ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname PlotGelman
#' @name PlotGelman
#' @order 4
#' @author Ahmed El-Gabbas

PlotGelman_Omega <- function(CodaObj, NOmega = 1000L, PlottingAlpha = 0.25) {

  # # ..................................................................... ###

  if (is.null(CodaObj)) {
    stop("CodaObj cannot be empty", call. = FALSE)
  }

  if (!is.numeric(NOmega) || NOmega <= 0) {
    stop("NOmega must be a positive integer.", call. = FALSE)
  }

  if (!inherits(CodaObj, "mcmc.list")) {
    stop("CodaObj has to be of class mcmc.list", call. = FALSE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Iter <- Type <- Sp_comb <- ShrinkFactor <- group <- NULL

  # # ..................................................................... ###

  Gelman_OmegaDT <- magrittr::extract2(CodaObj, 1) %>%
    attr("dimnames") %>%
    magrittr::extract2(2) %>%
    sample(min(NOmega, length(.))) %>%
    sort() %>%
    purrr::map_dfr(
      .f = function(x) {

        Omega_Preplot1 <- lapply(CodaObj, function(Y) {
          Y[, x, drop = TRUE]
        }) %>%
          coda::mcmc.list()

        Omega_Preplot2 <- try(
          gelman.preplot(
            x = Omega_Preplot1,
            bin.width = 10, max.bins = 50, confidence = 0.95,
            transform = FALSE, autoburnin = TRUE),
          silent = TRUE)

        if (inherits(Omega_Preplot2, "try-error")) {
          Omega_Preplot2 <- gelman.preplot(
            x = Omega_Preplot1,
            bin.width = 10, max.bins = 50, confidence = 0.95,
            transform = FALSE, autoburnin = FALSE)
        }

        Omega_Preplot2 %>%
          magrittr::extract2("shrink") %>%
          tibble::as_tibble(rownames = "Iter") %>%
          purrr::set_names(c("Iter", "Median", "Q97_5")) %>%
          dplyr::mutate(Iter = as.integer(Iter)) %>%
          tidyr::pivot_longer(
            cols = -Iter, names_to = "Type", values_to = "ShrinkFactor") %>%
          dplyr::arrange(Type, Iter) %>%
          dplyr::mutate(Type = factor(Type), Sp_comb = x)
      }) %>%
    dplyr::mutate(group = paste0(Sp_comb, "_", Type))

  # # ..................................................................... ###

  Gelman_Omega_Plot <- Gelman_OmegaDT %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(
      mapping = ggplot2::aes(
        x = Iter, y = ShrinkFactor, group = group, color = Type),
      alpha = PlottingAlpha) +
    ggplot2::scale_color_manual(
      values = c("Median" = "red", "Q97_5" = "black")) +
    ggplot2::geom_hline(
      yintercept = 1.1, linetype = "dashed", col = "darkgrey",
      linewidth = 0.8) +
    ggplot2::facet_grid(
      ~Type,
      labeller = ggplot2::as_labeller(
        c(`Median` = "Median", `Q97_5` = "97.5%"))) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::labs(
      title = paste0(
        "Gelman-Rubin-Brooks plot --- Omega --- ",
        NOmega, " species combination samples"),
      subtitle = NULL,
      caption = paste0(
        "This plot shows the evolution of Gelman and Rubin's shrink factor as ",
        "the number of iterations increases.")) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("Shrink factor") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5, 20, 5, 5),
      legend.position = "none",
      legend.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      axis.title = ggplot2::element_text(
        size = 12, colour = "darkgrey", face = "bold"),
      axis.text = ggplot2::element_text(size = 12),
      plot.title = ggtext::element_markdown(
        size = 18, face = "bold", color = "blue"),
      plot.subtitle = ggplot2::element_text(
        size = 12, face = "italic", color = "darkgrey"),
      panel.spacing = ggplot2::unit(0.85, "lines"),
      plot.caption = ggplot2::element_text(
        color = "darkgrey", face = "italic", size = 10))

  return(Gelman_Omega_Plot)
}


# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##
# # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ##

## |------------------------------------------------------------------------| #
# PlotGelman_Rho ----
## |------------------------------------------------------------------------| #

#' @export
#' @rdname PlotGelman
#' @name PlotGelman
#' @order 5
#' @author Ahmed El-Gabbas

PlotGelman_Rho <- function(CodaObj) {

  if (is.null(CodaObj)) {
    stop("CodaObj cannot be empty", call. = FALSE)
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Iter <- Type <- ShrinkFactor <- NULL

  Gelman_Rho_Plot <- try(
    gelman.preplot(
      x = CodaObj, bin.width = 10, max.bins = 50, confidence = 0.95,
      transform = FALSE, autoburnin = TRUE),
    silent = TRUE)

  if (inherits(Gelman_Rho_Plot, "try-error")) {
    Gelman_Rho_Plot <- gelman.preplot(
      x = CodaObj, bin.width = 10, max.bins = 50, confidence = 0.95,
      transform = FALSE, autoburnin = FALSE)
  }

  Gelman_Rho_Plot <- Gelman_Rho_Plot %>%
    magrittr::extract2("shrink") %>%
    tibble::as_tibble(rownames = "Iter") %>%
    purrr::set_names(c("Iter", "Median", "Q97_5")) %>%
    dplyr::mutate(Iter = as.integer(Iter)) %>%
    tidyr::pivot_longer(
      cols = -Iter, names_to = "Type", values_to = "ShrinkFactor") %>%
    dplyr::arrange(Type, Iter) %>%
    dplyr::mutate(Type = factor(Type)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(
      mapping = ggplot2::aes(x = Iter, y = ShrinkFactor, color = Type)) +
    ggplot2::scale_color_manual(
      values = c("Median" = "red", "Q97_5" = "black")) +
    ggplot2::geom_hline(
      yintercept = 1.1, linetype = "dashed", col = "darkgrey",
      linewidth = 0.8) +
    ggplot2::facet_grid(
      ~Type,
      labeller = ggplot2::as_labeller(
        c(`Median` = "Median", `Q97_5` = "97.5%"))) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::labs(
      title = "Gelman-Rubin-Brooks plot --- Rho",
      subtitle = NULL,
      caption = paste0(
        "This plot shows the evolution of Gelman and Rubin's ",
        "shrink factor as the number of iterations increases.")) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("Shrink factor") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5, 20, 5, 5),
      legend.position = "none",
      legend.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 14, face = "bold"),
      axis.title = ggplot2::element_text(
        size = 12, colour = "darkgrey", face = "bold"),
      axis.text = ggplot2::element_text(size = 12),
      plot.title = ggtext::element_markdown(
        size = 18, face = "bold", color = "blue"),
      plot.subtitle = ggplot2::element_text(
        size = 12, face = "italic", color = "darkgrey"),
      panel.spacing = ggplot2::unit(0.85, "lines"),
      plot.caption = ggplot2::element_text(
        color = "darkgrey", face = "italic", size = 10))

  return(Gelman_Rho_Plot)
}
