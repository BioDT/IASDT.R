## |------------------------------------------------------------------------| #
# PlotConvergence ----
## |------------------------------------------------------------------------| #

#' Plot model convergence of a selected model
#'
#' Plot model convergence of a selected model
#'
#' @param Path_Coda String. Path to the coda object.
#' @param Path_FittedModel String. Path to the fitted Hmsc model object
#' @param EnvFile String. Path to read the environment variables. Default value: `.env`
#' @param FromHPC Logical. Work from HPC? This is to adjust the file paths.
#' @param NChains Integer. Number of model chains
#' @param Title String. Title of the rho and alpha plots
#' @param NOmega Number of sample species interactions
#' @param NCores Number of parallel cores for parallelization
#' @param NRC Vector. Number of rows and columns
#' @param SavePlotData Logical. Save plots as RData file
#' @param Cols Colours for lines for each chain
#' @name PlotConvergence
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

PlotConvergence <- function(
    Path_Coda = NULL, Path_FittedModel = NULL, EnvFile = ".env",
    FromHPC = TRUE, NChains = 4, Title = " ", NOmega = 1000, NCores = NULL,
    NRC = c(2, 3), SavePlotData = FALSE,
    Cols = c("red", "blue", "darkgreen", "darkgrey")) {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  SpComb <- `25%` <- `75%` <- Class <- Order <- Family <- Plot <- DT <-
    IAS_ID <- Species <- Variable <- Chain <- Iter <- Value <- Var_Sp <-
    coda <- data <- dplyr <- ggExtra <- ggplot2 <- ggtext <- label <- y <-
    Var <- PlotFixedY <- Var <- NULL

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Check input arguments ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Check input arguments")

  AllArgs <- ls()
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "character",
    Args = c("Path_Coda", "Path_FittedModel", "EnvFile"))

  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "logical", Args = "FromHPC")
  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Type = "numeric",
    Args = c("NChains", "NOmega", "NCores", "NRC"))
  rm(AllArgs)

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Load environment variables ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Load environment variables")

  if (file.exists(EnvFile)) {
    readRenviron(EnvFile)
    Path_Scratch <- Sys.getenv("Path_LUMI_Scratch")
    if (FromHPC) {
      setwd(Path_Scratch)
    }
  } else {
    MSG <- paste0("Path for environment variables: ", EnvFile, " was not found")
    stop(MSG)
  }

  IASDT.R::CatTime("Create path")
  Path_Convergence <- Path_Coda %>%
    dirname() %>%
    dirname() %>%
    file.path("Model_Convergence")
  Path_Convergence_BySp <- file.path(Path_Convergence, "Beta_BySpecies")
  fs::dir_create(c(Path_Convergence, Path_Convergence_BySp))

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Prepare convergence data ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Prepare convergence data")

  IASDT.R::CatTime("  >>  Loading coda object")
  Coda_Obj <- IASDT.R::LoadAs(Path_Coda)

  IASDT.R::CatTime("  >>  Loading fitted model object")
  Model <- IASDT.R::LoadAs(Path_FittedModel)

  Obj_Omega <- Coda_Obj$Omega[[1]]
  Obj_Beta <- Coda_Obj$Beta

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Rho ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Rho")

  IASDT.R::CatTime("  >>  Prepare plot")
  PlotObj_Rho <- IASDT.R::PlotRho(
    Post = Coda_Obj, Model = Model, Title = Title, Cols = Cols)

  IASDT.R::CatTime("  >>  Save plot")
  ggplot2::ggsave(
    plot = PlotObj_Rho, dpi = 600, device = "pdf", width = 18, height = 12,
    filename = file.path(Path_Convergence, "Convergence_Rho.pdf"))

  if (SavePlotData) {
    IASDT.R::CatTime("  >>  Save plotting data")
    IASDT.R::SaveAs(
      InObj = PlotObj_Rho, OutObj = "Convergence_Rho",
      OutPath = file.path(Path_Convergence, "Convergence_Rho.RData"))
  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Alpha  ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Alpha")

  IASDT.R::CatTime("  >>  Prepare plot")
  PlotObj_Alpha <- IASDT.R::PlotAlpha(
    Post = Coda_Obj, Model = Model, Title = Title, NRC = NRC,
    AddFooter = FALSE, AddTitle = FALSE, Cols = Cols)

  IASDT.R::CatTime("  >>  Save plot")
  ggplot2::ggsave(
    plot = PlotObj_Alpha, dpi = 600, device = "pdf", width = 18, height = 12,
    filename = file.path(Path_Convergence, "Convergence_Alpha.pdf"))

  if (SavePlotData) {
    IASDT.R::CatTime("  >>  Save plotting data")
    IASDT.R::SaveAs(
      InObj = PlotObj_Alpha, OutObj = "Convergence_Alpha",
      OutPath = file.path(Path_Convergence, "Convergence_Alpha.RData"))
  }

  rm(Model, Coda_Obj, PlotObj_Alpha, PlotObj_Rho)
  invisible(gc())

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Omega  ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Omega")

  IASDT.R::CatTime("  >>  Coda to tibble")
  OmegaDF <- IASDT.R::Coda_to_tibble(
    CodaObj = Obj_Omega, Type = "omega", NOmega = NOmega, EnvFile = EnvFile)

  SelectedCombs <- unique(OmegaDF$SpComb)

  IASDT.R::CatTime("  >>  Prepare confidence interval data")
  CI <- Obj_Omega %>%
    purrr::map(.f = ~.x[, SelectedCombs]) %>%
    coda::mcmc.list() %>%
    summary(quantiles = c(0.25, 0.75)) %>%
    magrittr::extract2("quantiles") %>%
    as.data.frame() %>%
    tibble::as_tibble(rownames = "SpComb")

  IASDT.R::CatTime("  >>  Prepare working in parallel")
  c1 <- snow::makeSOCKcluster(NCores)
  future::plan(future::cluster, workers = c1, gc = TRUE)
  invisible(snow::clusterEvalQ(
    cl = c1, IASDT.R::LoadPackages(dplyr, coda, ggplot2, ggExtra, ggtext)))
  snow::clusterExport(
    cl = c1, list = c("NOmega", "CI", "SelectedCombs", "OmegaDF", "Cols"),
    envir = environment())

  IASDT.R::CatTime("  >>  Prepare trace plots data")
  OmegaTracePlots <- snow::parLapply(
    cl = c1, x = seq_len(NOmega),
    fun = function(x) {

      CI0 <- CI %>%
        dplyr::filter(SpComb == SelectedCombs[x]) %>%
        dplyr::select(-SpComb) %>%
        round(2) %>%
        unlist()

      CI1 <- CI0 %>%
        paste0(collapse = " - ") %>%
        paste0("<i>50% credible interval:</i> ", .) %>%
        data.frame(x = -Inf, y = -Inf, label = .)

      PanelTitle <- data.frame(
        x = Inf, y = Inf,
        label = paste0(OmegaDF$IAS1[[x]], " & ", OmegaDF$IAS2[[x]]))

      OmegaDF2 <- dplyr::filter(OmegaDF, SpComb == SelectedCombs[x]) %>%
        dplyr::pull("DT") %>%
        magrittr::extract2(1)

      ValRange <- range(OmegaDF2$Value)

      Plot <- OmegaDF2 %>%
        ggplot2::ggplot(
          mapping = ggplot2::aes(x = Iter, y = Value, color = factor(Chain))) +
        ggplot2::geom_line(linewidth = 0.15, alpha = 0.6) +
        ggplot2::geom_smooth(
          method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.8) +
        ggplot2::geom_point(alpha = 0) +
        ggplot2::geom_hline(
          yintercept = CI0, linetype = "dashed", color = "black", linewidth = 1) +
        ggplot2::scale_color_manual(values = Cols) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggtext::geom_richtext(
          mapping = ggplot2::aes(x = x, y = y, label = label),
          data = CI1, inherit.aes = FALSE, hjust = 0, vjust = 0, lineheight = 0,
          fill = NA, label.color = NA) +
        ggtext::geom_richtext(
          mapping = ggplot2::aes(x = x, y = y, label = label),
          data = PanelTitle, inherit.aes = FALSE, colour = "blue",
          hjust = 1, vjust = 1, lineheight = 0, fill = NA, label.color = NA) +
        ggplot2::theme_bw() +
        ggplot2::xlab(NULL) +
        ggplot2::ylab(NULL) +
        ggplot2::theme(legend.position = "none")

      if (dplyr::between(0, ValRange[1], ValRange[2])) {
        Plot <- Plot +
          ggplot2::geom_hline(
            yintercept = 0, linetype = "solid", color = "green",
            linewidth = 0.4)
      }

      Plot <- ggExtra::ggMarginal(
        p = Plot, type = "density", margins = "y", size = 5,
        color = "steelblue4")

      invisible(gc())
      return(Plot)
    })

  # Save plot
  IASDT.R::CatTime("  >>  Save plots")
  OmegaTracePlots %>%
    gridExtra::marrangeGrob(
      bottom = bquote(paste0("page ", g, " of ", npages)),
      top = grid::textGrob(
        label = paste0(
          "Convergence of the omega parameter - ", NOmega, " species pair sample"),
        gp = grid::gpar(fontface = "bold", fontsize = 20)),
      nrow = NRC[1], ncol = NRC[2]) %>%
    ggplot2::ggsave(
      dpi = 600, device = "pdf", width = 18, height = 12,
      filename = file.path(Path_Convergence, "Convergence_Omega.pdf"))

  if (SavePlotData) {
    IASDT.R::CatTime("  >>  Save plot data")
    IASDT.R::SaveAs(
      InObj = OmegaTracePlots, OutObj = "Convergence_OmegaTracePlots",
      OutPath = file.path(
        Path_Convergence, "Convergence_OmegaTracePlots.RData"))
  }

  snow::stopCluster(c1)
  closeAllConnections()
  rm(Obj_Omega, OmegaDF, SelectedCombs, CI, OmegaTracePlots)
  invisible(gc())

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Beta ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("Beta")

  IASDT.R::CatTime("  >>  Prepare trace plot data")

  IASDT.R::CatTime("  >> >> Prepare working in parallel")
  c1 <- snow::makeSOCKcluster(NCores)
  future::plan(future::cluster, workers = c1, gc = TRUE)

  IASDT.R::CatTime("  >> >> Prepare confidence interval data")
  CI <- summary(Obj_Beta, quantiles = c(0.25, 0.75))$quantiles %>%
    as.data.frame() %>%
    tibble::as_tibble(rownames = "Var_Sp") %>%
    dplyr::rename(CI_25 = `25%`, CI_75 = `75%`)

  IASDT.R::CatTime("  >> >> Coda to tibble")
  BetaTracePlots <- Obj_Beta %>%
    IASDT.R::Coda_to_tibble(Type = "beta", EnvFile = EnvFile) %>%
    dplyr::left_join(CI, by = "Var_Sp")
  rm(CI)

  VarRanges <- BetaTracePlots %>%
    dplyr::arrange(Variable, IAS_ID) %>%
    dplyr::select(Var = Variable, DT) %>%
    dplyr::mutate(
      Range = purrr::map(
        .x = DT,
        .f = ~{
          .x %>%
            dplyr::pull(Value) %>%
            range()
        })) %>%
    dplyr::select(-DT) %>%
    tidyr::nest(data = c("Range")) %>%
    dplyr::mutate(
      Range = purrr::map(
        .x = data,
        .f = ~{
          .x %>%
            pull(Range) %>%
            as.vector() %>%
            range() %>%
            purrr::set_names(c("Min", "Max"))
        })) %>%
    dplyr::select(-data) %>%
    tidyr::unnest_wider("Range")

  IASDT.R::CatTime("  >> >> Export objects to cores")
  invisible(snow::clusterEvalQ(
    cl = c1, IASDT.R::LoadPackages(dplyr, coda, ggplot2, ggExtra, ggtext)))
  snow::clusterExport(
    cl = c1, list = c("BetaTracePlots", "Cols", "VarRanges"),
    envir = environment())

  IASDT.R::CatTime("  >> >> Prepare data in parallel")
  BetaTracePlots <- snow::parLapply(
    cl = c1, x = BetaTracePlots$Var_Sp,
    fun = function(x) {

      SubDT <- dplyr::filter(BetaTracePlots, Var_Sp == x)
      Species <- SubDT$Species
      IAS_ID <- SubDT$IAS_ID
      DT <- SubDT$DT[[1]]

      Variable <- SubDT$Variable
      VariableLabel <- dplyr::case_when(
        Variable == "RoadRailLog" ~ "Road + Rail intensity",
        Variable == "BiasLog" ~ "Sampling intensity",
        .default = Variable)

      ValRange <- range(DT$Value)

      ValRangeByVar <- VarRanges %>%
        dplyr::filter(Var == Variable) %>%
        dplyr::select(-Var) %>%
        unlist()

      CI0 <- round(c(SubDT$CI_25, SubDT$CI_75), 2)
      CI1 <- CI0 %>%
        paste0(collapse = " - ") %>%
        paste0("<i>50% credible interval:</i> ", .) %>%
        data.frame(x = -Inf, y = -Inf, label = .)

      PanelTitle <- data.frame(
        x = Inf, y = Inf, label = paste0("<b><i>", Species, "</i></b>"))
      PanelTitle2 <- IASDT.R::GetSpeciesName(
        EnvFile = EnvFile, SpID = IAS_ID) %>%
        dplyr::select(Class, Order, Family) %>%
        unlist() %>%
        paste0(collapse = " | ") %>%
        paste0("<b>", ., "</b>") %>%
        paste0("<br>", IAS_ID) %>%
        data.frame(x = -Inf, y = Inf, label = .)
      PanelTitle3 <- data.frame(x = Inf, y = -Inf, label = VariableLabel)

      Plot <- DT %>%
        ggplot2::ggplot(mapping = ggplot2::aes(
          x = Iter, y = Value, color = factor(Chain))) +
        ggplot2::geom_line(linewidth = 0.15, alpha = 0.6) +
        ggplot2::geom_smooth(
          method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.8) +
        ggplot2::geom_point(alpha = 0) +
        ggplot2::geom_hline(
          yintercept = CI0, linetype = "dashed", color = "black",
          linewidth = 1) +
        ggplot2::scale_color_manual(values = Cols) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(limits = ValRange) +
        ggtext::geom_richtext(
          mapping = ggplot2::aes(x = x, y = y, label = label), data = CI1,
          inherit.aes = FALSE, hjust = 0, vjust = 0, lineheight = 0, fill = NA,
          label.color = NA) +
        ggtext::geom_richtext(
          mapping = ggplot2::aes(x = x, y = y, label = label),
          data = PanelTitle, inherit.aes = FALSE, colour = "blue",
          hjust = 1, vjust = 1, lineheight = 0, fill = NA, label.color = NA) +
        ggtext::geom_richtext(
          mapping = ggplot2::aes(x = x, y = y, label = label),
          data = PanelTitle2, inherit.aes = FALSE,
          hjust = 0, vjust = 1, lineheight = 0, fill = NA, label.color = NA) +
        ggtext::geom_richtext(
          mapping = ggplot2::aes(x = x, y = y, label = label), size = 6,
          data = PanelTitle3, inherit.aes = FALSE, colour = "darkgrey",
          hjust = 1, vjust = 0, lineheight = 0, fill = NA, label.color = NA) +
        ggplot2::theme_bw() +
        ggplot2::xlab(NULL) +
        ggplot2::ylab(NULL) +
        ggplot2::theme(
          legend.position = "none",
          axis.text = ggplot2::element_text(size = 10))

      suppressMessages({
        Plot2 <- Plot +
          ggplot2::scale_y_continuous(limits = ValRangeByVar)
      })

      if (dplyr::between(0, ValRangeByVar[1], ValRangeByVar[2])) {
        Plot2 <- Plot2 +
          ggplot2::geom_hline(
            yintercept = 0, linetype = "solid", color = "green",
            linewidth = 0.4)
      }

      if (dplyr::between(0, ValRange[1], ValRange[2])) {
        Plot <- Plot +
          ggplot2::geom_hline(
            yintercept = 0, linetype = "solid", color = "green",
            linewidth = 0.4)
      }

      Plot <- ggExtra::ggMarginal(
        p = Plot, type = "density", margins = "y", size = 5, color = "steelblue4")
      Plot2 <- ggExtra::ggMarginal(
        p = Plot2, type = "density", margins = "y", size = 5, color = "steelblue4")

      invisible(gc())
      return(tibble::tibble(Var_Sp = x, Plot = list(Plot), PlotFixedY = list(Plot2)))
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::right_join(BetaTracePlots, by = "Var_Sp")

  snow::stopCluster(c1)
  closeAllConnections()
  invisible(gc())

  IASDT.R::CatTime("  >>  Save trace plot data")
  if (SavePlotData) {
    IASDT.R::SaveAs(
      InObj = BetaTracePlots, OutObj = "Convergence_BetaTracePlots",
      OutPath = file.path(
        Path_Convergence, "Convergence_BetaTracePlots.RData"))
  }

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  ## Beta - by variable - Free y ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("  >>  Trace plots, grouped by variables - Free Y axis")

  IASDT.R::CatTime("  >>  >>  Preparing data")
  BetaTracePlots_ByVar <- BetaTracePlots %>%
    dplyr::arrange(Variable, IAS_ID) %>%
    dplyr::select(Variable, Plot) %>%
    tidyr::nest(data = c("Plot"))

  IASDT.R::CatTime("  >>  >>  Preparing working in parallel")
  c1 <- snow::makeSOCKcluster(NCores)
  future::plan(future::cluster, workers = c1, gc = TRUE)

  invisible(snow::clusterEvalQ(
    cl = c1, IASDT.R::LoadPackages(dplyr, coda, ggplot2, ggExtra, ggtext)))
  snow::clusterExport(
    cl = c1, list = c("BetaTracePlots_ByVar", "Path_Convergence", "NRC", "Cols"),
    envir = environment())

  IASDT.R::CatTime("  >>  >>  save plots in parallel")
  BetaTracePlots_ByVar0 <- snow::parLapply(
    cl = c1, x = BetaTracePlots_ByVar$Variable,
    fun = function(x) {

      VarName <- dplyr::case_when(
        x == "RoadRailLog" ~ "Road + Rail intensity",
        x == "BiasLog" ~ "Sampling intensity",
        .default = x)

      Path_Beta <- file.path(
        Path_Convergence, paste0("Convergence_Beta_", x, "_FreeY.pdf"))

      Plot <- dplyr::filter(BetaTracePlots_ByVar, Variable == x) %>%
        dplyr::pull(data) %>%
        magrittr::extract2(1) %>%
        dplyr::pull(1) %>%
        gridExtra::marrangeGrob(
          bottom = bquote(paste0("page ", g, " of ", npages)),
          top = grid::textGrob(
            label = paste0("Convergence of the beta parameter - ", VarName),
            gp = grid::gpar(fontface = "bold", fontsize = 20)),
          nrow = NRC[1], ncol = NRC[2])

      ggplot2::ggsave(
        plot = Plot, dpi = 600, device = "pdf", width = 18, height = 12,
        filename = Path_Beta)

      return(Path_Beta)
    })

  snow::stopCluster(c1)
  closeAllConnections()
  rm(BetaTracePlots_ByVar0, BetaTracePlots_ByVar)
  invisible(gc())

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  ## Beta - by variable - Fixed y ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("  >>  Trace plots, grouped by variables - Fixed Y axis")

  IASDT.R::CatTime("  >>  >>  Preparing data")
  BetaTracePlots_ByVar <- BetaTracePlots %>%
    dplyr::arrange(Variable, IAS_ID) %>%
    dplyr::select(Variable, PlotFixedY) %>%
    tidyr::nest(data = c("PlotFixedY"))

  IASDT.R::CatTime("  >>  >>  Preparing working in parallel")
  c1 <- snow::makeSOCKcluster(NCores)
  future::plan(future::cluster, workers = c1, gc = TRUE)

  invisible(snow::clusterEvalQ(
    cl = c1, IASDT.R::LoadPackages(dplyr, coda, ggplot2, ggExtra, ggtext)))
  snow::clusterExport(
    cl = c1, list = c("BetaTracePlots_ByVar", "Path_Convergence", "NRC", "Cols"),
    envir = environment())

  IASDT.R::CatTime("  >>  >>  save plots in parallel")
  BetaTracePlots_ByVar0 <- snow::parLapply(
    cl = c1, x = BetaTracePlots_ByVar$Variable,
    fun = function(x) {

      VarName <- dplyr::case_when(
        x == "RoadRailLog" ~ "Road + Rail intensity",
        x == "BiasLog" ~ "Sampling intensity",
        .default = x)

      Path_Beta <- file.path(
        Path_Convergence, paste0("Convergence_Beta_", x, "_FixedY.pdf"))

      Plot <- dplyr::filter(BetaTracePlots_ByVar, Variable == x) %>%
        dplyr::pull(data) %>%
        magrittr::extract2(1) %>%
        dplyr::pull(1) %>%
        gridExtra::marrangeGrob(
          bottom = bquote(paste0("page ", g, " of ", npages)),
          top = grid::textGrob(
            label = paste0("Convergence of the beta parameter - ", VarName, " - Fixed y-axis"),
            gp = grid::gpar(fontface = "bold", fontsize = 20)),
          nrow = NRC[1], ncol = NRC[2])

      ggplot2::ggsave(
        plot = Plot, dpi = 600, device = "pdf", width = 18, height = 12,
        filename = Path_Beta)

      return(Path_Beta)
    })

  snow::stopCluster(c1)
  closeAllConnections()
  rm(BetaTracePlots_ByVar0, BetaTracePlots_ByVar)
  invisible(gc())

  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  # Beta - by species ------
  # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::CatTime("  >>  Trace plots, grouped by species")

  IASDT.R::CatTime("  >>  >>  Preparing data")
  Order <- c("Intercept", "bio4", "bio6", "bio8", "bio12", "bio15",
             "bio18", "RoadRailLog", "BiasLog")
  BetaTracePlots_BySp <- BetaTracePlots %>%
    dplyr::slice(order(factor(Variable, levels = Order))) %>%
    dplyr::select(Species, IAS_ID, Plot) %>%
    tidyr::nest(data = -c("Species", "IAS_ID"))

  IASDT.R::CatTime("  >>  >>  Preparing working in parallel")
  c1 <- snow::makeSOCKcluster(NCores)
  future::plan(future::cluster, workers = c1, gc = TRUE)

  invisible(snow::clusterEvalQ(
    cl = c1, IASDT.R::LoadPackages(dplyr, coda, ggplot2, ggExtra, ggtext)))
  snow::clusterExport(
    cl = c1, list = c("BetaTracePlots_BySp", "Path_Convergence_BySp"),
    envir = environment())

  IASDT.R::CatTime("  >>  >>  save plots in parallel")
  BetaTracePlots_BySp0 <- snow::parLapply(
    cl = c1, x = BetaTracePlots_BySp$Species,
    fun = function(x) {

      VarName <- dplyr::case_when(
        x == "RoadRailLog" ~ "Road + Rail intensity",
        x == "BiasLog" ~ "Sampling intensity",
        .default = x)

      DT <- dplyr::filter(BetaTracePlots_BySp, Species == x)

      Path_Beta <- file.path(
        Path_Convergence_BySp,
        paste0("Convergence_Beta_", DT$IAS_ID, "_",
               IASDT.R::ReplaceSpace(x), ".pdf"))

      Plot <- DT %>%
        dplyr::pull(data) %>%
        magrittr::extract2(1) %>%
        dplyr::pull(1) %>%
        gridExtra::marrangeGrob(
          top = grid::textGrob(
            label = paste0("Convergence of the beta parameter - ", x),
            gp = grid::gpar(fontface = "bold", fontsize = 20)),
          nrow = 3, ncol = 3)

      ggplot2::ggsave(
        plot = Plot, dpi = 600, device = "pdf", width = 18, height = 12,
        filename = Path_Beta, create.dir = TRUE)

      return(Path_Beta)
    })
  rm(BetaTracePlots_BySp0)

  return(invisible(NULL))

}
