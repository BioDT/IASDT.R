# |---------------------------------------------------| #
# PlotBeta_byVar ----
# |---------------------------------------------------| #

#' PlotBeta_byVar
#'
#' PlotBeta_byVar
#' @param Beta Beta
#' @param NRC NRC
#' @param PdfPath PdfPath
#' @param FilePrefix FilePrefix
#' @param Verbose Verbose
#' @param VerboseSp VerboseSp
#' @name PlotBeta_byVar
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

PlotBeta_byVar <- function(
    Beta, NRC = c(3, 2), PdfPath = getwd(), FilePrefix = "Posterior_Sp",
    Verbose = TRUE, VerboseSp = FALSE) {

  fs::dir_create(PdfPath)

  PostID <- dimnames(Beta[[1]])[[2]] %>%
    stringr::str_split(",", simplify = TRUE) %>%
    as.data.frame() %>%
    stats::setNames(c("Var", "Sp")) %>%
    tibble::tibble() %>%
    # dplyr::arrange(Var) %>%
    dplyr::mutate(
      Var = stringr::str_remove_all(Var, "B\\[| \\(.+\\)|\\(|\\)"),
      Sp = stringr::str_trim(Sp),
      Sp = stringr::str_remove_all(Sp, "^Sp_| \\(S.+\\)\\]|\\]"),
      ID = seq_len(dplyr::n())) %>%
    dplyr::group_by(Var) %>%
    dplyr::group_split() %>%
    stats::setNames(purrr::map(., ~paste0("Var_", .x$Var[1])))

  load("Data/IAS_Final_Data/Sp_PA_Summary_DF.RData")

  PostID %>%
    purrr::map(
      .f = ~{
        if (Verbose) cat(paste0(.x$Var[1], "\n"))

        # NRC2 <- NRC[1]*NRC[2]
        NPlots <- nrow(.x)
        TitleIDs <- .x$ID[seq(1, NPlots, by = NRC[1])]

        paste0(FilePrefix, .x$Var[1], ".pdf") %>%
          file.path(PdfPath, .) %>%
          grDevices::pdf(height = 12, width = 9, pointsize = 12)

        graphics::par(mfrow = NRC, oma = c(0.05, 0.05, 0.05, 0.05), mar = c(4, 3, 3, 1))

        seq_len(nrow(.x)) %>%
          purrr::walk(
            .f = function(Y) {

              PlotIndex <- .x$ID[Y]
              SpNum <- .x$Sp[Y]
              SpName <- TaxaList %>%
                dplyr::filter(IAS_ID == as.numeric(SpNum)) %>%
                dplyr::pull(Species_name)

              if (VerboseSp) {
                cat(paste0("  >>>  " , .x$Sp[Y], " (", SpName, ")\n"))
              }

              Title <- paste0(as.numeric(SpNum), " (", SpName, ")")

              MCMC_List <- coda::mcmc.list(
                Beta[[1]][, PlotIndex], Beta[[2]][, PlotIndex],
                Beta[[3]][, PlotIndex], Beta[[4]][, PlotIndex])

              MCMC_List %>%
                coda::traceplot(
                  smooth = TRUE, lwd = 0.2, cex.axis = 1.25, cex.lab = 1.25,
                  col = c("black", "blue", "red", "darkgreen"), lty = 1)
              mtext(Title, cex = 1.1, adj = 0, col = "black", font = 2, padj = 0)

              MCMC_List %>%
                coda::densplot(
                  show.obs = TRUE, lwd = 0.75, main = NA, cex.axis = 1.25,
                  cex.lab = 1.25)

              if( PlotIndex %in% TitleIDs) {
                mtext(
                  paste0(.x$Var[1], "  "), cex = 1.75, adj = 1, col = "blue",
                  outer = TRUE, line = -0.75, font = 2, padj = 1)
              }

            })
        grDevices::dev.off()
      }) %>%
    invisible()
}
