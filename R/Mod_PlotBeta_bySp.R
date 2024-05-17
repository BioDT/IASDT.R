# |---------------------------------------------------| #
# PlotBeta_bySp ----
# |---------------------------------------------------| #

#' PlotBeta_bySp
#'
#' PlotBeta_bySp
#' @param Beta Beta
#' @param NRC NRC
#' @param PdfPath PdfPath
#' @param FilePrefix FilePrefix
#' @param Verbose Verbose
#' @name PlotBeta_bySp
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

PlotBeta_bySp <- function(
    Beta, NRC = c(3, 2), PdfPath = getwd(), FilePrefix = "Posterior_Sp",
    Verbose = TRUE) {

  fs::dir_create(PdfPath)

  PostID <- dimnames(Beta[[1]])[[2]] %>%
    stringr::str_split(",", simplify = TRUE) %>%
    as.data.frame() %>%
    stats::setNames(c("Var", "Sp")) %>%
    tibble::tibble() %>%
    dplyr::mutate(
      Var = stringr::str_remove_all(Var, "B\\[| \\(.+\\)|\\(|\\)"),
      Sp = stringr::str_trim(Sp),
      Sp = stringr::str_remove_all(Sp, "^Sp_| \\(S.+\\)\\]|\\]"),
      ID = seq_len(dplyr::n())) %>%
    dplyr::group_by(Sp) %>%
    dplyr::group_split() %>%
    stats::setNames(purrr::map(., ~paste0("Sp_", .x$Sp[1])))

  PostID %>%
    purrr::map(
      .f = ~{
        SpNum <- .x$Sp[1]
        SpName <- TaxaList %>%
          dplyr::filter(IAS_ID == as.numeric(SpNum)) %>%
          dplyr::pull(Species_name)

        if (Verbose) cat(paste0(as.numeric(SpNum), " (", SpName, ")\n"))

        paste0(FilePrefix, SpNum, ".pdf") %>%
          file.path(PdfPath, .) %>%
          grDevices::pdf(height = 12, width = 9, pointsize = 12)

        par(mfrow = NRC, oma = c(0.05, 0.05, 0.05, 0.05), mar = c(4, 3, 3, 1))

        seq_len(nrow(.x)) %>%
          purrr::walk(
            .f = function(Y) {

              Title <- paste0(as.numeric(SpNum), " (", SpName, ") / ", .x$Var[Y])
              PlotIndex <- .x$ID[Y]

              MCMC_List <- coda::mcmc.list(
                Beta[[1]][, PlotIndex], Beta[[2]][, PlotIndex])

              MCMC_List %>%
                coda::traceplot(
                  smooth = TRUE, lwd = 0.4, cex.axis = 1.25, cex.lab = 1.25,
                  col = c("black", "blue"))
              mtext(Title, cex = 1.1, adj = 0, col = "black", font = 2, padj = 0)

              MCMC_List %>%
                coda::densplot(
                  show.obs = TRUE, lwd = 0.75, main = NA, cex.axis = 1.25,
                  cex.lab = 1.25)
            })
        grDevices::dev.off()
      }) %>%
    invisible()
}
