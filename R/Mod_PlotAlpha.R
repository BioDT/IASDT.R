## |------------------------------------------------------------------------| #
# PlotAlpha ----
## |------------------------------------------------------------------------| #

#' Plot convergence traceplots for the alpha parameter
#'
#' Plot convergence traceplots for the alpha parameter
#'
#' @param Post Coda object or path to it
#' @param Model fitted model object or path to it
#' @param Title String. Plotting title
#' @param NRC Vector for the number of rows and columns per plot page
#' @param AddFooter Add footer to the plot for the page number
#' @param AddTitle Add main title to the plot
#' @param Cols Colours for lines for each chain
#' @name PlotAlpha
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

PlotAlpha <- function(
    Post = NULL, Model = NULL, Title = NULL, NRC = NULL, AddFooter = TRUE,
    AddTitle = TRUE, Cols = c("red", "blue", "darkgreen", "darkgrey")) {

  AllArgs <- ls()
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)
  IASDT.R::CheckArgs(AllArgs = AllArgs, Type = "character", Args = "Title")

  x <- Factor <- NULL

  # Load coda object
  if (inherits(Post, "character")) {
    Post <- IASDT.R::LoadAs(Post)
  }
  Post <- Post$Alpha[[1]]

  # Load model object
  if (inherits(Model, "character")) {
    Model <- IASDT.R::LoadAs(Model)
  }

  SampleSize <- Model$samples
  NChains <- length(Model$postList)
  rm(Model)
  invisible(gc())

  NLV <- ncol(Post[[1]])

  ## Gelman convergence diagnostic
  Gelman <- Post %>%
    coda::gelman.diag(multivariate = FALSE) %>%
    magrittr::extract2("psrf") %>%
    as.data.frame() %>%
    dplyr::pull(1) %>%
    round(2)

  ## Effective sample size
  ESS <- Post %>%
    coda::effectiveSize() %>%
    round() %>%
    magrittr::divide_by(NChains)

  ## quantiles
  CI <- Post %>%
    summary(quantiles = c(0.25, 0.75)) %>%
    magrittr::extract2("quantiles") %>%
    matrix(ncol = 2) %>%
    round(3)

  AlphaDF <- Post %>%
    IASDT.R::Coda_to_tibble(Type = "alpha") %>%
    dplyr::mutate(
      Factor2 = purrr::map_int(
        .x = Factor,
        .f = ~{
          as.character(.x) %>%
            stringr::str_remove("factor") %>%
            as.integer()
        }))

  if (is.null(NRC)) {
    NRC <- dplyr::case_when(
      NLV == 1 ~ c(1, 1),
      NLV == 2 ~ c(1, 2),
      NLV == 3 ~ c(1, 3),
      NLV == 4 ~ c(2, 2),
      .default = c(2, 3))
  }


  Plots <- purrr::map(
    .x = seq_len(NLV),
    .f = ~{

      ESS0 <- paste0(
        "<b><i>Mean effective sample size:</i></b> ", ESS[.x],
        " / ", SampleSize, " samples")

      CI0 <- CI[.x, ] %>%
        round(2) %>%
        paste0(collapse = "-") %>%
        paste0("<b><i>50% credible interval:</i></b> ", .)

      ESS_CI <- data.frame(
        x = -Inf, y = -Inf, label = paste0(ESS0, "<br>", CI0))

      Gelman0 <- paste0("<b><i>Gelman convergence diagnostic:</i></b> ", Gelman[.x])
      Title2 <- data.frame(x = Inf, y = Inf, label = Gelman0)
      Title3 <- data.frame(x = -Inf, y = Inf, label = paste0("Factor", .x))

      Plot <- AlphaDF %>%
        dplyr::filter(Factor2 == .x) %>%
        ggplot2::ggplot(
          mapping = ggplot2::aes(x = Iter, y = Value, color = Chain)) +
        ggplot2::geom_line(linewidth = 0.15, alpha = 0.6) +
        ggplot2::geom_smooth(
          method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.8) +
        ggplot2::geom_point(alpha = 0) +
        ggplot2::geom_hline(
          yintercept = CI[.x, ], linetype = "dashed", color = "black",
          linewidth = 1) +
        ggplot2::scale_color_manual(values = Cols) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::theme_bw() +
        ggplot2::xlab(NULL) +
        ggplot2::ylab(NULL) +
        ggtext::geom_richtext(
          mapping = ggplot2::aes(x = x, y = y, label = label),
          data = Title2, inherit.aes = FALSE, size = 4, hjust = 1, vjust = 1,
          lineheight = 0, fill = NA, label.color = NA) +
        ggtext::geom_richtext(
          mapping = ggplot2::aes(x = x, y = y, label = label),
          data = Title3, inherit.aes = FALSE, size = 4,
          hjust = -0.1, vjust = 1, color = "blue",
          lineheight  = 0, fill = NA, label.color = NA) +
        ggtext::geom_richtext(
          mapping = ggplot2::aes(x = x, y = y, label = label),
          data = ESS_CI, inherit.aes = FALSE, size = 4,
          hjust = 0, vjust = 0, lineheight  = 0, fill = NA, label.color = NA) +
        ggplot2::theme(
          legend.position = "none", axis.text = ggplot2::element_text(size = 12),
          axis.title = ggplot2::element_blank())

      Plot <- ggExtra::ggMarginal(
        p = Plot, type = "density", margins = "y", size = 4, color = "steelblue4")
      return(Plot)
    }) %>%

    if (AddTitle) {
      gridExtra::marrangeGrob(
        bottom = dplyr::if_else(
          condition = AddFooter,
          true = bquote(paste0("page ", g, " of ", npages)), false = NULL),
        top = grid::textGrob(
          label = Title, gp = grid::gpar(fontface = "bold", fontsize = 20)),
        nrow = NRC[1], ncol = NRC[2])
    } else {
      gridExtra::marrangeGrob(
        bottom = dplyr::if_else(
          condition = AddFooter,
          true = bquote(paste0("page ", g, " of ", npages)), false = NULL),
        top = NULL, nrow = NRC[1], ncol = NRC[2])
    }

  return(Plots)
}
