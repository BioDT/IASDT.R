# |---------------------------------------------------| #
# Range2NewVal ----
# |---------------------------------------------------| #
#
#' Change values of `vector`/`data.frame`/`raster` within/outside a certain range to another value
#'
#' Change values of `vector`/`data.frame`/`raster` within/outside a certain range to another value
#'
#' @name Range2NewVal
#' @author Ahmed El-Gabbas
#' @param x vector / data.frame / raster object.
#' @param Between Range of values to be changed or kept.
#' @param MoreThan The value above which the original values will be changed. Only applied if `Between` is not used.
#' @param LessThan The value below which the original values will be changed. Only applied if `Between` and `MoreThan` are not used.
#' @param NewVal The new value to be assigned. Default: `NA`.
#' @param InvertSelection Should the selection be inverted? Only valid if `Between` used. If `InvertSelection = TRUE`, then the values not in the range of `Between` argument will be changed. Default: `FALSE`.
#'
#' @export
#' @examples
#' # Vector
#'
#' Range2NewVal(x =  1:10, Between = c(5, 8), NewVal = NA)
#'
#' Range2NewVal(x =  1:10, Between = c(5, 8), NewVal = NA, InvertSelection = TRUE)
#'
#' # ---------------------------------------------
#'
#' # tibble
#'
#' iris %>%
#'  tibble::as_tibble() %>%
#'  dplyr::slice_head(n = 50) %>%
#'  dplyr::select(-Sepal.Length, -Petal.Length, -Petal.Width) %>%
#'  dplyr::mutate(
#'    Sepal.Width.New = Range2NewVal(
#'        x = Sepal.Width, Between = c(3, 3.5), NewVal = NA, InvertSelection = FALSE),
#'    Sepal.Width.Rev = Range2NewVal(
#'        x = Sepal.Width, Between = c(3, 3.5), NewVal = NA, InvertSelection = TRUE)) %>%
#'  dplyr::arrange(-Sepal.Width) %>%
#'  print(n = 50)
#'
#' # ---------------------------------------------
#'
#' # raster
#'
#' LoadPackages(raster)
#'
#' RRR <- system.file("external/test.grd", package = "raster") %>%
#'     raster::raster()
#'
#' RRR2 <- Range2NewVal(x = RRR, LessThan = 500, NewVal = NA)
#' RRR3 <- Range2NewVal(x = RRR, MoreThan = 500, NewVal = NA)
#' par(mar = c(0.5, 0.5, 3, 3))
#' plot(raster::stack(RRR, RRR2, RRR3), nr = 1, main = c("Original", "<500 to NA", ">500 to NA"))
#'
#' RRR2 <- Range2NewVal(x = RRR, Between = c(1000, 1800), NewVal = 1800, InvertSelection = FALSE)
#' RRR3 <- Range2NewVal(x = RRR, Between = c(1000, 1800), NewVal = 1800, InvertSelection = TRUE)
#' plot(stack(RRR>=1000, RRR2, RRR3), nr = 1, main = c(">1000 ?", "<500 to NA", ">500 to NA"))

Range2NewVal <- function(
    x = NULL, Between = NULL,
    MoreThan = NULL, LessThan = NULL,
    NewVal = NA, InvertSelection = FALSE) {

  if (!is.null(Between)) {
    if (length(Between) != 2) stop()
    Min <- Between[1]
    Max <- Between[2]
    if (inherits(x, "RasterLayer")) {
      X1 <- X2 <- x
      X1[X1 >= Max] <- NA
      X1[!is.na(X1)] <- 1
      X2[X2 <= Min] <- NA
      X2[!is.na(X2)] <- 1
      X3 <- sum(X1, X2, na.rm = TRUE)

      if (InvertSelection) {
        x[X3 == 1] <- NewVal
      } else {
        x[X3 == 2] <- NewVal
      }
    } else {

      if (is.null(Max)) Max <- max(x)
      if (is.null(Min)) Min <- min(x)
      if (Max <= Min) stop()

      if (InvertSelection) {
        x[!(x >= Min & x <= Max)] <- NewVal
      } else {
        x[x >= Min & x <= Max] <- NewVal
      }
    }
  } else {

    if (!is.null(MoreThan)) {
      x[x > MoreThan] <- NewVal
    } else {
      if (!is.null(LessThan)) x[x < LessThan] <- NewVal
    }
  }
  x
}
