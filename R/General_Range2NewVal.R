## |------------------------------------------------------------------------| #
# Range2NewVal ----
## |------------------------------------------------------------------------| #
#
#' Changes values within a specified range, `≥`, or `≤`  specific value to a new
#' value in a vector, data.frame, or raster.
#'
#' This function modifies values in the input object `x` based on the specified
#' conditions. It can operate on vectors, data.frames, or RasterLayer objects.
#' The function allows for changing values within a specified range (`Between`),
#' `≥` a specified value (`MoreThan`), or `≤` a specified value (`LessThan`) to
#' a new value (`NewVal`). An option to invert the selection is also available
#' for ranges.
#' @name Range2NewVal
#' @author Ahmed El-Gabbas
#' @param x A vector, data.frame, or RasterLayer object whose values are to be
#'   modified.
#' @param Between A numeric vector of length 2 specifying the range of values to
#'   be changed or kept. If specified, `MoreThan` and `LessThan` are ignored.
#' @param MoreThan,LessThan A numeric value specifying the threshold above/below
#'   which values in `x` will be changed to `NewVal`. Only applied if `Between`
#'   is not specified.
#' @param NewVal The new value to assign to the selected elements in `x`.
#' @param InvertSelection A logical value indicating whether to invert the
#'   selection specified by `Between`. If `TRUE`, values outside the specified
#'   range are changed to `NewVal`. Default is `FALSE`.
#' @return The modified object `x` with values changed according to the
#'   specified conditions.
#' @export
#' @examples
#' # Vector
#'
#' Range2NewVal(x =  1:10, Between = c(5, 8), NewVal = NA)
#'
#'
#' Range2NewVal(
#' x =  1:10, Between = c(5, 8), NewVal = NA, InvertSelection = TRUE)
#'
#' Range2NewVal(x =  1:10, Between = c(5, 8), NewVal = NA, MoreThan = 4)
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
#'        x = Sepal.Width, Between = c(3, 3.5), NewVal = NA,
#'         InvertSelection = FALSE),
#'    Sepal.Width.Rev = Range2NewVal(
#'        x = Sepal.Width, Between = c(3, 3.5), NewVal = NA,
#'         InvertSelection = TRUE)) %>%
#'  dplyr::arrange(-Sepal.Width) %>%
#'  print(n = 50)
#'
#' # ---------------------------------------------
#'
#' # raster
#'
#' library(raster)
#'
#' RRR <- raster::raster(system.file("external/test.grd", package = "raster"))
#'
#' RRR2 <- Range2NewVal(x = RRR, LessThan = 500, NewVal = NA)
#' RRR3 <- Range2NewVal(x = RRR, MoreThan = 500, NewVal = NA)
#' par(mar = c(0.5, 0.5, 3, 3))
#' plot(
#'    raster::stack(RRR, RRR2, RRR3), nr = 1,
#'    main = c("Original", "<500 to NA", ">500 to NA"))
#'
#' RRR2 <- Range2NewVal(
#'    x = RRR, Between = c(1000, 1800), NewVal = 1800, InvertSelection = FALSE)
#' RRR3 <- Range2NewVal(
#'    x = RRR, Between = c(1000, 1800), NewVal = 1800, InvertSelection = TRUE)
#' plot(
#'    raster::stack(RRR>=1000, RRR2, RRR3), nr = 1,
#'    main = c(">1000 ?", "<500 to NA", ">500 to NA"))

Range2NewVal <- function(
    x, Between = NULL, MoreThan = NULL, LessThan = NULL,
    NewVal, InvertSelection = FALSE) {

  if (is.null(x) || is.null(NewVal)) {
    stop("x and NewVal cannot be NULL", .call = FALSE)
  }

  if (all(is.null(MoreThan), is.null(LessThan), is.null(Between))) {
      stop(
        "At least one of MoreThan, LessThan, and Between should be not NULL",
        .call = FALSE)
  }

  if (!is.null(Between)) {

    if (length(Between) != 2) {
      stop(
        "Between should have exactly two values: a minimum and a maximum.", 
        .call = FALSE)
    }

    Min <- Between[1]
    Max <- Between[2]

    if (Max <= Min) {
      stop("Max must be greater than Min.", .call = FALSE)
    }


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

      if (InvertSelection) {
        x[!(x >= Min & x <= Max)] <- NewVal
      } else {
        x[x >= Min & x <= Max] <- NewVal
      }
    }
  } else {

    if (!is.null(MoreThan)) {
      x[x >= MoreThan] <- NewVal
    }

    if (!is.null(LessThan)) {
      x[x <= LessThan] <- NewVal
    }
  }
  return(x)
}
