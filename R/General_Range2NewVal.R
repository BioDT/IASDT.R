## |------------------------------------------------------------------------| #
# Range2NewVal ----
## |------------------------------------------------------------------------| #

#' Changes values within a specified range, or greater than or less than a
#' specific value to a new value in a vector, data.frame, or raster
#'
#' This function modifies values in the input object `x` based on the specified
#' conditions. It can operate on vectors, data.frames, or RasterLayer objects.
#' The function allows for changing values within a specified range (`Between`),
#' greater than or equals to (`MoreThan`) or  less than or equals to
#' (`LessThan`) a specified value to a new value (`NewVal`). An option to invert
#' the selection is also available for ranges.
#' @name Range2NewVal
#' @author Ahmed El-Gabbas
#' @param x A numeric `vector`, `data.frame`, `RasterLayer`, or `SpatRaster`
#'   object whose values are to be modified.
#' @param Between Numeric. A numeric vector of length 2 specifying the range of
#'   values to be changed or kept. If specified, `MoreThan` and `LessThan` are
#'   ignored.
#' @param MoreThan,LessThan Numeric. Threshold larger than or equal to/less than
#'   or equal to which values in `x` will be changed to `NewVal`. Only applied
#'   if `Between` is not specified.
#' @param NewVal The new value to assign to the selected elements in `x`.
#' @param Invert Logical. Whether to invert the selection specified by
#'   `Between`. If `TRUE`, values outside the specified range are changed to
#'   `NewVal`. Default is `FALSE`.
#' @return The modified object `x` with values changed according to the
#'   specified conditions.
#' @export
#' @examples
#' library(raster)
#' library(terra)
#' par(mar = c(0.5, 0.5, 1, 2.5), oma = c(0.5, 0.5, 0.5, 1))
#'
#' # ---------------------------------------------
#'
#' # Vector
#'
#' (VV <- seq_len(10))
#'
#' Range2NewVal(x = VV, Between = c(5, 8), NewVal = NA)
#'
#' Range2NewVal(x = VV, Between = c(5, 8), NewVal = NA, Invert = TRUE)
#'
#' # MoreThan is ignored as Between is specified
#' Range2NewVal(x = VV, Between = c(5, 8), NewVal = NA, MoreThan = 4)
#'
#' Range2NewVal(x = VV, NewVal = NA, MoreThan = 4)
#'
#' Range2NewVal(x = VV, NewVal = NA, LessThan = 4)
#'
#' # ---------------------------------------------
#'
#' # tibble
#'
#' iris2 <- iris %>%
#'   tibble::as_tibble() %>%
#'   dplyr::slice_head(n = 50) %>%
#'   dplyr::select(-Sepal.Length, -Petal.Length, -Petal.Width) %>%
#'   dplyr::arrange(-Sepal.Width)
#'
#' iris2 %>%
#'  dplyr::mutate(
#'    Sepal.Width.New = Range2NewVal(
#'       x = Sepal.Width, Between = c(3, 3.5), NewVal = NA, Invert = FALSE),
#'    Sepal.Width.Rev = Range2NewVal(
#'       x = Sepal.Width, Between = c(3, 3.5), NewVal = NA, Invert = TRUE)) %>%
#'  print(n = 50)
#'
#' # ---------------------------------------------
#'
#' # RasterLayer / SpatRaster
#'
#' grd_file <- system.file("external/test.grd", package = "raster")
#' R_raster <- raster::raster(grd_file)
#' R_terra <- terra::rast(grd_file)
#'
#' # Convert values less than 500 to NA
#' R_raster2 <- Range2NewVal(x = R_raster, LessThan = 500, NewVal = NA)
#' plot(
#'    raster::stack(R_raster, R_raster2), nr = 1,
#'    main = c("\nOriginal", "\n<500 to NA"),
#'    box = FALSE, axes = FALSE, legend.width = 2, colNA = "lightgrey",
#'    xaxs = "i", yaxs = "i")
#'
#' R_terra2 <- Range2NewVal(x = R_terra, LessThan = 500, NewVal = NA)
#' plot(
#'    c(R_terra, R_terra2), nr = 1, main = c("\nOriginal", "\n<500 to NA"),
#'    box = FALSE, axes = FALSE, colNA = "lightgrey", xaxs = "i", yaxs = "i")
#'
#'
#' # Convert values greater than 700 to NA
#' R_raster2 <- Range2NewVal(x = R_raster, MoreThan = 700, NewVal = NA)
#' plot(
#'    raster::stack(R_raster, R_raster2), nr = 1,
#'    main = c("\nOriginal", "\n>700 to NA"),
#'    box = FALSE, axes = FALSE, legend.width = 2, colNA = "lightgrey",
#'    xaxs = "i", yaxs = "i")
#'
#' R_terra2 <- Range2NewVal(x = R_terra, MoreThan = 700, NewVal = NA)
#' plot(
#'    c(R_terra, R_terra2), nr = 1, main = c("\nOriginal", "\n>700 to NA"),
#'    box = FALSE, axes = FALSE, colNA = "lightgrey", xaxs = "i", yaxs = "i")

Range2NewVal <- function(
    x = NULL, Between = NULL, MoreThan = NULL, LessThan = NULL,
    NewVal = NULL, Invert = FALSE) {

  if (is.null(x) || is.null(NewVal)) {
    stop("x and NewVal cannot be NULL", call. = FALSE)
  }

  if (all(is.null(MoreThan), is.null(LessThan), is.null(Between))) {
      stop(
        "At least one of MoreThan, LessThan, and Between should be not NULL",
        call. = FALSE)
  }

  if (!is.null(Between)) {

    if (length(Between) != 2) {
      stop(
        "Between should have exactly two values: a minimum and a maximum.",
        call. = FALSE)
    }

    Min <- Between[1]
    Max <- Between[2]

    if (Max <= Min) {
      stop("Max must be greater than Min.", call. = FALSE)
    }


    if (inherits(x, "RasterLayer")) {
      X1 <- X2 <- x
      X1[X1 >= Max] <- NA
      X1[!is.na(X1)] <- 1
      X2[X2 <= Min] <- NA
      X2[!is.na(X2)] <- 1
      X3 <- sum(X1, X2, na.rm = TRUE)

      if (Invert) {
        x[X3 == 1] <- NewVal
      } else {
        x[X3 == 2] <- NewVal
      }

    } else {

      if (Invert) {
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
