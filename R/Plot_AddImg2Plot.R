## |------------------------------------------------------------------------| #
# AddImg2Plot ----
## |------------------------------------------------------------------------| #
#
#' Add an image to an existing plot in R
#'
#' This function allows the user to add an image to an existing plot in R by
#' specifying the image object, its position, and its size. The function
#' calculates the necessary dimensions and places the image accordingly. The
#' function uses the existing plot's coordinate system and accounts for the
#' current plot dimensions to ensure accurate placement of the image. It also
#' allows for interpolation, which can improve the visual quality of the image.
#' @name AddImg2Plot
#' @source The source code of this function was taken from this
#'   [stackoverflow](https://stackoverflow.com/questions/27800307/) question.
#' @export
#' @param obj The image object to be added to the plot, expected to be an
#'   array-like structure (e.g., as read by [png::readPNG] or [jpeg::readJPEG]).
#' @param x,y Numeric, the x-coordinate or y-coordinate (in plot units) at which
#'   the center of the image should be placed.
#' @param width Numeric, the desired width of the image in plot units (not
#'   pixels or inches). The function will calculate the corresponding height to
#'   preserve the image's aspect ratio.
#' @param interpolate Logical, whether to apply linear interpolation to the
#'   image when drawing. Defaults to `TRUE`. Passed directly to
#'   [graphics::rasterImage]. Interpolation can improve image quality but may
#'   take longer to render.
#' @return This function does not return a value but modifies the current plot
#'   by adding an image.
#' @note The function will stop with an error message if any of the required
#'   arguments (`obj`, `x`, `y`, `width`) are `NULL`.
#' @examples
#' library(png)
#' myurl <- paste0("https://upload.wikimedia.org/wikipedia/commons/",
#'     "e/e1/Jupiter_%28transparent%29.png")
#' z <- tempfile()
#' download.file(myurl, z, mode="wb", quiet = TRUE)
#' pic <- png::readPNG(z)
#' file.remove(z) # cleanup
#'
#' image(volcano)
#' AddImg2Plot(pic, x = 0.3, y = 0.5, width = 0.2)
#' AddImg2Plot(pic, x = 0.7, y = 0.7, width = 0.2)
#' AddImg2Plot(pic, x = 0.7, y = 0.2, width = 0.1)

AddImg2Plot <- function(obj, x, y, width, interpolate = TRUE) {

  if (is.null(obj) || is.null(x) || is.null(y) || is.null(width)) {
    stop("Must provide args 'obj', 'x', 'y', and 'width'")
  }

  # A vector of the form c(x1, x2, y1, y2) giving the extremes of the user
  # coordinates of the plotting region
  USR <- graphics::par()$usr

  # The current plot dimensions, (width, height), in inches
  PIN <- graphics::par()$pin

  # number of x-y pixels for the image
  DIM <- dim(obj)

  # pixel aspect ratio (y/x)
  ARp <- DIM[1] / DIM[2]

  # convert width units to inches
  WIDi <- width / (USR[2] - USR[1]) * PIN[1]

  # height in inches
  HEIi <- WIDi * ARp

  # height in units
  HEIu <- HEIi / PIN[2] * (USR[4] - USR[3])

  graphics::rasterImage(
    image = obj, xleft = x - (width / 2), xright = x + (width / 2),
    ybottom = y - (HEIu / 2), ytop = y + (HEIu / 2),
    interpolate = interpolate)
}
