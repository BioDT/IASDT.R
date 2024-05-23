## |------------------------------------------------------------------------| #
# AddImg2Plot ----
## |------------------------------------------------------------------------| #
#
#' Add image to plot
#'
#' Add image to plot
#'
#' @name AddImg2Plot
#' @references [Click here](https://stackoverflow.com/questions/27800307/)
#' @export
#' @param obj an image file imported as an array (e.g. `png::readPNG`, `jpeg::readJPEG`)
#' @param x mid x coordinate for image
#' @param y mid y coordinate for image
#' @param width width of image (in x coordinate units)
#' @param interpolate (passed to `graphics::rasterImage`) A logical vector (or scalar) indicating whether to apply linear interpolation to the image when drawing. Default = `TRUE`.
#' @examples
#' LoadPackages(png)
#' myurl <- "https://upload.wikimedia.org/wikipedia/commons/e/e1/Jupiter_%28transparent%29.png"
#' z <- tempfile()
#' download.file(myurl, z, mode="wb")
#' pic <- readPNG(z)
#' file.remove(z) # cleanup
#'
#' image(volcano)
#' AddImg2Plot(pic, x = 0.3, y = 0.5, width = 0.2)
#' AddImg2Plot(pic, x = 0.7, y = 0.7, width = 0.2)
#' AddImg2Plot(pic, x = 0.7, y = 0.2, width = 0.1)

AddImg2Plot <- function(
    obj, x = NULL, y = NULL, width = NULL, interpolate = TRUE) {

  if (any(is.null(x), is.null(y), is.null(width))) {
    stop("Must provide args 'x', 'y', and 'width'")
  }

  USR <- graphics::par()$usr # A vector of the form c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region
  PIN <- graphics::par()$pin # The current plot dimensions, (width, height), in inches
  DIM <- dim(obj) # number of x-y pixels for the image
  ARp <- DIM[1] / DIM[2] # pixel aspect ratio (y/x)
  WIDi <- width / (USR[2] - USR[1]) * PIN[1] # convert width units to inches
  HEIi <- WIDi * ARp # height in inches
  HEIu <- HEIi / PIN[2] * (USR[4] - USR[3]) # height in units
  graphics::rasterImage(
    image = obj, xleft = x - (width / 2), xright = x + (width / 2),
    ybottom = y - (HEIu / 2), ytop = y + (HEIu / 2), interpolate = interpolate)
}
