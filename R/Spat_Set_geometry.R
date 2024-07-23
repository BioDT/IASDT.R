## |------------------------------------------------------------------------| #
# Set_geometry ----
## |------------------------------------------------------------------------| #

#' Set the geometry column of a simple feature (sf) data frame in the pipe pipeline.
#'
#' Set the geometry column of a simple feature (sf) data frame, making it particularly useful in data processing pipelines. By specifying the name of the geometry column, users can ensure that spatial operations utilize the correct data.
#'
#' @param x A simple feature (sf) data frame. This is the data frame whose geometry column will be set or changed.
#' @param Name A string specifying the name of the geometry column to be used or set in the `x` data frame.
#' @name Set_geometry
#' @return The modified simple feature (sf) data frame with the updated geometry column. The function returns the original data frame `x` with its geometry column set to `Name`.
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

Set_geometry <- function(x, Name) {
  sf::st_geometry(x) <- Name
  return(x)
}
