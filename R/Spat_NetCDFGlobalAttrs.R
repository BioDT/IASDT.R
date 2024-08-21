## |------------------------------------------------------------------------| #
# NetCDFGlobalAttrs ------
## |------------------------------------------------------------------------| #

#' Get global attributes for NetCDF files
#'
#' This function opens a NetCDF file, extracts all global attributes, and
#' returns them as a character vector where each element is an attribute
#' name-value pair.
#' @name NetCDFGlobalAttrs
#' @param nc A character string specifying the path to the NetCDF file. If
#'   `NULL`, the function will stop with an error message.
#' @return A character vector where each element is a global attribute.
#' @references [Click here](https://github.com/rspatial/terra/issues/1443)
#' @export

NetCDFGlobalAttrs <- function(nc = NULL) {

  # Input Validation
  if (is.null(nc)) {
    stop("Input file cannot be NULL", .call = FALSE)
  }

  # Open the NetCDF File
  nc <- RNetCDF::open.nc(nc)

  # Extracting Global Attributes
  GlobAttrs <- purrr::map_chr(
    .x = (seq_len(RNetCDF::file.inq.nc(nc)$ngatt) - 1),
    .f = ~{
      attn <- RNetCDF::att.inq.nc(nc, "NC_GLOBAL", .x)$name
      attv <- RNetCDF::att.get.nc(nc, "NC_GLOBAL", .x)
      paste0(attn, "=", attv)
    })

  # Closing the NetCDF File
  RNetCDF::close.nc(nc)

  return(GlobAttrs)
}
