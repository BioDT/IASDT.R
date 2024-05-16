# |---------------------------------------------------| #
# NetCDFGlobalAttrs ------
# |---------------------------------------------------| #

#' Get global attributes for NetCDF files
#'
#' Get global attributes for NetCDF files
#'
#' @name NetCDFGlobalAttrs
#' @param nc path to the NetCDF file
#' @export
#' @description
#' https://github.com/rspatial/terra/issues/1443

NetCDFGlobalAttrs <- function(nc) {
  nc <- RNetCDF::open.nc(nc)
  GlobAttrs <- purrr::map_chr(
    seq_len(RNetCDF::file.inq.nc(nc)$ngatt) - 1,
    ~{
      attn <- RNetCDF::att.inq.nc(nc, "NC_GLOBAL", .x)$name
      attv <- RNetCDF::att.get.nc(nc, "NC_GLOBAL", .x)
      paste0(attn, "=", attv)
    })
  RNetCDF::close.nc(nc)
  GlobAttrs
}
