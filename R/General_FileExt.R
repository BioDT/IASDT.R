## |------------------------------------------------------------------------| #
# FileExt ----
## |------------------------------------------------------------------------| #
#
#' Get file extension
#'
#' Get file extension
#' @name FileExt
#' @author Ahmed El-Gabbas
#' @return NULL
#' @param Path File path
#' @examples
#' FileExt("File.doc")
#' FileExt("D:/File.doc")
#'
#' @export

FileExt <- function(Path) {
  Ext <- Path %>%
    basename() %>%
    stringr::str_split("\\.", simplify = TRUE)
  Ext[length(Ext)]
}
