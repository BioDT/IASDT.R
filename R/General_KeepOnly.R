## |------------------------------------------------------------------------| #
# KeepOnly ----
## |------------------------------------------------------------------------| #

#' Keep only certain objects in memory, all other objects will be removed
#'
#' Keep only certain objects in memory, all other objects will be removed
#'
#' @name KeepOnly
#' @param Obj character vector for objects to be kept in memory
#' @param Verbose Should the names of kept and removed variables printed? Default: `TRUE`.
#' @author Ahmed El-Gabbas
#' @export
#' @examples
#' A <- B <- C <- 15
#' ls()
#'
#' KeepOnly("A")
#'
#' ls()
#' rm(list = ls())
#'
#'
#' A <- B <- C <- 15
#' KeepOnly(c("A","B"))
#' ls()

KeepOnly <- function(Obj = NULL, Verbose = TRUE) {
  if (is.null(Obj) || length(Obj) == 0) {
    stop()
  }
  AllObjects <- ls(pos = parent.frame())
  RemObjects <- setdiff(AllObjects, Obj)
  if (Verbose) {
    cat(crayon::red(paste0("Removed Variables (", length(RemObjects), "): ")), crayon::blue(paste0(seq_along(RemObjects), ":", RemObjects, collapse = " ||  ")), sep = "")
  }
  rm(list = RemObjects, pos = parent.frame())
  if (Verbose) {
    cat(crayon::red(paste0("\nKept Variables (", length(Obj), "): ")), crayon::blue(paste0(seq_along(Obj), ":", Obj, collapse = " ||  ")), "\n", sep = "")
  }
}
