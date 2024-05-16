# |---------------------------------------------------| #
# SetChainName ----
# |---------------------------------------------------| #

#' SetChainName
#'
#' SetChainName
#' @param Obj Obj
#' @param Chain Chain
#' @name SetChainName
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

SetChainName <- function(Obj, Chain) {
  Obj %>%
    unlist() %>%
    as.vector() %>%
    stats::setNames(paste0("Chain", unlist(Chain)))
}
