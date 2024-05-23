## |------------------------------------------------------------------------| #
# SetChainName ----
## |------------------------------------------------------------------------| #

#' Set name of the chain to a vector of strings
#'
#' Set name of the chain to a vector of strings
#' 
#' @param Obj Input object name
#' @param Chain Integer. Chain number
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
