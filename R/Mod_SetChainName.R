## |------------------------------------------------------------------------| #
# SetChainName ----
## |------------------------------------------------------------------------| #

#' Set name of the chain to a vector of strings
#'
#' This function takes an input object and a chain number, then assigns names to the elements of the object. #' The names are created by prefixing 'Chain' to the chain number for each element.
#'
#' @param Obj An object to which the chain names will be assigned.
#' @param Chain An integer or vector of integers representing the chain number(s).
#' @name SetChainName
#' @author Ahmed El-Gabbas
#' @return Returns the input object with elements named according to the specified chain number(s). The function modifies the names attribute of the object.
#' @export

SetChainName <- function(Obj, Chain) {

  if (is.null(Obj) || is.null(Chain)) {
    stop("Obj and Chain cannot be empty")
  }

  Obj %>%
    unlist() %>%
    as.vector() %>%
    stats::setNames(paste0("Chain", unlist(Chain)))
}
