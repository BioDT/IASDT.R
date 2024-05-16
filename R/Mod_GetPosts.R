# |---------------------------------------------------| #
# GetPosts ----
# |---------------------------------------------------| #

#' Combining posteriors into Hmsc object
#'
#' Combining posteriors into Hmsc object
#' @param FilePath vector of `rds` files created by Hmsc-HPC
#' @name GetPosts
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

GetPosts <- function(FilePath) {
  readRDS(file = FilePath) %>%
    magrittr::extract2(1) %>%
    jsonify::from_json() %>%
    magrittr::extract2(1)
}
