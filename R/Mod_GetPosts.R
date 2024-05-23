## |------------------------------------------------------------------------| #
# GetPosts ----
## |------------------------------------------------------------------------| #

#' Combining posteriors exported by `Hmsc-HPC` into Hmsc object
#'
#' Combining posteriors exported by `Hmsc-HPC` into Hmsc object
#' @param FilePath String. Vector of the file paths of the `rds` files created by `Hmsc-HPC`
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
