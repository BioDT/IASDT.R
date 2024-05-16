# |---------------------------------------------------| #
# List2RData ----
# |---------------------------------------------------| #

#' Split list items to separate RData files
#'
#' Split list items to separate RData files
#' @param List list object
#' @param Prefix prefix string
#' @param Dir output directory
#' @param Overwrite overwrite existing files?
#' @name List2RData
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export
#' @examples
#' # split iris data by species name
#' iris2 <- iris %>%
#'   tibble::tibble() %>%
#'   split(~Species)
#'
#' str(iris2, 1)
#'
#' (TMP_Folder <- file.path(tempdir(), stringi::stri_rand_strings(1, 5)))
#' list.files(TMP_Folder)
#'
#' List2RData(iris2, Dir = TMP_Folder)
#' list.files(TMP_Folder)

List2RData <- function(List, Prefix = "", Dir = getwd(), Overwrite = FALSE) {
  if (is.null(names(List))) {
    stop()
  } else {

    seq_along(List) %>%
      lapply(
        function(x) {

          if (Prefix != "") {
            Prefix <- paste0(Prefix, "_")
          }
          Name <- names(List[x])
          fs::dir_create(Dir)
          File <- paste0(Dir, "/", Prefix, Name, ".RData")

          if (file.exists(File)) {

            if (Overwrite) {
              assign(x = Name, value = List[[x]])
              save(list = Name, file = File)
            } else {
              "\n\nFile: {File} already exists" %>%
                stringr::str_glue() %>%
                cat()
            }
          } else {
            assign(x = Name, value = List[[x]])
            save(list = Name, file = File)
          }
        }
      ) %>%
      invisible()
  }
}
