## |------------------------------------------------------------------------| #
# List2RData ----
## |------------------------------------------------------------------------| #

#' Split list items into separate RData files
#'
#' This function takes a named list and saves each element of the list as a separate RData file. The names of the list elements are used as the base for the filenames, optionally prefixed. Files are saved in the specified directory, with an option to overwrite existing files.
#'
#' @param List A named list object to be split into separate RData files.
#' @param Prefix A character string to prefix to each filename. If empty (default), no prefix is added.
#' @param Dir The directory where the RData files will be saved. Defaults to the current working directory.
#' @param Overwrite A logical indicating whether to overwrite existing files.  Defaults to `FALSE`, in which case files that already exist will not be overwritten, and a message will be printed for each such file.
#' @name List2RData
#' @author Ahmed El-Gabbas
#' @return Invisible NULL. The function is called for its side effect of saving files and does not return a value.
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
#' List2RData(List = iris2, Dir = TMP_Folder)
#' list.files(TMP_Folder)

List2RData <- function(List, Prefix = "", Dir = getwd(), Overwrite = FALSE) {

  # Validation Checks
  if (is.null(List) || length(List) == 0) {
    stop("List cannot be NULL or empty.")
  }

  if (is.null(names(List))) {
    stop("List names cannot be NULL.")
  }

  # Directory Creation
  fs::dir_create(Dir)

  # File Saving Loop --- iterates over each element in the list and writes it to a separate RData file. The filename is the element's name

  IASDT.R::lapply_(
    X = seq_along(List),
    FUN = function(x) {

      # construct filename using the element's name and the optional prefix
      FileName <- if (Prefix == "") {
        names(List)[x]
      } else {
        paste0(Prefix, "_", names(List)[x])
      }
      File <- file.path(Dir, paste0(FileName, ".RData"))

      # check if the file already exists. If it does and Overwrite is FALSE, it prints a message indicating that the file already exists and will not be overwritten.
      if (file.exists(File) && !Overwrite) {
        stringr::str_glue("\n\nFile: {File} already exists. No files were created.") %>%
          IASDT.R::CatTime()
      } else {
        # If the file does not exist or Overwrite is TRUE, saves the list element as an RData file
        IASDT.R::SaveAs(InObj = List[[x]], OutObj = FileName,  OutPath = File)
      }
    })

  return(invisible(NULL))
}
