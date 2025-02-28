## |------------------------------------------------------------------------| #
# SaveAs ----
## |------------------------------------------------------------------------| #

#' Save an object to a file with a new name
#'
#' This function saves an R object to a specified file path with a potentially
#' new name. It is useful for renaming objects during the save process. The
#' function supports saving objects in `RData`, `qs2`, `feather`, and `rds`
#' formats. The format is determined by the extension of the file path.
#'
#' @param InObj The input object to be saved. This can be an actual R object or
#'   a character string representing the name of an object.
#' @param OutObj Character. The new name for the saved `RData` object. This name
#'   is used when the object is loaded back into R. Default is `NULL`. This is
#'   required when saving `RData` files.
#' @param OutPath Character. File path (ends with either `*.RData`, `*.qs2`,
#'   `feather`, and `rds`) where the object be saved. This includes the
#'   directory and the file name.
#' @param nthreads Character. Number of threads to use when compressing data.
#'   See [qs2::qs_save].
#' @param feather_compression Character. The compression algorithm to use when
#'   saving the object in the `feather` format. The default is "zstd". See
#'   [arrow::write_feather].
#' @param ... Additional arguments to be passed to the respective save
#'   functions. [base::save] for `RData` files; [qs2::qs_save] for `qs2` files;
#'   [arrow::write_feather] for `feather` files; and [base::saveRDS] for `rds`
#'   files.
#' @name SaveAs
#' @author Ahmed El-Gabbas
#' @return The function does not return a value but saves an object to the
#'   specified file path.
#' @export
#' @examples
#' TMP_Folder <- IASDT.R::Path(tempdir(), stringi::stri_rand_strings(1, 5))
#' fs::dir_create(TMP_Folder)
#' list.files(TMP_Folder)
#'
#' # save iris data in `iris2.RData` with `iris2` object name
#' SaveAs(iris, "iris2", IASDT.R::Path(TMP_Folder, "iris2.RData"))
#' list.files(TMP_Folder, pattern = "^.+.RData")
#'
#' (load(IASDT.R::Path(TMP_Folder, "iris2.RData")))
#'
#' tibble::tibble(iris2)

SaveAs <- function(
    InObj = NULL, OutObj = NULL, OutPath = NULL, nthreads = 5L,
    feather_compression = "zstd", ...) {

  if (is.null(InObj) || is.null(OutPath)) {
    stop("`InObj` and `OutPath` cannot be NULL", call. = FALSE)
  }

  if (inherits(InObj, "character")) {
    InObj <- get(InObj)
  }

  Extension <- stringr::str_to_lower(tools::file_ext(OutPath))

  if (!Extension %in% c("qs2", "rdata", "feather", "rds")) {
    stop(
      "Extension of `OutPath` must be either 'qs2', ",
      "'rdata', 'feather', or 'rds' (case-insensitive).", call. = FALSE)
  }

  # Create directory if not available
  fs::dir_create(dirname(OutPath))

  switch(
    Extension,
    qs2 = {
      qs2::qs_save(object = InObj, file = OutPath, nthreads = nthreads, ...)
    },
    rdata = {
      if (is.null(OutObj)) {
        stop("`OutObj` cannot be `NULL` for saving RData files", call. = FALSE)
      }
      OutObj <- eval(OutObj)
      assign(OutObj, InObj)
      save(list = OutObj, file = OutPath, ...)
    },
    feather = {
      arrow::write_feather(
        x = InObj, sink = OutPath, compression = feather_compression, ...)
    },
    rds = {
      saveRDS(object = InObj, file = OutPath, ...)
    },
    stop("Invalid file extension", call. = FALSE)
  )

  return(invisible())
}
