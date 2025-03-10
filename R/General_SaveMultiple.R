## |------------------------------------------------------------------------| #
# SaveMultiple ----
## |------------------------------------------------------------------------| #

#' Save multiple objects to their respective `.RData` files
#'
#' This function saves specified variables from the global environment to
#' separate `.RData` files. It allows for optional file prefixing and
#' overwriting of existing files.
#' @param Vars Character vector. Names of the variables to be saved. If `NULL`
#'   or any specified variable does not exist in the global environment, the
#'   function will stop with an error.
#' @param OutFolder Character. Path to the output folder where the `.RData`
#'   files will be saved. Defaults to the current working directory.
#' @param Overwrite Logical. Whether existing `.RData` files should be
#'   overwritten. If `FALSE` (Default) and files exist, the function will stop
#'   with an error message.
#' @param Prefix Character. Prefix of each output file name. Useful for
#'   organizing saved files or avoiding name conflicts. Defaults to an empty
#'   string.
#' @param Verbose Logical. Whether to print a message upon successful saving of
#'   files. Defaults to `FALSE`.
#' @name SaveMultiple
#' @author Ahmed El-Gabbas
#' @return The function is used for its side effect of saving files and does not
#'   return a value.
#' @export
#' @examples
#' \dontrun{
#'   TMP_Folder <- IASDT.R::Path(tempdir(), stringi::stri_rand_strings(1, 5))
#'   fs::dir_create(TMP_Folder)
#'
#'   # ----------------------------------------------
#'   # Save x1 and x2 to disk
#'   # ----------------------------------------------
#'   x1 = 10; x2 = 20
#'
#'   SaveMultiple(Vars = c("x1", "x2"), OutFolder = TMP_Folder)
#'
#'   list.files(path = TMP_Folder, pattern = "^.+.RData")
#'
#'   (x1Contents <- IASDT.R::LoadAs(IASDT.R::Path(TMP_Folder, "x1.RData")))
#'
#'   (x2Contents <- IASDT.R::LoadAs(IASDT.R::Path(TMP_Folder, "x2.RData")))
#'
#'   # ----------------------------------------------
#'   # Use Prefix
#'   # ----------------------------------------------
#'
#'   SaveMultiple(Vars = c("x1", "x2"), OutFolder = TMP_Folder, Prefix = "A_")
#'
#'   list.files(path = TMP_Folder, pattern = "^.+.RData")
#'
#'   # ----------------------------------------------
#'   # File exists, no save
#'   # ----------------------------------------------
#'   try(SaveMultiple(Vars = c("x1", "x2"), OutFolder = TMP_Folder))
#'
#'   # ----------------------------------------------
#'   # overwrite existing file
#'   # ----------------------------------------------
#'   x1 = 100; x2 = 200; x3 = 300
#'
#'   SaveMultiple(Vars = c("x1", "x2", "x3"),
#'      OutFolder = TMP_Folder, Overwrite = TRUE)
#'
#'   (x1Contents <- IASDT.R::LoadAs(IASDT.R::Path(TMP_Folder, "x1.RData")))
#'
#'   (x2Contents <- IASDT.R::LoadAs(IASDT.R::Path(TMP_Folder, "x2.RData")))
#'
#'      (x3Contents <- IASDT.R::LoadAs(IASDT.R::Path(TMP_Folder, "x3.RData")))
#' }

SaveMultiple <- function(
    Vars = NULL, OutFolder = getwd(),
    Overwrite = FALSE, Prefix = "", Verbose = FALSE) {

  # Check if Vars is provided correctly
  if (is.null(Vars) || !is.character(Vars)) {
    stop(
      "Vars should be a character vector for names of objects", call. = FALSE)
  }

  env <- rlang::new_environment()

  purrr::walk(
    .x = Vars, .f = ~assign(.x, get(.x, envir = parent.frame()), envir = env))

  # Check if all specified Vars are available in the caller environment
  missing_vars <- setdiff(Vars, ls(envir = env))

  if (length(missing_vars) > 0) {
    stop(
      "Variable(s) ", paste(missing_vars, collapse = " & "),
      " do not exist in the caller environment.\n", call. = FALSE)
  }
  fs::dir_create(OutFolder)

  # Check if files already exist
  FilesExist <- purrr::map_lgl(
    .x = IASDT.R::Path(OutFolder, paste0(Prefix, Vars, ".RData")),
    .f = file.exists) %>%
    any()

  if (FilesExist && isFALSE(Overwrite)) {
    message(
      "Some files already exist. No files are saved. ",
      "Please use overwrite = TRUE")
  } else {

    purrr::walk(
      .x = Vars,
      .f = ~{
        IASDT.R::SaveAs(
          InObj = get(.x, envir = env), OutObj = .x,
          OutPath = IASDT.R::Path(OutFolder, paste0(Prefix, .x, ".RData")))
      })

    AllExist <- all(
      file.exists(IASDT.R::Path(OutFolder, paste0(Prefix, Vars, ".RData"))))

    if (AllExist) {
      if (Verbose) {
        message("All files are saved to disk in ", OutFolder, " successfully.")
      }
    } else {
      message("Some files were not saved to disk! Please check again.")
    }
  }
  return(invisible(NULL))
}
