## |------------------------------------------------------------------------| #
# SaveMultiple ----
## |------------------------------------------------------------------------| #

#' Save multiple objects to their respective `.RData` files
#'
#' This function saves specified variables from the global environment to separate `.RData` files. It allows for optional file prefixing and overwriting of existing files.
#'
#' @param Vars A character vector specifying the names of the variables to be saved. If `NULL` or any specified variable does not exist in the global environment, the function will stop with an error.
#' @param OutFolder  A string specifying the path to the output folder where the `.RData` files will be saved. Defaults to the current working directory.
#' @param Overwrite A logical value indicating whether existing `.RData` files should be overwritten. If `FALSE` (Default) and files exist, the function will stop with an error message.
#' @param Prefix A string to be prefixed to each output file name. Useful for organizing saved files or avoiding name conflicts. Defaults to an empty string.
#' @param Verbose A logical value indicating whether to print a message upon successful saving of files. Defaults to `FALSE`.
#' @name SaveMultiple
#' @author Ahmed El-Gabbas
#' @return The function is used for its side effect of saving files and does not return a value.
#' @export
#' @examples
#' TMP_Folder <- file.path(tempdir(), stringi::stri_rand_strings(1, 5))
#' fs::dir_create(TMP_Folder)
#'
#' # ----------------------------------------------
#' # Save x1 and x2 to disk
#' # ----------------------------------------------
#' x1 = 10; x2 = 20
#'
#' SaveMultiple(Vars = c("x1", "x2"), OutFolder = TMP_Folder)
#'
#' list.files(path = TMP_Folder, pattern = "^.+.RData")
#'
#' (x1Contents <- IASDT.R::LoadAs(file.path(TMP_Folder, "x1.RData")))
#'
#' (x2Contents <- IASDT.R::LoadAs(file.path(TMP_Folder, "x2.RData")))
#'
#' # ----------------------------------------------
#' # Use Prefix
#' # ----------------------------------------------
#'
#' SaveMultiple(Vars = c("x1", "x2"), OutFolder = TMP_Folder, Prefix = "A_")
#'
#' list.files(path = TMP_Folder, pattern = "^.+.RData")
#'
#' # ----------------------------------------------
#' # File exists, no save
#' # ----------------------------------------------
#' try(SaveMultiple(Vars = c("x1", "x2"), OutFolder = TMP_Folder))
#'
#' # ----------------------------------------------
#' # overwrite existing file
#' # ----------------------------------------------
#' x1 = 100; x2 = 200; x3 = 300
#'
#' SaveMultiple(Vars = c("x1", "x2", "x3"), OutFolder = TMP_Folder, Overwrite = TRUE)
#'
#' (x1Contents <- IASDT.R::LoadAs(file.path(TMP_Folder, "x1.RData")))
#'
#' (x2Contents <- IASDT.R::LoadAs(file.path(TMP_Folder, "x2.RData")))
#'
#' (x3Contents <- IASDT.R::LoadAs(file.path(TMP_Folder, "x3.RData")))

SaveMultiple <- function(
    Vars = NULL, OutFolder = getwd(),
    Overwrite = FALSE, Prefix = "", Verbose = FALSE) {

  # Check if Vars is provided correctly
  if (is.null(Vars) || !is.character(Vars) || any(!Vars %in% ls(envir = rlang::caller_env()))) {
    stop("Vars should be a character vector for names of objects available in the specified environment")
  }

  # Get variables in the caller environment
  GlobalVars <- ls(rlang::caller_env())

  # Check if all specified Vars are available in the caller environment
  missing_vars <- setdiff(Vars, GlobalVars)

  if (length(missing_vars) > 0) {
    stop(paste0("Variable(s) ", paste0(missing_vars, collapse = " & "), " do not exist in the caller environment.\n"))
  } else {
    fs::dir_create(OutFolder)
  }

  # Check if files already exist
  FilesExist <- sapply(
    X = file.path(OutFolder, paste0(Prefix, Vars, ".RData")), FUN = file.exists) %>%
    any()

  if (FilesExist && !Overwrite) {
    message("Some files already exist. No files are saved. Please use overwrite = TRUE")
  } else {

    purrr::walk(
      .x = Vars,
      .f = ~{
        # Val <- get(.x, envir = parent.frame(n = 4))
        # assign(rlang::eval_tidy(.x), Val)
        # save(
        #   list = paste0(.x),
        #   file = file.path(OutFolder, paste0(Prefix, .x, ".RData")))
        Val <- get(.x, envir = rlang::caller_env())
        save(
          list = .x,
          file = file.path(OutFolder, paste0(Prefix, .x, ".RData")))
      })

    if (all(
      file.exists(file.path(OutFolder, paste0(Prefix, Vars, ".RData"))))) {
      if (Verbose) {
        message(paste0("All files are saved to disk in ", OutFolder, " successfully."))
      }
    } else {
      message("Some files were not saved to disk! Please check again.")
    }
  }
  return(invisible(NULL))
}
