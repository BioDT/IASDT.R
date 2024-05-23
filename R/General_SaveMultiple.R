## |------------------------------------------------------------------------| #
# SaveMultiple ----
## |------------------------------------------------------------------------| #

#' Save multiple objects to their respective `.RData` files
#'
#' Save multiple objects to their respective `.RData` files
#' @param Vars variables to save
#' @param OutFolder output path
#' @param Overwrite overwrite existing files? If file already exist, no files are saved unless `Overwrite` argument is set as `TRUE`
#' @param Prefix String prefix of the output file
#' @param Verbose Show message after saving files
#' @name SaveMultiple
#' @author Ahmed El-Gabbas
#' @return NULL
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
#' SaveMultiple(Vars = cc(x1, x2), OutFolder = TMP_Folder)
#'
#' list.files(path = TMP_Folder, pattern = "^.+.RData")
#'
#' (x1Contents <- LoadAs(file.path(TMP_Folder, "x1.RData")))
#'
#' (x2Contents <- LoadAs(file.path(TMP_Folder, "x2.RData")))
#'
#' # ----------------------------------------------
#' # Use Prefix
#' # ----------------------------------------------
#'
#' SaveMultiple(Vars = cc(x1, x2), OutFolder = TMP_Folder, Prefix = "A_")
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
#' SaveMultiple(Vars = cc(x1, x2, x3), OutFolder = TMP_Folder, Overwrite = TRUE)
#'
#' (x1Contents <- LoadAs(file.path(TMP_Folder, "x1.RData")))
#'
#' (x2Contents <- LoadAs(file.path(TMP_Folder, "x2.RData")))
#'
#' (x3Contents <- LoadAs(file.path(TMP_Folder, "x3.RData")))

SaveMultiple <- function(
    Vars = NULL, OutFolder = getwd(),
    Overwrite = FALSE, Prefix = "", Verbose = FALSE) {

  SkipFun <- any(
    is.null(Vars), !inherits(Vars, "character"),
    any(!Vars %in% ls(envir = rlang::caller_env())))

  if (SkipFun) {
    cat(paste0(crayon::red("Vars should be a character vector for names of objects available in the global environment"), "\n"))
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }

  GlobalVars <- ls(rlang::caller_env())
  VarsInGlobal <- Vars %>%
    purrr::map_lgl(~.x %in% GlobalVars) %>%
    all() %>%
    magrittr::not()

  if (VarsInGlobal) {
    "Please check that all variables to be saved exist at your global environment" %>%
      crayon::red() %>%
      cat(sep = "\n")

    MissedVars <- setdiff(Vars, GlobalVars)
    cat(crayon::red(paste0("Variable ", MissedVars, " does not exist.\n")))
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  } else {
    if (is.null(OutFolder)) {
      OutFolder <- getwd()
    }
    fs::dir_create(OutFolder)
  }

  FilesExist <- file.path(OutFolder, paste0(Prefix, Vars, ".RData")) %>%
    file.exists() %>%
    any()


  if (FilesExist && !Overwrite) {
    "Some files already exist.\nNo files are saved. Please use overwrite = TRUE" %>%
      crayon::red() %>%
      cat(sep = "\n")
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }

  purrr::walk(
    .x = Vars,
    .f = ~{
      Val <- get(.x, envir = parent.frame(n = 4))
      assign(rlang::eval_tidy(.x), Val)
      save(
        list = paste0(.x),
        file = file.path(OutFolder, paste0(Prefix, .x, ".RData")))
    })


  if (all(file.exists(file.path(OutFolder, paste0(Prefix, Vars, ".RData"))))) {
    if (Verbose) {
      cat(paste0(crayon::blue("All files are saved to disk", crayon::red(OutFolder), "successfully."), "\n"))
    }
  } else {
    cat(paste0(crayon::red("Some files were not saved to disk!\nplease check again"), "\n"))
  }
}
