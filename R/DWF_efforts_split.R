## |------------------------------------------------------------------------| #
# efforts_split ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @name efforts_data
#' @rdname efforts_data
#' @order 5
#' @export

efforts_split <- function(
  path_zip = NULL, env_file = ".env", chunk_size = 100000L) {

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  ID <- Col <- Path_Interim <- NULL

  # # ..................................................................... ###

  # Check if chunk_size is valid (greater than zero)
  if (!is.numeric(chunk_size) || chunk_size <= 0) {
    stop("chunk_size must be a positive number.", call. = FALSE)
  }

  # Check if `path_zip` is a character of length 1 and not empty. Also Check
  # if file exists
  if (!is.character(path_zip) || length(path_zip) != 1 || path_zip == "" ||
    is.null(path_zip)) {
    stop(
      "`path_zip` must be a character of length 1 and not empty.",
      call. = FALSE)
  }

  # Check if `path_zip` is a valid path
  if (!file.exists(path_zip)) {
    stop("`path_zip` is not a valid path.", call. = FALSE)
  }


  # Checking required bash tools
  Commands <- c("unzip", "cut", "sed", "split")
  CommandsAvail <- purrr::map_lgl(Commands, IASDT.R::check_system_command)
  if (!all(CommandsAvail)) {
    Missing <- paste(Commands[!CommandsAvail], collapse = " + ")
    stop("The following command(s) are not available: ", Missing, call. = FALSE)
  }

  # ensure that `chunk_size` is not formatted in scientific notation
  chunk_size <- format(chunk_size, scientific = FALSE)

  # # ..................................................................... ###

  # Environment variables ----
  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_Interim", "DP_R_Efforts_interim", FALSE, FALSE)
  # Assign environment variables and check file and paths
  IASDT.R::assign_env_vars(
    env_file = env_file, env_variables_data = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  OutPrefix <- stringr::str_replace(basename(path_zip), ".zip$", "_") %>%
    IASDT.R::path(Path_Interim, .)

  # nolint start
  CSV_File <- stringr::str_replace(basename(path_zip), ".zip$", ".csv")

  # extract column names and their numbers from the zipped file without
  # extraction read first line
  SelectedColNames <- c(
    "taxonRank", "decimalLatitude", "decimalLongitude",
    "coordinateUncertaintyInMeters", "speciesKey")

  SelectedCols <- "unzip -p {path_zip} {CSV_File} | head -n 1" %>%
    stringr::str_glue() %>%
    IASDT.R::system_command() %>%
    # Split the first row into column names. Data is tab-separated
    stringr::str_split("\t") %>%
    magrittr::extract2(1) %>%
    dplyr::tibble(Col = .) %>%
    # column number in the original data
    dplyr::mutate(ID = seq_len(dplyr::n())) %>%
    # Only keep selected columns
    dplyr::filter(Col %in% SelectedColNames) %>%
    dplyr::pull(ID) %>%
    paste0(collapse = ",")
  # nolint end

  Command <- stringr::str_glue(
    'unzip -p {path_zip} {CSV_File} | cut -f{SelectedCols} -d "\t" | \\
    sed -n "1!p" | split -l {chunk_size} -a 4 -d - {OutPrefix} \\
    --additional-suffix=.txt')

  Path_Chunks <- tryCatch(
    IASDT.R::system_command(Command, R_object = FALSE),
    error = function(e) {
      stop("Failed to execute system command: ", e$message, call. = FALSE)
    }
  )
  rm(Path_Chunks, envir = environment())

  return(
    list.files(
      Path_Interim, full.names = TRUE,
      pattern = paste0(basename(OutPrefix), ".+txt"))
  )
}
