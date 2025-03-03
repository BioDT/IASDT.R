## |------------------------------------------------------------------------| #
# Efforts_Split ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @name Efforts_data
#' @rdname Efforts_data
#' @order 5
#' @export

Efforts_Split <- function(
  Path_Zip = NULL, EnvFile = ".env", ChunkSize = 100000L) {

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  ID <- Col <- Path_Interim <- NULL

  # # ..................................................................... ###

  # Check if ChunkSize is valid (greater than zero)
  if (!is.numeric(ChunkSize) || ChunkSize <= 0) {
    stop("ChunkSize must be a positive number.", call. = FALSE)
  }

  # Check if `Path_Zip` is a character of length 1 and not empty. Also Check
  # if file exists
  if (!is.character(Path_Zip) || length(Path_Zip) != 1 || Path_Zip == "" ||
    is.null(Path_Zip)) {
    stop(
      "`Path_Zip` must be a character of length 1 and not empty.",
      call. = FALSE)
  }

  # Check if `Path_Zip` is a valid path
  if (!file.exists(Path_Zip)) {
    stop("`Path_Zip` is not a valid path.", call. = FALSE)
  }


  # Checking required bash tools
  Commands <- c("unzip", "cut", "sed", "split")
  CommandsAvail <- purrr::map_lgl(Commands, IASDT.R::CheckCommands)
  if (!all(CommandsAvail)) {
    Missing <- paste(Commands[!CommandsAvail], collapse = " + ")
    stop("The following command(s) are not available: ", Missing, call. = FALSE)
  }

  # ensure that `ChunkSize` is not formatted in scientific notation
  ChunkSize <- format(ChunkSize, scientific = FALSE)

  # # ..................................................................... ###

  # Environment variables ----
  EnvVars2Read <- tibble::tribble(
    ~VarName, ~Value, ~CheckDir, ~CheckFile,
    "Path_Interim", "DP_R_Efforts_interim", FALSE, FALSE)
  # Assign environment variables and check file and paths
  IASDT.R::AssignEnvVars(EnvFile = EnvFile, EnvVarDT = EnvVars2Read)
  rm(EnvVars2Read, envir = environment())

  # # ..................................................................... ###

  OutPrefix <- stringr::str_replace(basename(Path_Zip), ".zip$", "_") %>%
    IASDT.R::Path(Path_Interim, .)

  # nolint start
  CSV_File <- stringr::str_replace(basename(Path_Zip), ".zip$", ".csv")

  # extract column names and their numbers from the zipped file without
  # extraction read first line
  SelectedColNames <- c(
    "taxonRank", "decimalLatitude", "decimalLongitude",
    "coordinateUncertaintyInMeters", "speciesKey")

  SelectedCols <- "unzip -p {Path_Zip} {CSV_File} | head -n 1" %>%
    stringr::str_glue() %>%
    IASDT.R::System() %>%
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
    'unzip -p {Path_Zip} {CSV_File} | cut -f{SelectedCols} -d "\t" | \\
    sed -n "1!p" | split -l {ChunkSize} -a 4 -d - {OutPrefix} \\
    --additional-suffix=.txt')

  Path_Chunks <- tryCatch(
    IASDT.R::System(Command, RObj = FALSE),
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
