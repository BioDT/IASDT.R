## |------------------------------------------------------------------------| #
# Efforts_Split ----
## |------------------------------------------------------------------------| #

#' Split Order data into smaller chunks
#'
#' This function extracts Order data without extracting the zipped archive. The
#' function reads CSV files inside the zipped file, selects specified columns,
#' and splits the data into smaller chunks of specified row size, saving each
#' chunk as a separate file.
#' @param Path_Zip Character. The file path to the zip archive containing the
#'   CSV file. The file must be a ZIP archive containing a single CSV file.
#' @param Path_Output Character. The directory where the split files will be
#'   saved. The directory must exist.
#' @param ChunkSize Integer. The number of rows per chunk file. Default:
#'   `100,000`. Note: Larger chunk sizes may require significant memory and
#'   processing power.
#' @return A character vector of file path(s) to the created chunk files.
#' @name Efforts_Split
#' @author Ahmed El-Gabbas
#' @note This function is not intended to be used directly by the user or in the
#'   IAS-pDT, but only used inside the [Sampling_Efforts] function.
#' @export

Efforts_Split <- function(Path_Zip, Path_Output, ChunkSize = 100000) {
  ID <- Col <- NULL

  # Check if ChunkSize is valid (greater than zero)
  if (!is.numeric(ChunkSize) || ChunkSize <= 0) {
    stop("ChunkSize must be a positive number.")
  }

  # Checking required bash tools
  IASDT.R::CheckCommands(c("unzip", "cut", "sed", "split"))

  # ensure that `ChunkSize` is not formatted in scientific notation
  ChunkSize <- format(ChunkSize, scientific = FALSE)

  CSV_File <- stringr::str_replace(basename(Path_Zip), ".zip$", ".csv")
  OutPrefix <- stringr::str_replace(basename(Path_Zip), ".zip$", "_") %>%
    file.path(Path_Output, .)

  # extract column names and their numbers from the zipped file without
  # extraction read first line
  SelectedColNames <- c(
    "taxonRank", "decimalLatitude", "decimalLongitude",
    "coordinateUncertaintyInMeters", "speciesKey"
  )

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

  Command <- stringr::str_glue(
    'unzip -p {Path_Zip} {CSV_File} | cut -f{SelectedCols} -d "\t" | ',
    'sed -n "1!p" | split -l {ChunkSize} ',
    "-a 4 -d - {OutPrefix} --additional-suffix=.txt"
  )

  Path_Chunks <- tryCatch(
    IASDT.R::System(Command, RObj = FALSE),
    error = function(e) {
      stop("Failed to execute system command: ", e$message)
    }
  )

  return(
    list.files(
      Path_Output,
      full.names = TRUE,
      pattern = paste0(basename(OutPrefix), ".+txt")
    )
  )
}
