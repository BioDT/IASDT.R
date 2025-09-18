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
  column_id <- column <- path_interim <- NULL

  # # ..................................................................... ###

  # Check if chunk_size is valid (greater than zero)
  if (!is.numeric(chunk_size) || chunk_size <= 0) {
    ecokit::stop_ctx(
      "chunk_size must be a positive number.", chunk_size = chunk_size,
      include_backtrace = TRUE)
  }

  # Check if `path_zip` is a character of length 1 and not empty. Also Check
  # if file exists
  if (!is.character(path_zip) || length(path_zip) != 1 || path_zip == "" ||
      is.null(path_zip)) {
    ecokit::stop_ctx(
      "`path_zip` must be a character of length 1 and not empty.",
      path_zip = path_zip, include_backtrace = TRUE)
  }

  # Check if `path_zip` is a valid path
  if (!file.exists(path_zip)) {
    ecokit::stop_ctx(
      "`path_zip` is not a valid path.", path_zip = path_zip,
      include_backtrace = TRUE)
  }

  # Checking required bash tools
  commands <- c("unzip", "cut", "sed", "split")
  available_commands <- purrr::map_lgl(commands, ecokit::check_system_command)
  if (!all(available_commands)) {
    missing_commands <- paste(commands[!available_commands], collapse = " + ")
    ecokit::stop_ctx(
      "Missing commands", missing_commands = missing_commands,
      include_backtrace = TRUE)
  }

  # ensure that `chunk_size` is not formatted in scientific notation
  chunk_size <- format(chunk_size, scientific = FALSE)

  # # ..................................................................... ###

  # Environment variables ----

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "path_interim", "DP_R_efforts_interim", FALSE, FALSE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  # # ..................................................................... ###

  out_prefix <- stringr::str_replace(basename(path_zip), ".zip$", "_") %>%
    fs::path(path_interim, .)

  # nolint start
  csv_file <- stringr::str_replace(basename(path_zip), ".zip$", ".csv")

  # extract column names and their numbers from the zipped file without
  # extraction read first line
  selected_col_names <- c(
    "taxonRank", "decimalLatitude", "decimalLongitude",
    "coordinateUncertaintyInMeters", "speciesKey")

  selected_columns <- "unzip -p {path_zip} {csv_file} | head -n 1" %>%
    stringr::str_glue() %>%
    ecokit::system_command() %>%
    # Split the first row into column names. Data is tab-separated
    stringr::str_split("\t") %>%
    magrittr::extract2(1) %>%
    dplyr::tibble(column = .) %>%
    # column number in the original data
    dplyr::mutate(column_id = seq_len(dplyr::n())) %>%
    # Only keep selected columns
    dplyr::filter(column %in% selected_col_names) %>%
    dplyr::pull(column_id) %>%
    paste0(collapse = ",")
  # nolint end

  extract_command <- stringr::str_glue(
    'unzip -p {path_zip} {csv_file} | cut -f{selected_columns} -d "\t" | \\
    sed -n "1!p" | split -l {chunk_size} -a 4 -d - {out_prefix} \\
    --additional-suffix=.txt')

  path_chunks <- tryCatch(
    ecokit::system_command(extract_command, r_object = TRUE),
    error = function(e) {
      ecokit::stop_ctx(
        paste0("Failed to execute system command: ", e$message),
        include_backtrace = TRUE)
    }
  )
  rm(path_chunks, envir = environment())

  return(
    list.files(
      path_interim, full.names = TRUE,
      pattern = paste0(basename(out_prefix), ".+txt"))
  )
}
