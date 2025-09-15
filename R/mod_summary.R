## |------------------------------------------------------------------------| #
# mod_summary ----
## |------------------------------------------------------------------------| #

#' Summary of Hmsc model parameters
#'
#' This function provides a comprehensive summary of Hmsc model parameters,
#' including `Alpha`, `Beta`, `Rho`, and `Omega`. It processes the model's
#' output, performs statistical summaries, and optionally returns the summarised
#' data.
#' @param path_coda Character. Path to the `.qs2` / `.RData` file containing the
#'   coda object.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param return_data Logical. Whether the summarised data should be returned as
#'   an R object. If `TRUE`, the function returns a list containing summaries of
#'   `Alpha`, `Beta`, `Rho`, and `Omega` parameters. The default value is
#'   `FALSE`, which means the function will not return any data but will save
#'   the summaries to a specified directory.
#' @param spatial_model Logical. Whether the model is a spatial model. If `TRUE`
#'   (default), the function will also process the `Alpha` parameter.
#' @author Ahmed El-Gabbas
#' @return If `return_data` is `FALSE` (default), the function does not return
#'   anything and saves the summaries to a directory. If `return_data` is
#'   `TRUE`, it also returns the data as R object.
#' @export
#' @name mod_summary

mod_summary <- function(
    path_coda = NULL, env_file = ".env", return_data = FALSE,
    spatial_model = TRUE) {

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  if (is.null(path_coda)) {
    ecokit::stop_ctx(
      "`path_coda` cannot be empty", path_coda = path_coda,
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  VarSp <- Sp1_abb <- Sp2_abb <- ias_id <- Val <- taxon_name <- species_name <-
    Sp1_Species_name <- Sp2_Species_name <- Q2_5 <- Q97_5 <- Sp1_Sp_abb <-
    Sp2_Sp_abb <- taxa_info_file <- NULL

  # # ..................................................................... ###

  # Prepare Species list -----
  ecokit::cat_time("Prepare Species list")

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  env_vars_to_read <- tibble::tribble(
    ~var_name, ~value, ~check_dir, ~check_file,
    "taxa_info_file", "DP_R_taxa_info", FALSE, TRUE)
  # Assign environment variables and check file and paths
  ecokit::assign_env_vars(
    env_file = env_file, env_variables_data = env_vars_to_read)
  rm(env_vars_to_read, envir = environment())

  # # ..................................................................... ###

  # DataPrep helper function -------
  DataPrep <- function(DT) {
    as.data.frame(DT) %>%
      tibble::rownames_to_column(var = "VarSp") %>%
      tibble::as_tibble() %>%
      dplyr::rename(
        dplyr::any_of(c(
          Naive_SE = "Naive SE", TimeSeries_SE = "Time-series SE",
          Q2_5 = "2.5%", Q25 = "25%", Q50 = "50%", Q75 = "75%",
          Q97_5 = "97.5%")))
  }

  # # ..................................................................... ###

  # Loading coda object ------
  ecokit::cat_time("Loading coda object")
  Coda <- ecokit::load_as(path_coda)
  coda_names <- names(Coda)

  ecokit::cat_time("Beta summary", level = 1L)
  Beta_Summary <- summary(Coda$Beta)

  if (("Alpha" %in% coda_names) && spatial_model) {
    ecokit::cat_time("Alpha summary", level = 1L)
    Alpha_Summary <- summary(Coda$Alpha[[1]])
  } else {
    Alpha_Summary <- NULL
  }

  if ("Omega" %in% coda_names) {
    ecokit::cat_time("Omega summary", level = 1L)
    Omega_Summary <- summary(Coda$Omega[[1]])
  } else {
    Omega_Summary <- NULL
  }

  if ("Rho" %in% coda_names) {
    ecokit::cat_time("Rho summary", level = 1L)
    Rho_Summary <- summary(Coda$Rho)
  } else {
    Rho_Summary <- NULL
  }

  rm(Coda, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  ecokit::cat_time("Extracting summary data")

  # Beta ------
  ecokit::cat_time("Beta", level = 1L)
  Beta_Summary <- DataPrep(Beta_Summary$statistics) %>%
    dplyr::full_join(DataPrep(Beta_Summary$quantiles), by = "VarSp") %>%
    dplyr::mutate(
      VarSp = purrr::map(
        .x = VarSp,
        .f = ~{
          .x %>%
            stringr::str_remove_all("B\\[|B\\[|\\(|\\)|\\]|stats::poly") %>%
            stringr::str_replace_all(", degree = 2, raw = TRUE", "_") %>%
            stringr::str_split(",", simplify = TRUE) %>%
            as.data.frame() %>%
            tibble::as_tibble() %>%
            stats::setNames(c("Variable", "ias_id")) %>%
            dplyr::mutate_all(stringr::str_trim) %>%
            dplyr::mutate(
              Variable = purrr::map_chr(
                .x = Variable,
                .f = function(V) {
                  V %>%
                    stringr::str_replace_all("_1$", "_L") %>%
                    stringr::str_replace_all("_2$", "_Q")
                }))
        })) %>%
    tidyr::unnest_wider(col = VarSp) %>%
    dplyr::mutate(
      CI_Overlap_0 = purrr::map2_lgl(
        .x = Q2_5, .y = Q97_5, ~dplyr::between(x = 0, left = .x, right = .y)))

  # Alpha -----
  if (("Alpha" %in% coda_names) && spatial_model) {
    ecokit::cat_time("Alpha", level = 1L)
    Alpha_Summary <- DataPrep(Alpha_Summary$statistics) %>%
      dplyr::full_join(DataPrep(Alpha_Summary$quantiles), by = "VarSp") %>%
      dplyr::mutate(
        VarSp = purrr::map(
          .x = VarSp,
          .f = ~{
            stringr::str_split(
              string = .x, pattern = "\\[", simplify = TRUE) %>%
              stringr::str_remove(pattern = "\\]") %>%
              as.vector() %>%
              stringr::str_trim() %>%
              stats::setNames(c("Alpha", "Factor"))
          })) %>%
      tidyr::unnest_wider(col = VarSp) %>%
      dplyr::mutate(
        CI_Overlap_0 = purrr::map2_lgl(
          .x = Q2_5, .y = Q97_5, ~dplyr::between(x = 0, left = .x, right = .y)))
  } else {
    Alpha_Summary <- tibble::tibble()
  }

  # Rho ----
  if ("Rho" %in% coda_names) {
    ecokit::cat_time("Rho", level = 1L)
    Rho_Summary <- dplyr::bind_rows(
      as.data.frame(as.matrix(Rho_Summary$statistics)),
      as.data.frame(as.matrix(Rho_Summary$quantiles))) %>%
      stats::setNames("Val") %>%
      tibble::rownames_to_column(var = "Var") %>%
      tibble::as_tibble() %>%
      tidyr::pivot_wider(names_from = "Var", values_from = Val) %>%
      dplyr::rename(
        dplyr::any_of(c(
          Naive_SE = "Naive SE", TimeSeries_SE = "Time-series SE",
          Q2_5 = "2.5%", Q25 = "25%", Q50 = "50%", Q75 = "75%", Q97_5 = "97.5%")
        )) %>%
      dplyr::mutate(Rho = "Taxonomy", .before = 1) %>%
      dplyr::mutate(  # nolint: consecutive_mutate_linter
        CI_Overlap_0 = purrr::map2_lgl(
          .x = Q2_5, .y = Q97_5, ~dplyr::between(x = 0, left = .x, right = .y)))
  } else {
    Rho_Summary <- tibble::tibble()
  }

  # Omega ------

  if ("Omega" %in% coda_names) {
    ecokit::cat_time("Omega", level = 1L)

    ListSp <- utils::read.delim(taxa_info_file, sep = "\t") %>%
      tibble::tibble() %>%
      dplyr::mutate(
        ias_id = stringr::str_pad(ias_id, pad = "0", width = 4),
        ias_id = paste0("Sp_", ias_id)) %>%
      dplyr::select(
        Sp_abb = ias_id, taxon_name, species_name, tidyselect::everything())

    # Prepare summary for the omega parameter
    Omega_Summary <- DataPrep(Omega_Summary$statistics) %>%
      dplyr::full_join(DataPrep(Omega_Summary$quantiles), by = "VarSp") %>%
      dplyr::mutate(
        VarSp = purrr::map(
          .x = VarSp,
          .f = ~{
            stringr::str_remove_all(.x, pattern = "\\]|Omega1|\\[") %>%
              stringr::str_split(pattern = ",", simplify = TRUE) %>%
              as.data.frame() %>%
              tibble::tibble() %>%
              stats::setNames(c("Sp1_abb", "Sp2_abb")) %>%
              dplyr::mutate_all(stringr::str_trim)
          })) %>%
      tidyr::unnest_wider(col = VarSp) %>%
      dplyr::left_join(
        dplyr::rename_all(ListSp, ~paste0("Sp1_", .x)),
        by = dplyr::join_by(Sp1_abb == Sp1_Sp_abb)) %>%
      dplyr::left_join(
        dplyr::rename_all(ListSp, ~paste0("Sp2_", .x)),
        by = dplyr::join_by(Sp2_abb == Sp2_Sp_abb)) %>%
      dplyr::select(
        Sp1_abb, Sp2_abb, tidyselect::everything(),
        tidyselect::starts_with("Sp1"), tidyselect::starts_with("Sp2"),
        -tidyselect::ends_with("_name2"), -tidyselect::ends_with("_File")) %>%
      dplyr::relocate(Sp1 = Sp1_Species_name, .after = "Sp1_abb") %>%
      dplyr::relocate(Sp2 = Sp2_Species_name, .after = "Sp2_abb") %>%
      dplyr::mutate(
        CI_Overlap_0 = purrr::map2_lgl(
          .x = Q2_5, .y = Q97_5, ~dplyr::between(x = 0, left = .x, right = .y)))
  } else {
    Omega_Summary <- tibble::tibble()
  }

  # # ..................................................................... ###

  # Saving ------
  ecokit::cat_time("Saving")
  Path_Out <- dirname(dirname(path_coda)) %>%
    fs::path("Model_Postprocessing", "Parameters_Summary")
  fs::dir_create(Path_Out)

  ecokit::save_as(
    object = list(
      Alpha = Alpha_Summary, Beta = Beta_Summary,
      Rho = Rho_Summary, Omega = Omega_Summary),
    object_name = "Parameters_Summary",
    out_path = fs::path(Path_Out, "Parameters_Summary.RData"))

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time, prefix = "Extracting model summary took ")

  # # ..................................................................... ###

  if (return_data) {
    return(
      list(
        Alpha = Alpha_Summary, Beta = Beta_Summary,
        Rho = Rho_Summary, Omega = Omega_Summary
      ))
  } else {
    return(invisible(NULL))
  }
}
