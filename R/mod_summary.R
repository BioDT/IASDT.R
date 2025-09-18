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
  var_species <- sp1_abb <- sp2_abb <- ias_id <- value <- taxon_name <-
    species_name <- sp1_species_name <- sp1_species_name <- q_2_5 <- q_97_5 <-
    sp1_sp_abb <- sp2_sp_abb <- taxa_info_file <- NULL

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

  # data_prep helper function -------
  data_prep <- function(data) {
    as.data.frame(data) %>%
      tibble::rownames_to_column(var = "var_species") %>%
      tibble::as_tibble() %>%
      dplyr::rename(
        dplyr::any_of(c(
          Naive_SE = "Naive SE", TimeSeries_SE = "Time-series SE",
          q_2_5 = "2.5%", q_25 = "25%", q_50 = "50%", q75 = "75%",
          q_97_5 = "97.5%")))
  }

  # # ..................................................................... ###

  # Loading coda object ------
  ecokit::cat_time("Loading coda object")
  coda_obj <- ecokit::load_as(path_coda)
  coda_names <- names(coda_obj)

  ecokit::cat_time("Beta summary", level = 1L)
  beta_summary <- summary(coda_obj$Beta)

  if (("Alpha" %in% coda_names) && spatial_model) {
    ecokit::cat_time("Alpha summary", level = 1L)
    alpha_summary <- summary(coda_obj$Alpha[[1]])
  } else {
    alpha_summary <- NULL
  }

  if ("Omega" %in% coda_names) {
    ecokit::cat_time("Omega summary", level = 1L)
    omega_summary <- summary(coda_obj$Omega[[1]])
  } else {
    omega_summary <- NULL
  }

  if ("Rho" %in% coda_names) {
    ecokit::cat_time("Rho summary", level = 1L)
    rho_summary <- summary(coda_obj$Rho)
  } else {
    rho_summary <- NULL
  }

  rm(coda_obj, envir = environment())
  invisible(gc())

  # # ..................................................................... ###

  ecokit::cat_time("Extracting summary data")

  # Beta ------
  ecokit::cat_time("Beta", level = 1L)
  beta_summary <- data_prep(beta_summary$statistics) %>%
    dplyr::full_join(data_prep(beta_summary$quantiles), by = "var_species") %>%
    dplyr::mutate(
      var_species = purrr::map(
        .x = var_species,
        .f = ~{
          .x %>%
            stringr::str_remove_all("B\\[|B\\[|\\(|\\)|\\]|stats::poly") %>%
            stringr::str_replace_all(", degree = 2, raw = TRUE", "_") %>%
            stringr::str_split(",", simplify = TRUE) %>%
            as.data.frame() %>%
            tibble::as_tibble() %>%
            stats::setNames(c("variable", "ias_id")) %>%
            dplyr::mutate_all(stringr::str_trim) %>%
            dplyr::mutate(
              variable = purrr::map_chr(
                .x = variable,
                .f = function(v) {
                  v %>%
                    stringr::str_replace_all("_1$", "_L") %>%
                    stringr::str_replace_all("_2$", "_Q")
                }))
        })) %>%
    tidyr::unnest_wider(col = var_species) %>%
    dplyr::mutate(
      CI_Overlap_0 = purrr::map2_lgl(
        .x = q_2_5, .y = q_97_5, ~dplyr::between(x = 0, left = .x, right = .y)))

  # Alpha -----
  if (("Alpha" %in% coda_names) && spatial_model) {
    ecokit::cat_time("Alpha", level = 1L)
    alpha_summary <- data_prep(alpha_summary$statistics) %>%
      dplyr::full_join(
        data_prep(alpha_summary$quantiles), by = "var_species") %>%
      dplyr::mutate(
        var_species = purrr::map(
          .x = var_species,
          .f = ~{
            stringr::str_split(
              string = .x, pattern = "\\[", simplify = TRUE) %>%
              stringr::str_remove(pattern = "\\]") %>%
              as.vector() %>%
              stringr::str_trim() %>%
              stats::setNames(c("Alpha", "Factor"))
          })) %>%
      tidyr::unnest_wider(col = var_species) %>%
      dplyr::mutate(
        CI_Overlap_0 = purrr::map2_lgl(
          .x = q_2_5, .y = q_97_5,
          ~dplyr::between(x = 0, left = .x, right = .y)))
  } else {
    alpha_summary <- tibble::tibble()
  }

  # Rho ----
  if ("Rho" %in% coda_names) {
    ecokit::cat_time("Rho", level = 1L)
    rho_summary <- dplyr::bind_rows(
      as.data.frame(as.matrix(rho_summary$statistics)),
      as.data.frame(as.matrix(rho_summary$quantiles))) %>%
      stats::setNames("value") %>%
      tibble::rownames_to_column(var = "Var") %>%
      tibble::as_tibble() %>%
      tidyr::pivot_wider(names_from = "Var", values_from = value) %>%
      dplyr::rename(
        dplyr::any_of(c(
          Naive_SE = "Naive SE", TimeSeries_SE = "Time-series SE",
          q_2_5 = "2.5%", q_25 = "25%", q_50 = "50%",
          q75 = "75%", q_97_5 = "97.5%")
        )) %>%
      dplyr::mutate(Rho = "Taxonomy", .before = 1) %>%
      dplyr::mutate(  # nolint: consecutive_mutate_linter
        CI_Overlap_0 = purrr::map2_lgl(
          .x = q_2_5, .y = q_97_5,
          ~dplyr::between(x = 0, left = .x, right = .y)))
  } else {
    rho_summary <- tibble::tibble()
  }

  # Omega ------

  if ("Omega" %in% coda_names) {
    ecokit::cat_time("Omega", level = 1L)

    list_sp <- utils::read.delim(taxa_info_file, sep = "\t") %>%
      tibble::tibble() %>%
      dplyr::mutate(
        ias_id = stringr::str_pad(ias_id, pad = "0", width = 4),
        ias_id = paste0("sp_", ias_id)) %>%
      dplyr::select(
        sp_abb = ias_id, taxon_name, species_name, tidyselect::everything())

    # Prepare summary for the omega parameter
    omega_summary <- data_prep(omega_summary$statistics) %>%
      dplyr::full_join(
        data_prep(omega_summary$quantiles), by = "var_species") %>%
      dplyr::mutate(
        var_species = purrr::map(
          .x = var_species,
          .f = ~{
            stringr::str_remove_all(.x, pattern = "\\]|Omega1|\\[") %>%
              stringr::str_split(pattern = ",", simplify = TRUE) %>%
              as.data.frame() %>%
              tibble::tibble() %>%
              stats::setNames(c("sp1_abb", "sp2_abb")) %>%
              dplyr::mutate_all(stringr::str_trim)
          })) %>%
      tidyr::unnest_wider(col = var_species) %>%
      dplyr::left_join(
        dplyr::rename_all(list_sp, ~paste0("sp1_", .x)),
        by = dplyr::join_by(sp1_abb == sp1_sp_abb)) %>%
      dplyr::left_join(
        dplyr::rename_all(list_sp, ~paste0("sp2_", .x)),
        by = dplyr::join_by(sp2_abb == sp2_sp_abb)) %>%
      dplyr::select(
        sp1_abb, sp2_abb, tidyselect::everything(),
        tidyselect::starts_with("sp1"), tidyselect::starts_with("sp2"),
        -tidyselect::ends_with("_name2"), -tidyselect::ends_with("_File")) %>%
      dplyr::relocate(sp1 = sp1_species_name, .after = "sp1_abb") %>%
      dplyr::relocate(sp2 = sp1_species_name, .after = "sp2_abb") %>%
      dplyr::mutate(
        CI_Overlap_0 = purrr::map2_lgl(
          .x = q_2_5, .y = q_97_5,
          ~dplyr::between(x = 0, left = .x, right = .y)))
  } else {
    omega_summary <- tibble::tibble()
  }

  # # ..................................................................... ###

  # Saving ------
  ecokit::cat_time("Saving")
  path_out <- dirname(dirname(path_coda)) %>%
    fs::path("model_postprocessing", "parameters_summary")
  fs::dir_create(path_out)

  ecokit::save_as(
    object = list(
      Alpha = alpha_summary, Beta = beta_summary,
      Rho = rho_summary, Omega = omega_summary),
    object_name = "parameters_summary",
    out_path = fs::path(path_out, "parameters_summary.RData"))

  # # ..................................................................... ###

  ecokit::cat_diff(
    init_time = .start_time, prefix = "Extracting model summary took ")

  # # ..................................................................... ###

  if (return_data) {
    return(
      list(
        Alpha = alpha_summary, Beta = beta_summary,
        Rho = rho_summary, Omega = omega_summary
      ))
  } else {
    return(invisible(NULL))
  }
}
