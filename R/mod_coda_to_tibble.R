## |------------------------------------------------------------------------| #
# coda_to_tibble ----
## |------------------------------------------------------------------------| #

#' Convert a Coda object to a tibble with specified parameter transformations
#'
#' This function converts a Coda object (`mcmc.list` or `mcmc`) into a tibble
#' format, facilitating further analysis and visualisation. It supports
#' transformation for specific parameter types: `rho`, `alpha`, `omega`, and
#' `beta`.
#' @param coda_object An object of class `mcmc.list` or `mcmc`, representing the
#'   MCMC output.
#' @param posterior_type Character. The parameter type to transform and extract.
#'   Must be one of `rho`, `alpha`, `omega`, or `beta`.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param n_omega Integer. The number of species to be sampled for the `Omega`
#'   parameter transformation. Defaults to 100.
#' @name coda_to_tibble
#' @author Ahmed El-Gabbas
#' @return A tibble containing the transformed parameters based on the specified
#'   `posterior_type`. The structure of the returned tibble varies depending on
#'   the `posterior_type` parameter.
#' @export
#' @examples
#' ecokit::load_packages(Hmsc, coda, dplyr, data.table)
#'
#' coda_object <- Hmsc::convertToCodaObject(Hmsc::TD$m)
#' ecokit::ht(
#'   IASDT.R::coda_to_tibble(
#'     coda_object = coda_object$Alpha[[1]], posterior_type = "Alpha"))

coda_to_tibble <- function(
    coda_object = NULL, posterior_type = NULL,
    env_file = ".env", n_omega = 100) {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  chain <- iter <- Alpha <- alpha_num <- factor <- CHAIN <- ITER <-
    value <- ias_id <- species_name <- var_sp <- variable <- species <-
    species_combs <- sp1 <- naps_1 <- sp2 <- naps_2 <- taxa_info_file <- NULL

  # # |||||||||||||||||||||||||||||||||||||||
  # Check missing arguments ----
  # # |||||||||||||||||||||||||||||||||||||||

  if (is.null(coda_object) || is.null(posterior_type)) {
    ecokit::stop_ctx(
      "None of `coda_object` or `posterior_type` cannot be empty",
      coda_object = coda_object, posterior_type = posterior_type,
      include_backtrace = TRUE)
  }

  # # |||||||||||||||||||||||||||||||||||||||
  # Check class of input coda object -----
  # # |||||||||||||||||||||||||||||||||||||||

  posterior_type <- tolower(posterior_type)

  if (!posterior_type %in% c("rho", "alpha", "omega", "beta")) {
    ecokit::stop_ctx(
      "posterior_type has to be one of rho, alpha, omega, or beta",
      posterior_type = posterior_type, include_backtrace = TRUE)
  }

  if (!(inherits(coda_object, "mcmc.list") || inherits(coda_object, "mcmc"))) {
    ecokit::stop_ctx(
      "Input Coda object has to be of class mcmc.list or mcmc",
      class_coda_obj = class(coda_object), include_backtrace = TRUE)
  }

  # # |||||||||||||||||||||||||||||||||||||||
  # Convert to tibble ----
  # # |||||||||||||||||||||||||||||||||||||||

  if (posterior_type == "omega") {
    CombSample <- sample.int(n = dim(coda_object[[1]])[2], size = n_omega)      # nolint: object_name_linter

    coda_data <- purrr::map(coda_object, ~.x[, CombSample]) %>%
      coda::as.mcmc.list() %>%
      as.matrix(iter = TRUE, chain = TRUE) %>%
      tibble::as_tibble() %>%
      dplyr::rename(chain = CHAIN, iter = ITER) %>%
      dplyr::mutate(chain = factor(chain), iter = as.integer(iter)) %>%
      dplyr::arrange(chain, iter) %>%
      tidyr::pivot_longer(
        -c(chain, iter), names_to = "species_combs", values_to = "value")
  } else {
    coda_data <- as.matrix(coda_object, iter = TRUE, chain = TRUE) %>%
      tibble::as_tibble() %>%
      dplyr::rename(chain = CHAIN, iter = ITER) %>%
      dplyr::mutate(chain = factor(chain), iter = as.integer(iter)) %>%
      dplyr::arrange(chain, iter)
  }


  # # |||||||||||||||||||||||||||||||||||||||
  # Rho parameter ----
  # # |||||||||||||||||||||||||||||||||||||||

  if (posterior_type == "rho") {
    coda_data <- dplyr::rename(coda_data, dplyr::any_of(c(value = "var1")))
  }

  # # |||||||||||||||||||||||||||||||||||||||
  # Alpha parameter ----
  # # |||||||||||||||||||||||||||||||||||||||

  if (posterior_type == "alpha") {

    coda_data <- tidyr::pivot_longer(
      data = coda_data,
      cols = -c(chain, iter), names_to = "Alpha", values_to = "value")

    coda_data <- dplyr::distinct(coda_data, Alpha) %>%
      dplyr::mutate(
        Alpha2 = purrr::map(
          .x = Alpha,
          .f = ~{
            stringr::str_split(.x, "\\[", simplify = TRUE) %>%
              as.character() %>%
              stringr::str_remove("\\]")
          })) %>%
      tidyr::unnest_wider("Alpha2", names_sep = "_") %>%
      purrr::set_names(c("Alpha", "alpha_num", "factor")) %>%
      dplyr::right_join(coda_data, by = "Alpha") %>%
      dplyr::mutate(
        Alpha = factor(Alpha), alpha_num = factor(alpha_num),
        factor = factor(factor)) %>%
      dplyr::select(
        Alpha, alpha_num, factor, chain, iter, value, dplyr::everything()) %>%
      dplyr::arrange(Alpha, chain, iter)
  }

  # # |||||||||||||||||||||||||||||||||||||||
  # Beta parameters ------
  # # |||||||||||||||||||||||||||||||||||||||

  if (posterior_type == "beta") {
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

    sp_names <- readr::read_tsv(
      file = taxa_info_file, show_col_types = FALSE, progress = FALSE) %>%
      dplyr::select(ias_id, species = species_name) %>%
      dplyr::mutate(
        ias_id = stringr::str_pad(ias_id, pad = "0", width = 4),
        ias_id = paste0("sp_", ias_id))

    coda_data <- tidyr::pivot_longer(
      coda_data, -c(chain, iter), names_to = "var_sp", values_to = "value")


    var_species <- coda_data %>%
      dplyr::distinct(var_sp) %>%
      dplyr::mutate(
        species = purrr::map(
          .x = var_sp,
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
      tidyr::unnest_wider("species") %>%
      dplyr::left_join(sp_names, by = "ias_id")

    coda_data <- var_species %>%
      dplyr::right_join(coda_data, by = "var_sp") %>%
      dplyr::select(
        var_sp, variable, species, ias_id, chain, iter, value,
        dplyr::everything()) %>%
      tidyr::nest(data = -c(variable, ias_id, species, var_sp)) %>%
      dplyr::arrange(variable, ias_id)
  }

  # # |||||||||||||||||||||||||||||||||||||||
  # Omega parameters ------
  # # |||||||||||||||||||||||||||||||||||||||

  if (posterior_type == "omega") {

    naps <- IASDT.R::get_species_name(env_file = env_file) %>%      # nolint: object_name_linter
      dplyr::select(ias_id, species_name)

    coda_data <- tibble::tibble(
      species_combs = unique(coda_data$species_combs)) %>%
      dplyr::mutate(
        SP = purrr::map(
          .x = species_combs,
          .f = ~{
            naps_n <- stringr::str_remove_all(.x, "Omega1\\[|\\]") %>%
              stringr::str_split(",| ", simplify = TRUE) %>%
              as.character() %>%
              stringr::str_subset("^sp_") %>%
              purrr::set_names(c("sp1", "sp2"))
            naps_1 <- dplyr::filter(naps, ias_id == naps_n[1])$species_name
            naps_2 <- dplyr::filter(naps, ias_id == naps_n[2])$species_name

            return(c(naps_n, naps_1 = naps_1, naps_2 = naps_2))
          })) %>%
      tidyr::unnest_wider("SP") %>%
      dplyr::right_join(coda_data, by = "species_combs") %>%
      dplyr::select(
        species_combs, sp1, naps_1, sp2, naps_2, chain, iter, value,
        dplyr::everything()) %>%
      tidyr::nest(data = -c(species_combs, sp1, naps_1, sp2, naps_2)) %>%
      dplyr::arrange(species_combs)

  }
  return(coda_data)
}
