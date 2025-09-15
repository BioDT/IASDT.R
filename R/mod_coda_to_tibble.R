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
#' ecokit::load_packages(Hmsc, coda, dplyr)
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
  Chain <- Iter <- Alpha <- AlphaNum <- Factor <- CHAIN <- ITER <-
    Value <- ias_id <- species_name <- Var_Sp <- Variable <- Species <-
    species_combs <- Sp1 <- IAS1 <- Sp2 <- IAS2 <- taxa_info_file <- NULL

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

    Coda <- purrr::map(coda_object, ~.x[, CombSample]) %>%
      coda::as.mcmc.list() %>%
      as.matrix(iter = TRUE, chain = TRUE) %>%
      tibble::as_tibble() %>%
      dplyr::rename(Chain = CHAIN, Iter = ITER) %>%
      dplyr::mutate(Chain = factor(Chain), Iter = as.integer(Iter)) %>%
      dplyr::arrange(Chain, Iter) %>%
      tidyr::pivot_longer(
        -c(Chain, Iter), names_to = "species_combs", values_to = "Value")
  } else {
    Coda <- as.matrix(coda_object, iter = TRUE, chain = TRUE) %>%
      tibble::as_tibble() %>%
      dplyr::rename(Chain = CHAIN, Iter = ITER) %>%
      dplyr::mutate(Chain = factor(Chain), Iter = as.integer(Iter)) %>%
      dplyr::arrange(Chain, Iter)
  }


  # # |||||||||||||||||||||||||||||||||||||||
  # Rho parameter ----
  # # |||||||||||||||||||||||||||||||||||||||

  if (posterior_type == "rho") {
    Coda <- dplyr::rename(Coda, dplyr::any_of(c(Value = "var1")))
  }

  # # |||||||||||||||||||||||||||||||||||||||
  # Alpha parameter ----
  # # |||||||||||||||||||||||||||||||||||||||

  if (posterior_type == "alpha") {

    Coda <- tidyr::pivot_longer(
      data = Coda,
      cols = -c(Chain, Iter), names_to = "Alpha", values_to = "Value")

    Coda <- dplyr::distinct(Coda, Alpha) %>%
      dplyr::mutate(
        Alpha2 = purrr::map(
          .x = Alpha,
          .f = ~{
            stringr::str_split(.x, "\\[", simplify = TRUE) %>%
              as.character() %>%
              stringr::str_remove("\\]")
          })) %>%
      tidyr::unnest_wider("Alpha2", names_sep = "_") %>%
      purrr::set_names(c("Alpha", "AlphaNum", "Factor")) %>%
      dplyr::right_join(Coda, by = "Alpha") %>%
      dplyr::mutate(
        Alpha = factor(Alpha), AlphaNum = factor(AlphaNum),
        Factor = factor(Factor)) %>%
      dplyr::select(
        Alpha, AlphaNum, Factor, Chain, Iter, Value, dplyr::everything()) %>%
      dplyr::arrange(Alpha, Chain, Iter)
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

    SpeciesNames <- readr::read_tsv(
      file = taxa_info_file, show_col_types = FALSE, progress = FALSE) %>%
      dplyr::select(ias_id, Species = species_name) %>%
      dplyr::mutate(
        ias_id = stringr::str_pad(ias_id, pad = "0", width = 4),
        ias_id = paste0("Sp_", ias_id))

    Coda <- tidyr::pivot_longer(
      Coda, -c(Chain, Iter), names_to = "Var_Sp", values_to = "Value")


    VarSp <- Coda %>%
      dplyr::distinct(Var_Sp) %>%
      dplyr::mutate(
        Species = purrr::map(
          .x = Var_Sp,
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
      tidyr::unnest_wider("Species") %>%
      dplyr::left_join(SpeciesNames, by = "ias_id")



    Coda <- VarSp %>%
      dplyr::right_join(Coda, by = "Var_Sp") %>%
      dplyr::select(
        Var_Sp, Variable, Species, ias_id, Chain, Iter, Value,
        dplyr::everything()) %>%
      tidyr::nest(DT = -c(Variable, ias_id, Species, Var_Sp)) %>%
      dplyr::arrange(Variable, ias_id)
  }

  # # |||||||||||||||||||||||||||||||||||||||
  # Omega parameters ------
  # # |||||||||||||||||||||||||||||||||||||||

  if (posterior_type == "omega") {

    IAS <- IASDT.R::get_species_name(env_file = env_file) %>%      # nolint: object_name_linter
      dplyr::select(ias_id, species_name)

    Coda <- tibble::tibble(species_combs = unique(Coda$species_combs)) %>%
      dplyr::mutate(
        SP = purrr::map(
          .x = species_combs,
          .f = ~{
            IAS_N <- stringr::str_remove_all(.x, "Omega1\\[|\\]") %>%
              stringr::str_split(",| ", simplify = TRUE) %>%
              as.character() %>%
              stringr::str_subset("^Sp_") %>%
              purrr::set_names(c("Sp1", "Sp2"))
            IAS1 <- dplyr::filter(IAS, ias_id == IAS_N[1])$species_name
            IAS2 <- dplyr::filter(IAS, ias_id == IAS_N[2])$species_name

            return(c(IAS_N, IAS1 = IAS1, IAS2 = IAS2))
          })) %>%
      tidyr::unnest_wider("SP") %>%
      dplyr::right_join(Coda, by = "species_combs") %>%
      dplyr::select(
        species_combs, Sp1, IAS1, Sp2, IAS2, Chain, Iter, Value,
        dplyr::everything()) %>%
      tidyr::nest(DT = -c(species_combs, Sp1, IAS1, Sp2, IAS2)) %>%
      dplyr::arrange(species_combs)

  }
  return(Coda)
}
