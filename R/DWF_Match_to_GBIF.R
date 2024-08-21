## |------------------------------------------------------------------------| #
# Match_to_GBIF ----
## |------------------------------------------------------------------------| #

#' Match taxonomy names with GBIF and return standardized results.
#'
#' This function matches given taxonomy names with the GBIF database,
#' potentially returning more than one match per name. It can operate in
#' parallel.
#' @param taxon_name A character vector of taxonomy names to be matched.
#' @param taxon_id An optional numeric or character vector of local identifiers
#'   corresponding to `taxon_name`. If `NULL`, a sequence along `taxon_name` is
#'   used. Defaults to `NULL`.
#' @param include_genus Logical, whether to include matches at the genus level
#'   in the results. Defaults to `FALSE`.
#' @param Parallel Logical, whether to perform the matching in parallel.
#'   Defaults to `FALSE`.
#' @param Progress Logical, whether to display a progress bar during the
#'   operation. Defaults to `FALSE`.
#' @name Match_to_GBIF
#' @importFrom rlang .data
#' @author Marina Golivets
#' @return A `tibble` containing the standardized taxonomy results, including
#'   both the best matches and alternatives, filtered to include only vascular
#'   plants (Tracheophyta) and excluding non-matches and higher rank matches.
#'   The results are further refined based on the confidence score and the
#'   status (ACCEPTED, SYNONYM, DOUBTFUL) of the matches.
#' @export
#' @details as input, provide a vector of verbatim taxon names (preferably with
#' authorship) and a vector of existing local identifiers for those names

Match_to_GBIF <- function(
    taxon_name, taxon_id = NULL, include_genus = FALSE,
    Parallel = FALSE, Progress = FALSE) {

  if (is.null(taxon_name)) {
    stop("taxon_name cannot be NULL", .call = FALSE)
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  phylum <- NULL

  if (is.null(taxon_id)) taxon_id <- seq_along(taxon_name)

  if (Parallel) {
    all_matches <- furrr::future_map(
      taxon_name, rgbif::name_backbone_verbose,
      kingdom = "plants", strict = TRUE,
      .progress = Progress, .options = furrr::furrr_options(seed = TRUE))
  } else {
    if (Progress) {
      ProgrOptns <- list(
        type = "iterator", clear = TRUE,
        format = "{cli::pb_bar} {cli::pb_percent} [{cli::pb_elapsed}]")
    } else {
      ProgrOptns <- FALSE
    }

    all_matches <- purrr::map(
      taxon_name, rgbif::name_backbone_verbose,
      kingdom = "plants", strict = TRUE, .progress = ProgrOptns)
  }

  # retrieve alternative matches
  alternative_matches <- lapply(
    all_matches,
    function(x) {
      y <- x$alternatives
      if (nrow(y) == 0) {
        y[1, 1] <- NA
        colnames(y) <- "usageKey"
      } else {
        y <- y
      }
      return(y)
    }
  ) %>%
    mapply(
      cbind, .,
      taxon_name = taxon_name, taxon_id = taxon_id,
      stringsAsFactors = FALSE, SIMPLIFY = FALSE
    ) %>%
    data.table::rbindlist(fill = TRUE) %>%
    dplyr::filter(!is.na(.data$usageKey)) %>%
    dplyr::distinct() %>%
    # filter only if phylum column exists
    dplyr::filter(
      (if ("phylum" %in% names(.)) phylum else NULL) == "Tracheophyta")

  # retrieve best matches
  best_matches <- lapply(all_matches, function(x) x$data) %>%
    mapply(
      cbind, .,
      taxon_name = taxon_name, taxon_id = taxon_id,
      stringsAsFactors = FALSE, SIMPLIFY = FALSE
    ) %>%
    data.table::rbindlist(fill = TRUE) %>%
    dplyr::distinct() %>%
    dplyr::filter(
      (if ("phylum" %in% names(.)) phylum else NULL) == "Tracheophyta")

  matched <- best_matches %>%
    dplyr::filter(!(.data$matchType %in% c("NONE", "HIGHERRANK")))

  matched_alternative <- try(
    alternative_matches %>%
      # use only vascular plants
      # filter only if phylum column exists
      dplyr::filter(
        (if ("phylum" %in% names(.)) phylum else NULL) == "Tracheophyta") %>%
      dplyr::filter(.data$confidence >= 0) %>%
      dplyr::filter(!taxon_id %in% matched$taxon_id)
  )
  if (class(matched_alternative)[1] == "try-error") {
    taxon_list <- matched
  } else {
    taxon_list <- dplyr::bind_rows(matched, matched_alternative)
  }

  if (include_genus == FALSE) {
    taxon_list <- taxon_list %>%
      dplyr::filter(rank != "GENUS")
  }

  # get names that were matched as accepted
  accepted <- taxon_list %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::filter(.data$status == "ACCEPTED")
  if (nrow(accepted) > 0) {
    accepted <- accepted %>%
      dplyr::filter(.data$confidence == max(.data$confidence)) %>%
      dplyr::ungroup()
  } else {
    accepted <- dplyr::ungroup(accepted)
  }

  # get names that were matched as synonyms only
  synonyms <- taxon_list %>%
    dplyr::group_by(.data$taxon_id) %>%
    dplyr::summarise(
      has_accepted = dplyr::n_distinct(.data$status == "ACCEPTED") > 1) %>%
    dplyr::full_join(taxon_list, by = "taxon_id") %>%
    dplyr::filter(.data$has_accepted == FALSE) %>%
    dplyr::filter(.data$status == "SYNONYM")
  if (nrow(synonyms) > 0) {
    synonyms <- synonyms %>%
      dplyr::group_by(taxon_id) %>%
      dplyr::filter(.data$confidence == max(.data$confidence)) %>%
      dplyr::ungroup()
  } else {
    synonyms <- dplyr::ungroup(synonyms)
  }

  # get names that were matched as doubtful only
  doubtful <- taxon_list %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::summarise(
      has_accepted = dplyr::n_distinct(.data$status == "ACCEPTED") > 1) %>%
    dplyr::full_join(taxon_list, by = "taxon_id") %>%
    dplyr::filter(.data$has_accepted == FALSE) %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::filter(.data$status == "DOUBTFUL")
  if (nrow(doubtful) > 0) {
    doubtful <- doubtful %>%
      dplyr::filter(.data$confidence == max(.data$confidence)) %>%
      dplyr::ungroup()
  } else {
    doubtful <- dplyr::ungroup(doubtful)
  }

  # combine all names
  taxon_list_final <- dplyr::bind_rows(accepted, synonyms, doubtful) %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::filter(.data$confidence == max(.data$confidence)) %>%
    dplyr::filter(.data$status != "NONE") %>% # exclude non-matched names
    dplyr::select(-"has_accepted") %>%
    dplyr::ungroup() %>%
    dplyr::relocate(taxon_name, taxon_id)

  return(taxon_list_final)
}
