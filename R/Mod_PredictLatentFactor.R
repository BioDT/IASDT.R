## |------------------------------------------------------------------------| #
# predictLatentFactor ----
## |------------------------------------------------------------------------| #

#' predictLatentFactor
#'
#' Draws samples from the conditional predictive distribution of latent factors.
#' This function is adapted from [Hmsc::predictLatentFactor].
#'
#' @param unitsPred a factor vector with random level units for which
#'   predictions are to be made
#' @param modelunits a factor vector with random level units that are
#'   conditioned on
#' @param postEta a list containing samples of random factors at conditioned
#'   units
#' @param postAlpha a list containing samples of range (lengthscale) parameters
#'   for latent factors
#' @param NCores Integer specifying the number of cores to use for parallel
#'   processing. Defaults to 8.
#' @param TempDir Character string specifying the path for temporary storage of
#'   intermediate files.
#' @param ModelName Character string used as a prefix for temporary file names.
#' @param rL a HmscRandomLevel-class object that describes the random level
#'   structure
#' @export
#' @name predictLatentFactor

predictLatentFactor <- function(
    unitsPred, modelunits, postEta, postAlpha, rL, NCores = 8,
    TempDir = "TEMP2Pred", ModelName = NULL) {

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/

  SampleID <- Unit_ID <- LF <- LF_ID <- etaPred <- Sample_IDs <- File <-
    Alpha_ID <- NULL

  # # ..................................................................... ###

  fs::dir_create(TempDir)

  if (inherits(postEta, "character")) {
    postEta <- qs::qread(postEta, nthreads = 5)
  }

  indOld <- (unitsPred %in% modelunits)
  indNew <- !(indOld)

  # # ..................................................................... ###

  # Calculate D11 and D12 only once

  alphapw <- rL$alphapw
  s1 <- rL$s[modelunits, , drop = FALSE]
  s2 <- rL$s[unitsPred[indNew], , drop = FALSE]
  D11 <- Rfast::Dist(s1)
  D12 <- Rfast::dista(s1, s2)

  Path_D11 <- file.path(TempDir, paste0(ModelName, "_D11.qs"))
  qs::qsave(D11, file = Path_D11, preset = "fast")
  Path_D12 <- file.path(TempDir, paste0(ModelName, "_D12.qs"))
  qs::qsave(D12, file = Path_D12, preset = "fast")
  rm(D11, D12)
  invisible(gc())

  on.exit(
    try(fs::file_delete(c(Path_D11, Path_D12)), silent = TRUE),
    add = TRUE)

  # # ..................................................................... ###

  postAlpha_tibble <- do.call(rbind, postAlpha) %>%
    as.data.frame() %>%
    tibble::tibble() %>%
    stats::setNames(paste0("LF_", seq_len(ncol(.))))

  postAlpha_unique <- postAlpha_tibble %>%
    # ID column represents the original row number
    dplyr::mutate(Sample_IDs = dplyr::row_number()) %>%
    tidyr::nest(Sample_IDs = Sample_IDs) %>%
    dplyr::mutate(Sample_IDs = purrr::map(Sample_IDs, unlist))

  Unique_Alpha <- postAlpha_unique %>%
    dplyr::select(-Sample_IDs) %>%
    tidyr::pivot_longer(
      cols = names(.), values_to = "Alpha_ID", names_to = "LF") %>%
    dplyr::mutate(
      LF_ID = as.integer(stringr::str_remove(LF, "LF_")),
      .before = "Alpha_ID") %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      Denom = purrr::map_dbl(Alpha_ID, ~ alphapw[.x, 1]),
      SampleID = purrr::map2(
        .x = LF, .y = Alpha_ID,
        .f = ~ {
          postAlpha_unique %>%
            dplyr::filter(.[[.x]] == .y) %>%
            dplyr::pull(Sample_IDs) %>%
            unlist() %>%
            as.vector()
        }),
      SampleID = purrr::map(SampleID, unlist)) %>%
    dplyr::mutate(
      File = file.path(
        TempDir,
        paste0(ModelName, "_postEta_ch", dplyr::row_number(), ".qs")),
      Export = purrr::map2(
        .x = SampleID, .y = File,
        .f = ~ {
          Out <- postEta[.x]
          qs::qsave(Out, file = .y, preset = "fast")
        }
      ),
      Export = NULL
    )

  # # ..................................................................... ###

  calc_eta_pred <- function(RowNum) {
    Denom <- Unique_Alpha$Denom[[RowNum]]
    LF_ID <- Unique_Alpha$LF_ID[[RowNum]]
    SampleID <- Unique_Alpha$SampleID[[RowNum]]
    File <- Unique_Alpha$File[[RowNum]]

    postEta <- qs::qread(File, nthreads = 5)

    D11 <- qs::qread(Path_D11, nthreads = 5)
    D12 <- qs::qread(Path_D12, nthreads = 5)

    # Calculate K11 and K12 if Denom > 0, otherwise set `eta_indNew` to zero
    if (Denom > 0) {
      K11 <- IASDT.R::exp_neg_div(D11, Denom)
      K12 <- IASDT.R::exp_neg_div(D12, Denom)

      etaPred <- purrr::map_dfr(
        .x = seq_len(length(SampleID)),
        .f = function(ID) {
          eta <- postEta[[ID]][, LF_ID]
          eta_indNew <- IASDT.R::Solve2vect(K11, eta) %>%
            as.vector() %>%
            Matrix::crossprod(K12, .)

          etaPred <- matrix(NA, length(unitsPred), 1)

          etaPred[indOld, ] <- eta[match(unitsPred[indOld], modelunits)]
          etaPred[indNew, ] <- eta_indNew

          tibble::tibble(
            SampleID = SampleID[ID], etaPred = as.vector(etaPred),
            Unit_ID = seq_len(length(etaPred)))
        }) %>%
        tidyr::unnest("etaPred")
    } else {
      etaPred <- purrr::map_dfr(
        .x = seq_len(length(SampleID)),
        .f = function(ID) {
          eta <- postEta[[ID]][, LF_ID]
          etaPred <- matrix(NA, length(unitsPred), 1)

          etaPred[indOld, ] <- eta[match(unitsPred[indOld], modelunits)]
          etaPred[indNew, ] <- 0

          tibble::tibble(
            SampleID = SampleID[ID],
            etaPred = as.vector(etaPred),
            Unit_ID = seq_len(length(etaPred)))
        }) %>%
        tidyr::unnest("etaPred")
    }

    fs::file_delete(File)
    return(etaPred)
  }

  # # ..................................................................... ###

  withr::local_options(
    future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE)
  c1 <- snow::makeSOCKcluster(NCores)
  on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
  future::plan("future::cluster", workers = c1, gc = TRUE)
  on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)

  etaPreds <- future.apply::future_lapply(
    X = seq_len(nrow(Unique_Alpha)),
    FUN = calc_eta_pred,
    future.seed = TRUE,
    future.globals = c(
      "Unique_Alpha", "Path_D11", "Path_D12", "indNew", "unitsPred",
      "indOld", "modelunits"),
    future.packages = c(
      "Rcpp", "RcppArmadillo", "dplyr", "tidyr", "tibble",
      "Matrix", "Hmsc", "qs", "fs", "purrr"))

  postEtaPred <- Unique_Alpha %>%
    dplyr::mutate(etaPred = etaPreds) %>%
    dplyr::select(LF, LF_ID, etaPred) %>%
    tidyr::unnest("etaPred") %>%
    tidyr::pivot_wider(
      id_cols = c(SampleID, Unit_ID),
      names_from = LF, values_from = etaPred) %>%
    dplyr::arrange(SampleID, Unit_ID) %>%
    dplyr::group_split(SampleID) %>%
    purrr::map(
      .f = ~ {
        Mat <- dplyr::select(.x, tidyselect::starts_with("LF_")) %>%
          as.matrix()
        rownames(Mat) <- unitsPred
        colnames(Mat) <- NULL
        return(Mat)
      })

  return(postEtaPred)
}
