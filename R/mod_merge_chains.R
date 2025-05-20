## |------------------------------------------------------------------------| #
# Merge model chains ----
## |------------------------------------------------------------------------| #

#' Merge model chains into `Hmsc` and `coda` objects
#'
#' These functions merge posterior chains from multiple runs of `Hmsc` models
#' into unified `Hmsc` and `coda` objects, facilitating further analysis. They
#' check for missing or incomplete chains, optionally report these issues, and
#' save the processed results to disk. `mod_merge_chains` handles regular
#' models, while `mod_merge_chains_CV` is designed for cross-validated models.
#'
#' @param model_dir Character. Path to the root directory where the model was
#'   fitted. For `mod_merge_chains`, subdirectories `Model_Fitted` and
#'   `Model_Coda` are created within this path to store the merged `Hmsc` and
#'   `coda` objects, respectively. For `mod_merge_chains_CV`, merged objects are
#'   stored under `Model_Fitting_CV/Model_Fitted`.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#'   Defaults to 8L.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "sequential", "multisession" (default), "multicore", and
#'   "cluster". See [future::plan()] and [ecokit::set_parallel()] for details.
#' @param model_info_name Character. Name of the file (without extension) where
#'   updated model information is saved. If `NULL`, overwrites the existing
#'   `Model_Info.RData` file in `model_dir` directory. If specified, creates a
#'   new `.RData` file with this name in `model_dir` directory. Applies only to
#'   `mod_merge_chains`.
#' @param print_incomplete Logical. If `TRUE`, prints the names of model
#'   variants that failed to merge due to missing or incomplete chains. Defaults
#'   to `TRUE`.
#' @param from_JSON Logical. Whether to convert loaded models from JSON format
#'   before reading. Defaults to `FALSE`.
#' @param out_extension Character. File extension (without dot) for output files
#'   containing merged `Hmsc` and `coda` objects. Options are `qs2` (faster
#'   read/write via the `qs2` package) or `RData` (standard R format). Defaults
#'   to `qs2`.
#' @param CV_names Character vector. Names of cross-validation strategies to
#'   merge, matching those used during model setup. Defaults to `c("CV_Dist",
#'   "CV_Large")`. The names should be one of `CV_Dist`, `CV_Large`, or
#'   `CV_SAC`. Applies only to `mod_merge_chains_CV`.
#' @export
#' @rdname mod_merge_chains
#' @name mod_merge_chains
#' @order 1
#' @author Ahmed El-Gabbas
#' @return Both functions return `invisible(NULL)` and save processed model
#'   information and merged objects to disk in the specified locations.
#' @details
#' - `mod_merge_chains` merges posterior chains from multiple runs of an `Hmsc`
#' model fitted without cross-validation. It checks for missing or incomplete
#' chains, aligns posteriors (using `alignPost = TRUE`, falling back to `FALSE`
#' if alignment fails), and saves a merged `Hmsc` object and a `coda` object for
#' MCMC diagnostics. It also records fitting time and memory usage from progress
#' files.
#' - `mod_merge_chains_CV` performs a similar merging process for
#' cross-validated `Hmsc` models, processing each fold of the specified
#' `CV_names` separately. It saves merged `Hmsc` objects per fold but does not
#' generate `coda` objects.

## |------------------------------------------------------------------------| #
# mod_merge_chains ----
## |------------------------------------------------------------------------| #

mod_merge_chains <- function(
    model_dir = NULL, n_cores = 8L, strategy = "multisession",
    model_info_name = NULL, print_incomplete = TRUE, from_JSON = FALSE,
    out_extension = "qs2") {

  # # ..................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Post_Path <- Post_Path <- M_Name_Fit <- Path_FittedMod <- Path_Coda <-
    NMissingChains <- MissingModels <- Model_Finished <- Path_ModProg <-
    Post_Aligned2 <- FittingTime <- FittingMemory <- NULL

  # # ..................................................................... ###

  # Checking arguments ----

  if (is.null(model_dir)) {
    ecokit::stop_ctx(
      "`model_dir` cannot be empty", model_dir = model_dir,
      include_backtrace = TRUE)
  }

  if (is.null(n_cores)) {
    ecokit::stop_ctx(
      "`n_cores` cannot be empty", n_cores = n_cores, include_backtrace = TRUE)
  }

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("model_dir", "out_extension", "strategy"))

  ecokit::check_args(
    args_all = AllArgs, args_to_check = "n_cores", args_type = "numeric")

  ecokit::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c("print_incomplete", "from_JSON"))

  rm(AllArgs, envir = environment())

  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores <= 0) {
    ecokit::stop_ctx(
      "n_cores must be a single positive integer.", n_cores = n_cores,
      include_backtrace = TRUE)
  }

  if (strategy == "sequential") {
    n_cores <- 1L
  }
  if (length(strategy) != 1L) {
    ecokit::stop_ctx(
      "`strategy` must be a character vector of length 1",
      strategy = strategy, length_strategy = length(strategy))
  }
  valid_strategy <- c("sequential", "multisession", "multicore", "cluster")
  if (!strategy %in% valid_strategy) {
    ecokit::stop_ctx("Invalid `strategy` value", strategy = strategy)
  }

  if (!dir.exists(model_dir)) {
    ecokit::stop_ctx(
      "`model_dir` directory does not exist", model_dir = model_dir,
      include_backtrace = TRUE)
  }

  if (length(out_extension) != 1) {
    ecokit::stop_ctx(
      "`out_extension` must be a single string.",
      out_extension = out_extension,
      length_out_extension = length(out_extension), include_backtrace = TRUE)
  }

  if (!out_extension %in% c("qs2", "RData")) {
    ecokit::stop_ctx(
      "`out_extension` must be either 'qs2' or 'RData'.",
      out_extension = out_extension, include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c("Hmsc", "coda", "purrr", "ecokit", "dplyr", "fs"),
    strategy = strategy)

  # # ..................................................................... ###

  # Creating paths -----

  Path_Fitted_Models <- fs::path(model_dir, "Model_Fitted")
  Path_Coda <- fs::path(model_dir, "Model_Coda")
  fs::dir_create(c(Path_Fitted_Models, Path_Coda))

  # # ..................................................................... ###

  # Loading model info ----

  Path_ModInfo <- fs::path(model_dir, "Model_Info.RData")

  if (!file.exists(Path_ModInfo)) {
    ecokit::stop_ctx(
      "ModInfo file does not exist", Path_ModInfo = Path_ModInfo,
      include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Remove temp files and incomplete RDs files ----

  Path_Model_Fit <- fs::path(model_dir, "Model_Fitting_HPC")
  tempFiles <- list.files(
    path = Path_Model_Fit, pattern = ".rds_temp$", full.names = TRUE)

  if (length(tempFiles) > 0) {
    ecokit::cat_time(
      paste0(
        "There are ", length(tempFiles),
        " unsuccessful model variants to be removed"))

    tempFilesRDs <- stringr::str_replace_all(tempFiles, ".rds_temp$", ".rds")
    tempFilesProgress <- stringr::str_replace_all(
      tempFiles, "_post.rds_temp", "_Progress.txt")

    purrr::walk(
      .x = c(tempFilesRDs, tempFiles, tempFilesProgress),
      .f = ~{
        if (file.exists(.x)) {
          file.remove(.x)
        }
      })
  }

  # # ..................................................................... ###

  Model_Info2 <- ecokit::load_as(Path_ModInfo)

  # Prepare working in parallel -----
  if (n_cores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = min(n_cores, nrow(Model_Info2)),
      future_max_size = 800L, strategy = strategy, level = 1L)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  ecokit::cat_time("Check if any posterior files is missing")
  # Check if any posterior files is missing
  Model_Info2 <- Model_Info2 %>%
    dplyr::mutate(
      Post_Missing = furrr::future_map_lgl(
        .x = Post_Path,
        .f = function(x) {

          purrr::map_lgl(
            .x = as.character(x),
            .f = function(y) {

              if (isFALSE(fs::file_exists(y))) {
                return(TRUE)
              } else if (ecokit::check_rds(y)) {
                return(FALSE)
              } else {
                fs::file_delete(y)
                return(TRUE)
              }

            }) %>%
            any()
        }),
      # delete these columns if already exist from previous function execution
      Path_FittedMod = NULL, Path_Coda = NULL)

  # # ..................................................................... ###

  # Merge posteriors and save as Hmsc model / coda object
  ecokit::cat_time("Merge posteriors and save as Hmsc model / coda object")

  Model_Info3 <- future.apply::future_lapply(
    X = seq_len(nrow(Model_Info2)),
    FUN = function(x) {

      if (Model_Info2$Post_Missing[[x]]) {
        return(
          list(
            Path_FittedMod = NA_character_,
            Path_Coda = NA_character_, Post_Aligned2 = NA))
      }

      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # Merge fitted models
      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      M_Name_Fit <- Model_Info2$M_Name_Fit[[x]]

      Path_FittedMod <- fs::path(
        Path_Fitted_Models, paste0(M_Name_Fit, "_Model.", out_extension))

      # Check if model exists and is valid
      ModelFileOkay <- ecokit::check_data(Path_FittedMod, warning = FALSE)

      if (isFALSE(ModelFileOkay)) {

        # delete corrupted file
        if (file.exists(Path_FittedMod)) {
          fs::file_delete(Path_FittedMod)
        }

        # Get posteriors
        Posts <- purrr::map(
          .x = as.character(Model_Info2$Post_Path[[x]]),
          .f = IASDT.R::mod_get_posteriors, from_JSON = from_JSON)

        # Convert to Hmsc object
        # Try with `alignPost = TRUE`
        Model_Fit <- Hmsc::importPosteriorFromHPC(
          m = ecokit::load_as(Model_Info2$M_Init_Path[x]),
          postList = Posts, nSamples = Model_Info2$M_samples[x],
          thin = Model_Info2$M_thin[x],
          transient = Model_Info2$M_transient[x],
          alignPost = TRUE) %>%
          try(silent = TRUE)

        # If failed, use `alignPost = FALSE`
        if (inherits(Model_Fit, "try-error")) {
          Model_Fit <- try(
            Hmsc::importPosteriorFromHPC(
              m = ecokit::load_as(Model_Info2$M_Init_Path[x]),
              postList = Posts, nSamples = Model_Info2$M_samples[x],
              thin = Model_Info2$M_thin[x],
              transient = Model_Info2$M_transient[x], alignPost = FALSE),
            silent = TRUE)
          Post_Aligned2 <- FALSE
        } else {
          Post_Aligned2 <- TRUE
        }

        rm(Posts, envir = environment())
        invisible(gc())

        if (inherits(Model_Fit, "try-error")) {
          return(
            list(
              Path_FittedMod = NA_character_,
              Path_Coda = NA_character_,
              Post_Aligned2 = NA))
        }
        ecokit::save_as(
          object = Model_Fit, object_name = paste0(M_Name_Fit, "_Model"),
          out_path = Path_FittedMod)

      } else {
        Post_Aligned2 <- Model_Info2$Post_Aligned[[x]]
      }

      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # Convert to Coda object
      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      Path_Coda <- fs::path(
        Path_Coda, paste0(M_Name_Fit, "_Coda.", out_extension))

      if (!exists("Model_Fit") && !file.exists(Path_FittedMod)) {
        return(
          list(
            Path_FittedMod = NA_character_, Path_Coda = NA_character_,
            Post_Aligned2 = NA))
      }

      if (!exists("Model_Fit")) {
        Model_Fit <- ecokit::load_as(Path_FittedMod)
      }

      # Check if coda file exists and is valid
      if (isFALSE(ecokit::check_data(Path_Coda, warning = FALSE))) {

        if (file.exists(Path_Coda)) {
          fs::file_delete(Path_Coda)
        }

        Mod_Coda <- Hmsc::convertToCodaObject(
          Model_Fit, spNamesNumbers = c(TRUE, FALSE),
          covNamesNumbers = c(TRUE, FALSE))

        ecokit::save_as(
          object = Mod_Coda, object_name = paste0(M_Name_Fit, "_Coda"),
          out_path = Path_Coda)

        rm(Mod_Coda, envir = environment())
      }

      rm(Model_Fit, envir = environment())
      invisible(gc())

      # Return list of objects
      return(
        list(
          Path_FittedMod = Path_FittedMod, Path_Coda = Path_Coda,
          Post_Aligned2 = Post_Aligned2))
    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = pkg_to_export,
    future.globals = c(
      "out_extension", "Model_Info2", "Path_Fitted_Models",
      "from_JSON", "Path_Coda"))

  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 1L)
    future::plan("sequential", gc = TRUE)
  }

  # # ..................................................................... ###

  ecokit::cat_time("Extract information on elapsed time and memory usage")

  Model_Info2 <- dplyr::mutate(Model_Info2, ModelPosts = Model_Info3) %>%
    tidyr::unnest_wider("ModelPosts") %>%
    dplyr::mutate(Post_Aligned = dplyr::coalesce(Post_Aligned2)) %>%
    dplyr::select(-Post_Aligned2) %>%
    dplyr::mutate(

      # Check if both merged fitted model and coda file exist
      Model_Finished = purrr::map2_lgl(
        .x = Path_FittedMod, .y = Path_Coda,
        .f = ~all(file.exists(c(.x, .y)))),

      # Extract fitting time from the progress file
      FittingTime = purrr::map(
        .x = Path_ModProg,
        .f = ~{
          purrr::map(
            .x = .x,
            .f = function(File) {
              if (file.exists(File)) {
                readr::read_lines(file = File, progress = FALSE) %>%
                  stringr::str_subset("Whole Gibbs sampler elapsed") %>%
                  stringr::str_remove("Whole Gibbs sampler elapsed") %>%
                  stringr::str_trim() %>%
                  as.numeric() %>%
                  magrittr::divide_by(60) %>%
                  round(1)
              } else {
                NA_real_
              }
            }) %>%
            unlist()
        }),

      # Mean elapsed time
      FittingTimeMean = purrr::map2_dbl(
        .x = Model_Finished, .y = FittingTime,
        .f = ~{
          if (.x) {
            mean(.y)
          } else {
            NA
          }
        }),

      # Maximum memory
      FittingMemory = purrr::map(
        .x = Path_ModProg,
        .f = ~{
          purrr::map(
            .x = .x,
            .f = function(File) {
              if (file.exists(File)) {
                readr::read_lines(file = File, progress = FALSE) %>%
                  stringr::str_subset("Maximum resident set size") %>%
                  stringr::str_remove_all(
                    "\t|:|Maximum resident set size \\(kbytes\\)") %>%
                  stringr::str_trim() %>%
                  as.numeric() %>%
                  magrittr::divide_by(1024 * 1024) %>%
                  round(2)
              } else {
                NA_real_
              }
            }) %>%
            unlist()
        }),


      # Mean memory
      FittingMemoryMean = purrr::map2_dbl(
        .x = Model_Finished, .y = FittingMemory,
        .f = ~{
          if (.x) {
            mean(.y)
          } else {
            NA
          }
        }))

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Print to the console the name of failed models and number of missing chain
  # files

  if (print_incomplete) {
    MissingModelVars <- Model_Info2 %>%
      dplyr::filter(!Model_Finished) %>%
      dplyr::mutate(
        NMissingChains = purrr::map_int(
          .x = Post_Path, .f = ~sum(!file.exists(.x))),
        MissingModels = paste0(M_Name_Fit, " (", NMissingChains, " chains)")
      ) %>%
      dplyr::pull(MissingModels) %>%
      gtools::mixedsort()

    if (length(MissingModelVars) > 0) {
      ecokit::cat_time("Unsuccessful models")
      purrr::walk(
        .x = MissingModelVars, .f = ecokit::cat_time,
        cat_timestamp = FALSE, level = 1L)
    }
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Save Model_Info to disk

  if (is.null(model_info_name)) {
    ecokit::save_as(
      object = Model_Info2, object_name = "Model_Info", out_path = Path_ModInfo)
  } else {
    ecokit::save_as(
      object = Model_Info2, object_name = model_info_name,
      out_path = fs::path(model_dir, paste0(model_info_name, ".RData")))
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_diff(init_time = .start_time, prefix = "Merging chains took ")

  return(invisible(NULL))
}

# # ========================================================================== #
# # ========================================================================== #

## |------------------------------------------------------------------------| #
# mod_merge_chains_CV ----
# Merge chains for cross-validated models
## |------------------------------------------------------------------------| #

#' @export
#' @rdname mod_merge_chains
#' @name mod_merge_chains
#' @order 2
#' @author Ahmed El-Gabbas

mod_merge_chains_CV <- function(
    model_dir = NULL, n_cores = 8L, strategy = "multisession",
    CV_names = c("CV_Dist", "CV_Large"), from_JSON = FALSE,
    out_extension = "qs2") {

  # # ..................................................................... ###

  ecokit::cat_time("Merge chains for cross-validated models")
  .start_time <- lubridate::now(tzone = "CET")

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Model_Finished <- FittingMemory <- Path_Post <- Path_ModFitted <-
    FittingTime <- Path_ModProg <- CV_name <- NULL

  # # ..................................................................... ###

  # Checking arguments ----

  if (is.null(model_dir)) {
    ecokit::stop_ctx(
      "`model_dir` cannot be empty", model_dir = model_dir,
      include_backtrace = TRUE)
  }

  if (is.null(n_cores)) {
    ecokit::stop_ctx(
      "`n_cores` cannot be empty", n_cores = n_cores, include_backtrace = TRUE)
  }

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  ecokit::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("model_dir", "out_extension", "strategy"))
  ecokit::check_args(
    args_all = AllArgs, args_to_check = "n_cores", args_type = "numeric")
  ecokit::check_args(
    args_all = AllArgs, args_type = "logical", args_to_check = "from_JSON")
  rm(AllArgs, envir = environment())

  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores <= 0) {
    ecokit::stop_ctx(
      "n_cores must be a single positive integer.", n_cores = n_cores,
      include_backtrace = TRUE)
  }

  if (strategy == "sequential") {
    n_cores <- 1L
  }
  if (length(strategy) != 1L) {
    ecokit::stop_ctx(
      "`strategy` must be a character vector of length 1",
      strategy = strategy, length_strategy = length(strategy))
  }
  valid_strategy <- c("sequential", "multisession", "multicore", "cluster")
  if (!strategy %in% valid_strategy) {
    ecokit::stop_ctx("Invalid `strategy` value", strategy = strategy)
  }

  if (!dir.exists(model_dir)) {
    ecokit::stop_ctx(
      "`model_dir` directory does not exist", model_dir = model_dir,
      include_backtrace = TRUE)
  }

  if (length(out_extension) != 1) {
    ecokit::stop_ctx(
      "`out_extension` must be a single string.",
      out_extension = out_extension, length_out_extension = out_extension,
      include_backtrace = TRUE)
  }

  if (!out_extension %in% c("qs2", "RData")) {
    ecokit::stop_ctx(
      "`out_extension` must be either 'qs2' or 'RData'.",
      out_extension = out_extension, include_backtrace = TRUE)
  }

  if (!all(CV_names %in% c("CV_Dist", "CV_Large", "CV_SAC"))) {
    ecokit::stop_ctx(
      paste0(
        "Invalid value for CV_names argument. Valid values ",
        "are: 'CV_Dist', 'CV_Large', or `CV_SAC`"),
      CV_names = CV_names, include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c(
      "Hmsc", "purrr", "ecokit", "IASDT.R", "dplyr", "stringr", "fs"),
    strategy = strategy)

  # # ..................................................................... ###

  # Creating paths -----

  Path_Fitted_Models <- fs::path(model_dir, "Model_Fitting_CV", "Model_Fitted")
  fs::dir_create(Path_Fitted_Models)

  # # ..................................................................... ###

  # Loading CV model info -----

  Path_CV_DT <- fs::path(model_dir, "Model_Fitting_CV", "CV_DT.RData")
  if (!file.exists(Path_CV_DT)) {
    ecokit::stop_ctx(
      "CV_DT file does not exist", Path_CV_DT = Path_CV_DT,
      include_backtrace = TRUE)
  }
  if (isFALSE(ecokit::check_data(Path_CV_DT, warning = FALSE))) {
    ecokit::stop_ctx(
      "CV_DT file is not a valid file", Path_CV_DT = Path_CV_DT,
      include_backtrace = TRUE)
  }

  CV_DT <- ecokit::load_as(Path_CV_DT) %>%
    # filter only selected cross-validation strategies
    dplyr::filter(CV_name %in% stringr::str_remove(CV_names, "CV_"))

  # # ..................................................................... ###

  # Prepare working in parallel
  if (n_cores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    ecokit::set_parallel(
      n_cores = min(n_cores, nrow(CV_DT)), level = 1L, future_max_size = 800L,
      strategy = strategy)
    withr::defer(future::plan("sequential", gc = TRUE))
  }

  # # ..................................................................... ###

  CV_DT <- CV_DT %>%
    dplyr::mutate(
      Post_Missing = furrr::future_map_lgl(
        .x = Path_Post,
        .f = function(x) {
          purrr::map_lgl(
            .x = as.character(x),
            .f = function(y) {

              if (isFALSE(fs::file_exists(y))) {
                return(TRUE)
              } else if (ecokit::check_rds(y)) {
                return(FALSE)
              } else {
                fs::file_delete(y)
                return(TRUE)
              }
            }) %>%
            any()
        },
        .options = furrr::furrr_options(seed = TRUE, packages = "fs")),

      Path_ModFitted = stringr::str_replace_all(
        Path_ModFitted, "RData$|qs2$", out_extension))

  invisible(gc())

  # # ..................................................................... ###

  # Merge posteriors and save as Hmsc object
  ecokit::cat_time("Merge posteriors and save as Hmsc object", level = 1L)

  CV_DT2 <- future.apply::future_lapply(
    X = seq_len(nrow(CV_DT)),
    FUN = function(x) {

      if (CV_DT$Post_Missing[[x]]) {
        return(NA)
      }

      Path_Fitted <- CV_DT$Path_ModFitted[[x]]

      # Check if model exists and is valid
      ModelFileOkay <- ecokit::check_data(Path_Fitted, warning = FALSE)

      if (isFALSE(ModelFileOkay)) {

        # delete corrupted file
        if (file.exists(Path_Fitted)) {
          fs::file_delete(Path_Fitted)
        }

        # Get posteriors
        Posts <- purrr::map(
          .x = as.character(CV_DT$Path_Post[[x]]),
          .f = IASDT.R::mod_get_posteriors, from_JSON = from_JSON)

        # Convert to Hmsc object
        Model_init_rds <- ecokit::load_as(CV_DT$Path_ModInit_rds[x])
        Model_init <- ecokit::load_as(CV_DT$Path_ModInit[x])
        Model_Fit <- Hmsc::importPosteriorFromHPC(
          m = Model_init, postList = Posts, nSamples = Model_init_rds$samples,
          thin = Model_init_rds$thin, transient = Model_init_rds$transient,
          alignPost = TRUE) %>%
          try(silent = TRUE)

        # If failed, use `alignPost = FALSE`
        if (inherits(Model_Fit, "try-error")) {
          Model_Fit <- try(
            Hmsc::importPosteriorFromHPC(
              m = ecokit::load_as(CV_DT$M_Init_Path[x]),
              postList = Posts, nSamples = CV_DT$M_samples[x],
              thin = CV_DT$M_thin[x],
              transient = CV_DT$M_transient[x], alignPost = FALSE),
            silent = TRUE)
          Post_Aligned <- FALSE
        } else {
          Post_Aligned <- TRUE
        }

        rm(Posts, envir = environment())

        if (inherits(Model_Fit, "try-error")) {
          return(NA)
        }

        ecokit::save_as(
          object = Model_Fit,
          object_name = stringr::str_remove(
            basename(Path_Fitted), ".RData$|.qs2"),
          out_path = Path_Fitted)

      } else {
        Post_Aligned <- NA
      }

      invisible(gc())

      # Return list of objects
      return(Post_Aligned)
    },
    future.scheduling = Inf, future.seed = TRUE,
    future.packages = pkg_to_export,
    future.globals = c("out_extension", "CV_DT", "from_JSON"))

  invisible(gc())

  # # ..................................................................... ###

  # Check saved Hmsc object and extract info on model fitting
  ecokit::cat_time(
    "Check saved Hmsc object and extract info on model fitting", level = 1L)

  # Check if both merged fitted model file exist
  CV_DT <- CV_DT %>%
    dplyr::mutate(

      Post_Aligned = unlist(CV_DT2),

      Model_Finished = furrr::future_map_lgl(
        .x = Path_ModFitted, .f = ecokit::check_data, warning = FALSE,
        .options = furrr::furrr_options(seed = TRUE, packages = pkg_to_export)),

      # Extract fitting time from the progress file
      FittingTime = purrr::map(
        .x = Path_ModProg,
        .f = ~{
          purrr::map(
            .x = .x,
            .f = function(File) {
              if (file.exists(File)) {
                readr::read_lines(file = File, progress = FALSE) %>%
                  stringr::str_subset("Whole Gibbs sampler elapsed") %>%
                  stringr::str_remove("Whole Gibbs sampler elapsed") %>%
                  stringr::str_trim() %>%
                  as.numeric() %>%
                  magrittr::divide_by(60) %>%
                  round(1)
              } else {
                NA_real_
              }
            }) %>%
            unlist()
        }),

      # Mean elapsed time
      FittingTimeMean = purrr::map2_dbl(
        .x = Model_Finished, .y = FittingTime,
        .f = ~{
          if (.x) {
            mean(.y)
          } else {
            NA
          }
        }),

      # Maximum memory
      FittingMemory = purrr::map(
        .x = Path_ModProg,
        .f = ~{
          purrr::map(
            .x = .x,
            .f = function(File) {
              if (file.exists(File)) {
                readr::read_lines(file = File, progress = FALSE) %>%
                  stringr::str_subset("Maximum resident set size") %>%
                  stringr::str_remove_all(
                    "\t|:|Maximum resident set size \\(kbytes\\)") %>%
                  stringr::str_trim() %>%
                  as.numeric() %>%
                  magrittr::divide_by(1024 * 1024) %>%
                  round(2)
              } else {
                NA_real_
              }
            }) %>%
            unlist()
        }),

      # Mean memory
      FittingMemoryMean = purrr::map2_dbl(
        .x = Model_Finished, .y = FittingMemory,
        .f = ~{
          if (.x) {
            mean(.y)
          } else {
            NA
          }
        }))

  invisible(gc())

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # stopping the cluster
  if (n_cores > 1) {
    ecokit::set_parallel(stop_cluster = TRUE, level = 2L)
    future::plan("sequential", gc = TRUE)
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Save Model_Info to disk

  ecokit::save_as(
    object = CV_DT, object_name = "CV_DT_fitted",
    out_path = fs::path(model_dir, "Model_Fitting_CV", "CV_DT_fitted.RData"))

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  ecokit::cat_diff(init_time = .start_time, prefix = "Merging chains took ")

  return(invisible(NULL))
}
