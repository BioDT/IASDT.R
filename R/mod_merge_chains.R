## |------------------------------------------------------------------------| #
# mod_merge_chains ----
## |------------------------------------------------------------------------| #

#' Merge model chains
#'
#' This function merges model chains, checks for missing or incomplete model
#' runs, and optionally prints information about incomplete models. It saves the
#' processed model objects to disk.
#' @param model_dir Character. Path to the root directory of the fitted model.
#'   Two folders will be created `Model_Fitted` and `Model_Coda` to store merged
#'   model and coda objects, respectively.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#' @param model_info_name Character. Name of the file (without extension) where
#'   the processed model information will be saved. If `NULL`, it overwrites the
#'   `Model_Info.RData` file in the `model_dir` directory. If provided, a new
#'   `.RData` file will be created with this name.
#' @param print_incomplete Logical. Whether to print the names of models that
#'   were not successfully fitted to the console. Defaults to `TRUE`.
#' @param from_JSON Logical. Whether to convert loaded models from JSON format
#'   before reading. Defaults to `FALSE`.
#' @param out_extension Character. File extension (without dot) for the files
#'   containing the fitted model and coda objects. It can be either `qs2` or
#'   `RData`. Defaults to `qs2`.
#' @name mod_merge_chains
#' @author Ahmed El-Gabbas
#' @return The function does not return anything but saves the processed model
#'   information to disk.
#' @export

mod_merge_chains <- function(
    model_dir = NULL, n_cores = NULL, model_info_name = NULL,
    print_incomplete = TRUE, from_JSON = FALSE, out_extension = "qs2") {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Post_Path <- Post_Path <- M_Name_Fit <- Path_FittedMod <- Path_Coda <-
    NMissingChains <- MissingModels <- Model_Finished <- Path_ModProg <-
    Post_Aligned2 <- FittingTime <- FittingMemory <- NULL

  # # ..................................................................... ###

  # Checking arguments ----

  if (is.null(model_dir) || is.null(n_cores)) {
    stop("`model_dir` and `n_cores` cannot be empty", call. = FALSE)
  }

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::check_args(
    args_all = AllArgs, args_type = "character",
    args_to_check = c("model_dir", "out_extension"))

  IASDT.R::check_args(
    args_all = AllArgs, args_to_check = "n_cores", args_type = "numeric")

  IASDT.R::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c("print_incomplete", "PlotVP"))

  rm(AllArgs, envir = environment())

  if (!dir.exists(model_dir)) {
    stop(
      "`model_dir` directory (`", model_dir, "`) does not exist", call. = FALSE)
  }

  if (length(out_extension) != 1) {
    stop("`out_extension` must be a single string.", call. = FALSE)
  }

  if (!out_extension %in% c("qs2", "RData")) {
    stop("`out_extension` must be either 'qs2' or 'RData'.", call. = FALSE)
  }

  # # ..................................................................... ###

  # Creating paths -----

  Path_Fitted_Models <- IASDT.R::path(model_dir, "Model_Fitted")
  Path_Coda <- IASDT.R::path(model_dir, "Model_Coda")
  fs::dir_create(c(Path_Fitted_Models, Path_Coda))

  # # ..................................................................... ###

  # Loading model info ----

  Path_ModInfo <- IASDT.R::path(model_dir, "Model_Info.RData")

  if (!file.exists(Path_ModInfo)) {
    stop("ModInfo file `", Path_ModInfo, "` does not exist", call. = FALSE)
  }

  # # ..................................................................... ###

  # Remove temp files and incomplete RDs files ----

  Path_Model_Fit <- IASDT.R::path(model_dir, "Model_Fitting_HPC")
  tempFiles <- list.files(
    path = Path_Model_Fit, pattern = ".rds_temp$", full.names = TRUE)

  if (length(tempFiles) > 0) {
    IASDT.R::cat_time(
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

  Model_Info2 <- IASDT.R::load_as(Path_ModInfo)

  # Prepare working on parallel -----
  if (n_cores == 1) {
    future::plan("future::sequential", gc = TRUE)
  } else {
    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE,
      future.seed = TRUE)
    c1 <- snow::makeSOCKcluster(min(n_cores, nrow(Model_Info2)))
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("future::cluster", workers = c1, gc = TRUE)
    withr::defer(future::plan("future::sequential", gc = TRUE))
  }


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
              } else if (IASDT.R::check_rds(y)) {
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

      Path_FittedMod <- IASDT.R::path(
        Path_Fitted_Models, paste0(M_Name_Fit, "_Model.", out_extension))

      # Check if model exists and is valid
      ModelFileOkay <- IASDT.R::check_data(Path_FittedMod, warning = FALSE)

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
          m = IASDT.R::load_as(Model_Info2$M_Init_Path[x]),
          postList = Posts, nSamples = Model_Info2$M_samples[x],
          thin = Model_Info2$M_thin[x],
          transient = Model_Info2$M_transient[x],
          alignPost = TRUE) %>%
          try(silent = TRUE)

        # If failed, use `alignPost = FALSE`
        if (inherits(Model_Fit, "try-error")) {
          Model_Fit <- try(
            Hmsc::importPosteriorFromHPC(
              m = IASDT.R::load_as(Model_Info2$M_Init_Path[x]),
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
        IASDT.R::save_as(
          object = Model_Fit, object_name = paste0(M_Name_Fit, "_Model"),
          out_path = Path_FittedMod)

      } else {
        Post_Aligned2 <- Model_Info2$Post_Aligned[[x]]
      }

      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # Convert to Coda object
      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      Path_Coda <- IASDT.R::path(
        Path_Coda, paste0(M_Name_Fit, "_Coda.", out_extension))

      if (!exists("Model_Fit") && !file.exists(Path_FittedMod)) {
        return(
          list(
            Path_FittedMod = NA_character_, Path_Coda = NA_character_,
            Post_Aligned2 = NA))
      }

      if (!exists("Model_Fit")) {
        Model_Fit <- IASDT.R::load_as(Path_FittedMod)
      }

      # Check if coda file exists and is valid
      if (isFALSE(IASDT.R::check_data(Path_Coda, warning = FALSE))) {

        if (file.exists(Path_Coda)) {
          fs::file_delete(Path_Coda)
        }

        Mod_Coda <- Hmsc::convertToCodaObject(
          Model_Fit, spNamesNumbers = c(TRUE, FALSE),
          covNamesNumbers = c(TRUE, FALSE))

        IASDT.R::save_as(
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
    future.globals = c(
      "out_extension", "Model_Info2", "Path_Fitted_Models",
      "from_JSON", "Path_Coda"),
    future.packages = c("Hmsc", "coda", "purrr", "IASDT.R", "dplyr"))

  if (n_cores > 1) {
    snow::stopCluster(c1)
    future::plan("future::sequential", gc = TRUE)
  }

  # # ..................................................................... ###

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
                readr::read_lines(file = File) %>%
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
        .x =  Model_Finished, .y = FittingTime,
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
                readr::read_lines(file = File) %>%
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
        .x =  Model_Finished, .y = FittingMemory,
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
      IASDT.R::sort2()

    if (length(MissingModelVars) > 0) {
      IASDT.R::cat_time("Unsuccessful models")
      purrr::walk(MissingModelVars, IASDT.R::cat_time, level = 1)
    }
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Save Model_Info to disk

  if (is.null(model_info_name)) {
    IASDT.R::save_as(
      object = Model_Info2, object_name = "Model_Info", out_path = Path_ModInfo)
  } else {
    IASDT.R::save_as(
      object = Model_Info2, object_name = model_info_name,
      out_path = IASDT.R::path(model_dir, paste0(model_info_name, ".RData")))
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  IASDT.R::cat_diff(init_time = .StartTime, prefix = "Merging chains took ")

  return(invisible(NULL))
}
