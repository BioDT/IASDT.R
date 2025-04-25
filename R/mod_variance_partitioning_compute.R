## |------------------------------------------------------------------------| #
# variance_partitioning_compute ----

#' Computes and visualize variance partitioning of Hmsc models
#'
#' The **`variance_partitioning_compute()`** function computes variance
#' components with respect to given grouping of fixed effects and levels of
#' random effects. This function inherits the main functionality from the
#' `Hmsc::computeVariancePartitioning` function, but with the added
#' functionality of parallel computation and using TensorFlow.<br/>The
#' **`variance_partitioning_plot()`** function generates plots for variance
#' partitioning as JPEG files. It allows for sorting the predictors and species;
#' e.g., by the mean value per predictor; and by original species order. It also
#' plots the raw variance partitioning (relative variance partitioning
#' multiplied by the Tjur-R<sup>2</sup> value).
#' @param path_model Character. Path to fitted `Hmsc` model object.
#' @param group vector of numeric values corresponding to group identifiers in
#'   groupnames. If the model was defined with `XData` and `XFormula`, the
#'   default is to use model terms.
#' @param group_names vector of names for each group of fixed effect. Should
#'   match `group`. If the model was defined with `XData` and `XFormula`, the
#'   default is to use the labels of model terms.
#' @param start index of first MCMC sample included
#' @param na.ignore Logical. If `TRUE`, covariates are ignored for sites where
#'   the focal species is NA when computing variance-covariance matrices for
#'   each species.
#' @param n_cores Integer. Number of CPU cores to use for computing variance
#'   partitioning using TensorFlow. This is only effective when `use_TF` is
#'   `TRUE`. Default: `1`.
#' @param chunk_size Integer. Size of each chunk of samples to process in
#'   parallel. Only relevant for TensorFlow. Default: `50`.
#' @param verbose Logical. Whether to print progress messages. Default: `TRUE`.
#' @param temp_cleanup Logical. Whether to delete temporary files after
#'   processing. Default: `TRUE`.
#' @param VP_commands_only Logical. If `TRUE`, returns the commands to run the
#'   Python script. Default is `FALSE`. Only relevant when `use_TF` is `TRUE`.
#' @param VP_file Character. Name of the output file to save the results.
#'   Default: `VarPar`.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param width,height Numeric. Width and height of the output plot in
#'   centimeters. Default: `30` and `15`, respectively.
#' @param Axis_text Numeric. Size of the axis text. Default: `4`.
#' @author Ahmed El-Gabbas
#' @inheritParams predict_hmsc
#' @inheritParams predict_latent_factor
#' @name variance_partitioning
#' @rdname variance_partitioning
#' @order 1
#' @export

variance_partitioning_compute <- function(
    path_model, group = NULL, group_names = NULL, start = 1L, na.ignore = FALSE,
    n_cores = 8L, use_TF = TRUE, TF_environ = NULL, TF_use_single = FALSE,
    temp_cleanup = TRUE, chunk_size = 50L, verbose = TRUE,
    VP_file = "VarPar", VP_commands_only = FALSE) {

  # # .................................................................... ###

  if (isFALSE(verbose)) {
    sink(file = nullfile())
    on.exit(sink(), add = TRUE)
  }

  .StartTime <- lubridate::now(tzone = "CET")

  # # .................................................................... ###

  # Check if the virtual environment and Python scripts exist

  if (use_TF) {

    # Check python virtual environment
    if (isFALSE(VP_commands_only) && .Platform$OS.type == "windows" &&
        (is.null(TF_environ) || !dir.exists(TF_environ))) {
      IASDT.R::stop_ctx(
        paste0(
          "When running on Windows and `use_TF` is TRUE, `TF_environ` must ",
          "be specified and point to an existing directory with a ",
          "Python virtual environment"),
        VP_commands_only = VP_commands_only, OS = .Platform$OS.type,
        TF_environ = TF_environ)
    }

    # Determine the Python executable path

    # On Windows, the TF calculations has to be done through a valid virtual
    # environment; the path to the virtual environment must be specified in
    # `TF_environ`. On LUMI, this is not needed as the compatible python
    # installation is loaded automatically when loading tensorflow module. When
    # using another HPC system, the function needs to be adapted accordingly to
    # point to a valid python virtual environment.

    if (.Platform$OS.type == "windows") {
      python_executable <- IASDT.R::path(TF_environ, "Scripts", "python.exe")

      if (isFALSE(VP_commands_only) && !file.exists(python_executable)) {
        IASDT.R::stop_ctx(
          "Python executable not found in the virtual environment",
          python_executable = python_executable)
      }
    } else {
      python_executable <- "/usr/bin/time -v python3"
    }

    # Check GPU availability
    if (isFALSE(VP_commands_only) && .Platform$OS.type == "windows") {
      result <- system(
        paste0(
          python_executable,
          " -c \"import tensorflow as tf; print(len(",
          "tf.config.list_physical_devices('GPU')))\""),
        intern = TRUE)
      N_GPU <- result[length(result)]
      if (N_GPU == 0) {
        IASDT.R::cat_time(
          "No GPU found; Calculations will use CPU.",
          cat_timestamp = FALSE, cat_bold = TRUE, cat_red = TRUE)
      } else {
        IASDT.R::cat_time(
          paste0(N_GPU, " GPUs were found. Calculations will use GPU."),
          cat_timestamp = FALSE, cat_bold = TRUE, cat_red = TRUE)
      }
    }

    # Paths to the Python scripts
    Script_geta <- system.file("VP_geta.py", package = "IASDT.R")
    Script_getf <- system.file("VP_getf.py", package = "IASDT.R")
    Script_gemu <- system.file("VP_gemu.py", package = "IASDT.R")
    if (!all(file.exists(c(Script_geta, Script_getf, Script_gemu)))) {
      IASDT.R::stop_ctx(
        "Necessary python scripts do not exist",
        Script_geta = Script_geta, Script_getf = Script_getf,
        Script_gemu = Script_gemu)
    }
  }

  # # .................................................................... ###

  # Load model object ------
  IASDT.R::cat_time("Load model object")

  if (is.null(path_model) || !file.exists(path_model)) {
    IASDT.R::stop_ctx(
      "Model path is NULL or does not exist", path_model = path_model)
  }

  Model <- IASDT.R::load_as(path_model)

  # # .................................................................... ###

  # Create folder for variance partitioning results
  Path_VarPar <- IASDT.R::path(
    dirname(dirname(path_model)),
    "Model_Postprocessing", "Variance_Partitioning")
  fs::dir_create(Path_VarPar)

  # # .................................................................... ###

  ny <- Model$ny
  nc <- Model$nc
  ns <- Model$ns
  nr <- Model$nr

  if (is.null(group)) {
    # names of variables used in the model
    group_names <- names(Model$XData)

    # actual variables used in the model, including quadratic terms
    ModelVars <- dimnames(Model$X)[[2]][-1]

    # group variable to combine variable and its quadratic terms together
    group <- purrr::map(
      ModelVars, ~ which(stringr::str_detect(.x, group_names))) %>%
      unlist() %>%
      # add intercept to the first group
      c(.[1], .)
  }

  # If na.ignore=T, convert XData to a list
  if (na.ignore) {
    xl <- list()
    for (s in seq_len(ns)) {
      xl[[s]] <- Model$X
    }
    Model$X <- xl
  }

  switch(
    class(Model$X)[1L],
    matrix = {
      cMA <- stats::cov(Model$X)
    },
    list = {
      if (na.ignore) {
        cMA <- list()
        for (s in seq_len(ns)) {
          cMA[[s]] <- stats::cov(Model$X[[s]][which(Model$Y[, s] > -Inf), ])
        }
      } else {
        cMA <- lapply(Model$X, stats::cov)
      }
    }
  )

  # # .................................................................... ###

  # Prepare postList-----

  IASDT.R::cat_time("Prepare postList")
  postList <- Hmsc::poolMcmcChains(Model$postList, start = start) %>%
    purrr::map(
      ~ {
        Items2Delete <- c(
          "Eta", "Psi", "V", "sigma", "Delta", "Alpha",
          "rho", "wRRR", "PsiRRR", "DeltaRRR")
        .x[Items2Delete] <- NULL
        return(.x)
      })

  # Remove unnecessary elements from the model object
  names_to_remove <- c(
    "postList", "Y", "XScaled", "rL", "ranLevels", "XData", "dfPi",
    "studyDesign", "C", "Pi", "phyloTree", "XFormula", "XScalePar",
    "aSigma", "bSigma", "TrScaled", "YScalePar", "call", "rhopw",
    "distr", "V0", "UGamma", "YScaled")
  Model[names_to_remove] <- NULL

  poolN <- length(postList)
  ngroups <- max(group)

  invisible(gc())

  # # .................................................................... ###
  # # .................................................................... ###

  # Prepare `la`/`lf`/`lmu` lists -----

  if (use_TF) {

    # Create the temporary directory
    Path_Temp <- IASDT.R::path(dirname(dirname(path_model)), "TEMP_VP")
    fs::dir_create(Path_Temp)

    FileSuffix <- stringr::str_pad(
      string = seq_len(poolN), pad = "0", width = 4)

    # List of feather files resulted from `geta` function
    Files_la <- IASDT.R::path(
      Path_Temp, paste0("VP_A_", FileSuffix, ".feather"))
    Files_la_Exist <- all(file.exists(Files_la))

    # List of feather files resulted from `getf` function
    Files_lf <- IASDT.R::path(
      Path_Temp, paste0("VP_F_", FileSuffix, ".feather"))
    Files_lf_Exist <- all(file.exists(Files_lf))

    # List of feather files resulted from `gemu` function
    Files_lmu <- IASDT.R::path(
      Path_Temp, paste0("VP_Mu_", FileSuffix, ".feather"))
    Files_lmu_Exist <- all(file.exists(Files_lmu))

    Beta_Files <- IASDT.R::path(
      Path_Temp, paste0("VP_Beta_", FileSuffix, ".feather"))

    if (all(c(Files_la_Exist, Files_lf_Exist, Files_lmu_Exist))) {

      # All feather data are already processed and available on disk
      IASDT.R::cat_time(
        "Data for `la`/`lf`/`lmu` lists were already processed")

    } else {

      ## Prepare la/lf/lmu lists using TensorFlow ----
      IASDT.R::cat_time("Prepare la/lf/lmu lists using TensorFlow")

      ### Prepare data for TensorFlow ----
      IASDT.R::cat_time("Prepare data for TensorFlow", level = 1)

      #### X data -----
      # needed only to calculate `geta` and `getf` functions
      if (!all(c(Files_la_Exist, Files_lf_Exist))) {
        IASDT.R::cat_time("X", level = 2)
        Path_X <- IASDT.R::path(Path_Temp, "VP_X.feather")
        if (!file.exists(Path_X)) {
          arrow::write_feather(as.data.frame(Model$X), Path_X)
        }
      }

      #### Tr / Gamma ------
      # needed only to calculate `geta` and `gemu` functions
      if (!all(c(Files_la_Exist, Files_lmu_Exist))) {
        IASDT.R::cat_time("Tr", level = 2)
        Path_Tr <- IASDT.R::path(Path_Temp, "VP_Tr.feather")
        if (!file.exists(Path_Tr)) {
          arrow::write_feather(as.data.frame(Model$Tr), Path_Tr)
        }

        # Gamma - convert each list item into a column in a data frame
        IASDT.R::cat_time("Gamma", level = 2)
        Path_Gamma <- IASDT.R::path(Path_Temp, "VP_Gamma.feather")
        if (!file.exists(Path_Gamma)) {
          Gamma_data <- postList %>%
            purrr::map(~as.vector(.x[["Gamma"]])) %>%
            as.data.frame() %>%
            stats::setNames(paste0("Sample_", seq_len(ncol(.))))
          arrow::write_feather(Gamma_data, Path_Gamma)
        }
      }

      #### Beta -----
      # only needed for `getf` function
      if (!Files_lf_Exist) {
        # Beta -- Each element of Beta is a matrix, so each list item is saved
        # to separate feather file
        IASDT.R::cat_time("Beta", level = 2)

        if (!all(file.exists(Beta_Files))) {

          IASDT.R::cat_time("Prepare working in parallel", level = 3)
          c1 <- parallel::makePSOCKcluster(n_cores)
          on.exit(try(parallel::stopCluster(c1), silent = TRUE), add = TRUE)

          IASDT.R::cat_time("Export necessary objects to cores", level = 3)
          parallel::clusterExport(
            cl = c1, varlist = c("Beta_Files", "postList"),
            envir = environment())

          IASDT.R::cat_time("Load libraries on each core", level = 3)
          invisible(
            parallel::clusterEvalQ(
              cl = c1, expr = sapply("arrow", library, character.only = TRUE)))

          IASDT.R::cat_time("Processing beta in parallel", level = 3)
          Beta0 <- parallel::parLapply(
            cl = c1,
            X = seq_along(postList),
            fun = function(x) {
              Beta_File <- Beta_Files[x]
              if (!file.exists(Beta_File)) {
                Beta <- as.data.frame(postList[[x]][["Beta"]])
                arrow::write_feather(x = Beta, sink = Beta_File)
              }
              return(NULL)
            })
          rm(Beta0)

          IASDT.R::cat_time("stop cluster", level = 3)
          parallel::stopCluster(c1)
        }
      }
      invisible(gc())

      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||

      ### Processing geta -----

      if (Files_la_Exist) {

        IASDT.R::cat_time(
          "All `la` data were already available on disk", level = 1)

      } else {

        IASDT.R::cat_time("Processing `geta` function", level = 1)
        Path_Out_a <- IASDT.R::path(Path_Temp, "VP_A.feather") %>%
          IASDT.R::normalize_path()

        cmd_a <- paste(
          python_executable, Script_geta,
          "--tr", IASDT.R::normalize_path(Path_Tr, must_work = TRUE),
          "--x", IASDT.R::normalize_path(Path_X, must_work = TRUE),
          "--gamma", IASDT.R::normalize_path(Path_Gamma, must_work = TRUE),
          "--output", Path_Out_a,
          "--ncores", n_cores,
          "--chunk_size", chunk_size)

        if (TF_use_single) {
          cmd_a <- paste0(cmd_a, " --use_single")
        }

        if (VP_commands_only) {

          # Save command to file
          # Redirect results of time to a log file
          path_log_a <- stringr::str_replace(Path_Out_a, ".feather", ".log")
          cmd_a <- paste0(cmd_a, paste0(" >> ", path_log_a, " 2>&1"))
          readr::write_lines(
            x = cmd_a,
            file = IASDT.R::path(Path_Temp, "VP_A_Command.txt"), append = FALSE)

        } else {

          # Run the command using system
          la <- system(cmd_a, wait = TRUE, intern  = TRUE)

          # Check for errors
          if (inherits(la, "error") || la[length(la)] != "Done") {
            IASDT.R::stop_ctx(
              "Error in computing geta", la = la, class_la = class(la))
          }

          if (length(la) != 1) {
            cat(la, sep = "\n")
          }
        }
      }

      invisible(gc())

      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||

      ### Processing getf ----

      if (Files_la_Exist) {

        IASDT.R::cat_time(
          "All `lf` data were already available on disk", level = 1)

      } else {

        IASDT.R::cat_time("Processing `getf` function", level = 1)
        Path_Out_f <- IASDT.R::path(Path_Temp, "VP_F.feather") %>%
          IASDT.R::normalize_path()

        cmd_f <- paste(
          python_executable, Script_getf,
          "--x", IASDT.R::normalize_path(Path_X, must_work = TRUE),
          "--beta_dir", IASDT.R::normalize_path(Path_Temp, must_work = TRUE),
          "--output", Path_Out_f,
          "--ncores", n_cores)

        if (TF_use_single) {
          cmd_f <- paste0(cmd_f, " --use_single")
        }

        if (VP_commands_only) {

          # Save command to file
          # Redirect results of time to a log file
          path_log_f <- stringr::str_replace(Path_Out_f, ".feather", ".log")
          cmd_f <- paste0(cmd_f, paste0(" >> ", path_log_f, " 2>&1"))
          readr::write_lines(
            x = cmd_f,
            file = IASDT.R::path(Path_Temp, "VP_F_Command.txt"),
            append = FALSE)

        } else {

          # Run the command using system
          lf <- system(cmd_f, wait = TRUE, intern  = TRUE)

          # Check for errors
          if (inherits(lf, "error") || lf[length(lf)] != "Done") {
            IASDT.R::stop_ctx(
              "Error in computing geta", lf = lf, class_lf = class(lf))
          }

          if (length(lf) != 1) {
            cat(lf, sep = "\n")
          }
        }
      }

      invisible(gc())

      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||

      ### Processing gemu ----

      if (Files_lmu_Exist) {

        IASDT.R::cat_time(
          "All `lmu` data were already available on disk", level = 1)

      } else {

        IASDT.R::cat_time("Processing `gemu` function", level = 1)
        Path_Out_mu <- IASDT.R::path(Path_Temp, "VP_Mu.feather") %>%
          IASDT.R::normalize_path()

        cmd_mu <- paste(
          python_executable, Script_gemu,
          "--tr", IASDT.R::normalize_path(Path_Tr, must_work = TRUE),
          "--gamma", IASDT.R::normalize_path(Path_Gamma, must_work = TRUE),
          "--output", Path_Out_mu,
          "--ncores", n_cores,
          "--chunk_size", chunk_size)

        if (TF_use_single) {
          cmd_mu <- paste0(cmd_mu, " --use_single")
        }

        if (VP_commands_only) {

          # Save command to file
          # Redirect results of time to a log file
          path_log_mu <- stringr::str_replace(Path_Out_mu, ".feather", ".log")
          cmd_mu <- paste0(cmd_mu, paste0(" >> ", path_log_mu, " 2>&1"))
          readr::write_lines(
            x = cmd_mu,
            file = IASDT.R::path(Path_Temp, "VP_mu_Command.txt"),
            append = FALSE)

        } else {

          # Run the command using system
          lmu <- system(cmd_mu, wait = TRUE, intern  = TRUE)

          # Check for errors
          if (inherits(lmu, "error") || lmu[length(lmu)] != "Done") {
            IASDT.R::stop_ctx(
              "Error in computing geta", lmu = lmu, class_lmu = class(lmu))
          }

          if (length(lmu) != 1) {
            cat(lmu, sep = "\n")
          }
        }
      }
    }

    invisible(gc())

    if (VP_commands_only) {
      return(invisible(NULL))
    }

  } else {

    # Prepare la/lf/lmu lists using original R code -----

    IASDT.R::cat_time("Prepare la/lf/lmu lists using original R code")

    ## geta -----
    IASDT.R::cat_time("Running geta", level = 1)
    geta <- function(a) {
      switch(
        class(Model$X)[1L],
        matrix = {
          res <- Model$X %*% (t(Model$Tr %*% t(a$Gamma)))
        },
        list = {
          res <- matrix(NA, ny, ns)
          for (j in seq_len(ns)) {
            res[, j] <- Model$X[[j]] %*% (t(Model$Tr[j, ] %*% t(a$Gamma)))
          }
        })
      return(res)
    }

    la <- lapply(postList, geta)

    invisible(gc())

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||

    ## getf ------

    IASDT.R::cat_time("Running getf", level = 1)
    getf <- function(a) {
      switch(
        class(Model$X)[1L],
        matrix = {
          res <- Model$X %*% (a$Beta)
        },
        list = {
          res <- matrix(NA, ny, ns)
          for (j in seq_len(ns)) res[, j] <- Model$X[[j]] %*% a$Beta[, j]
        })
      return(res)
    }

    lf <- lapply(postList, getf)

    invisible(gc())

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||

    ## gemu ------

    IASDT.R::cat_time("Running gemu", level = 1)
    gemu <- function(a) {
      res <- t(Model$Tr %*% t(a$Gamma))
      return(res)
    }

    lmu <- lapply(postList, gemu)

    invisible(gc())

  }

  # # .................................................................... ###

  # Running gebeta -------
  IASDT.R::cat_time("Running gebeta")
  gebeta <- function(a) {
    res <- a$Beta
    return(res)
  }
  lbeta <- lapply(postList, gebeta)

  # # .................................................................... ###

  # Remove Gamma from postList ------
  IASDT.R::cat_time("Remove Gamma from postList")
  postList <- purrr::map(
    postList,
    ~ {
      .x["Gamma"] <- NULL
      .x
    })

  invisible(gc())

  # # .................................................................... ###
  # # .................................................................... ###

  # Computing variance partitioning -------

  if (use_TF) {

    IASDT.R::cat_time("Computing variance partitioning in parallel")

    IASDT.R::cat_time("Split `lbeta` list into small qs2 files", level = 1)
    path_lbeta <- IASDT.R::path(Path_Temp, "lbeta")
    purrr::walk(
      .x = seq_along(lbeta),
      .f = ~ {
        IASDT.R::save_as(
          object = lbeta[[.x]],
          out_path = IASDT.R::path(
            path_lbeta,
            paste0(
              "lbeta_", stringr::str_pad(.x, width = 4, pad = "0"), ".qs2")))
      })

    IASDT.R::cat_time("Split `postList` list into small qs2 files", level = 1)
    path_postList <- IASDT.R::path(Path_Temp, "postList")
    purrr::walk(
      .x = seq_along(postList),
      .f = ~{
        IASDT.R::save_as(
          object = postList[[.x]],
          out_path = IASDT.R::path(
            path_postList,
            paste0(
              "post_", stringr::str_pad(.x, width = 4, pad = "0"), ".qs2")))
      })

    IASDT.R::cat_time("removing `postList` and `lbeta` list objects", level = 1)
    n_postList <- length(postList)
    rm(postList, lbeta, envir = environment())
    invisible(gc())

    IASDT.R::cat_time("Prepare working in parallel", level = 1)
    IASDT.R::set_parallel(n_cores = n_cores, level = 2)
    withr::defer(future::plan("future::sequential", gc = TRUE))

    IASDT.R::cat_time("Processing in parallel", level = 1)
    Res <- future.apply::future_lapply(
      X = seq_len(n_postList),
      FUN = function(i) {

        mm <- methods::getMethod("%*%", "Matrix")

        curr_postList <- IASDT.R::path(
          path_postList,
          paste0(
            "post_", stringr::str_pad(i, width = 4, pad = "0"), ".qs2")) %>%
          IASDT.R::load_as()
        Beta <- curr_postList$Beta
        Lambdas <- curr_postList$Lambda

        # Suppress warnings when no trait information is used in the models
        # cor(Beta[k, ], lmu[k, ]) : the standard deviation is zero
        DT_lmu <- as.matrix(arrow::read_feather(Files_lmu[i]))
        curr_lbeta <- IASDT.R::path(
          path_lbeta,
          paste0(
            "lbeta_", stringr::str_pad(i, width = 4, pad = "0"), ".qs2")) %>%
          IASDT.R::load_as()
        R2T.Beta <- purrr::map_dbl(
          .x = seq_len(nc),
          .f = ~ {
            suppressWarnings(stats::cor(curr_lbeta[.x, ], DT_lmu[.x, ])^2)
          })

        fixed1 <- matrix(0, nrow = ns, ncol = 1)
        fixedsplit1 <- matrix(0, nrow = ns, ncol = ngroups)
        random1 <- matrix(0, nrow = ns, ncol = nr)

        DT_la <- as.matrix(arrow::read_feather(Files_la[i]))
        DT_lf <- as.matrix(arrow::read_feather(Files_lf[i]))
        # a <- DT_la - matrix(rep(rowMeans(DT_la), ns), ncol = ns)
        # f <- DT_lf - matrix(rep(rowMeans(DT_lf), ns), ncol = ns)
        a <- DT_la - rowMeans(DT_la)
        f <- DT_lf - rowMeans(DT_lf)

        res1 <- sum((rowSums((a * f)) / (ns - 1))^2)
        res2 <- sum((rowSums((a * a)) / (ns - 1)) *
                      (rowSums((f * f)) / (ns - 1)))
        R2T.Y <- res1 / res2


        for (j in seq_len(ns)) {
          switch(
            class(Model$X)[1L],
            matrix = {
              cM <- cMA
            },
            list = {
              cM <- cMA[[j]]
            })

          ftotal <- Matrix::crossprod(Beta[, j], mm(cM, Beta[, j]))
          fixed1[j] <- fixed1[j] + ftotal

          for (k in seq_len(ngroups)) {
            sel <- (group == k)
            fpart <- Matrix::crossprod(
              Beta[sel, j], mm(cM[sel, sel], Beta[sel, j]))
            fixedsplit1[j, k] <- fixedsplit1[j, k] + fpart
          }

        }

        for (level in seq_len(nr)) {
          Lambda <- Lambdas[[level]]
          nf <- dim(Lambda)[[1]]
          for (factor in seq_len(nf)) {
            random1[, level] <- random1[, level] +
              t(Lambda[factor, ]) * Lambda[factor, ]
          }
        }


        if (nr > 0) {
          tot <- fixed1 + rowSums(random1)
          fixed <- fixed1 / tot
          for (level in seq_len(nr)) {
            random <- random1[, level] / tot
          }
        } else {
          fixed <- matrix(1, nrow = ns, ncol = 1)
          random <- random1
        }

        # fixedsplit <- matrix(0, nrow = ns, ncol = ngroups)
        # for (k in seq_len(ngroups)) {
        #   fixedsplit[, k] <- fixedsplit1[, k] / rowSums(fixedsplit1)
        # }
        fixedsplit <- fixedsplit1 / rowSums(fixedsplit1)

        return(
          list(
            fixed = fixed, random = random, fixedsplit = fixedsplit,
            R2T.Y = R2T.Y, R2T.Beta = R2T.Beta))

      },
      future.scheduling = Inf, future.seed = TRUE,
      future.packages = c(
        "Matrix", "dplyr", "arrow", "purrr", "IASDT.R", "qs2", "methods"),
      future.globals = c(
        "ngroups", "Files_la", "Files_lf", "Files_lmu", "path_lbeta", "nc",
        "Model", "path_postList", "poolN", "ns", "nr", "cMA", "group"))

    # stopping the cluster
    IASDT.R::set_parallel(stop_cluster = TRUE, level = 2)

    # Summarize the results
    IASDT.R::cat_time("Summarize the results", level = 1)
    fixed <- Reduce("+", purrr::map(Res, ~ .x$fixed)) / poolN
    random <- Reduce("+", purrr::map(Res, ~ .x$random)) / poolN
    fixedsplit <- Reduce("+", purrr::map(Res, ~ .x$fixedsplit)) / poolN
    R2T.Y <- Reduce("+", purrr::map(Res, ~ .x$R2T.Y)) / poolN
    R2T.Beta <- Reduce("+", purrr::map(Res, ~ .x$R2T.Beta)) / poolN

  } else {

    IASDT.R::cat_time("Computing variance partitioning sequentially")

    mm <- methods::getMethod("%*%", "Matrix")

    fixed <- matrix(0, nrow = ns, ncol = 1)
    fixedsplit <- matrix(0, nrow = ns, ncol = ngroups)
    random <- matrix(0, nrow = ns, ncol = nr)
    R2T.Y <- 0
    R2T.Beta <- rep(0, nc)

    for (i in seq_len(poolN)) {

      if (i %% 200 == 0) {
        IASDT.R::cat_time(sprintf("Processing iteration %d of %d", i, poolN))
      }

      DT_la <- la[[i]]
      DT_lf <- lf[[i]]
      DT_lmu <- lmu[[i]]

      # Suppress warnings when no trait information is used in the models
      # cor(Beta[k, ], lmu[k, ]) : the standard deviation is zero
      suppressWarnings({
        for (k in seq_len(nc)) {
          R2T.Beta[k] <- R2T.Beta[k] +
            stats::cor(lbeta[[i]][k, ], DT_lmu[k, ])^2
        }
      })

      fixed1 <- matrix(0, nrow = ns, ncol = 1)
      fixedsplit1 <- matrix(0, nrow = ns, ncol = ngroups)
      random1 <- matrix(0, nrow = ns, ncol = nr)
      Beta <- postList[[i]]$Beta
      Lambdas <- postList[[i]]$Lambda

      a <- DT_la - matrix(rep(rowMeans(DT_la), ns), ncol = ns)
      f <- DT_lf - matrix(rep(rowMeans(DT_lf), ns), ncol = ns)

      res1 <- sum((rowSums((a * f)) / (ns - 1))^2)
      res2 <- sum((rowSums((a * a)) / (ns - 1)) *
                    (rowSums((f * f)) / (ns - 1)))
      R2T.Y <- R2T.Y + res1 / res2

      for (j in seq_len(ns)) {
        switch(
          class(Model$X)[1L],
          matrix = {
            cM <- cMA
          },
          list = {
            cM <- cMA[[j]]
          })

        ftotal <- Matrix::crossprod(Beta[, j], mm(cM, Beta[, j]))
        fixed1[j] <- fixed1[j] + ftotal

        for (k in seq_len(ngroups)) {
          sel <- (group == k)
          fpart <- Matrix::crossprod(
            Beta[sel, j], mm(cM[sel, sel], Beta[sel, j]))
          fixedsplit1[j, k] <- fixedsplit1[j, k] + fpart
        }
      }

      for (level in seq_len(nr)) {
        Lambda <- Lambdas[[level]]
        nf <- dim(Lambda)[[1]]
        for (factor in seq_len(nf)) {
          random1[, level] <- random1[, level] +
            t(Lambda[factor, ]) * Lambda[factor, ]
        }
      }

      if (nr > 0) {
        tot <- fixed1 + rowSums(random1)
        fixed <- fixed + fixed1 / tot
        for (level in seq_len(nr)) {
          random[, level] <- random[, level] + random1[, level] / tot
        }
      } else {
        fixed <- fixed + matrix(1, nrow = ns, ncol = 1)
      }

      for (k in seq_len(ngroups)) {
        fixedsplit[, k] <- fixedsplit[, k] +
          fixedsplit1[, k] / rowSums(fixedsplit1)
      }
    }

    fixed <- fixed / poolN
    random <- random / poolN
    fixedsplit <- fixedsplit / poolN
    R2T.Y <- R2T.Y / poolN
    R2T.Beta <- R2T.Beta / poolN
  }

  # # .................................................................... ###

  vals <- matrix(0, nrow = ngroups + nr, ncol = ns)
  for (i in seq_len(ngroups)) {
    vals[i, ] <- fixed * fixedsplit[, i]
  }
  for (i in seq_len(nr)) {
    vals[ngroups + i, ] <- random[, i]
  }

  names(R2T.Beta) <- Model$covNames
  colnames(vals) <- Model$spNames
  leg <- group_names
  for (r in seq_len(nr)) {
    leg <- c(leg, paste0("Random: ", Model$rLNames[r]))
  }
  rownames(vals) <- leg

  VP <- list()
  VP$vals <- vals
  VP$R2T <- list(Beta = R2T.Beta, Y = R2T.Y)
  VP$group <- group
  VP$groupnames <- group_names

  # # .................................................................... ###

  # Save the results
  IASDT.R::cat_time("Save the variance partitioning results")

  File_VarPar <- IASDT.R::path(Path_VarPar, paste0(VP_file, ".RData"))
  IASDT.R::save_as(object = VP, object_name = VP_file, out_path = File_VarPar)

  VP$File <- File_VarPar

  # # .................................................................... ###

  if (temp_cleanup) {

    IASDT.R::cat_time("Clean up temporary files")

    if (use_TF) {

      if (fs::dir_exists(path_lbeta)) {
        try({
          IASDT.R::system_command(
            paste0("rm -rf ", IASDT.R::normalize_path(path_lbeta)),
            ignore.stderr = TRUE, ignore.stdout = TRUE)
        },
        silent = TRUE)
      }
      if (fs::dir_exists(path_postList)) {
        try({
          IASDT.R::system_command(
            paste0("rm -rf ", IASDT.R::normalize_path(path_postList)),
            ignore.stderr = TRUE, ignore.stdout = TRUE)
        },
        silent = TRUE)
      }
    }

    Path_Temp <- IASDT.R::normalize_path(Path_Temp)
    try(
      expr = {
        file_paths <- list.files(
          path = Path_Temp,
          pattern = "(VP_).+(feather|log)$", full.names = TRUE)

        fs::file_delete(file_paths)
      },
      silent = TRUE)

    if (fs::dir_exists(Path_Temp) && length(fs::dir_ls(Path_Temp)) == 0) {
      try({
        IASDT.R::system_command(
          paste0("rm -rf ", Path_Temp),
          ignore.stderr = TRUE, ignore.stdout = TRUE)
      },
      silent = TRUE)
    }
  }

  # # .................................................................... ###

  IASDT.R::cat_diff(init_time = .StartTime)

  return(VP)
}
