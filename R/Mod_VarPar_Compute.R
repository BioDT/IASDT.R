## |------------------------------------------------------------------------| #
# VarPar_Compute ----
## |------------------------------------------------------------------------| #

#' Computes variance partitioning of Hmsc models
#'
#' Computes variance components with respect to given grouping of fixed effects
#' and levels of random effects. This function is a wrapper around the
#' `Hmsc::computeVariancePartitioning`, but with the added functionality of
#' parallel computation using TensorFlow.
#'
#' @param Path_Model Character. Path to fitted `Hmsc` model object.
#' @param group vector of numeric values corresponding to group identifiers in
#'   groupnames. If the model was defined with `XData` and `XFormula`, the
#'   default is to use model terms.
#' @param groupnames vector of names for each group of fixed effect. Should
#'   match `group`. If the model was defined with `XData` and `XFormula`, the
#'   default is to use the labels of model terms.
#' @param start index of first MCMC sample included
#' @param na.ignore logical. If TRUE, covariates are ignored for sites where the
#'   focal species is NA when computing variance-covariance matrices for each
#'   species
#' @param NCores Integer specifying the number of parallel cores for
#'   parallelization. This is only used when `UseTF` is TRUE.
#' @param Chunk_size Integer. Size of each chunk of samples to process in
#'   parallel. Only relevant for TensorFlow. Default: `50`.
#' @param Verbose Logical. Indicates whether progress messages should be
#'   displayed. Defaults to `TRUE`.
#' @param Temp_Cleanup Logical. Whether to delete temporary files after
#'   processing. Default: `TRUE`.
#' @param VP_Commands_Only logical. If `TRUE`, returns the commands to run the
#'   Python script. Default is `FALSE`. Only relevant when `UseTF` is `TRUE`.
#' @param VarParFile Character. Name of the output file to save the results.
#' @export
#' @name VarPar_Compute
#' @author Ahmed El-Gabbas
#' @inheritParams Predict_Hmsc
#' @inheritParams Predict_LF
#' @export

VarPar_Compute <- function(
    Path_Model, group = NULL, groupnames = NULL, start = 1L, na.ignore = FALSE,
    NCores = 8L, UseTF = TRUE, TF_Environ = NULL, TF_use_single = FALSE,
    Temp_Cleanup = TRUE, Chunk_size = 50L, Verbose = TRUE,
    VarParFile = "VarPar", VP_Commands_Only = FALSE) {

  # # .................................................................... ###

  if (isFALSE(Verbose)) {
    sink(file = nullfile())
    on.exit(sink(), add = TRUE)
  }

  .StartTime <- lubridate::now(tzone = "CET")

  # # .................................................................... ###

  # Check if the virtual environment and Python scripts exist

  if (UseTF) {
    if (isFALSE(VP_Commands_Only) &&
        (is.null(TF_Environ) || !dir.exists(TF_Environ))) {
      stop(
        paste0(
          "When `UseTF` is TRUE, `TF_Environ` must be specified and should ",
          "point to an existing directory with a Python environment"),
        call. = FALSE)
    }

    # Determine the Python executable path
    python_executable <- if (.Platform$OS.type == "windows") {
      if (isFALSE(VP_Commands_Only) && !file.exists(python_executable)) {
        stop(
          "Python executable not found in the virtual environment.",
          call. = FALSE)
      }
      file.path(TF_Environ, "Scripts", "python.exe")
    } else {
      "/usr/bin/time -v python3"
    }

    # Check GPU availability
    if (isFALSE(VP_Commands_Only) && .Platform$OS.type == "windows") {
      result <- system(
        paste0(
          python_executable,
          " -c \"import tensorflow as tf; print(len(",
          "tf.config.list_physical_devices('GPU')))\""),
        intern = TRUE)
      N_GPU <- result[length(result)]
      if (N_GPU == 0) {
        IASDT.R::CatTime(
          "No GPU found; Calculations will use CPU.",
          Time = FALSE, Bold = TRUE, Red = TRUE)
      } else {
        IASDT.R::CatTime(
          paste0(N_GPU, " GPUs were found. Calculations will use GPU."),
          Time = FALSE, Bold = TRUE, Red = TRUE)
      }
    }

    # Paths to the Python scripts
    Script_geta <- system.file("VP_geta.py", package = "IASDT.R")
    Script_getf <- system.file("VP_getf.py", package = "IASDT.R")
    Script_gemu <- system.file("VP_gemu.py", package = "IASDT.R")
    if (!all(file.exists(c(Script_geta, Script_getf, Script_gemu)))) {
      stop(
        "Necessary python scripts do not exist", call. = FALSE)
    }
  }

  # # .................................................................... ###

  # Load model object ------
  IASDT.R::CatTime("Load model object")

  if (is.null(Path_Model) || !file.exists(Path_Model)) {
    stop(
      paste0("Model path is NULL or does not exist: ", Path_Model),
      call. = FALSE)
  }

  Model <- IASDT.R::LoadAs(Path_Model)

  # # .................................................................... ###

  # Create folder for variance partitioning results
  Path_VarPar <- file.path(
    dirname(dirname(Path_Model)), "Model_Postprocessing/Variance_Partitioning")
  fs::dir_create(Path_VarPar)

  # # .................................................................... ###

  ny <- Model$ny
  nc <- Model$nc
  ns <- Model$ns
  nr <- Model$nr

  if (is.null(group)) {
    # names of variables used in the model
    groupnames <- names(Model$XData)

    # actual variables used in the model, including the intercept and quadratic
    # terms
    ModelVars <- dimnames(Model$X)[[2]][-1]

    # group variable to combine variable and its quadratic terms together
    group <- purrr::map(
      ModelVars, ~ which(stringr::str_detect(.x, groupnames))) %>%
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

  IASDT.R::CatTime("Prepare postList")
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

  if (UseTF) {

    # Create the temporary directory
    Path_Temp <- file.path(dirname(dirname(Path_Model)), "TEMP_VP")
    fs::dir_create(Path_Temp)

    FileSuffix <- stringr::str_pad(
      string = seq_len(poolN), pad = "0", width = 4)

    # List of feather files resulted from `geta` function
    Files_la <- file.path(Path_Temp, paste0("VP_A_", FileSuffix, ".feather"))
    Files_la_Exist <- all(file.exists(Files_la))

    # List of feather files resulted from `getf` function
    Files_lf <- file.path(Path_Temp, paste0("VP_F_", FileSuffix, ".feather"))
    Files_lf_Exist <- all(file.exists(Files_lf))

    # List of feather files resulted from `gemu` function
    Files_lmu <- file.path(Path_Temp, paste0("VP_Mu_", FileSuffix, ".feather"))
    Files_lmu_Exist <- all(file.exists(Files_lmu))

    Beta_Files <- file.path(
      Path_Temp, paste0("VP_Beta_", FileSuffix, ".feather"))

    if (all(c(Files_la_Exist, Files_lf_Exist, Files_lmu_Exist))) {

      # All feather data are already processed and available on disk
      IASDT.R::CatTime(
        "Data for `la`/`lf`/`lmu` lists were already processed")

    } else {

      ## Prepare la/lf/lmu lists using TensorFlow ----
      IASDT.R::CatTime("Prepare la/lf/lmu lists using TensorFlow")

      ### Prepare data for TensorFlow ----
      IASDT.R::CatTime("Prepare data for TensorFlow", Level = 1)

      #### X data -----
      # needed only to calculate `geta` and `getf` functions
      if (!all(c(Files_la_Exist, Files_lf_Exist))) {
        IASDT.R::CatTime("X", Level = 2)
        Path_X <- file.path(Path_Temp, "VP_X.feather")
        if (!file.exists(Path_X)) {
          arrow::write_feather(as.data.frame(Model$X), Path_X)
        }
      }

      #### Tr / Gamma ------
      # needed only to calculate `geta` and `gemu` functions
      if (!all(c(Files_la_Exist, Files_lmu_Exist))) {
        IASDT.R::CatTime("Tr", Level = 2)
        Path_Tr <- file.path(Path_Temp, "VP_Tr.feather")
        if (!file.exists(Path_Tr)) {
          arrow::write_feather(as.data.frame(Model$Tr), Path_Tr)
        }

        # Gamma - convert each list item into a column in a data frame
        IASDT.R::CatTime("Gamma", Level = 2)
        Path_Gamma <- file.path(Path_Temp, "VP_Gamma.feather")
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
        IASDT.R::CatTime("Beta", Level = 2)

        if (!all(file.exists(Beta_Files))) {



          IASDT.R::CatTime("Prepare working on parallel", Level = 3)
          c1 <- parallel::makePSOCKcluster(NCores)
          on.exit(try(parallel::stopCluster(c1), silent = TRUE), add = TRUE)

          IASDT.R::CatTime("Export necessary objects to cores", Level = 3)
          parallel::clusterExport(
            cl = c1, varlist = c("Beta_Files", "postList"),
            envir = environment())

          IASDT.R::CatTime("Load libraries on each core", Level = 3)
          invisible(
            parallel::clusterEvalQ(
              cl = c1, expr = sapply("arrow", library, character.only = TRUE)))

          IASDT.R::CatTime("Processing beta on parallel", Level = 3)
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

          IASDT.R::CatTime("stop cluster", Level = 3)
          parallel::stopCluster(c1)
        }
      }
      invisible(gc())
      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||

      ### Processing geta -----

      if (!Files_la_Exist) {

        IASDT.R::CatTime("Processing `geta` function", Level = 1)
        Path_Out_a <- file.path(Path_Temp, "VP_A.feather") %>%
          normalizePath(winslash = "/", mustWork = FALSE)

        cmd_a <- paste(
          python_executable, Script_geta,
          "--tr", normalizePath(Path_Tr, winslash = "/", mustWork = TRUE),
          "--x", normalizePath(Path_X, winslash = "/", mustWork = TRUE),
          "--gamma", normalizePath(Path_Gamma, winslash = "/", mustWork = TRUE),
          "--output", Path_Out_a,
          "--ncores", NCores,
          "--chunk_size", Chunk_size)

        if (TF_use_single) {
          cmd_a <- paste0(cmd_a, " --use_single")
        }

        if (VP_Commands_Only) {

          # Save command to file
          # Redirect results of time to log file
          path_log_a <- stringr::str_replace(Path_Out_a, ".feather", ".log")
          cmd_a <- paste0(cmd_a, paste0(" >> ", path_log_a, " 2>&1"))
          readr::write_lines(
            x = cmd_a,
            file = file.path(Path_Temp, "VP_A_Command.txt"), append = FALSE)

        } else {
          # Run the command using system
          la <- system(cmd_a, wait = TRUE, intern  = TRUE)

          # Check for errors
          if (inherits(la, "error") || la[length(la)] != "Done") {
            stop(paste0("Error in computing geta: ", la), call. = FALSE)
          }

          if (length(la) != 1) {
            cat(la, sep = "\n")
          }
        }

      } else {
        IASDT.R::CatTime(
          "All `la` data were already available on disk", Level = 1)
      }

      invisible(gc())

      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||

      ### Processing getf ----

      if (!Files_la_Exist) {

        IASDT.R::CatTime("Processing `getf` function", Level = 1)
        Path_Out_f <- file.path(Path_Temp, "VP_F.feather") %>%
          normalizePath(winslash = "/", mustWork = FALSE)

        cmd_f <- paste(
          python_executable, Script_getf,
          "--x", normalizePath(Path_X, winslash = "/", mustWork = TRUE),
          "--beta_dir",
          normalizePath(Path_Temp, winslash = "/", mustWork = TRUE),
          "--output", Path_Out_f,
          "--ncores", NCores)

        if (TF_use_single) {
          cmd_f <- paste0(cmd_f, " --use_single")
        }

        if (VP_Commands_Only) {

          # Save command to file
          # Redirect results of time to log file
          path_log_f <- stringr::str_replace(Path_Out_f, ".feather", ".log")
          cmd_f <- paste0(cmd_f, paste0(" >> ", path_log_f, " 2>&1"))
          readr::write_lines(
            x = cmd_f,
            file = file.path(Path_Temp, "VP_F_Command.txt"), append = FALSE)

        } else {

          # Run the command using system
          lf <- system(cmd_f, wait = TRUE, intern  = TRUE)

          # Check for errors
          if (inherits(lf, "error") || lf[length(lf)] != "Done") {
            stop(paste0("Error in computing geta: ", lf), call. = FALSE)
          }

          if (length(lf) != 1) {
            cat(lf, sep = "\n")
          }
        }

      } else {
        IASDT.R::CatTime(
          "All `lf` data were already available on disk", Level = 1)
      }

      invisible(gc())

      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||

      ### Processing gemu ----

      if (!Files_lmu_Exist) {

        IASDT.R::CatTime("Processing `gemu` function", Level = 1)
        Path_Out_mu <- file.path(Path_Temp, "VP_Mu.feather") %>%
          normalizePath(winslash = "/", mustWork = FALSE)

        cmd_mu <- paste(
          python_executable, Script_gemu,
          "--tr", normalizePath(Path_Tr, winslash = "/", mustWork = TRUE),
          "--gamma", normalizePath(Path_Gamma, winslash = "/", mustWork = TRUE),
          "--output", Path_Out_mu,
          "--ncores", NCores,
          "--chunk_size", Chunk_size)

        if (TF_use_single) {
          cmd_mu <- paste0(cmd_mu, " --use_single")
        }

        if (VP_Commands_Only) {

          # Save command to file
          # Redirect results of time to log file
          path_log_mu <- stringr::str_replace(Path_Out_mu, ".feather", ".log")
          cmd_mu <- paste0(cmd_mu, paste0(" >> ", path_log_mu, " 2>&1"))
          readr::write_lines(
            x = cmd_mu,
            file = file.path(Path_Temp, "VP_mu_Command.txt"), append = FALSE)

        } else {

          # Run the command using system
          lmu <- system(cmd_mu, wait = TRUE, intern  = TRUE)

          # Check for errors
          if (inherits(lmu, "error") || lmu[length(lmu)] != "Done") {
            stop(paste0("Error in computing geta: ", lmu), call. = FALSE)
          }

          if (length(lmu) != 1) {
            cat(lmu, sep = "\n")
          }
        }
      } else {
        IASDT.R::CatTime(
          "All `lmu` data were already available on disk", Level = 1)
      }
    }

    invisible(gc())

    if (VP_Commands_Only) {
      return(invisible(NULL))
    }

  } else {

    # Prepare la/lf/lmu lists using original R code -----

    IASDT.R::CatTime("Prepare la/lf/lmu lists using original R code")

    ## geta -----
    IASDT.R::CatTime("Running geta", Level = 1)
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

    IASDT.R::CatTime("Running getf", Level = 1)
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

    IASDT.R::CatTime("Running gemu", Level = 1)
    gemu <- function(a) {
      res <- t(Model$Tr %*% t(a$Gamma))
      return(res)
    }

    lmu <- lapply(postList, gemu)

    invisible(gc())

  }

  # # .................................................................... ###

  # Running gebeta -------
  IASDT.R::CatTime("Running gebeta")
  gebeta <- function(a) {
    res <- a$Beta
    return(res)
  }
  lbeta <- lapply(postList, gebeta)

  # # .................................................................... ###

  # Remove Gamma from postList ------
  IASDT.R::CatTime("Remove Gamma from postList")
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


  mm <- methods::getMethod("%*%", "Matrix")

  if (UseTF) {

    IASDT.R::CatTime("Computing variance partitioning in parallel")

    IASDT.R::CatTime("Prepare working on parallel", Level = 1)
    c1 <- parallel::makePSOCKcluster(NCores)
    on.exit(try(parallel::stopCluster(c1), silent = TRUE), add = TRUE)

    IASDT.R::CatTime("Load libraries", Level = 1)
    invisible(
      parallel::clusterEvalQ(
        cl = c1,
        expr = sapply(
          c("Matrix", "dplyr", "arrow"), library, character.only = TRUE)))

    IASDT.R::CatTime("Export objects to cores", Level = 1)
    parallel::clusterExport(
      cl = c1, varlist = c(
        "ngroups", "Files_la", "Files_lf", "Files_lmu", "lbeta",
        "postList", "poolN", "ns", "nr"),
      envir = environment())

    IASDT.R::CatTime("Processing on parallel", Level = 1)
    Res <- parallel::parLapply(
      cl = c1,
      X = seq_along(postList),
      fun = function(i) {

        DT_la <- as.matrix(arrow::read_feather(Files_la[i]))
        DT_lf <- as.matrix(arrow::read_feather(Files_lf[i]))
        DT_lmu <- as.matrix(arrow::read_feather(Files_lmu[i]))

        # Suppress warnings when no trait information is used in the models
        # cor(Beta[k, ], lmu[k, ]) : the standard deviation is zero
        R2T.Beta <- purrr::map_dbl(
          .x = seq_len(nc),
          .f = ~ {
            suppressWarnings(stats::cor(lbeta[[i]][.x, ], DT_lmu[.x, ])^2)
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
        R2T.Y <- res1 / res2

        for (j in seq_len(ns)) {
          switch(
            class(Model$X)[1L],
            matrix = {
              cM <- cMA
            },
            list = {
              cM <- cMA[[j]]
            }
          )

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

        fixedsplit <- matrix(0, nrow = ns, ncol = ngroups)
        for (k in seq_len(ngroups)) {
          fixedsplit[, k] <- fixedsplit1[, k] / rowSums(fixedsplit1)
        }

        list(
          fixed = fixed, random = random, fixedsplit = fixedsplit,
          R2T.Y = R2T.Y, R2T.Beta = R2T.Beta) %>%
          return()
      }
    )

    # Summarize the results
    IASDT.R::CatTime("Summarize the results", Level = 1)
    fixed <- Reduce("+", purrr::map(Res, ~ .x$fixed)) / poolN
    random <- Reduce("+", purrr::map(Res, ~ .x$random)) / poolN
    fixedsplit <- Reduce("+", purrr::map(Res, ~ .x$fixedsplit)) / poolN
    R2T.Y <- Reduce("+", purrr::map(Res, ~ .x$R2T.Y)) / poolN
    R2T.Beta <- Reduce("+", purrr::map(Res, ~ .x$R2T.Beta)) / poolN

  } else {

    IASDT.R::CatTime("Computing variance partitioning sequentially")

    fixed <- matrix(0, nrow = ns, ncol = 1)
    fixedsplit <- matrix(0, nrow = ns, ncol = ngroups)
    random <- matrix(0, nrow = ns, ncol = nr)
    R2T.Y <- 0
    R2T.Beta <- rep(0, nc)

    for (i in seq_len(poolN)) {

      if (i %% 200 == 0) {
        IASDT.R::CatTime(sprintf("Processing iteration %d of %d", i, poolN))
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
          }
        )

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
  leg <- groupnames
  for (r in seq_len(nr)) {
    leg <- c(leg, paste("Random: ", Model$rLNames[r], sep = ""))
  }
  rownames(vals) <- leg

  VP <- list()
  VP$vals <- vals
  VP$R2T <- list(Beta = R2T.Beta, Y = R2T.Y)
  VP$group <- group
  VP$groupnames <- groupnames

  # # .................................................................... ###

  # Save the results
  IASDT.R::CatTime("Save the variance partitioning results")

  File_VarPar <- file.path(Path_VarPar, paste0(VarParFile, ".RData"))
  IASDT.R::SaveAs(InObj = VP, OutObj = VarParFile, OutPath = File_VarPar)

  VP$File <- File_VarPar

  # # .................................................................... ###

  if (UseTF && Temp_Cleanup) {

    IASDT.R::CatTime("Clean up temporary files")

    Path_Temp <- normalizePath(Path_Temp, winslash = "/", mustWork = FALSE)

    if (dir.exists(Path_Temp)) {
      # Function to delete files or directories
      delete_files <- function(file_path) {
        if (.Platform$OS.type == "windows") {
          # Use rmdir command on Windows
          system2(
            "cmd", c("/c", "rmdir /s /q", shQuote(file_path)),
            stdout = NULL, stderr = NULL)
        } else {
          # Use rm command on Linux/macOS
          system2(
            "rm", c("-rf", shQuote(file_path)),
            stdout = NULL, stderr = NULL)
        }
      }

      # List all files to delete
      files_to_delete <- list.files(Path_Temp, full.names = TRUE)

      if (UseTF) {
        try(parallel::parLapply(c1, files_to_delete, delete_files))

        # stop the cluster
        parallel::stopCluster(c1)
      } else {
        try(purrr::walk(files_to_delete, delete_files))
      }
      fs::dir_delete(Path_Temp)
    }
  }

  # # .................................................................... ###

  IASDT.R::CatDiff(InitTime = .StartTime)

  return(VP)
}
