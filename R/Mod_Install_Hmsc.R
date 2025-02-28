## |------------------------------------------------------------------------| #
# Install_Hmsc ----
## |------------------------------------------------------------------------| #

#' Install Hmsc-HPC in a python virtual environment on Windows
#'
#' This function sets up a Python virtual environment and installs the
#' [`Hmsc-HPC`](https://github.com/aniskhan25/hmsc-hpc) package from a specified
#' Git repository. It also installs `TensorFlow` and performs checks to to
#' verify the installation of Python, the virtual environment, and the packages.
#' @param Path_Python Character. Path to the Python executable.
#' @param Path_VE Character. Path where the virtual environment will be created.
#'   This can not be an existing folder.
#' @param URL_Hmsc Character. URL of the Git repository (and branch name) of the
#'   `Hmsc-HPC` package.
#' @param URL_rdata Character. URL of the Git repository (and branch name) of
#'   the `rdata` package. This is temporary to allow `rdata` package to write
#'   rds file. In the near future, this functionality will be pushed to the
#'   main branch of `rdata`.
#' @param Force Logical. Whether to force the installation of `Hmsc-HPC`; i.e.
#'   using `--force-reinstall` as a suffix to the `pip install` command
#' @return The function performs installation steps and returns NULL invisibly.
#' @export
#' @name Install_Hmsc
#' @author Ahmed El-Gabbas
#' @details The function performs the following steps:
#' - Checks if the virtual environment directory already exists and stops with
#'   an error if it does
#' - Verifies the Python version and installation
#' - Creates a new Python virtual environment
#' - Upgrades `pip` in the virtual environment
#' - Installs the `Hmsc-HPC` package from the provided Git URL
#' - Installs `TensorFlow` version 2.15
#' - Installs `rdata`
#' - Checks `TensorFlow` and `Hmsc-HPC` package installations
#' @examples
#' \dontrun{
#' Install_Hmsc(
#'    Path_Python = "C:/Python/Python310/python.exe",
#'    Path_VE = "D:/Hmsc-HPC")
#' }

Install_Hmsc <- function(
    Path_Python, Path_VE,
    URL_Hmsc = "https://github.com/trossi/hmsc-hpc.git@simplify-io",
    URL_rdata = "https://github.com/trossi/rdata.git@test-for-hmsc-v2",
    Force = FALSE) {

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Check if the virtual environment directory already exists
  if (fs::dir_exists(Path_VE)) {
    stop(
      "Path to the virtual environment already exists:\n", Path_VE,
      call. = FALSE)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Checking Python installation

  IASDT.R::InfoChunk("Checking Python")
  PythonVersion <- system2(
    Path_Python, "--version", stdout = TRUE, stderr = TRUE)

  if (!stringr::str_detect(PythonVersion, "^Python")) {
    stop(
      "Python was not installed or not found at the provided path",
      Path_Python, call. = FALSE)
  }

  cat(paste0("  >>  ", PythonVersion, "\n"))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Creating the virtual environment

  IASDT.R::InfoChunk("Creating virtual environment")
  system2(Path_Python, paste0("-m venv ", Path_VE))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Checking Python version in the virtual environment

  IASDT.R::InfoChunk("Checking virtual environment")
  python <- IASDT.R::Path(Path_VE, "Scripts", "python.exe")
  PythonVersion <- system2(python, "--version", stdout = TRUE, stderr = TRUE)

  if (!stringr::str_detect(PythonVersion, "^Python")) {
    stop("Python was not installed in the virtual environment", call. = FALSE)
  }

  cat(paste0("  >>  ", PythonVersion, "\n"))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Upgrading pip

  IASDT.R::InfoChunk("Upgrading pip")
  system2(python, "-m pip install --upgrade pip")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Installing Hmsc-HPC package from Git

  IASDT.R::InfoChunk(
    "Installing Hmsc-HPC package from Git")
  if (Force) {
    Command <- paste0("-m pip install git+", URL_Hmsc, " --force-reinstall")
  } else {
    Command <- paste0("-m pip install git+", URL_Hmsc)
  }
  system2(python, Command)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Install TensorFlow v2.15

  # This is to avoid the following error: ImportError: This version of
  # TensorFlow Probability requires TensorFlow version >= 2.15; Detected an
  # installation of version 2.13.1. Please upgrade TensorFlow to proceed.

  IASDT.R::InfoChunk("Installing TensorFlow v2.15")
  system2(python, "-m pip install tensorflow==2.15")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Checking TensorFlow
  IASDT.R::InfoChunk("Checking TensorFlow")
  system2(python, '-c "import tensorflow as tf; print(tf.constant(1))"')

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Install `rdata`
  IASDT.R::InfoChunk("Installing rdata")
  Command <- paste0("-m pip install git+", URL_rdata)
  system2(python, Command)

  # It is necessary to downgrade the version of numpy to work with tensorflow

  # ERROR: pip's dependency resolver does not currently take into account all
  # the packages that are installed. This behaviour is the source of the
  # following dependency conflicts. tensorflow-intel 2.15.0 requires
  # numpy<2.0.0, >=1.23.5, but you have numpy 2.0.2 which is incompatible.

  system2(python, '-m pip install numpy==1.26.4"')

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Checking Hmsc
  IASDT.R::InfoChunk("Checking Hmsc")
  system2(python, '-c "import hmsc"')

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  return(invisible(NULL))
}
