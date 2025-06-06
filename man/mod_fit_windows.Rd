% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mod_fit_windows.R
\name{mod_fit_windows}
\alias{mod_fit_windows}
\title{Fit Hmsc-HPC models on UFZ Windows Server}
\usage{
mod_fit_windows(
  path_model = NULL,
  python_VE = NULL,
  n_cores = NULL,
  strategy = "multisession"
)
}
\arguments{
\item{path_model}{Character. Path to the model files. This argument can not
be empty.}

\item{python_VE}{Character. Path to a valid Python virtual environment.
Defaults to \code{NULL}. This argument can not be empty.}

\item{n_cores}{Integer. Number of CPU cores to use for parallel processing.}

\item{strategy}{Character. The parallel processing strategy to use. Valid
options are "sequential", "multisession" (default), "multicore", and
"cluster". See \code{\link[future:plan]{future::plan()}} and \code{\link[ecokit:set_parallel]{ecokit::set_parallel()}} for details.}
}
\value{
The function does not return anything but prints messages to the
console regarding the progress and completion of model fitting.
}
\description{
This function fits Hmsc models on a UFZ Windows Server. It reads model
configurations from a specified path, loads environment variables, checks
input arguments for validity, and executes model fitting in parallel if
required.
}
\author{
Ahmed El-Gabbas
}
