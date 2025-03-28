% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mod_RespCurv_PrepData.R,
%   R/Mod_RespCurv_PlotSp.R, R/Mod_RespCurv_PlotSpAll.R, R/Mod_RespCurv_PlotSR.R
\name{Response_curves}
\alias{Response_curves}
\alias{RespCurv_PrepData}
\alias{RespCurv_PlotSp}
\alias{RespCurv_PlotSpAll}
\alias{RespCurv_PlotSR}
\title{Prepare and plot response curve data for Hmsc models}
\usage{
RespCurv_PrepData(
  Path_Model = NULL,
  N_Grid = 50,
  NCores = 8,
  ReturnData = FALSE,
  Probabilities = c(0.025, 0.5, 0.975),
  UseTF = TRUE,
  TF_Environ = NULL,
  TF_use_single = FALSE,
  LF_NCores = NCores,
  LF_Check = FALSE,
  LF_Temp_Cleanup = TRUE,
  LF_Commands_Only = FALSE,
  Temp_Dir = "TEMP_Pred",
  Temp_Cleanup = TRUE,
  Verbose = TRUE
)

RespCurv_PlotSp(
  ModelDir = NULL,
  NCores = 20,
  EnvFile = ".env",
  ReturnData = FALSE
)

RespCurv_PlotSpAll(
  ModelDir = NULL,
  NCores = 8L,
  ReturnData = FALSE,
  PlottingAlpha = 0.3
)

RespCurv_PlotSR(ModelDir, Verbose = TRUE, NCores = 8L)
}
\arguments{
\item{Path_Model}{Character. Path to the file containing the fitted Hmsc
model.}

\item{N_Grid}{Integer. Number of points along the gradient for continuous
focal variables. Higher values result in smoother curves. Default: 50. See
\link[Hmsc:constructGradient]{Hmsc::constructGradient} for details.}

\item{NCores}{Integer. Number of CPU cores to use for parallel processing.
Defaults to 8 for all functions, except for \code{RespCurv_PlotSp}, in which it
defaults to 20.}

\item{ReturnData}{Logical. If \code{TRUE}, the function returns processed data as
an R object. Default: \code{FALSE}.}

\item{Probabilities}{Numeric vector. Quantiles to calculate in response curve
predictions. Default: \code{c(0.025, 0.5, 0.975)}. See \link[stats:quantile]{stats::quantile} for
details.}

\item{UseTF}{Logical. Whether to use TensorFlow for calculations. Defaults to
\code{TRUE}.}

\item{TF_Environ}{Character. Path to the Python environment. This argument is
required if \code{UseTF} is \code{TRUE} under Windows. Defaults to \code{NULL}.}

\item{TF_use_single}{Logical. Whether to use single precision for the
TensorFlow calculations. Defaults to \code{FALSE}.}

\item{LF_NCores}{Integer. Number of cores to use for parallel processing of
latent factor prediction. Defaults to 8L.}

\item{LF_Check}{Logical. If \code{TRUE}, the function checks if the output files
are already created and valid. If \code{FALSE}, the function will only check if
the files exist without checking their integrity. Default is \code{FALSE}.}

\item{LF_Temp_Cleanup}{Logical. Whether to delete temporary files in the
\code{Temp_Dir} directory after finishing the LF predictions.}

\item{LF_Commands_Only}{Logical. If \code{TRUE}, returns the command to run the
Python script. Default is \code{FALSE}.}

\item{Temp_Dir}{Character. Path for temporary storage of intermediate files.}

\item{Temp_Cleanup}{Logical. Whether to clean up temporary files.
Defaults to \code{TRUE}.}

\item{Verbose}{Logical. Whether to print a message upon
successful saving of files. Defaults to \code{FALSE}.}

\item{ModelDir}{Character. Path to the root directory containing fitted
models. The function reads data from the \code{RespCurv_DT} subdirectory, which
is created by \code{RespCurv_PrepData}.}

\item{EnvFile}{Character. Path to the environment file containing paths to
data sources. Defaults to \code{.env}.}

\item{PlottingAlpha}{Numeric. Opacity level for response curve lines (0 =
fully transparent, 1 = fully opaque). Default: 0.3.}
}
\description{
The \verb{RespCurv_*()} functions process and visualize response curves for Hmsc
models. They support parallel computation and optionally return processed
data. There are four functions in this group:
\itemize{
\item \code{RespCurv_PrepData()}: Prepares response curve data for analysis
\item \code{RespCurv_PlotSp()}: Generates response curve plots for individual species
\item \code{RespCurv_PlotSpAll()}: Generates response curves for all species together
in a single plot
\item \code{RespCurv_PlotSR()}: Plots response curves for species richness.
}
}
\author{
Ahmed El-Gabbas
}
