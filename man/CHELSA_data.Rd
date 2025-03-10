% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DWF_CHELSA_Process.R, R/DWF_CHELSA_Prepare.R,
%   R/DWF_CHELSA_Project.R
\name{CHELSA_data}
\alias{CHELSA_data}
\alias{CHELSA_Process}
\alias{CHELSA_Prepare}
\alias{CHELSA_Project}
\title{Process CHELSA Climate Data for the \code{IAS-pDT}}
\usage{
CHELSA_Process(
  EnvFile = ".env",
  NCores = 8L,
  Download = FALSE,
  Overwrite = FALSE,
  Download_Attempts = 10L,
  Sleep = 5L,
  OtherVars = "npp",
  Download_NCores = 4,
  CompressLevel = 5,
  OverwriteProcessed = FALSE
)

CHELSA_Prepare(
  EnvFile = ".env",
  Download = FALSE,
  NCores = 8L,
  Overwrite = FALSE,
  Download_Attempts = 10L,
  Sleep = 5L,
  OtherVars = "npp"
)

CHELSA_Project(Metadata = NULL, EnvFile = ".env", CompressLevel = 5)
}
\arguments{
\item{EnvFile}{Character. Path to the environment file containing paths to
data sources. Defaults to \code{.env}.}

\item{NCores}{Integer. Number of CPU cores to use for parallel processing.
Default: 8.}

\item{Download}{Logical. If \code{TRUE}, downloads CHELSA files. Default: \code{FALSE}.}

\item{Overwrite}{Logical. If \code{TRUE}, re-downloads existing files. Default:
\code{FALSE}.}

\item{Download_Attempts}{Integer. Maximum download retries. Default: \code{10}.}

\item{Sleep}{Integer. Seconds to wait between download attempts. Default:
\code{5}.}

\item{OtherVars}{Character. Additional variables to process (e.g., \code{"npp"}
for Net Primary Productivity alongside 19 bioclimatic variables
bio1-bio19). Use \code{""} for bioclimatic only. See \link{CHELSA_Vars} for details.
Default: \code{"npp"}.}

\item{Download_NCores}{Integer. Number of CPU cores to use for parallel
downloading of CHELSA data. Only valid if Download = \code{TRUE}. Defaults to 4.}

\item{CompressLevel}{Integer. NetCDF compression level (1 = least, 9 = most).
Default: \code{5}.}

\item{OverwriteProcessed}{Logical. If \code{TRUE}, overwrites processed files.
Default: \code{FALSE}.}

\item{Metadata}{Tibble. Single-row metadata for input files, prepared by
\code{CHELSA_Prepare()}}
}
\description{
Downloads, processes, and projects \href{https://chelsa-climate.org/}{CHELSA}
climate data at the European scale for the Invasive Alien Species prototype
Digital Twin (\code{IAS-pDT}). Supports multiple climate scenarios, outputting
data in TIFF and NetCDF formats. Orchestrated by \code{CHELSA_Process()}, with
helper functions \code{CHELSA_Prepare()} and \code{CHELSA_Project()}.
}
\note{
\itemize{
\item \code{CHELSA_Prepare()} and \code{CHELSA_Project()} are internal helpers, not for
direct use.
\item Processes 19 bioclimatic variables (bio1–bio19) plus optional variables
(e.g., NPP) for 46 scenarios (1 current, 45 future).
\item Time-intensive; depends on file size and compute resources.
}
}
\section{Functions details}{

\itemize{
\item \strong{\code{CHELSA_Process()}}: Main function; optionally downloads CHELSA data,
processes it to the European scale and reference grid, and saves TIFF and
NetCDF outputs for 46 climate scenarios.
\item \strong{\code{CHELSA_Prepare()}}: Extracts metadata from local URL text files and
manages optional downloads.
\item \strong{\code{CHELSA_Project()}}: Projects data to the IAS-pDT reference grid with
optional transformations.
}
}

\author{
Ahmed El-Gabbas
}
