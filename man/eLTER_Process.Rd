% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DWF_elTER_Process.R
\name{eLTER_Process}
\alias{eLTER_Process}
\title{Process eLTER data for the \code{IAS-pDT}}
\usage{
eLTER_Process(EnvFile = ".env", StartYear = 1981)
}
\arguments{
\item{EnvFile}{Character. Path to the environment file containing paths to
data sources. Defaults to \code{.env}.}

\item{StartYear}{Numeric. The starting year for the occurrence data. Only
records from this year onward will be processed. Default is \code{1981}, which
matches the year ranges of CHELSA current climate data.}
}
\value{
Returns \code{NULL} invisibly after saving the processed data.
}
\description{
This function processes pre-cleaned and pre-standardized Integrated European
Long-Term Ecosystem, critical zone and socio-ecological Research
(\href{https://elter-ri.eu/}{eLTER}) data.
}
\note{
This function processes pre-cleaned vascular plants data from eLTER
sites, harmonized by Ahmed El-Gabbas. The original eLTER biodiversity data
were highly heterogeneous in format and structure, requiring
standardization and cleaning before use. Taxonomic standardization with the
GBIF backbone was performed by Marina Golivets (Feb. 2024).
}
\author{
Ahmed El-Gabbas
}
