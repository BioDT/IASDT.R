% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/General_AssignEnvVars.R
\name{AssignEnvVars}
\alias{AssignEnvVars}
\title{Assign environment variables from a .env file}
\usage{
AssignEnvVars(EnvFile = ".env", EnvVarDT = NULL)
}
\arguments{
\item{EnvFile}{Character. Path to the environment file containing paths to
data sources. Defaults to \code{.env}.}

\item{EnvVarDT}{\code{data.frame}. A data frame or tibble containing the columns
\code{VarName}, \code{Value}, \code{CheckDir}, and \code{CheckFile}. Each row specifies an
environment variable, the name to assign it to, and whether to check if it
is a directory or file. This structure allows for additional validation
checks on the variables being imported.}
}
\value{
This function is used for its side effects of setting environment
variables in the R environment. It assigns each variable from the \code{.env}
file to the R environment with the name specified in the \code{EnvVarDT} data
frame.
}
\description{
This function reads environment variables from a specified \code{.env} file and
assigns them to variables in the R environment based on a data frame
containing variable names, values, and checks for directories and files. It
is designed to facilitate the management of environment variables in a
structured and reproducible manner.
}
\author{
Ahmed El-Gabbas
}
