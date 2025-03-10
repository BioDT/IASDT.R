% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/General_KeepOnly.R
\name{KeepOnly}
\alias{KeepOnly}
\title{Keep only specified objects in the environment, removing all others.}
\usage{
KeepOnly(Obj, Verbose = TRUE)
}
\arguments{
\item{Obj}{Character vector. Names of the objects to be kept in the
environment.}

\item{Verbose}{Logical. Whether to print the names of kept and removed
variables. Default to \code{TRUE}.}
}
\value{
No return value, called for side effects.
}
\description{
This function selectively retains the objects specified in the \code{Obj}
parameter in the current environment, removing all other objects. It is
useful for memory management by clearing unnecessary objects from the
environment. The function also provides an option to print the names of the
kept and removed variables.
}
\examples{
A <- B <- C <- 15
ls()

KeepOnly("A")

ls()
rm(list = ls())


A <- B <- C <- 15
KeepOnly(c("A","B"))
ls()
}
\author{
Ahmed El-Gabbas
}
