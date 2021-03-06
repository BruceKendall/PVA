% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractVRs.R
\name{extractVRs}
\alias{extractVRs}
\title{Extract vital rates from time series of Comadre/Compadre matrices}
\usage{
extractVRs(cadre)
}
\arguments{
\item{cadre}{A Compadre or Comadre object, subset to contain two or more "Individual" matrices for a single population}
}
\value{
A list with two elements
\describe{
\item{vrmat}{A matrix with a row for each vital rate and a column for each matrix from cadre}
\item{matstruct}{A list with lots of structural information about the vital rate decomposition of the matrix, which can be used by \code{\link{vrtomat}} to reconstruct the matrix.}
}
}
\description{
Extract vital rates from time series of Comadre/Compadre matrices
}
