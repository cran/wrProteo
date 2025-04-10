% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readProteomeDiscovererFile.R
\name{.commonSpecies}
\alias{.commonSpecies}
\title{Get Matrix With UniProt Abbreviations For Selected Species As Well As Simple Names}
\usage{
.commonSpecies()
}
\value{
This function returns a 2-column matrix with species names
}
\description{
This (low-level) function allows accessing matrix with UniProt abbreviations for species frequently used in research.
This information may be used to harmonize species descriptions or extract species information out of protein-names.
}
\examples{
.commonSpecies()
}
\seealso{
used eg in \code{readProtDiscovererFile},  \code{\link{readMaxQuantFile}}, \code{\link{readProlineFile}}, \code{\link{readFragpipeFile}}
}
