% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exportAsWombatP.R
\name{exportAsWombatP}
\alias{exportAsWombatP}
\title{Export As Wombat-P Set Of Files}
\usage{
exportAsWombatP(
  wrProtObj,
  path = ".",
  combineFractions = "mean",
  silent = FALSE,
  debug = FALSE,
  callFrom = NULL
)
}
\arguments{
\item{wrProtObj}{(list produced by any import-function from wrProteo) object which will be exported as Wombat-P format}

\item{path}{(character) the location where the data should be exorted to}

\item{combineFractions}{(\code{NULL} or character (length=1)) if not \code{NULL} this assigns the method how multiple farctions should be combined
(at this point only the method 'mean' is implemented)}

\item{silent}{(logical) suppress messages}

\item{debug}{(logical) display additional messages for debugging}

\item{callFrom}{(character) allows easier tracking of messages produced}
}
\value{
This function creates a set of files (\code{README.md}, \code{test_params.yml}), plus a sud-directory containig file(s) (\code{stand_prot_quant_method.csv}); finally the function returns  (\code{NULL}),
}
\description{
This function allows exporting objects created from wrProteo to the format of Wombat-P \href{https://github.com/wombat-p}{Wombat-P}.
}
\examples{
path1 <- system.file("extdata", package="wrProteo")
fiNa <- "proteinGroupsMaxQuant1.txt.gz"
specPr <- c(conta="conta|CON_|LYSC_CHICK", mainSpecies="YEAST", spike="HUMAN_UPS")
dataMQ <- readMaxQuantFile(path1, file=fiNa, specPref=specPr, tit="tiny MaxQuant")

exportAsWombatP(dataMQ, path=tempdir())
}
\seealso{
\code{\link{readMaxQuantFile}}, \code{\link{readProteomeDiscovererFile}}; \code{\link[wrMisc]{moderTestXgrp}} or \code{\link[wrMisc]{moderTest2grp}}
}
