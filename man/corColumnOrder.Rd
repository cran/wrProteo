% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corColumnOrder.R
\name{corColumnOrder}
\alias{corColumnOrder}
\title{Order columns in list of matrixes (or matrix)}
\usage{
corColumnOrder(
  dat,
  sampNames,
  useListElem = c("quant", "raw"),
  silent = FALSE,
  callFrom = NULL
)
}
\arguments{
\item{dat}{(matrix, list or MArrayLM-object from limma) main input of which columns should get re-ordered, may be output from \code{\link{moderTestXgrp}} or \code{\link{moderTest2grp}}.}

\item{sampNames}{(character) column-names in desired order for output}

\item{useListElem}{(character) in case \code{dat} is list, all list-elements who's columns should get (re-)ordered}

\item{silent}{(logical) suppress messages}

\item{callFrom}{(character) allow easier tracking of message(s) produced}
}
\value{
object of same class as input \code{dat}
}
\description{
This function orders columns in list of matrixes (or matrix) according to argument \code{sampNames}. 
This function can be used to adjut/correct the order of samples after reading data using \code{readMaxQuantFile()}, \code{readPDExport()} etc. 
The input may also be MArrayLM-type object from package \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{limma} or from \code{\link{moderTestXgrp}} or \code{\link{moderTest2grp}}.
}
\examples{
grp <- factor(rep(LETTERS[c(3,1,4)], c(2,3,3)))
dat1 <- matrix(1:15, ncol=5, dimnames=list(NULL,c("D","A","C","E","B")))
corColumnOrder(dat1,sampNames=LETTERS[1:5])  

dat1 <- list(quant=dat1,raw=dat1)
  dat1
corColumnOrder(dat1,sampNames=LETTERS[1:5]) 
}
\seealso{
\code{\link{moderTestXgrp}} for single comparisons, \code{\link[base]{order}}
}
