% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/countNoOfCommonPeptides.R
\name{countNoOfCommonPeptides}
\alias{countNoOfCommonPeptides}
\title{Compare in-silico digested proteomes for unique and shared peptides, counts per protein or as peptides
  
Compare in-silico digested proteomes for unique and shared peptides, counts per protein or as peptides.
The in-silico digestion may be performed separately using the package \href{https://bioconductor.org/packages/release/bioc/html/cleaver.html}{cleaver}.
Note: input must be list (or multiple names lists) of proteins with their respective peptides (eg by in-silico digestion).}
\usage{
countNoOfCommonPeptides(
  ...,
  prefix = c("Hs", "Sc", "Ec"),
  sep = "_",
  silent = FALSE,
  debug = FALSE,
  callFrom = NULL
)
}
\arguments{
\item{...}{(list) multiple lists of (ini-silico) digested proteins (typically protein ID as names) with their respectice peptides (AA sequence), one entry for each species}

\item{prefix}{(character) optional (species-) prefix for entries in '...', will be only considered if '...' has no names}

\item{sep}{(character) concatenation symbol}

\item{silent}{(logical) suppress messages}

\item{debug}{(logical) display additional messages for debugging}

\item{callFrom}{(character) allow easier tracking of message(s) produced}
}
\value{
This function returns a list with $byPep as list of logical matrixes for each peptide (as line) and unique/shared/etc for each species; $byProt as list of matrixes with count data per proten (as line) for each species; $tab with simple summary-type count data
}
\description{
Compare in-silico digested proteomes for unique and shared peptides, counts per protein or as peptides
  
Compare in-silico digested proteomes for unique and shared peptides, counts per protein or as peptides.
The in-silico digestion may be performed separately using the package \href{https://bioconductor.org/packages/release/bioc/html/cleaver.html}{cleaver}.
Note: input must be list (or multiple names lists) of proteins with their respective peptides (eg by in-silico digestion).
}
\examples{
## The example mimics a proteomics experiment where extracts form E coli and 
## Saccharomyces cerevisiae were mixed, thus not all peptdes may occur unique.  
(mi2 = countNoOfCommonPeptides(Ec=list(E1=letters[1:4],E2=letters[c(3:7)],
  E3=letters[c(4,8,13)],E4=letters[9]),Sc=list(S1=letters[c(2:3,6)], 
  S2=letters[10:13],S3=letters[c(5,6,11)],S4=letters[c(11)],S5="n")))
##  a .. uni E, b .. inteR, c .. inteR(+intra E), d .. intra E  (no4), e .. inteR, 
##  f .. inteR +intra E   (no6), g .. uni E, h .. uni E  no 8), i .. uni E, 
##  j .. uni S (no10), k .. intra S  (no11), l .. uni S (no12), m .. inteR  (no13)
lapply(mi2$byProt,head)
mi2$tab
}
\seealso{
\code{\link{readFasta2}} and/or \code{cleave-methods} in package \href{https://bioconductor.org/packages/release/bioc/html/cleaver.html}{cleaver}
}
