% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readMaxQuantFile.R
\name{readMaxQuantFile}
\alias{readMaxQuantFile}
\title{Read proteinGroups.txt files exported from MaxQuant}
\usage{
readMaxQuantFile(
  path,
  fileName = "proteinGroups.txt",
  normalizeMeth = "median",
  quantCol = "LFQ.intensity",
  contamCol = "Potential.contaminant",
  pepCountCol = c("Razor + unique peptides", "Unique peptides", "MS.MS.count"),
  uniqPepPat = NULL,
  refLi = NULL,
  extrColNames = c("Majority.protein.IDs", "Fasta.headers", "Number.of.proteins"),
  specPref = c(conta = "conta|CON_|LYSC_CHICK", mainSpecies = "OS=Homo sapiens"),
  remRev = TRUE,
  separateAnnot = TRUE,
  tit = NULL,
  wex = 1.6,
  plotGraph = TRUE,
  silent = FALSE,
  callFrom = NULL
)
}
\arguments{
\item{path}{(character) path of file to be read}

\item{fileName}{(character) name of file to be read (default 'proteinGroups.txt' as typically generated by MaxQuant in txt folder). Gz-compressed files can be read, too.}

\item{normalizeMeth}{(character) normalization method (for details see \code{\link[wrMisc]{normalizeThis}})}

\item{quantCol}{(character or integer) exact col-names, or if length=1 content of \code{quantCol} will be used as pattern to search among column-names for $quant using \code{grep}}

\item{contamCol}{(character or integer, length=1) which columns should be used for contaminants marked by ProteomeDiscoverer}

\item{pepCountCol}{(character) pattern to search among column-names for count data (1st entry for 'Razor + unique peptides', 2nd fro 'Unique peptides', 3rd for 'MS.MS.count' (PSM))}

\item{uniqPepPat}{(character, length=1) depreciated, please use \code{pepCountCol} instead}

\item{refLi}{(character or integer) custom specify which line of data is main species, if character (eg 'mainSpe'), the column 'SpecType' in $annot will be searched for exact match of the (single) term given}

\item{extrColNames}{(character) column names to be read (1: prefix for LFQ quantitation, default 'LFQ.intensity'; 2: column name for protein-IDs, default 'Majority.protein.IDs'; 3: column names of fasta-headers, default 'Fasta.headers', 4: column name for number of protein IDs matching, default 'Number.of.proteins')}

\item{specPref}{(character) prefix to identifiers allowing to separate i) recognize contamination database, ii) species of main identifications and iii) spike-in species}

\item{remRev}{(logical) option to remove all protein-identifications based on reverse-peptides}

\item{separateAnnot}{(logical) if \code{TRUE} output will be organized as list with \code{$annot}, \code{$abund} for initial/raw abundance values and \code{$quant} with final normalized quantitations}

\item{tit}{(character) custom title to plot}

\item{wex}{(numeric)  relative expansion factor of the violin in plot}

\item{plotGraph}{(logical) optional plot vioplot of initial and normalized data (using \code{normalizeMeth}); alternatively the argument may contain numeric details that will be passed to \code{layout} when plotting}

\item{silent}{(logical) suppress messages}

\item{callFrom}{(character) allow easier tracking of message produced}
}
\value{
list with  \code{$raw} (initial/raw abundance values), \code{$quant} with final normalized quantitations, \code{$annot} (columns ), \code{$counts} an array with 'PSM' and 'NoOfRazorPeptides', \code{$quantNotes} and \code{$notes}; or a data.frame with quantitation and annotation if \code{separateAnnot=FALSE}
}
\description{
Quantification results form \href{https://www.maxquant.org/}{MaxQuant} can be read using this function and relevant information extracted. 
Innput files compressed as .gz can be read as well. Besides protein abundance values (XIC) peptide counting information like number of unique razor-peptides or PSM values can be extracted, too.
The protein abundance values mat be normalized using multiple methods (median normalization is default), the determination of normalization values can be restricted to specific proteins
(normalization to bait protein(s), or to matrix in UPS1 spike-in experiments).
Besides, a graphical display of the distruibution of protein abundance values may be generated.
}
\details{
This function has been developed using MaxQuant versions 1.6.10.x to 1.6.17.x, the format of resulting file 'proteinGroups.txt' is typically well conserved.  
The final output is a list containing these elements: \code{$raw}, \code{$quant}, \code{$annot}, \code{$counts}, \code{$quantNotes}, \code{$notes}, or (if \code{separateAnnot=FALSE}) data.frame 
  with annotation- and main quantification-content.
}
\examples{
path1 <- system.file("extdata", package="wrProteo")
# Here we'll load a short/trimmed example file (thus not the MaxQuant default name) 
fiNa <- "proteinGroupsMaxQuant1.txt.gz"
specPr <- c(conta="conta|CON_|LYSC_CHICK", mainSpecies="YEAST",spike="HUMAN_UPS")
dataMQ <- readMaxQuantFile(path1, file=fiNa, specPref=specPr, tit="tiny MaxQuant")
summary(dataMQ$quant)
matrixNAinspect(dataMQ$quant, gr=gl(3,3)) 
}
\seealso{
\code{\link[utils]{read.table}}, \code{\link[wrMisc]{normalizeThis}}) , \code{\link{readProlineFile}}
}
