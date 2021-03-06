% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readProtDiscovFile.R
\name{readProtDiscovFile}
\alias{readProtDiscovFile}
\title{Read tabulated files exported by ProteomeDiscoverer}
\usage{
readProtDiscovFile(
  fileName,
  path = NULL,
  normalizeMeth = "median",
  sampleNames = NULL,
  infoFile = TRUE,
  read0asNA = TRUE,
  quantCol = "^Abundances*",
  annotCol = NULL,
  contamCol = "Contaminant",
  refLi = NULL,
  separateAnnot = TRUE,
  FDRCol = list(c("^Protein.FDR.Confidence", "High"), c("^Found.in.Sample.", "High")),
  plotGraph = TRUE,
  tit = "Proteome Discoverer",
  graphTit = NULL,
  wex = 1.6,
  specPref = c(conta = "CON_|LYSC_CHICK", mainSpecies = "OS=Homo sapiens"),
  silent = FALSE,
  debug = FALSE,
  callFrom = NULL
)
}
\arguments{
\item{fileName}{(character) name of file to be read}

\item{path}{(character) path of file to be read}

\item{normalizeMeth}{(character) normalization method (will be sent to  \code{\link[wrMisc]{normalizeThis}})}

\item{sampleNames}{(character) new column-names for quantification data (ProteomeDiscoverer does not automatically use file-names from spectra); this argument has priority over \code{infoFile}}

\item{infoFile}{(character or logical) filename containing additional information about MS-samples (produced by ProteomeDiscoverer default '*.InputFiles.txt'), 
if \code{TRUE} the first file in \code{path} containing the default name will be used. If no specific \code{sampleNames} given, the filenames will will be trimmed to remove redundant text and used as \code{sampleNames}. 
Besides, ProteomeDiscoverer version number and full raw-file path will be extracted for $notes in fial output.}

\item{read0asNA}{(logical) decide if initial quntifications at 0 should be transformed to NA}

\item{quantCol}{(character or integer) exact col-names, or if length=1 content of \code{quantCol} will be used as pattern to search among column-names for $quant using \code{grep}}

\item{annotCol}{(character) column names to be read/extracted for the annotation section (default  c("Accession","Description","Gene","Contaminant","Sum.PEP.Score","Coverage....","X..Peptides","X..PSMs","X..Unique.Peptides", "X..AAs","MW..kDa.") )}

\item{contamCol}{(character or integer, length=1) which columns should be used for contaminants marked by ProteomeDiscoverer. 
If a column named \code{contamCol} is found, the data will be lateron filtered to remove all contaminants, set to \code{NULL} for keeping all contaminants}

\item{refLi}{(character or integer) custom specify which line of data is main species, if character (eg 'mainSpe'), the column 'SpecType' in $annot will be searched for exact match of the (single) term given}

\item{separateAnnot}{(logical) if \code{TRUE} output will be organized as list with \code{$annot}, \code{$abund} for initial/raw abundance values and \code{$quant} with final normalized quantitations}

\item{FDRCol}{(list) optional indication to search for protein FDR information}

\item{plotGraph}{(logical) optional plot of type vioplot of initial and normalized data (using \code{normalizeMeth}); if integer, it will be passed to \code{layout} when plotting}

\item{tit}{(character) custom title to plot}

\item{graphTit}{(character) depreciated custom title to plot, please use 'tit'}

\item{wex}{(integer) relative expansion factor of the violin-plot (will be passed to \code{\link[wrGraph]{vioplotW}})}

\item{specPref}{(character or list) define characteristic text for recognizing (main) groups of species (1st for comtaminants - will be marked as 'conta', 2nd for main species- marked as 'mainSpe', 
and optional following ones for supplemental tags/species - maked as 'species2','species3',...); 
if list and list-element has multiple values they will be used for exact matching of accessions (ie 2nd of argument \code{annotCol})}

\item{silent}{(logical) suppress messages}

\item{debug}{(logical) additional messages for debugging}

\item{callFrom}{(character) allow easier tracking of messages produced}
}
\value{
list with \code{$raw} (initial/raw abundance values), \code{$quant} with final normalized quantitations, \code{$annot}, \code{$counts} an array with number of peptides, \code{$quantNotes} 
 and \code{$notes}; or if \code{separateAnnot=FALSE} the function returns a data.frame with annotation and quantitation only
}
\description{
Protein quantification results form \href{https://www.thermofisher.com/order/catalog/product/OPTON-30812}{Thermo ProteomeDiscoverer} 
which were exported as tabulated text can be imported and relevant information extracted. 
The final output is a list containing 3 elements: \code{$annot}, \code{$raw} and optional \code{$quant}, or returns data.frame with entire content of file if \code{separateAnnot=FALSE}.
}
\details{
This function has been developed using Thermo ProteomeDiscoverer versions 2.2 to 2.5.
The format of resulting files at export also depends which columns are chosen as visible inside ProteomeDiscoverer and subsequently get chosen for export.
Using the argument \code{infoFile} it is possible to specify a specific file (or search for default file) to read for extracting file-names as sample-names and other experiment realted information.  
If a column named \code{contamCol} is found, the data will be lateron filtered to remove all contaminants, set to \code{NULL} for keeping all contaminants 
This function replaces the depreciated function \code{readPDExport}.
}
\examples{
path1 <- system.file("extdata", package="wrProteo")
fiNa <- "tinyPD_allProteins.txt.gz"
dataPD <- readProtDiscovFile(file=fiNa, path=path1)
summary(dataPD$quant)

}
\seealso{
\code{\link[utils]{read.table}}, \code{\link[wrMisc]{normalizeThis}}) , \code{\link{readMaxQuantFile}}, \code{\link{readProlineFile}}
}
