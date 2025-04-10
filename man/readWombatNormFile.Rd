% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readWombatNormFile.R
\name{readWombatNormFile}
\alias{readWombatNormFile}
\title{Read (Normalized) Quantitation Data Files Produced By Wombat At Protein Level}
\usage{
readWombatNormFile(
  fileName,
  path = NULL,
  quantSoft = "(quant software not specified)",
  fasta = NULL,
  isLog2 = TRUE,
  normalizeMeth = "none",
  quantCol = "abundance_",
  contamCol = NULL,
  pepCountCol = c("number_of_peptides"),
  read0asNA = TRUE,
  refLi = NULL,
  sampleNames = NULL,
  extrColNames = c("protein_group"),
  specPref = NULL,
  remRev = TRUE,
  remConta = FALSE,
  separateAnnot = TRUE,
  gr = NULL,
  sdrf = NULL,
  suplAnnotFile = NULL,
  groupPref = list(lowNumberOfGroups = TRUE, chUnit = TRUE),
  titGraph = NULL,
  wex = 1.6,
  plotGraph = TRUE,
  silent = FALSE,
  debug = FALSE,
  callFrom = NULL
)
}
\arguments{
\item{fileName}{(character) name of file to be read (default 'proteinGroups.txt' as typically generated by Compomics in txt folder). Gz-compressed files can be read, too.}

\item{path}{(character) path of file to be read}

\item{quantSoft}{(character) qunatification-software used inside Wombat-P}

\item{fasta}{(logical or character) if \code{TRUE} the (first) fasta from one direcory higher than \code{fileName} will be read as fasta-file to extract further protein annotation;
if \code{character} a fasta-file at this location will be read/used/}

\item{isLog2}{(logical) typically data read from Wombat are expected to be \code{isLog2=TRUE}}

\item{normalizeMeth}{(character) normalization method, defaults to \code{median}, for more details see \code{\link[wrMisc]{normalizeThis}})}

\item{quantCol}{(character or integer) exact col-names, or if length=1 content of \code{quantCol} will be used as pattern to search among column-names for $quant using \code{grep}}

\item{contamCol}{(character or integer, length=1) which columns should be used for contaminants}

\item{pepCountCol}{(character) pattern to search among column-names for count data (1st entry for 'Razor + unique peptides', 2nd fro 'Unique peptides', 3rd for 'MS.MS.count' (PSM))}

\item{read0asNA}{(logical) decide if initial quntifications at 0 should be transformed to NA (thus avoid -Inf in log2 results)}

\item{refLi}{(character or integer) custom specify which line of data should be used for normalization, ie which line is main species; if character (eg 'mainSpe'), the column 'SpecType' in $annot will be searched for exact match of the (single) term given}

\item{sampleNames}{(character) custom column-names for quantification data; this argument has priority over \code{suplAnnotFile}}

\item{extrColNames}{(character) column names to be read (1st position: prefix for LFQ quantitation, default 'LFQ.intensity'; 2nd: column name for protein-IDs, default 'Majority.protein.IDs'; 3rd: column names of fasta-headers, default 'Fasta.headers', 4th: column name for number of protein IDs matching, default 'Number.of.proteins')}

\item{specPref}{(character) prefix to identifiers allowing to separate i) recognize contamination database, ii) species of main identifications and iii) spike-in species}

\item{remRev}{(logical) option to remove all protein-identifications based on reverse-peptides}

\item{remConta}{(logical) option to remove all proteins identified as contaminants}

\item{separateAnnot}{(logical) if \code{TRUE} output will be organized as list with \code{$annot}, \code{$abund} for initial/raw abundance values and \code{$quant} with final normalized quantitations}

\item{gr}{(character or factor) custom defined pattern of replicate association, will override final grouping of replicates from \code{sdrf} and/or \code{suplAnnotFile} (if provided)   \code{}}

\item{sdrf}{(logical, character, list or data.frame) optional extraction and adding of experimenal meta-data:
if \code{sdrf=TRUE} the 1st sdrf in the directory above \code{fileName} will be used
if character, this may be the ID at ProteomeExchange,
the second element may give futher indicatations for automatic organization of groups of replicates.
Besides, the output from \code{readSdrf} or a list from \code{defineSamples} may be provided; if \code{gr} is provided, \code{gr} gets priority for grouping of replicates}

\item{suplAnnotFile}{(logical or character) optional reading of supplemental files produced by Compomics; if \code{gr} is provided, it gets priority for grouping of replicates
if \code{TRUE} default to files 'summary.txt' (needed to match information of \code{sdrf}) and 'parameters.txt' which can be found in the same folder as the main quantitation results;
if \code{character} the respective file-names (relative ro absolute path), 1st is expected to correspond to 'summary.txt' (tabulated text, the samples as given to Compomics) and 2nd to 'parameters.txt' (tabulated text, all parameters given to Compomics)}

\item{groupPref}{(list) additional parameters for interpreting meta-data to identify structure of groups (replicates), will be passed to \code{readSampleMetaData}.
May contain \code{lowNumberOfGroups=FALSE} for automatically choosing a rather elevated number of groups if possible (defaults to low number of groups, ie higher number of samples per group)
May contain \code{chUnit} (logical or character) to be passed to \code{readSampleMetaData()} for (optional) adjustig of unit-prefixes in meta-data group labels, in case multiple different unit-prefixes 
are used (eg '100pMol' and '1nMol').}

\item{titGraph}{(character) custom title to plot of distribution of quantitation values}

\item{wex}{(numeric)  relative expansion factor of the violin in plot}

\item{plotGraph}{(logical) optional plot vioplot of initial and normalized data (using \code{normalizeMeth}); alternatively the argument may contain numeric details that will be passed to \code{layout} when plotting}

\item{silent}{(logical) suppress messages}

\item{debug}{(logical) additional messages for debugging}

\item{callFrom}{(character) allow easier tracking of messages produced}
}
\value{
This function returns a list with  \code{$raw} (initial/raw abundance values), \code{$quant} with final normalized quantitations, \code{$annot} (columns ), \code{$counts} an array with 'PSM' and 'NoOfRazorPeptides',
  \code{$quantNotes}, \code{$notes} and optional \code{setup} for meta-data from \code{sdrf}; or a data.frame with quantitation and annotation if \code{separateAnnot=FALSE}
}
\description{
Protein quantification results from \href{https://github.com/wombat-p}{Wombat-P} using the Bioconductor package Normalizer can be read using this function and relevant information extracted.
Input files compressed as .gz can be read as well.
The protein abundance values (XIC), peptide counting get extracted. Since protein annotation is not very extensive with this format of data, the function allows reading the
initial fasta files (from the directory above the quantitation-results) allowing to extract more protein-annotation (like species).
Sample-annotation (if available) can be extracted from  sdrf files, which are typically part of the Wombat output, too.
The protein abundance values may be normalized using multiple methods (median normalization as default), the determination of normalization factors can be restricted to specific proteins
(normalization to bait protein(s), or to invariable matrix of spike-in experiments).
The protein annotation data gets parsed to extract specific fields (ID, name, description, species ...).
Besides, a graphical display of the distribution of protein abundance values may be generated before and after normalization.
}
\details{
By standard workflow of Wombat-P writes the results of each analysis-method/quantification-algorithm as .csv files
Meta-data describing the proteins may be available from two sources :
a) The 1st column of the Wombat/normalizer output.
b) Form the .fasta file in the directory above the analysis/quantiication results of the Wombar-workflow


Meta-data describing the samples and experimental setup may be available from a sdrf-file (from the directory above the analysis/quantiication results)
If available, the meta-data will be examined for determining groups of replicates and
the results thereof can be found in $sampleSetup$levels.
Alternatively, a dataframe formatted like sdrf-files (ie for each sample a separate line, see also function \code{readSdrf}) may be given, too.

This import-function has been developed using Wombat-P version 1.x.
The final output is a list containing these elements: \code{$raw}, \code{$quant}, \code{$annot}, \code{$counts}, \code{$sampleSetup}, \code{$quantNotes}, \code{$notes}, or (if \code{separateAnnot=FALSE}) data.frame
  with annotation- and main quantification-content. If \code{sdrf} information has been found, an add-tional list-element \code{setup}
will be added containg the entire meta-data as \code{setup$meta} and the suggested organization as \code{setup$lev}.
}
\examples{
path1 <- system.file("extdata", package="wrProteo")
# Here we'll load a short/trimmed example file (originating from Compomics)
fiNa <- "tinyWombCompo1.csv.gz"
dataWB <- readWombatNormFile(file=fiNa, path=path1, tit="tiny Wombat/Compomics, Normalized ")
summary(dataWB$quant)
}
\seealso{
\code{\link[utils]{read.table}}, \code{\link[wrMisc]{normalizeThis}}) , \code{\link{readProteomeDiscovererFile}}; \code{\link{readProlineFile}} (and other import-functions), \code{\link{matrixNAinspect}}
}
