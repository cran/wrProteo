% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readMaxQuantPeptides.R
\name{readMaxQuantPeptides}
\alias{readMaxQuantPeptides}
\title{Read Peptide Identification and Quantitation Data-Files (peptides.txt) Produced By MaxQuant}
\usage{
readMaxQuantPeptides(
  path,
  fileName = "peptides.txt",
  normalizeMeth = "median",
  quantCol = "Intensity",
  contamCol = "Potential.contaminant",
  pepCountCol = "Experiment",
  refLi = NULL,
  sampleNames = NULL,
  extrColNames = c("Sequence", "Proteins", "Leading.razor.protein", "Start.position",
    "End.position", "Mass", "Missed.cleavages", "Unique..Groups.", "Unique..Proteins.",
    "Charges"),
  specPref = c(conta = "conta|CON_|LYSC_CHICK", mainSpecies = "HUMAN"),
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
\item{path}{(character) path of file to be read}

\item{fileName}{(character) name of file to be read (default 'peptides.txt' as typically generated by MaxQuant in txt folder). Gz-compressed files can be read, too.}

\item{normalizeMeth}{(character) normalization method (for details see \code{\link[wrMisc]{normalizeThis}})}

\item{quantCol}{(character or integer) exact col-names, or if length=1 content of \code{quantCol} will be used as pattern to search among column-names for $quant using \code{grep}}

\item{contamCol}{(character or integer, length=1) which columns should be used for contaminants}

\item{pepCountCol}{(character) pattern to search among column-names for count data (defaults to 'Experiment')}

\item{refLi}{(character or integer) custom specify which line of data should be used for normalization, ie which line is main species; if character (eg 'mainSpe'), the column 'SpecType' in $annot will be searched for exact match of the (single) term given}

\item{sampleNames}{(character) custom column-names for quantification data; this argument has priority over \code{suplAnnotFile}}

\item{extrColNames}{(character) column names to be read (1st position: prefix for quantitation, default 'intensity';
2nd: column name for peptide-IDs, default )}

\item{specPref}{(character) prefix to identifiers allowing to separate i) recognize contamination database,
ii) species of main identifications and iii) spike-in species}

\item{remRev}{(logical) option to remove all peptide-identifications based on reverse-peptides}

\item{remConta}{(logical) option to remove all peptides identified as contaminants}

\item{separateAnnot}{(logical) if \code{TRUE} output will be organized as list with \code{$annot}, \code{$abund}
for initial/raw abundance values and \code{$quant} with final normalized quantitations}

\item{gr}{(character or factor) custom defined pattern of replicate association, will override final grouping of
replicates from \code{sdrf} and/or \code{suplAnnotFile} (if provided)   \code{}}

\item{sdrf}{(character, list or data.frame) optional extraction and adding of experimenal meta-data: if character, this may be the ID at ProteomeExchange,
the second & third elements may give futher indicatations for automatic organization of groups of replicates.
Besides, the output from \code{readSdrf} or a list from \code{defineSamples} may be provided;
if \code{gr} is provided, \code{gr} gets priority for grouping of replicates;
if \code{sdrfOrder=TRUE} the output will be put in order of sdrf}

\item{suplAnnotFile}{(logical or character) optional reading of supplemental files produced by MaxQuant; if \code{gr} is provided, it gets priority for grouping of replicates
if \code{TRUE} default to files 'summary.txt' (needed to match information of \code{sdrf}) and 'parameters.txt' which can be found in the same folder as the main quantitation results;
if \code{character} the respective file-names (relative ro absolute path), 1st is expected to correspond to 'summary.txt' (tabulated text, the samples as given to MaxQuant) and 2nd to 'parameters.txt' (tabulated text, all parameters given to MaxQuant)}

\item{groupPref}{(list) additional parameters for interpreting meta-data to identify structure of groups (replicates), will be passed to \code{readSampleMetaData}.
May contain \code{lowNumberOfGroups=FALSE} for automatically choosing a rather elevated number of groups if possible (defaults to low number of groups, ie higher number of samples per group)
May contain \code{chUnit} (logical or character) to be passed to \code{readSampleMetaData()} for (optional) adjustig of unit-prefixes in meta-data group labels, in case multiple different unit-prefixes 
are used (eg '100pMol' and '1nMol').}

\item{titGraph}{(character) custom title to plot}

\item{wex}{(numeric)  relative expansion factor of the violin in plot}

\item{plotGraph}{(logical) optional plot vioplot of initial and normalized data (using \code{normalizeMeth}); alternatively the argument may contain numeric details that will be passed to \code{layout} when plotting}

\item{silent}{(logical) suppress messages}

\item{debug}{(logical) additional messages for debugging}

\item{callFrom}{(character) allows easier tracking of messages produced}
}
\value{
This function returns a list with  \code{$raw} (initial/raw abundance values), \code{$quant} with final normalized quantitations, \code{$annot} (columns ), \code{$counts} an array with 'PSM' and 'NoOfRazorPeptides',
  \code{$quantNotes}, \code{$notes} and optional \code{setup} for meta-data from \code{sdrf}; or a data.frame with quantitation and annotation if \code{separateAnnot=FALSE}
}
\description{
Peptide level identification and quantification data produced by \href{https://www.maxquant.org/}{MaxQuant} can be read using
this function and relevant information extracted.
Input files compressed as .gz can be read as well.
The peptide abundance values (XIC), peptide counting information and sample-annotation (if available) can be extracted, too.
}
\details{
The peptide annotation data gets parsed to extract specific fields (ID, name, description, species ...).
Besides, a graphical display of the distribution of peptide abundance values may be generated before and after normalization.


\href{https://www.maxquant.org/}{MaxQuant} is proteomics quantification software provided by the MaxPlanck institute.
By default MaxQuant write the results of each run to the path \code{combined/txt}, from there (only) the files
 'peptides.txt' (main quantitation at peptide level), 'summary.txt' and 'parameters.txt' will be used for this function.

Meta-data describing the samples and experimental setup may be available from two sources :
a) The file \code{summary.txt} which gets produced by MaxQuant in the same folder as the main quantification data.
b) Furthermore, meta-data deposited as \code{sdrf} at Pride can be retreived (via the respective github page) when giving
the accession number in argument \code{sdrf}.
Then, the meta-data will be examined for determining groups of replicates and
the results thereof can be found in $sampleSetup$levels.
Alternatively, a dataframe formatted like sdrf-files (ie for each sample a separate line, see also function \code{readSdrf}) may be given.
In tricky cases it is also possible to precise the column-name to use for defining the groups of replicates or the method for automatically choosing
 the most suited column via the 2nd value of the argument \code{sdrf}, see also the function \code{defineSamples} (which gets used internally).
Please note, that sdrf is still experimental and only a small fraction of proteomics-data on Pride have been annotated accordingly.
If a valid sdrf is furnished, it's information has priority over the information extracted from the MaxQuant produced file summary.txt.

This function has been developed using MaxQuant versions 1.6.10.x to 2.0.x, the format of the resulting file 'peptides.txt'
is typically well conserved between versions.
The final output is a list containing these elements: \code{$raw}, \code{$quant}, \code{$annot}, \code{$counts}, \code{$sampleSetup},
\code{$quantNotes}, \code{$notes}, or (if \code{separateAnnot=FALSE}) data.frame
  with annotation- and main quantification-content. If \code{sdrf} information has been found, an add-tional list-element \code{setup}
will be added containg the entire meta-data as \code{setup$meta} and the suggested organization as \code{setup$lev}.
}
\examples{
# Here we'll load a short/trimmed example file (thus not the MaxQuant default name)
MQpepFi1 <- "peptides_tinyMQ.txt.gz"
path1 <- system.file("extdata", package="wrProteo")
specPref1 <- c(conta="conta|CON_|LYSC_CHICK", mainSpecies="YEAST", spec2="HUMAN")
dataMQpep <- readMaxQuantPeptides(path1, file=MQpepFi1, specPref=specPref1,
  tit="Tiny MaxQuant Peptides")
summary(dataMQpep$quant)
}
\seealso{
\code{\link[utils]{read.table}}, \code{\link[wrMisc]{normalizeThis}}), for reading protein level \code{\link{readMaxQuantFile}}, \code{\link{readProlineFile}}
}
