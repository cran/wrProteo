#' Read Tabulated Files Exported By ProteomeDiscoverer At Protein Level, Deprecated
#'
#' Deprecated old version of Protein identification and quantification results from
#' \href{https://www.thermofisher.com/order/catalog/product/OPTON-30812}{Thermo ProteomeDiscoverer}
#' which were exported as tabulated text can be imported and relevant information extracted.
#' The final output is a list containing 3 elements: \code{$annot}, \code{$raw} and optional \code{$quant},
#' or returns data.frame with entire content of file if \code{separateAnnot=FALSE}.
#' Please use readProteomeDiscovererFile() from the same package instead !
#'
#'
#' @param fileName (character) name of file to be read
#' @param path (character) path of file to be read
#' @param normalizeMeth (character) normalization method, defaults to \code{median}, for more details see \code{\link[wrMisc]{normalizeThis}})
#' @param sampleNames (character) custom column-names for quantification data (ProteomeDiscoverer does not automatically use file-names from spectra); this argument has priority over \code{suplAnnotFile}
#' @param read0asNA (logical) decide if initial quntifications at 0 should be transformed to NA
#' @param quantCol (character or integer) exact col-names, or if length=1 content of \code{quantCol} will be used as pattern to search among column-names for $quant using \code{grep}
#' @param contamCol (character or integer, length=1) which columns should be used for contaminants marked by ProteomeDiscoverer.
#'  If a column named \code{contamCol} is found, the data will be mater on filtered to remove all contaminants, set to \code{NULL} for keeping all contaminants
#' @param refLi (character or integer) custom specify which line of data is main species, if character (eg 'mainSpe'), the column 'SpecType' in $annot will be searched for exact match of the (single) term given
#' @param separateAnnot (logical) if \code{TRUE} output will be organized as list with \code{$annot}, \code{$abund} for initial/raw abundance values and \code{$quant} with final log2 (normalized) quantitations
#' @param annotCol (character) column names to be read/extracted for the annotation section (default  c("Accession","Description","Gene","Contaminant","Sum.PEP.Score","Coverage....","X..Peptides","X..PSMs","X..Unique.Peptides", "X..AAs","MW..kDa.") )
#' @param FDRCol (list) optional indication to search for protein FDR information
#' @param specPref (character or list) define characteristic text for recognizing (main) groups of species (1st for contaminants - will be marked as 'conta', 2nd for main species- marked as 'mainSpe',
#'  and optional following ones for supplemental tags/species - marked as 'species2','species3',...);
#'  if list and list-element has multiple values they will be used for exact matching of accessions (ie 2nd of argument \code{annotCol})
#' @param gr (character or factor) custom defined pattern of replicate association, will override final grouping of replicates from \code{sdrf} and/or \code{suplAnnotFile} (if provided)   \code{}
#' @param sdrf (character, list or data.frame) optional extraction and adding of experimenal meta-data: if character, this may be the ID at ProteomeExchange,
#'   the second element may give further indications for automatic organization of groups of replicates.
#'   Besides, the output from \code{readSdrf} or a list from \code{defineSamples} may be provided; if \code{gr} is provided, \code{gr} gets priority for grouping of replicates
#' @param suplAnnotFile (logical or character) optional reading of supplemental files produced by ProteomeDiscoverer; however, if \code{gr} is provided, \code{gr} gets priority for grouping of replicates;
#'  if \code{TRUE} defaults to file '*InputFiles.txt' (needed to match information of \code{sdrf}) which can be exported next to main quantitation results;
#'  if \code{character} the respective file-name (relative or absolute path)
#' @param groupPref (list) additional parameters for interpreting meta-data to identify structure of groups (replicates), will be passed to \code{readSampleMetaData}.
#'   May contain \code{lowNumberOfGroups=FALSE} for automatically choosing a rather elevated number of groups if possible (defaults to low number of groups, ie higher number of samples per group)
#' @param plotGraph (logical or integer) optional plot of type vioplot of initial and normalized data (using \code{normalizeMeth}); if integer, it will be passed to \code{layout} when plotting
#' @param titGraph (character) custom title to plot of distribution of quantitation values
#' @param wex (integer) relative expansion factor of the violin-plot (will be passed to \code{\link[wrGraph]{vioplotW}})
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a list with \code{$raw} (initial/raw abundance values), \code{$quant} with final normalized quantitations, \code{$annot}, \code{$counts} an array with number of peptides, \code{$quantNotes}
#'  and \code{$notes}; or if \code{separateAnnot=FALSE} the function returns a data.frame with annotation and quantitation only
#' @seealso \code{\link[utils]{read.table}}, \code{\link[wrMisc]{normalizeThis}}) , \code{\link{readMaxQuantFile}}, \code{\link{readProlineFile}}, \code{\link{readFragpipeFile}}
#' @examples
#' path1 <- system.file("extdata", package="wrProteo")
#' fiNa <- "tinyPD_allProteins.txt.gz"
#' ## Please use the function readProteinDiscovererFile(), as shown below (same syntax)
#' dataPD <- readProteomeDiscovererFile(file=fiNa, path=path1, suplAnnotFile=FALSE)
#' summary(dataPD$quant)
#'
#' @export
readProtDiscovFile <- function(fileName, path=NULL, normalizeMeth="median", sampleNames=NULL, read0asNA=TRUE, quantCol="^Abundances*", annotCol=NULL, contamCol="Contaminant",
  refLi=NULL, separateAnnot=TRUE, FDRCol=list(c("^Protein.FDR.Confidence","High"), c("^Found.in.Sample.","High")), gr=NULL, sdrf=NULL, suplAnnotFile=TRUE,
  groupPref=list(lowNumberOfGroups=TRUE), specPref=c(conta="CON_|LYSC_CHICK", mainSpecies="OS=Homo sapiens"), plotGraph=TRUE, wex=1.6, titGraph="Proteome Discoverer",
  silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## read ProteomeDiscoverer exported txt

  .Deprecated(new="readProtDiscovererFile", package="wrProteo", msg="The function readProtDiscovFile() has been deprecated and replaced by readProtDiscovererFile() from the same package
     \n +++ Synthax and structure of output remain the same ! \n")

}
  
