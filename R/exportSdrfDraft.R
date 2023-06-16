#' Export Sample Meta-data from Quantification-Software as Sdrf-draft
#'
#' Sample/experimental annotation meta-data form \href{https://www.maxquant.org/}{MaxQuant} that was previously import can now be formatted in sdrf-style and exported
#' using this function to write a draft-sdrf-file.
#' Sdrf-files provide additional meta-information about samles and MS-runs in a standardized format, they may also be part of submissions to \href{https://www.ebi.ac.uk/pride/}{Pride}.
#'
#'
#' @details
#' Gathering information about samples and MS-runs requires that at data-import the additioal files created from software, like MaxQuant, was present and imported.
#' After exporting the draft sdrf the user is advised to check and complete the information in the resulting file.
#' Unfortunately, not all information present in a standard sdrf-file (like on \href{https://www.ebi.ac.uk/pride/}{Pride}) cannot be gathered automatically,
#' but key columns are already present and thus may facilitate completing.
#' Please note, that the file-format has been defined as \code{.tsv}, thus columns/fields shoud be separated by tabs.
#' At manual editing and completion, some editing- or tabulator-software may change the file-extesion to \code{.tsv.txt},
#' in this case the final files should be renamed as \code{.tsv} to remain compatible with Pride.
#'
#' At this point only the import of data from MaxQuant via \code{\link{readMaxQuantFile}} has been developed to extract information for creating a draft-sdrf.
#' Other data/file-import functions will be further developed to gather equivalent information in the future.
#'
#' @param lst (list) object created by import-function (MaxQuant)
#' @param fileName (character) file-name (and path) to be used when exprting
#' @param correctFileExtension (logical) if \code{TRUE} the fileName will get a \code{.tsv}-extension if not already present
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function writes an Sdrf draft to file
#'
#' @seealso This function may be used after reading/importig data by \code{\link{readMaxQuantFile}} in absence of sdrf
#' @examples
#' path1 <- system.file("extdata", package="wrProteo")
#' fiNaMQ <- "proteinGroups.txt.gz"
#' dataMQ <- readMaxQuantFile(path1, file=fiNaMQ, refLi="mainSpe", sdrf=FALSE, suplAnnotFile=TRUE)
#' ## Here we'll write simply in the current temporary directory of this R-session
#' exportSdrfDraft(dataMQ, file.path(tempdir(),"testSdrf.tsv"))
#'
#' @export
exportSdrfDraft <- function(lst, fileName="sdrfDraft.tsv", correctFileExtension=TRUE, silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## export sdrf-like
  ## Note : lst$sampleSetup$sdrfExport  is created only if import-function was run with sdrf=FALSE !!
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="exportSdrfDraft")
  if(isTRUE(debug)) silent <- FALSE
  if(!isTRUE(silent)) silent <- FALSE
  datOK <- TRUE
  msg <- "Invalid entry - nothing to do.  Argument 'lst' should be object like output of import-function from wrProteo"
  if(!is.list(lst) || length(lst) <1) { datOK <- FALSE
    if(!silent) message(fxNa, msg) }
  if(datOK) {
    expSdrf <- lst$sampleSetup$sdrfExport
    datDim <- dim(lst$quant)
    datOK <- length(expSdrf) >0 && isTRUE(datDim[2] >0)}
  msg <- "Nothing to do; invalid entry or no data for exporting available !!"
  if(!datOK) { if(!silent) message(fxNa, msg)
  } else {
    ## prepare for writing file
    fileName <- if(length(fileName) <1) "sdrfDraft.tsv" else fileName[1]
    if(is.na(fileName)) fileName <- "sdrfDraft.tsv"
    plCompl <- "please complete"
    na <- "not available"
    organisms <- table(lst$annot[which(lst$annot[,"Contam"] != "TRUE"), "Species"])
    organisms <- names(sort(organisms, decreasing=TRUE))
    organisms <- matrix(rep(organisms, each=datDim[2]), nrow=datDim[2], dimnames=list(NULL, rep("char_organism", length(organisms))))
    modMatr <- c(if(length(lst$sampleSetup$summaryD$Variable.modifications) > 0) lst$sampleSetup$summaryD$Variable.modifications[1],
      if(length(lst$sampleSetup$summaryD$Fixed.modifications) >0) lst$sampleSetup$summaryD$Fixed.modifications[1] )
    if(length(wrMisc::naOmit(modMatr)) >0) {
      modMatr <- unlist(lapply(modMatr, function(x) strsplit(as.character(x), ";") ))
      modMatr <- wrMisc::naOmit(unique(modMatr))
      modMatr <- matrix(rep(modMatr, each=datDim[2]), nrow=datDim[2], dimnames=list(NULL, rep("com_modification_parameters", length(modMatr))))
    } else modMatr <- NULL
    expo <- cbind(source.name=lst$sampleSetup$summaryD$ref, organisms, char_organism_part=na, char_disease=na, char_cell_type=na,
      char_mass=plCompl, assay.name=paste0("run",1:(datDim[2])), char_spiked_compound=plCompl,
      com_label=expSdrf["label"], com_instrument=plCompl, modMatr)
    expo <- cbind( expo, com_precursor_mass_tolerance=expSdrf["precMassTol"], com_fragment_mass_tolerance=expSdrf["fragMassTol"],
      com_technical_replicate=plCompl, com_fraction_identifyer=plCompl,
      com_file_uri="not available", com_data_file=lst$sampleSetup$summaryD[,1] , material_type=plCompl )
    ## chat_ for   characteristics[..]; com_ for comment[..]
    colnames(expo) <- gsub("_"," ", sub("com_","comment\\[", sub("char_","characteristics\\[", colnames(expo))))
    colnames(expo)[grep("\\[", colnames(expo))] <- paste0(colnames(expo)[grep("\\[", colnames(expo))],"]")   # complete brackets
    ## main write/export
    chWr <- try(utils::write.table(expo, fileName, quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE))
    if(inherits(chWr, "try-error")) {  if(!silent) message(fxNa,"FAILED to write file  (check if rights to write)") 
    } else if(!silent) message(fxNa,"Successfully exported sdrf-draft to file '",fileName,"'") 
  }

}

