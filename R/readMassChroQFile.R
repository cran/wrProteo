#' Read tabulated files imported from MassChroQ
#'
#' Quantification results using MassChroQ should be initially treated using the R-package MassChroqR (both distributed by the PAPPSO at http://pappso.inrae.fr/)
#'  for initial normalization on peptide-level and combination of peptide values into protein abundances.
#'
#' The final output of this fucntion is a list containing 3 elements: \code{$annot}, \code{$raw}, \code{$quant} and  \code{$notes}, or returns data.frame with entire content of file if \code{separateAnnot=FALSE}. Other list-elements remain empty to keep format compatible to other import functions.
#'
#' @details
#' This function has been developed using MassChroQ version 2.2 and R-package MassChroqR version 0.4.0. Both are distributed by the PAPPSO (http://pappso.inrae.fr/).
#' When saving quantifications generated in R as RData (with extension .rdata or .rda) using the R-packages associated with MassChroq, the ABUNDANCE_TABLE produced by  mcq.get.compar(XICAB) should be used.
#'
#' After import data get (re-)normalized according to \code{normalizeMeth} and \code{refLi}, and boxplots or vioplots drawn.
#'
#'
#' @param fileName (character) name of file to be read (may be tsv, csv, rda or rdata); both US and European csv formats are supported
#' @param path (character) path of file to be read
#' @param normalizeMeth (character) normalization method (will be sent to  \code{\link[wrMisc]{normalizeThis}})
#' @param sampleNames (character) custom column-names for quantification data; this argument has priority over \code{suplAnnotFile}
#' @param refLi (character or integer) custom specify which line of data is main species, if character (eg 'mainSpe'), the column 'SpecType' in $annot will be searched for exact match of the (single) term given
#' @param separateAnnot (logical) if \code{TRUE} output will be organized as list with \code{$annot}, \code{$abund} for initial/raw abundance values and \code{$quant} with final normalized quantitations
#' @param titGraph (character) custom title to plot of distribution of quantitation values
#' @param wex (integer) relative expansion factor of the violin-plot (will be passed to \code{\link[wrGraph]{vioplotW}})
#' @param specPref (character or list) define characteristic text for recognizing (main) groups of species (1st for comtaminants - will be marked as 'conta', 2nd for main species- marked as 'mainSpe',
#'  and optional following ones for supplemental tags/species - maked as 'species2','species3',...);
#'  if list and list-element has multiple values they will be used for exact matching of accessions (ie 2nd of argument \code{annotCol})
#' @param separateAnnot (logical) if \code{TRUE} output will be organized as list with \code{$annot}, \code{$abund} for initial/raw abundance values and \code{$quant} with final normalized quantitations
#' @param gr (character or factor) custom defined pattern of replicate association, will override final grouping of replicates from \code{sdrf} and/or \code{suplAnnotFile} (if provided)   \code{}
#' @param sdrf (character, list or data.frame) optional extraction and adding of experimenal meta-data: if character, this may be the ID at ProteomeExchange,
#'   the second & third elements may give futher indicatations for automatic organization of groups of replicates.
#'   Besides, the output from \code{readSdrf} or a list from \code{defineSamples} may be provided;
#'   if \code{gr} is provided, \code{gr} gets priority for grouping of replicates;
#'   if \code{sdrfOrder=TRUE} the output will be put in order of sdrf
#' @param suplAnnotFile (logical or character) optional reading of supplemental files produced by ProteomeDiscoverer; however, if \code{gr} is provided, \code{gr} gets priority for grouping of replicates;
#'  if \code{TRUE} defaults to file '*InputFiles.txt' (needed to match information of \code{sdrf}) which can be exported next to main quantitation results;
#'  if \code{character} the respective file-name (relative or absolute path)
#' @param groupPref (list) additional parameters for interpreting meta-data to identify structure of groups (replicates), will be passed to \code{readSampleMetaData}.
#'   May contain \code{lowNumberOfGroups=FALSE} for automatically choosing a rather elevated number of groups if possible (defaults to low number of groups, ie higher number of samples per group)
#'   May contain \code{chUnit} (logical or character) to be passed to \code{readSampleMetaData()} for (optional) adjustig of unit-prefixes in meta-data group labels, in case multiple different unit-prefixes 
#'   are used (eg '100pMol' and '1nMol').
#' @param plotGraph (logical) optional plot of type vioplot of initial and normalized data (using \code{normalizeMeth}); if integer, it will be passed to \code{layout} when plotting
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns list with \code{$raw} (initial/raw abundance values), \code{$quant} with final normalized quantitations, \code{$annot}, \code{$counts} an array with number of peptides, \code{$quantNotes} and \code{$notes}; or if \code{separateAnnot=FALSE} the function returns a data.frame with annotation and quantitation only
#' @seealso \code{\link[utils]{read.table}}, \code{\link[wrMisc]{normalizeThis}}) , \code{\link{readProlineFile}}
#' @examples
#' path1 <- system.file("extdata", package="wrProteo")
#' fiNa <- "tinyMC.RData"
#' dataMC <- readMassChroQFile(file=fiNa, path=path1)
#' @export
readMassChroQFile <- function(fileName, path=NULL, normalizeMeth="median", sampleNames=NULL, refLi=NULL, separateAnnot=TRUE, titGraph="MassChroQ", wex=NULL,
  specPref=c(conta="CON_|LYSC_CHICK", mainSpecies="OS=Homo sapiens"), gr=NULL, sdrf=NULL, suplAnnotFile=FALSE,
  groupPref=list(lowNumberOfGroups=TRUE, chUnit=TRUE), plotGraph=TRUE, silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## read MassChroQ (pre-)treated data
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readMassChroQFile")

  oparMar <- graphics::par("mar")                     # old margins, for rest after figure
  oparLayout <- graphics::par("mfcol")                # old layout, for rest after figure
  if(plotGraph) on.exit(graphics::par(mar=oparMar, mfcol=oparLayout))            # restore old mar settings

  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  if(!requireNamespace("utils", quietly=TRUE)) stop("package 'utils' not found ! Please install first from CRAN")
  tmp <- list()                                                        # initialize
  infoDat <- infoFi <- setupSd <- parametersD <- NULL        # initialize for sdrf annotation
  counts <- NULL                                             # so far PSM data are not accessible

  ## check & read file
  paFi <- wrMisc::checkFilePath(fileName, path, expectExt=NULL, compressedOption=TRUE, stopIfNothing=TRUE, callFrom=fxNa, silent=silent,debug=debug)

  if(length(c(grep("\\.txt$",paFi), grep("\\.txt\\.gz$",paFi))) >0) tmp[[1]] <- try(utils::read.delim(paFi, stringsAsFactors=FALSE), silent=TRUE)             # read tabulated text-file
  if(length(c(grep("\\.csv$",paFi), grep("\\.csv\\.gz$",paFi))) >0) tmp[[2]] <- try(utils::read.csv(paFi, stringsAsFactors=FALSE), silent=TRUE)               # read US csv-file
  if(length(c(grep("\\.csv$",paFi), grep("\\.csv\\.gz$",paFi))) >0) tmp[[3]] <- try(utils::read.csv2(paFi, stringsAsFactors=FALSE), silent=TRUE)              # read Euro csv-file
  if(length(c(grep("\\.tsv$",paFi), grep("\\.tsv\\.gz$",paFi))) >0) tmp[[4]] <- try(utils::read.csv(file=paFi, stringsAsFactors=FALSE, sep='\t', header=TRUE)) # read US comma tsv-file
  if(length(c(grep("\\.rda$",paFi), grep("\\.rdata$",tolower(paFi)))) >0) {
    ls1 <- ls()
    tmp[[5]] <- try(load(paFi))
    if(!inherits(tmp[[5]], "try-error")) {          # dont know under which name the object was saved in RData..
      if(length(ls1) +2 ==length(ls())) {
        tmp[[5]] <- get(ls()[which(!ls() %in% ls1 & ls() != "ls1")])            # found no way of removing initial object
        if(!silent) message(fxNa,"Loading R-object '",ls()[which(!ls() %in% ls1 & ls() != "ls1")],"' as quantification data out of ",paFi)
      } else stop(" Either .RData is empty or element loaded has name of one of the arguments of this function and can't be recognized as such")
    } else stop("Failed to load .RData") }
  if(debug) {message(fxNa,"mc1")}

  if(length(tmp) <1) stop("Failed to recognize file extensions of input data (unknown format)")
  chCl <- sapply(tmp, inherits, "try-error")
  if(all(chCl)) stop(" Failed to extract data from '",fileName,"'  (check format & rights to read)")
  nCol <- sapply(tmp, function(x) if(length(x) >0) {if(!inherits(x, "try-error")) ncol(x) else NA} else NA)
  bestT <- which.max(nCol)
  fiType <- c("txt","UScsv","EURcsv","tsv","RData")[bestT]
  tmp <- tmp[[bestT]]
  ## tibble colnames may include/start with '#' ... adopt to rest
  corColNa <- grep("^#",colnames(tmp))
  if(length(corColNa) >0) colnames(tmp)[which(corColNa)] <- sub("^#","X.",colnames(tmp)[which(corColNa)])    # make colnames alike
  if(debug) {message(fxNa,"mc2")}
  

  ## recover OS
  tmp2 <- sub("^[[:alpha:]]+\\|","", rownames(tmp))               # trim heading database-name
  annot <- cbind(Accession=sub("\\|[[:upper:]]+[[:digit:]]*_{0,1}[[:upper:]][[:print:]]*","",tmp2),
    EntryName=sub("^[[:upper:]]+[[:digit:]]*\\|","",tmp2), GeneName=NA, Species=NA, Contam=NA, SpecType=NA)      #

  ## extract species out of EntryName
  commonSpec <- .commonSpecies()
  spec <- apply(commonSpec, 1, function(x) grep(paste0(x[1],"$"), tmp2))
  chLe <- sapply(spec,length) >0
  if(any(chLe)) for(i in which(chLe)) annot[spec[[i]],"Species"] <- commonSpec[i,2]
  ## SpecType
  if(length(specPref) >0) {
    for(i in 1:length(specPref)) { ch1 <- grep(specPref[i], rownames(tmp))
      if(length(ch1) >0) annot[which(ch1),"SpecType"] <- names(specPref)[i] } }
  ## checke for unique rownames
  ch1 <- duplicated(annot[,1])
  if(all(!ch1)) rownames(annot) <- annot[,1]
  if(debug) {message(fxNa,"mc3")}  

  ## colnames for quantitative data
  if(length(sampleNames) >0) if(length(sampleNames)==ncol(tmp)) colnames(tmp) <- sampleNames else message(fxNa,"invalid entry of 'sampleNames' (incorrect length)")
  
  ## check for reference for normalization
  refLiIni <- refLi
  if(is.character(refLi) && length(refLi)==1) { 
    refLi <- which(annot[,"SpecType"]==refLi)
    if(length(refLi) <1 ) { refLi <- 1:nrow(tmp)
      if(!silent) message(fxNa,"Could not find any proteins matching argument 'refLi=",refLiIni,"', ignoring ...")
    } else {
      if(!silent) message(fxNa,"Normalize using (custom) subset of ",length(refLi)," lines specified as '",refLiIni,"'")}}    # may be "mainSpe"

  ## normalize
  abund <- if(!is.matrix(tmp)) as.matrix(tmp) else tmp
  quant <- wrMisc::normalizeThis(abund, method=normalizeMeth, mode="additive", refLines=refLi, callFrom=fxNa)       #
  if(debug) {message(fxNa,"mc4")}

  ### GROUPING OF REPLICATES AND SAMPLE META-DATA
  if(length(suplAnnotFile) >0) {
    if(length(sampleNames) %in% c(1, ncol(abund))) groupPref$sampleNames <- sampleNames
    if(length(gr) %in% c(1, ncol(abund))) groupPref$gr <- gr
    setupSd <- readSampleMetaData(sdrf=sdrf, suplAnnotFile=suplAnnotFile, quantMeth="MC", path=path, abund=utils::head(quant), chUnit=isTRUE(groupPref$chUnit), groupPref=groupPref, silent=silent, debug=debug, callFrom=fxNa)
  }
  if(debug) {message(fxNa,"rmc13b .."); rmc13b <- list(sdrf=sdrf,gr=gr,suplAnnotFile=suplAnnotFile,abund=abund, refLi=refLi,annot=annot,setupSd=setupSd,sampleNames=sampleNames)}

  ## finish groups of replicates & annotation setupSd
    setupSd <- .checkSetupGroups(abund=abund, setupSd=setupSd, gr=gr, sampleNames=sampleNames, quantMeth="MC", silent=silent, debug=debug, callFrom=fxNa)

    ## option : set order of samples as sdrf
    if("sdrfOrder" %in% names(sdrf) && isTRUE(as.logical(sdrf["sdrfOrder"])) && length(setupSd$iniSdrfOrder)==ncol(abund) && ncol(abund) >1) {  # set order according to sdrf (only if >1 samples)
      nOrd <- order(setupSd$iniSdrfOrder)
      ## rename columns according to sdrf and set order of quant and abund ..
      abund <- abund[,nOrd]
      if(length(quant) >0) quant <- quant[,nOrd]
      if(length(setupSd$sampleNames)==ncol(quant)) {
        colNa <- colnames(abund) <- setupSd$sampleNames <- setupSd$sampleNaSdrf[nOrd]  #old# setupSd$sampleNames[nOrd]         ## take sample names from sdrf via  setupSd$sampleNaSdrf
        if(length(quant) >0) colnames(quant) <- setupSd$sampleNaSdrf[nOrd]  #old# setupSd$sampleNames[nOrd]
      } else colNa <- colnames(abund)
      ## now adapt order of setupSd, incl init Sdrf
      if(length(setupSd) >0) { 
        is2dim <- sapply(setupSd, function(x,le) length(dim(x))==2 && nrow(x)==le, le=length(nOrd))    # look for matr or df to change order of lines
        if(any(is2dim) >0) for(i in which(is2dim)) setupSd[[i]] <- setupSd[[i]][nOrd,]
        isVe <- sapply(setupSd, function(x,le) length(x)==le && length(dim(x)) <1, le=length(nOrd))    # look for vector to change order in setupSd
        if(any(isVe) >0) for(i in which(isVe)) setupSd[[i]] <- setupSd[[i]][nOrd] }
      gr <- gr[nOrd]

      if(length(counts) >0 && length(dim(counts))==3) counts <- array(counts[,nOrd,], dim=c(nrow(counts), length(nOrd), dim(counts)[3]), 
        dimnames=list(rownames(counts), colnames(counts)[nOrd], dimnames(counts)[[3]]))
      ## try re-adjusting levels
      tm1 <- sub("^[[:alpha:]]+( |_|-|\\.)+[[:alpha:]]+","", colnames(abund))  # remove heading text
      if(all(grepl("^[[:digit:]]", tm1))) {
        tm1 <- try(as.numeric(sub("( |_|-|\\.)*[[:alpha:]].*","", tm1)), silent=TRUE)   # remove tailing text and try converting to numeric
        if(!inherits(tm1, "try-error")) {
          setupSd$level <- match(tm1, sort(unique(tm1)))
          names(setupSd$level) <- tm1
          if(!silent) message(fxNa,"Sucessfully re-adjusted levels after bringing in order of Sdrf")}
      }     
    } else {

      ## harmonize sample-names/2
      colNa <- colnames(abund)
      chGr <- grepl("^X[[:digit:]]", colNa)                                                # check & remove heading 'X' from initial column-names starting with digits
      if(any(chGr)) colNa[which(chGr)] <- sub("^X","", colNa[which(chGr)])                 #
      colnames(quant) <- colNa
      if(length(abund) >0) colnames(abund) <- colNa  
    }  
    if(length(setupSd$sampleNames)==ncol(abund)) setupSd$sampleNames <- colNa #no#else setupSd$groups <- colNa
    if(length(dim(counts)) >1 && length(counts) >0) colnames(counts) <- colNa
  if(debug) {message(fxNa,"Read sample-meta data, rmc14"); rmc14 <- list()}

  ## main plotting of distribution of intensities
  custLay <- NULL
  if(is.numeric(plotGraph) && length(plotGraph) >0) {custLay <- as.integer(plotGraph); plotGraph <- TRUE} else {
      if(!isTRUE(plotGraph)) plotGraph <- FALSE}
  if(plotGraph) .plotQuantDistr(abund=abund, quant=quant, custLay=custLay, normalizeMeth=normalizeMeth, softNa="MassChroQ",
    refLi=refLi, refLiIni=refLiIni, tit=titGraph, silent=silent, callFrom=fxNa, debug=debug)


  ## meta-data
  notes <- c(inpFile=paFi, qmethod="MassChroQ", qMethVersion=if(length(infoDat) >0) unique(infoDat$Software.Revision) else NA,
    identType="protein", rawFilePath= if(length(infoDat) >0) infoDat$File.Name[1] else NA, normalizeMeth=normalizeMeth, call=deparse(match.call()),
    created=as.character(Sys.time()), wrProteo.version=paste(utils::packageVersion("wrProteo"), collapse="."), machine=Sys.info()["nodename"])
  ## prepare for final output
  if(isTRUE(separateAnnot)) list(raw=abund, quant=quant, annot=annot, counts=NULL, quantNotes=NULL, notes=notes) else data.frame(abund, annot)
}
  
