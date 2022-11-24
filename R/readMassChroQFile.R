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
#' @param sampleNames (character) new column-names for quantification data (ProteomeDiscoverer does not automatically use file-names from spectra)
#' @param refLi (character or integer) custom specify which line of data is main species, if character (eg 'mainSpe'), the column 'SpecType' in $annot will be searched for exact match of the (single) term given
#' @param separateAnnot (logical) if \code{TRUE} output will be organized as list with \code{$annot}, \code{$abund} for initial/raw abundance values and \code{$quant} with final normalized quantitations
#' @param tit (character) custom title to plot
#' @param graphTit (character) depreciated custom title to plot, please use 'tit'
#' @param wex (integer) relative expansion factor of the violin-plot (will be passed to \code{\link[wrGraph]{vioplotW}})
#' @param specPref (character or list) define characteristic text for recognizing (main) groups of species (1st for comtaminants - will be marked as 'conta', 2nd for main species- marked as 'mainSpe',
#'  and optional following ones for supplemental tags/species - maked as 'species2','species3',...);
#'  if list and list-element has multiple values they will be used for exact matching of accessions (ie 2nd of argument \code{annotCol})
#' @param separateAnnot (logical) if \code{TRUE} output will be organized as list with \code{$annot}, \code{$abund} for initial/raw abundance values and \code{$quant} with final normalized quantitations
#' @param gr (character or factor) custom defined pattern of replicate association, will override final grouping of replicates from \code{sdrf} and/or \code{suplAnnotFile} (if provided)   \code{}
#' @param sdrf (character, list or data.frame) optional extraction and adding of experimenal meta-data: if character, this may be the ID at ProteomeExchange. Besides, the output from \code{readSdrf} or a list from \code{defineSamples} may be provided; if \code{gr} is provided, it gets priority for grouping of replicates
#' @param suplAnnotFile (logical or character) optional reading of supplemental files; however, if \code{gr} is provided, \code{gr} gets priority for grouping of replicates;
#'  if \code{character} the respective file-name (relative or absolute path)
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
#'
#' @export
readMassChroQFile <- function(fileName, path=NULL, normalizeMeth="median", sampleNames=NULL, refLi=NULL, separateAnnot=TRUE, tit="MassChroQ", graphTit=NULL, wex=NULL,
  specPref=c(conta="CON_|LYSC_CHICK", mainSpecies="OS=Homo sapiens"), gr=NULL, sdrf=NULL, suplAnnotFile=FALSE, plotGraph=TRUE, silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## read MassChroQ (pre-)treated data
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readMassChroQFile")
  oparMar <- if(plotGraph) graphics::par("mar") else NULL       # only if figure might be drawn

  if(!isTRUE(silent)) silent <- FALSE
  if(!requireNamespace("utils", quietly=TRUE)) stop("package 'utils' not found ! Please install first")

  ## check & read file
  infoDat <- infoFi <- setupSd <- parametersD <- NULL        # initialize for sdrf annotation
  counts <- NULL                                             # so far PSM data are not accessible

  chPa <- length(grep("/",fileName)) >0 | length(grep("\\\\",fileName)) >0       # check for path already in fileName "
  if(length(path) <1) path <- "."
  paFi <- if(!chPa) file.path(path[1],fileName[1]) else fileName[1]              # use path only when no path joined to fileName
  chFi <- file.exists(paFi)
  if(!chFi) stop(" file ",fileName," was NOT found ",if(chPa) paste(" in path ",path)," !")
  tmp <- list()                                                        # initialize
  if(length(c(grep("\\.txt$",fileName), grep("\\.txt\\.z$",fileName))) >0) tmp[[1]] <- try(utils::read.delim(paFi, stringsAsFactors=FALSE), silent=TRUE)             # read tabulated text-file
  if(length(c(grep("\\.csv$",fileName), grep("\\.csv\\.gz$",fileName))) >0) tmp[[2]] <- try(utils::read.csv(paFi, stringsAsFactors=FALSE), silent=TRUE)               # read US csv-file
  if(length(c(grep("\\.csv$",fileName), grep("\\.csv\\.gz$",fileName))) >0) tmp[[3]] <- try(utils::read.csv2(paFi, stringsAsFactors=FALSE), silent=TRUE)              # read Euro csv-file
  if(length(c(grep("\\.tsv$",fileName), grep("\\.tsv\\.gz$",fileName))) >0) tmp[[4]] <- try(utils::read.csv(file=paFi, stringsAsFactors=FALSE, sep='\t', header=TRUE)) # read US comma tsv-file
  if(length(c(grep("\\.rda$",fileName), grep("\\.rdata$",tolower(fileName)))) >0) {
    ls1 <- ls()
    tmp[[5]] <- try(load(paFi))
    if(!inherits(tmp[[5]], "try-error")) {          # dont know under which name the object was saved in RData..
      if(length(ls1) +2 ==length(ls())) {
        tmp[[5]] <- get(ls()[which(!ls() %in% ls1 & ls() != "ls1")])            # found no way of removing initial object
        if(!silent) message(fxNa,"Loading object '",ls()[which(!ls() %in% ls1 & ls() != "ls1")],"' as quantification data out of ",fileName)
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

  ## normalize
  if(!is.matrix(tmp)) tmp <- as.matrix(tmp)
  abund <- wrMisc::normalizeThis(tmp, method=normalizeMeth, mode="additive", refLines=refLi, callFrom=fxNa)       #
  if(debug) {message(fxNa,"mc4")}

  ### GROUPING OF REPLICATES AND SAMPLE META-DATA
  if(any((!isFALSE(suplAnnotFile) & length(suplAnnotFile) >0), length(sdrf) >0)) {
    setupSd <- readSampleMetaData(sdrf=sdrf, suplAnnotFile=suplAnnotFile, quantMeth="MC", path=path, abund=utils::head(tmp), silent=silent, debug=debug, callFrom=fxNa)
  }
  if(debug) {message(fxNa,"mc13b")}

  ## finish groups of replicates & annotation setupSd
  if(length(setupSd) >0 & length(setupSd$groups) <1) {                       # if nothing found/matching from sdrf & file, try getting sample-setup from colnames (use if at least 1 replicate found)
    if(length(gr) ==ncol(abund)) setupSd$groups <- gr else {
      if(debug) {message(fxNa,"mc13e  Note: setupSd$groups is still empty ! ")}
      if(length(setupSd$lev) ==ncol(abund)) setupSd$groups <- setupSd$lev else {
        ## try defining groups based on colnames
        if(debug) {message(fxNa,"mc13f  Note: setupSd is still empty !  .. try getting sample-setup from colnames")}
        delPat <- "_[[:digit:]]+$|\\ [[:digit:]]+$|\\-[[:digit:]]+$"       # remove enumerators, ie trailing numbers after separator
        grou <- sub(delPat,"", colnames(abund))
        if(length(unique(grou)) >1 & length(unique(grou)) < ncol(abund)) setupSd$groups <- grou else {
          grou <- sub("[[:digit:]]+$","", colnames(abund))
          if(length(unique(grou)) >1 & length(unique(grou)) < ncol(abund)) setupSd$groups <- grou }
      } }
  }

  ## finish sample-names: use file-names from meta-data if no custom 'sampleNames' furnished
  ## One more check for colnames & sampleNames
  if(length(sampleNames) == ncol(abund)) {
    colnames(abund) <- colnames(abund) <- sampleNames     # use custom provided sampleNames
    if(length(dim(counts)) >1 & length(counts) >0) colnames(counts) <- sampleNames
  } else {
    if(debug) { message(fxNa,"mc13c") }
    sampleNames <- if(length(setupSd$annotBySoft$File.Name)==ncol(abund)) {    # try taking sampleNames from sdrf
      basename(sub("\\.raw$|\\.Raw$|\\.RAW$","", setupSd$annotBySoft$File.Name)) } else colnames(abund)     # otherwise remain with colnames from main file
    sampleNames <- wrMisc::trimRedundText(sampleNames, minNchar=2, spaceElim=TRUE, silent=silent, callFrom=fxNa, debug=debug)
    colnames(abund) <- colnames(abund) <- sampleNames
    if(length(dim(counts)) >1 & length(counts) >0) colnames(counts) <- sampleNames

  }
  if(debug) {message(fxNa,"Read sample-meta data, mc14"); mc14 <- list(sdrf=sdrf,suplAnnotFile=suplAnnotFile,abund=abund, abund=abund,refLi=refLi,annot=annot,setupSd=setupSd,sampleNames=sampleNames)}

  ## plot distribution of intensities
  custLay <- NULL
  if(length(plotGraph) >0) {if(is.numeric(plotGraph)) {custLay <- plotGraph; plotGraph <- TRUE
    } else  {plotGraph <- isTRUE(plotGraph[1]) }}
  if(plotGraph) {
    if(length(custLay) >0) graphics::layout(custLay) else graphics::layout(1:2)
    graphics::par(mar=c(3, 3, 3, 1))                           # mar: bot,le,top,ri
    misPa <- c(requireNamespace("wrGraph", quietly=TRUE), requireNamespace("sm", quietly=TRUE))
    if(is.null(tit)) tit <- "MassChroQ Quantification "
    titSu <- if(length(refLi) >0) paste(" by",length(refLi),"selected lines")  else NULL

    if(any(misPa)) {
      if(!silent) message(fxNa," missing package ",wrMisc::pasteC(c("wrGraph","sm")[which(misPa)],quoteC="'")," for drawing vioplots")
      ## wrGraph not available : simple boxplot
      graphics::boxplot(tmp, main=paste(tit," (initial)"), las=1, outline=FALSE)
      graphics::abline(h=round(stats::median(tmp,na.rm=TRUE)) +(-2:2), lty=2, col=grDevices::grey(0.6))
      ## now normalized
      graphics::boxplot(abund, main=paste(tit," (",normalizeMeth,"-normalized",titSu,")"), las=1, outline=FALSE)
      graphics::abline(h=round(stats::median(abund,na.rm=TRUE)) +(-2:2), lty=2, col=grDevices::grey(0.6))
    } else {                                                  # wrGraph and sm are available
      wrGraph::vioplotW(tmp, tit=paste(tit," (initial)"), wex=wex)
      graphics::abline(h=round(stats::median(tmp,na.rm=TRUE)) +(-2:2), lty=2, col=grDevices::grey(0.6))
      ## now normalized
      wrGraph::vioplotW((abund), tit=paste(tit," , ",normalizeMeth,"-normalized",titSu), wex=wex)
      graphics::abline(h=round(stats::median(abund,na.rm=TRUE)) +(-2:2), lty=2, col=grDevices::grey(0.6))
    }
    on.exit(graphics::par(mar=oparMar)) }   # restaure old settings

  ## meta-data
  notes <- c(inpFile=paFi, qmethod="MassChroQ", qMethVersion=if(length(infoDat) >0) unique(infoDat$Software.Revision) else NA,
   	rawFilePath= if(length(infoDat) >0) infoDat$File.Name[1] else NA, normalizeMeth=normalizeMeth, call=match.call(),
    created=as.character(Sys.time()), wrProteo.version=utils::packageVersion("wrProteo"), machine=Sys.info()["nodename"])
  ## prepare for final output
  if(isTRUE(separateAnnot)) list(raw=tmp, quant=abund, annot=annot, counts=NULL, quantNotes=NULL, notes=notes) else data.frame(abund, annot)
}
  
