#' Read csv or txt files exported from MS-Angel and Proline
#'
#' Quantification results form MS-Angel and Proline \href{http://www.profiproteomics.fr/proline/}{Proline} should be first saved via Excel or LibreOffice as csv or tabulated txt. 
#' Such files can be read by this function and relevant information be extracted. 
#' The final output is a list containing 3 elements: \code{$annot}, \code{$abund} and optional \code{$quant}, or returns data.frame with entire content of file if \code{separateAnnot=FALSE}.
#' Note: There is no normalization by default since quite frequently data produced by Proline are already sufficiently normalized. 
#' In case of doubt the figure prouced using the argument \code{plotGraph=TRUE} may help judging if distribtions are aligned suffiently well.  
#' 
#' @param fileName (character) name of file to read
#' @param path (character) optional path (note: Windows backslash sould be protected or written as '/')
#' @param logConvert (logical) convert numeric data as log2, will be placed in $quant
#' @param quantCol (character or integer) exact col-names, or if length=1 content of \code{quantCol} will be used as pattern to search among column-names for $quant using \code{grep} 
#' @param annotCol (character) (character) exact col-names or if length=1 pattern to search among column-names for $annot
#' @param separateAnnot (logical) separate annotation form numeric data (quantCol and annotCol must be defined)
#' @param refLi (integer) custom decide which line of data is main species, if single character entry it will be used to choose a group of species (eg 'mainSpe')
#' @param plotGraph (logical or matrix of integer) optional plot vioplot of initial data; if integer, it will be passed to \code{layout} when plotting
#' @param tit (character) custom title to plot
#' @param graphTit (character) (depreciated custom title to plot), please use 'tit'
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @return list with \code{$annot}, \code{$raw} and optional \code{$quant}, or returns data.frame with entire content of file if \code{separateAnnot=FALSE}
#' @seealso \code{\link[utils]{read.table}} 
#' @examples
#' path1 <- system.file("extdata", package="wrProteo")
#' fiNa <- "exampleProlineABC.csv"
#' dataABC <- readProlineFile(file.path(path1,fiNa))
#' summary(dataABC$quant)
#' matrixNAinspect(dataABC$quant, gr=as.factor(substr(colnames(dataABC$quant),1,1))) 
#' @export
readProlineFile <- function(fileName,path=NULL,logConvert=TRUE,quantCol="^abundance_",annotCol=c("accession","description","is_validated","coverage","X.sequences","X.peptides","protein_set.score"), 
  refLi=NULL,separateAnnot=TRUE,plotGraph=TRUE,tit=NULL,graphTit=NULL,silent=FALSE,callFrom=NULL){
  ## 'quantCol', 'annotCol' (character) exact col-names or if length=1 pattern to search among col-names for $quant or $annot
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readProlineFile")
  opar <- graphics::par(no.readonly=TRUE)
  ## check & read file
  chPa <- length(grep("/",fileName)) >0 | length(grep("\\\\",fileName)) >0      # check for path already in fileName "
  if(length(path) <1) path <- "."
  paFi <- if(!chPa) file.path(path[1],fileName[1]) else fileName[1]             # use path only when no path joined to fileName
  chFi <- file.exists(paFi)
  if(!chFi) stop(" file ",fileName," was NOT found ",if(length(path) >0) paste(" in path ",path)," !")
  if(length(grep("\\.xlsx$",fileName)) >0) message(fxNa," Trouble ahead, extracting out of Excel should be done via saving as csv or txt !!")
  ## try to find out which input format; read 3x and choose the one with most columns 
  tmp <- list()
  if(length(grep("\\.txt$",fileName)) >0) tmp[[1]] <- try(utils::read.delim(paFi, stringsAsFactors=FALSE), silent=TRUE)             # read tabulated text-file
  if(length(grep("\\.csv$",fileName)) >0) tmp[[2]] <- try(utils::read.csv(paFi, stringsAsFactors=FALSE), silent=TRUE)               # read US csv-file
  if(length(grep("\\.csv$",fileName)) >0) tmp[[3]] <- try(utils::read.csv2(paFi, stringsAsFactors=FALSE), silent=TRUE)              # read Euro csv-file
  chCl <- sapply(tmp,class) =="try-error" 
  if(length(chCl) <1) stop("Failed to recognize file extensions of input data (unknown format)")
  if(any(chCl)) {if(all(chCl)) stop(" Failed to extract data (unknown format) from ",fileName)}
  nCol <- sapply(tmp, function(x) if(length(x) >0) {if(class(x) != "try-error") ncol(x) else NA} else NA)
  bestT <- which.max(nCol)
  out <- tmp[[bestT]]
  if(any(c(length(quantCol), length(annotCol)) <1)) separateAnnot <- FALSE else {
    if(any(all(is.na(quantCol)), all(is.na(annotCol)))) separateAnnot <- FALSE
  }
  
  metaCo <- wrMisc::naOmit(if(length(annotCol) >1) wrMisc::extrColsDeX(out, extrCol=annotCol, doExtractCols=FALSE, callFrom=fxNa) else grep(annotCol,colnames(out)))
  ## locate & extract abundance/quantitation data
  if(length(quantCol) >1) { abund <- as.matrix(wrMisc::extrColsDeX(out, extrCol=quantCol, doExtractCols=TRUE, callFrom=fxNa))
  } else {
    quantCol <- grep(quantCol, colnames(out)) 
    chNa <- is.na(quantCol)
    if(all(chNa)) stop("Could not find any of of the columns specified in argument 'quantCol' !")
    if(any(chNa)) { 
      if(!silent) message(fxNa," Could not find columns ",wrMisc::pasteC(quantCol[which(chNa)],quote="'")," .. omit")
      quantCol <- wrMisc::naOmit(quantCol)} 
    abund <- as.matrix(out[,quantCol]) }           # abundance val
  ## check abundance/quantitation data
  chNum <- is.numeric(abund)
  if(!chNum) {abund <- apply(out[,quantCol], 2, wrMisc::convToNum, convert="allChar", callFrom=fxNa)}    
  if(is.character(refLi) & length(refLi)==1) refLi <- which(out[,"Spec"]==refLi)   # may be "mainSpe"
  ## don't normalize here by default since data are typically sufficiently normalized  
  colnames(abund) <- wrMisc::.trimFromStart(wrMisc::.trimFromEnd(colnames(abund)))
  ## plot distribution of intensities
  custLay <- NULL
  if(length(plotGraph) >0) {if(is.numeric(plotGraph)) {custLay <- plotGraph; plotGraph <- TRUE
    } else { plotGraph <- as.logical(plotGraph[1])}}
  if(plotGraph){
    if(length(custLay) >0) graphics::layout(custLay)
    graphics::par(mar=c(3, 3, 3, 1))                          # mar: bot,le,top,ri
    if(length(graphTit) >0) message(fxNa,"argument 'graphTit' is depreciated, please rather use 'tit' ")
    if(is.null(tit) & !is.null(graphTit)) tit <- graphTit
    if(is.null(tit)) tit <- "Distribution of quantification values"
    chGr <- try(find.package("wrGraph"), silent=TRUE)
    chSm <- try(find.package("sm"), silent=TRUE)
    misPa <- c("try-error" %in% class(chGr),"try-error" %in% class(chSm))
    if(any(misPa)) { 
      if(!silent) message(fxNa," missing package ",wrMisc::pasteC(c("wrGraph","sm")[which(misPa)],quoteC="'")," for drawing vioplots")
      ## wrGraph not available : simple boxplot  
      graphics::boxplot(log2(abund), main=tit, las=1, outline=FALSE) 
      graphics::abline(h=round(log2(stats::median(abund, na.rm=TRUE))) +c(-1:1), lty=2, col=grDevices::grey(0.6))      
    } else {                                          # wrGraph & sm are available
      wrGraph::vioplotW(log2(abund), tit=paste(tit," (initial)",sep=" ")) 
      graphics::abline(h=round(stats::median(abund, na.rm=TRUE)) +(-1:1), lty=2, col=grDevices::grey(0.6))           
    }
    on.exit(graphics::par(opar)) }                         #
  ## meta-data
  notes <- c(qmethod="Proline", normalizeMeth="none", call=match.call(), created=as.character(Sys.time()), 
    wrProteo.version=utils::packageVersion("wrProteo"), machine=Sys.info()["nodename"])
  ##
  if(separateAnnot) {
    if(!is.numeric(abund) & logConvert) {message(fxNa," Problem: Abundance data seem not numeric, can't transform log2 !")}
    out <- list(raw=abund, quant=abund, annot=as.matrix(out[,metaCo]), notes=notes)
    if(logConvert) {
      ch0 <- which(out$quant <=0)
      out$quant <- abund
      if(length(ch0) >0) { out$quant[ch0] <- NA
        if(!silent) message(fxNa," NOTE : ",length(ch0)," elements of '",fileName,"' are 0 or negative, will be transformed to NA for log2 transformation")}
      out$quant <- log2(out$quant) }
    }
  ## final  
  out }
   
