#' Read csv or txt files exported from MS-Angel and Proline
#'
#' Quantification results form MS-Angel and Proline \href{http://proline.profiproteomics.fr/}{Proline} should be first saved via Excel or LibreOffice as csv or tabulated txt. 
#' Such files can be read by this function and relevant information be extracted. 
#' The final output is a list containing 3 elements: \code{$annot}, \code{$abund} and optional \code{$quant}, or returns data.frame with entire content of file if \code{separateAnnot=FALSE}. 
#' 
#' @param fileNa (character) name of file to read
#' @param wdir (character) optional path (note: Windows backslash sould be protected or written as '/')
#' @param logConvert (logical) convert numeric data as log2, will be placed in $quant
#' @param quantCol (character) (character) exact col-names or if length=1 pattern to search among column-names for $quant 
#' @param annotCol (character) (character) exact col-names or if length=1 pattern to search among column-names for $annot
#' @param separateAnnot (logical) separate annotation form numeric data (quantCol and annotCol must be defined)
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message produced
#' @return list with \code{$annot}, \code{$abund} and optional \code{$quant}, or returns data.frame with entire content of file if \code{separateAnnot=FALSE}
#' @seealso \code{\link[utils]{read.table}} 
#' @examples
#' path1 <- system.file("extdata",package="wrProteo")
#' fiNa <- "exampleProlineABC.csv"
#' dataABC <- readProlineFile(file.path(path1,fiNa))
#' summary(dataABC$abund)
#' matrixNAinspect(dataABC$quant,gr=as.factor(substr(colnames(dataABC$abund),1,1))) 
#' @export
readProlineFile <- function(fileNa,wdir=NULL,logConvert=TRUE,quantCol="^abundance_",annotCol=c("accession","description","is_validated","coverage","X.sequences","X.peptides","protein_set.score"), separateAnnot=TRUE,silent=FALSE,callFrom=NULL){
  ## 'quantCol', 'annotCol' (character) exact col-names or if length=1 pattern to search among col-names for $quant or $annot
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="readProlineFile")
  chPa <- length(grep("/",fileNa)) >0 | length(grep("\\\\",fileNa)) >0
  paFi <- if(length(wdir) >0) file.path(wdir,fileNa) else fileNa
  chFi <- file.exists(paFi)
  if(!chFi) stop(" file ",fileNa," was NOT found ",if(length(wdir) >0) paste(" in path ",wdir)," !")
  if(length(grep("\\.xlsx$",fileNa)) >0) message(fxNa," Trouble ahead, extracting out of Excel should be done via saving as csv or txt !!")
  tmp <- list()
  if(length(grep("\\.txt$",fileNa)) >0) tmp[[1]] <- try(utils::read.delim(paFi,stringsAsFactors=FALSE),silent=TRUE)             # read tabulated text-file
  if(length(grep("\\.csv$",fileNa)) >0) tmp[[2]] <- try(utils::read.csv(paFi,stringsAsFactors=FALSE),silent=TRUE)               # read US csv-file
  if(length(grep("\\.csv$",fileNa)) >0) tmp[[3]] <- try(utils::read.csv2(paFi,stringsAsFactors=FALSE),silent=TRUE)              # read Euro csv-file
  chCl <- sapply(tmp,class) =="try-error" 
  if(length(chCl) <1) stop("Failed to recognize file extesions of inout data (unknown format)")
  if(any(chCl)) {if(all(chCl)) stop(" Failed to extract data (unknown format) from ",fileNa)}
  nCol <- sapply(tmp,function(x) if(length(x) >0) {if(class(x) != "try-error") ncol(x) else NA} else NA)
  bestT <- which.max(nCol)
  out <- tmp[[bestT]]
  if(any(c(length(quantCol),length(annotCol)) <1)) separateAnnot <- FALSE else {
    if(any(all(is.na(quantCol)),all(is.na(annotCol)))) separateAnnot <- FALSE
  }
  ##
  if(separateAnnot) {
    metaCo <- wrMisc::naOmit(if(length(annotCol) >1) match(annotCol,colnames(out)) else grep(annotCol,colnames(out)))
    quantCo <- wrMisc::naOmit(if(length(quantCol) >1) match(quantCol,colnames(out)) else grep(quantCol,colnames(out)))
    out <- list(abund=as.matrix(out[,quantCo]),annot=as.matrix(out[,metaCo]))
    colnames(out$abund) <- wrMisc::.trimFromStart(wrMisc::.trimFromEnd(colnames(out[[1]])))
    if(logConvert) {
      ch0 <- which(out$abund <=0)
      out$quant <- out$abund
      if(length(ch0) >0) { out$quant[ch0] <- NA
        if(!silent) message(fxNa," NOTE : ",length(ch0)," elements of '",fileNa,"' are 0 or negative, will be transformed to NA for log2 transformation")}
      out$quant <- log2(out$quant) }
    }
  out }
   
