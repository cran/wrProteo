#' Read proteomics meta-data as sdrf file 
#'  
#' This function allows reading proteomics meta-data from sdrf file, as they are provided on https://github.com/bigbio/proteomics-metadata-standard. 
#' Then, a data.frame with all annotation data will be returned. To stay conform with the (non-obligatory) recommendations, column-names will be shown as lower caps.  
#' The package utils must be installed.
#' 
#' @param fi (character) main input; may be full path or url to the file with meta-annotation. If a short project-name is given, 
#'   it will be searched based at the location of \code{urlPrefix}  
#' @param chCol (character, length=1) optional checking of column-names
#' @param urlPrefix (character, length=1) prefix to add to search when no complete path or url is given on \code{fi}, defaults to proteomics-metadata-standard on github
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of message(s) produced
#' @param debug (logical) display additional messages for debugging
#' @return This function returns the content of Sdrf-file as data.frame 
#' @seealso  in \code{\link[utils]{read.table}}
#' @examples
#' 
#' pxd001819 <- readSdrf("PXD001819")
#' str(pxd001819)
#' 
#'  
#' @export  
readSdrf <- function(fi, chCol="auto", urlPrefix="github", silent=FALSE, callFrom=NULL, debug=FALSE) {
  ## read proteomics meta-data as sdrf file 
  ## see https://github.com/bigbio/proteomics-metadata-standard
  ##  return data.frame, testing for 
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readSdrf")
  datOK <- if(length(fi) <1) FALSE else TRUE
  reqPa <- c("utils","wrMisc")
  chPa <- sapply(reqPa, requireNamespace, quietly=TRUE)
  if(any(!chPa)) { datOK <- FALSE; if(!silent) message(fxNa,"package(s) '",paste(reqPa[which(!chPa)], collapse="','"),"' not found ! Please install first from CRAN")}
  if(datOK) { if(length(fi) >1) fi <- fi[1]  
    if(is.na(fi)) { datOK <- FALSE
      if(!silent) message(fxNa," argument 'fi' is NA - nothing to do")} }
  if(datOK) {
    chFi <- file.exists(fi)
    if(!chFi & length(urlPrefix)==1 & !grepl("^https?://",fi)) {
      if(debug) {message(fxNa,"could not find as file, content of 'fi' expanded to url  rs1")
        rs1 <- list(fi=fi,chCol=chCol,urlPrefix=urlPrefix,datOK=datOK,chFi=chFi) }
      if(identical(urlPrefix,"github")) urlPrefix <- "https://github.com/bigbio/proteomics-metadata-standard/blob/master/annotated-projects/"
      fi <- paste0(urlPrefix, if(grepl("/",fi)) fi else paste0(fi,"/",sub("\\.sdrf\\.tsv$","",fi),".sdrf.tsv") )
  } }
  ## Main reading
  if(datOK) { 
    if(debug & !grepl("\\.sdrf",fi)) message(fxNa,"Trouble ahead, '",fi,"' does not contain '.sdrf' ...")
    out <- suppressWarnings(try(utils::read.delim(wrMisc::gitDataUrl(fi), sep='\t', header=TRUE, fill=TRUE), silent=!isTRUE(debug)))  
    if(inherits(out, "try-error")) { message(fxNa," failed reading '",fi,"'  (possibly bad url/path ?)"); return(NULL)
    } else {
      fi2 <- sub("[[:print:]]+/","", sub("\\.sdrf\\.tsv","",fi))
      if(any(length(dim(out)) !=2, dim(out) < 2:1)) { message(fxNa," data in bad format")
      } else {
        ## diagnostic for expected column-names
        colnames(out) <- tolower(colnames(out))
        if(any(sapply(c("auto","def","default"), identical, chCol))) chCol <- c("source.name", "assay.name",
          "characteristics.biological.replicate.","characteristics.organism.", "comment.data.file.","comment.file.uri." )
          #"characteristics.spiked.compound.","factor.value.spiked.compound.")
        if(length(chCol) >0) locCol <- match(chCol, colnames(out))
        if(any(is.na(locCol))) message(fxNa,"File ",fi2,"  Can't find columns ",wrMisc::pasteC(chCol[which(is.na(locCol))], quoteC="'"))    
        return(out) } }
  } 
} 
  
