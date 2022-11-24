#' Read proteomics meta-data as sdrf file
#'
#' @description This function allows reading proteomics meta-data from sdrf file, as they are provided on https://github.com/bigbio/proteomics-metadata-standard. 
#' A data.frame containing all annotation data will be returned. To stay conform with the (non-obligatory) recommendations, columnnames are shown as lower caps.  
#' 
#' @details The packages utils and wrMisc must be installed.
#' Please note that reading sdrf files (if not provided as local copy) will take a few seconds, depending on the responsiveness of github.
#' 
#' @param fi (character) main input; may be full path or url to the file with meta-annotation. If a short project-name is given, 
#'   it will be searched based at the location of \code{urlPrefix}  
#' @param chCol (character, length=1) optional checking of column-names
#' @param urlPrefix (character, length=1) prefix to add to search when no complete path or url is given on \code{fi}, defaults to proteomics-metadata-standard on github
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of messages produced
#' @param debug (logical) display additional messages for debugging
#' @return This function returns the content of sdrf-file as data.frame (or \code{NULL} if the corresponding file was not found) 
#' @seealso  in \code{\link[utils]{read.table}}
#' @examples 
#' ## This may take a few sconds...
#' sdrf001819 <- readSdrf("PXD001819") 
#' str(sdrf001819)
#' 
#'  
#' @export  
readSdrf <- function(fi, chCol="auto", urlPrefix="github", silent=FALSE, callFrom=NULL, debug=FALSE) {
  ## read proteomics meta-data as sdrf file 
  ## see https://github.com/bigbio/proteomics-metadata-standard
  ##  return data.frame, testing for 
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readSdrf")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  
  datOK <- if(length(fi) <1) FALSE else TRUE
  reqPa <- c("utils","wrMisc")
  chPa <- sapply(reqPa, requireNamespace, quietly=TRUE)
  if(any(!chPa)) { datOK <- FALSE; if(!silent) message(fxNa,"package(s) '",paste(reqPa[which(!chPa)], collapse="','"),"' not found ! Please install first from CRAN")}
  if(datOK) { if(length(fi) >1) fi <- fi[1]                   # make length=1
    if(is.na(fi)) { datOK <- FALSE
      if(!silent) message(fxNa," argument 'fi' is NA - nothing to do")} }
  if(datOK) {
    chFi <- file.exists(fi)
    if(!chFi & length(urlPrefix)==1 & !grepl("^https?://",fi)) {
      if(debug) {message(fxNa," rs1"); rs1 <- list(fi=fi,chCol=chCol,urlPrefix=urlPrefix,datOK=datOK,chFi=chFi)}
      if(identical(urlPrefix,"github")) urlPrefix <- "https://github.com/bigbio/proteomics-metadata-standard/blob/master/annotated-projects/"
      fi <- paste0(urlPrefix, if(grepl("/",fi)) fi else paste0(fi,"/",sub("\\.sdrf\\.tsv$","",fi),".sdrf.tsv") )
      if(debug) message(fxNa,"Could not find as file, content of 'fi' expanded to url ",fi," ")
  } }


  ## Main reading
  if(datOK) { 
    if(debug & !grepl("\\.sdrf",fi)) message(fxNa,"Trouble ahead ?  '",fi,"' does not contain '.sdrf' ...")
    out <- suppressWarnings(try(utils::read.delim(wrMisc::gitDataUrl(fi), sep='\t', header=TRUE, fill=TRUE), silent=!isTRUE(debug)))
    if(inherits(out, "try-error") & grepl("pxd", fi)) { 
      if(!silent) message(fxNa,"First try not successful; trying rather as '", gsub("pxd","PXD", fi)," (instead of '",fi,"')")
      fi <- gsub("pxd","PXD", fi)	
      out <- suppressWarnings(try(utils::read.delim(wrMisc::gitDataUrl(fi), sep='\t', header=TRUE, fill=TRUE), silent=!isTRUE(debug)))    
    }
    if(debug) {message(fxNa," rs2"); rs2 <- list(fi=fi,out=out,chCol=chCol,urlPrefix=urlPrefix,datOK=datOK,chFi=chFi)}
    if(inherits(out, "try-error")) { message(fxNa," FAILED reading '",fi,"'  (possibly bad url/path ?)"); return(NULL)
    } else {
      fi2 <- sub("[[:print:]]+/","", sub("\\.sdrf\\.tsv","",fi))
      if(debug) {message(fxNa," rs3"); rs3 <- list(fi=fi,out=out,fi2=fi2,chCol=chCol,urlPrefix=urlPrefix,datOK=datOK,chFi=chFi)}

      if(any(length(dim(out)) !=2, dim(out) < 1, na.rm=TRUE)) { message(fxNa," data in bad format")           # min 1 line, 1 col
      } else {
        ## diagnostic for expected column-names
        colnames(out) <- tolower(colnames(out))
        if(any(sapply(c("auto","def","default"), identical, chCol), na.rm=TRUE)) chCol <- c("source.name", "assay.name",
          "characteristics.biological.replicate.","characteristics.organism.", "comment.data.file.","comment.file.uri." )
        if(length(chCol) >0) locCol <- match(chCol, colnames(out))
        if(any(is.na(locCol))) message(fxNa,"File ",fi2,"  Can't find columns ",wrMisc::pasteC(chCol[which(is.na(locCol))], quoteC="'"))    
        if(debug) message(fxNa,"Successfully read ",ncol(out)," annotation columns for ",nrow(out)," samples")
        return(out) } }
  } 
} 
    
