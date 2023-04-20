#' Read proteomics meta-data as sdrf file
#'
#' @description This function allows reading proteomics meta-data from sdrf file, as they are provided on https://github.com/bigbio/proteomics-sample-metadata.
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
  ## see https://github.com/bigbio/proteomics-sample-metadata
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
    flexFileNa <- TRUE
    chFi <- file.exists(fi)
    if(!chFi && length(urlPrefix)==1 && !grepl("^https?://",fi)) {
      if(debug) {message(fxNa," rs1"); rs1 <- list(fi=fi,chCol=chCol,urlPrefix=urlPrefix,datOK=datOK,chFi=chFi)}
      if(identical(urlPrefix,"github")) urlPrefix <- "https://github.com/bigbio/proteomics-sample-metadata/blob/master/annotated-projects/"
      if(grepl("/",fi)) {        # if abs or relative path- do not adjust lower/upper case
        if(flexFileNa) fi <- c(fi, file.path(dirname(fi),"sdrf.tsv"))
      } else {            # simple name, need to add folder of pxd project
        projN <- toupper(sub("\\.sdrf\\.tsv$","",tolower(fi)))
        fi <- paste0(projN,".sdrf.tsv")    #  set 'pxd'-part as upper case
        if(isTRUE(flexFileNa)) fi <- c(fi,"sdrf.tsv")
        fi <- paste0(urlPrefix,"/",projN,"/", fi )
      }
      if(debug) message(fxNa,"Could not find as file, content of 'fi' expanded to url ",wrMisc::pasteC(fi)," ")
  } }
  if(debug) {message(fxNa," rs1b"); rs1b <- list(fi=fi,chCol=chCol,urlPrefix=urlPrefix,datOK=datOK,chFi=chFi)}

  ## Main reading
  if(datOK) {
    if(debug && !grepl("\\.sdrf",fi[1])) message(fxNa,"Trouble ahead ?  '",fi,"' does not contain '.sdrf' ...")
    out <- suppressWarnings(try(utils::read.delim(wrMisc::gitDataUrl(fi[1]), sep='\t', header=TRUE, fill=TRUE), silent=!isTRUE(debug)))
    if(inherits(out, "try-error")) {
      if(any(grepl("pxd",fi))) {
        if(debug) message(fxNa,"First try not successful; trying rather as '", gsub("pxd","PXD", fi[1])," (instead of '",fi[1],"')")
        fi[1] <- gsub("pxd","PXD", fi[1])
        out <- suppressWarnings(try(utils::read.delim(wrMisc::gitDataUrl(fi[1]), sep='\t', header=TRUE, fill=TRUE), silent=!isTRUE(debug)))}
      if(inherits(out, "try-error") && length(fi) >1) {
        if(debug) message(fxNa,"So far not successful; trying rather as '",fi[2],"')")
        out <- suppressWarnings(try(utils::read.delim(wrMisc::gitDataUrl(fi[2]), sep='\t', header=TRUE, fill=TRUE), silent=!isTRUE(debug)))
        if(!inherits(out, "try-error")) fi <- fi[2] }
    }
    if(debug) {message(fxNa," rs2"); rs2 <- list(fi=fi,out=out,chCol=chCol,urlPrefix=urlPrefix,datOK=datOK,chFi=chFi)}
    if(inherits(out, "try-error")) { message(fxNa," FAILED reading '",fi[1],"'  (possibly bad url/path ?)"); return(NULL)
    } else {
      fi2 <- sub("[[:print:]]+/","", sub("\\.sdrf\\.tsv","",fi))   # pxd part
      if(debug) {message(fxNa," rs3"); rs3 <- list(fi=fi,out=out,fi2=fi2,chCol=chCol,urlPrefix=urlPrefix,datOK=datOK,chFi=chFi)}

      if(any(length(dim(out)) !=2, dim(out) < 1, na.rm=TRUE)) { message(fxNa," data in bad format")           # min 1 line, 1 col
      } else {
        ## diagnostic for expected column-names
        colnames(out) <- tolower(colnames(out))
        if(any(sapply(c("auto","def","default"), identical, chCol), na.rm=TRUE)) chCol <- c("source.name", "assay.name",
          "characteristics.biological.replicate.","characteristics.organism.", "comment.data.file.","comment.file.uri." )
        if(length(chCol) >0) locCol <- match(chCol, colnames(out))
        if(any(is.na(locCol)) && !silent) message(fxNa,"Data-Annotation ",fi2,"  Can't find column(s) ",wrMisc::pasteC(chCol[which(is.na(locCol))], quoteC="'"))
        if(!silent) message(fxNa,"Successfully read ",ncol(out)," annotation columns for ",nrow(out)," samples")
        return(out) } }
  }
}
  
