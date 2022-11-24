#' Complement Missing EntryNames In Annotation 
#'
#' This function helps replacing missing EntryNames (in $annot) after reading quantification results. 
#' To do so the comumn-names of \code{annCol} will be used : 
#' The content of 2nd element (and optional 3rd element) will be used to replace missing content in column defined by 1st element.
#' 
#' @param x (list) output of \code{readMaxQuantFile}, \code{readProtDiscovFile} or \code{readProlineFile}. 
#'   This list must be a matrix and contain $annot with the columns designated in \code{annCol}.
#' @param annCol (character) the column-names form \code{x$annot}) which will be used : The first column designs the
#'   column where empty fields are searched and the 2nd and (optional) 3rd will be used to fill the empty spots in the st column
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @return This function returns a list (like as input), but with missing elments of $annot completed (if available in other columns) 
#' @seealso \code{\link{readMaxQuantFile}}, \code{\link{readProtDiscovFile}}, \code{\link{readProlineFile}} 
#' @examples
#' dat <- list(quant=matrix(sample(11:99,9,replace=TRUE), ncol=3), annot=cbind(EntryName=c(
#'   "YP010_YEAST","",""),Accession=c("A5Z2X5","P01966","P35900"), SpecType=c("Yeast",NA,NA)))
#' replMissingProtNames(dat)
#' @export
replMissingProtNames <- function(x, annCol=c("EntryName","Accession","SpecType"), silent=FALSE, callFrom=NULL) {
  ## replace in $annot missing EntryNames by concatenating Accession + SpecType (ie 2nd & 3rd of annCol)
  ## move to wrProteo ?
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="replreplMissingProtNames")
  msg <- "argument 'x' should be list containg list-element called 'annot' (matrix), as produced by readMaxQuantFile(), readProtDiscovFile() etc"
  if(!is.list(x) | length(x) <1) stop(msg)
  if(!"annot" %in% names(x)) stop(msg)  
  if(any(length(dim(x$annot)) !=2, dim(x$annot) < 2:3)) stop("x$annot must be matrix or data.frame with min 2 lines and 3 cols")
  ## main  
  chCol <- annCol %in% colnames(x$annot)
  if(all(chCol[1:2])) {
    if(length(chCol) >2) {if(!chCol[3]) chCol <- chCol[1:2] }   # omit 3rd column-name if not present in x$annot
    chNA1 <- is.na(x$annot[,annCol[1]]) | x$annot[,annCol[1]]==""
    chNA2 <- is.na(x$annot[,annCol[2]]) | x$annot[,annCol[2]]==""
    if(length(chCol) >2) {chNA3 <- is.na(x$annot[,annCol[3]]) | x$annot[,annCol[3]]==""
      chNA <- chNA1 & (!chNA2 | !chNA3)
    } else chNA <-  chNA1 & !chNA2  
    if(any(chNA)) { 
      if(!silent) message(fxNa," ..trying to replace ",sum(chNA)," '",annCol[1],"'")
        if(length(chCol) >2) {
          x$annot[which(chNA),annCol[1]] <- paste(x$annot[which(chNA),annCol[2]], x$annot[which(chNA),annCol[3]], sep="_")
          x$annot[which(chNA),annCol[1]] <- sub("_","", sub("_NA","", sub("NA_","", sub("NA NA","NA", x$annot[which(chNA),annCol[1]]))))
        } else x$annot[which(chNA),annCol[1]] <- x$annot[which(chNA),annCol[2]] } 
  } else message(fxNa," Nothing to do. Column-names ",wrMisc::pasteC(annCol[which(!chCol)],quoteC="'")," not found in x$annot !")
  x }
   
