#' Remove Samples/Columns From list of matrixes
#'   
#' @description
#' Remove samples (ie columns) from every instance of list of matrixes.
#' Note: This function assumes same order of columns in list-elements 'listElem' !
#'  
#' @param dat (list) main input to be filtered
#' @param remSamp (integer) column number to exclude
#' @param listElem (character) names of list-elements where columns indicated with 'remSamp' should be removed
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This function returns a matrix including imputed values or list of final and matrix with number of imputed by group (plus optional plot)
#' @seealso \code{\link{testRobustToNAimputation}}  
#' @examples
#' set.seed(2019)
#' datT6 <- matrix(round(rnorm(300)+3,1), ncol=6, dimnames=list(paste("li",1:50,sep=""),
#'   letters[19:24]))
#' datL <- list(abund=datT6, quant=datT6, annot=matrix(nrow=nrow(datT6), ncol=2)) 
#' datDelta2 <- removeSampleInList(datL, remSam=2)
#' @export
removeSampleInList <- function(dat, remSamp, listElem=c("abund","quant"), silent=FALSE, callFrom=NULL) {
  ##  
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="removeSampleInList")
  msg <- "'dat' should be list or S3-object with $abund, $quant, $annot; invalid entry - can't do anything"
  if(!isTRUE(silent)) silent <- FALSE
  datOK <- FALSE
  if(length(dat) >0 & is.list(dat)) datOK <- TRUE
  if(datOK) {  
    chLst <- listElem %in% names(dat)
    if(sum(chLst) <1) stop("Can't find any of the list-elements defined in 'listElem' - nothing to do")
    listElem <- listElem[which(chLst)]
    remSamp <- wrMisc::convToNum(remSamp)
    chRm <- 1:ncol(dat$quant) %in% remSamp
    if(any(chRm)) {
      for(i in listElem) dat[[i]] <- if(sum(!chRm) >1) dat[[i]][,-1*which(chRm)] else matrix(dat[[i]][,-1*which(chRm)],
        nrow=nrow(dat[[i]]), dimnames=list(rownames(dat[[i]]),colnames(dat[[i]])[-1*which(chRm)])) }
    } else if(!silent) message(fxNa,msg)      
  dat }
   
