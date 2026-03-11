#' Remove Samples/Columns From List Of Matrixes 
#'
#' @description
#' Remove samples (ie columns) from every instance of list of matrixes, like from objects created with proteomics import-functions from this package.
#' Note: This function assumes same order of columns in list-elements 'listElem' !
#'
#' @param dat (list) main input to be filtered
#' @param remSamp (integer) column number to exclude
#' @param listElem (character) names of list-elements where columns indicated with 'remSamp' should be removed
#' @param silent (logical) suppress messages
#' @param debug (logical) display additional messages for debugging
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @return This function returns a matrix including imputed values or list of final and matrix with number of imputed by group (plus optional plot)
#' @seealso \code{\link{testRobustToNAimputation}}
#' @examples
#' set.seed(2019)
#' datT6 <- matrix(round(rnorm(300)+3,1), ncol=6, dimnames=list(paste("li",1:50,sep=""),
#'   letters[19:24]))
#' datL <- list(raw=datT6, quant=datT6, annot=matrix(nrow=nrow(datT6), ncol=2))
#' datDelta2 <- removeSampleInList(datL, remSam=2)
#' @export
removeSampleInList <- function(dat, remSamp, listElem=c("raw","quant","counts","sampleSetup"), silent=FALSE, debug=FALSE, callFrom=NULL) {
  ##
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="removeSampleInList")
  msg <- c("'dat' should be list or S3-object with $raw, $quant, $annot","; invalid entry - can't do anything ...","'remSamp' should be index of columns to remove")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  nSamp <- transpEl <- NULL  #out
  datOK <- TRUE
  if(length(dat) <1 || !is.list(dat)) { datOK <- FALSE
    if(!silent) message(fxNa, msg[1:2])}
  if(length(remSamp) <1) { datOK <- FALSE
    if(!silent) message(fxNa, msg[3:2])}
  ## main
  if(datOK) {
    if(any(is.na(listElem))) listElem <- wrMisc::naOmit(listElem)
    chLst <- if(length(listElem) ==0) FALSE else listElem %in% names(dat)
    if(debug) {message(fxNa, "rSIL1"); rSIL1 <- list(dat=dat,remSamp=remSamp,listElem=listElem,chLst=chLst)}
    if(all(!chLst)) {
      warning(fxNa,"Can't find any of the list-elements defined via 'listElem' - nothing to do")
      datOK <- chRm <- FALSE
    } else {                # (some) list-elements to consider for removing samples
      listElemI <- wrMisc::naOmit(match(listElem, names(dat)))    # remove non-existing list-elements
      if(debug) { message(fxNa," rSIL2"); rSIL2 <- list(dat=dat,remSamp=remSamp,listElem=listElem,listElemI=listElemI)}
      ## check for list-elements eligible
      if(length(listElemI) >1) {    ## special vase : $sampleSetup may be list
        liDim <- lapply(listElemI, function(x) dim(dat[[x]]) )
        chLiDim <- sapply(liDim, length)
        if(debug) { message(fxNa," rSIL3"); rSIL3 <- list(dat=dat,remSamp=remSamp,listElem=listElem,listElemI=listElemI,liDim=liDim,chLiDim=chLiDim)}
        if(any(chLiDim >1)) {
          ## check list-elements matching names more in detail (for same number of cols/rows to qualify) ; uses 1st listElemI ie (dat[[listElemI[1]]]) as reference number of columns to consider !
          nSamp <- ncol(dat[[listElemI[1]]])
          byDim <- list(byCol=sapply(liDim, function(x) which(x[2]==ncol(dat[[listElemI[1]]]))),
            byRow=sapply(liDim, function(x) which(x[1]==ncol(dat[[listElemI[1]]]))))
          byDim <- lapply(byDim, function(x) unlist(listElemI[which(sapply(x, length) >0)]))  
          if(debug) { message(fxNa," rSIL4"); rSIL4 <- list(dat=dat,remSamp=remSamp,listElem=listElem,listElemI=listElemI, nSamp=nSamp,byDim=byDim,liDim=liDim )}
          if(length(byDim$byCol) >0) {
            for(i in unlist(byDim$byCol)) {
              if(all(grepl("^[[:digit:]]+$", remSamp))) { remSamp <- try(as.numeric(remSamp), silent=TRUE)       # check removing colmns
                if(debug) message(fxNa,"Removing  col(s) ",wrMisc::pasteC(remSamp)," from ",i) 
                remSamp <- remSamp[which(remSamp >0 & remSamp <= ncol(dat[[i]]))]
              } else { remSamp <- wrMisc::naOmit(match(remSamp, colnames(dat[[i]])))}         
              if(length(remSamp) >0)  dat[[i]] <- if(length(dim(dat[[i]]))==2) dat[[i]][,-remSamp, drop=FALSE] else dat[[i]][,-remSamp, , drop=FALSE]}
          }
          if(debug) { message(fxNa," rSIL5"); rSIL5 <- list()}
          if(length(byDim$byRow) >0) {
            for(i in unlist(byDim$byRow)) {                                                         # check removing rows
              if(all(grepl("^[[:digit:]]+$", remSamp))) { remSamp <- try(as.numeric(remSamp)); 
                remSamp <- remSamp[which(remSamp >0 & remSamp <= nrow(dat[[i]]))]
              } else { remSamp <- wrMisc::naOmit(match(remSamp, colnames(dat[[i]])))}
              if(length(remSamp) >0)  dat[[i]] <- dat[[i]][-remSamp,, drop=FALSE] }
          }
          chVect <- sapply(liDim, length) ==1
          if(any(chVect) && length(nSamp) >0) for(i in which(chVect)) {if(length(dat[[listElemI[i]]])==nSamp) dat[[listElemI[i]]] <- dat[[listElemI[i]]][-remSamp] }    
          
        }
        if(debug) { message(fxNa," rSIL6"); rSIL6 <- list(dat=dat,remSamp=remSamp,listElem=listElem,listElemI=listElemI, nSamp=nSamp,byDim=byDim,liDim=liDim)} #,chSamAn=chSamAn
        if("sampleSetup" %in% names(dat) && "sampleSetup" %in% listElem) {
          for(i in 1 :length(dat$sampleSetup)) { chDim <- dim(dat$sampleSetup[[i]])
            if(length(chDim) ==2 && nrow(dat$sampleSetup[[i]]) >= remSamp) dat$sampleSetup[[i]] <- dat$sampleSetup[[i]][-1*remSamp,,drop=FALSE]
            if(length(chDim) ==3 && nrow(dat$sampleSetup[[i]]) >= remSamp) dat$sampleSetup[[i]] <- dat$sampleSetup[[i]][-1*remSamp,,,drop=FALSE]
            if(length(chDim) ==0 && length(dat$sampleSetup[[i]]) >=remSamp) dat$sampleSetup[[i]] <- dat$sampleSetup[[i]][-1*remSamp]            
          }   
        }
      }
    }
    dat
  } else NULL  
}        
          
