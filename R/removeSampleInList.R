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
  out <-  transpEl <- NULL
  datOK <- TRUE
  if(length(dat) <1 || !is.list(dat)) { datOK <- FALSE
    if(!silent) message(fxNa, msg[1:2])}
  if(length(remSamp) <1) {datOK <- FALSE
    if(!silent) message(fxNa, msg[3:2])}
  ## main
  if(datOK) {
    chLst <- listElem %in% names(dat)
    if(sum(chLst) <1) {
      warning("Can't find any of the list-elements defined via 'listElem' - nothing to do")
      datOK <- chRm <- FALSE
    } else {                # (some) list-elements fit dat
      listElemI <- wrMisc::naOmit(match(listElem, names(dat)))    # remove non-existing list-elements
      if(debug) { message(fxNa," rSIL1"); rSIL1 <- list(dat=dat,remSamp=remSamp,listElem=listElem,listElemI=listElemI)}
      if(length(listElemI) >1) {
        liDim <- lapply(dat[listElemI], dim)
        chLiDim <- sapply(liDim, length)
        if(any(chLiDim <1)) {    # remove linear (or lists like $sampleSetup)
          if(debug) message(fxNa,"Don't treat list-elements ",wrMisc::pasteC(names(dat)[which(chLiDim <2)], quoteC="'")," (not matrix, data.frame or array)")
          listElemI <- listElemI[-which(chLiDim <2)]
          liDim <- liDim[-which(chLiDim <2)]
      } }
      if(debug) { message(fxNa," rSIL2"); rSIL2 <- list()}
      if(length(listElemI) >1) {
        nCol <- sapply(liDim, function(x) x[2])
        chDupNC <- duplicated(nCol, fromLast=FALSE)
        if(any(!chDupNC[-1])) {
          useNCol <- names(which.max(table(nCol)))   # assume most frequent is good one !
          if(debug) message(fxNa,"Variable number of columns, using most frequent : ",useNCol," columns in ",wrMisc::pasteC(names(dat)[listElemI[which(nCol==useNCol)]], quoteC="'"))
          nRow <- sapply(liDim, function(x) x[1])
          chTra <- nCol == nRow & nCol != useNCol
          if(any(chTra)) {  # transpose case
            transpEl <- which(chTra)
            listElemI <- listElemI[-which(chTra)]
            if(debug) message(fxNa,"Matrix(es) ",wrMisc::pasteC(names(dat)[which(chTra)], quoteC="'"," to be trated as transposed ") )
          } else listElemI <- listElemI[which(nCol==useNCol)]
        }
      }
      if(debug) { message(fxNa," rSIL3"); rSIL3 <- list(dat=dat,remSamp=remSamp,listElem=listElem,listElemI=listElemI,transpEl=transpEl)}
      if(length(listElemI) >1) {      ## remove columns
        ## convert text-entries to index ?

        if(!is.integer(remSamp)) remSamp <- try(as.integer(remSamp), silent=TRUE)
        if(inherits(remSamp, "try-error")) stop("Invalid argument 'remSamp' (must be integer to design column(s) to be removed) !")
        chRm <- remSamp %in% 1:ncol(dat[[listElemI[1]]])
        if(all(!chRm)) stop("Invalid columns selected !")
        if(any(!chRm)) {
          if(!silent) message(fxNa,"Removing column(s) ",wrMisc::pasteC(remSamp[which(!chRm)], quoteC="'"," from 'remSamp' (not existing in data !)"))
          remSamp <- remSamp[which(chRm)]
        }
        if(debug) { message(fxNa," rSIL4"); rSIL4 <- list()}
        ## main removing of columns
        for(i in listElemI) { dat[[i]] <- if(length(dim(dat[[i]])) ==2) {if(length(remSamp) < ncol(dat[[i]]) -1) dat[[i]][,-remSamp] else {
            matrix(dat[[i]][,-remSamp], ncol=1, dimnames=list(rownames(dat[[i]]), colnames(dat[[i]])[-remSamp])) }
          } else if(length(remSamp) < ncol(dat[[i]]) -1) dat[[i]][,-remSamp,] else array(dat[[i]][,-remSamp,], dim=c(nrow(dat[[i]]),1, dim(dat[[i]])[3]))   # not fully tested ?          }
        }
        if(debug) message(fxNa,"Removed column(s) number ",wrMisc::pasteC(remSamp)," from matrix content")
        if(length(transpEl) >0) for(i in transpEl) dat[[i]] <- dat[[i]][-remSamp,]
        if(debug) { message(fxNa," rSIL5"); rSIL5 <- list(dat=dat,remSamp=remSamp,listElem=listElem,listElemI=listElemI,transpEl=transpEl,chRm=chRm)}

        if("sampleSetup" %in% listElem) {     # special case...
          if(debug) message(fxNa,"list-element $sampleSetup found, treat ..")
          for(i in c("lev","level","groups","sampleNames")) if(i %in% names(dat$sampleSetup)) dat$sampleSetup[[i]] <- dat$sampleSetup[[i]][-remSamp]
          if(length(dat$sampleSetup$sdrfDat) >0) dat$sampleSetup$sdrfDat <- if(length(remSamp) < nrow(dat$sampleSetup$sdrfDat) -1) dat$sampleSetup$sdrfDat[-remSamp,] else {
            matrix(dat$sampleSetup$sdrfDat[-remSamp,], nrow=1, dimnames=list(rownames(dat$sampleSetup$sdrfDat)[-remSamp], colnames(dat$sampleSetup$sdrfDat))) }
          if(length(dat$sampleSetup$annotBySoft) >0) dat$sampleSetup$annotBySoft <- if(length(remSamp) < nrow(dat$sampleSetup$annotBySoft) -1) dat$sampleSetup$annotBySoft[-remSamp,] else {
            matrix(dat$sampleSetup$annotBySoft[-remSamp,], nrow=1, dimnames=list(rownames(dat$sampleSetup$annotBySoft)[-remSamp], colnames(dat$sampleSetup$annotBySoft))) }
        }
      }
    }
  }
  if(datOK) dat else NULL }
   
