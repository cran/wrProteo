#' Find On-Off Situations In Pairwise Comparisons
#'
#' This function aims to find instances from pairwise comparisons where one group of samples is entirely absent (at least \code{lim} percent of \code{NA}) 
#' and all present (no more than 1 - \code{lim} persent \code{NA}) in the second group.
#'  
#' @details
#' Cases where all samples of a given group all contain a given protein while all samples of another group are all absent below detction limit (ie \code{NA}) 
#' are great interest for interpretation but are difficult to identify with t-test from imputed data.
#' This function allows counting of such constellations with a given tolerance (\code{lim}) for a given percentage of \code{NA}s.
#' 
#' Resulting counts can be returned as matrix or as list-element $onOff to the input object \code{stat} (typically from function \code{testRobustToNAimputation}).
#'
#' @param stat ('MArrayLM'-object or list) designed for the output from \code{moderTest2grp} or \code{moderTestXgrp}
#' @param lim (numeric) the threshold to be applied on \code{statTy} for the result of the statistcal testing (after multiple testing correction)
#' @param compNo (integer) so far without any effect, in future to select the comparison(s) name/number/index to be used, defaults to all comparisons
#' @param asFullObj (logical) return entire 'MArrayLM'-object or list with additional $onOff (default) or direcly the onOff results
#' @param debug (logical) supplemental messages for debugging
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This function returns (depending on argument \code{asFullObj}) an 'MArrayLM'-object or list wth additional matrix named $onOff or a matrix with on-off instances 
#' @seealso  \code{\link{testRobustToNAimputation}}
#' @examples
#' dat1 <- list(raw=matrix(round(runif(120), 2), ncol=6))
#' dat1$raw[upper.tri(dat1$raw)] <- NA
#' dat1$raw[7:8,1:3] <- NA
#' dat1$isNA <- is.na(dat1$raw)
#' dat1$setup <- list(grp=rep(LETTERS[4:3], each=3), grpNa =c(LETTERS[4:3]), 
#'   pwGrpNa = matrix(LETTERS[4:3], nrow=1, dimnames=list("D  /  C",c("samp","ref"))) )
#' dat1$means <- matrix(runif(nrow(dat1$raw)*2), ncol=2, dimnames=list(NULL,dat1$setup$grpNa))
#' dat1$filter <- matrix(TRUE, nrow=nrow(dat1$raw), ncol=1)
#' 
#' head(findOnOff(dat1, asFullObj=FALSE), 9)
#' 
#' @export
findOnOff <- function(stat, lim=0.9, compNo=NULL, asFullObj=TRUE, silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## find on/off (and reverse) settings
  ## locate on/off situations for all pairwise comparisons
  ## stat .. ('MArrayLM'-object or list) from the output from \code{testRobustToNAimputation} or \code{moderTestXgrp}
  ## lim ..(numeric, length=1) percent ratio cutoff for considering (all) on or (all) off
  ## addTheo (logical) allows adding $countOnOff to list

  ## move to wrProteo 
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="findOnOff")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) { silent <- FALSE } else { debug <- FALSE }
  if(!"list" %in% mode(stat) || length(stat) <1) stop(fxNa,"Argument 'stat' must be a list or 'MArrayLM'-object from limma")
  argNa <- deparse(substitute(stat))
  if(!("MArrayLM" %in% class(stat)) && !silent) message(fxNa," Caution, '",argNa,"' is not 'MArrayLM'-object as expected")
  if(length(lim) !=1 || !is.finite(lim) || lim >1) {lim <- 0.9; message(fxNa,"Invalid input for 'lim', setting to default (lim=0.9)")}
  if(debug) { message(fxNa,"foo1")}

  ## main extracting
  #if("nNArat" %in% names(stat) && length(stat$setup$pwIndex) !=0) {

  if("isNA" %in% names(stat) && length(stat$setup$pwGrpNa) !=0) {
    if(debug) { message(fxNa,"foo1a")}
      #onOff <- NULL
      #if(!"onOff" %in% names(stat))
      #if(!"nNArat" %in% names(stat)) stat$nNArat <- countOnOff(stat, lim=lim, addTy="table") 
    if(debug) {message(fxNa,"foo1b"); foo1b <- list(stat=stat,lim=lim,asFullObj=asFullObj)}
    if("nNArat" %in% names(stat) && length(stat$nNArat) !=0 && !all(is.na(stat$nNArat))) {
      nNA <- matrix(as.integer(sub("/[0-9]+","", stat$nNArat)), ncol=ncol(stat$nNArat), dimnames=dimnames(stat$nNArat))  # count of NAs
      nSa <- as.integer(sub("[0-9]+/", "", stat$nNArat[1,]))
      if(debug) { message(fxNa,"foo2a")}

    } else {
        #(stat$setup$grpNa, function(x) which(stat$setup$grp==x) )
      nNA <- sapply(stat$setup$grpNa, function(x) rowSums(stat$isNA[, which(stat$setup$grp==x), drop=FALSE]) )
      nSa <- table(stat$setup$grp)       
      stat$nNArat <- nNA/matrix(rep(nSa, each=nrow(nNA)), nrow=nrow(nNA))
      if(debug) { message(fxNa,"foo2b")}
    }
    
    nNA <- nNA[, match( colnames(stat$means), sub("^nNA\\.","", colnames(nNA))), drop=FALSE]                           # in order of means
    if(debug) { message(fxNa,"foo2c"); foo2c <- list(stat=stat, nNA=nNA, nSa=nSa, lim=lim, asFullObj=asFullObj)}
    if(all(is.na(nNA))) stop(fxNa,"Problem extracting number of NA or matching colnames")
    stat$setup$nSa <- nSa[match( colnames(stat$means), names(nSa))]      # in order of stat$means
    stat$nNArat <- stat$nNArat[, match( colnames(stat$means), colnames(stat$nNArat)), drop=FALSE]   # in order of means
    if(debug) {message(fxNa,"foo3"); foo3 <- list(stat=stat,nNA=nNA,nSa=nSa, lim=lim)}

    onOff <- apply(stat$setup$pwGrpNa, 1, function(x, li=lim) { 
      ((stat$nNArat[, grep(paste0("^",x[1],"$"), colnames(stat$nNArat))] >= li) & (stat$nNArat[, grep(paste0("^",x[2],"$"), colnames(stat$nNArat))] <= 1 - li)) |  
      ((stat$nNArat[, grep(paste0("^",x[1],"$"), colnames(stat$nNArat))] <= 1 - li) & (stat$nNArat[, grep(paste0("^",x[2],"$"), colnames(stat$nNArat))] >= li))  })
    if(all(colnames(onOff)==colnames(stat$filter), dim(onOff)==dim(stat$filter))) onOff <- onOff & stat$filter else warning(fxNa,"Can't integrate filtering (no $filter or bad dimensions)")
  }
  if(isTRUE(asFullObj)) { stat$onOff <- onOff; stat } else onOff
}  

