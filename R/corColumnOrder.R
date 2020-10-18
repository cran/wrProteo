#' Order columns in list of matrixes (or matrix) 
#'
#' This function orders columns in list of matrixes (or matrix) according to argument \code{sampNames}. 
#' This function can be used to adjut/correct the order of samples after reading data using \code{readMaxQuantFile()}, \code{readPDExport()} etc. 
#' The input may also be MArrayLM-type object from package \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{limma} or from \code{\link{moderTestXgrp}} or \code{\link{moderTest2grp}}.
#'
#' @param dat (matrix, list or MArrayLM-object from limma) main input of which columns should get re-ordered, may be output from \code{\link{moderTestXgrp}} or \code{\link{moderTest2grp}}.
#' @param sampNames (character) column-names in desired order for output
#' @param useListElem (character) in case \code{dat} is list, all list-elements who's columns should get (re-)ordered
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @return object of same class as input \code{dat}
#' @seealso \code{\link{moderTestXgrp}} for single comparisons, \code{\link[base]{order}}  
#' @examples
#' grp <- factor(rep(LETTERS[c(3,1,4)], c(2,3,3)))
#' dat1 <- matrix(1:15, ncol=5, dimnames=list(NULL,c("D","A","C","E","B")))
#' corColumnOrder(dat1,sampNames=LETTERS[1:5])  
#' 
#' dat1 <- list(quant=dat1,raw=dat1)
#'   dat1
#' corColumnOrder(dat1,sampNames=LETTERS[1:5]) 
#' @export
corColumnOrder <- function(dat, sampNames, useListElem=c("quant","raw"), silent=FALSE, callFrom=NULL) {
  ## order columns in list of matrixes (or matrix) according to 'sampNames'
  ## This function can be used to adjut/correct the order of samples after reading data using \code{readMaxQuantFile()}, \code{readPDExport()} etc.
  ## dat (list or matrix) main input of which columns should get re-ordered
  ## sampNames (character) column-names in desired order for output
  ## useListElem (character) all names of list-elements where the reordering should be performed
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="corColumnOrder")
  chLE <- useListElem %in% names(dat)
  if(any(chLE) & is.list(dat)) {
    ## this is list of multiple matrixes
    useListElem <- useListElem[which(chLE)]
    newO <- lapply(useListElem, function(x) match(sampNames, colnames(dat[[x]])))
    names(newO) <- useListElem
    if(all(newO[[1]] == 1:ncol(dat[[useListElem[1]]]))) message(fxNa," order already correct !") else {
      if(any(is.na(newO[[1]])) & !silent) message(fxNa," Note : ",sum(is.na(newO[[1]])),
        " (out of ",ncol(dat[[useListElem[1]]]),") column-names not found !  (ignoring)")    
      names(newO) <- useListElem
      for(i in useListElem) dat[[i]] <- dat[[i]][, newO[[i]]] }
  } else {
    ## for simple matrix
    newO <- match(sampNames, colnames(dat)) 
    if(all(newO == 1:ncol(dat))) message(fxNa," order already correct !") else {
      if(any(is.na(newO)) & !silent) message(fxNa," Note : ",sum(is.na(newO))," (out of ",ncol(dat),") column-names not found !  (ignoring)")
      dat <- dat[,wrMisc::naOmit(newO)]  }      
    } 
  dat }
#example
 