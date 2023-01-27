#' Order Columns in list of matrixes
#'
#' @description
#' This function orders columns in list of matrixes (or matrix) according to argument \code{sampNames}.
#' This function can be used to adjust/correct the order of samples after reading data using \code{readMaxQuantFile()}, \code{readPDExport()} etc.
#' The input may also be MArrayLM-type object from package \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{limma} or from \code{\link{moderTestXgrp}} or \code{\link{moderTest2grp}}.
#'
#' @param dat (matrix, list or MArrayLM-object from limma) main input of which columns should get re-ordered, may be output from \code{\link{moderTestXgrp}} or \code{\link{moderTest2grp}}.
#' @param replNames (character) new column-names (in order as input from \code{dat}), allows renaming colnames before defining new order
#' @param sampNames (character) column-names in desired order for output (must match colnames of \code{dat} or \code{replNames}, if used)
#' @param newNames depreciated, plese use \code{replNames} instead
#' @param useListElem (character) in case \code{dat} is list, all list-elements who's columns should get (re-)ordered
#' @param annotElem (character) name of list-element of \code{dat} with annotation data to get in new order
#' @param silent (logical) suppress messages
#' @param debug (logical) display additional messages for debugging
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @return This function returns an object of same class as input \code{dat}  (ie matrix, list or MArrayLM-object from limma)
#' @seealso \code{\link{moderTestXgrp}} for single comparisons; \code{\link[base]{order}}
#' @examples
#' grp <- factor(rep(LETTERS[c(3,1,4)], c(2,3,3)))
#' dat1 <- matrix(1:15, ncol=5, dimnames=list(NULL,c("D","A","C","E","B")))
#' corColumnOrder(dat1, sampNames=LETTERS[1:5])
#'
#' dat1 <- list(quant=dat1,raw=dat1)
#'   dat1
#' corColumnOrder(dat1, sampNames=LETTERS[1:5])
#' @export
corColumnOrder <- function(dat, replNames=NULL, sampNames, useListElem=c("quant","raw","counts"), annotElem="sampleSetup", newNames=NULL, silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## order columns in list of matrixes (or matrix) according to 'sampNames'
  ## This function can be used to adjust/correct the order of samples after reading data using \code{readMaxQuantFile()}, \code{readPDExport()} etc.
  ## dat (list or matrix) main input of which columns should get re-ordered
  ## sampNames (character) column-names in desired order for output
  ## useListElem (character) all names of list-elements where the reordering should be performed
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="corColumnOrder")
  .corPathW <- function(x) gsub("\\\\", "/", x)
  .corEnum <- function(colNa, repl=c("Samp","samp","Rep","rep","Re","re","R","r","Number","number","No","no","N","n"), sep=c("_","-")) {
    ## function to match to remove enumeration characters in colNa; check for 'abc_Rep123' and correct to 'abc_123'
    ## colNa (character)
    ##
    out1 <- lapply(repl, function(x) {ch2 <- grep(paste0(".\\",sep[1],x,"[[:digit:]]+$"), colNa)
      if(identical(ch2, 1:length(colNa))) {le <- sub(paste0("\\",sep[1],x,"[[:digit:]]+$"), paste0("\\",sep[1]), colNa); paste0(le, substr(colNa, nchar(le)+2, nchar(colNa)))}})
    chLe1 <- sapply(out1, length)
    out2 <- lapply(repl, function(x) {ch2 <- grep(paste0(".\\",sep[2],x,"[[:digit:]]+$"), colNa)
      if(identical(ch2, 1:length(colNa))) {le <- sub(paste0("\\",sep[2],x,"[[:digit:]]+$"), paste0("\\",sep[2]), colNa); paste0(le, substr(colNa, nchar(le)+2, nchar(colNa)))}})
    chLe2 <- sapply(out2, length)
    ch3 <- chLe1[which.max(chLe1)] > chLe2[which.max(chLe2)]
    out <- if(ch3) {if(chLe1[which.max(chLe1)] >0) out1[[which.max(chLe1)]] else colNa} else {if(chLe2[which.max(chLe2)] >0) out1[[which.max(chLe2)]] else colNa}
    out }

  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  datOK <- TRUE
  alreadyOK <- FALSE
  newO <- NA            # initialize
  chSetupNa <- c("groups","level","lev", "sdrfDat", "annotBySoft")
  ## checks
  if(length(dat) <0) { datOK <- FALSE
    msg <- "'dat' is empty, nothing to do" }
  if(datOK & length(names(dat)) <0) { datOK <- FALSE
    msg <- "'dat' has no names, nothing to do" }
  if(datOK) { ch1 <- useListElem %in% names(dat)    # check for useListElem
    if(any(ch1, na.rm=TRUE)) { if(any(!ch1, na.rm=TRUE)) useListElem <- useListElem[which(ch1)]   # update
      } else {datOK <- FALSE; if(!silent) message(fxNa,"list-elements ",wrMisc::pasteC(useListElem),"  not found in 'dat'")}}
  if(debug) {message(fxNa,"cCO1")}

  if(length(useListElem) >0) {
    chEl <- useListElem %in% names(dat)
    if(any(!chEl, na.rm=TRUE)) useListElem <- useListElem[which(chEl)]
  }

  ## main
  if(all(datOK, length(replNames) >0, length(replNames)==length(sampNames))) {
    ## replace colnames (if needed)
    if(debug) {message(fxNa,"replace colnames    cCO1b")}
    if(is.list(dat)) { for(i in useListElem) colnames(dat[[i]]) <- replNames
    } else if(is.matrix(dat)) {
      if(ncol(dat) != length(sampNames)) warning(fxNa,"'dat' has different number of columns as length of 'sampNames' !!  The function might be using the wrong ones !")
      for(i in useListElem) colnames(dat[[i]])[1:length(sampNames)] <- replNames }
    if(debug) {message(fxNa,"cCO1b"); cCO1b <- list(at=dat,sampNames=sampNames,replNames=replNames,useListElem=useListElem,datOK=datOK) }
  }

  if(datOK) {
    if(debug) {message(fxNa,"cCO2"); cCO2 <- list(dat=dat,sampNames=sampNames,replNames=replNames, datOK=datOK,useListElem=useListElem) }
    ## compare sampNames & colnames of $quant
    ## Note : comparing by $sampleSetup$groups won't work well due to repeated levels
    newO <- match(sampNames, colnames(dat[[useListElem[1]]]))
    if(any(is.na(newO))) { if(!silent) message(fxNa,"Colnames of 'dat$quant' differ from 'sampNames', trying to adjust .. ")
      ch3 <- .corEnum(colnames(dat[[useListElem[1]]]))
      if(length(ch3)==length(sampNames)) colnames(dat[[useListElem[1]]]) <- ch3
      newO <- match(sampNames, colnames(dat[[useListElem[1]]]))  # update
    }
    if(any(is.na(newO))) { datOK <- FALSE
      if(!silent) message(fxNa,"Failed : Unable to match ",sum(is.na(newO))," suggested sampNames ( ",wrMisc::pasteC(sampNames[which(is.na(newO))], quoteC="'")," )")
    } else {
      ## finally adjust order of $quant etc based on $quant
      if(debug) {message(fxNa,"cCO2b"); cCO2b <- list(dat=dat,newO=newO,sampNames=sampNames,newNames=newNames,useListElem=useListElem,datOK=datOK) }
      if(identical(newO, 1:length(sampNames))) {
        if(!silent) message(fxNa,"Quant/counting data already in good order ..")
        alreadyOK <- TRUE
      } else {
        if(is.list(dat)) { for(i in wrMisc::naOmit(match(useListElem, names(dat)))) { dat[[i]] <- if(length(dim(dat[[i]])) ==2) {
            if(any(dim(dat[[i]])==1)) matrix(dat[[i]][,newO], nrow=nrow(dat[[i]]), dimnames=dimnames(dat[[i]][,newO])) else dat[[i]][,newO]
          } else { if(length(dim(dat[[i]])) ==3) array(as.numeric(dat[[i]][,newO,]), dim=c(nrow(dat[[i]]), length(newO), dim(dat[[i]])[3]),
            dimnames=list(rownames(dat[[i]]), colnames(dat[[i]])[newO], dimnames(dat[[i]])[[3]])) }
        if(debug) message(fxNa,"Sucessfully adjusted quantitation data to new order") }}
      }
    }
  } else if(!silent) message(fxNa,"Failed to adjust quantitative/count data")
  if(debug) {message(fxNa,"cCO3"); cCO3 <- list(dat=dat,sampNames=sampNames,useListElem=useListElem,annotElem=annotElem,newO=newO)}
  if(datOK & !alreadyOK) {
    ## try adjusting order of $quant, $raw and $counts
    ## Continue adjusting order, now sample annotation
    ## check if $sampleSetup present, => look for filenames in $sampleSetup : $sampleSetup$sdrfDat$comment.file.uri. or $comment.data.file.  OR   $sampleSetup$annotBySoft$File.Name 
    ##  .. this part won't work with PL ??
    setupOK <- length(annotElem) >0
    if(setupOK) { if(length(annotElem) >1) annotElem <- annotElem[1]
      setupOK <- annotElem[1] %in% names(dat)}
    if(setupOK) setupOK <- chSetupNa %in% names(dat[[annotElem]])
    if(any(setupOK, na.rm=TRUE)) {
      ## what to use from sampleSetup for good comparison ?
      if(debug) {message(fxNa,"cCO3b"); cCO3b <- list(dat=dat,newO=newO,sampNames=sampNames,newNames=newNames,useListElem=useListElem,datOK=datOK,annotElem=annotElem) }
      if("annotBySoft" %in% names(dat[[annotElem]])) {
        newO <- apply(dat[[annotElem]][["annotBySoft"]], 2, function(x) match(sampNames, x))
        ch4 <- colSums(is.na(newO))
        if(any(ch4 <1)) { newO <- newO[,which(ch4 <1)[1]]
        } else {                                     # no column with direct matches
          ##
          useCol <- match(c("file.name","raw.file","file"), tolower(colnames(dat[[annotElem]][["annotBySoft"]])))
          useCol <- if(any(!is.na(useCol))) useCol[which(!is.na(useCol))[1]] else NA
          if(!is.na(useCol)) {     ## no drect matches
            ch5 <- .corEnum(sub("\\.RAW$|\\.Raw$|\\.raw$","", basename(.corPathW(dat[[annotElem]][["annotBySoft"]][,useCol]))))
            newO <- match(sampNames, ch5)  # update
            if(any(is.na(newO))) newO <- match(wrMisc::trimRedundText(sampNames,silent=silent,callFrom=fxNa), wrMisc::trimRedundText(ch5,silent=silent, callFrom=fxNa))  # update
        }
      }
      if(debug) {message(fxNa,"cCO3c"); cCO3c <- list(dat=dat,newO=newO,sampNames=sampNames,newNames=newNames,useListElem=useListElem,datOK=datOK) }

      } else {
        if("sdrfDat" %in% names(dat[[annotElem]])) {
          if(debug) message(fxNa," Since $annotBySoft not found, trying to check order based on sdrfDat")
          newO <- apply(dat[[annotElem]][["sdrfDat"]], 2, function(x) match(sampNames, x))
          ch4 <- colSums(!is.na(newO))
          if(any(ch4 <1)) { newO <- newO[,which(ch4 <1)[1]]
          } else {   # no column with direct matches
            ##
            useCol <- match(c("comment.data.file."), colnames(dat[[annotElem]][["sdrfDat"]]))
            useCol <- if(any(is.na(useCol)) | length(useCol) >1) useCol[which(!is.na(useCol))[1]] else NA
            if(!is.na(useCol)) {     ## no drect matches
              ch5 <- .corEnum(sub("\\.RAW$|\\.Raw$|\\.raw$","", basename(.corPathW(dat[[annotElem]][["sdrfDat"]][,useCol]))))
              newO <- match(sampNames, ch5)  # update
                     #match(wrMisc::trimRedundText(sampNames), wrMisc::trimRedundText(ch5))
              if(any(is.na(newO))) newO <- match(wrMisc::trimRedundText(sampNames,silent=silent,callFrom=fxNa), wrMisc::trimRedundText(ch5,silent=silent, callFrom=fxNa))  # update
          } }
      } }
      if(debug) {message(fxNa,"cCO4") }
      if(all(!is.na(newO)) & any(setupOK, na.rm=TRUE)) {
        ## apply new order to sample annotation
        if(identical(newO,1:length(sampNames))) { if(!silent) message(fxNa,"Sample annottaion data alreay in correct order ...")
        } else {
          for(i in 1:length(dat[[annotElem]])) { if(length(dat[[annotElem]][[i]]) >1 & length(dim(dat[[annotElem]][[1]]) >1)) {
              if(ncol(dat[[annotElem]][[i]])==length(newO)) dat[[annotElem]][[i]] <- dat[[annotElem]][[i]][,newO]
            } else if(length(dat[[annotElem]][[i]])==length(newO)) dat[[annotElem]][[i]] <- dat[[annotElem]][[i]][newO]
        } }
      if(debug) message(fxNa,"Sucessfully adjusted sample annotation to new order")
    } else { if(!silent & !alreadyOK) message(fxNa,"Failed to adjust order of sample annotation")} }
  }                    ## end corColumnOrder
  dat }
  
