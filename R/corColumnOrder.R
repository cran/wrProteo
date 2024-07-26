#' Order Columns In List Of Matrixes And Vectors
#'
#' @description
#' This function orders columns in list of matrixes (or matrix) according to argument \code{sampNames}.
#' It can be used to adjust/correct the order of samples after reading data using \code{readMaxQuantFile()}, \code{readProteomeDiscovererFile()} etc.
#' The input may also be MArrayLM-type object from package \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{limma} or 
#' from functions \code{moderTestXgrp} or \code{moderTest2grp}.
#'
#' @param dat (matrix, list or MArrayLM-object from limma) main input of which columns should get re-ordered, may be output from \code{moderTestXgrp} or \code{moderTest2grp}.
#' @param replNames (character) new column-names (in order as input from \code{dat}), allows renaming colnames before defining new order
#' @param sampNames (character) column-names in desired order for output (must match colnames of \code{dat} or \code{replNames}, if used)
#' @param newNames depreciated, plese use \code{replNames} instead
#' @param useListElem (character) in case \code{dat} is list, all list-elements who's columns should get (re-)ordered
#' @param annotElem (character) name of list-element of \code{dat} with annotation data to get in new order
#' @param silent (logical) suppress messages
#' @param debug (logical) display additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns an object of same class as input \code{dat}  (ie matrix, list or MArrayLM-object from limma)
#' @seealso  \code{\link{readMaxQuantFile}}, \code{\link{readProteomeDiscovererFile}}; \code{\link[wrMisc]{moderTestXgrp}} or \code{\link[wrMisc]{moderTest2grp}}
#' @examples
#' grp <- factor(rep(LETTERS[c(3,1,4)], c(2,3,3)))
#' dat1 <- matrix(1:15, ncol=5, dimnames=list(NULL,c("D","A","C","E","B")))
#' corColumnOrder(dat1, sampNames=LETTERS[1:5])
#'
#' dat2 <- list(quant=dat1, raw=dat1)
#' dat2
#' corColumnOrder(dat2, sampNames=LETTERS[1:5])
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
      if(identical(ch2, 1:length(colNa))) {le <- sub(paste0("\\",sep[1],x,"[[:digit:]]+$"), paste0("\\",sep[1]), colNa); paste0(le, substr(colNa, nchar(le)+2, nchar(colNa)))}}) # nolint # nolint: line_length_linter.
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
  if(datOK && length(names(dat)) <0) { datOK <- FALSE
    msg <- "'dat' has no names, nothing to do" }
  if(datOK) { ch1 <- useListElem %in% names(dat)    # check for useListElem
    if(any(ch1, na.rm=TRUE)) { if(any(!ch1, na.rm=TRUE)) useListElem <- useListElem[which(ch1)]   # update
      ch2 <- colnames(dat[[useListElem[1]]])
      if(length(ch2) != ncol(dat[[useListElem[1]]])) { datOK <- FALSE; if(!silent) message(fxNa,"Checked list-element ",useListElem[1]," : has no or no differentiating colnames")}
    } else { datOK <- FALSE; if(!silent) message(fxNa,"list-elements ",wrMisc::pasteC(useListElem),"  not found in 'dat'")}}
  if(debug) {message(fxNa,"cCO1"); cCO1 <- list(dat=dat,sampNames=sampNames,replNames=replNames,useListElem=useListElem,annotElem=annotElem,newNames=newNames,datOK=datOK) }

  if(length(useListElem) >0) {
    chEl <- useListElem %in% names(dat)
    if(any(!chEl, na.rm=TRUE)) useListElem <- useListElem[which(chEl)]
    if(length(useListElem) <1) {warning(fxNa,"After cleaning 'useListElem' nothing remains !?!");  datOK <- FALSE}
  }

  ## main
  if(datOK && length(replNames) >0 && length(replNames)==length(sampNames)) {
    ## replace colnames (if valid 'replNames' given) - replNames must be in same order as  colnames(dat$quant) !!
    if(debug) {message(fxNa,"Replace colnames    cCO1b")}
    if(is.list(dat)) { for(i in useListElem) colnames(dat[[i]]) <- replNames
    } else if(is.matrix(dat)) {
      if(ncol(dat) != length(sampNames)) warning(fxNa,"'dat' has different number of columns as length of 'sampNames' !!  The function might be using the wrong ones !")
      for(i in useListElem) colnames(dat[[i]])[1:length(sampNames)] <- replNames }
    if(debug) { message(fxNa,"cCO1b"); cCO1b <- list() }
  }

  if(datOK) {
    ## try adjusting order of $quant, $raw and $counts
    if(debug) {message(fxNa,"cCO2"); cCO2 <- list(dat=dat,sampNames=sampNames,replNames=replNames, datOK=datOK,useListElem=useListElem) }
    ## Adjust order of columns in $quant etc: compare sampNames & colnames of $quant
    if(length(sampNames) ==ncol(dat[[useListElem[1]]])) {
      ## Note : comparing by $sampleSetup$groups won't work well due to repeated levels
      newO <- match(sampNames, colnames(dat[[useListElem[1]]]))
      if(any(is.na(newO))) { if(!silent) message(fxNa,"Colnames of 'dat$quant' differ from 'sampNames', trying to adjust .. ")
        ## modify colnames of $abund to remove enumerator-names
        ch3 <- wrMisc::rmEnumeratorName(colnames(dat[[useListElem[1]]]), sepEnum=c(""," ","-","_"), newSep="_", incl=c("anyCase","trim1"))
        if(length(ch3) ==length(sampNames)) colnames(dat[[useListElem[1]]]) <- ch3
        ## same treatment to sampleNames  for higher chances of matching
        sampNa2 <- wrMisc::rmEnumeratorName(sampNames, sepEnum=c(""," ","-","_"), newSep="_", incl=c("anyCase","trim1"))
        newO <- match(sampNa2, ch3)  # update
        ## this could be made after trimming if still not successful, see also wrMisc::matchMatrixLinesToRef()
      }
      if(any(is.na(newO))) { datOK <- FALSE
        if(!silent) message(fxNa,"Failed : Unable to match ",sum(is.na(newO))," suggested sampNames ( ",wrMisc::pasteC(sampNames[which(is.na(newO))], quoteC="'")," )")
      } else {
        ## finally adjust order of $quant etc based on $quant
        if(debug) {message(fxNa,"cCO2b"); cCO2b <- list() }
        if(identical(newO, 1:length(sampNames))) {
          if(!silent) message(fxNa,"Quant/counting data already in good order ..")
          alreadyOK <- TRUE
        } else {
          if(is.list(dat)) { for(i in wrMisc::naOmit(match(useListElem, names(dat)))) { dat[[i]] <- if(length(dim(dat[[i]])) ==2) {
              if(any(dim(dat[[i]]) ==1)) matrix(dat[[i]][,newO], nrow=nrow(dat[[i]]), dimnames=dimnames(dat[[i]][,newO])) else dat[[i]][,newO]
            } else { if(length(dim(dat[[i]])) ==3) array(as.numeric(dat[[i]][,newO,]), dim=c(nrow(dat[[i]]), length(newO), dim(dat[[i]])[3]),
              dimnames=list(rownames(dat[[i]]), colnames(dat[[i]])[newO], dimnames(dat[[i]])[[3]])) }
          if(debug) message(fxNa,"Sucessfully adjusted quantitation data to new order") }}
        }
      }
    } else { newO <- NA; datOK <- FALSE
      if(!silent) message(fxNa,"Failed to adjust quantitative/count data")}

    }

  if(debug) {message(fxNa,"cCO3"); cCO3 <- list(dat=dat,sampNames=sampNames,useListElem=useListElem,annotElem=annotElem,newO=newO,alreadyOK=alreadyOK)}
  if(datOK && !alreadyOK && length(dat[[annotElem]]) >0) {
    ## Continue adjusting order, now in $sampleSetup
    ## presume $sampleSetup is in same order as $quant & $raw
    ch1 <- sapply(dat[[annotElem]], function(x) if(length(dim(x)) >0) c(NA, dim(x)==length(sampNames)) else c(length(x)==length(sampNames),NA,NA))
    if(length(dim(ch1)) <2) ch1 <- matrix(ch1, ncol=length(dat[[annotElem]]), dimnames=list(NULL,names(dat[[annotElem]])))
    if(any(ch1[1,], na.rm=TRUE)) for(i in which(ch1[1,])) dat[[annotElem]][[i]] <- dat[[annotElem]][[i]][newO]
    if(any(ch1[2,], na.rm=TRUE)) for(i in which(ch1[2,])) dat[[annotElem]][[i]] <- dat[[annotElem]][[i]][newO,]
    if(any(ch1[3,], na.rm=TRUE)) for(i in which(ch1[3,])) dat[[annotElem]][[i]] <- dat[[annotElem]][[i]][,newO]

  }                    ## end corColumnOrder
  dat } 
    
