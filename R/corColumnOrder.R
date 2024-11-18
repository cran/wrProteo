#' Order Columns In List Of Matrixes, Data.frames And Vectors
#'
#' @description
#' This function orders columns in list of matrixes (or matrix) according to argument \code{sampNames} and also offers an option for changing names of columns.
#' It was (initially) designed to adjust/correct the order of samples after import using \code{readMaxQuantFile()}, \code{readProteomeDiscovererFile()} etc.
#' The input may also be MArrayLM-type object from package \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{limma} or 
#' from functions \code{moderTestXgrp} or \code{moderTest2grp}.
#'
#' @param dat (matrix, list or MArrayLM-object from limma) main input of which columns should get re-ordered, may be output from \code{moderTestXgrp} or \code{moderTest2grp}.
#' @param sampNames (character) column-names in desired order for output (its content must match colnames of \code{dat} or \code{replNames}, if used)
#' @param replNames (character) option for replacing column-names by new/different colnames; should be vector of NEW column-names (in order as input from \code{dat} !), allows renaming colnames before defining new order
#' @param newNames depreciated, pleqse use \code{replNames} instead
#' @param useListElem (character) in case \code{dat} is list, all list-elements who's columns should get (re-)ordered
#' @param annotElem (character) name of list-element of \code{dat} with annotation data to get in new order
#' @param silent (logical) suppress messages
#' @param debug (logical) display additional messages for debugging
#' @param callFrom (character) allows easier tracking of messages produced
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
#' corColumnOrder(dat2, sampNames=LETTERS[1:5], replNames=c("Dd","Aa","Cc","Ee","Bb"))
#' @export
corColumnOrder <- function(dat, sampNames, replNames=NULL, useListElem=c("quant","raw","counts"), annotElem="sampleSetup", newNames=NULL, silent=FALSE, debug=FALSE, callFrom=NULL) {
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
  setupNa <- "sampleSetup"     # name of list-element for sample-setup; so far treated separately to useListElem since may contain vectors and matrix or data.frame
  
  ## checks
  if(debug) {message(fxNa,"cCO0"); cCO0 <- list(dat=dat,sampNames=sampNames,replNames=replNames,useListElem=useListElem,annotElem=annotElem,newNames=newNames,datOK=datOK) }
  if(length(dat) <0) { datOK <- FALSE
    msg <- "'dat' is empty, nothing to do" }
  if(datOK && length(names(dat)) <0) { datOK <- FALSE
    msg <- "'dat' has no names, nothing to do" }
  if(datOK) { 
    if(is.list(dat) && length(dim(dat)) !=2  ) {            # need to add include limma_type obj ?
      ch1 <- useListElem %in% names(dat)    # check for useListElem
      if(any(ch1, na.rm=TRUE)) { if(any(!ch1, na.rm=TRUE)) useListElem <- useListElem[which(ch1)]   # update
        ch2 <- colnames(dat[[useListElem[1]]])
        if(length(ch2) != ncol(dat[[useListElem[1]]])) { datOK <- FALSE; if(!silent) message(fxNa,"Checked list-element ",useListElem[1]," : has no or no differentiating colnames")}
      } else { datOK <- FALSE; if(!silent) message(fxNa,"list-elements ",wrMisc::pasteC(useListElem),"  not found in 'dat'")}
    } else if(length(dim(dat)) ==2) {  # otherwise 'dat' may be matrix or data.frame
      useListElem <- NULL
      if(length(sampNames) <1) { sampNames <- colnames(dat); newNa <- 1:ncol(dat)
      } else newNa <- match(sampNames, colnames(dat))
      if(any(is.na(newNa))) stop(fxNa," ",sum(is.na(newNa))," 'sampNames' do NOT match colnames of dat")
      ## replace names (only in data.frame & matrix)
      if(all(newNa == 1:ncol(dat))) colnames(dat) <- replNames else {        
        if(length(replNames) ==ncol(dat)) colnames(dat) <- replNames }       ## different order
      dat <- dat[,newNa]
      replNames <- NULL       
    } else stop(fxNa, "Unknown type of object as 'dat'") 
  }    
  if(debug) {message(fxNa,"cCO1"); cCO1 <- list(dat=dat,sampNames=sampNames,replNames=replNames,useListElem=useListElem,annotElem=annotElem,newNames=newNames,datOK=datOK) }

  if(datOK && length(useListElem) >0) {
    ## check useListElem  (& remove invalid entries) 
    chEl <- useListElem %in% names(dat) 
    if(any(!chEl, na.rm=TRUE)) useListElem <- useListElem[which(chEl)]
    if(length(useListElem) >0) { chEl <- length(dim(dat[[useListElem[1]]])) >1
      if(any(!chEl)) useListElem <- useListElem[which(chEl)]}
    if(length(useListElem) <1) {warning(fxNa,"After cleaning 'useListElem' nothing remains !?!");  datOK <- FALSE}
  }
  if(debug) message(fxNa,"cCO1b") 

  ## main
  ## extract init colnames
  if(datOK) {
    iniColNa <- if(length(useListElem) >0) colnames(dat[[useListElem[1]]]) else colnames(dat) 
    if(length(iniColNa) != length(sampNames)) warning(fxNa,"'dat' has ",length(iniColNa)," columns while length of 'sampNames'=",length(sampNames)," => DIFFERENT !!  The input seem erroneous !")
    if(length(iniColNa) <1) datOK <- FALSE } else iniColNa <- NULL 

  ## replace colnames (if valid replNames) ..
  if(datOK && length(replNames)==length(iniColNa)) {  # ncol(dat[[useListElem[1]]])
    iniColNa2 <- replNames
    if(length(useListElem) >0) {
      ## replace colnames (if valid 'replNames' given) - replNames must be in same order as  colnames(dat$quant) !!
      if(debug) {message(fxNa,"Replace colnames    cCO1c"); cCO1c <- list(dat=dat,sampNames=sampNames,replNames=replNames,useListElem=useListElem,annotElem=annotElem,newNames=newNames,datOK=datOK)}
      for(i in useListElem) { if(length(dat[[i]]) >0 && length(dim(dat[[i]])) ==2) colnames(dat[[i]]) <- replNames } 
    } else colnames(dat) <- replNames    
  } else iniColNa2 <- NULL 
  if(debug) { message(fxNa,"cCO1d"); cCO1d <- list(dat=dat,sampNames=sampNames,replNames=replNames,iniColNa=iniColNa, useListElem=useListElem,annotElem=annotElem,newNames=newNames,datOK=datOK) }  

  ## find matches for adjusting order
  if(datOK && length(sampNames)==length(iniColNa)) {  # ncol(dat[[useListElem[1]]])
    ## try adjusting order of $quant, $raw and $counts
    if(debug) {message(fxNa,"cCO2"); cCO2 <- list(dat=dat,sampNames=sampNames,replNames=replNames, datOK=datOK,useListElem=useListElem) }
    ## transform to list in case of matrix or data.frame :
    if(length(useListElem) <1 && is.matrix(dat) || is.data.frame(dat)) { is2dim <- TRUE; dat <- list(dat=dat); useListElem <- "dat"} else is2dim <- FALSE

    ## Adjust order of columns in $quant etc: compare sampNames & colnames of $quant
      ## Note : comparing by $sampleSetup$groups won't work well due to repeated levels
    newO <- match(sampNames, iniColNa)      
    if(length(sampNames) != length(iniColNa) && !silent) message(fxNa,"Note : Output assigned to contain fewer columns")

    ## try recuperating names if any not found/matched
    if(any(is.na(newO))) { if(!silent) message(fxNa, if(debug) "cCO2b  ","Colnames of 'dat$quant' differ from 'sampNames', trying to adjust .. ")
      if(length(iniColNa2) >0) {newO2 <- match(sampNames, iniColNa2)
        if(sum(is.na(newO) > sum(is.na(newO2)))) {newO <- newO2; iniColNa <- iniColNa2 }} }
    if(any(is.na(newO))) {       
      ## modify colnames of $abund to remove enumerator-names
      colNa3 <- wrMisc::rmEnumeratorName(iniColNa, sepEnum=c(""," ","-","_"), newSep="_", incl=c("anyCase","trim1"))
       if(debug) cCO2b <- list(newO=newO,colNa3=colNa3,dat=dat,sampNames=sampNames,replNames=replNames,iniColNa=iniColNa, useListElem=useListElem,annotElem=annotElem,newNames=newNames,datOK=datOK)
      newO2 <- match(sampNames, colNa3)         # check
      if(sum(is.na(newO) > sum(is.na(newO2)))) {newO <- newO2} } # update
 
    if(any(is.na(newO))) {
      sampNa2 <- wrMisc::rmEnumeratorName(sampNames, sepEnum=c(""," ","-","_"), newSep="_", incl=c("anyCase","trim1"))
      newO2 <- match(sampNa2, colNa3)           # check, use if improvement in no of matches           
      if(sum(is.na(newO2)) < sum(is.na(newO))) {newO <- newO2; sampNames <- sampNa2} # colnames(dat[[useListElem[1]]]) <- colNa3
    }
    if(debug) { message(fxNa,"cCO2c"); cCO2c <- list(dat=dat,sampNames=sampNames,replNames=replNames,newO=newO, useListElem=useListElem,annotElem=annotElem,newNames=newNames,datOK=datOK) }
    if(sum(!is.na(newO)) >1) {
      ## apply newO
      if(any(is.na(newO)) && !silent) message(fxNa,"Note : ",sum(is.na(newO))," colnames NOT found - will be omitted")
      alreadyOK <- identical(newO, 1:length(sampNames))
      sameNames <- identical(iniColNa, sampNames)
      if(length(newO) <2) {
        if(!silent) message(fxNa,"Quant/counting data already in good order ..")
       } else {
        for(i in wrMisc::naOmit(match(useListElem, names(dat)))) { dat[[i]] <- if(length(dim(dat[[i]])) ==2) {
            if(any(dim(dat[[i]]) ==1)) matrix(dat[[i]][,newO], nrow=nrow(dat[[i]]), dimnames=dimnames(dat[[i]][,newO])) else dat[[i]][,newO] 
          } else { if(length(dim(dat[[i]])) ==3) array(as.numeric(dat[[i]][,newO,]), dim=c(nrow(dat[[i]]), length(newO), dim(dat[[i]])[3]),
            dimnames=list(rownames(dat[[i]]), colnames(dat[[i]])[newO], dimnames(dat[[i]])[[3]])) }
          if(!sameNames && length(dim(dat[[i]])) >1)  colnames(dat[[i]]) <- sampNames  
        }
        ## Continue adjusting order, now in $sampleSetup :  list of vectors, matr and/or df
        ## presume $sampleSetup is in same order as $quant & $raw
        if(length(dat) >1 && setupNa %in% names(dat) && !setupNa %in% useListElem) {   # special treatment of useListElem
          ch1 <- sapply(dat[[setupNa]], function(x) { if(length(dim(x)) >0) c(NA, dim(x)==ncol(dat[[useListElem[1]]])) else c(length(x)==ncol(dat[[useListElem[1]]]), NA,NA) })  
          if(length(dim(ch1)) ==0) ch1 <- as.matrix(ch1)
          liEl2 <- c("sdrfDat", "annotBySof"); liEl3 <- c("level", "iniSdrfOrder","sampleNames","sampleNaSdrf")    # for explicit - not used any more
          if(any(ch1[1,], na.rm=TRUE)) for(i in which(ch1[1,])) dat[[setupNa]][[i]] <- dat[[setupNa]][[i]][newO]
          if(any(ch1[2,], na.rm=TRUE)) for(i in which(ch1[2,])) dat[[setupNa]][[i]] <- dat[[setupNa]][[i]][newO,]
          if(any(ch1[3,], na.rm=TRUE)) for(i in which(ch1[3,])) dat[[setupNa]][[i]] <- dat[[setupNa]][[i]][,newO]
          #if(debug)
        }
        if(debug) message(fxNa,"Successfully adjusted quantitation data to new order")
      }
      if(is2dim) dat <- dat$dat         # return to initial type/level of object    
    } else { datOK <- FALSE
      if(!silent) message(fxNa,"Failed : Unable to match ",sum(is.na(newO))," suggested sampNames ( ",wrMisc::pasteC(sampNames[which(is.na(newO))], quoteC="'")," )")}
   ## end corColumnOrder
  } else {
    if(!silent) message(fxNa,"Failed to adjust quantitative/count data")}
  if(debug) {message(fxNa,"cCO3"); cCO3 <- list(dat=dat,sampNames=sampNames,useListElem=useListElem,annotElem=annotElem,newO=newO,alreadyOK=alreadyOK)}                
  dat } 
      
