#' Extract results from moderated t-tests
#'
#' This function allows convenient access to results produced using the functions \code{\link{moderTest2grp}} or \code{moderTestXgrp}.
#' The user can define the threshold which type of multiple testing correction should be used
#'  (as long as the  multiple testing correction method was actually performed as part of testing).
#'  
#' @param stat ('MArrayLM'-object or list) Designed for the output from \code{moderTest2grp} or \code{moderTestXgrp}
#' @param compNo (integer) the comparison number/index to be used 
#' @param statTy (character) the multiple-testing correction type to be considered when looking for significant changes  with threshold \code{thrsh} (depends on which have been run initially with \code{moderTest2grp} or \code{moderTestXgrp})
#' @param thrsh (numeric) the threshold to be applied on \code{statTy} for the result of the statistcal testing (after multiple testing correction)
#' @param FCthrs (numeric) Fold-Change threshold given as Fold-change and NOT log2(FC), default at 1.5 (for filtering at M-value =0.585)
#' @param annotCol (character) column-names from the annotation to be included
#' @param nSign (integer) number of significant digits whe returning results
#' @param addTy (character) additional groups to add (so far only "allMeans" available) in addition to the means used in the pairwise comparison
#' @param filename (character) optional (path and) file-name for exporting results to csv-file
#' @param fileTy (character) file-type to be used with argument \code{filename}, may be 'csvEur' or 'csvUS'
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @return limma-type MA-object (list)
#' @seealso \code{\link{moderTest2grp}} for single comparisons, \code{\link{moderTestXgrp}} for multiple comparisons, \code{\link[limma]{lmFit}} and the \code{eBayes}-family of functions in package \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{limma}
#' @examples
#' grp <- factor(rep(LETTERS[c(3,1,4)],c(2,3,3)))
#' set.seed(2017); t8 <- matrix(round(rnorm(208*8,10,0.4),2), ncol=8,
#'   dimnames=list(paste(letters[],rep(1:8,each=26),sep=""), paste(grp,c(1:2,1:3,1:3),sep="")))
#' t8[3:6,1:2] <- t8[3:6,1:2] +3                    # augment lines 3:6 (c-f) 
#' t8[5:8,c(1:2,6:8)] <- t8[5:8,c(1:2,6:8)] -1.5    # lower lines 
#' t8[6:7,3:5] <- t8[6:7,3:5] +2.2                  # augment lines 
#' ## expect to find C/A in c,d,g, (h)
#' ## expect to find C/D in c,d,e,f
#' ## expect to find A/D in f,g,(h) 
#' library(wrMisc)     # for testing we'll use this package
#' test8 <- moderTestXgrp(t8, grp) 
#' extractTestingResults(test8)
#' @export
extractTestingResults <- function(stat, compNo=1, statTy="BH",thrsh=0.05, FCthrs=1.5, annotCol=c("Accession","EntryName","GeneName"), 
  nSign=6, addTy=c("allMeans"), filename=NULL, fileTy="csvUS", silent=FALSE, callFrom=NULL) {
  ##
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="extractTestingResults")
  argNa <- deparse(substitute(stat))
  if(!isTRUE(silent)) silent <- FALSE
  if(!"list" %in% mode(stat) | length(stat) <1) stop("'stat' must be a list or 'MArrayLM'-object from limma")
  if(!("MArrayLM" %in% class(stat)) & !silent) message(fxNa," Caution, '",argNa,"' is not 'MArrayLM'-object as expected")
  if(length(statTy) <1) { statTy <- "BH"
    if(!silent) message(fxNa," argument 'statTy' empty, setting to default 'BH'")}
  useFdrTy <- if(identical(statTy,"BH") & "FDR" %in% names(stat)) "FDR" else statTy    # maybe not needed to force to FDR
  chLstEl <- c(useFdrTy,"annot") %in% names(stat) 
  if(!chLstEl[1]) stop("Cannot find list-element '",useFdrTy,"' in 'stat'")
  if(!chLstEl[2]) {
    if(!silent & length(annotCol) >0) message(fxNa,"Cannot find list-element 'annot' in ",argNa)
    annotCol <- NULL }
  if(length(compNo) != 1) stop("'compNo' must be numeric and of length=1")
  if(compNo > ncol(stat[[useFdrTy]])) { compNo <- 1
    message(fxNa," Invalid entry of 'compNo', setting to defaul compNo=1")}   
  if(is.na(FCthrs) | !is.numeric(FCthrs)) FCthrs <- NULL
  if(length(FCthrs) >0 & is.numeric(FCthrs)) FCthrs <- log2(FCthrs) else FCthrs <- NULL
  groupSep <- "-"                                    # used to separate comparison groups  
  ## main extracting
  ## redo sample-pair assoc
  avCol <- wrMisc::sampNoDeMArrayLM(stat, compNo, lstP=useFdrTy)          # ultimately switch to function in wrMisc

  ## filtering ?  normally already taken care of during testing
  
  ## need FC values
  logFC <- stat$means[,avCol[2]] - stat$means[,avCol[1]]
  fcOk <- if(length(FCthrs) ==1) abs(logFC) > FCthrs else rep(TRUE, length(logFC))
  chNa <- is.na(fcOk)
  if(any(chNa)) fcOk[which(chNa)] <- FALSE
  ## check for Fdr results
  fdrOk <- stat[[useFdrTy]][,compNo] < thrsh
  chNa <- is.na(fdrOk)
  if(any(chNa)) fdrOk[which(chNa)] <- FALSE
  if(any(fcOk & fdrOk)) {
    extrLi <- if(TRUE) which(fcOk & fdrOk) else which(fcOk | fdrOk)
    useCompNo <- c(compNo,(1:ncol(stat[[useFdrTy]]))[-compNo])
    ## some results .. continue
    if(any(c("all","allFDR") %in% addTy) & length(avCol) >2) {
      out <- stat[[useFdrTy]][extrLi,]
      colnames(out) <- paste0(useFdrTy,".",colnames(out))
      out <- if(nrow(out) >1) out[,useCompNo] else matrix(out[,useCompNo], nrow=1, dimnames=list(names(extrLi),colnames(out)[useCompNo]))   # place comparison of interest first      
    } else out <- data.frame(FDR=stat[[useFdrTy]][extrLi,compNo])  
    
    ## prepare FC
    if(any(c("all","allFC") %in% addTy) & length(avCol) >2) {
      outX <- sapply(useCompNo, function(x) wrMisc::sampNoDeMArrayLM(stat, x, lstP=useFdrTy))
      out2 <- (stat$means[,outX[2,]] - stat$means[,avCol[1,]])[extrLi,useCompNo]
      colnames(out2) <- if(length(useCompNo) >1) paste0("logFC.",apply(outX,2, function(x) paste(colnames(stat$means)[x],collapse="-"))) else "logFC"
    } else {out2 <- as.matrix(logFC[extrLi]); colnames(out2) <- paste0("logFC.",colnames(stat[[useFdrTy]])[compNo]) }
    ## prepare means
    ch1 <- any(c("all","allMeans") %in% addTy) | length(avCol) <3
    out3 <- if(ch1) stat$means[extrLi,] else stat$means[extrLi,avCol]
    if(length(extrLi)==1) out3 <- matrix(out3, nrow=1, dimnames=list(names(extrLi),  if(ch1) colnames(stat$means) else colnames(stat$means)[avCol] ))
    if(length(dim(out3)) >1) colnames(out3) <- paste0("av.",colnames(out3))
    out <- signif(cbind(out, out2, out3), nSign)
    if(length(annotCol) >0) out <- cbind(stat$annot[extrLi,wrMisc::naOmit(match(annotCol,colnames(stat$annot)))], out)  # add annotation

    ## optional writing to file
    if(length(filename)==1) {
      digits <- min(nSign, 12)
      tmp <- if(identical(fileTy,"csvEur")) {
        try(utils::write.csv2(as.matrix(format(out, digits=digits)), filename, row.names=FALSE,quote=FALSE),silent=silent)
      } else try(utils::write.csv(as.matrix(format(out, digits=digits)), filename, row.names=FALSE,quote=FALSE),silent=silent)
      if("try-error" %in% class(tmp)) message(fxNa," Note: Did not manage to write results to file '",filename,"', check for rights to write  ...") else {
        if(!silent) message(fxNa," Wrote results successfully to file '",filename,"'")}
    }
    out
  } else {if(!silent) message(fxNa,"No results pass thresholds"); return(NULL)} }
  
#' @export
.extractTestingResultsOld <- function(stat, compNo=1, statTy="BH",thrsh=0.05, FCthrs=1.5, annotCol=c("Accession","EntryName","GeneName"), nSign=6,add="allMeans",filename=NULL,fileTy="csvUS",silent=FALSE,callFrom=NULL) {
  ##
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="extractTestingResults")
  argNa <- deparse(substitute(stat))
  if(!"list" %in% mode(stat) | length(stat) <1) stop("'stat' must be a list or 'MArrayLM'-object from limma")
  if(!("MArrayLM" %in% class(stat)) & !silent) message(fxNa," Caution, '",argNa,"' is not 'MArrayLM'-object as expected")
  if(length(statTy) <1) { statTy <- "BH"
    if(!silent) message(fxNa," argument 'statTy' empty, setting to default 'BH'")}
  useFdrTy <- if(identical(statTy,"BH") & "FDR" %in% names(stat)) "FDR" else statTy    # maybe not needed to force to FDR
  chLstEl <- c(useFdrTy,"annot") %in% names(stat) 
  if(!chLstEl[1]) stop("Cannot find list-element '",useFdrTy,"' in 'stat'")
  if(!chLstEl[2]) {
    if(!silent & length(annotCol) >0) message(fxNa,"Cannot find list-element 'annot' in ",argNa)
    annotCol <- NULL }
  if(length(compNo) !=1) stop("'compNo' must be numeric and of length=1")
  if(compNo > ncol(stat[[useFdrTy]])) { compNo <- 1
    message(fxNa," Invalid entry of 'compNo', setting to defaul compNo=1")}   
  if(is.na(FCthrs)) FCthrs <- NULL
  if(length(FCthrs) >0 & is.numeric(FCthrs)) FCthrs <- log2(FCthrs) else FCthrs <- NULL
  groupSep <- "-"                                    # used to separate comparison groups  
  ## main extracting
  ## redo sample-pair assoc
  avCol <- wrMisc::sampNoDeMArrayLM(stat, compNo, lstP=useFdrTy)          # ultimately switch to function in wrMisc
  
  ## filtering ?
  
  ## extract means
  sig1 <- order(stat[[useFdrTy]][,compNo])[1: sum(stat[[useFdrTy]][,compNo] < thrsh, na.rm=TRUE)]
  out <- data.frame(FDR=signif(stat[[useFdrTy]][sig1,compNo],nSign)) 
  colnames(out) <- statTy
  if(length(annotCol) >0) out <- data.frame(stat$annot[sig1,annotCol], out) 
  if("means" %in% names(stat)) {
    logFC <- stat$means[sig1,avCol[2]] - stat$means[sig1,avCol[1]]
    out <- cbind(out, signif(logFC,nSign), signif(stat$means[sig1,],nSign) )
    colnames(out)[(length(annotCol) +2) :ncol(out)] <- c("logFC", paste0(colnames(stat$means),".av")) }
    if(length(FCthrs) ==1) { chFC <- abs(out[,"logFC"]) >= FCthrs
      if(all(!chFC)) message(fxNa," NOTE :  NONE of the ",nrow(out)," proteins pass FCthrs !") 
      out <- out[chFC,]   
    } else chFC <- rep(TRUE, nrow(out))
  if("allMeans" %in% add) {                           # add all grp-means used in pairwise test
    avD <- signif(stat$means[sig1,avCol], nSign)[chFC,]
    colnames(avD) <- paste0("av.",colnames(avD))
    out <- cbind(out[,1:(length(annotCol) +2)], avD)
  } else {                                            # stay with grp-means used in pairwise test only, simple colnames
    colnames(out)[length(annotCol) +(1:2)] <- paste(c(statTy,"logFC"),colnames(stat[[useFdrTy]])[compNo],sep=".")}
  ## optional writing to file
  if(length(filename)==1) {
    digits <- min(nSign, 12)
    tmp <- if(identical(fileTy,"csvEur")) {
     try(utils::write.csv2(as.matrix(format(if(length(stat)>0) out else stat, digits=digits)), filename, row.names=FALSE,quote=FALSE),silent=silent)
    } else try(utils::write.csv(as.matrix(format(out, digits=digits)), filename, row.names=FALSE,quote=FALSE),silent=silent)
    if("try-error" %in% class(tmp)) message(fxNa," Note: Did not manage to write results to file '",filename,"', check for rights to write  ...") else {
      if(!silent) message(fxNa," Wrote results successfully to file '",filename,"'")}
  }
  out
}
   
#' @export
.sampNoDeMArrayLM3 <- function(MArrayObj, useComp, groupSep="-",lstMeans="means",lstP="BH",silent=FALSE,callFrom=NULL) {
  ## old version, not used any more
  ## locate sample index from index or name of pair-wise comparisons in list or MArrayLM-object
  fxNa <- wrMisc::.composeCallName(callFrom, newNa=".sampNoDeMArrayLM")
  errMsg <- c("argument 'MArrayObj' is ","empty","doesn't contain the list-element needed  ('",lstMeans,"') !")
  if(length(MArrayObj) <1) stop(errMsg[1:2])
  if(length(MArrayObj[[lstMeans]]) <1) stop(errMsg[-2])
  if(length(colnames(MArrayObj[[lstMeans]])) <1)  stop(" problem with 'MArrayObj$lstMeans' (does not contain matrix of means)")
  if(ncol(MArrayObj[[lstMeans]])==2) {             # only 2 mean-values, no other choice, don't need to try matching anything
    if(!identical(as.character(useComp),"1") & !silent) message(fxNa,"Only 2 columns of mean-values available, can't interpret properly 'useComp=",useComp,"'")
    out <- 1:2
  } else {
    if(length(lstP) <0) stop(" 'lstP' is empty !")
    if(length(colnames(MArrayObj[[lstP]])) <1)  stop(" problem with 'MArrayObj' (does not contain matrix of p-values)")
    ## convert/locate names to index 
    if(is.character(useComp) & length(grep("[[:alpha:]]",useComp)) >0) useComp <- wrMisc::naOmit(match(useComp, MArrayObj[[lstP]] ))
    if(length(useComp) <1) stop("argument 'useComp' is empty or can't locate in comparison-names")
    ## main
    out <- if(ncol(MArrayObj[[lstP]])==2 & colnames(MArrayObj[[lstP]])[1] =="(Intercept)") 1:2 else {
      wrMisc::matchSampToPairw(grpNa=colnames(MArrayObj[[lstMeans]]), pairwNa=colnames(MArrayObj[[lstP]])[useComp], sep=groupSep,silent=silent,callFrom=fxNa)}
    if(length(useComp)==1) out <- as.integer(out)}
  out }
    
