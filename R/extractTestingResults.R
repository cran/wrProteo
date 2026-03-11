#' Extract Results From Moderated t-tests
#'
#' This function allows convenient access to results produced using the functions \code{\link[wrMisc]{moderTest2grp}} or \code{moderTestXgrp}.
#' The user can define the threshold which type of multiple testing correction should be used
#'  (as long as the  multiple testing correction method was actually performed as part of testing).
#'  
#' @details
#' The number of columns returned can be modulated using the argument \code{addTy}. 
#' \code{addTy='slim'} gives the fewest columns, only FDR and FC of current comparison specified in argument \code{compNo} will be returned.
#' \code{addTy='all'} results in a rather exaustive collection of columns (incluing annotation data, if present), while default 'allMeans' gives a compromise
#' 
#' 
#' @param stat ('MArrayLM'-object or list) Designed for the output from \code{moderTest2grp} or \code{moderTestXgrp}
#' @param compNo (integer) the comparison number/index to be used 
#' @param statTy (character) the multiple-testing correction type to be considered when looking for significant changes  with threshold \code{thrsh} (depends on which have been run initially with \code{moderTest2grp} or \code{moderTestXgrp})
#' @param thrsh (numeric) the threshold to be applied on \code{statTy} for the result of the statistcal testing (after multiple testing correction)
#' @param FCthrs (numeric) Fold-Change threshold given as Fold-change and NOT log2(FC), default at 1.5 (for filtering at M-value =0.585)
#' @param annotCol (character) column-names from the annotation to be included
#' @param nSign (integer) number of significant digits whe returning results
#' @param addTy (character) additional groups of columns to add ("allMeans", "all" or "slim") in addition;
#'   if \code{addTy="slim"} only the FDR- and FC-columns directly concerned will be extracted/displayed
#' @param filename (character) optional (path and) file-name for exporting results to csv-file
#' @param fileTy (character) file-type to be used with argument \code{filename}, may be 'csvEur' or 'csvUS'
#' @param silent (logical) suppress messages
#' @param debug (logical) display additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a data.frame or matrix (if no annotation added) with values (and annotation) conform to fiter criteria
#' @seealso \code{\link{testRobustToNAimputation}}, \code{\link[wrMisc]{moderTestXgrp}} or \code{\link[wrMisc]{moderTest2grp}}
#' @examples
#' grp <- factor(rep(LETTERS[4:2],c(2,3,3)))
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
#' extractTestingResults(test8, compNo=c("B-C"))
#' 
#' extractTestingResults(test8, addTy="all")
#' extractTestingResults(test8, addTy="slim")
#' @export
extractTestingResults <- function(stat, compNo=1, statTy="BH", thrsh=0.05, FCthrs=1.5, annotCol=c("Accession","EntryName","GeneName"), 
  nSign=6, addTy=c("allMeans"), filename=NULL, fileTy="csvUS", silent=FALSE, debug=FALSE, callFrom=NULL) {
  ##
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="extractTestingResults")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  ouy <- outX <- NULL

  argNa <- deparse(substitute(stat))
  if(!isTRUE(silent)) silent <- FALSE
  if(!"list" %in% mode(stat) || length(stat) <1) stop("Argument 'stat' must be a list or 'MArrayLM'-object from limma")
  if(!("MArrayLM" %in% class(stat)) && !silent) message(fxNa," Caution, '",argNa,"' is not 'MArrayLM'-object as expected")
  if(length(statTy) <1) { statTy <- "BH"
    if(!silent) message(fxNa,"Argument 'statTy' empty, setting to default 'BH'")}
  useFdrTy <- if(identical(statTy,"BH") & "FDR" %in% names(stat)) "FDR" else statTy    # maybe not needed to force to FDR
  chLstEl <- c(useFdrTy,"annot") %in% names(stat) 
  if(!chLstEl[1]) stop(fxNa,"Cannot find list-element '",useFdrTy,"' in 'stat'")
  if(!chLstEl[2]) {
    if(!silent && length(annotCol) >0) message(fxNa,"Cannot find list-element 'annot' in ",argNa)
    annotCol <- NULL }
  if(length(compNo) != 1) stop("Argument 'compNo' must be numeric and of length=1")
  if(compNo > ncol(stat[[useFdrTy]])) { compNo <- 1
    message(fxNa,"Invalid entry of argument 'compNo', setting to defaul compNo=1")}   
  if(is.na(FCthrs) || !is.numeric(FCthrs)) FCthrs <- NULL
  if(length(FCthrs) >0 && is.numeric(FCthrs)) FCthrs <- log2(FCthrs) else FCthrs <- NULL
  if(length(addTy) ==0 || any(is.na(addTy))) addTy <- "slim"
  groupSep <- "-"                                    # used to separate comparison groups  
  ## main extracting
  ## redo sample-pair assoc
  avCol <- wrMisc::sampNoDeMArrayLM(stat, compNo, outTy=useFdrTy)          # 
  if(debug) {message(fxNa,"eTR0"); eTR0 <- list(stat=stat,compNo=compNo,statTy=statTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,addTy=addTy,filename=filename, avCol=avCol,groupSep=groupSep)}      

  ## filtering ?  normally already taken care of during testing
  
  ## need FC values
  if(!"means" %in% names(stat)) stop(fxNa,"Need list-elment stat$means - but is ABSENT !") else {
    if(!"ind" %in% names(avCol)) stop(fxNa,"Need list-elment avCol$ind - but is ABSENT !")}
  logFC <- stat$means[,avCol$ind[1]] - stat$means[,avCol$ind[2]]
  fcOk <- if(length(FCthrs) ==1) abs(logFC) > FCthrs else rep(TRUE, length(logFC))
  chNa <- is.na(fcOk)
  if(any(chNa)) fcOk[which(chNa)] <- FALSE
  ## check for Fdr results
  fdrOk <- stat[[useFdrTy]][,compNo] < thrsh
  chNa <- is.na(fdrOk)
  if(any(chNa)) fdrOk[which(chNa)] <- FALSE
  if(debug) {message(fxNa,"eTR0a"); eTR0a <- list(stat=stat,compNo=compNo,statTy=statTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,addTy=addTy,filename=filename, avCol=avCol,groupSep=groupSep,useFdrTy=useFdrTy,fcOk=fcOk,fdrOk=fdrOk)}      
  if(any(fcOk & fdrOk)) {
    extrLi <- if(TRUE) which(fcOk & fdrOk) else which(fcOk | fdrOk)
    useCompNo <- if("slim" %in% addTy) compNo else c(compNo,(1:ncol(stat[[useFdrTy]]))[-compNo])   # index of targert comparison and remaining comparisons
    if(debug) {message(fxNa,"eTR0b"); eTR0b <- list(stat=stat,compNo=compNo,statTy=statTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,addTy=addTy,filename=filename, avCol=avCol,groupSep=groupSep,useFdrTy=useFdrTy,fcOk=fcOk,fdrOk=fdrOk)}      
    ## some results .. continue
    if(any(c("all","allFDR","slim") %in% addTy) && ncol(stat[[useFdrTy]]) >1) {    # possibility of extracting multiple FDR
      out <- stat[[useFdrTy]][extrLi,, drop=FALSE]
      out <- out[,if("all" %in% addTy)  unique(c(compNo, 1:ncol(out))) else useCompNo, drop=FALSE]            # place comparison of interest first (if many retained)
      if(debug) {message(fxNa,"eTR0c"); eTR0c <- list(stat=stat,compNo=compNo,statTy=statTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,addTy=addTy,filename=filename, avCol=avCol,logFC=logFC,extrLi=extrLi,useCompNo=useCompNo,useFdrTy=useFdrTy,chNa=chNa,out=out)}  
    } else out <- data.frame(FDR=stat[[useFdrTy]][extrLi, compNo, drop=FALSE])
    colnames(out) <- paste0(useFdrTy,".", colnames(out))
    if(debug) {message(fxNa,"eTR1"); eTR1 <- list(stat=stat,compNo=compNo,statTy=statTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,addTy=addTy,filename=filename, avCol=avCol,logFC=logFC,extrLi=extrLi,useCompNo=useCompNo,useFdrTy=useFdrTy,chNa=chNa,out=out)}  
    
    
    if(any(c("all","allFC","slim") %in% addTy)) {     #  && length(avCol) >2
      message(fxNa, "Functionality 'addTy=' disabled temporariliy, please check for package updates (also for wrMisc)")
      addTy <- c("allMeans")
    }  


    ## prepare FC
    if(any(c("all","allFC","slim") %in% addTy)) {     #  && length(avCol) >2
      ## index of grp-means cited in comparison
      #curComp <- colnames(stat$lods)[compNo]
      outX <- sapply(useCompNo, function(x) wrMisc::sampNoDeMArrayLM(stat, x, outT=useFdrTy))        # need to finish !!
      out2 <- (stat$means[,outX[[1]][2,], drop=FALSE] - stat$means[,outX[[1]][1,], drop=FALSE]) [extrLi, useCompNo, drop=FALSE]    ## creat M-values (FC)      
      colnames(out2) <- if(length(useCompNo) >1) paste0("logFC.",apply(outX,2, function(x) paste(colnames(stat$means)[x],collapse="-"))) else "logFC"
    } else {out2 <- as.matrix(logFC[extrLi]); colnames(out2) <- paste0("logFC.",colnames(stat[[useFdrTy]])[compNo]) }
    if(debug) {message(fxNa,"eTR1b"); eTR1b <- list(stat=stat,compNo=compNo,statTy=statTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,addTy=addTy,filename=filename, avCol=avCol,logFC=logFC,extrLi=extrLi,useFdrTy=useFdrTy,chNa=chNa,out=out,out2=out2,outX=outX,avCol=avCol)}  

    ## prepare means
    ch2 <- any(c("slim") %in% addTy)
    ch1 <- any(c("all","allMeans") %in% addTy)  && !ch2  # || length(avCol) <3)
    out3 <- if(ch1) stat$means[extrLi,] else stat$means[extrLi,avCol, drop=FALSE]
    if(length(dim(out3)) >1) colnames(out3) <- paste0("av.",colnames(out3))
    out <- if("annot" %in% names(stat)) signif(data.frame(out, out2, out3), nSign) else signif(cbind(out, out2, out3), nSign)                 # combine FDR, logFC & means

    ## adjust annotCol accoring to addTy
    if("all" %in% addTy && "annot" %in% names(stat) && "ProteinName" %in% colnames(stat$annot)) annotCol <- unique(c(annotCol, "EntryName","ProteinName"))
    outAnn <- stat$annot[extrLi, wrMisc::naOmit(match(annotCol,colnames(stat$annot))), drop=FALSE]
    ## trim to max char size if slim ?
    
    
    if(length(annotCol) >0) out <- cbind(outAnn, out)  # add annotation
    if(debug) {message(fxNa,"eTR2"); eTR2 <- list(stat=stat,compNo=compNo,statTy=statTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,addTy=addTy,filename=filename, logFC=logFC,chNa=chNa,out=out,useFdrTy=useFdrTy)}  
    
    ## optional writing to file
    if(length(filename)==1) {
      digits <- min(nSign, 12)
      tmp <- if(identical(fileTy,"csvEur")) {
        try(utils::write.csv2(as.matrix(format(out, digits=digits)), filename, row.names=FALSE, quote=FALSE), silent=silent)
      } else try(utils::write.csv(as.matrix(format(out, digits=digits)), filename, row.names=FALSE, quote=FALSE), silent=silent)
      if(inherits(tmp, "try-error")) message(fxNa,"Note: Did not manage to write results to file '",filename,"', check for rights to write  ...") else {
        if(!silent) message(fxNa,"Wrote results successfully to file '",filename,"'")}
    }
    out
  } else {if(!silent) message(fxNa,"No results pass thresholds"); return(NULL)} }
           
