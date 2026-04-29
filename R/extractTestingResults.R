#' Extract Results From Moderated t-tests
#'
#' This function allows convenient access to results produced using the functions \code{\link[wrMisc]{moderTestXgrp}} or \code{\link[wrMisc]{moderTest2grp}}.
#' The user can define the threshold and which type of multiple testing correction should be used
#'  (as long as the multiple testing correction method cited was actually performed as part of testing).
#'  
#' @details
#' One single given comparison has to be selected by the user via argument \code{compNo}, multiple options exist :
#' 1) The use may give and index which refers to the i'th comparison stored in the object \code{stat}
#' 2) It may be more convenient to directy cite a goven comparison (eg \code{compNo=c("B-C")} in the examples below). 
#'   Such comparions may be 'inversed' to the ones orginally performed, the M-values will be adjusted accordingly (eg \code{compNo=c("C-B")})
#' 3) It is also possible to give a vector of 2 integers which will be uses as index to the levels of the initial group-assignment, beware that levels are sorted.
#' 
#' The number of columns returned can be modulated using the argument \code{addTy}. 
#' \code{addTy='slim'} gives the fewer columns, only FDR and FC of current comparison specified in argument \code{compNo} will be returned.
#' With \code{addTy='slim2'} the column-names for statistical testing and log Fold-Change will be generic, all colnames will be trimmed to max 25 characters
#' \code{addTy='all'} results in a rather exaustive collection of columns (including annotation data, if present), 
#' The option \code{addTy='slim3'} is like slim2 but without column ProteinName
#' while default \code{addTy='allMeans'} gives a compromise for the number of columns returned.
#' 
#' @param stat ('MArrayLM'-object or list) designed for the output from \code{moderTest2grp} or \code{moderTestXgrp}
#' @param compNo (integer) the comparison name/number/index to be used
#' @param statTy (character) the multiple-testing correction type to be considered when looking for significant changes  with threshold \code{thrsh} (depends on which have been run initially with \code{moderTest2grp} or \code{moderTestXgrp})
#' @param thrsh (numeric) the threshold to be applied on \code{statTy} for the result of the statistcal testing (after multiple testing correction)
#' @param FCthrs (numeric) Fold-Change threshold given as Fold-change and NOT log2(FC), default at 1.5 (for filtering at M-value =0.585)
#' @param annotCol (character) column-names from the annotation to be included
#' @param nSign (integer) number of significant digits whe returning results
#' @param addTy (character) additional groups of columns to add ("allMeans", "all" or "slim") in addition;
#'   if \code{addTy="slim"} only the FDR- and FC-columns directly concerned will be extracted/displayed, column-names will be trimmed to max 26 characters;
#'   if \code{addTy="slim2"} only the FDR- and FC-columns directly concerned will be extracted/displayed, FDR- and FC- columns will have generic names (without indicating the comparison), all text will be trimmed to max 38 characters;
#'   if \code{addTy="slim3"} like slim2 but without ProteinName
#' @param sortBy (character or logical) optional sorting of results, must specify a column-name of output, if \code{sortBy=TRUE} the fisrt column of argument \code{statTy} will be used
#' @param rowNames (character or logical) optional custom rownames, or if \code{rowNames=FALSE} no rownames
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
extractTestingResults <- function(stat, compNo=1, statTy="BH", thrsh=0.05, FCthrs=1.5, annotCol=c("Accession","EntryName","GeneName","ProteinName"), 
  nSign=6, addTy=c("allMeans"), sortBy=NULL, rowNames=TRUE, filename=NULL, fileTy="csvUS", silent=FALSE, debug=FALSE, callFrom=NULL) {
  ##
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="extractTestingResults")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  out <- outX <- NULL
  fdrOk <- FALSE
  argNa <- deparse(substitute(stat))
  
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
  if(debug) {message(fxNa,"eTR00"); eTR00 <- list(stat=stat,compNo=compNo,statTy=statTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,addTy=addTy,filename=filename,groupSep=groupSep)}      

  ## main extracting
  ## redo sample-pair assoc
  avCol <- wrMisc::sampNoDeMArrayLM(stat, compNo, outTy="list", silent=silent, callFrom=fxNa, debug=debug)          # 
  if(debug) {message(fxNa,"eTR0"); eTR0 <- list(stat=stat,compNo=compNo,statTy=statTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,addTy=addTy,avCol=avCol,filename=filename, groupSep=groupSep)}      

  ## filtering ?  normally already taken care of during testing
  
  ## need FC values
  if(!"means" %in% names(stat)) stop(fxNa,"Need list-elment stat$means - but is ABSENT !") else {
    if(!any(c("pwGrpNa","pwInd") %in% names(avCol))) stop(fxNa,"Need list-element  avCol$pwGrpNa  (or avCol$avInd) - but is ABSENT !")}
  logFC <- if(length(avCol$pwGrpNa) >1) try(stat$means[,avCol$pwGrpNa[1]] - stat$means[,avCol$pwGrpNa[2]]) else try(stat$means[,avCol$avInd[1]] - stat$means[,avCol$avInd[2]])  ## note that order of columns may not be alphabetically
  if(inherits(logFC, "try-error") || length(logFC)==0) { fcOk <- FALSE
    if(!silent) message(fxNa,"Having PROBLEMS calulating M-values (FC-values), possible problem of indexing or non-numeric values !  ... nothing to do")    
  } else {
    fcOk <- if(length(FCthrs) ==1) abs(logFC) > FCthrs else rep(TRUE, length(logFC))
    chNa <- is.na(fcOk)
    if(any(chNa)) fcOk[which(chNa)] <- FALSE
    ## check for Fdr results
    fdrOk <- stat[[useFdrTy]][,compNo] < thrsh
    chNa <- is.na(fdrOk)
    if(any(chNa)) fdrOk[which(chNa)] <- FALSE
    if(debug) {message(fxNa,"eTR0a"); eTR0a <- list(stat=stat,compNo=compNo,statTy=statTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,addTy=addTy,filename=filename, avCol=avCol,groupSep=groupSep,useFdrTy=useFdrTy,fcOk=fcOk,fdrOk=fdrOk)}        
  }
  
  if(any(fcOk & fdrOk)) {
    extrLi <- if(TRUE) which(fcOk & fdrOk) else which(fcOk | fdrOk)
    if("slim3" %in%  addTy && isTRUE("ProteinName" %in% annotCol)) annotCol <- annotCol[-which(annotCol =="ProteinName")]
    
    useCompNo <- if(any(grepl("^slim", addTy))) compNo else c(compNo,(1:ncol(stat[[useFdrTy]]))[-compNo])   # index of targert comparison and remaining comparisons
    if(debug) {message(fxNa,"eTR0b"); eTR0b <- list(stat=stat,compNo=compNo,statTy=statTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,addTy=addTy,filename=filename, avCol=avCol,groupSep=groupSep,useFdrTy=useFdrTy,fcOk=fcOk,fdrOk=fdrOk)}      
    ## some results .. continue
    if(any(c("all","allFDR","slim","slim2","slim3") %in% addTy) && ncol(stat[[useFdrTy]]) >1) {    # possibility of extracting multiple FDR
      out <- stat[[useFdrTy]][extrLi,, drop=FALSE]
      out <- out[,if("all" %in% addTy)  unique(c(compNo, 1:ncol(out))) else useCompNo, drop=FALSE]            # place comparison of interest first (if many retained)
      if(debug) {message(fxNa,"eTR0c"); eTR0c <- list(stat=stat,compNo=compNo,statTy=statTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,addTy=addTy,filename=filename, avCol=avCol,logFC=logFC,extrLi=extrLi,useCompNo=useCompNo,useFdrTy=useFdrTy,chNa=chNa,out=out)}  
    } else out <- data.frame(FDR=stat[[useFdrTy]][extrLi, compNo, drop=FALSE])
    colnames(out) <- paste0(useFdrTy,".", colnames(out))
    if(debug) {message(fxNa,"eTR1"); eTR1 <- list(stat=stat,compNo=compNo,statTy=statTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,addTy=addTy,filename=filename, avCol=avCol,logFC=logFC,extrLi=extrLi,useCompNo=useCompNo,useFdrTy=useFdrTy,chNa=chNa,out=out)}  
    
  
    if(any(c("all","allFC","slim","slim2","slim3") %in% addTy)) addTy <- c(addTy,"allMeans")         # Do we realy want to report all means when using slim ??

    ## prepare FC
    #if(any(c("all","allFC","slim","slim2") %in% addTy)) {     #  && length(avCol) >2
    #  ## index of grp-means cited in comparison
    #  #curComp <- colnames(stat$lods)[compNo]
    #  outX <- sapply(useCompNo, function(x) wrMisc::sampNoDeMArrayLM(stat, x, outT=useFdrTy))        # need to finish !!
    #  out2 <- (stat$means[,outX[[1]][2,], drop=FALSE] - stat$means[,outX[[1]][1,], drop=FALSE]) [extrLi, useCompNo, drop=FALSE]    ## creat M-values (FC)      
    #  colnames(out2) <- if(length(useCompNo) >1) paste0("logFC.",apply(outX,2, function(x) paste(colnames(stat$means)[x],collapse="-"))) else "logFC"
    #} else {out2 <- as.matrix(logFC[extrLi]); colnames(out2) <- paste0("logFC.",colnames(stat[[useFdrTy]])[compNo]) }
    out2 <- as.matrix(logFC[extrLi])
    colnames(out2) <- paste0("logFC.",colnames(stat[[useFdrTy]])[compNo])
    if(debug) {message(fxNa,"eTR1b"); eTR1b <- list(stat=stat,compNo=compNo,statTy=statTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,addTy=addTy,filename=filename, avCol=avCol,logFC=logFC,extrLi=extrLi,useFdrTy=useFdrTy,chNa=chNa,out=out,out2=out2,outX=outX,avCol=avCol)}  

    ## prepare means
    ch2 <- any(grepl("^slim", addTy))
    ch1 <- any(c("all","allMeans") %in% addTy)  && !ch2  # || length(avCol) <3)
    out3 <- if(ch1) stat$means[extrLi,] else stat$means[extrLi, which(colnames(stat$means) %in% avCol$pwGrpNa), drop=FALSE]         # the group-means
    if(length(dim(out3)) >1) colnames(out3) <- paste0("av.",colnames(out3))
    if(any(grepl("^slim[2-3]", addTy))) {
      colnames(out) <- if(any(c("BH","BY") %in% statTy)) paste0(statTy[1],".FDR") else statTy[1] 
      colnames(out2) <- "logFC" 
    }
    out <- if("annot" %in% names(stat)) signif(data.frame(out, out2, out3), nSign) else signif(cbind(out, out2, out3), nSign)                 # combine FDR, logFC & means

    ## adjust annotCol accoring to addTy
    if("all" %in% addTy && "annot" %in% names(stat) && "ProteinName" %in% colnames(stat$annot)) annotCol <- unique(c(annotCol, "EntryName","ProteinName"))
    outAnn <- stat$annot[extrLi, wrMisc::naOmit(match(annotCol,colnames(stat$annot))), drop=FALSE]
    ## trim to max char size if slim ?
    
    if(length(annotCol) >0) out <- cbind(outAnn, out)  # add annotation
    if(debug) {message(fxNa,"eTR2"); eTR2 <- list(stat=stat,out=out,out2=out2,compNo=compNo,statTy=statTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,addTy=addTy,filename=filename,outAnn=outAnn,sortBy=sortBy,rowNames=rowNames, logFC=logFC,chNa=chNa,out=out,useFdrTy=useFdrTy)}  
    
    ## optional writing to file
    if(length(filename)==1) {
      digits <- min(nSign, 12)
      tmp <- if(identical(fileTy,"csvEur")) {
        try(utils::write.csv2(as.matrix(format(out, digits=digits)), filename, row.names=FALSE, quote=FALSE), silent=silent)
      } else try(utils::write.csv(as.matrix(format(out, digits=digits)), filename, row.names=FALSE, quote=FALSE), silent=silent)
      if(inherits(tmp, "try-error")) message(fxNa,"Note: Did not manage to write results to file '",filename,"', check for rights to write  ...") else {
        if(!silent) message(fxNa,"Wrote results successfully to file '",filename,"'")}
    }
    ## final optional trimming of text
    if(any(grepl("^slim", addTy))) {
      ## remove all empty cols
      chNA <- colSums(is.na(out))==nrow(out)
      if(any(chNA)) out <- out[,which(!chNA), drop=FALSE]
      ## trim colnames
      nColChar <- nchar(colnames(out)) >26
      if(any(nColChar)) colnames(out)[which(nColChar)] <- substr(colnames(out)[which(nColChar)], 1, 26)
      ## trim  text to max length
      if(any(grepl("^slim[2-3]", addTy))) {
        if(length(outAnn) >0) {
          nChar <- nchar(outAnn) >38
          if(any(nChar)) outAnn[which(nChar)] <- paste0(substr(outAnn[which(nChar)], 1, 38),"..") }
      }
    }  
    ## Custom Sorting
    if(length(sortBy)==1) {               # prepare custom sorting
      if(isTRUE(sortBy)) if(length(statTy)==0) {
          sortBy <- NULL
          if(!silent) message(fxNa,"Unable to sort with argument 'statTy' being empty")
        } else sortBy <- statTy[1]
    }  
    if(length(sortBy)==1 && is.character(sortBy)) {               # custom sorting    
      sortBy1 <- which(colnames(out) %in% sortBy)
      if(length(sortBy1) ==0) sortBy1 <- grep(paste0("^",sortBy,"\\."), colnames(out))  # may be composed colname
      if(length(sortBy1) >0) out <- out[order(out[,sortBy1[1]]),,drop=FALSE] else {
        if(!silent) message(fxNa,"Unable to understand argument 'sortBy' (not found in colnames of output)")
      } 
    }
    ## customize rownames
    if(length(rowNames) !=0 && !any(is.na(rowNames))) {                 
      if(isFALSE(rowNames)) rownames(out) <- NULL else {
        if(length(rowNames) ==nrow(out) && sum(duplicated(rowNames))==0) rownames(out) <- rowNames
      }
    }
    out
  } else {if(!silent) message(fxNa,"No results pass thresholds"); return(NULL)} }
              

