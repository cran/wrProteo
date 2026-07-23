#' Extract Results From Moderated t-tests
#'
#' This function allows convenient access to results produced using the functions \code{\link[wrMisc]{moderTestXgrp}} or \code{\link[wrMisc]{moderTest2grp}}.
#' The user can define the threshold and which type of multiple testing correction should be used
#'  (as long as the multiple testing correction method cited was actually performed as part of testing).
#'  
#' @details
#' One single given comparison has to be selected by the user via argument \code{compNo}, multiple options exist :
#' 1) The use may give and index which refers to the i'th comparison stored in the object \code{stat$setup}
#' 2) It may be more convenient to directy cite a given comparison (eg \code{compNo=c("B-C")} in the examples below). 
#'   Such comparions may be 'inversed' to the ones orginally performed, the M-values will be adjusted accordingly (eg \code{compNo=c("C-B")})
#' 3) It is also possible to give a vector of 2 integers which will be uses as index to the levels of the initial group-assignment, beware that levels are sorted.
#' 
#' The number and type of columns returned can be modulated using the argument \code{outTy}. 
#' \code{outTy='slim'} gives the fewer columns, only FDR and FC of current comparison specified in argument \code{compNo} will be returned.
#' With \code{outTy='slim2'} the column-names for statistical testing and log Fold-Change will be generic, all colnames will be trimmed to max 25 characters
#' \code{outTy='all'} results in a rather exaustive collection of columns (including annotation data, if present), 
#' The option \code{outTy='slim3'} is like slim2 but without column ProteinName
#' while default \code{outTy='allMeans'} gives a compromise for the number of columns returned.
#' 
#' @param stat ('MArrayLM'-object or list) designed for the output from \code{moderTest2grp} or \code{moderTestXgrp}
#' @param compNo (integer) the comparison name/number/index to be used
#' @param statTy (character) the multiple-testing correction type to be considered when looking for significant changes  with threshold \code{thrsh} (depends on which have been run initially with \code{moderTest2grp} or \code{moderTestXgrp})
#' @param thrsh (numeric) the threshold to be applied on \code{statTy} for the result of the statistcal testing (after multiple testing correction)
#' @param FCthrs (numeric) Fold-Change threshold given as Fold-change and NOT log2(FC), default at 1.5 (for filtering at M-value =0.585)
#' @param annotCol (character) column-names from the annotation to be included
#' @param nSign (integer) number of significant digits whe returning results
#' @param outTy (character) additional groups of columns to add ("allMeans", "all" or "slim") in addition;
#'   if \code{outTy="slim"} only the FDR- and FC-columns directly concerned will be extracted/displayed, column-names will be trimmed to max 26 characters;
#'   if \code{outTy="slim2"} only the FDR- and FC-columns directly concerned will be extracted/displayed, FDR- and FC- columns will have generic names (without indicating the comparison), all text will be trimmed to max 38 characters;
#'   if \code{outTy="slim3"} like slim2 but without ProteinName
#' @param sep (character of length=1 or \code{NULL}) custom choice for separator for column-names of pairwise comparisons, otherwise based on \code{stat} or (if \code{outTy="slim3"} or \code{outTy="slim4"}) '_vs_' (or '__vs__' if '-' or '_' appear in groups)
#' @param addTy  this argument has been deprecated and replaced by \code{outTy} (content of addTy will be ignored) 
#' @param sortBy (character or logical) optional sorting of results according to argument \code{compNo}, may be name of mutiple-testing correction (eg 'BH', see also argument \code{statTy}), if \code{sortBy=TRUE} the first column of argument \code{statTy} will be used
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
#' extractTestingResults(test8, outTy="all")
#' extractTestingResults(test8, outTy="slim")
#' @export
extractTestingResults <- function(stat, compNo=1, statTy="BH", thrsh=0.05, FCthrs=1.5, annotCol=c("Accession","EntryName","GeneName","ProteinName"), 
  nSign=6, outTy="slim", sep=NULL, addTy=NULL, sortBy=NULL, rowNames=TRUE, filename=NULL, fileTy="csvUS", silent=FALSE, debug=FALSE, callFrom=NULL) {
  ##
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="extractTestingResults")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) { silent <- FALSE } else { debug <- FALSE }
  shortColNa <- FALSE     # in case of extracting single comparison shorten colnames for FDR & logFC (ie not citing groups)
  #later#newSep <- if(any("slim2" %in% outTy)) "_vs_" else "__vs__"                          # note cbind will change colnames with '  /  ' to '.....'
  
  out <- outX <- NULL
  fdrOk <- FALSE
  argNa <- deparse(substitute(stat))
  
  if(!"list" %in% mode(stat) || length(stat) <1) stop("Argument 'stat' must be a list or 'MArrayLM'-object from limma")
  if(!("MArrayLM" %in% class(stat)) && !silent) message(fxNa," Caution, '",argNa,"' is not 'MArrayLM'-object as expected")
  if(length(statTy) ==0) { statTy <- "BH"
    if(!silent) message(fxNa,"Argument 'statTy' empty, setting to default 'BH'")}
  useFdrTy <- if(identical(statTy,"BH") && "FDR" %in% names(stat)) "FDR" else statTy    # maybe not needed to force to FDR
  datOK <- all(c("means", useFdrTy) %in% names(stat) )
  
  hasAnn <- "annot" %in% names(stat) && length(stat$annot) != 0 
  if(!hasAnn && !silent) message(fxNa,"Note : no list-element 'annot' found in '",argNa,"' (would be useful for better annotation of results)")
  chLstEl <- c(useFdrTy) %in% names(stat)   
  if(!chLstEl[1]) stop(fxNa,"Cannot find list-element '",useFdrTy,"' in 'stat'")

  if(length(compNo) != 1) stop("Argument 'compNo' must be numeric and of length=1")
  if(compNo > ncol(stat[[useFdrTy[1]]])) { compNo <- 1
    message(fxNa,"Invalid entry of argument 'compNo', setting to defaul compNo=1")}   
  if(is.na(FCthrs) || !is.numeric(FCthrs)) FCthrs <- NULL
  if(length(FCthrs) >0 && is.numeric(FCthrs)) FCthrs <- log2(FCthrs) else FCthrs <- NULL
  if(length(outTy) != 1 || any(is.na(outTy))) outTy <- "slim"
  
  ## main extracting
  ## redo sample-pair assoc
  avCol <- wrMisc::sampNoDeMArrayLM(stat, compNo, outTy="list", silent=silent, callFrom=fxNa, debug=debug)          # 
  if(debug) {message(fxNa,"eTR0"); eTR0 <- list(stat=stat,compNo=compNo,statTy=statTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,outTy=outTy,avCol=avCol,filename=filename,sep=sep, hasAnn=hasAnn)}      
  
  ## define separator for different types of extracting
  newSep <- if(length(sep) !=1 || identical(sep,NA)) { 
    if(grepl("^slim[3-4]", outTy)) { if(any(grepl("_|-", stat$setup$grpNa))) "__vs__" else "_vs_"} else stat$setup$sep    # special sep for slim3 & slim4
  } else sep
  #newSep <- if("slim" %in% outTy) else stat$setup$sep if(grepl("^slim[3-4]", outTy)) "_vs_" else "__vs__"                          # note:  cbind will change colnames with '  /  ' to '.....' 
  if(any(grepl(newSep, stat$setup$grpNa), na.rm=TRUE)) warning(fxNa,"NOTE : (New) separator 'sep' for pairwise groups seems not well chosen - occuring inside some of the group-names !")
  shortColNa2 <- grepl("^slim4", outTy)    # give only 'FDR' or 'logFC' if single column/comparison
      
  useCompNo <- if(any(grepl("^slim", outTy)) || ncol(stat$setup$pwGrpNa) ==1) compNo else c(compNo,(1:ncol(stat[[useFdrTy[1]]]))[-compNo])   # index of targert comparison and remaining comparisons
  if(debug) {message(fxNa,"eTR0a"); eTR0a <- list(stat=stat,compNo=compNo,useFdrTy=useFdrTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,outTy=outTy,avCol=avCol,filename=filename, hasAnn=hasAnn,newSep=newSep,useCompNo=useCompNo)}      
    
  if(datOK) {      
    ## prepare names & co
    if(is.null(names((stat$setup$grpNa)))) { 
      shortGrp <- if(length(stat$setup$grpNa) <27) LETTERS else {
        if(length(stat$setup$grpNa) <27) paste0(rep(LETTERS, each=26), rep(LETTERS, 26)) else paste0(rep(LETTERS, each=676), rep(LETTERS, each=26), rep(LETTERS, 676))}
      shortGrp <- shortGrp[1:length(stat$setup$grpNa)]
      names(shortGrp) <- stat$setup$grpNa
      names(stat$setup$grpNa) <- shortGrp
    } else { shortGrp <- names(stat$setup$grpNa) }

    useGrpNaFull <- stat$setup$pwGrpNa
    if(any(c("slim3","slim4") %in% outTy)) {                                         
      useGrpNa <- matrix(shortGrp[stat$setup$pwIndex], nrow=nrow(stat$setup$pwIndex))
    } else { 
      useGrpNa <- stat$setup$pwGrpNa }
        
    usePwNa <- paste0(useGrpNa[,1], newSep, useGrpNa[,2]) 
    usePwNaFull <- paste0(stat$setup$pwGrpNa[,1], newSep, stat$setup$pwGrpNa[,2]) 
    if(any(grepl("^slim", outTy)) && nrow(stat$setup$pwIndex) >1) {     ## adjust order to bring question of interest to front
      useGrpNa <- useGrpNa[useCompNo,, drop=FALSE]
      useGrpNaFull <- useGrpNaFull[useCompNo,, drop=FALSE]
      usePwNaFull <- usePwNaFull[useCompNo]
      usePwNa <- usePwNa[useCompNo] }
    if(debug) {message(fxNa,"eTR0b"); eTR0b <- list(stat=stat,compNo=compNo,useFdrTy=useFdrTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,outTy=outTy,avCol=avCol,useFdrTy=useFdrTy,filename=filename, shortColNa2=shortColNa2,hasAnn=hasAnn,newSep=newSep,useCompNo=useCompNo,useGrpNa=useGrpNa,usePwNa=usePwNa,usePwNaFull=usePwNaFull)}
               
    ## FDR results
    out <- as.data.frame(stat[[useFdrTy[1]]][ , useCompNo, drop=FALSE] )   
    names(out) <- if(shortColNa2 && ncol(out) ==1) "FDR" else { paste0(useFdrTy,".", if(!"lfdr" %in% useFdrTy) "FDR.", usePwNa)}
    if(length(useFdrTy) >1) {
      for(i in 2:length(useFdrTy)) out <- cbind(out, stat[[useFdrTy[i]]][ , useCompNo, drop=FALSE])
    }
      
    ## FC results                            if(shortColNa2 && ncol(out)==1)
    logFC <- as.matrix(apply(stat$setup$pwGrpNa[useCompNo,, drop=FALSE], 1, function(grN) stat$means[,grN[1]] - stat$means[,grN[2]]))
    if(length(useCompNo) !=1) logFC <- logFC[ , useCompNo, drop=FALSE]       # place main question/comparison first
    colnames(logFC) <- if(shortColNa2 && ncol(logFC) ==1) "logFC" else paste0("logFC.",usePwNa)

    ## prepare filter 
    fcOK <- if(length(FCthrs) ==1) abs(logFC[,1]) > FCthrs else rep(TRUE, nrow(logFC))
    fdrOK <- out[,1] < thrsh
    extrLi <- if("all" %in% outTy) rep(TRUE, nrow(out)) else { fcOK & fdrOK }
    chNA <- is.na(extrLi)
    if(any(chNA)) extrLi[which(chNA)] <- FALSE
    if(debug) { message(fxNa,"eTR0c"); eTR0c <- list(stat=stat,compNo=compNo,useFdrTy=useFdrTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,outTy=outTy,avCol=avCol,useFdrTy=useFdrTy,filename=filename, shortColNa2=shortColNa2,hasAnn=hasAnn,newSep=newSep,useCompNo=useCompNo,useGrpNa=useGrpNa,usePwNa=usePwNa,out=out,logFC=logFC,fcOK=fcOK,fdrOK=fdrOK,extrLi=extrLi,shortGrp=shortGrp)}
  }  
    
  if(datOK && sum(extrLi) !=0) {
    extrLi <- which(extrLi)
    
    ## group means
    useMeanCol <- if(grepl("^slim[[:digit:]]",outTy)) match(useGrpNaFull[1,], colnames(stat$means)) else 1:ncol(stat$means)
    outD <- stat$means[, useMeanCol ]
    colnames(outD) <- paste0("av.", if(grepl("^slim[3-4]",outTy)) shortGrp[match(colnames(outD),names(shortGrp))] else colnames(outD))
    
    out <- signif(cbind(out, logFC, outD)[extrLi,, drop=FALSE], nSign)                 # combine FDR, logFC & means
    if(debug) { message(fxNa,"eTR2"); eTR2 <- list(stat=stat,compNo=compNo,statTy=statTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,outTy=outTy,avCol=avCol,useFdrTy=useFdrTy,filename=filename, shortColNa2=shortColNa2,hasAnn=hasAnn,newSep=newSep,useCompNo=useCompNo,useGrpNa=useGrpNa,usePwNa=usePwNa,out=out,outD=outD,logFC=logFC,fcOK=fcOK,fdrOK=fdrOK,extrLi=extrLi,shortGrp=shortGrp)}
    
    ## number of NAs per group of replicates & protein
    if(length(stat$isNA) ==0 && length(stat$nNA) ==0) {                # && packageVersion("wrMisc") < 2.1.0     # primarily needed with older wrMisc 
      chDat <- c("raw","quant") %in% names(stat)
      if(any(chDat)) stat$isNA <- is.na(stat[[which(chDat)[1]]])  
    }
    if(length(stat$nNArat) ==0 && length(stat$isNA) !=0) {
      stat$nNArat <- wrMisc::rowGrpNA(stat, mode="complete")  # gives '0/3' etc
      if(all(grepl(".nNA$", colnames(stat$nNArat)))) colnames(stat$nNArat) <- paste0("nNA.",sub("X{0,1}\\.nNA$","", colnames(stat$nNArat) ))     # shuffle 'nNA' to beginning of colnames to remain similar
    } 
    if(length(stat$nNArat) !=0 && grepl("slim[3-4]", outTy)) {                                      # adjust names for nNA to shortGrp (eg nNA.A)  etc
      colnames(stat$nNArat) <- paste0("nNA.",shortGrp[match(names(shortGrp),colnames(stat$means))])      
    }           
    if(debug) { message(fxNa,"eTR3"); eTR3 <- list(stat=stat,compNo=compNo,useFdrTy=useFdrTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,outTy=outTy,avCol=avCol,filename=filename, hasAnn=hasAnn,newSep=newSep,useCompNo=useCompNo,useGrpNa=useGrpNa,usePwNa=usePwNa,out=out,outD=outD,logFC=logFC,fcOK=fcOK,fdrOK=fdrOK,extrLi=extrLi,shortGrp=shortGrp,useMeanCol=useMeanCol,sortBy=sortBy)}
    
    ## annotation
    if(hasAnn) {
      if("slim3" %in% outTy && isTRUE("ProteinName" %in% annotCol)) {annotCol <- annotCol[-which(annotCol =="ProteinName")] }     
      if(grepl("^slim", outTy)) {                   # remove all NA columns from annot (unless all cols would ghet removed)
        isNA <- colSums(!is.na(stat$annot[extrLi, annotCol, drop=FALSE])) ==0
        if(any(isNA) && !all(isNA)) annotCol <- annotCol[which(!isNA)]                 # adjust if some cols maissing
      }
      useAnnLi <- match( rownames(out), stat$annot[,"EntryName"])
      if(all(is.na(useAnnLi))) useAnnLi <- match( rownames(out), stat$annot[,"Accession"])
      if(debug && any(is.na(useAnnLi))) message(fxNa,"Note : ","Adding Annotation : Mapping of gave ",sum(is.na(useAnnLi)),"/",length(useAnnLi)," NAs !")

      out <- data.frame(stat$annot[useAnnLi, annotCol, drop=FALSE], stat$nNArat[extrLi, useMeanCol, drop=FALSE], out)   ## combine &  apply filter
    }     
     
    ## sort
    newOrd <- 1:nrow(out)    # init default
    if(length(sortBy) ==1) { if(isTRUE(sortBy)) { 
        newOrd <- order(stat[[useFdrTy[1]]][ extrLi, useCompNo], decreasing=FALSE ) 
      } else { if(sortBy %in% names(stat) && ncol(stat[[sortBy]]) >= useCompNo) { newOrd <- order(stat[[sortBy]][ extrLi, useCompNo], decreasing=FALSE ) 
        } else if(!isFALSE(sortBy)) warning(fxNa,"Invalid entry for argument 'sortBy' (The corresponding element doesn't have sufficient columns compared to argument 'compNo')")
      }
    } 

    out <- out[newOrd,, drop=FALSE]
    if(length(rowNames)==nrow(out) && length(rowNames) !=0 && !any(duplicated(rowNames))) rownames(out) <- rowNames
    if(isFALSE(rowNames)) rownames(out) <- NULL
    if(debug) { message(fxNa,"eTR4"); eTR4 <- list(stat=stat,compNo=compNo,useFdrTy=useFdrTy,thrsh=thrsh,FCthrs=FCthrs,annotCol=annotCol,outTy=outTy,avCol=avCol,filename=filename, hasAnn=hasAnn,newSep=newSep,useCompNo=useCompNo,useGrpNa=useGrpNa,usePwNa=usePwNa,out=out,outD=outD,logFC=logFC,fcOK=fcOK,fdrOK=fdrOK,extrLi=extrLi,shortGrp=shortGrp,useMeanCol=useMeanCol,newOrd=newOrd)}
      
    ## export
    if(length(filename)==1 && length(out) !=0 && length(dim(out))==2) {
      if(any(c("csvUS","UScsv") %in% fileTy)) {
        ch1 <- try(utils::write.csv(out, file=filename, row.names=FALSE, quote=FALSE))
      } else {
        ch1 <- try(utils::write.csv2(out, file=filename, row.names=FALSE, quote=FALSE))
      }
      if(inherits(ch1, "try-error")) { warning("FAILED to export as file '",filename,"'") } }
  }
      
   ## return results
   out
}
      
   
