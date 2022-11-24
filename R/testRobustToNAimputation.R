#' Pair-wise testing robust to NA-imputation
#'
#' \code{testRobustToNAimputation} replaces \code{NA} values based on group neighbours (based on grouping of columns in argument \code{gr}), following overall assumption of close to Gaussian distribution.
#' Furthermore, it is assumed that \code{NA}-values originate from experimental settings where measurements at or below detection limit are recoreded as \code{NA}.
#' In  such cases (eg in proteomics) it is current practice to replace \code{NA}-values by very low (random) values in order to be able to perform t-tests.
#' However, random normal values used for replacing may in rare cases deviate from the average (the 'assumed' value) and in particular, if multiple \code{NA} replacements are above the average, 
#' may look like induced biological data and be misinterpreted as so.      
#' The statistical testing uses \code{eBayes} from Bioconductor package \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{limma} for robust testing in the context of small numbers of replicates. 
#' By repeating multiple times the process of replacing \code{NA}-values and subsequent testing the results can be sumarized afterwards by median over all repeated runs to remmove the stochastic effect of individual NA-imputation.
#' Thus, one may gain stability towards random-character of \code{NA} imputations by repeating imputation & test 'nLoop' times and summarize p-values by median (results stabilized at 50-100 rounds).
#' It is necessary to define all groups of replicates in \code{gr} to obtain all possible pair-wise testing (multiple columns in $BH, $lfdr etc). 
#' The modified testing-procedure of Bioconductor package \href{https://bioconductor.org/packages/release/bioc/html/ROTS.html}{ROTS} may optionaly be included, if desired.
#' This function returns a \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{limma}-like S3 list-object further enriched by additional fields/elements.

#' @details
#' The argument \code{multCorMeth} allows to choose which multiple correction algorimths will be used and included to the final results.
#' Possible options are 'lfdr','BH','BY','tValTab', ROTSn='100' (name to element necessary) or 'noLimma' (to add initial p.values and BH to limma-results). By default 'lfdr' (local false discovery rate from package 'fdrtools') and 'BH' (Benjamini-Hochberg FDR) are chosen.
#' The option 'BY' referrs to Benjamini-Yakuteli FDR, 'tValTab' allows exporting all individual t-values from the repeated NA-substitution and subsequent testing.
#' 
#' This function is compatible with automatic extraction of experimental setup based on sdrf or other quantitation-specific sample annotation.
#' In this case, the results of automated importing and mining of sample annotation should be stored as \code{$sampleSetup$groups} or \code{$sampleSetup$lev}  
#' 
#' 
#' @param dat (matrix or data.frame) main data (may contain \code{NA}); if \code{dat} is list containing $quant and $annot as matrix, the element $quant will be used
#' @param gr (character or factor) replicate association; if \code{dat} contains a list-element \code{$sampleSetup$groups} or \code{$sampleSetup$lev} this may be used in case \code{gr=NULL}
#' @param annot (matrix or data.frame) annotation (lines must match lines of data !), if \code{annot} is \code{NULL} and argument \code{dat} is a list containing both $quant and $annot, the element $annot will be used 
#' @param retnNA (logical) retain and report number of \code{NA}
#' @param avSdH (numeric) population characteristics (mean and sd) for >1 \code{NA} neighbours 'high' (per line)
#' @param avSdL  depreciated argument, no longer used 
#' @param plotHist (logical) additional histogram of original, imputed and resultant distribution (made using \code{\link{matrixNAneighbourImpute}} )
#' @param xLab (character) custom x-axis label 
#' @param tit (character) custom title
#' @param imputMethod (character) choose the imputation method (may be 'mode2'(default), 'mode1', 'datQuant', 'modeAdopt', 'informed' or 'none', for details see \code{\link{matrixNAneighbourImpute}} )
#' @param seedNo (integer) seed-value for normal random values
#' @param multCorMeth (character) define which method(s) for correction of multipl testing should be run (for choice : 'BH','lfdr','BY','tValTab', choosing several is possible)
#' @param nLoop (integer) number of runs of independent \code{NA}-imputation
#' @param lfdrInclude (logical) depreciated, please used \code{multCorMeth} instead (include lfdr estimations, may cause warning message(s) concerning convergence if few too lines/proteins in dataset tested).
#' @param ROTSn (integer) depreciated, please used \code{multCorMeth} instead (number of repeats by \code{ROTS}, if \code{NULL} \code{ROTS} will not be called)
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) This function allows easier tracking of messages produced
#' @return This function returns a limma-type S3 object of class 'MArrayLM' (which can be accessed lika a list); multiple results of testing or multiple testing correction types may get included ('p.value','FDR','BY','lfdr' or 'ROTS.BH')
#' @seealso \code{\link[wrMisc]{moderTest2grp}}, \code{\link[wrMisc]{pVal2lfdr}}, \code{eBayes} in Bioconductor package \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{limma}, \code{\link[stats]{t.test}},\code{ROTS} of Bioconductor package \href{https://bioconductor.org/packages/release/bioc/html/ROTS.html}{ROTS}   
#' @examples
#' set.seed(2015); rand1 <- round(runif(600) +rnorm(600,1,2),3)
#' dat1 <- matrix(rand1,ncol=6) + matrix(rep((1:100)/20,6),ncol=6)
#' dat1[13:16,1:3] <- dat1[13:16,1:3] +2      # augment lines 13:16 
#' dat1[19:20,1:3] <- dat1[19:20,1:3] +3      # augment lines 19:20
#' dat1[15:18,4:6] <- dat1[15:18,4:6] +1.4    # augment lines 15:18 
#' dat1[dat1 <1] <- NA                        # mimick some NAs for low abundance
#' ## normalize data
#' boxplot(dat1, main="data before normalization")
#' dat1 <- wrMisc::normalizeThis(as.matrix(dat1), meth="median")
#' ## designate replicate relationships in samples ...  
#' grp1 <- gl(2, 3, labels=LETTERS[1:2])                   
#' ## moderated t-test with repeated inputations (may take >10 sec,  >60 sec if ROTSn >0 !) 
#' PLtestR1 <- testRobustToNAimputation(dat=dat1, gr=grp1, retnNA=TRUE, nLoop=70)
#' names(PLtestR1)
#' @export
testRobustToNAimputation <- function(dat, gr=NULL, annot=NULL, retnNA=TRUE, avSdH=c(0.15,0.5), avSdL=NULL, plotHist=FALSE, xLab=NULL, tit=NULL, imputMethod="mode2", 
  seedNo=NULL,  multCorMeth=NULL, nLoop=100, lfdrInclude=NULL, ROTSn=NULL, silent=FALSE, debug=FALSE, callFrom=NULL) {
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="testRobustToNAimputation")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  if(!isTRUE(plotHist)) plotHist <- FALSE
  datOK <- TRUE
  msg <- grIni <- NULL
  ## start testing input
  if(is.list(dat)) { if(all(c("quant","annot") %in% names(dat))) {
    if(length(dim(dat$annot)) ==2 & length(annot) <1) annot <- dat$annot else if(!silent) message(fxNa,"Invalid '$annot'") # recover$annot if not given separately
    if("sampleSetup" %in% names(dat) & length(gr) <1) {
      gr <- if("groups" %in% names(dat$sampleSetup)) dat$sampleSetup$groups else dat$sampleSetup$lev
      if(all(match(gr, unique(gr)) ==match(names(gr), unique(names(gr))))) gr <- names(gr)              # rather use name instead of index
      if(all(length(gr) >1, length(dat$sampleSetup$col) <2, length(names(gr)) <1)) names(gr) <- dat$sampleSetup$sdrf[,dat$sampleSetup$col]       # in case not names provided
    } else grIni <- gr
    dat <- dat$quant } else { datOK <- FALSE; msg <- "Invalid 'dat' : does NOT contain both '$quant' and '$annot !"} } 
  if(datOK) { if(length(unique(gr))==length(gr)) { datOK <- FALSE
    msg <- "Argument 'gr' is empty or does not design any replicates !  (nothing to do)"} } 
  if(datOK) if(any(length(dim(dat)) !=2, dim(dat) < 1:2, na.rm=TRUE)) { datOK <- FALSE
    msg <- "'dat' must be matrix or data.frame with >1 columns"}  
  if(datOK) {
    if(is.data.frame(dat)) dat <- as.matrix(dat)
    if(length(gr) != ncol(dat)) { datOK <- FALSE
      msg <- "Number of columns in 'dat' and number of (group-)elements in 'gr' do NOT match !"} }
  if(datOK) {
    if(!is.factor(gr)) gr <- as.factor(gr)
    if(is.null(xLab)) xLab <- "values"            
    if(length(annot) <1) annot <- matrix(NA, nrow=nrow(dat), ncol=1, dimnames=list(rownames(dat),"rowNa"))
    if(length(ROTSn) >0) message(fxNa,"Argument 'ROTSn' is depreciated, please used argument 'multCorMeth' instead (like multCorMeth=c(ROTSn='10'))")
    if(length(lfdrInclude) >0) message(fxNa,"Argument 'lfdrInclude' is depreciated, please used argument 'multCorMeth' instead (like multCorMeth='lfdrInclude')")
    ## get compatible to old arguments lfdrInclude & ROTSn
    ROTSn <- NULL
    multCorMeth <- if(length(multCorMeth) <1) c("lfdr","FDR","means") else unique(c(multCorMeth, "means"))
    if(length(multCorMeth) ==1 & is.numeric(multCorMeth)) { 
       multCorMeth <- if(multCorMeth >1) c("lfdr", ROTSn=as.integer(multCorMeth), "means") else "lfdr"}
    
    if("ROTSn" %in% names(multCorMeth)) { ROTSn <- try(as.integer(multCorMeth["ROTSn"]), silent=TRUE)
      if(inherits(ROTSn, "try-error")) {ROTSn <- NULL; comp<- NULL }} else {ROTSn <- NULL; comp<- NULL }
  
    if("lfdr" %in% multCorMeth) { lfdrInclude <- TRUE
    } else if("lfdr" %in% names(multCorMeth)) { lfdrInclude <- try(as.logical(multCorMeth["lfdr"]), silent=TRUE)
        if(inherits(lfdrInclude, "try-error")) { lfdrInclude <- FALSE; multCorMeth <- multCorMeth[-which(names(multCorMeth)== "lfdr")] } }
    if(length(lfdrInclude) <1) lfdrInclude <- FALSE            # if for some reason whatsoever ...
        
    ## main
    isNA <- is.na(dat)
    chNA <- any(isNA)
    nNAmat <- matrix(0, nrow=nrow(dat), ncol=length(levels(gr)), dimnames=list(NULL,levels(gr)))
    seedNo <- as.integer(seedNo)[1]
    gr <- try(as.factor(gr)) 
    if(inherits(gr, "try-error")) message("+++++\n",fxNa," MAJOR PROBLEM with argument 'gr' !!  (possibly not sufficient level-names ?) \n+++++")
    
    
    ## 1st pass
    if(debug) message(fxNa," imputMethod=",imputMethod,"  tRN1")
    if("none" %in% tolower(imputMethod)) {
      if(debug) {message(fxNa,"Starting 1st pass,  no of NA: ",sum(isNA)); tt1 <- list(dat=dat,gr=gr,imputMethod=imputMethod,annot=annot,seedNo=seedNo, retnNA=retnNA, avSdH=avSdH,ROTSn=ROTSn,lfdrInclude=lfdrInclude,multCorMeth=multCorMeth)}
      nLoop <- 1
      chFin <- is.finite(dat)
      if(any(!chFin, na.rm=TRUE)) dat[which(!chFin)] <- NA   # need to replace Inf & -Inf by NA to avoid problems at limma::lmFit() 
      datI <- list(data=dat)
      datFi <- combineMultFilterNAimput(dat=dat, imputed=datI, grp=gr, annDat=annot, abundThr=stats::quantile(if(is.list(dat)) dat$quant else dat, 0.02,na.rm=TRUE), silent=silent, callFrom=fxNa)  # number of unique peptides unknown !
       
    } else {
      datI <- matrixNAneighbourImpute(dat, gr, imputMethod=imputMethod, retnNA=retnNA ,avSdH=avSdH, plotHist=plotHist, xLab=xLab, tit=tit, seedNo=seedNo, silent=silent, callFrom=fxNa,debug=debug)
      if(debug) {message(fxNa,"Start combineMultFilterNAimput   tt2a")}
      datFi <- combineMultFilterNAimput(dat=dat, imputed=datI, grp=gr, annDat=annot, abundThr=stats::quantile(if(is.list(dat)) dat$quant else dat, 0.02,na.rm=TRUE), silent=silent, callFrom=fxNa)  # number of unique peptides unknown !
      if(debug) {message(fxNa,"Done combineMultFilterNAimput   tt2b")}     # done combineMultFilterNAimput
    }
    ## prepare for testing
    if(lfdrInclude) {
      chLfdr <- try(find.package("fdrtool"), silent=TRUE)
      if(inherits(chLfdr, "try-error")) { 
        message(fxNa,"Package 'fdrtool' NOT found ! Please install first from CRAN for calculating lfdr-values. Omitting (defaut) 'lfdr' option from argument 'multCorMeth' ..")
        lfdrInclude <- FALSE } }
    if(debug) {message(fxNa,"  tRN2"); tRN2 <- list(dat=dat,gr=gr,datFi=datFi,multCorMeth=multCorMeth,datI=datI,datFi=datFi)}
    
    pwComb <- wrMisc::triCoord(length(levels(gr)))
    if(debug) message(fxNa,"Start 1st moderTestXgrp()")
    out <- wrMisc::moderTestXgrp(datFi$data, grp=gr, limmaOutput=TRUE, addResults=multCorMeth, silent=silent, callFrom=fxNa)   # can't do question specific filtering w/o explicit loop
   if(debug) {message(fxNa,"  tRN3"); tRN3 <- list(dat=dat,gr=gr,datFi=datFi,multCorMeth=multCorMeth,datI=datI,datFi=datFi,out=out,pwComb=pwComb)}
    
    chFDR <- names(out) =="FDR"
    if(any(chFDR, na.rm=TRUE)) names(out)[which(chFDR)] <- "BH"                 # rename $FDR to $BH
    rownames(pwComb) <- colnames(out$t) 
    ## need to add $ROTS.p
    if(length(ROTSn)==1) if(ROTSn >0 & !is.na(ROTSn)) {  
      chPa <- requireNamespace("ROTS", quietly=TRUE)
      if(!chPa) { message(fxNa,"Package 'RORS' not found/installed (please install from Bioconductor), omitting argument 'ROTSn'")
        ROTSn <- 0 }
    } else ROTSn <- NULL
    if(length(ROTSn)==1) if(ROTSn >0) {
      ## this requires package ROTS
      if(debug) message(fxNa,"Start ROTS   n=",ROTSn)
      comp <- wrMisc::triCoord(length(levels(gr)))
      rownames(comp) <- paste(levels(gr)[comp[,1]], levels(gr)[comp[,2]],sep="-")
      tmRO <- matrix(nrow=nrow(datFi$data), ncol=nrow(comp))
      comPair <- matrix(unlist(strsplit(rownames(comp),"-")), ncol=nrow(comp))
      useCol <- apply(comPair, 2, function(x) gr %in% x)
      for(i in 1:nrow(comp)) tmRO[which(datFi$filt[,i]),i] <- ROTS::ROTS(datFi$data[which(datFi$filt[,i]), 
        which(useCol[,i])], groups=as.numeric(as.factor(gr[which(useCol[,i])]))-1, B=ROTSn)$pvalue        # K=500  
      out$ROTS.p <- tmRO
      out$ROTS.BH <- apply(tmRO, 2, stats::p.adjust, method="BH") 
      if(lfdrInclude) out$ROTS.lfdr <- apply(tmRO, 2, wrMisc::pVal2lfdr) 
    } 
    if(debug) message(fxNa,"tRN4")
          
    ## subsequent rounds of NA-imputation  
    if(chNA & nLoop >1) { 
      if(debug) message(fxNa,"Subsequent rounds of NA-imputation   nLoop=",nLoop)
      pValTab <- tValTab <- array(NA, dim=c(nrow(dat), nrow(pwComb), nLoop))
      datIm <- array(NA, dim=c(nrow(dat), ncol(dat), nLoop))
      datIm[,,1] <- datFi$data
      pValTab[,,1] <- out$p.value
      tValTab[,,1] <- out$t
      if(length(ROTSn)==1) if(ROTSn >0) {
        pVaRotsTab <- array(NA, dim=c(nrow(dat), nrow(pwComb), min(10,nLoop)))
        pVaRotsTab[,,1] <- out$ROTS.p }
      for(i in 2:nLoop) {
        ## the repeated NA-imputation & testing   
        if(length(seedNo)==1) seedNo <- seedNo +i 
        if(debug) {message("tRN5"); tRN5 <- list(dat=dat,gr=gr,seedNo=seedNo, retnNA=retnNA, avSdH=avSdH,datI=datI, pValTab=pValTab,datIm=datIm,ROTSn=ROTSn)}
      
        datX <- matrixNAneighbourImpute(dat, gr, imputMethod=imputMethod, seedNo=seedNo, retnNA=retnNA, avSdH=avSdH, NAneigLst=datI$NAneigLst, plotHist=FALSE, silent=TRUE, callFrom=fxNa)$data
       if(debug) message(fxNa,"Passed matrixNAneighbourImpute()   in loop no ",i,"  tRN5b")
  
    #1st round# datI <- matrixNAneighbourImpute(dat, gr, seedNo=seedNo, retnNA=retnNA ,avSdH=avSdH, plotHist=plotHist, xLab=xLab, tit=tit, silent=silent, callFrom=fxNa)
    #1st round# datFi <- combineMultFilterNAimput(dat=dat, imputed=datI, grp=gr, annDat=annot, abundThr=stats::quantile(dat,0.02,na.rm=TRUE), silent=silent, callFrom=fxNa) 
    #1st round# out <- wrMisc::moderTestXgrp(datFi$data, grp=gr, limmaOutput=TRUE, addResults="", silent=silent, callFrom=fxNa)       
        fitX <- limma::eBayes(limma::contrasts.fit(limma::lmFit(datX[,], out$design), contrasts=out$contrasts))
        datIm[,,i] <- datX
        pValTab[,,i] <- fitX$p.value
        tValTab[,,i] <- fitX$t
        if(length(ROTSn)==1) if(ROTSn >0 & i < min(10, nLoop)) {       # test using ROTS (TAKES MUCH TIME !!)
          for(i in 1:nrow(comp)) tmRO[which(datFi$filt[,i]),i] <- ROTS::ROTS(datFi$data[which(datFi$filt[,i]),which(useCol[,i])], groups=as.numeric(as.factor(gr[which(useCol[,i])])), B=ROTSn)$pvalue   # ,K=500  
          pVaRotsTab[,,i] <- tmRO } }
      if(debug) message("tRN6")
      
      ## propagate filtering results to p-values (disqualify to NA)
      if(any(!datFi$filt, na.rm=TRUE)) {
        fiAr <- rep(datFi$filt, nLoop)    
        pValTab[which(!fiAr)] <- NA }   
          
      ## resume indiv rounds of imputation, optinal return details 
      out$datImp <- as.matrix(apply(datIm, 1:2, mean, na.rm=TRUE))
      rownames(out$datImp) <- if(is.null(rownames(dat))) rownames(annot) else rownames(dat)
      if("tValTab" %in% multCorMeth) { out$tValArr <- tValTab
        out$pValArr <- pValTab }
      if("noLimma" %in% multCorMeth) out$simple.p.value <- out$p.value
      out$p.value <- as.matrix(apply(pValTab, 1:2, stats::median, na.rm=FALSE))
      out$t <- as.matrix(apply(tValTab, 1:2, stats::median, na.rm=FALSE))
      colnames(out$p.value) <- colnames(out$t) <- rownames(pwComb)
      ## when converting t-value to p how to consider n due to nLoop ??
    } else out$datImp <- datFi$data
    out$annot <- annot
    out$filter <- datFi$filt
    if(debug) message("tRN7")     
  
    ## update dimnames of out$datImp
    dimnames(out$datImp) <- list(if(is.null(rownames(out$lods))) rownames(out$annot) else rownames(out$lods), colnames(dat))
    rownames(out$t) <- rownames(out$p.value) <- rownames(out$datImp)
    ## integrate column specific filtering
    if(any(!datFi$filt, na.rm=TRUE)) out$p.value[which(!datFi$filt)] <- NA 
    if(lfdrInclude) {out$lfdr <- as.matrix(apply(out$p.value, 2, wrMisc::pVal2lfdr, callFrom=fxNa)) 
      dimnames(out$lfdr) <- list(rownames(out$lods), colnames(out$contrasts))
      if("noLimma" %in% multCorMeth) out$simple.lfdr <- if(ncol(out$simple.p.value) >1) apply(out$simple.p.value, 2, wrMisc::pVal2lfdr) else wrMisc::pVal2lfdr(out$simple.p.value)
    } 
    if(any(c("FDR","BH") %in% multCorMeth, na.rm=TRUE)) { out$BH <- as.matrix(apply(out$p.value, 2, stats::p.adjust, method="BH"))
      dimnames(out$BH) <- list(rownames(out$lods), colnames(out$contrasts))
      if("noLimma" %in% multCorMeth) out$simple.BH <- if(ncol(out$simple.p.value) >1) apply(out$simple.p.value, 2, stats::p.adjust, method="BH") else stats::p.adjust(out$simple.p.value, method="BH")
      } 
    if("BY" %in% multCorMeth) {out$BY <- as.matrix(apply(out$p.value, 2, stats::p.adjust, method="BY"))
      dimnames(out$BY) <- list(rownames(out$lods), colnames(out$contrasts))}
    if(length(ROTSn)==1) if(ROTSn >0 & chNA & nLoop >1) {
      out$ROTS.p <- apply(pVaRotsTab, 1:2, stats::median, na.rm=TRUE)
      if(any(!datFi$filt, na.rm=TRUE)) out$ROTS.p[which(!datFi$filt)] <- NA    
      out$ROTS.BH <- as.matrix(as.matrix(apply(out$ROTS.p, 2, stats::p.adjust, method="BH")))
      dimnames(out$ROTS.BH) <- list(rownames(out$lods), colnames(out$contrasts) )
      if(lfdrInclude) {out$ROTS.lfdr <- as.matrix(as.matrix(apply(out$ROTS.p, 2, wrMisc::pVal2lfdr)))
        dimnames(out$ROTS.lfdr) <- list(rownames(out$lods), colnames(out$contrasts))} }
      out 
  } else { warning(fxNa,msg)
    return(NULL) }
  }
   
