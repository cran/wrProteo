#' Test robust to NA-imputation
#'
#' \code{testRobustToNAimputation} replaces \code{NA} values based on group neigbours (based on grouping of columns in argument \code{gr}), following overall assumption of close to Gaussian distribution.
#' Furthermore, it is assumed that \code{NA}-values originate from experimental settings where measurements at or below detection limit are recoreded as \code{NA}.
#' In  such cases (eg in proteomics) it is current practice to replace \code{NA}-values by very low (random) values in order to be able to perform t-tests.
#' However, random normal values used for replacing may in rare cases deviate from the average (the 'assumed' value) and in particular, if multiple \code{NA} replacements are above the average, 
#' may look like induced biological data and be misinterpreted as so.      
#' By repeating multiple times the process of replacing \code{NA}-values and subsequent testing the results can be sumarized afterwards by median over all repeated runs to remmove the stochastic effect of individual NA-imputation.
#' Thus, one may gain stability towards random-character of \code{NA} imputations by repeating imputation & test 'nLoop' times and summarize p-values by median (results stabilized at 50-100 rounds).
#' It is necessary to define all groups of replicates in \code{gr} to obtain all possible pair-wise testing (multiple columns in $BH, $lfdr etc). 
#' Testing by package \href{https://bioconductor.org/packages/release/bioc/html/ROTS.html}{ROTS} may optionaly be included.
#' This function returns limma-like S3 list-object further enriched by additional fields/elements.
#' @param dat (matrix or data.frame) main data (may contain \code{NA}); if \code{dat} is list containing $quant and $annot as matrix, the element $quant will be used
#' @param gr (character or factor) replicate association
#' @param annot (matrix or data.frame) annotation (lines must match lines of data !), if \code{annot} is \code{NULL} and argument \code{dat} is a list containing both $quant and $annot, the element $annot will be used 
#' @param retnNA (logical) retain and report number of \code{NA}
#' @param avSdH (numeric) population characteristics (mean and sd) for >1 \code{NA} neighbours 'high' (per line)
#' @param avSdL (numeric) population characteristics (mean and sd) for >0 \code{NA} neighbours 'low' (per line)
#' @param plotHist (logical) additional histogram of original, imputed and resultant distribution (made using \code{\link{matrixNAneighbourImpute}} )
#' @param xLab (character) custom x-axis label 
#' @param tit (character) custom title
#' @param seedNo (integer) seed-value for normal random values
#' @param nLoop (integer) number of runs of independent \code{NA}-imputation
#' @param lfdrInclude (logical) include lfdr estimations (may cause warning message(s) concerning convergence if few too lines/proteins in dataset tested).
#' @param ROTSn (integer) number of repeats by \code{ROTS}, if \code{NULL} \code{ROTS} will not be called
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of messages produced
#' @return matrix including imputed values or list of final and matrix with number of imputed by group. Various options of multiple testing ccorrection are implemented ('BY','lfdr','p.value' and 'ROTS.BH'
#' @seealso  \code{\link[wrMisc]{moderTest2grp}}, \code{\link[wrMisc]{pVal2lfdr}}, \code{\link[limma]{eBayes}}, \code{\link[stats]{t.test}},\code{\link[ROTS]{ROTS}}  
#' @examples
#' set.seed(2015); rand1 <- round(runif(600)+rnorm(600,1,2),3)
#' dat1 <- matrix(rand1,ncol=6) + matrix(rep((1:100)/20,6),ncol=6)
#' dat1[13:16,1:3] <- dat1[13:16,1:3]+2   # augment lines 13:16 
#' dat1[19:20,1:3] <- dat1[19:20,1:3]+3   # augment lines 19:20
#' dat1[15:18,4:6] <- dat1[15:18,4:6]+1.4   # augment lines 15:18 
#' dat1[dat1 <1] <- NA                    # mimick some NAs for low abundance
#' ## normalize data
#' boxplot(dat1,main="data before normalization")
#' dat1 <- wrMisc::normalizeThis(as.matrix(dat1),meth="median")
#' ## designate replicate relationships in samples ...  
#' grp1 <- gl(2,3,labels=LETTERS[1:2])                   
#' ## moderated t-test with repeated inputations (may take >10 sec,  >60 sec if ROTSn >0 !) 
#' PLtestR1 <- testRobustToNAimputation(dat=dat1,gr=grp1,retnNA=TRUE,nLoop=100,ROTSn=0,lfdr=FALSE)
#' names(PLtestR1)
#' @export
testRobustToNAimputation <- function(dat,gr,annot=NULL,retnNA=TRUE,avSdH=c(0.18,0.5),avSdL=c(0.1,0.5),plotHist=FALSE,xLab=NULL,tit=NULL,seedNo=2018,nLoop=20,lfdrInclude=TRUE,ROTSn=NULL,silent=FALSE,callFrom=NULL){
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="testRobustToNAimputation")
  if(is.list(dat)) { if(all(c("quant","annot") %in% names(dat))) if(length(dim(dat$quant))==2 & length(dim(dat$annot)) ==2) {
    if(length(annot) <1) { annot <- dat$annot} 
    dat <- dat$quant }}
  if(length(dim(dat)) !=2) stop("'dat' must be matrix or data.frame with >1 columns")
  if(is.data.frame(dat)) dat <- as.matrix(dat)
  if(length(gr) != ncol(dat)) stop("Number of columns in 'dat' and number of (group-)elements in 'gr' do not match !")
  if(!is.factor(gr)) gr <- as.factor(gr)
  if(is.null(xLab)) xLab <- "values"            
  if(length(annot) <1) annot <- matrix(NA,nrow=nrow(dat),ncol=1,dimnames=list(rownames(dat),"rowNa"))
  ## main
  isNA <- is.na(dat)
  chNA <- any(isNA)
  nNAmat <- matrix(0,nrow=nrow(dat),ncol=length(levels(gr)),dimnames=list(NULL,levels(gr)))
  seedNo <- as.integer(seedNo)[1]
  gr <- as.factor(gr)
  callFro <- try(as.factor(gr)) 
  if(class(callFro) == "try-error") message("+++++\n",fxNa," MAJOR PROBLEM with argument 'gr' !!  (possibly not sufficient level-names ?) \n+++++")
  ## 1st pass
  datI <- matrixNAneighbourImpute(dat,gr,seedNo=seedNo,retnNA=retnNA,avSdH=avSdH,avSdL=avSdL,plotHist=plotHist,xLab=xLab,tit=tit)
  datFi <- combineMultFilterNAimput(dat=dat,imputed=datI,grp=gr,annDat=annot,abundThr=stats::quantile(dat,0.02,na.rm=TRUE))  # number of unique peptides unknown !
  ## done combineMultFilterNAimput
  ## prepare for testing
  if(lfdrInclude) {
    chPa <- requireNamespace("fdrtool", quietly=TRUE)
    if(!chPa) { message(fxNa,": package 'fdrtool' not installed, omit argument 'lfdrInclude'")
      lfdrInclude <- FALSE }}  
  pwComb <- wrMisc::triCoord(length(levels(gr)))
  out <- wrMisc::moderTestXgrp(datFi$data,grp=gr,limmaOutput=TRUE,addResults="means")        #   ,anno=annot  # can't do question specific filtering w/o explicit loop 
  ## need to add $ROTS.p
  if(length(stats::na.omit(ROTSn))==1) if(ROTSn >0) {  
    chPa <- requireNamespace("ROTS", quietly=TRUE)
    if(!chPa) { message(fxNa,": package 'RORS' not installed, omit argument 'ROTSn'")
      ROTSn <- NULL }}  
  if(length(stats::na.omit(ROTSn))==1) if(ROTSn >0) {
    ## this requires package ROTS
    comp <- wrMisc::triCoord(length(levels(gr)))
    rownames(comp) <- paste(levels(gr)[comp[,1]],levels(gr)[comp[,2]],sep="-")
    tmRO <- matrix(nrow=nrow(datFi$data),ncol=nrow(comp))
    comPair <- matrix(unlist(strsplit(rownames(comp),"-")),ncol=nrow(comp))
    useCol <- apply(comPair,2,function(x) gr %in% x)
    for(i in 1:nrow(comp)) tmRO[which(datFi$filt[,i]),i] <- ROTS::ROTS(datFi$data[which(datFi$filt[,i]),which(useCol[,i])], groups=as.numeric(as.factor(gr[which(useCol[,i])]))-1, B=ROTSn)$pvalue   # ,K=500  
    out$ROTS.p <- tmRO
    out$ROTS.BH <- apply(tmRO,2,stats::p.adjust,method="BH") 
    if(lfdrInclude) out$ROTS.lfdr <- apply(tmRO,2,wrMisc::pVal2lfdr) 
  }   
  ## subsequent rounds of NA-imputation  
  if(chNA & nLoop >1) { 
    pValTab <- tValTab <- array(NA,dim=c(nrow(dat),nrow(pwComb),nLoop))
    datIm <- array(NA,dim=c(nrow(dat),ncol(dat),nLoop))
    datIm[,,1] <- datFi$data
    pValTab[,,1] <- out$p.value
    tValTab[,,1] <- out$t
    if(length(stats::na.omit(ROTSn))==1) if(ROTSn >0) {
      pVaRotsTab <- array(NA,dim=c(nrow(dat),nrow(pwComb),min(10,nLoop)))
      pVaRotsTab[,,1] <- out$ROTS.p}
    for(i in 2:nLoop) {
      datX <- dat
      ## idea 17 oct change seed intiation ?
      datX <- .imputeNA(dat,gr=gr,impParam=datI$randParam +c(0,0,0,i),exclNeg=TRUE)$data          
      fitX <- limma::eBayes(limma::contrasts.fit(limma::lmFit(datX[,],out$design),contrasts=out$contrasts))
      datIm[,,i] <- datX
      pValTab[,,i] <- fitX$p.value
      tValTab[,,i] <- fitX$t
      if(length(stats::na.omit(ROTSn))==1) if(ROTSn >0 & i < min(10,nLoop)) {    # test using ROTS (TAKES MUCH TIME !!)
        for(i in 1:nrow(comp)) tmRO[which(datFi$filt[,i]),i] <- ROTS::ROTS(datFi$data[which(datFi$filt[,i]),which(useCol[,i])], groups=as.numeric(as.factor(gr[which(useCol[,i])])), B=ROTSn)$pvalue   # ,K=500  
        pVaRotsTab[,,i] <- tmRO } }
    out$datImp <- as.matrix(apply(datIm,1:2,mean,na.rm=TRUE))
    out$p.value <- as.matrix(apply(pValTab,1:2,stats::median,na.rm=TRUE))
    out$t.value <- as.matrix(apply(tValTab,1:2,stats::median,na.rm=TRUE))
    ## when converting t-value to p how to consider n due to nLoop ??
  } else out$datImp <- datFi$data
  out$annot <- annot
    ## integrate column specific filtering
  if(any(!datFi$filt)) out$p.value[which(!datFi$filt)] <- NA 
  if(lfdrInclude) {out$lfdr <- as.matrix(apply(out$p.value,2,wrMisc::pVal2lfdr)) 
    colnames(out$lfdr) <- colnames(out$contrasts)} 
  if(length(dat) >0) {out$BH <- as.matrix(apply(out$p.value,2,stats::p.adjust,method="BH"))
    colnames(out$BH) <- colnames(out$contrasts)} 
  if(length(dat) >0) {out$BY <- as.matrix(apply(out$p.value,2,stats::p.adjust,method="BY"))
    colnames(out$BY) <- colnames(out$contrasts)} 
  if(length(stats::na.omit(ROTSn))==1) if(ROTSn >0 & chNA & nLoop >1){
    out$ROTS.p <- apply(pVaRotsTab,1:2,stats::median,na.rm=TRUE)
    if(any(!datFi$filt)) out$ROTS.p[which(!datFi$filt)] <- NA    
    out$ROTS.BH <- as.matrix(as.matrix(apply(out$ROTS.p,2,stats::p.adjust,method="BH")))
    colnames(out$ROTS.BH) <- colnames(out$contrasts) 
    if(lfdrInclude) {out$ROTS.lfdr <- as.matrix(as.matrix(apply(out$ROTS.p,2,wrMisc::pVal2lfdr)))
      colnames(out$ROTS.lfdr) <-colnames(out$contrasts)} }
  out }

