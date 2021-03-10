#' Summarize statistical test result for plotting ROC-curves
#'
#' \code{summarizeForROC} takes statistical testing results (obtained using \code{\link{testRobustToNAimputation}} or \code{\link[wrMisc]{moderTest2grp}}, based on \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{limma}) and calculates specifcity and sensitivity values for plotting ROC-curves along a panel of thresholds.  
#' Based on column from test$annot and argument 'spec' TP,FP,FN and TN are determined. Special consideration is made to 3 species mix samples as found in proteomics benchmark-tests.  
#' See also \href{https://en.wikipedia.org/wiki/Receiver_operating_characteristic}{ROC on Wkipedia} for explanations of TP,FP,FN and TN as well as examples.
#' An optional plot may be produced, too.
#' Return matrix with TP,FP,FN,TN,spec,sens,prec,accur and FDR count values along the various thrsholds specified in column 'alph'.
#' Note that numerous other packages also provide support for building and plotting ROC-curves : Eg \href{https://CRAN.R-project.org/package=dlstats}{rocPkgShort}, 
#'  \href{https://CRAN.R-project.org/package=ROCR}{ROCR}, \href{https://CRAN.R-project.org/package=pROC}{pROC} or \href{https://CRAN.R-project.org/package=ROCit}{ROCit} 
#'  
#' @param test (class \code{MArrayLM}, S3-object from limma) from testing (eg \code{\link{testRobustToNAimputation}} or \code{\link{test2grp}}
#' @param thr (numeric) p-value threshold, if \code{NULL} a panel of 108 p-value threshold-levels values will be used for calculating specifcity and sensitivity 
#' @param tyThr (character,length=1) type of test-result to be used for sensitivity and specificity calculations (eg 'BH','lfdr' or 'p.value'), must be list-element of 'test'
#' @param useComp (character or integer) in case multiple comparisons (ie multiple columns 'test$tyThr'); which pairwise comparison to used
#' @param columnTest depreciated, pleas use 'useComp' instead
#' @param FCthrs (numeric) Fold-Change threshold (display as line) give as Fold-change and NOT log2(FC), default at 1.5, set to \code{NA} for omitting


#' @param spec (character) labels for species will be matched to column 'spec' of test$annot and used for sensitivity and specificity calculations. Important : 1st label for matrix (expected as constant) and subsequent labels for spike-ins (variable)
#' @param annotCol (character) column name of \code{test$annot} to use to separate species
#' @param tit (character) optinal custom title in graph 
#' @param color (character or integer) color in graph
#' @param plotROC (logical) toogle plot on or off  
#' @param pch (integer) type of symbol to be used (see \code{\link[graphics]{par}})
#' @param bg (character) backgroud in plot (see \code{\link[graphics]{par}})
#' @param overlPlot (logical) overlay to existing plot if \code{TRUE} 
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of message(s) produced
#' @return matrix including imputed values or list of final and matrix with number of imputed by group (plus optional plot)
#' @seealso replot the figure \code{\link{plotROC}}, calculate AUC using \code{\link{AucROC}}, robust test for preparing tables \code{\link{testRobustToNAimputation}}, \code{\link[wrMisc]{moderTest2grp}}, \code{\link{test2grp}}, \code{eBayes} in package \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{limma}, \code{\link[stats]{t.test}}  
#' @examples
#' set.seed(2019); test1 <- list(annot=cbind(spec=c(rep("b",35),letters[sample.int(n=3,
#'   size=150,replace=TRUE)])),BH=matrix(c(runif(35,0,0.01),runif(150)),ncol=1))
#' tail(roc1 <- summarizeForROC(test1,spec=c("a","b","c")))
#' 
#' @export
summarizeForROC <- function(test,thr=NULL,tyThr="BH",useComp=1,columnTest=NULL,FCthrs=NULL,spec=c("H","E","S"),annotCol="spec",tit=NULL,color=1,plotROC=TRUE,pch=1,bg=NULL,overlPlot=FALSE,silent=FALSE,callFrom=NULL) {
  ## summarize esting result by species (3rd is supposed as reference)
  argN <- deparse(substitute(test))
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="summarizeForROC")
  inclFilter <- TRUE     # use $filtFin
  badFCtoNA <- FALSE     # how to disqualify FC not passing filter
  
  chLst <- tyThr %in% names(test)
  if(any(!chLst)) stop("Don't know what kind of test-results to use.  Can't find element '",tyThr[!chLst],"' in elements of 'test' !!")
  if(!"annot" %in% names(test)) stop("test$annot is needed to map content of 'spec'")
  if(is.null(thr)) thr <- signif(c(as.numeric(sapply((1:4)*2, function(x) x*c(1e-4,1e-5,1e-6,1e-7))),
    seq(0,1,length.out=50)^5,seq(0,1,length.out=50)^2,seq(0,1,length.out=61),4^(-2:-10),1.01),2)
  thr <- sort(unique(abs(wrMisc::naOmit(thr))))                                    # 151 -> 108 values for default
  pp <- matrix(nrow=length(thr), ncol=4, dimnames=list(NULL,c("TP","FP","FN","TN")))
  oriKeep <- matrix(rep(TRUE,nrow(test[[tyThr]])), nrow=nrow(test[[tyThr]]), ncol=2,dimnames=list(NULL,c("filtFin","passFC")))     # the lines passing various filtering (filtFin, FCthrs)  
  chSpec <- spec %in% test$annot[,annotCol]
  if(all(!chSpec)) stop(" None of the elements of argument 'spec' found in column '",annotCol,"' !!")
  if(any(!chSpec)) message(fxNa," Trouble ahead ?  Elements ",wrMisc::pasteC(spec[which(!chSpec)],quoteC="'")," not found !!")  
  ## look for (global) filtering (from test$filtFin)
  if("filtFin" %in% names(test) & inclFilter) {
    chLe <- if(length(dim(test$filtFin)) >1) nrow(test$filtFin) else length(test$filtFin) ==nrow(test[[tyThr]]) 
    if(chLe) {
      if(length(dim(test$filtFin)) >1) test$filtFin <- test$filtFin[,useComp]
      if(sum(test$filtFin) <1) {test$filtFin <- rep(TRUE,length(test$filtFin)); message(fxNa," nothing passing filtering, ignoring filter !")
      } else if(!silent) message(fxNa," filtering: ",sum(test$filtFin)," out of ",nrow(test$filtFin)," pass filtering")      
      ## apply filtFin
      if(any(!test$filtFin)) {oriKeep[which(!test$filtFin),] <- FALSE
       test[[tyThr]][which(!test$filtFin),useComp] <- NA }
    } else if(!silent) message(fxNa,"  ",argN,"$filtFin does not have same length as ",argN,"$",tyThr," , ignoring filter")
  }  
  ## FC-filtering
  if(length(FCthrs) ==1) if(is.numeric(FCthrs) & !is.na(FCthrs)) {
    ## FC-threshold, need to locate means to construct FC
    chM <- "means" %in% names(test)
    if("means" %in% names(test)) {
      if(nrow(test$means) ==nrow(test[[tyThr]])) {
        ## identify sample-groups to comparison(s) - needed lateron
        pairwCol <- .sampNoDeMArrayLM2(test, useComp, lstMeans="means", lstP=which(names(test)==tyThr),callFrom=fxNa,silent=silent) 
        grpMeans <- cbind(mean1=test$means[,pairwCol[1]], mean2=test$means[,pairwCol[2]])
        FCval <- grpMeans[,2] - grpMeans[,1]  
        ## FC-filter
        chFC <- abs(FCval) >= log2(FCthrs)
        if(any(!chFC))  {
          oriKeep[which(!chFC),2] <- FALSE    # update oriKeep
          ## disqualify FDR for those not passing FCthrs as FDR=1.0 so they won't get counted (appear only at end)
          if(is.logical(badFCtoNA)) test[[tyThr]][which(!chFC ),useComp] <- if(badFCtoNA) NA else 1
          ## need also to explore other ways of dynamic combining FCthrs to FDRthrs
          
        }
        if(!silent) message(fxNa," FC-filter: ",sum(chFC)," out of ",length(FCval)," pass 'FCthrs'=",FCthrs )
      } else warning(argN,"$means has not correct number of rows, ignoring") 
    } else warning("Could not find suitable field '$means' in '",argN,"'")    
  } else { FCval <- NULL }   # needed ?
  if(sum(oriKeep[,2]) <2) warning("TROUBLE AHEAD: Only ",sum(oriKeep)," out of ",length(oriKeep)," elements pass filtering")

  spiSpec <- if(ncol(test$annot) >1) {test$annot[,annotCol] %in% spec[-1]} else {test$annot[,1] %in% spec[-1]} # locate who is not 1st/ref species
  spiSpec[which(!oriKeep)] <- NA
  if(all(!spiSpec)) message(fxNa," Trouble ahead ?  Could not find any element annotated as species to search for (ie, TP will always remain 0)")
  if(length(dim(test[[tyThr]])) ==2) if(ncol(test[[tyThr]])==2 & identical(colnames(test[[tyThr]])[1],"(Intercept)") & useComp !=2) {
    ## single comparison, thus 2 cols of p-val => use 2nd
    useComp <- 2 }
  
  tmp1 <- if(length(dim(test[[tyThr]])) ==2) test[[tyThr]][which(spiSpec),useComp] else test[[tyThr]][which(spiSpec)]
  tmp2 <- if(length(dim(test[[tyThr]])) ==2) test[[tyThr]][which(!spiSpec),useComp] else test[[tyThr]][which(!spiSpec)]
  if(all(is.na(tmp1))) {
    message(fxNa," PROBLEM :\n  ***  None of the elements annotated as 'positive' species to search for has any valid testing results ! Unable to construct TP ! ***")
    keyVal <- NULL
  } else {
    pp <- cbind(
      TP= sapply(thr,function(x) sum(tmp1 <=x,na.rm=TRUE)),
      FP= sapply(thr,function(x) sum(tmp2 <=x,na.rm=TRUE)),
      FN= sapply(thr,function(x) sum(tmp1 >x,na.rm=TRUE)),
      TN= sapply(thr,function(x) sum(tmp2 >x,na.rm=TRUE)))
  
    keyVal <- cbind(alph=thr,
      spec=as.numeric(pp[,"TN"]/(pp[,"TN"]+pp[,"FP"])), sens=as.numeric(c(pp[,"TP"]/(pp[,"TP"]+pp[,"FN"]))),
      prec=as.numeric(pp[,"TP"]/(pp[,"FP"]+pp[,"TP"])), accur=as.numeric((pp[,"TP"]+pp[,"TN"])/rowSums(pp)), 
      FDR=as.numeric(pp[,"FP"]/(pp[,"TP"]+pp[,"FP"])))
    chNaN <- colSums(is.nan(keyVal)) ==nrow(keyVal)
    if(any(chNaN)) keyVal[,which(chNaN)] <- 0
    ## add no of lines/prot retained for each species
    tmp3 <- test$annot[,annotCol] %in% spec[2]
    tmp3 <- if(length(dim(test[[tyThr]])) >1) test[[tyThr]][which(tmp3),useComp] else test[[tyThr]][which(tmp3)]        # pvalues for 1st spike-in species (eg E)
  
    if(length(spec) >2) {
      tmp4 <- test$annot[,annotCol] %in% spec[3]
      tmp4 <- if(length(dim(test[[tyThr]])) >1) test[[tyThr]][which(tmp4),useComp] else test[[tyThr]][which(tmp4)] }    # pvalues for 2nd spike-in species (eg S)
    tmp3 <- cbind(Sp1Pos=pp[,"FP"], Sp2Pos=sapply(thr, function(x) sum(tmp3 <=x,na.rm=TRUE)),
      Sp3Pos=if(length(spec) >2) sapply(thr,function(x) sum(tmp4 <=x,na.rm=TRUE)) else NULL)
    colnames(tmp3) <- paste("n.pos",spec,sep=".")  
    keyVal <- cbind(keyVal,tmp3)  
    if(plotROC) {
      if(is.null(tit)) tit <- paste("ROC of ",argN)
      xLab <- "1 - Specificity"
      yLab <- "Sensitivity"
      if(overlPlot) graphics::points(1-keyVal[,"spec"], keyVal[,"sens"], col=color,pch=pch,bg=NULL,type="s") else {
        graphics::plot(1-keyVal[,"spec"], keyVal[,"sens"], col=color,pch=pch,bg=bg,type="s",main=tit,xlab=xLab,ylab=yLab,xlim=c(0,1), ylim=c(0,1), las=1)}
      cutP <- keyVal[which(keyVal[,1]==0.05),-1]
      newPch <- cbind(c(1,16,2,17, 7,15,5,6), new=c(21,21,24,24,22,22,23,25))                                   # transform open or plain filled points to color-filled
      if(pch %in% newPch[,1]){ pch2 <- newPch[which(newPch[,1]==pch),2]; bg <- color; col2 <- grDevices::grey(0.2)} else {pch2 <- pch; col2 <- color}   
      col2 <- wrMisc::convColorToTransp(col2, alph=0.9)
      graphics::points(1-cutP["spec"],cutP["sens"], col=col2, pch=pch2, bg=bg, cex=1.4) 
      graphics::mtext(paste("AUC = ",signif(AucROC(keyVal),3)), side=3, cex=0.8,adj=0)
  } }
  keyVal }     
    
