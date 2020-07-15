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
#' @param thr (numeric) threshold, if \code{NULL} a panel of 108 values will be used for calculating specifcity and sensitivity 
#' @param tyThr (character,length=1) type of test-result to be used for sensitivity and specificity calculations (eg 'BH','lfdr' or 'p.value'), must be list-element of 'test'
#' @param columnTest (character or integer) only in case 'tyThr' is matrix (as typically the case after \code{testRobustToNAimputation}) : which column of 'test$tyThr' should be used as test-result 
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
#' @seealso replot the figure \code{\link[wrProteo]{plotROC}}, robust test for preparing tables \code{\link{testRobustToNAimputation}}, \code{\link[wrMisc]{moderTest2grp}}, \code{\link{test2grp}}, \code{eBayes} in package \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{limma}, \code{\link[stats]{t.test}}  
#' @examples
#' set.seed(2019); test1 <- list(annot=cbind(spec=c(rep("b",35),letters[sample.int(n=3,
#'   size=150,replace=TRUE)])),BH=matrix(c(runif(35,0,0.01),runif(150)),ncol=1))
#' tail(roc1 <- summarizeForROC(test1,spec=c("a","b","c")))
#' 
#' @export
summarizeForROC <- function(test,thr=NULL,tyThr="BH",columnTest=1,spec=c("H","E","S"),annotCol="spec",tit=NULL,color=1,plotROC=TRUE,pch=1,bg=NULL,overlPlot=FALSE,silent=FALSE,callFrom=NULL) {
  ## summarize esting result by species (3rd is supposed as reference)
  argN <- deparse(substitute(test))
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="summarizeForROC")
  chLst <- tyThr %in% names(test)
  if(!any(chLst)) stop("Don't know what kind of test-results to use.  Can't find element '",tyThr[!chLst],"' in elements of 'test' !!")
  if(!"annot" %in% names(test)) stop("test$annot is needed to map content of 'spec'")
  if(is.null(thr)) thr <- signif(c(as.numeric(sapply((1:4)*2,function(x) x*c(1e-4,1e-5,1e-6,1e-7))),
    seq(0,1,length.out=50)^5,seq(0,1,length.out=50)^2,seq(0,1,length.out=61),4^(-2:-10),1.01),2)
  thr <- sort(unique(abs(wrMisc::naOmit(thr))))                                    # 151 -> 108 values for default
  pp <- matrix(nrow=length(thr),ncol=4,dimnames=list(NULL,c("TP","FP","FN","TN")))
  spiSpec <- if(ncol(test$annot) >1) {test$annot[,annotCol] %in% spec[-1]} else {test$annot[,1] %in% spec[-1]} 
  if(length(dim(test[[tyThr]])) ==2) if(ncol(test[[tyThr]])==2 & identical(colnames(test[[tyThr]])[1],"(Intercept)") & columnTest !=2) {
    message(fxNa," Value of argument 'columnTest' seems bizzare, setting to =2 !  (to avoid testing '(Intercept)')")
    columnTest <- 2 }
  tmp1 <- if(length(dim(test[[tyThr]])) ==2) test[[tyThr]][which(spiSpec),columnTest] else test[[tyThr]][which(spiSpec)]
  tmp2 <- if(length(dim(test[[tyThr]])) ==2) test[[tyThr]][which(!spiSpec),columnTest] else test[[tyThr]][which(!spiSpec)]
  pp <- cbind(
    TP= sapply(thr,function(x) sum(tmp1 <=x,na.rm=TRUE)),
    FP= sapply(thr,function(x) sum(tmp2 <=x,na.rm=TRUE)),
    FN= sapply(thr,function(x) sum(tmp1 >x,na.rm=TRUE)),
    TN= sapply(thr,function(x) sum(tmp2 >x,na.rm=TRUE)))
  keyVal <- cbind(alph=thr,
    spec=as.numeric(pp[,"TN"]/(pp[,"TN"]+pp[,"FP"])), sens=as.numeric(c(pp[,"TP"]/(pp[,"TP"]+pp[,"FN"]))),
    prec=as.numeric(pp[,"TP"]/(pp[,"FP"]+pp[,"TP"])), accur=as.numeric((pp[,"TP"]+pp[,"TN"])/rowSums(pp)), 
    FDR=as.numeric(pp[,"FP"]/(pp[,"TP"]+pp[,"FP"])))
  ## add no of lines/prot retained for each species
  tmp3 <- test$annot[,annotCol] %in% spec[2]
  tmp3 <- if(length(dim(test[[tyThr]])) >1) test[[tyThr]][which(tmp3),columnTest] else test[[tyThr]][which(tmp3)]       # pvalues for 1st spike-in species (eg E)
  if(length(spec) >2) {
    tmp4 <- test$annot[,annotCol] %in% spec[3]
    tmp4 <- if(length(dim(test[[tyThr]])) >1) test[[tyThr]][which(tmp4),columnTest] else test[[tyThr]][which(tmp4)] }      # pvalues for 2nd spike-in species (eg S)
  tmp3 <- cbind(Sp1Pos=pp[,"FP"], Sp2Pos=sapply(thr,function(x) sum(tmp3 <=x,na.rm=TRUE)),
    Sp3Pos=if(length(spec) >2) sapply(thr,function(x) sum(tmp4 <=x,na.rm=TRUE)) else NULL)
  colnames(tmp3) <- paste("n.pos",spec,sep=".")  
  keyVal <- cbind(keyVal,tmp3)  
  if(plotROC) {
    if(is.null(tit)) tit <- paste("ROC of ",argN)
    xLab <- "1 - Specificity"
    yLab <- "Sensitivity"
    if(overlPlot) graphics::points(1-keyVal[,annotCol],keyVal[,"sens"],col=color,pch=pch,bg=NULL,type="S") else {
      graphics::plot(1-keyVal[,annotCol],keyVal[,"sens"],col=color,pch=pch,bg=bg,type="S",main=tit,xlab=xLab,ylab=yLab,xlim=c(0,1),ylim=c(0,1))}
    cutP <- keyVal[which(keyVal[,1]==0.05),-1]
    newPch <- cbind(c(1,16,2,17, 7,15,5,6),new=c(21,21,24,24,22,22,23,25))                                   # transform open or plain filled points to color-filled
    if(pch %in% newPch[,1]){ pch2 <- newPch[which(newPch[,1]==pch),2]; bg <- color; col2 <- grDevices::grey(0.2)} else {pch2 <- pch; col2 <- color}   
    col2 <- wrMisc::convColorToTransp(col2,alph=0.9)
    graphics::points(1-cutP[annotCol],cutP["sens"],col=col2,pch=pch2,bg=bg,cex=1.4)  
  }
  keyVal }     
     
