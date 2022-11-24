#' Summarize statistical test result for plotting ROC-curves
#'
#' This function takes statistical testing results (obtained using \code{\link{testRobustToNAimputation}} or \code{\link[wrMisc]{moderTest2grp}}, 
#' based on \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{limma}) and calculates specifcity and sensitivity values for plotting ROC-curves along a panel of thresholds.  
#' Based on annotation (from test$annot) with the user-defined column for species (argument 'spec') the counts of TP (true positives), FP (false positves), FN (false negatives) and TN are determined. 
#' In addition, an optional plot may be produced.
#'
#' @details
#' Determining TP and FP counts requires 'ground trouth' experiments, where it is known in advance which proteins are expected to change abundance between two groups of samples.
#' Typically this is done by mixing proteins of different species origin, the first species noted by argument 'spec' designes the species to be considered constant (expected as FN in statistical tests).
#' Then, one or mutiple additional spike-in species can be defined. As the spike-in cocentration should have been altered between different gruops of samples, they are expected as TP. 
#'  
#' The main aim of this function consists in providing specifcity and sensitivity values, plus counts of TP (true positives), FP (false positves), FN (false negatives) and TN (true negatives),  
#' along various thrsholds (specified in column 'alph') for statistical tests preformed prior to calling this function.
#'  
#' Note, that the choice of species-annotation plays a crucial role who the counting results are obtained. 
#' In case of multiple spike-in species the user should pay attention if they all are expected to change abundance at the same ratio.
#' If not, it is advised to run this function multiple times sperately only with the subset of those species expected to change at same ratio.
#'   
#' The dot on the plotted curve shows the results at the level of the single threshold alpha=0.05.   
#' For plotting multiple ROC curves as overlay and additional graphical parameters/options you may use \code{\link{plotROC}}.
#'  
#' See also \href{https://en.wikipedia.org/wiki/Receiver_operating_characteristic}{ROC on Wkipedia} for explanations of TP,FP,FN and TN as well as examples.
#' Note that numerous other packages also provide support for building and plotting ROC-curves : Eg \href{https://CRAN.R-project.org/package=dlstats}{rocPkgShort}, 
#'  \href{https://CRAN.R-project.org/package=ROCR}{ROCR}, \href{https://CRAN.R-project.org/package=pROC}{pROC} or \href{https://CRAN.R-project.org/package=ROCit}{ROCit} 
#'  
#'  
#' @param test (list or class \code{MArrayLM}, S3-object from limma) from testing (eg \code{\link{testRobustToNAimputation}} or \code{\link{test2grp}}
#' @param useComp (character or integer) in case multiple comparisons (ie multiple columns 'test$tyThr'); which pairwise comparison to used
#' @param tyThr (character,length=1) type of statistical test-result to be used for sensitivity and specificity calculations (eg 'BH','lfdr' or 'p.value'), must be list-element of 'test'
#' @param thr (numeric) stat test (FDR/p-value) threshold, if \code{NULL} a panel of 108 p-value threshold-levels values will be used for calculating specifcity and sensitivity 
#' @param columnTest depreciated, please use 'useComp' instead
#' @param FCthrs (numeric) Fold-Change threshold (display as line) give as Fold-change and NOT as log2(FC), default at 1.5, set to \code{NA} for omitting
#' @param spec (character) labels for those species which should be matched to column \code{annotCol} ('spec') of test$annot and used for sensitivity and specificity calculations. Important : 1st entry for species designed as constant (ie matrix) and subsequent labels for spike-ins (expected variable)
#' @param annotCol (character, length=1) column name of \code{test$annot} to use to separate species
#' @param filterMat (character) name (or index) of element of \code{test} containing matrix or vector of logical filtering results
#' @param tit (character) optinal custom title in graph 
#' @param color (character or integer) color in graph
#' @param plotROC (logical) toogle plot on or off  
#' @param pch (integer) type of symbol to be used (see \code{\link[graphics]{par}})
#' @param bg (character) backgroud in plot (see \code{\link[graphics]{par}})
#' @param overlPlot (logical) overlay to existing plot if \code{TRUE} 
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging 
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This function returns a numeric matrix containing the columns 'alph', 'spec', 'sens', 'prec', 'accur', 'FD' plus two columns with absolute numbers of lines (genes/proteins) passing the current threshold level alpha (1st species, all other species) 
#' @seealso replot the figure using \code{\link{plotROC}}, calculate AUC using \code{\link{AucROC}}, robust test for preparing tables \code{\link{testRobustToNAimputation}}, \code{\link[wrMisc]{moderTest2grp}}, \code{\link{test2grp}}, \code{eBayes} in package \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{limma}, \code{\link[stats]{t.test}}  
#' @examples
#' set.seed(2019); test1 <- list(annot=cbind(Species=c(rep("b",35), letters[sample.int(n=3,
#'   size=150, replace=TRUE)])), BH=matrix(c(runif(35,0,0.01), runif(150)), ncol=1))
#' tail(roc1 <- summarizeForROC(test1, spec=c("a","b","c"), annotCol="Species"))
#' 
#' @export
summarizeForROC <- function(test, useComp=1, tyThr="BH", thr=NULL, columnTest=NULL, FCthrs=NULL, spec=c("H","E","S"), annotCol="Species", filterMat="filter", 
  tit=NULL, color=1, plotROC=TRUE, pch=1, bg=NULL, overlPlot=FALSE, silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## summarize esting result by species (3rd is supposed as reference)
  argN <- deparse(substitute(test))
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="summarizeForROC")
  inclFilter <- TRUE     # use $filtFin
  badFCtoNA <- FALSE     # how to disqualify FC not passing filter
  ## checking
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  if(!isFALSE(plotROC)) plotROC <- TRUE
  if(!isTRUE(overlPlot)) overlPlot  <- FALSE
  
  if(!"annot" %in% names(test)) stop("test$annot is needed to map content of 'spec'")
  chLst <- tyThr %in% names(test)
  if(any(!chLst, na.rm=TRUE)) stop("Don't know what kind of test-results to use.  Can't find element '",tyThr[!chLst],"' in elements of 'test' !!")
  if(length(dim(test[[tyThr]])) !=2) stop("Can't find data.frame or matrix as list-elment '",annotCol,"' in 'test'")
  if(length(annotCol) <1) stop("Argument 'annotCol' must be specified (should be name of column of 'test$annot' to check for species from argument 'spec')")
  if(length(spec) <1) stop("Argument 'spec' must be specified (should contain species names present in tst$annot$",annotCol,")")
  if(is.numeric(useComp)) { msg <- "Argument 'useComp' points to comparison not existing, resetting to 1"
    if(useComp > ncol(test[[tyThr]])) { warning(msg); useComp <- 1 } 
  } else if(!(useComp %in% test[[tyThr]])) { warning(msg); useComp <- 1 } 
  if(debug) {message(fxNa,"sROC0"); sROC0 <- list(test=test, chLst=chLst,tyThr=tyThr,thr=thr,useComp=useComp,annotCol=annotCol,spec=spec)}
  if(length(annotCol) >1) { annotCol <- annotCol[1]
    if(!silent) message(fxNa,"NOTE : Argument 'annotCol' may not be longer than 1, using only 1st")}
  
  if(!(annotCol %in% colnames(test$annot))) stop("Can't find column '",annotCol,"' from argument 'annotCol' in test$annot")  
  chSpec <- spec %in% test$annot[,annotCol]
  if(!any(chSpec, na.rm=TRUE)) stop("None of the elements of argument 'spec' found in column '",annotCol,"' !!")
  if(any(!chSpec, na.rm=TRUE)) message(fxNa,"Note : Species-types ",wrMisc::pasteC(unique(spec[which(!chSpec)]),quoteC="'")," will be ignored")  # check for species tags not specified, ie ignored
  if(is.null(thr)) thr <- signif(c(as.numeric(sapply((1:4)*2, function(x) x*c(1e-4,1e-5,1e-6,1e-7))),
    seq(0,1,length.out=50)^5, seq(0,1,length.out=50)^2, seq(0,1,length.out=61), 4^(-2:-10),1.01), 2)
  thr <- sort(unique(abs(wrMisc::naOmit(thr))))                                    # 151 -> 108 values for default
  pp <- matrix(nrow=length(thr), ncol=4, dimnames=list(NULL,c("TP","FP","FN","TN")))
  oriKeep <- matrix(rep(TRUE,nrow(test[[tyThr]])), nrow=nrow(test[[tyThr]]), ncol=2, dimnames=list(NULL,c("filtFin","passFC")))     # the lines passing various filtering (filtFin, FCthrs) 
  if(debug) {message(fxNa,"sROC1")}
  
  ## look for (global) filtering (from test[[filterMat]])
  if(length(filterMat) ==1) if(filterMat %in% names(test) & inclFilter) {
    chLe <- (if(length(dim(test[[filterMat]])) >1) nrow(test[[filterMat]]) else length(test[[filterMat]])) ==nrow(test[[tyThr]])   # check if filtFin seems to fit to testing results 
    if(chLe) {
      if(length(dim(test[[filterMat]])) >1) test[[filterMat]] <- test[[filterMat]][,useComp]
      if(sum(test[[filterMat]]) <1) {test[[filterMat]] <- rep(TRUE,length(test[[filterMat]])); message(fxNa,"Nothing passing filtering, ignoring filter !")
      } else if(!silent) message(fxNa,"Filtering: ",sum(test[[filterMat]])," out of ",length(test[[filterMat]])," pass filtering")      
      ## apply filtFin
      if(any(!test[[filterMat]], na.rm=TRUE)) { oriKeep[which(!test[[filterMat]]),] <- FALSE
        test[[tyThr]][which(!test[[filterMat]]),useComp] <- NA }
      if(debug) {
        message(fxNa,"Filter away ",sum(is.na(test[[tyThr]]))," instances (",round(100*sum(is.na(test[[tyThr]]))/prod(dim(test[[tyThr]])),1),"%);  sROC1b")}
    } else if(!silent) message(fxNa,"  ",argN,"[[",filterMat,"]] does not have same length as ",argN,"$",tyThr," , ignoring filter")
  }  
  if(debug) {message(fxNa,"sROC2"); sROC2 <- list(test=test, chLst=chLst,tyThr=tyThr,thr=thr,useComp=useComp,chSpec=chSpec,pp=pp,oriKeep=oriKeep,filterMat=filterMat,inclFilter=inclFilter)}
  
  ## FC-filtering
  if(length(FCthrs) ==1) if(is.numeric(FCthrs) & !is.na(FCthrs)) {
    ## FC-threshold, need to locate means to construct FC
    chM <- "means" %in% names(test)
    if("means" %in% names(test)) {
      if(nrow(test$means) ==nrow(test[[tyThr]])) {
        ## identify sample-groups to comparison(s) - needed lateron
        pairwCol <- wrMisc::sampNoDeMArrayLM(test, useComp, lstMeans="means", lstP=which(names(test)==tyThr),callFrom=fxNa,silent=silent) 
        grpMeans <- cbind(mean1=test$means[,pairwCol[1]], mean2=test$means[,pairwCol[2]])
        FCval <- grpMeans[,2] - grpMeans[,1]  
        ## FC-filter
        chFC <- abs(FCval) >= log2(FCthrs)
        if(any(!chFC, na.rm=TRUE))  {
          oriKeep[which(!chFC),2] <- FALSE    # update oriKeep
          ## disqualify FDR for those not passing FCthrs as FDR=1.0 so they won't get counted (appear only at end)
          if(is.logical(badFCtoNA)) test[[tyThr]][which(!chFC ),useComp] <- if(badFCtoNA) NA else 1
          ## need also to explore other ways of dynamic combining FCthrs to FDRthrs          
        }
        if(debug) {message(fxNa,"sROC2b")}
  
        if(!silent) message(fxNa,"FC-filter: ",sum(chFC)," out of ",length(FCval)," pass 'FCthrs'=",FCthrs )
      } else warning(argN,"$means has not correct number of rows, ignoring") 
    } else warning("Could not find suitable field '$means' in '",argN,"'")    
  } else { FCval <- NULL }   # needed ?
  if(sum(oriKeep[,2]) <2) warning("TROUBLE AHEAD:  ",sum(oriKeep)," out of ",length(oriKeep)," elements pass filtering")

  spiSpec <- if(ncol(test$annot) >1) {test$annot[,annotCol] %in% spec[-1]} else {test$annot[,1] %in% spec[-1]} # locate who NOT is 1st species (ie Spike and NOT matrix)
  matrSpec <- if(ncol(test$annot) >1) {test$annot[,annotCol] %in% spec[1]} else {test$annot[,1] %in% spec[1]}  # locate who 1st/ref species  (ie matrix)
  if(debug) {message(fxNa,"sROC4")}
  
  if(any(!oriKeep[,2], na.rm=TRUE)) spiSpec[which(!oriKeep[,2])] <- matrSpec[which(!oriKeep[,2])] <- NA
  msg2 <- " Trouble ahead ?  Could not find any element annotated as "
  if(!any(wrMisc::naOmit(matrSpec), na.rm=TRUE)) message(fxNa, msg2, "matrix species to search for (ie, TN will always remain 0)")
  if(!any(wrMisc::naOmit(spiSpec), na.rm=TRUE)) message(fxNa, msg2, "spike-species to search for (ie, TP will always remain 0)")
  if(length(dim(test[[tyThr]])) ==2) if(ncol(test[[tyThr]])==2 & identical(colnames(test[[tyThr]])[1],"(Intercept)") & useComp !=2) {
    ## single comparison, thus 2 cols of p-val => use 2nd
    useComp <- 2 }
  
  teSpi <- if(length(dim(test[[tyThr]])) ==2) test[[tyThr]][which(spiSpec & oriKeep[,2]), useComp] else test[[tyThr]][which(spiSpec & oriKeep[,2])]    # spike,ie non-ref species
  teMat <- if(length(dim(test[[tyThr]])) ==2) test[[tyThr]][which(matrSpec & oriKeep[,2]), useComp] else test[[tyThr]][which(matrSpec & oriKeep[,2])]  # ref-species
  if(debug) {message(fxNa,"sROC5")}

  if(all(is.na(teSpi))) {
    message(fxNa," PROBLEM :\n  ***  None of the elements annotated as 'positive' species to search for has any valid testing results ! Unable to construct TP ! ***")
    keyVal <- NULL
  } else {
    pp <- cbind(
      TP= sapply(thr, function(x) sum(teSpi <=x, na.rm=TRUE)),
      FP= sapply(thr, function(x) sum(teMat <=x, na.rm=TRUE)),
      FN= sapply(thr, function(x) sum(teSpi >x, na.rm=TRUE)),
      TN= sapply(thr, function(x) sum(teMat >x, na.rm=TRUE)) )
    rownames(pp) <- thr   
  
    keyVal <- cbind(alph=thr,
      spec=as.numeric(pp[,"TN"]/(pp[,"TN"] +pp[,"FP"])), sens=as.numeric(c(pp[,"TP"]/(pp[,"TP"] +pp[,"FN"]))),
      prec=as.numeric(pp[,"TP"]/(pp[,"FP"] +pp[,"TP"])), accur=as.numeric((pp[,"TP"] +pp[,"TN"])/rowSums(pp)), 
      FDR=as.numeric(pp[,"FP"]/(pp[,"TP"] +pp[,"FP"])))
    chNaN <- colSums(is.nan(keyVal)) ==nrow(keyVal)
    if(any(chNaN, na.rm=TRUE)) keyVal[,which(chNaN)] <- 0

    ## add no of lines/prot retained for each species
    tmp3 <- test$annot[,annotCol] %in% spec[2:length(spec)]
    tmp3 <- if(length(dim(test[[tyThr]])) >1) test[[tyThr]][which(tmp3),useComp] else test[[tyThr]][which(tmp3)]        # pvalues for 1st spike-in species (eg E)

    if(debug) {message(fxNa,"sROC5b")}
  
    tmp3 <- cbind(Sp1Pos=pp[,"FP"], Sp2Pos=sapply(thr, function(x) sum(tmp3 <=x,na.rm=TRUE)))
    colnames(tmp3) <- paste("n.pos",c(spec[1], if(length(spec) >2) paste(spec[-1],collapse="+") else spec[-1]), sep=".")  
    keyVal <- cbind(keyVal, tmp3)  

    if(plotROC) {
      if(is.null(tit)) tit <- paste("ROC of ",argN)
      xLab <- "1 - Specificity"
      yLab <- "Sensitivity"
      if(overlPlot) graphics::points(1-keyVal[,"spec"], keyVal[,"sens"], col=color, pch=pch, bg=NULL, type="s") else {
        graphics::plot(1-keyVal[,"spec"], keyVal[,"sens"], col=color, pch=pch, bg=bg, type="s",main=tit,xlab=xLab,ylab=yLab,xlim=c(0,1), ylim=c(0,1), las=1)}
      cutP <- keyVal[which(keyVal[,1]==0.05),-1]
      newPch <- cbind(c(1,16,2,17, 7,15,5,6), new=c(21,21,24,24,22,22,23,25))                                   # transform open or plain filled points to color-filled
      if(pch %in% newPch[,1]){ pch2 <- newPch[which(newPch[,1]==pch),2]; bg <- color; col2 <- grDevices::grey(0.2)} else {pch2 <- pch; col2 <- color}   
      col2 <- wrMisc::convColorToTransp(col2, alph=0.9)
      graphics::points(1-cutP["spec"], cutP["sens"], col=col2, pch=pch2, bg=bg, cex=1.4) 
      graphics::mtext(paste("AUC = ",signif(AucROC(keyVal),3)), side=3, cex=0.8,adj=0)
  } }
  keyVal }     
    
