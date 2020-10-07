#' Impute random values to NAs in matrix based on replicates (neighbour) values
#'
#' It is assumed that \code{NA}-values appear the data when quantitation values are very low, as this appears eg in proteomics. 
#' Thus the remaining lowest values may be used to guide imputation.
#' Here, groups of replicate samples (grouping defined via \code{gr} of columns of \code{dat}) are inspected for each line to gather NA-neighbour values.
#' Eg, if a given line contains for a set of 4 replicates 2 \code{NA}-values, the remaining 2 non-\code{NA}-values will be considered as NA-neighbours.
#' Then, this function replaces \code{NA}-values based the sub-population of all NA-neighbours (across all groups of replicates and all lines), assuming a Gaussian distribution.
#' Indeed, in a number of experimental settings some actual measurements may not meet an arbitrary defined baseline (as 'zero') or may be too low to be distinguishable from  
#' noise that associated measures were initially recorded as \code{NA}. In several types of (quantitative) measurments in proteomics and transcriptomics this is known to happen.
#' So this function allows to model and subsequently replace all \code{NA}-values by Gaussian random values based on the characteristics of \code{NA}-neighbours in the same data-set.
#' However, defining these characteristics (via the arguments \code{avSdH} and \code{avSdL}) may be very delicate and visual verification of the plots produced is highly encouraged ! 
#' If more than 300 \code{NA}-neighbours were detected, the imputation will be based on a more restricted sub-set of data with >1 \code{NA} values (ie via the argument \code{avSdH}). 
#' Optionally a histogram may be plotted showing the initial, imputed and final distribution to check if the global hypothesis that \code{NA}-values arose 
#' from very low measurements and to appreciate the impact of the imputed values to the overall final distribution.
#' Of course, all decisions to replace \code{values} do have a strong impact on further steps of data-analysis and should be performed with care. 
#' Please note, that no distinction is made if values seem totally absent (all values of given line and group) as \code{NA} or partially absent (mixture of \code{NA} and real quantitations).
#' Thus, truly absent groups may be over-estimated.
#'  
#' @param dat (matrix or data.frame) main data (may contain \code{NA})
#' @param gr (character or factor) grouping of columns of 'dat', replicate association
#' @param retnNA (logical) decide if NA values should be removed or retained
#' @param avSdH (numerical,length=2) population characteristics 'high' (mean and sd) for >1 \code{NA}-neighbours (per line)
#' @param avSdL (numerical,length=2) population characteristics 'low' (mean and sd) for >0 \code{NA}-neighbours
#' @param plotHist (logical) decide if supplemental figure with histogram shoud be drawn
#' @param xLab (character) label on x-axis on plot
#' @param tit (character) title on plot
#' @param addImputDetail (logical) display details about data (number of NAs) and imputation in graph (min number of NA-neighbours per protein and group, quantile to model, mean and sd of imputed)
#' @param seedNo (integer) seed-value for normal random values
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @return list with \code{$data} .. matrix of data where \code{NA} are replaced by imputed values, \code{$nNA} .. number of \code{NA} by group, \code{$randParam} .. parameters used for making random data 
#' @seealso  \code{\link[graphics]{hist}}, \code{\link[stats]{na.fail}},  \code{\link[wrMisc]{naOmit}} 
#' @examples
#' set.seed(2013)
#' datT6 <- matrix(round(rnorm(300)+3,1),ncol=6,dimnames=list(paste("li",1:50,sep=""),
#'   letters[19:24]))
#' datT6 <- datT6 +matrix(rep(1:nrow(datT6),ncol(datT6)),ncol=ncol(datT6))
#' datT6[6:7,c(1,3,6)] <- NA
#' datT6[which(datT6 < 11 & datT6 > 10.5)] <- NA
#' datT6[which(datT6 < 6 & datT6 > 5)] <- NA
#' datT6[which(datT6 < 4.6 & datT6 > 4)] <- NA
#' datT6b <- matrixNAneighbourImpute(datT6,gr=gl(2,3))
#' head(datT6b$data)
#' @export
matrixNAneighbourImpute <- function(dat,gr,retnNA=TRUE,avSdH=c(0.18,0.5),avSdL=c(0.1,0.5),plotHist=TRUE,xLab=NULL,tit=NULL,addImputDetail=TRUE,seedNo=2018,silent=FALSE,callFrom=NULL){
  ## replace NA values based on group neigbours (based on grouping of columns in gr), overall assumption of close to Gaussian distrib
  ## return matrix including imputed values or list of final & matrix with number of imputed by group
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="matrixNAneighbourImpute")
  if(length(dim(dat)) !=2) stop("'dat' must be matrix or data.frame with >1 columns")
  if(is.data.frame(dat)) dat <- as.matrix(dat)
  if(length(gr) != ncol(dat)) stop("Number of columns in 'dat' and number of (group-)elements in 'gr' do not match !")
  if(!is.factor(gr)) gr <- as.factor(gr)
  if(is.null(xLab)) xLab <- "values"            
  ## main
  isNA <- is.na(dat)
  chNA <- any(isNA)
  nNAmat <- matrix(0, nrow=nrow(dat), ncol=length(levels(gr)), dimnames=list(NULL,unique(wrMisc::naOmit(gr))))
  if(!chNA) {
    if(plotHist) {graphics::hist(dat, br="FD", border=grDevices::grey(0.85), col=grDevices::grey(0.92), xlab=xLab, las=1, main=tit)
      graphics::mtext("no NA-replacement needed  ",adj=1,cex=0.6,line=-0.3)    
      graphics::mtext(paste("  n=",length(dat)),side=3,line=-0.3,cex=0.55,adj=0,col=grDevices::grey(0.3))}
    return(if(retnNA) list(data=dat,nNA=nNAmat) else dat)
  } else {
    NAneig <- NAneig2 <- numeric() 
    grLev <- unique(wrMisc::naOmit(gr))
    for(i in c(1:length(grLev))) {
      curCol <- which(gr==grLev[i]) 
      nNAmat[,i] <- rowSums(isNA[,curCol])     
      maxCol <- length(curCol)
      useLi <- which(nNAmat[,i] >0 & nNAmat[,i] < maxCol-1)           # 1 or more NA
      useL2 <- which(nNAmat[,i] >1 & nNAmat[,i] < maxCol-1)           # min 2 NAs
      if(length(useLi) >0) NAneig <- c(NAneig, wrMisc::naOmit(as.numeric(dat[useLi,curCol])))
      if(length(useL2) >0) NAneig2 <- c(NAneig2, wrMisc::naOmit(as.numeric(dat[useL2,curCol])))
    }
    ## need to optimize : the higher NA% the higher aver to model  
    randParam <- if(length(NAneig2) >300) c(stats::quantile(NAneig2,avSdH[1]), stats::sd(NAneig2)*avSdH[2]) else c(stats::quantile(NAneig,avSdL[1]), stats::sd(NAneig)*avSdL[2])
    randParam <- c(mean=signif(randParam[1],4), sd=signif(randParam[2],4), impQuant=NA,seed=seedNo)
    randParam[3] <- signif(which.min(abs(dat -randParam[1]))/sum(!isNA),4)       # recalculate quantile
    dat1 <- .imputeNA(dat, gr=gr, impParam=randParam, exclNeg=TRUE)
    lowValMod <- dat1$lowValMod
    dat1 <- dat1$data
    ranPar2 <- if(length(NAneig2) >300) c(n=length(NAneig2), avSdH, 2) else c(n=length(NAneig), avSdL, 1)
    msg <- list(c(" n.woNA=",sum(!isNA)," n.NA =",sum(isNA)), 
      c("model",100*ranPar2[2],"%-tile of (min",ranPar2[4],"NA/grp)",ranPar2[1],"NA-neighbour values"),
      c("imputation: mean=",signif(randParam[1],3),"  sd=",signif(randParam[2],2)))
    if(!silent) message(fxNa,paste(sapply(msg,paste,collapse=" "),collapse="\n    "))  
    if(plotHist) {
      ranPar2 <- if(length(NAneig2) >300) c(n=length(NAneig2), avSdH,2) else c(n=length(NAneig), avSdL, 1)
      hi1 <- graphics::hist(dat1,br="FD",col=grDevices::grey(0.9),border=grDevices::grey(0.8),xlab=xLab,las=1,main=paste(tit,"at NA-replacement"))  # grey cols (final distr)
      if(FALSE) graphics::abline(v=stats::quantile(dat,0.03,na.rm=TRUE),lty=2,col="red")                                                     # 3% quantile as red line
      graphics::hist(dat,breaks=hi1$breaks,border=grDevices::grey(0.75),col=grDevices::rgb(0.1,1,0.1,0.15),add=TRUE)                          # orig data in green
      graphics::hist(lowValMod[1:sum(isNA)],br=hi1$breaks,border=grDevices::grey(0.75),col=grDevices::rgb(0,0,0.7,0.2),add=TRUE)           # add purple hist to plot
      colPanel <- c(grDevices::grey(0.6),grDevices::rgb(0,0.7,0,0.6),grDevices::rgb(0.15,0.15,0.7,0.7))
      if(addImputDetail) graphics::mtext(paste(sapply(msg,paste,collapse=" "),collapse="\n "),side=3,line=-1,cex=0.75,adj=0,col=grDevices::grey(0.3))
      graphics::legend("topright",c("final","initial","imputed"),col=colPanel,text.col=colPanel,cex=0.9,seg.len=0.3,lwd=4)}
    if(retnNA) list(data=dat1,nNA=nNAmat,randParam=randParam) else dat1 }}

#' @export  
.imputeNA <- function(dat,gr,impParam,exclNeg=TRUE,inclLowValMod=TRUE) {
  ## basic NA imputation for 'dat' using 'impParam' (mean, sd and seed)
  ## 'impParam' .. (numeric) 1st for mean; 2nd for sd; 3rd for seed
  isNa <- is.na(dat)
  if(length(impParam) >3) set.seed(as.integer(impParam[4]))
  lowValMod <- stats::rnorm(round(1.5*sum(isNa)), mean=impParam[1], sd=impParam[2])
  if(exclNeg) lowValMod <- lowValMod[which(lowValMod >0)]
  dat[which(isNa)] <- lowValMod[1:sum(isNa)]
  if(inclLowValMod) list(data=dat, lowValMod=lowValMod) else dat }
   
