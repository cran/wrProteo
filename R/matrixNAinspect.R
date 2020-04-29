#' Histogram of content of NAs in matrix
#'
#' \code{matrixNAinspect} makes histograms of the full data and shows sub-population of \code{NA}-neighbour values.  
#' The aim of this function is to investigate the nature of \code{NA} values in matrix (of experimental measures) where replicate measurements are available.
#' If a given element was measured twice, and one of these measurements revealed a \code{NA} while the other one gave a (finite) numeric value, the non-NA-value is considered a \code{NA}-neighbour.  
#' The subpopulation of these \code{NA}-neighbour values will then be highlighted in the resulting histogram.
#' In a number of experimental settiongs some actual measurements may not meet an arbitrary defined baseline (as 'zero') or may be too low to be distinguishable from noise that 
#' associated measures were initially recorded as \code{NA}. In several types of measurments in proteomics and transcriptomics this may happen.
#' So this fucntion allows to collect all \code{NA}-neighbour values and compare them to the global distribution of the data to investigate if \code{NA}-neighbours are typically very low values.
#' In case of data with multiple replicates \code{NA}-neighbour values may be distinguished for the case of 2 \code{NA} per group/replicate-set.
#' The resulting plots are typically used to decide if and how \code{NA} values may get replaced by imputed random values or wether measues containing \code{NA}-values should rather me omitted.
#' Of course, such decisions do have a strong impact on further steps of data-analysis and should be performed with care. 
#' 
#' @param dat (matrix or data.frame) main numeric data
#' @param gr (charcter or factor) grouping of columns of dat indicating who is a replicate of whom (ie the length of 'gr' must be equivalent to the number of columns in 'dat')
#' @param retnNA (logical) report number of NAs in graphic
#' @param xLab (character) custom x-label
#' @param tit (character) custom title
#' @param xLim (numerical,length=2) custom x-axis limits
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of messages produced
#' @return graphic only
#' @seealso  \code{\link[graphics]{hist}}, \code{\link[stats]{na.fail}}, \code{\link[wrMisc]{naOmit}} 
#' @examples
#' set.seed(2013)
#' datT6 <- matrix(round(rnorm(300)+3,1),ncol=6,dimnames=list(paste("li",1:50,sep=""),letters[19:24]))
#' datT6 <- datT6 +matrix(rep(1:nrow(datT6),ncol(datT6)),ncol=ncol(datT6))
#' datT6[6:7,c(1,3,6)] <- NA
#' datT6[which(datT6 < 11 & datT6 > 10.5)] <- NA
#' datT6[which(datT6 < 6 & datT6 > 5)] <- NA
#' datT6[which(datT6 < 4.6 & datT6 > 4)] <- NA
#' matrixNAinspect(datT6,gr=gl(2,3)) 
#' @export
matrixNAinspect <- function(dat,gr,retnNA=TRUE,xLab=NULL,tit=NULL,xLim=NULL,silent=FALSE,callFrom=NULL) {
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="matrixNAinspect")
  if(length(dim(dat)) !=2) stop("'dat' must be matrix or data.frame with >1 columns")
  if(is.data.frame(dat)) dat <- as.matrix(dat)
  if(length(gr) != ncol(dat)) stop("Number of columns in 'dat' and number of (group-)elements in 'gr' do not match !")
  if(!is.factor(gr)) gr <- as.factor(gr)
  if(is.null(xLab)) xLab <- "values"            
  quaCol <- if(requireNamespace("RColorBrewer", quietly=TRUE))  RColorBrewer::brewer.pal(4,"Set1")[c(3,2,4)] else c(3:4,2) 
  if(is.null(tit)) tit <- "Distribution of values and NA-neighbours"
  cexMain <- if(nchar(tit) < 25) 1.8 else 1.2
  ## main
  NAneig <- NAneig2 <- numeric() 
  isNA <- is.na(dat)
  chNA <- any(isNA)
  nNAmat <- matrix(0,nrow=nrow(dat),ncol=length(levels(gr)),dimnames=list(NULL,levels(gr)))
  colPanel <- c(grDevices::grey(0.6),grDevices::rgb(0,0.7,0,0.6),grDevices::rgb(0.15,0.15,0.7,0.7))
  if(chNA) {
    for(i in c(1:length(levels(gr)))) {
      curCol <- which(gr==levels(gr)[i]) 
      nNAmat[,i] <- if(length(curCol) >1) rowSums(isNA[,curCol]) else (isNA[,curCol])    
      maxCol <- length(curCol)
      useLi <- which(nNAmat[,i] >0 & nNAmat[,i] < maxCol)           # 1 or 2 NAs  (but not all)
      useL2 <- which(nNAmat[,i] >1 & nNAmat[,i] < maxCol)           # just 2 NAs  (but not all)
      if(length(useLi) >0) NAneig <- c(NAneig,wrMisc::naOmit(as.numeric(dat[useLi,curCol])))
      if(length(useL2) >0) NAneig2 <- c(NAneig2,wrMisc::naOmit(as.numeric(dat[useL2,curCol])))
      }
    n <- c(sum(!is.na(dat)),length(NAneig),length(NAneig2))
    perc <- c("",paste(" (",round(100*n[2:3]/n[1],1),"%)"))
    hi1 <- graphics::hist(dat,breaks="FD",plot=FALSE)
    if(is.null(xLim)) graphics::plot(hi1,border=grDevices::grey(0.85),col=grDevices::grey(0.92),xlab=xLab,las=1,main=tit,cex.main=cexMain) else {
      graphics::plot(hi1,border=grDevices::grey(0.85),col=grDevices::grey(0.92),xlab=xLab,las=1,main=tit,xlim=xLim,cex.main=cexMain)}
    graphics::abline(v=stats::quantile(dat,c(0.03,0.05,0.1),na.rm=TRUE),col=quaCol,lty=2)
    graphics::mtext(paste(c(" (bar) all data",paste(" (box) ",c("any","min 2")," NA-neighbour values"))," n=",n,perc),col=colPanel[1:3],cex=0.65,adj=0,line=c(0.6,-0.1,-0.7),side=3)
    graphics::mtext(paste(" - -",c(3,5,10),"%-quantile"),col=quaCol,cex=0.6,adj=0,line=c(-1.4,-2,-2.6),side=3)
    if(chNA) graphics::hist(NAneig,breaks=hi1$breaks,border=grDevices::grey(0.75),col=grDevices::rgb(0.1,1,0.1,0.15),add=TRUE);                        # in green
    if(chNA) graphics::hist(NAneig2,breaks=hi1$breaks,border=grDevices::grey(0.75),col=grDevices::rgb(0,0,0.7,0.2),add=TRUE);                          # in purple
  } else {graphics::hist(dat,breaks="FD",border=grDevices::grey(0.85),col=grDevices::grey(0.92),xlab=xLab,las=1,main=tit,cex.main=cexMain)
    graphics::mtext(paste(" (bar) all data  n=",length(dat)),col=colPanel[1],cex=0.7,adj=0,line=0.6,side=3)
    graphics::abline(v=stats::quantile(dat,c(0.03,0.05,0.1),na.rm=TRUE),col=quaCol,lty=2) 
    graphics::mtext(paste(" - -",c(3,5,10),"%-quantile"),col=quaCol,cex=0.6,adj=0,line=c(-1.4,-2,-2.6),side=3)}
}
   
