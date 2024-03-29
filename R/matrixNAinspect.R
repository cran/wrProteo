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
#' @param debug (logical) additional messages for debugging 
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function produces a graphic (to the current graphical device)
#' @seealso  \code{\link[graphics]{hist}}, \code{\link[stats]{na.fail}}, \code{\link[wrMisc]{naOmit}} 
#' @examples
#' set.seed(2013)
#' datT6 <- matrix(round(rnorm(300)+3,1), ncol=6, 
#'   dimnames=list(paste("li",1:50,sep=""), letters[19:24]))
#' datT6 <- datT6 +matrix(rep(1:nrow(datT6),ncol(datT6)), ncol=ncol(datT6))
#' datT6[6:7,c(1,3,6)] <- NA
#' datT6[which(datT6 < 11 & datT6 > 10.5)] <- NA
#' datT6[which(datT6 < 6 & datT6 > 5)] <- NA
#' datT6[which(datT6 < 4.6 & datT6 > 4)] <- NA
#' matrixNAinspect(datT6, gr=gl(2,3)) 
#' @export
matrixNAinspect <- function(dat, gr=NULL, retnNA=TRUE, xLab=NULL, tit=NULL, xLim=NULL, silent=FALSE, debug=FALSE, callFrom=NULL) {
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="matrixNAinspect")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  ## extract if object
  if(is.list(dat)) {
    if("sampleSetup" %in% names(dat) & length(gr) <1) {
      gr <- dat$sampleSetup$lev
      if(length(gr) >1 && length(dat$sampleSetup$col) <2) names(gr) <- dat$sampleSetup$meta[,dat$sampleSetup$col]  # in case names are not provided
    }
    if(!silent) message(fxNa,"Trying to extract quantitation data to use as 'dat' out of list ..")
    dat <- dat$quant }
  
  if(any(length(dim(dat)) !=2, dim(dat) < 2, na.rm=TRUE)) stop("Invalid argument 'dat'; must be matrix or data.frame with min 2 lines and 2 cols")
  if(is.data.frame(dat)) dat <- as.matrix(dat)
  chGr <- FALSE
  if(length(gr) != ncol(dat)) stop("Number of columns in 'dat' and number of (group-)elements in 'gr' do not match !")
  if(length(gr)==length(unique(gr))) { hasNaNeigh <- FALSE
    if(!silent) message(fxNa,"NOTE : The argument 'gr' does not designate any replicates, can't determine NA-neighbours !")
  } else hasNaNeigh <- TRUE
  
  if(!is.factor(gr)) gr <- as.factor(gr)
  if(is.null(xLab)) xLab <- "(log2) Abundance"
  chRColB <- requireNamespace("RColorBrewer", quietly=TRUE) 
  if(!chRColB) message(fxNa,"More/better colors may be displayed with package 'RColorBrewer' installed; consider installing it !")
  quaCol <- if(chRColB) RColorBrewer::brewer.pal(4,"Set1")[c(3,2,4)] else c(3:4,2) 
  if(is.null(tit)) tit <- "Distribution of values and NA-neighbours"
  cexMain <- if(nchar(tit) < 25) 1.4 else 1.1
  ## main
  NAneig <- NAneig2 <- numeric() 
  isNA <- is.na(dat)
  chNA <- any(isNA)
  nNAmat <- matrix(0, nrow=nrow(dat), ncol=length(levels(gr)), dimnames=list(NULL,levels(gr)))
  colPanel <- c(grDevices::grey(0.6), grDevices::rgb(0,0.7,0,0.6), grDevices::rgb(0.15,0.15,0.7,0.7))
  if(debug) {message(fxNa,"Ffound ",sum(isNA,na.rm=TRUE)," NAs (out of ",prod(dim(dat))," values); chNA=",chNA,"   mMNi0"); 
    mMNi0 <- list(dat=dat,gr=gr,isNA=isNA,chNA=chNA,nNAmat=nNAmat,hasNaNeigh=hasNaNeigh)}
  if(chNA & hasNaNeigh) {
    ## extract NA-neighbours
    for(i in c(1:length(levels(gr)))) {
      curCol <- which(gr==levels(gr)[i]) 
      nNAmat[,i] <- if(length(curCol) >1) rowSums(isNA[,curCol]) else (isNA[,curCol])    
      maxCol <- length(curCol)
      useLi <- which(nNAmat[,i] >0 & nNAmat[,i] < maxCol)           # 1 or 2 NAs  (but not all)
      useL2 <- which(nNAmat[,i] >1 & nNAmat[,i] < maxCol)           # just 2 NAs  (but not all)
      if(length(useLi) >0) NAneig <- c(NAneig, wrMisc::naOmit(as.numeric(dat[useLi,curCol])))
      if(length(useL2) >0) NAneig2 <- c(NAneig2, wrMisc::naOmit(as.numeric(dat[useL2,curCol])))
      }
    n <- c(sum(!is.na(dat)), length(NAneig), length(NAneig2))
    perc <- c("",paste(" (",round(100*n[2:3]/n[1],1),"%)"))
    if(debug) {message(fxNa,"mMNi1  n=",wrMisc::pasteC(n))}
    hi1 <- graphics::hist(dat, breaks="FD", plot=FALSE)
    if(is.null(xLim)) graphics::plot(hi1, border=grDevices::grey(0.85), col=grDevices::grey(0.92), xlab=xLab, las=1, main=tit, cex.main=cexMain) else {
      graphics::plot(hi1, border=grDevices::grey(0.85), col=grDevices::grey(0.92), xlab=xLab, las=1, main=tit, xlim=xLim,cex.main=cexMain)}
    graphics::abline(v=stats::quantile(dat,c(0.05,0.1,0.15),na.rm=TRUE), col=c(quaCol[-1],"tomato3"), lty=2)
    graphics::mtext(paste(c(" (bar) all data",paste(" (box) ",c("any","min 2")," NA-neighbour values"))," n=",n,perc), col=colPanel[1:3],cex=0.65,adj=0,line=c(0.6,-0.1,-0.7),side=3)
    graphics::mtext(paste(" - -",c(5,10,15),"%-quantile (all data)"), col=c(quaCol[-1],"tomato3"), cex=0.6, adj=0, line=c(-1.5,-2.1,-2.7), side=3)
    if(length(NAneig) >10) {    # display mode
      yLim <- signif(graphics::par("usr")[3:4], 3)             # current y-limits
      mod <- signif(wrMisc::stableMode(if(length(NAneig2) >300) NAneig2 else NAneig, method="density"), 3)
      graphics::mtext(paste(" (arrow) mode of",if(length(NAneig2) >300) "2-"," NA-neighbours :",signif(mod,3)), col="sienna2", cex=0.7, adj=0, line=-3.4, side=3)
      graphics::arrows(mod, yLim[1]+(yLim[2]-yLim[1])*0.4, mod, yLim[1]+(yLim[2]-yLim[1])/4, length=0.1, col="sienna2", lwd=2)
    }    
    graphics::hist(NAneig, breaks=hi1$breaks, border=grDevices::grey(0.75), col=grDevices::rgb(0.1,1,0.1,0.15), add=TRUE);                    # in green
    graphics::hist(NAneig2, breaks=hi1$breaks, border=grDevices::grey(0.75), col=grDevices::rgb(0,0,0.7,0.2), add=TRUE);                      # in purple
    
  } else {
    graphics::hist(dat, breaks="FD", border=grDevices::grey(0.85), col=grDevices::grey(0.92), xlab=xLab, las=1,main=tit,cex.main=cexMain)
    graphics::mtext(paste(" (bar) all data  n=",length(dat)), col=colPanel[1], cex=0.7, adj=0, line=0.6, side=3)
    graphics::abline(v=stats::quantile(dat,c(0.05,0.1,0.15),na.rm=TRUE), col=c(quaCol[-1],"tomato3"), lty=2) 
    graphics::mtext(paste(" - -",c(5,10,15),"%-quantile (all data)"), col=c(quaCol[-1],"tomato3"), cex=0.6, adj=0, line=c(-1.5,-2.1,-2.7), side=3)}
}
    
