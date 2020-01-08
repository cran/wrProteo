#' Plot ROC curves
#'
#' \code{plotROC} plots ROC curves based on results from \code{\link{summarizeForROC}}. 
#' Does not return any data, plot only.
#'  
#' @param dat (matrix) from testing (eg  \code{\link{summarizeForROC}} )
#' @param ...  additional input (must be of same type of format as 'dat')
#' @param useCol (integer or character) colors to be used
#' @param methNames (character) names of methods (data-sets) to be displayed
#' @param col (character) custom color
#' @param pch (integer) type of symbol to be used (see \code{\link[graphics]{par}})
#' @param bg (character) background color in plot (see \code{\link[graphics]{par}})
#' @param tit (character) custom title
#' @param point05 (numeric) specific point to highlight in plot (typically at alpha=0.05) 
#' @param pointSi (numeric) size of points (as expansion factor \code{cex}) 
#' @param nByMeth (integer) value of n to display
#' @param colPanel (character) custom colors
#' @param speciesOrder (integer) custom order of species in legend
#  @param speciesOrder (integer) optional custom order for counts per species (eg number of proteins) in legend (eg 'n.H/S/E')
#' @param txtLoc (numeric) location for text
#' @param legCex (numeric) cex expansion factor for legend
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of messages produced
#' @return plot only
#' @seealso  \code{\link{summarizeForROC}}, \code{\link[wrMisc]{moderTest2grp}}  
#' @examples
#' set.seed(2019); test1 <- list(annot=cbind(spec=c(rep("b",35),letters[sample.int(n=3,
#'   size=150,replace=TRUE)])),BH=matrix(c(runif(35,0,0.01),runif(150)),ncol=1))
#' tail(roc1 <- summarizeForROC(test1,spec=c("a","b","c"),plotROC=FALSE))
#' plotROC(roc1)
#' @export
plotROC <- function(dat,...,useCol=2:3,methNames=NULL,col=NULL,pch=1,bg=NULL,tit=NULL,point05=0.05,pointSi=0.85,nByMeth=NULL,
  colPanel=c(grDevices::grey(0.4),2:5),speciesOrder=NULL,txtLoc=c(0.4,0.3,0.04),legCex=0.72,silent=FALSE,callFrom=NULL) {
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="plotROC")
  inpS <- list(...)
  addSupl <- TRUE
  chInp <- lapply(c("dat","useCol","methNames","col","pch","bg","tit","point05","pointSi","nByMeth","txtLoc","legCex"),wrMisc::.seqCutStr,startFr=2,reverse=TRUE)
  chAr <- names(inpS) %in% unlist(chInp)
  if(any(chAr)) {
    ## if argument names changed/ not complete need to change/adjust code here !!
    inpS <- inpS[which(!chAr)]  
  }
  if(is.null(tit)) tit <- paste("ROC")
  if(is.null(col)) col <- 1:(1+length(inpS))
  xLab <- "1 - Specificity"
  yLab <- "Sensitivity"
  graphics::plot(1-dat[,useCol[1]],dat[,useCol[2]],type="n",col=col[1],pch=pch,bg=bg,main=tit,xlab=xLab,ylab=yLab,xlim=c(0,1),ylim=c(0,1))
  col2 <- col
  cutP <- dat[which(dat[,1]==point05),]
  if(length(stats::na.omit(point05))==1) { 
    newPch <- cbind(c(1,16,2,17, 7,15,5,6),new=c(21,21,24,24,22,22,23,25))                                   # transform open or plain filled points to color-filled
    if(pch %in% newPch[,1]){ pch2 <- newPch[which(newPch[,1]==pch),2]; bg <- wrMisc::.convColorToTransp(col,0.1); col2 <- grDevices::grey(0.4)} else {pch2 <- pch; col2 <- col}   
    graphics::points(1-cutP[useCol[1]],cutP[useCol[2]],col=col2,pch=pch2,bg=col[1],cex=pointSi) }
  graphics::points(1-dat[,useCol[1]],dat[,useCol[2]],type="s",col=col[1],pch=pch,bg=bg)    # main curve
  coColN <- colnames(dat)[wrMisc::naOmit(grep("n\\.pos\\.",colnames(dat)))]            # more flexible (also to number of species/tags)
  if(length(coColN) <2) { coColN <- colnames(dat)[(ncol(dat)-2):ncol(dat)]
    if(!silent) message(fxNa," Can't find 'n.pos.' tag among colnames of 'dat', assuming last 3 columns") }
  if(length(speciesOrder) <length(coColN)) speciesOrder <- c(1:length(coColN))
  coColN <- coColN[speciesOrder]
  coColN1 <- sub("n\\.pos\\.","",coColN)
  coColN2 <- paste(" n.",paste(coColN1,sep="",collapse="/")," ",sep="")
  if(addSupl) {       # add legend-like method-name/descr
    txt <- if(is.null(nByMeth)) methNames[1] else paste(methNames[1]," (n.test=",nByMeth[1],") ",sep="")
    graphics::text(txtLoc[1]-txtLoc[3],txtLoc[2]+txtLoc[3],paste("Values at threshold of",point05,":"),cex=0.75,col=grDevices::grey(0.4),adj=0)   
    graphics::text(txtLoc[1],txtLoc[2],txt,cex=legCex+0.02,col=col[1],adj=1)    
    graphics::text(txtLoc[1]+0.02,txtLoc[2],paste(paste(paste(c("prec=","accur=","FDR="),        #"n.E/S/H="
      signif(cutP[c("prec","accur","FDR")],2)),collapse="  "),coColN2 ,cutP[coColN[1]],cutP[coColN[2]],
      if(length(coColN)>2) cutP[coColN[3]]),cex=legCex,col=col[1],adj=0) }  
  if(length(inpS) >0) {
    for(i in 1:length(inpS)) {
      cutP <- inpS[[i]][which(inpS[[i]][,1]==point05),]
      if(point05) graphics::points(1-cutP[useCol[1]],cutP[useCol[2]],col=col2,pch=pch2[1],bg=col[i+1],cex=pointSi)
      graphics::points(1-inpS[[i]][,useCol[1]],inpS[[i]][,useCol[2]],type="s",col=col[i+1],pch=pch,bg=bg[i+1])
      if(TRUE) {
        txt <- if(is.null(nByMeth)) methNames[i+1] else paste(methNames[i+1]," (n.test=",nByMeth[i+1],") ",sep="")
        graphics::text(txtLoc[1],txtLoc[2]-txtLoc[3]*i,txt,cex=legCex+0.02,col=col[i+1],adj=1)  
        graphics::text(txtLoc[1]+0.02,txtLoc[2]-txtLoc[3]*i,paste(paste(paste(c("prec=","accur=","FDR="),
          signif(cutP[c("prec","accur","FDR")],2)),collapse="  "),coColN2,cutP[coColN[1]],cutP[coColN[2]],
          if(length(coColN)>2) cutP[coColN[3]]),cex=legCex,col=col[i+1],adj=0)
        }   
    } }
}
     
