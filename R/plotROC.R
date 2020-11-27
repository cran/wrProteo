#' Plot ROC curves
#'
#' \code{plotROC} plots ROC curves based on results from \code{\link{summarizeForROC}}. 
#' Does not return any data, plot only. Allows printing simultaneously multiple ROC curves from different studies. 
#' Was made for special consideration of 3 species mix as in proteomics benchmark

#' In the simplest case data were prepared using \code{\link[wrMisc]{moderTest2grp}} 
#'  
#' @param dat (matrix) from testing (eg  \code{\link{summarizeForROC}} )
#' @param ... optional additional data-sets to include as seprate ROC-curves to same plot (must be of same type of format as 'dat')
#' @param useCol (integer or character, length=2) columns from \code{dat} to be used for pecificity and sensitivity
#' @param methNames (character) names of methods (data-sets) to be displayed
#' @param col (character) custom colors for lines and text (choose one color for each different data-set)
#' @param pch (integer) type of symbol to be used (see also \code{\link[graphics]{par}})
#' @param bg (character) background color in plot (see also \code{\link[graphics]{par}})
#' @param tit (character) custom title
#' @param xlim (numeric, length=2) custom x-axis limits
#' @param ylim (numeric, length=2) custom y-axis limits
#' @param point05 (numeric) specific point to highlight in plot (typically at alpha=0.05) 
#' @param pointSi (numeric) size of points (as expansion factor \code{cex}) 
#' @param nByMeth (integer) value of n to display
#' @param speciesOrder (integer) custom order of species in legend
#  @param speciesOrder (integer) optional custom order for counts per species (eg number of proteins) in legend (eg 'n.H/S/E')
#' @param txtLoc (numeric, length=3) location for text (x, y and proportional factor for line-offset, default is c(0.4,0.3,0.04))
#' @param legCex (numeric) cex expansion factor for legend (see also \code{\link[graphics]{par}})
#' @param las (numeric) factor for text-orientation (see also \code{\link[graphics]{par}})
#' @param addSuplT (logical) add text with information about precision,accuracy and FDR
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of messages produced
#' @return plot only
#' @seealso \code{\link[wrProteo]{summarizeForROC}}, \code{\link[wrMisc]{moderTest2grp}}  
#' @examples
#' roc0 <- cbind(alph=c(2e-6,4e-5,4e-4,2.7e-3,1.6e-2,4.2e-2,8.3e-2,1.7e-1,2.7e-1,4.1e-1,5.3e-1,
#' 	 6.8e-1,8.3e-1,9.7e-1), spec=c(1,1,1,1,0.957,0.915,0.915,0.809,0.702,0.489,0.362,0.234,
#'   0.128,0.0426), sens=c(0,0,0.145,0.942,2.54,2.68,3.33,3.99,4.71,5.87,6.67,8.04,8.77,
#'   9.93)/10, n.pos.a=c(0,0,0,0,2,4,4,9,14,24,36,41) )
#' plotROC(roc0)
#' @export
plotROC <- function(dat,...,useCol=2:3, methNames=NULL, col=NULL, pch=1, bg=NULL, tit=NULL, xlim=NULL, ylim=NULL, point05=0.05, pointSi=0.85, nByMeth=NULL,
  speciesOrder=NULL, txtLoc=NULL, legCex=0.72, las=1, addSuplT=TRUE, silent=FALSE, callFrom=NULL) {
  ##
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="plotROC")
  inpS <- list(...)
  chInp <- lapply(c("dat","useCol","methNames","col","pch","bg","tit","point05","pointSi","nByMeth","txtLoc","legCex"),wrMisc::.seqCutStr,startFr=2,reverse=TRUE)
  chAr <- names(inpS) %in% unlist(chInp)
  if(any(chAr)) {
    ## if argument names changed/ not complete need to change/adjust code here !!
    inpS <- inpS[which(!chAr)]  
  }
  if(is.null(tit)) tit <- paste("ROC")
  if(is.null(col)) col <- if(length(inpS)==1) 1 else c(grDevices::grey(0.4), 2:(1+length(inpS)))
  xLab <- "1 - Specificity"
  yLab <- "Sensitivity"
  if(!is.numeric(xlim) | length(xlim) !=2) xlim <- c(0,1)
  graphics::plot(1 -dat[,useCol[1]], dat[,useCol[2]], type="n", col=col[1], pch=pch, bg=bg, main=tit, xlab=xLab, ylab=yLab, xlim=xlim, ylim=if(length(ylim)==2) ylim else c(0,1),las=las)
  col2 <- col
  cutP <- dat[which(dat[,1]==point05),]
  if(length(stats::na.omit(point05))==1) { 
    newPch <- cbind(c(1,16,2,17, 7,15,5,6), new=c(21,21,24,24,22,22,23,25))                                   # transform open or plain filled points to color-filled
    if(pch %in% newPch[,1]) { pch2 <- newPch[which(newPch[,1]==pch),2]; bg <- wrMisc::convColorToTransp(col,0.1);
      col2 <- grDevices::grey(0.4) } else {pch2 <- pch; col2 <- col}   
    graphics::points(1 -cutP[useCol[1]], cutP[useCol[2]], col=col2, pch=pch2, bg=col[1], cex=pointSi) }
  graphics::points(1 -dat[,useCol[1]], dat[,useCol[2]], type="s", col=col[1], pch=pch, bg=bg)    # main curve
  coColN <- colnames(dat)[wrMisc::naOmit(grep("n\\.pos\\.",colnames(dat)))]            # more flexible (also to number of species/tags)
  if(length(coColN) <2) { coColN <- colnames(dat)[(ncol(dat) -2):ncol(dat)]
    if(!silent) message(fxNa," Can't find 'n.pos.' tag among colnames of 'dat', assuming last 3 columns") }
  if(length(speciesOrder) <length(coColN)) speciesOrder <- c(1:length(coColN))
  coColN <- coColN[speciesOrder]
  coColN1 <- sub("n\\.pos\\.","",coColN)
  coColN2 <- paste(" n.",paste(coColN1,sep="",collapse="/")," ",sep="")
  if(length(txtLoc) !=3) txtLoc <- graphics::par("usr")
  if(length(txtLoc) !=3) { figDim <- signif(graphics::par("usr"),3)
    txtLoc <- c(x=figDim[1] +0.42*(figDim[2] -figDim[1]), y=figDim[3] + (0.3 +length(inpS))*(figDim[4] -figDim[3])/30, fac=0.035*(figDim[4] -figDim[3])) }
   
  if(addSuplT) {            # add legend-like method-name/descr
    txt <- if(is.null(nByMeth)) methNames[1] else paste(methNames[1]," (n.test=",nByMeth[1],") ",sep="")
    graphics::text(txtLoc[1] -txtLoc[3], txtLoc[2] +txtLoc[3], paste("Values at threshold of",point05,":"), cex=0.75, col=grDevices::grey(0.4), adj=0)   
    graphics::text(txtLoc[1], txtLoc[2], txt, cex=legCex+0.02, col=col[1], adj=1)    
    if(addSuplT) graphics::text(txtLoc[1] +0.02, txtLoc[2], paste(paste(paste(c("prec=","accur=","FDR="),        #"n.E/S/H="
      signif(cutP[c("prec","accur","FDR")],2)),collapse="  "),coColN2 ,cutP[coColN[1]],cutP[coColN[2]],
      if(length(coColN)>2) cutP[coColN[3]]), cex=legCex, col=col[1], adj=0) }  
  if(length(inpS) >0) {
    for(i in 1:length(inpS)) {
      if(length(inpS[[i]]) >0) if(nrow(inpS[[i]]) >0) {
        cutP <- inpS[[i]][which(inpS[[i]][,1]==point05),]
        if(point05) graphics::points(1 -cutP[useCol[1]], cutP[useCol[2]], col=col2, pch=pch2[1], bg=col[i+1], cex=pointSi)      # new point at alpha
        graphics::points(1-inpS[[i]][,useCol[1]],inpS[[i]][,useCol[2]], type="s",col=col[i+1], pch=pch, bg=bg[i+1])             # new ROC curve
        if(addSuplT) {
          txt <- if(is.null(nByMeth)) methNames[i+1] else paste(methNames[i+1]," (n.test=",nByMeth[i+1],") ",sep="")
          graphics::text(txtLoc[1], txtLoc[2]-txtLoc[3]*i,txt, cex=legCex+0.02, col=col[i+1], adj=1)                     # first block
          graphics::text(txtLoc[1]+0.02, txtLoc[2]-txtLoc[3]*i, paste(paste(paste(c("prec=","accur=","FDR="),
            signif(cutP[c("prec","accur","FDR")],2)),collapse="  "),coColN2,cutP[coColN[1]],cutP[coColN[2]],            # first block (with counting data)
            if(length(coColN)>2) cutP[coColN[3]]), cex=legCex, col=col[i+1], adj=0)
        }   
    } } }
}
     
