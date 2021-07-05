#' Add arrow for expected Fold-Change to VolcanoPlot or MA-plot  
#'
#' NOTE : This function is deprecated, please use \code{\link[wrGraph]{foldChangeArrow}} instead !!
#' This function was made for adding an arrow indicating a fold-change to MA- or Volcano-plots. 
#' When comparing mutiple concentratios of standards in benchmark-tests it may be useful to indicate the expected ratio in a pair-wise comparison.
#' In case of main input as list or MArrayLM-object (as generated from limma), the colum-names of multiple pairwise comparisons can be used 
#' for extracting a numeric content (supposed as concentrations in sample-names) which will be used to determine the expected ratio used for plotting. 
#' Optionally the ratio used for plotting can be returned as numeric value. 
#' 
#' @param FC (numeric, list or MArrayLM-object) main information for drawing arrow : either numeric value for fold-change/log2-ratio of object to search for colnames of statistical testing for extracting numeric part
#' @param useComp (integer) only used in case FC is list or MArrayLM-object an has multiple pairwise-comparisons  
#' @param isLin (logical) inidicate if \code{FC} is log2 or not
#' @param asX (logical) indicate if arrow should be on x-axis 
#' @param col (integer or character) custom color 
#' @param arr (numeric, length=2) start- and end-points of arrow (as relative to entire plot)
#' @param lwd (numeric) line-width of arrow
#' @param addText (logical or named vector) indicate if text explaining arrow should be displayed, use \code{TRUE} for default (on top right of plot), 
#'    or any combination of 'loc','line','cex','side','adj','col','text' (or 'txt') for customizing specific elements
#' @param returnRatio (logical) return ratio
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @param silent (logical) suppress messages
#' @return plots arrow only (and explicative text), if \code{returnRatio=TRUE} also returns numeric value for extracted ratio
#'
#' @details The argument \code{addText} also allows specifying a fixed position when using \code{addText=c(loc="bottomleft")}, also bottomright, topleft, topright, toleft and toright may be used.
#'   In this case the elemts \code{side} and \code{adjust} will be redefined to accomodate the text in the corner specified. 
#'
#'  Ultimately this function will be integated to the package wrGraph. 
#'
#' @seealso new version : \code{\link[wrGraph]{foldChangeArrow}}; used with \code{\link[wrGraph]{MAplotW}}, \code{\link[wrGraph]{VolcanoPlotW}}
#' @examples
#' plot(rnorm(20,1.5,0.1),1:20)
#' #deprecated# foldChangeArrow2(FC=1.5) 
#' 
#' @export
foldChangeArrow2 <- function(FC, useComp=1, isLin=TRUE, asX=TRUE, col=2, arr=c(0.005,0.15), lwd=NULL, 
  addText=c(line=-0.9,cex=0.7,txt="expected",loc="toright"), returnRatio=FALSE,silent=FALSE, callFrom=NULL){
  ##
  .Deprecated("Please use the foldChangeArrow() function form the package wrGrpah instead !")
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="foldChangeArrow2")
  figCo <- graphics::par("usr")                         #  c(x1, x2, y1, y2)
  if(all(length(FC) >1, any(c("MArrayLM","list") %in% class(FC)) )) {
    ## try working based on MArrayLM-object or list
    ##  look for names of pairwise comparisons to extract numeric parts for calculating expected ratio 
    chNa <- names(FC) %in% c("t","BH","FDR","p.value")
    if(any(chNa)) {
      if(all(length(useComp)==1, length(dim(FC[[which(chNa)[1]]])) ==2, dim(FC[[which(chNa)[1]]]) > 0:1)) {
        colNa <- colnames(FC[[which(chNa)[1]]])
        ch2 <- colNa[1]=="(Intercept)" & length(colNa)==2        
      } else ch2 <- TRUE 
      if(!ch2) {
        regStr <-"[[:space:]]*[[:alpha:]]+[[:punct:]]*[[:alpha:]]*"
        colNa <- sub(paste0("^",regStr),"", sub(paste0(regStr,"$"), "", unlist(strsplit(colNa[useComp], "-"))) )
        chN2 <- try(as.numeric(colNa), silent=TRUE)
        if(!"try-error" %in% class(chN2) & length(colNa)==2) {
          FC <- chN2[2] / chN2[1]
        } else ch2 <- TRUE
        ## note: wrMisc::numPairDeColNames() sorts numeric values, can't use
        chN2 <- all(length(colNa)==2, nchar(sub("[[:digit:]]*\\.?[[:digit:]]*","",colNa)) <1)   # contains only usable digits
        isLin <- FALSE                # assume log2 when from testing result
      } else ch2 <- TRUE      
    } else ch2 <- TRUE
    if(ch2) FC <- NULL
  }  
  FC <- try(as.numeric(FC), silent=TRUE)
  if(!"try-error" %in% class(FC)) {
    if(!isLin) FC <- log2(FC)
    if(any(identical(addText, TRUE), c("line","cex","side","adj","col","text","txt","loc") %in% names(addText))) {
       cat(" .. FC",FC,"  addText:",addText,"\n")
      ## bottomleft, bottomright, topleft, topright, toleft and toright 
      mLi <- if("line" %in% names(addText)) try(as.numeric(addText["line"][1]),silent=TRUE) else -0.9
      mCex <- if("cex" %in% names(addText)) try(as.numeric(addText["cex"][1]),silent=TRUE) else 0.7
      mSide <- if("side" %in% names(addText)) try(as.integer(addText["side"][1]),silent=TRUE) else 1
      mAdj <- if("adj" %in% names(addText)) try(as.integer(addText["adj"][1]),silent=TRUE) else 1
      mCol <- if("col" %in% names(addText)) addText["col"][1] else col
      mTxt <- if("text" %in% names(addText)) addText["text"][1] else {if("txt" %in% names(addText)) addText["txt"][1] else "arrow: expected="} 
      if("loc" %in% names(addText)) {
        mTxt <- paste(mTxt, signif(if(isLin) 2^FC else FC,3))
        ## check for left/right/center
        chRi <- grep("right$",as.character(addText["loc"]))
        chLe <- grep("left$",as.character(addText["loc"]))
        chCe <- grep("center$",as.character(addText["loc"]))
        if(length(chLe) >0) { mAdj <- 0; mTxt <- paste0(" ",mTxt)                          # this is left.xxx
          if(arr[1] < 0.15 & FC < figCo[1] +diff(figCo[1:2])/3) arr[1] <- 0.015            # raise min starting hight to avoid crossing text
        } else { if(length(chRi) >0) {mAdj <- 1; mTxt <- paste0(mTxt," ")                  # this is right.xxx
          if(arr[1] < 0.15 & FC > figCo[2] -diff(figCo[1:2])/3) arr[1] <- 0.015            # raise min starting hight to avoid crossing text
          } else {
            if(length(chCe) >0) mAdj <- 0.5; if(arr[1] < 0.15) arr[1] <- 0.015 }} 
        ## check for top/bottom
        chTop <- grep("^top",as.character(addText["loc"]))
        chBot <- grep("^bottom",as.character(addText["loc"]))
        if(length(chTop) >0) mSide <- 3 else if(length(chBot) >0) mSide <- 1 
        chTo <- c(grep("tori",as.character(addText["loc"])), grep("tole",as.character(addText["loc"])))
        if(length(chTo) >0) graphics::mtext(mTxt, at=FC, side=mSide, adj=0, col=mCol, cex=mCex, line=mLi) else {
          graphics::mtext(mTxt, side=mSide, adj=mAdj, col=mCol, cex=mCex, line=mLi)
        } }
    }
    chArr <- (arr[2] -arr[1]) > 0.05
    if(!chArr) arr[2] <- arr[1] +0.04
    ## draw arrow
    if(asX) graphics::arrows(FC, figCo[3] + arr[1]*(figCo[4]-figCo[3]), FC, figCo[3] + arr[2]*(figCo[4]-figCo[3]), 
      col=col,lwd=lwd,length=0.1) else { graphics::arrows(figCo[3] + arr[1]*(figCo[4]-figCo[3]), FC, 
      figCo[3] + arr[2]*(figCo[4]-figCo[3]), FC, col=col, lwd=lwd,length=0.1) }
    if(returnRatio) return(FC)
  } else if(!silent) message("unable to extract usable values for drawing arrow")
}
  
