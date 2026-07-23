#' Add arrow for expected Fold-Change to VolcanoPlot or MA-plot  
#'
#' NOTE : This function is deprecated, please use \code{\link[wrGraph]{foldChangeArrow}} instead !!
#' This function was made for adding an arrow indicating a fold-change to MA- or Volcano-plots. 
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
#' @param silent (logical) suppress messages
#' @param debug (logical) display additional messages for debugging
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @return plots arrow only (and explicative text), if \code{returnRatio=TRUE} also returns numeric value for extracted ratio
#'
#' @details The argument \code{addText} also allows specifying a fixed position when using \code{addText=c(loc="bottomleft")}, also bottomright, topleft, topright, toleft and toright may be used.
#'   In this case the elemts \code{side} and \code{adjust} will be redefined to accomodate the text in the corner specified. 
#'
#'  Ultimately this function will be integated to the package wrGraph. 
#'
#' @seealso new version : \code{\link[wrGraph]{foldChangeArrow}}; used with \code{\link[wrGraph]{MAplotW}}, \code{\link[wrGraph]{VolcanoPlotW}}
#' @examples
#' #deprecated# foldChangeArrow2(FC=1.5) 
#' 
#' @export
foldChangeArrow2 <- function(FC, useComp=1, isLin=TRUE, asX=TRUE, col=1, arr=c(0.005,0.15), lwd=NULL, 
  addText=c(line=-0.9,cex=0.7,txt="expected",loc="toright"), returnRatio=FALSE, silent=FALSE, debug=FALSE, callFrom=NULL){
  ##
  .Deprecated(new="foldChangeArrow", package="wrGraph", msg="Please use the foldChangeArrow() function form the package wrGraph instead !")
}
  
