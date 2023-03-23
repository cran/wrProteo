#' Deprecialed Volcano-plot 
#'
#' Please use VolcanoPlotW() from package wrGraph.
#' This function does NOT produce a plot any more.
#'	 
#' @param Mvalue (numeric or matrix) data to plot; M-values are typically calculated as difference of log2-abundance values and 'pValue' the mean of log2-abundance values;
#'   M-values and p-values may be given as 2 columsn of a matrix, in this case the argument \code{pValue} should remain NULL 
#' @param pValue (numeric, list or data.frame) if \code{NULL} it is assumed that 2nd column of 'Mvalue' contains the p-values to be used
#' @param useComp (integer, length=1) choice of which of multiple comparisons to present in \code{Mvalue} (if generated using \code{moderTestXgrp()})  
#' @param filtFin (matrix or logical) The data may get filtered before plotting: If \code{FALSE} no filtering will get applied; if matrix of \code{TRUE}/\code{FALSE} it will be used as optional custom filter, otherwise (if \code{Mvalue} if an \code{MArrayLM}-object eg from limma) a default filtering based on the \code{filtFin} element will be applied 
#' @param ProjNa (character) custom title
#' @param FCthrs (numeric) Fold-Change threshold (display as line) give as Fold-change and NOT log2(FC), default at 1.5, set to \code{NA} for omitting
#' @param FdrList (numeric) FDR data or name of list-element
#' @param FdrThrs (numeric) FDR threshold (display as line), default at 0.05, set to \code{NA} for omitting 
#' @param FdrType (character) FDR-type to extract if \code{Mvalue} is 'MArrayLM'-object (eg produced by from \code{moderTest2grp} etc);
#'   if \code{NULL} it will search for suitable fields/values in this order : 'FDR','BH',"lfdr" and 'BY'
#' @param subTxt (character) custom sub-title
#' @param grayIncrem (logical) if \code{TRUE}, display overlay of points as increased shades of gray
#' @param col (character) custom color(s) for points of plot (see also \code{\link[graphics]{par}})
#' @param pch (integer) type of symbol(s) to plot (default=16) (see also \code{\link[graphics]{par}}) 
#' @param compNa (character) names of groups compared
#' @param batchFig (logical) if \code{TRUE} figure title and axes legends will be kept shorter for display on fewer splace  
#' @param cexMa (numeric) font-size of title, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}})
#' @param cexLa (numeric) size of axis-labels, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}}) 
#' @param limM (numeric, length=2) range of axis M-values 
#' @param limp (numeric, length=2) range of axis FDR / p-values
#' @param annotColumn (character) column names of annotation to be extracted (only if \code{Mvalue} is \code{MArrayLM}-object containing matrix $annot).
#'   The first entry (typically 'SpecType') is used for different symbols in figure, the second (typically 'GeneName') is used as prefered text for annotating the best points (if \code{namesNBest} allows to do so.)
#' @param annColor (character or integer) colors for specific groups of annoatation (only if \code{Mvalue} is \code{MArrayLM}-object containing matrix $annot)
#' @param cexPt (numeric) size of points, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}})
#' @param cexSub (numeric) size of subtitle, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}})
#' @param cexTxLab (numeric) size of text-labels for points, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}})
#' @param namesNBest (integer or character) number of best points to add names in figure; if 'passThr' all points passing FDR and FC-filtes will be selected; 
#'   if the initial object \code{Mvalue} contains a list-element called 'annot' the second of the column specified in argument \code{annotColumn} will be used as text
#' @param NbestCol (character or integer) colors for text-labels of best points
#' @param sortLeg (character) sorting of 'SpecType' annotation either ascending ('ascend') or descending ('descend'), no sorting if \code{NULL}
#' @param NaSpecTypeAsContam (logical) consider lines/proteins with \code{NA} in Mvalue$annot[,"SpecType"] as contaminants (if a 'SpecType' for contaminants already exits)
#' @param useMar (numeric,length=4) custom margings (see also \code{\link[graphics]{par}})
#' @param returnData (logical) optional returning data.frame with (ID, Mvalue, pValue, FDRvalue, passFilt) 
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @param debug (logical) additional messages for debugging 
#' @return deprecated - returns nothing
#' @seealso this function was replaced by \code{\link[wrGraph]{plotPCAw}})
#' @examples
#' set.seed(2005); mat <- matrix(round(runif(900),2), ncol=9)
#' @export
VolcanoPlotW2 <- function(Mvalue, pValue=NULL, useComp=1, filtFin=NULL, ProjNa=NULL, FCthrs=NULL, FdrList=NULL, FdrThrs=NULL, FdrType=NULL,
  subTxt=NULL, grayIncrem=TRUE, col=NULL, pch=16, compNa=NULL, batchFig=FALSE, cexMa=1.8, cexLa=1.1, limM=NULL, limp=NULL,
  annotColumn=NULL, annColor=NULL, cexPt=NULL, cexSub=NULL, 
  cexTxLab=0.7, namesNBest=NULL, NbestCol=1, sortLeg="descend", NaSpecTypeAsContam=TRUE, useMar=c(6.2,4,4,2), returnData=FALSE, callFrom=NULL, silent=FALSE,debug=FALSE) {
  ## MA plot
  ## optional arguments for explicit title in batch-mode
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="VolcanoPlotW2")
  .Deprecated("++ NOTE : OLD VERSION !! ++   Please use VolcanoPlotW() from package wrGraph")
 }
   
