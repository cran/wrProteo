#' Summarize On-Off Situations In Pairwise Comparisons Plus FDR
#'
#' This function aims to find instances from pairwise comparisons where one group of samples is entirely absent (at least \code{lim} percent of \code{NA}) 
#' and all present (no more than 1 - \code{lim} persent \code{NA}) in the second group, counts their frequency and calculates a False Discovery Rate (FDR) for this frequency of \code{NA}s.
#'  
#' @details
#' Cases where a given protein is found in all samples of a given group but not at all in another group is of great interest in pairwise comparisons, 
#' but typically difficult to compare based on statistical testing since missing values may render t-test impossible.
#' 
#' Here the frequency of all-missing and all-present is compared to the theoreticl expected frequency of such instances given the frequency of \code{NA}s in each sample.
#' This comparison can then be summarized as (global) FDR for the entire dataset. The FDR is given inside the '(' ')'. 
#' Please note that with a global FDR (like in this case with the Benjamini-Hochberg procedure used here) 
#' it is not possible to distinguish some specific lines as being more or less significant.
#' Thus, if the calculated FDR is below a given threshold (eg 0.05), on may conclude that significantly more on-off constellations were observed than expected solely on the frequency of \code{NA}s. 
#' 
#' This function uses internally \code{\link{testRobustToNAimputation}} if a $onOff is not yet present in the main object \code{stat}

#' @param stat ('MArrayLM'-object or list) designed for the output from \code{moderTest2grp} or \code{moderTestXgrp}
#' @param lim (numeric) the threshold to be applied on \code{statTy} for the result of the statistcal testing (after multiple testing correction)
#' @param addTy (character) influences output format : 'FDR' (default): concatenate real counts with FDR (n.theor / n.obs); 
#'   'paste' concatenate real counts with n.theor; 'table' table of observed and theor as separate columns
#' @param debug (logical) supplemental messages for debugging
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This function returns a character vector indicating the number of on-off constellations and the corresponding FDR.
#' @seealso  \code{\link{findOnOff}}, \code{\link{testRobustToNAimputation}}
#' @examples
#' dat1 <- list(raw=matrix(round(runif(120), 2), ncol=6))
#' dat1$raw[upper.tri(dat1$raw)] <- NA
#' dat1$raw[7:8,1:3] <- NA
#' dat1$isNA <- is.na(dat1$raw)
#' dat1$setup <- list(grp=rep(LETTERS[4:3], each=3), grpNa =c(LETTERS[4:3]), 
#'   pwGrpNa = matrix(LETTERS[4:3], nrow=1, dimnames=list("D  /  C",c("samp","ref"))) )
#' dat1$means <- matrix(runif(nrow(dat1$raw)*2), ncol=2, dimnames=list(NULL,dat1$setup$grpNa))
#' dat1$filter <- matrix(TRUE, nrow=nrow(dat1$raw), ncol=1)
#' 
#' head(resumeOnOff(dat1), 9)
#' 
#' @export
resumeOnOff <- function(stat, lim=0.9, addTy="FDR", silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## extract number of on/off events and calculate estimated (pure random) number for on/off 
  ## addTy (char) .. choose output format : 'FDR' : concatenate real counts with FDR (n.theor / n.obs); 'paste' concatenate real counts with n.theor; 'table' table of observed and theor as separate columns
  ##  buids on findOnOff()                       
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="findOnOff")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) { silent <- FALSE } else { debug <- FALSE }
  argNa <- deparse(substitute(stat))
  if(debug) { message(fxNa,"rOOff0")}
  if(!"list" %in% mode(stat) || length(stat) <1 && !silent) message(fxNa," Caution, '",argNa,"' must be a list or 'MArrayLM'-object from limma")
  if(length(lim) !=1 || !is.finite(lim) || lim >1) {lim <- 0.9; message(fxNa,"Invalid input for 'lim', setting to default (lim=0.9)")}
  ## main
  if(length(stat$onOff) ==0) stat <- findOnOff(stat, lim=lim, asFullObj=TRUE, silent=silent, debug=debug, callFrom=fxNa)      # observed on.off
  if(length(stat$onOff) !=0) {
    if(debug) { message(fxNa,"rOOff1")}
    frqNA <- colSums(stat$isNA) / nrow(stat$isNA) 
    pOff <- sapply(stat$setup$grpNa, function(x) prod(frqNA[which(stat$setup$grp %in% x)]))                             # p per sample pNA
    pOpw <- apply(stat$setup$pwGrpNa, 1, function(x) pOff[x[1]] *(1 - pOff[x[2]]) + (1 - pOff[x[1]]) *(pOff[x[2]]))     # p pairwise in both orientations
    if(debug) { message(fxNa,"rOOff2"); rOOff2 <- list(stat=stat,lim=lim,frqNA=frqNA,pOff=pOff,pOpw=pOpw ) }
    expOff <- round(pOpw * colSums(stat$filter), 1)                  # number of onOff expected at n.filtered repeats
    switch(addTy,
      paste = paste( colSums(stat$onOff), "/", expOff),                                                # for each pw question: absolute (/out of)
      FDR = paste0(colSums(stat$onOff), " (", round(expOff/(expOff + colSums(stat$onOff)), 3),")"),    # for each pw question: absolute (+FDR)
      table= cbind(observed=colSums(stat$onOff), theor=expOff))                                        # for each pw question: 2 col
  }    
}  

