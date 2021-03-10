#' Filter based on either number of total peptides and specific peptides or number of razor petides
#'
#' \code{razorNoFilter} filters based on either a) number of total peptides and specific peptides or b) numer of razor petides.
#' This function was designed for filtering using a mimimum number of (PSM-) count values following the common practice to consider results with 2 or more peptide counts as reliable. 
#' The function be (re-)run independently on each of various questions (comparisons).
#' Note: Non-integer data will be truncated to integer (equivalent to  \code{floor}). 
#'  
#' @param annot (matrix or data.frame) main data (may contain NAs) with (PSM-) count values for each protein
#' @param speNa (integer or character) indicate which column of 'annot' has number of specific peptides
#' @param totNa (integer or character) indicate which column of 'annot' has number of total peptides
#' @param minRazNa (integer or character) name of column with number of razor peptides, alternative to 'minSpeNo'& 'minTotNo'
#' @param minSpeNo (integer) minimum number of pecific peptides
#' @param minTotNo (integer) minimum total ie max razor number of peptides
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of messages produced
#' @return vector of logical values if corresponding line passes filter criteria  
#' @seealso \code{\link[wrMisc]{presenceFilt}} 
#' @examples
#' set.seed(2019); datT <- matrix(sample.int(20,60,replace=TRUE),ncol=6,
#'   dimnames=list(letters[1:10],LETTERS[1:6])) -3
#' datT[,2] <- datT[,2] +2
#' datT[which(datT <0)] <- 0
#' razorNoFilter(datT,speNa="A",totNa="B")
#' @export
razorNoFilter <- function(annot,speNa=NULL,totNa=NULL,minRazNa=NULL,minSpeNo=1,minTotNo=2,silent=FALSE,callFrom=NULL) {
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="razorNoFilter")
  if(is.null(minRazNa)) {
    specPe <- as.integer(annot[,speNa]) >= minSpeNo
    totPe <- as.integer(annot[,totNa]) >= minTotNo
    filt <- (specPe & totPe)  
  } else {
    filt <- as.integer(annot[,minRazNa]) >= minTotNo
  }
  filt } 

#' @export
.checkKnitrProt <- function(tryF=FALSE) {
  ## function for checking presence of knitr and rmarkdown
  ## needed to explicitely call functions of packages
  chPaR <- try(find.package("rmarkdown"), silent=TRUE)
  chPaK <- try(find.package("knitr"), silent=TRUE)
  if("try-error" %in% class(chPaR)) warning("package 'rmarkdown' not found ! Please install from CRAN") else {
    if(tryF) rmarkdown::pandoc_available() }
  if("try-error" %in% class(chPaK)) warning("package 'knitr' not found ! Please install from CRAN") else {
    if(tryF) knitr::kable(matrix(1:4, ncol=2)) }
  }  
   
