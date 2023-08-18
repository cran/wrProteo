#' AUC from ROC-curves
#'
#' This function calculates the AUC (area under the curve) from ROC data in matrix of specificity and sensitivity values,
#' as provided in the output from  \code{\link{summarizeForROC}}.
#'
#' @param dat (matrix or data.frame) main inut containig sensitivity and specificity data (from \code{summarizeForROC})
#' @param useCol (character or integer) column names to be used: 1st for specificity and 2nd for sensitivity count columns
#' @param returnIfInvalid (\code{NA} or \code{NULL}) what to return if data for calculating ROC is invalid or incomplete
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This function returns a matrix including imputed values or list of final and matrix with number of imputed by group (plus optional plot)
#' @seealso preparing ROC data \code{\link{summarizeForROC}}, (re)plot the ROC figure \code{\link{plotROC}};
#'   note that numerous other packages also provide support for working with ROC-curves : Eg \href{https://CRAN.R-project.org/package=dlstats}{rocPkgShort},
#'    \href{https://CRAN.R-project.org/package=ROCR}{ROCR}, \href{https://CRAN.R-project.org/package=pROC}{pROC} or \href{https://CRAN.R-project.org/package=ROCit}{ROCit}, etc.
#' @examples
#' set.seed(2019); test1 <- list(annot=cbind(Species=c(rep("b",35), letters[sample.int(n=3,
#'   size=150,replace=TRUE)])), BH=matrix(c(runif(35,0,0.01), runif(150)), ncol=1))
#' roc1 <- summarizeForROC(test1, spec=c("a","b","c"), annotCol="Species")
#' AucROC(roc1)
#' @export
AucROC <- function(dat, useCol=c("spec","sens"), returnIfInvalid=NA, silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## calculate AUC (area under the curve) from ROC data
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="AucROC")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  if(!is.null(returnIfInvalid)) returnIfInvalid <- NA
  ## basic checks
  datOK <- length(dat) >0 && length(useCol) >1
  if(datOK)  { chDim <- dim(dat)
    if(any(length(dim(dat)) !=2, dim(dat) < 1:2)) datOK <- FALSE }   # must be matrix or data.frame with min 1 lines and 2 cols")
  if(datOK) { if(is.numeric(useCol)) { if(any(useCol[1:2] <1) || any(useCol[1:2] > ncol(dat))) datOK <- FALSE
    } else datOK <- all(useCol %in% colnames(dat)) }
  if(debug) {message(fxNa,"aucR1   datOK : ",datOK); aucR1 <- list(dat=dat,useCol=useCol)}
  if(datOK) {
      if(is.character(useCol) && length(useCol) >1) {
      useCol <- which(colnames(dat) %in% useCol[1:2])
      dataOK <- !any(is.na(useCol)) }
  }
  out <- returnIfInvalid
  if(datOK) {
    ## check for NA
    chNa <- is.na(dat[,useCol])
    if(debug) {message(fxNa,"aucR2"); aucR2 <- list(dat=dat,useCol=useCol,chNa=chNa, datOK=datOK)}

    if(any(chNa)) { ch1 <- is.na(dat[1, useCol])
      if(!silent) message(fxNa,"NOTE : the data contain ",sum(chNa)," NAs, replacing by preceeding value")
      if(any(ch1)) dat[1,useCol] <- c(if(ch1[1]) 1 else dat[1,useCol[1]], if(ch1[2]) 0 else dat[1,useCol[2]])
      for(i in 1:2) {ch1 <- is.na(dat[,useCol[i]])
        if(any(ch1)) dat[which(ch1), useCol[i]] <- dat[which(ch1) -1, useCol[i]]} }
    if(debug) {message(fxNa,"aucR4")}
    ## Normalize (if needed)
    for(i in 1:2)  if(max(dat[,useCol[i]]) >1) dat[,useCol[i]] <- dat[,useCol[i]]/max(dat[,useCol[i]])
    ##
    if(dataOK) out <- sum(abs(diff(dat[,useCol[1]])) *dat[-nrow(dat), useCol[2]]) else {
      if(!silent) message(fxNa,"Invalid input / nothing to do") } }
  out }
  
