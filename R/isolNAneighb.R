#' Isolate NA-neighbours
#'
#' This functions extracts all replicate-values where at least one of the replicates is \code{NA} and sorts by number of \code{NA}s per group.
#' A list with all \code{NA}-neighbours organized by the number of \code{NA}s gets returned.
#'
#' @param mat (matrix or data.frame) main data (may contain \code{NA})
#' @param gr (character or factor) grouping of columns of 'mat', replicate association
#' @param silent (logical) suppress messages
#' @param debug (logical) display additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a list with NA-neighbours sorted by number of NAs in replicate group
#' @seealso This function gets used by \code{\link{matrixNAneighbourImpute}} and \code{\link{testRobustToNAimputation}}; estimation of mode \code{\link[wrMisc]{stableMode}}; detection of NAs \code{\link[stats]{na.fail}}
#' @examples
#' mat1 <- c(22.2, 22.5, 22.2, 22.2, 21.5, 22.0, 22.1, 21.7, 21.5, 22, 22.2, 22.7,
#'   NA, NA, NA, NA, NA, NA, NA, 21.2,   NA, NA, NA, NA,
#'   NA, 22.6, 23.2, 23.2,  22.4, 22.8, 22.8, NA,  23.3, 23.2, NA, 23.7,
#'   NA, 23.0, 23.1, 23.0,  23.2, 23.2, NA, 23.3,  NA, NA, 23.3, 23.8)
#' mat1 <- matrix(mat1, ncol=12, byrow=TRUE)
#' gr4 <- gl(3, 4)
#' isolNAneighb(mat1, gr4)
#' @export
isolNAneighb <- function(mat, gr, silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## isolate NA-neighbours
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="isolNAneighb")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  msg <- NULL
  datOK <- TRUE
  NAneig <- NULL
  if(any(length(mat) <1, length(dim(mat)) !=2, dim(mat) < c(2,1))) { datOK <- FALSE
    msg <- "'mat' should be matrix or data.frame with min 2 rows & 1 column, nothing to do, return NULL"}
  if(datOK && length(gr) !=ncol(mat)) { datOK <- FALSE
    msg <- "Length of 'gr' must match number of columns in 'mat', nothing to do, return NULL" }
  if(datOK && length(unique(gr)) ==length(gr)) { datOK <- FALSE
    msg <- "No replicates, can't isolate NA-neighbours" }
  if(!datOK && !silent) message(fxNa,msg)
  if(debug) {message(fxNa," iNN1"); iNN1 <- list(mat=mat, gr=gr, datOK=datOK)}

  if(datOK) {
    ## basic (optimized) extraction of NA-neighbours
    maxHi <- max(tapply(gr, gr, length)) -1             # max number of NA-neighbours (exclude group/line with all NA)
    NAneig <- lapply(1:maxHi, function(x) NULL)         # initialize output
    names(NAneig) <- paste0("n", 1:maxHi)
    ## need first to separate by groups of replicates
    matR <- lapply(unique(gr), function(x) {mat[, which(gr ==x)]})
    ## now separate NA-neighbours for each group/line
    nNA <- as.integer(sapply(matR, function(x) rowSums(is.na(x))))
    naNei <- wrMisc::partUnlist(lapply(matR, apply, 1, function(x) {chN <- is.na(x); if(sum(chN) ==0 || sum(chN)==length(x)) NULL else x[which(!chN)]}))
    ## combine according to number of NA-values in group/line
    for(i in 1:maxHi) {ch1 <- which(nNA==i); if(length(ch1) >0) NAneig[[i]] <- unlist(naNei[ch1])}}
  NAneig }
  

