#' Isolate NA-neighbours
#' 
#' This functions extracts all replicate-values where at least one of the replicates is \code{NA}. 
#' Then, the non-\code{NA} values are sorted by the number of \code{NA}s which occored in this group of replicates.
#' A list with all \code{NA}-neighbours organized by the number of \code{NA}s gets returned.           
#' 
#' @param mat (matrix or data.frame) main data (may contain \code{NA})
#' @param gr (character or factor) grouping of columns of 'mat', replicate association
#' @param maxHi (integer) maximum count of NAs to consider separately (higher ones will be counted/pooled as maxHi)
#' @param iniCheck (logical) check at beginning if executing this function is useful (presence any \code{NA})
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @return list with NA-neighbours sorted by number of NAs in replicate group
#' @seealso this function gets used by \code{\link{matrixNAneighbourImpute}} and \code{\link{testRobustToNAimputation}}; estimation of mode \code{\link[wrMisc]{stableMode}}; detection of NAs \code{\link[stats]{na.fail}}
#' @examples
#' mat1 <- c(22.2, 22.5, 22.2, 22.2, 21.5, 22.0, 22.1, 21.7, 21.5, 22, 22.2, 22.7,
#'   NA, NA, NA, NA, NA, NA, NA, 21.2,   NA, NA, NA, NA,
#'   NA, 22.6, 23.2, 23.2,  22.4, 22.8, 22.8, NA,  23.3, 23.2, NA, 23.7,
#'   NA, 23.0, 23.1, 23.0,  23.2, 23.2, NA, 23.3,  NA, NA, 23.3, 23.8)
#' mat1 <- matrix(mat1, ncol=12, byrow=TRUE)
#' gr4 <- gl(3, 4)
#' isolNAneighb(mat1, gr4)
#' 
#' @export
isolNAneighb <- function(mat, gr, maxHi=3, iniCheck=TRUE, silent=FALSE, callFrom=NULL) {
  ## isolate NA-neighbours
  ## 'maxHi' .. maximum count of NAs to consider separately (higher ones will be counted/pooled as maxHi)
  if(!isTRUE(silent)) silent <- FALSE
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="isolNAneighb")
  datOK <- TRUE
  msg <- NULL
  if(!isTRUE(silent)) silent <- FALSE
  if(any(length(mat) <1, length(dim(mat)) !=2, dim(mat) < c(2,1))) { datOK <- FALSE
    msg <- "'mat' should be matrix or data.frame with min 2 rows & 1 col , nothing to do, return NULL"}
  if(length(gr) !=ncol(mat)) { datOK <- FALSE
    msg <- "length of 'gr' must match number of columns in 'mat', nothing to do, return NULL" } 
  
  if(datOK) {
    NAneig <- lapply(1:maxHi, function(x) NULL)         # initialize
    names(NAneig) <- paste0("n", 1:maxHi)
    chNa <- if(!isFALSE(iniCheck)) is.na(mat) else TRUE
    if(!any(chNa)) { return(NAneig)
    } else { for(i in wrMisc::naOmit(unique(gr))) {
      tmp <- mat[,which(gr==i)]
      nNA <- rowSums(is.na(tmp))
      chHi <- nNA > maxHi & nNA < ncol(tmp)
      if(any(chHi)) nNA[which(chHi)] <- maxHi
      chExtr <- nNA %in% c(0,ncol(tmp))
      if(any(chExtr)) nNA[which(chExtr)] <- NA
      NAnei <- by(tmp,nNA,function(x) wrMisc::naOmit(as.numeric(as.matrix(x))))
      for(j in names(NAnei)) { x <- as.numeric(j); NAneig[[x]] <- c(NAneig[[x]], NAnei[[j]])}
    } }
    NAneig 
  } else { if(!silent) message(fxNa, msg)
    NULL} }
  
#' @export
.nNAbyGroup <- function(mat, gr) {
  ## get number of NAs per line & group of replicates
  ## replaced by wrMisc::rowGrpNA
  gr1 <- wrMisc::naOmit(unique(gr))
  nNA <- matrix(nrow=nrow(mat), ncol=length(gr1), dimnames=list(NULL,gr1))
  for(i in 1:length(gr1)) {
    nNA[,i] <- rowSums(is.na(mat[,which(gr==gr1[i])])) }
  nNA }  
  
