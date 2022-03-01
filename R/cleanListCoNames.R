#' Selective batch cleaning of sample- (ie column-) names in list  
#'
#' This function allows to manipulate sample-names (ie colnames) from data stored as multiple matrixes or data.frames in multiple sheets of a list in a batch-wise manner.
#' Import functions such as \code{readMaxQuantFile()} organize initial flat files into lists (of matrixes) of the different types of data.
#' Many times all column names in such lists carry long names including redundant information, like the overall experiment name or date, etc. 
#' The aim of this function is to facilitate 'cleaning' the sample- (ie column-) names to obtain short and concise names.
#' Character terms to be removed (via argument \code{rem}) and/or replaced/subsitituted (via argument \code{subst}) should be given as they are, characters with special behaviour in \code{grep} (like '.') will be protected internally.
#' Note, that the character substitution part will be done first, and the removal part (without character replacement) afterwards.
#'
#' @param dat (list) main input
#' @param rem (character) character string to be removed, may be named 'left' and 'right' for more specific exact pattern matching 
#'   (this part will be perfomed before character substitutions by \code{subst})
#' @param subst (character of length=2, or matrix with 2 columns) pair(s) of character-strings for replacement (1st as search-item and 2nd as replacement); this part is performed after character-removal via \code{rem} 
#' @param lstE (character, length=1) names of list-elements where colnames should be cleaned
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of messages produced
#' @seealso \code{\link[base]{grep}}  
#' @return This function returns a list (equivalent to input \code{dat})
#' @examples
#' dat1 <- matrix(1:12, ncol=4, dimnames=list(1:3, paste0("sample_R.",1:4)))
#' dat1 <- list(raw=dat1, quant=dat1, notes="other..") 
#' cleanListCoNames(dat1, rem=c(left="sample_"), c(".","-")) 
#' @export
cleanListCoNames <- function(dat, rem=NULL, subst=c("-","_"), lstE=c("raw","quant","counts"), silent=FALSE, callFrom=NULL) {
  ## clean/stratify columnames
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="cleanListCoNames")
  ch1 <- lstE %in% names(dat)
  if(!is.list(dat) | all(!ch1)) {if(!silent) message(fxNa," Nothing to do (verifiy your input)"); ok <- FALSE} else ok <- TRUE 
  if(ok & any(!ch1)) { if(!silent) message(fxNa," Term ",wrMisc::pasteC(lstE[which(!ch1)],quoteC="'")," not found in 'dat', ignoring ..")
    lstE <- lstE[which(ch1)]  }
    
  needToProt <- paste0("\\",c(".","+","*", "^","$","?", "(",")","\\"))                                     #" #for editor
  if(ok & length(rem) >0) {
    ## protect special characters (if needed)
    for(i in needToProt)  rem <- gsub(i, paste0("\\",i), rem)                                              #" #for editor
    ## look for left side removal
    ch1 <- names(rem) %in% c("l","le","left")
    if(any(ch1)) rem[which(ch1)] <- paste0("^",rem[which(ch1)])
    ## look for right side removal
    ch1 <- names(rem) %in% c("r","ri","right")
    if(any(ch1)) rem[which(ch1)] <- paste0(rem[which(ch1)],"$")
    for(i in lstE) for(j in rem) colnames(dat[[i]]) <- sub(j,"",colnames(dat[[i]]))
  }
  if(ok & length(subst) >0) {
    ## character substitution
    if(length(dim(subst)) >1) if(any(dim(subst) < 1:2)) { subst <- NULL
      message(fxNa," Invalid argument 'subst' (should be matrix with left column for term to search for and right column with replacement-term), ignoring") }
    if(length(dim(subst)) >1) {  
      for(j in needToProt) subst[,1] <- gsub(j, paste0("\\",i), subst[,1])                                              #" #for editor
      for(i in lstE) for(j in 1:nrow(subst)) colnames(dat[[i]]) <- sub(subst[j,1], subst[j,2], colnames(dat[[i]]))
    } else if(length(subst) >1) {
      for(j in needToProt) subst[1] <- gsub(j, paste0("\\",i), subst[1])                      #" #for editor
      for(i in lstE) 
      colnames(dat[[i]]) <- sub(subst[1], subst[2], colnames(dat[[i]]))}
  }
  dat }
  
