#' Extract species annotation
#'
#' \code{extrSpeciesAnnot} identifies species-related annotation (as suffix to identifyers) for data comnining multiple species and returns alternative (short) names.  
#' This function also suppresses extra heading or tailing space or punctuation characters.
#' In case multiple tags are found, the last tag is reported and a message of alert may be displayed.  
#' 
#' @param annot (character) vector with initial annotation
#' @param spec (character) the tags to be identified
#' @param shortNa (character) the final abbreviation used, order and lengt must fit to argument \code{annot}
#' @param silent (logical) suppress messages
#' @param debug (logical) display additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a character vector with single (last of multiple) term if found in argument \code{annot}
#' @seealso  \code{\link[base]{grep}} 
#' @examples
#' spec <- c("keratin_CONT","AB_HUMAN","CD_YEAST","EF_G_HUMAN","HI_HUMAN_ECOLI","_YEAST_012")
#' extrSpeciesAnnot(spec) 
#' @export
extrSpeciesAnnot <- function(annot, spec=c("_CONT","_HUMAN","_YEAST","_ECOLI"), shortNa=c("cont","H","S","E"), silent=FALSE, debug=FALSE, callFrom=NULL){
  ## extract species information for element of 'annot'
  ## return character vector with single (last of) term if found in 'annot'
  ## 'annot' .. character vector
  ## 'spec' ..  (character) term to search
  ## 'shortextrSpeciesAnnotNa' .. (character) term to code output
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="extrSpeciesAnnot")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE

  msg <- "Argument 'shortNa' doesn't fit to length of 'spec'"
  if(length(shortNa) < length(spec) && length(shortNa) >0) {
   if(!silent) message(fxNa,msg," ignoring")
   shortNa <- NULL }
  if(is.null(shortNa)) {
    shortNa <- sub("[[:punct:]]+[[:blank:]]+[[:punct:]]+|^[[:blank:]]+[[:punct:]]+[[:blank:]]+","",spec)
    trim <- substr(shortNa,1,1)
    if(length(unique(trim)) < length(spec)) {
      shortNa <- substr(shortNa,1,2)
    } else shortNa <- trim
    if(!silent) message(fxNa,"Constructing 'shortNa'.. replace by 1st alphanum-character : ",shortNa) }
  ## main
  tmp <- list()
  out <- rep(NA, length(annot))
  for(i in 1:length(spec)) {tmp[[i]] <- grep(spec[i], annot)
    out[tmp[[i]]] <- shortNa[i] }
  che <- table(table(unlist(tmp)))
  if(any(as.numeric(names(che)) >1) && !silent) message(fxNa,"Multiple/conflicting annotation in ",sum(che[which(as.numeric(names(che)) != 1)])," cases")  
  out }
    
