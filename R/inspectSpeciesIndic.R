#' Inspect Species Indictaion Or Group of Proteins
#'
#' This function inspects its main argument to convert a species indication to the scientific name or to return all protein-accession numbers for a name of a standard collection like UPS1.
#' 
#' 
#' 
#' @param x (character) species indication or name of collection of proteins (so far only UPS1 & UPS2)
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of messages produced
#' @param debug (logical) display additional messages for debugging
#' @return This function returns a character vector
#' @seealso  \code{\link{getUPS1acc}};
#' @examples
#' inspectSpeciesIndic("Human")
#' inspectSpeciesIndic("UPS1")
#' @export
inspectSpeciesIndic <- function(x, silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## inspect species indication and convert to scientific name, in case of UPS1 return all accessions
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="inspectSpeciesIndic")
  out <- NULL
  if(length(x) <1 || all(is.na(x))) warning("Invalid entry of 'x'") else {
    if(length(x) >1) {x <- as.character(x[1]); warning(fxNa,"'x' may only be of length=1, truncating")}
    commSpec <- .commonSpecies()
    commSpec[,1] <- sub("^_","", commSpec[,1])
    x2 <- tolower(x)
    # as is #
    ch1 <- x %in% commSpec
    if(!any(ch1, na.rm=TRUE)) ch1 <- tolower(x) %in% tolower(commSpec)  ## additional search with all as lower caps
    if(any(ch1, na.rm=TRUE)) { useLi <- which(commSpec==x, arr.ind=TRUE)[1]
      out <- commSpec[useLi, 2]           # convert to scientific name
      names(out) <- commSpec[useLi, 3]
    } else {
      if(any(c("UPS1","UPS-1","UPS2","UPS-2") ==x, na.rm=TRUE)) {out <- wrProteo::getUPS1acc()$ac; names(out) <- wrProteo::getUPS1acc()$uniProt}
    }
    if(length(out) <1) {if(!silent) message(fxNa,"Unknow species indication '",x,"', trying to use as is"); out <- x}
  }
  out }
