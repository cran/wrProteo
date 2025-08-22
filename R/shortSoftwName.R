#' Get Short Names of Proteomics Quantitation Software
#'
#' Get/convert short names of various proteomics quantitation software names for software results handeled by this package.
#' A 2-letter abbreviation will be returned
#' 
#' @details 
#' So far thuis function recognizes the following software names:
#' "DIA-NN", "ProteomeDiscoverer", "Compomics", "MaxQuant", "Proline", "TPP", "FragPipe", "MassChroQ", "OpenMS", "Ionbot" and "Sage"
#' 
#' @param x (character) software (full) name
#' @param tryAsLower (logical) include lower-caps writing to search 
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a vector with 2-letter abbreviation for the software
#' @seealso \code{\link{readMaxQuantFile}}
#' @examples
#' shortSoftwName(c("maxquant","DIANN"))
#' @export
shortSoftwName <- function(x, tryAsLower=TRUE, silent=FALSE, debug=FALSE, callFrom=NULL)  {
    ## convert software-algorith names to 2-letter appreviation
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="shortwSoftName")
  if(isTRUE(debug)) silent <- FALSE
  if(!isTRUE(silent)) silent <- FALSE
  y <- cbind(softna=c("DIA-NN","ProteomeDiscoverer","Compomics","MaxQuant","Proline","TPP","FragPipe","MassChroQ","OpenMS","Ionbot","Sage"),
    shortna= c("DN","PD","CP","MQ","PL","TP","FP","MC","OM","IB","SA")  )
  out <- y[match(x, y[,1]), 2]
  chNa <- is.na(out)
  if(any(chNa) && tryAsLower) out[which(chNa)] <- y[match(tolower(sub("\\-","",x[which(chNa)])), tolower(y[,1])), 2]
  out }
    
