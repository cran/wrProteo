#' Get Short Names of Proteomics Quantitation Software
#'
#' Get/convert short names of various proteomics quantitation software names.
#' A 2-letter abbreviation will be returned
#' 
#' @param x (character) 'mono' or 'average'
#' @param tryAsLower (logical) 
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a vector with masses for all amino-acids (argument 'massTy' to switch from mono-isotopic to average mass)
#' @seealso \code{\link{massDeFormula}}, \code{\link[wrMisc]{convToNum}}
#' @examples
#' shortSoftwName(c("maxquant","DIANN"))
#' @export
shortSoftwName <- function(x, tryAsLower=TRUE, silent=FALSE, debug=FALSE, callFrom=NULL)  {
    ## convert software-algorith names to 2-letter appreviation
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="shortwSoftName")
  if(isTRUE(debug)) silent <- FALSE
  if(!isTRUE(silent)) silent <- FALSE
  y <- cbind(softna=c("DIA-NN","ProteomeDiscoverer","Compomics","MaxQuant","Proline","TPP","FragPipe","MassChroQ","OpenMS"),
    shortna= c("DN","PD","CP","MQ","PL","TP","FP","MC","OM")  )
  out <- y[match(x, y[,1]), 2]
  chNa <- is.na(out)
  if(any(chNa) && tryAsLower) out[which(chNa)] <- y[match(tolower(sub("\\-","",x[which(chNa)])), tolower(y[,1])), 2]
  out }
    
