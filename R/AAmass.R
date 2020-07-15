#' Molecular mass for amino-acids 
#'
#' Calculate molecular mass based on atomic composition    
#' 
#' @param massTy (character) 'mono' or 'average'
#' @param inPept (logical) remove H20 corresponding to water loss at peptide bond formaton
#' @param inclSpecAA (logical) include ornithine O & selenocysteine U 
#' @return vector with masses for all amino-acids (argument 'massTy' to switch form mono-isotopic to average mass)
#' @seealso \code{\link{massDeFormula}}, \code{\link[wrMisc]{convToNum}}
#' @examples
#' massDeFormula(c("12H12O","HO"," 2H 1 Se, 6C 2N","HSeCN"," ","e"))
#' AAmass()
#' @export
AAmass <- function(massTy="mono",inPept=TRUE,inclSpecAA=FALSE) {
  ## return vector with masses for all amino-acids (argument 'massTy' to switch form mono-isotopic to average mass)
  ## 'inPept' will remove H20 corresponding to water loss at peptide bond formaton
  ## 'inclSpecAA' .. include ornithine O & selenocysteine U
  ##  so far all LETTERS exept B,J,X,Z (ie 2,10,24,26)   spec (used) O,U
  msg <- " argument 'massTy' must bei either  'mono' or 'average' !" 
  chTy <- length(massTy)
  if(chTy <1) stop(msg) else massTy <- massTy[1]
  chTy <- c("mono","average") %in% massTy
  if(!any(chTy)) stop(msg)  
  aaComp <- cbind(C=c(3,6,4,4,3,5,5,2,6,6,6,6,5,9,5,3,4,11,9,5,5,3),
    H=c(5,12,6,5,5,7,8,3,7,11,11,12,9,9,7,5,7,10,9,9,12,5),
    O=c(1,1,2,3,1,3,2,1,1,1,1,1,1,1,1,2,2,1,2,1,2,1),
    N=c(1,4,2,1,1,1,2,1,3,1,1,2,1,1,1,1,1,2,1,1,2,1),
    S=c(0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0),
    Se=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1))
  rownames(aaComp) <- c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V","O","U")
  if(!inPept) aaComp[,2:3] <- aaComp[,2:3] + matrix(rep(2:1, each=nrow(aaComp)), ncol=2)
  atoMass <- .atomicMasses()[,massTy]
  AAmass <- aaComp*matrix(rep(atoMass[match(colnames(aaComp),names(atoMass))], each=nrow(aaComp)), nrow=nrow(aaComp))
  AAmass <- rowSums(AAmass)
  if(!inclSpecAA) AAmass <- AAmass[1:20]  # so far exclude ornithine O & selenocysteine U
  AAmass }

