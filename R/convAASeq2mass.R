#' Molecular mass for amino-acids 
#'
#' This function calculates the molecular mass of one-letter code amion-acid sequences.    
#' 
#' @param x (character) aminoacid sequence (single upper case letters for describing a peptide/protein)
#' @param massTy (character) default 'mono' for mono-isotopic masses (alternative 'average')
#' @param seqName (logical) optional (alternative) names for the content of 'x' (ie aa seq) as name (always if 'x' has no names)
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @return This functions returns a vector with masses for all amino-acids (argument 'massTy' to switch form mono-isotopic to average mass)
#' @seealso \code{\link{massDeFormula}}, \code{\link{AAmass}}, \code{\link[wrMisc]{convToNum}}
#' @examples
#' convAASeq2mass(c("PEPTIDE","fPROTEINES"))
#' pep1 <- c(aa="AAAA", de="DEFDEF")
#' convAASeq2mass(pep1, seqN=FALSE)
#' @export
convAASeq2mass <- function(x, massTy="mono", seqName=TRUE, silent=FALSE, callFrom=NULL) {
  ## convert (character) aminoacid sequence vector (ie AA seq in single upper case letters) to mass with corresp modif
  ## 'seqName'  .. to use 'x' (aa seq) as name (always if 'x' has no names)
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="convAASeq2mass")
  AAmass1 <- AAmass(massTy=massTy ,inPept=TRUE)
  mH20 <- massDeFormula("2HO", massTy=massTy)
  if(length(names(x)) <1) seqName <- TRUE
  chNoLet <- which(!LETTERS %in% names(AAmass1))
  chNoLe2 <- nchar(x) == nchar(sapply(LETTERS[chNoLet],gsub,"",x))
  chNoLe2 <- if(length(dim(chNoLe2)) >1) colSums(!chNoLe2) >0 else !chNoLe2
  if(any(chNoLe2)) warning(fxNa,"Encountered/ignoring non-attributed sequence character: ",wrMisc::pasteC(LETTERS[chNoLet][which(chNoLe2)],quoteC="'")," !!")
  pep1 <- lapply(strsplit(x,""), match, names(AAmass1))                              # transform into indexes of AA-letters
  out <- sapply(pep1, function(x) sum(AAmass1[x], mH20, na.rm=TRUE))                 # basic mass (as sum of its AA)
  names(out) <- if(seqName) x else names(x)
  chNa <- names(out) %in% "0z" 
  if(any(chNa)) names(out)[which(chNa)] <- ""
  out }
  
