#' UniProt Accession-Numbers And Names Of UPS1 Proteins
#'
#' UPS1 (see https://www.sigmaaldrich.com/FR/en/product/sigma/ups1) and UPS2 are commerical products consisting of a mix of 48 human (purified) proteins.
#' They are frequently used as standard in spike-in experiments, available from Sigma-Aldrich (\code{https://www.sigmaaldrich.com/GB/en}).
#' This function allows accessing their protein accession numbers and associated names on \href{https://www.uniprot.org/}{UniProt}
#'  
#' @details
#' Please note that the UniProt accession 'P62988' for 'UBIQ_HUMAN' (as originally cited by Sigma-Aldrich)
#' has been withdrawn and replaced in 2010 by \href{https://www.uniprot.org/}{UniProt} by the accessions 'P0CG47', 'P0CG48', 'P62979', and 'P62987'.
#' This initial accession is available via \code{getUPS1acc()$acOld}, now \code{getUPS1acc()$ac} contains 'P0CG47'. 
#'  
#' 
#' @param updated (logical) return updated accession number (of UBB)
#' @return This function returns data.frame with accession-numbers as stated by the supplier (\code{$acFull}),
#'  trimmed accession-numbers, ie without version numbers (\code{$ac}), 
#'  and associated (\code{UniProt}) entry-names  (\code{$EntryName}) from \href{https://www.uniprot.org/}{UniProt} 
#'  as well as the species designation for the collection of 48 human UPS1 or UPS2 proteins.
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allows easier tracking of messages produced
#' 
#' @return This function returns a matrix including imputed values or list of final and matrix with number of imputed by group (plus optional plot)
#' @examples
#' head(getUPS1acc())
#' @export
getUPS1acc <- function(updated=TRUE, silent=FALSE, debug=FALSE, callFrom=NULL) {
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="getUPS1acc")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  ## The UniProt accession numbers for the UPS1 proteins
  UPS1 <- data.frame( ac=rep(NA,48),
  acFull=c("P00915", "P00918", "P01031", "P69905", "P68871", "P41159", "P02768", "P0CG47", "P04040", 
    "P00167", "P01133", "P02144", "P15559", "P62937", "Q06830", "P63165", "P00709", "P06732", 
    "P12081", "P61626", "Q15843", "P02753", "P16083", "P63279", "P01008", "P61769", "P55957", 
    "O76070", "P08263", "P01344", "P01127", "P10599", "P99999", "P06396", "P09211", "P01112", 
    "P01579", "P02787", "O00762", "P51965", "P08758", "P02741", "P05413", "P10145", "P02788", 
    "P10636-8", "P00441", "P01375"),
  acNew=NA,
  acOld=NA,
  Gene=c("CAH1", "CAH2", "CO5", "HBA", "HBB", "LEP", "ALBU", "UBB", "CATA",
    "CYB5", "EGF", "MYG", "NQO1", "PPIA", "PRDX1", "SUMO1", "LALBA", "KCRM",
    "HARS1", "LYSC", "NEDD8", "RBP4", "NQO2", "UBC9", "ANT3", "B2MG", "BID",
    "SYUG", "GSTA1", "IGF2", "PDGFB", "THIO", "CYC", "GELS", "GSTP1", "RASH",
    "IFNG", "TRFE", "UBE2C", "UB2E1", "ANXA5", "CRP", "FABPH", "IL8", "TRFL",
    "TAU", "SODC", "TNFA"),
  species=rep("Homo sapiens", 48),
  EntryName=c("CAH1_HUMAN","CAH2_HUMAN","CO5_HUMAN","HBA_HUMAN","HBB_HUMAN", "LEP_HUMAN","ALBU_HUMAN","UBB_HUMAN","CATA_HUMAN",
    "CYB5_HUMAN","EGF_HUMAN","MYG_HUMAN","NQO1_HUMAN","PPIA_HUMAN","PRDX1_HUMAN", "SUMO1_HUMAN","LALBA_HUMAN","KCRM_HUMAN",
    "HARS1_HUMAN","LYSC_HUMAN","NEDD8_HUMAN","RET4_HUMAN","NQO2_HUMAN","UBC9_HUMAN", "ANT3_HUMAN","B2MG_HUMAN","BID_HUMAN",
    "SYUG_HUMAN","GSTA1_HUMAN","IGF2_HUMAN","PDGFB_HUMAN","THIO_HUMAN","CYC_HUMAN","GELS_HUMAN", "GSTP1_HUMAN","RASH_HUMAN",
    "IFNG_HUMAN","TRFE_HUMAN","UBE2C_HUMAN",
    "UB2E1_HUMAN","ANXA5_HUMAN","CRP_HUMAN","FABPH_HUMAN","IL8_HUMAN", "TRFL_HUMAN","TAU_HUMAN","SODC_HUMAN","TNFA_HUMAN")
  #ProteinName=NA
  )
  if(debug) message(fxNa," gUP1")
  UPS1$acOld <- UPS1$acNew <- sub("\\-[[:digit:]]+", "", UPS1$acFull)
  UPS1$acOld[8] <- UPS1$acFull[8] <- "P62988"   # the accession cited by Sigma-Aldrich
  UPS1$ac <- if(isFALSE(updated)) UPS1$acOld else UPS1$acNew
  UPS1 }
   
