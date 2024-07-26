#' Accession-Numbers And Names Of UPS1 Proteins
#'
#' UPS1 (see https://www.sigmaaldrich.com/FR/en/product/sigma/ups1) and UPS2 are commerical products consisting of a mix of 48 human (purified) proteins.
#' They are frequently used as standard in spike-in experiments, available from \href{https://www.sigmaaldrich.com}{Sigma-Aldrich}.
#' This function allows accessing their protein accession numbers and associated names on \href{https://www.uniprot.org/}{UniProt}
#'  
#' @return This function returns data.frame with accession-numbers as stated by the supplier (\code{$acFull}),
#'  trimmed accession-numbers, ie without version numbers (\code{$ac}) 
#'  and associated (\code{UniProt}) names on \href{https://www.uniprot.org/}{UniProt} as well as the species designation for the collection of 48 human UPS-1 proteins.
#' @examples
#' head(getUPS1acc())
#' @export
getUPS1acc <- function() {
  ## The accession numbers for the UPS1 proteins
  UPS1 <- data.frame( ac=rep(NA,48),
  acFull=c("P00915", "P00918", "P01031", "P69905", "P68871", "P41159", "P02768", "P62988",
    "P04040", "P00167", "P01133", "P02144", "P15559", "P62937", "Q06830", "P63165",
    "P00709", "P06732", "P12081", "P61626", "Q15843", "P02753", "P16083", "P63279",
    "P01008", "P61769", "P55957", "O76070", "P08263", "P01344", "P01127", "P10599",
    "P99999", "P06396", "P09211", "P01112", "P01579", "P02787", "O00762", "P51965",
    "P08758", "P02741", "P05413", "P10145", "P02788", "P10636-8", "P00441", "P01375"),
  uniProt=c("CAH1", "CAH2", "CO5", "HBA", "HBB", "LEP", "ALBU", "UBIQ", "CATA",
    "CYB5", "EGF", "MYG", "NQO1", "PPIA", "PRDX1", "SUMO1", "LALBA", "KCRM",
    "SYHC", "LYSC", "NEDD8", "RETBP", "NQO2", "UBC9", "ANT3", "B2MG", "BID",
    "SYUG", "GSTA1", "IGF2", "PDGFB", "THIO", "CYC", "GELS", "GSTP1", "RASH",
    "IFNG", "TRFE", "UBE2C", "UB2E1", "ANXA5", "CRP", "FABPH", "IL8", "TRFL",
    "TAU", "SODC", "TNFA"),
  species=rep("Homo sapiens", 48),
    name=NA)
  UPS1$ac <- sub("\\-[[:digit:]]+","", UPS1$acFull)
  UPS1 }
  
