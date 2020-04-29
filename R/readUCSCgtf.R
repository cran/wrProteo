#' Read gtf annotation files from UCSC
#'
#' This function allos reading and importing genomic \href{https://genome.ucsc.edu/cgi-bin/hgTables}{gtf-annotation} data from UCSC.
#' In addition, a file with Ensrnot accessions only can be exported for further batch conversion on \href{https://www.uniprot.org/uploadlists/}{UniProt} by copy/pase to site,
#' since UCSC annotation data does not contain UniProt IDs. Subsequently, the resulting conversion by UniProt can be read and combined with UCSC annotation using \code{\link{readUniProtExport}}.
#' 
#' @param gtfFiName (character) name (and path) of file to read
#' @param exportFileNa (character) optional file-name to be exported, if \code{NULL} no file will be written
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of message(s) produced
#' @return data.frame (with columns  $chr, $source, $type, $start, $end, $score, $strand, $frame, $features, $Ensrnot) 
#' @seealso \code{\link{readUniProtExport}}, \code{\link{readPDExport}}, \code{\link{readMaxQuantFile}},
#' @examples
#' path1 <- system.file("extdata",package="wrProteo")
#' gtfFi <- file.path(path1,"UCSC_hg38_chr11extr.gtf")
#' # here we'll write the file for UniProt conversion to tempdir to keep things tidy
#' expFi <- file.path(tempdir(),"deUcscForUniProt2.txt")
#' UcscAnnot1 <- readUCSCgtf(gtfFi,exportFileNa=expFi)
#' ## results can be further combined with readUniProtExport() 
#' deUniProtFi <- file.path(path1,"deUniProt_hg38chr11extr.tab")
#' deUniPr1 <- readUniProtExport(deUniProtFi,deUcsc=UcscAnnot1,
#'   targRegion="chr11:1-135,086,622")  
#' deUniPr1[1:5,-5] 
#' @export
readUCSCgtf <- function(gtfFiName,exportFileNa="deUcsc2UniProt.dat",silent=FALSE,callFrom=NULL) {
  ## read & parse ensGene.gtf type file from UCSC, (optional) export to file for batch conversion on UniProt, return annotation (matrix)
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="readUCSCgtf")
  if(length(gtfFiName) >1) gtfFiName <- gtfFiName[1] else {if(length(gtfFiName) < 1) stop(" argument 'gtfFiName' seems empty")}
  chFi <- file.exists(gtfFiName)
  if(!chFi) stop(" file '",gtfFiName,"' not found !")
  chGz <- length(grep("\\.gz$",gtfFiName)) >0
  if(chGz & !silent) message(fxNa," Beware, reading of compressed files may pose problems" ) 
  ## main
  ensG1 <- try(utils::read.delim(if(chGz) unz(gtfFiName) else gtfFiName,header=FALSE,stringsAsFactors=FALSE))
  if("try-error" %in% class(chFi)) stop("Can't read file '",gtfFiName,"' - please check format")
  if(ncol(ensG1) != 9) stop("file seems not to be UCSC gtf format (does not contain sufficent number of columns) !")  
  colnames(ensG1) <- c("chr","source","type","start","end","score","strand","frame","features")
  if(!silent) message(" read ",nrow(ensG1)," lines and ",ncol(ensG1)," cols")
  EnsrnotAllChr <- tolower(sub("\\.[[:digit:]]+[[:print:]]+","",sub("^[[:print:]]*transcript_id ","",ensG1[,9])))
  ensG1 <- cbind(ensG1,Ensrnot=EnsrnotAllChr)
  isDupEnsAllChr <- duplicated(EnsrnotAllChr)
  ## better make unique now in respect to transcript_id
  ensG1 <- ensG1[which(!isDupEnsAllChr),]
  EnsrnotAllChr <- EnsrnotAllChr[which(!isDupEnsAllChr)]
  ## write to file  (for batch search/vonversion on UniProt site)
  if(length(exportFileNa) >0) { cat(unique(EnsrnotAllChr), file=exportFileNa[1], sep="\n") 
    if(!silent) message(fxNa," Exporting file  '",exportFileNa,"'  for conversion on https://www.uniprot.org/uploadlists") }                 # export to file for batch conversion on UniProt
  ensG1 }
   
