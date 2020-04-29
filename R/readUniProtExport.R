#' Read export from UniProt batch-conversion
#'
#' This function allows reading and importing protein-ID conversion results from \href{https://www.uniprot.org/uploadlists/}{UniProt}.
#' To do so, first copy/paste your query IDs into UniProt 'Retrieve/ID mapping' field called '1. Provide your identifiers' (or upload as file), verify '2. Select options'.
#' In a typical case of 'enst000xxx' IDs  you may leave default settings, ie 'Ensemble Transcript' as input and 'UniProt KB' as output. Then, 'Submit' your search and retreive results via 
#' 'Download', you need to specify a 'Tab-separated' format ! If you download as 'Compressed' you need to decompress the .gz file before running the function \code{readUCSCgtf} 
#' In addition, a file with UCSC annotation (Ensrnot accessions and chromosomic locations, obtained using \code{\link{readUCSCgtf}}) can be integrated.
#' @details
#' In a typicall use case, first chromosomic loacation annotation is extracted from UCSC for the species of interest and imported to R using  \code{\link{readUCSCgtf}} . 
#' However, the tables provided by UCSC don't contain Uniprot IDs. Thus, an additional (batch-)conversion step needs to get added. 
#' For this reason \code{\link{readUCSCgtf}} allows writing a file with Ensemble transcript IDs which can be converted tu UniProt IDs at the site of  \href{https://www.uniprot.org/uploadlists/}{UniProt}. 
#' Then, UniProt annotation (downloaded as tab-separated) can be imported and combined with the genomic annotation using this function.  
#' @param UniProtFileNa (character) name (and path) of file exported from Uniprot (tabulated text file inlcuding headers) 
#' @param deUcsc (data.frame) object produced by \code{readUCSCgtf} to be combined with data from \code{UniProtFileNa}
#' @param targRegion (character or list) optional marking of chromosomal locations to be part of a given chromosomal target region, may be given as character like \code{chr11:1-135,086,622} or as \code{list} with a firts component chracterizing the chromosome and a integer-vector with start- and end- sites 
#' @param useUniPrCol (character) optional declaration which colums from UniProt exported file should be used/imported (default 'EnsID','Entry','Entry.name','Status','Protein.names','Gene.names','Length').
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of message(s) produced
#' @return data.frame (with columns $EnsID, $Entry, $Entry.name, $Status, $Protein.names, $Gene.names, $Length; if \code{deUcsc} is integrated plus: $chr, $type, $start, $end, $score, $strand, $Ensrnot, $avPos) 
#' @seealso \code{\link{readUCSCgtf}}
#' @examples
#' path1 <- system.file("extdata",package="wrProteo")
#' deUniProtFi <- file.path(path1,"deUniProt_hg38chr11extr.tab")
#' deUniPr1a <- readUniProtExport(deUniProtFi) 
#' ## with including chromosomic location extracted via Ucsc gtf files
#' gtfFi <- file.path(path1,"UCSC_hg38_chr11extr.gtf")
#' ## Here we won't write the file for UniProt since the the results of the
#' ##   conversion at Uniprot are already vailable as file "deUniProt_hg38chr11extr.tab"
#' UcscAnnot1 <- readUCSCgtf(gtfFi,exportFileNa=NULL)
#' deUniPr1 <- readUniProtExport(deUniProtFi,deUcsc=UcscAnnot1,
#'   targRegion="chr11:1-135,086,622")  
#' @export
readUniProtExport <- function(UniProtFileNa,deUcsc=NULL,targRegion=NULL,useUniPrCol=NULL,silent=FALSE,callFrom=NULL) {         
  ## read annotation exported from https://www.uniprot.org/uploadlists/  upload  Ensemble Transcript => UniprotKB => export 
  ## targRegion : list('chr1',pos=c(198110001,198570000)) or 'chr11:1-135,086,622'
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="readUniProtExport")
  if(length(UniProtFileNa) >1) UniProtFileNa <- UniProtFileNa[1] else {if(length(UniProtFileNa) < 1) stop(" argument 'UniProtFileNa' seems empty")}
  chFi <- file.exists(UniProtFileNa)
  if(!chFi) stop(" file '",UniProtFileNa,"' not found !")
  chExt <- length(grep("\\.gz$",UniProtFileNa)) >0
  if(chExt & !silent) message(fxNa," Beware, reading of compressed files may pose problems" ) 
  ## main  
  deUniProt <- try(utils::read.delim(if(chExt) unz(UniProtFileNa) else UniProtFileNa,stringsAsFactors=FALSE))
  if("try-error" %in% class(deUniProt)) stop("Can't read file '",UniProtFileNa,"' - please check format !")
  if(ncol(deUniProt) <7) stop("file seems not to be in UniProt 'tab-separated' format (does not contain sufficent number of columns) !")  
  colnames(deUniProt)[1:2] <- c("EnsID",sapply(colnames(deUniProt)[2],function(x) sub("\\.[[:upper:]][[:alnum:]]+","",x)))
  if(nrow(deUniProt) <2 & !silent) message(fxNa," CAUTION, file '",UniProtFileNa,"' contains only ",nrow(deUniProt)," lines !") 
  ## combine with data initially/previously read from Ucsc
  multID <- NULL
  if(length(useUniPrCol) <1) useUniPrCol <- c("EnsID","Entry","Entry.name","Status","Protein.names","Gene.names","Length")
  useUniPrCo <- wrMisc::extrColsDeX(deUniProt,useUniPrCol,doExtractCols=FALSE,callFrom=fxNa,silent=silent)
  if(length(deUcsc) >0) {
    matchUniprInUcsc <- match(deUniProt[,1], deUcsc$Ensrnot)
    if(!silent) message(fxNa," intergrating genomic information for ",length(matchUniprInUcsc)," entries (",sum(is.na(matchUniprInUcsc))," not found)")
    ## add chrom Loc to deUniProt => combined DB
    combAllChrDB <- cbind(deUniProt[,useUniPrCo], deUcsc[matchUniprInUcsc,c(1,3:5,7,10)])       ## add Ensrnot
    if(!silent) message(fxNa," ",nrow(combAllChrDB)," IDs in output")
    combAllChrDB <- cbind(combAllChrDB,avPos=if(all(c("start","end") %in% colnames(combAllChrDB))) {
      round(rowMeans(combAllChrDB[,c("start","end")])) } else NA)    # add mean gene-position for easier sorting
    ## mark if genimic positions in targer region
    if(!all(c("chr","start") %in% colnames(combAllChrDB))) targRegion <- NULL
    if(length(targRegion) >0) if(is.character(targRegion) & length(targRegion)==1) {
      targRegion <- unlist(strsplit(targRegion,":"))
      targRegion <- list(targRegion[1],gsub(",","",unlist(strsplit(targRegion[2],"-")))) }
    combAllChrDB <- cbind(combAllChrDB,inTarg=if(length(targRegion) >0) {
      combAllChrDB[,"chr"]==targRegion[[1]] & as.integer(combAllChrDB[,"start"]) >targRegion[[2]][1] & as.integer(combAllChrDB[,"end"]) <targRegion[[2]][2]} else NA)
  } else combAllChrDB <- deUniProt[,useUniPrCo]
  ## convert factor-columns to character
  chFa <- rep(NA,ncol(combAllChrDB))
  for(i in 1:ncol(combAllChrDB)) chFa[i] <- is.factor(combAllChrDB[,i])
  if(any(chFa)) for(i in which(chFa)) combAllChrDB[,i] <- as.character(combAllChrDB[,i])
  combAllChrDB }    
  
