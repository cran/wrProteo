#' Read protein annotation as exported from UniProt batch-conversion
#'
#' This function allows reading and importing protein-ID conversion results from \href{https://www.uniprot.org/uploadlists/}{UniProt}.
#' To do so, first copy/paste your query IDs into \href{https://www.uniprot.org/uploadlists/}{UniProt} 'Retrieve/ID mapping' field called '1. Provide your identifiers' (or upload as file), verify '2. Select options'.
#' In a typical case of 'enst000xxx' IDs  you may leave default settings, ie 'Ensemble Transcript' as input and 'UniProt KB' as output. Then, 'Submit' your search and retreive results via 
#' 'Download', you need to specify a 'Tab-separated' format ! If you download as 'Compressed' you need to decompress the .gz file before running the function \code{readUCSCtable} 
#' In addition, a file with UCSC annotation (Ensrnot accessions and chromosomic locations, obtained using \code{\link{readUCSCtable}}) can be integrated.
#' @details
#' In a typicall use case, first chromosomic location annotation is extracted from UCSC for the species of interest and imported to R using  \code{\link{readUCSCtable}} . 
#' However, the tables provided by UCSC don't contain Uniprot IDs. Thus, an additional (batch-)conversion step needs to get added. 
#' For this reason \code{\link{readUCSCtable}} allows writing a file with Ensemble transcript IDs which can be converted tu UniProt IDs at the site of  \href{https://www.uniprot.org/uploadlists/}{UniProt}. 
#' Then, UniProt annotation (downloaded as tab-separated) can be imported and combined with the genomic annotation using this function.  
#' @param UniProtFileNa (character) name (and path) of file exported from Uniprot (tabulated text file inlcuding headers) 
#' @param deUcsc (data.frame) object produced by \code{readUCSCtable} to be combined with data from \code{UniProtFileNa}
#' @param targRegion (character or list) optional marking of chromosomal locations to be part of a given chromosomal target region, 
#'   may be given as character like \code{chr11:1-135,086,622} or as \code{list} with a first component characterizing the chromosome and a integer-vector with start- and end- sites 
#' @param useUniPrCol (character) optional declaration which colums from UniProt exported file should be used/imported (default 'EnsID','Entry','Entry.name','Status','Protein.names','Gene.names','Length').
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of message(s) produced
#' @return data.frame (with columns $EnsID, $Entry, $Entry.name, $Status, $Protein.names, $Gene.names, $Length; if \code{deUcsc} is integrated plus: $chr, $type, $start, $end, $score, $strand, $Ensrnot, $avPos) 
#' @seealso \code{\link{readUCSCtable}}
#' @examples
#' path1 <- system.file("extdata",package="wrProteo")
#' deUniProtFi <- file.path(path1,"deUniProt_hg38chr11extr.tab")
#' deUniPr1a <- readUniProtExport(deUniProtFi) 
#' str(deUniPr1a)
#' 
#' ## Workflow starting with UCSC annotation (gtf) files :
#' gtfFi <- file.path(path1,"UCSC_hg38_chr11extr.gtf.gz")
#' UcscAnnot1 <- readUCSCtable(gtfFi)
#' ## Results of conversion at UniProt are already available (file "deUniProt_hg38chr11extr.tab")
#' myTargRegion <- list("chr1", pos=c(198110001,198570000))
#' myTargRegion2 <-"chr11:1-135,086,622"      # works equally well
#' deUniPr1 <- readUniProtExport(deUniProtFi,deUcsc=UcscAnnot1,
#'   targRegion=myTargRegion)
#' ## Now UniProt IDs and genomic locations are both available :
#' str(deUniPr1)
#' @export
readUniProtExport <- function(UniProtFileNa, deUcsc=NULL, targRegion=NULL, useUniPrCol=NULL, silent=FALSE, callFrom=NULL) {         
  ## read annotation exported from https://www.uniprot.org/uploadlists/  upload  Ensemble Transcript => UniprotKB => export 
  ## targRegion : list('chr1',pos=c(198110001,198570000)) or 'chr11:1-135,086,622'
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="readUniProtExport")
  if(!isTRUE(silent)) silent <- FALSE  
  if(length(UniProtFileNa) >1) UniProtFileNa <- UniProtFileNa[1] else {if(length(UniProtFileNa) < 1) stop(" argument 'UniProtFileNa' seems empty")}
  chFi <- file.exists(UniProtFileNa)
  if(!chFi) stop(" file '",UniProtFileNa,"' not found !")
  chExt <- length(grep("\\.gz$", UniProtFileNa, fixed=FALSE, perl=FALSE)) >0  
  chPa <- try(find.package("utils"),silent=TRUE)
  if("try-error" %in% class(chPa)) stop("package 'utils' not found ! Please install first")   
  ## main  
  deUniProt <- try(utils::read.delim(UniProtFileNa,stringsAsFactors=FALSE), silent=TRUE)
  errMsg1 <- " seems not to be in UniProt 'tab-separated' format (does not contain sufficent number of columns) !"
  if("try-error" %in% class(deUniProt)) { 
    deUniProt <- try(wrMisc::readVarColumns(if(chExt) unz(UniProtFileNa) else UniProtFileNa,callFrom=fxNa), silent=TRUE)
    if("try-error" %in% class(deUniProt)) stop("Can't read file '",UniProtFileNa,"' - please check format !") else {
      if(!silent) message(fxNa," Managed to read file using readVarColumns()") }
    if(ncol(deUniProt) <9) stop("file ",UniProtFileNa,errMsg1)  
    colnames(deUniProt)[1:9] <- c("EnsTraID","xx","UniprotID",colnames(deUniProt)[c(2:7)])  # initial colnames by readVarColumns are shifted    
  } 
  if(ncol(deUniProt) <7) stop("file ",UniProtFileNa,errMsg1)     # check if (in)sufficient numer of columns
  if(nrow(deUniProt) <2 & !silent) message(fxNa," CAUTION, file '",UniProtFileNa,"' contains only ",nrow(deUniProt)," lines !") 
  ## correct colnames
  chCol <- c(grep("yourlist.",colnames(deUniProt)[1]) >0, grep("isomap.",colnames(deUniProt)[2]) >0, "Entry" %in% colnames(deUniProt))
  if(chCol[1]) colnames(deUniProt)[1] <- "EnsTraID"
  if(chCol[2]) colnames(deUniProt)[2] <- "xx"                 # this column contains almost no information
  colnames(deUniProt)[3] <- "UniProtID"  
  ## combine with data initially/previously read from Ucsc
  multID <- NULL
  colnames(deUniProt) <- sub(" ",".",colnames(deUniProt))
  if(length(useUniPrCol) <1) useUniPrCol <- c("EnsTraID","UniProtID","Entry.name","Status","Protein.names","Gene.names","Length")
  useUniPrCo <- wrMisc::extrColsDeX(deUniProt, useUniPrCol, doExtractCols=FALSE, callFrom=fxNa, silent=silent)
  ##  treat multi-Ensemble entries : need to replicate lines of table for multiple concatenated (eg ENSRNOT00000031808,ENSRNOT00000093745)
  splitExtendConcat <- function(mat,useCol=1,sep=",",sep2="[[:digit:]],[[:alpha:]]+"){
    ## extend matrix or data.frame by additional lines if column 'useCol' contains multiple concatenated terms (content of other columns will be duplicated)
    ## 'sep' used with strsplit() and grep() to identify lines and split, also used to construct (generic) term for keeping just first 
    ## 'sep2' optional custom pattern used with grep() to identify lines; will be used instead of 'generic' sep to identify entries to split lateron 
    ## main
    chMult <- grep(if(length(sep2) >0) sep2 else sep, mat[,useCol], fixed=FALSE, perl=FALSE)
    if(length(chMult) >0)  {
      ##
      spl1 <- strsplit(mat[chMult,useCol],sep, fixed=FALSE, perl=FALSE)
      spl2 <- unlist(lapply(spl1, function(x) x[-1]), use.names=FALSE)
      toLine <- rep(chMult, sapply(spl1,length) -1)
      mat[,useCol] <- sub(paste(sep,"[[:print:]]*$",sep=""),"",mat[,useCol], fixed=FALSE, perl=FALSE)
      mat2 <-  cbind(spl2,mat[c(toLine),-1])
      colnames(mat2)[1] <- colnames(mat)[1]
      mat <- rbind(mat,mat2)
    }
    mat }
  deUniProt <- splitExtendConcat(deUniProt, sep=",", sep2="[[:digit:]],[[:upper:]]+")   
  if(length(deUcsc) >0) {
    chGeneId <- which(colnames(deUcsc) =="gene_id")
    if(length(chGeneId) <1) stop("Invalid file-content: The file '",UniProtFileNa,"' does not conatain a column 'gene_id' ! Please check the input file")
    deUcsc[,"gene_id"] <- sub("\\.[[:digit:]]+$","",deUcsc[,"gene_id"])
    useUcCol <- wrMisc::naOmit(match(c("gene_id","chr","start","end","strand","frame"),colnames(deUcsc)))
    deUcsc <- wrMisc::convMatr2df(deUcsc[,useUcCol], addIniNa=FALSE, callFrom=fxNa,silent=silent)    
    matchUniprInUcsc <- match(deUniProt[,1], deUcsc[,"gene_id"])
    if(sum(!is.na(matchUniprInUcsc)) <4) { 
      if(!silent) message(fxNa," low yield matching ",wrMisc::pasteC(deUniProt[1:3,1],quoteC="'")," and ",
        wrMisc::pasteC(deUcsc[1:3,"gene_id"],quoteC="'"), " convert all to lower case and remove version numbers ('xxx.2') for better matching")    
      matchUniprInUcsc <- match(sub("\\.[[:digit:]]+$","", tolower(deUniProt[,1])), sub("\\.[[:digit:]]+$","", tolower(deUcsc[,"gene_id"])))
      if(sum(!is.na(matchUniprInUcsc)) <4) warning(fxNa," Matching failed : Very few or no matches between UniProtFile and deUcsc !")} 
    if(!silent) message(fxNa," intergrating genomic information for ",length(matchUniprInUcsc)," entries (",sum(is.na(matchUniprInUcsc))," not found)")
    ## add chrom Loc to deUniProt => combined DB
    combAllChrDB <- cbind(deUniProt[,useUniPrCo], deUcsc[matchUniprInUcsc,])       ## add Ensrnot   c(1,3:5,7,10)
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
  chEnsID <- "gene_id" %in% colnames(combAllChrDB)
  if(chEnsID)  combAllChrDB <- combAllChrDB[,-1*which(colnames(combAllChrDB)=="gene_id")]
  combAllChrDB }   
  
