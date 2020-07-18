#' Read annotation files from UCSC
#'
#' This function allows reading and importing genomic \href{https://genome.ucsc.edu/cgi-bin/hgTables}{UCSC-annotation} data.
#' Files can be read as default UCSC exprt or as GTF-format. 
#' In the context of proteomics we noticed that sometimes UniProt tables from UCSC are hard to match to identifiers from UniProt Fasta-files, ie many protein-identifiers won't match.
#' For this reason additional support is given to reading 'Genes and Gene Predictions': Since this table does not include protein-identifiers, a non-redundant list of ENSxxx transcript identifiers 
#' can be exprted as file for an additional stop of conversion, eg using a batch conversion tool at the site of \href{https://www.uniprot.org/uploadlists/}{UniProt}. 
#' The initial genomic annotation can then be complemented using \code{\link{readUniProtExport}}. 
#' Using this more elaborate route, we found higher coverage when trying to add genomic annotation to protein-identifiers to proteomics results with annnotation based on an initial Fasta-file. 
#' 
#' @param fiName (character) name (and path) of file to read
#' @param exportFileNa (character) optional file-name to be exported, if \code{NULL} no file will be written
#' @param gtf (logical) specify if file \code{fiName} in gtf-format (see \href{https://genome.ucsc.edu/cgi-bin/hgTables}{UCSC})
#' @param simplifyCols (character) optional list of column-names to be used for simplification  (if 6 column-headers are given) : the 1st value will be used to identify the column  
#'  used as refence to summarize all lines with this ID; for the 2nd (typically chromosome names) will be taken a representative value, 
#'  for the 3rd (typically gene start site) will be taken the minimum, 
#'  for the 4th (typically gene end site) will be taken the maximum, for the 5th and 6th a representative values will be reported;  
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of message(s) produced
#' @return matrix, optionally the file 'exportFileNa' may be written
#' @seealso \code{\link{readUniProtExport}}, \code{\link{readPDExport}}, \code{\link{readMaxQuantFile}},
#' @examples
#' path1 <- system.file("extdata",package="wrProteo")
#' gtfFi <- file.path(path1,"UCSC_hg38_chr11extr.gtf")
#' # here we'll write the file for UniProt conversion to tempdir() to keep things tidy
#' expFi <- file.path(tempdir(),"deUcscForUniProt2.txt")
#' UcscAnnot1 <- readUCSCtable(gtfFi,exportFileNa=expFi)
#' 
#' ## results can be further combined with readUniProtExport() 
#' deUniProtFi <- file.path(path1,"deUniProt_hg38chr11extr.tab")
#' deUniPr1 <- readUniProtExport(deUniProtFi,deUcsc=UcscAnnot1,
#'   targRegion="chr11:1-135,086,622")  
#' deUniPr1[1:5,-5] 
#' @export
readUCSCtable <- function(fiName, exportFileNa=NULL, gtf=NA, simplifyCols=c("gene_id","chr","start","end","strand","frame"), silent=FALSE, callFrom=NULL) {
  ## read & parse ensGene.gtf type file from UCSC, (optional) export to file for batch conversion on UniProt, return annotation (matrix)
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="readUCSCtable")
  if(length(fiName) >1) fiName <- fiName[1] else {if(length(fiName) < 1) stop(" argument 'fiName' seems empty")}
  chFi <- file.exists(fiName)
  if(!chFi) stop(" file '",fiName,"' not found ! (maybe you are not pointing to the correct direcory ?)")
  chPa1 <- try(find.package("utils"), silent=TRUE)
  chPa2 <- try(find.package("R.utils"), silent=TRUE)
  if("try-error" %in% class(chPa1)) stop("package 'utils' not found ! Please install first")     
  if("try-error" %in% class(chPa2)) stop("package 'R.utils' not found ! Please install first")     
  ## check/find out if file is gtf
  chGz <-  R.utils::isGzipped(fiName)  
  if(is.na(gtf)) {                           # try to guess/check if gtf=TRUE
    gtfColNa <- c("#bin	name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2",
      "cdsStartStat","cdsEndStat","exonFrame","#chrom","chromStart","chromEnd","score","thickStart","thickEnd","match","misMatch","repMatch",
      "nCount","acc","status","ensGene","ensTrans")
    tmp <- try(utils::read.delim(fiName, header=FALSE, stringsAsFactors=FALSE, nrows=7))
    if("try-error" %in% class(tmp)) {
      gtf <- length(grep("\\.gtf$|\\.gtf\\.gz$", tolower(basename(fiName)))) >0
      if(!silent) message(fxNa," Quick guess if compressed file '",basename(fiName),"' is gtf : ",gtf)
    } else {
      gtf <- sum(gtfColNa %in% tmp[1,],na.rm=TRUE) <7
      if(!silent) message(fxNa," File '",basename(fiName),"' is gtf : ",gtf)
    }
  }
  ## main reading
  ensG1 <- try(utils::read.delim(fiName, header=!gtf, stringsAsFactors=FALSE))
  if("try-error" %in% class(chFi)) {
    ## try other format ?
    errMsg <- "' - please check format or see if file is readable"
    if(chGz) {                                     # need to decompress first : copy to tempdir, uncompress, read      
      file.copy(fiName, file.path(tempdir(),basename(fiName)))
      ## decompress      
      tmp <- try(R.utils::gunzip(file.path(tempdir(),basename(fiName))),silent=TRUE)
      if("try-error" %in% class(tmp)) stop(" Failed to decompress file ",fiName)
      ## read 
      ensG1 <- try(wrMisc::readVarColumns(file.path(tempdir(),basename(fiName),header=!gtf,silent=silent,callFrom=fxNa), silent=TRUE,callFrom=fxNa))
      if("try-error" %in% class(ensG1)) stop("Can't read file '",fiName,errMsg)
      rmFi <- file.path(tempdir(), c(basename(fiName), sub("\\.gz$","",basename(fiName))))
      sapply(rmFi,function(x) if(file.exists(x)) file.remove(x))              # clean up files
    } else ensG1 <- try(wrMisc::readVarColumns(fiName, header=!gtf, silent=silent, callFrom=fxNa), silent=TRUE)
    if("try-error" %in% class(ensG1)) stop(" Could not succed to read file '",basename(fiName),errMsg)  
  }
  ## correct colnames, export ENSG if available  
  if(gtf) {
    ##  gtf-format, add standard colnames
    colnames(ensG1) <- c("chr","source","type","start","end","score","strand","frame","features")[1:ncol(ensG1)]
    ## now split last column to gene_id and transcript_id
    chTId <- grep("; transcript_id [[:alpha:]]+[[:digit:]]+" , as.character(ensG1[1:30,9]))
    chGId <- grep("^gene_id [[:alpha:]]+[[:digit:]]+", as.character(ensG1[1:30,9]))
    chId <- grep("transcript_id [[:punct:]]{0,1}ENS[[:upper:]]+[[:digit:]][[:print:]]+", as.character(ensG1[1:30,9]))
    if(length(chTId) >0 & length(chGId) >0) {
      ## need to split column
      newC <- matrix(unlist(strsplit(sub("^gene_id ","",as.character(ensG1[,9])),"; transcript_id ")),
        byrow=TRUE, ncol=2, dimnames=list(NULL,c("gene_id","transcript_id")))
      newC <- sub("; $","",newC)  
      ensG1 <- cbind(ensG1[,-9], newC)
      ## export (gtf case) 
    }
  }  
  ## option to summarize by first column of 'simplifyCols'

  if(length(simplifyCols) >5) { simplifyCols <- simplifyCols[which(simplifyCols %in% colnames(ensG1))]
    if(length(simplifyCols) <6 & !silent) message(fxNa," Cannot find sufficient column-names given in argument 'simplifyCols', ignoring ..") }  
  if(length(simplifyCols) >5) { 
    iniDim <- dim(ensG1)
    ensG1 <- matrix(unlist(by(ensG1[,simplifyCols[1:6]], ensG1[,simplifyCols[1]], 
      function(x) { if(length(dim(x)) <2) x <- matrix(x, ncol=6, dimnames=list(NULL,simplifyCols[1:6]))
      c(x[1,1],x[1,2], min(x[,3],na.rm=TRUE), max(x[,4],na.rm=TRUE), x[1,5], x[1,6])} )),
      byrow=TRUE, ncol=6, dimnames=list(NULL,simplifyCols))
    if(!silent) message(fxNa," simplifed from ",iniDim[1]," to ",nrow(ensG1)," non-redundant ",simplifyCols[1])  
    ## export ENSRNOT for onversion at UniProt site
    if(length(exportFileNa) >0) {         # also need to remove ENST-version-tags (since UniProt won't recognize Ensemble gene IDs with version tags)
      exportFileNa <- gsub("\\\\","/",exportFileNa)
      forFile <- unique(sub("\\.[[:digit:]]+$","",ensG1[,"gene_id"]))
      msg <- c("'  for conversion on https://www.uniprot.org/uploadlists")
      if(!silent) message(fxNa," Write to file : ",paste(utils::head(forFile,4),collapse=", ")," ...")
      if(!silent) if(file.exists(exportFileNa[1])) message(fxNa," Beware, file '",exportFileNa[1],"' will be overwritten !") else  message(fxNa,
        " Exporting file  '",exportFileNa,msg[1])          # export to file for batch conversion on UniProt
      tmp <- try(cat(forFile, file=exportFileNa[1], sep="\n"), silent=TRUE)
      if("try-error" %in% class(tmp)) warning(fxNa," Beware: Did not succed to write results to file")}       
      }
  ##
  ensG1 }  
    
