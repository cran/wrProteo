#' Read Tabulated Files Exported by ProteomeDiscoverer At Peptide Level
#'
#' Initials petide identificationa and quantification results form \href{https://www.thermofisher.com/order/catalog/product/OPTON-30812}{Thermo ProteomeDiscoverer}
#' which were exported as tabulated text can be imported and relevant information extracted.
#' The final output is a list containing 3 elements: \code{$annot}, \code{$raw} and optional \code{$quant}, or returns data.frame with entire content of file if \code{separateAnnot=FALSE}.
#'
#' @details
#' This function has been developed using Thermo ProteomeDiscoverer versions 2.2 to 2.5.
#' The format of resulting files at export also depends which columns are chosen as visible inside ProteomeDiscoverer and subsequently get chosen for export.
#' Using the argument \code{suplAnnotFile} it is possible to specify a specific file (or search for default file) to read for extracting file-names as sample-names and other experiment realted information.
#' Precedent and following aminoacids (relative to identified protease recognition sites) will be removed form peptide sequences and be displayed in $annot as columns 'prec' and 'foll'.
#' If a column named \code{contamCol} is found, the data will be lateron filtered to remove all contaminants, set to \code{NULL} for keeping all contaminants
#' This function replaces the depreciated function \code{readPDExport}.
#'
#'  Besides, ProteomeDiscoverer version number and full raw-file path will be extracted for $notes in final output.
#'
#' @param fileName (character) name of file to be read
#' @param path (character) path of file to be read
#' @param normalizeMeth (character) normalization method, defaults to \code{median}, for more details see \code{\link[wrMisc]{normalizeThis}})
#' @param sampleNames (character) new column-names for quantification data (ProteomeDiscoverer does not automatically use file-names from spectra); this argument has priority over \code{suplAnnotFile}
##' @param sdrf (character, list or data.frame) optional extraction and adding of experimenal meta-data: if character, this may be the ID at ProteomeExchange. Besides, the output from \code{readSdrf} or a list from \code{defineSamples} may be provided; if \code{gr} is provided, it gets priority for grouping of replicates
#' @param suplAnnotFile (logical or character) optional reading of supplemental files produced by MaxQuant; if \code{gr} is provided, it gets priority for grouping of replicates
#'  if \code{TRUE} default to files 'summary.txt' (needed to match information of \code{sdrf}) and 'parameters.txt' which can be found in the same folder as the main quantitation results;
#'  if \code{character} the respective file-names (relative ro absolute path), 1st is expected to correspond to 'summary.txt' (tabulated text, the samples as given to MaxQuant) and 2nd to 'parameters.txt' (tabulated text, all parameters given to MaxQuant)
#' @param read0asNA (logical) decide if initial quntifications at 0 should be transformed to NA
#' @param quantCol (character or integer) exact col-names, or if length=1 content of \code{quantCol} will be used as pattern to search among column-names for $quant using \code{grep}
#' @param contamCol (character or integer, length=1) which columns should be used for contaminants marked by ProteomeDiscoverer.
#'  If a column named \code{contamCol} is found, the data will be lateron filtered to remove all contaminants, set to \code{NULL} for keeping all contaminants
#' @param refLi (character or integer) custom specify which line of data is main species, if character (eg 'mainSpe'), the column 'SpecType' in $annot will be searched for exact match of the (single) term given
#' @param separateAnnot (logical) if \code{TRUE} output will be organized as list with \code{$annot}, \code{$abund} for initial/raw abundance values and \code{$quant} with final normalized quantitations
#' @param annotCol (character) column names to be read/extracted for the annotation section (default  c("Accession","Description","Gene","Contaminant","Sum.PEP.Score","Coverage....","X..Peptides","X..PSMs","X..Unique.Peptides", "X..AAs","MW..kDa.") )
#' @param FDRCol (list) optional indication to search for protein FDR information
#' @param tit (character) custom title to plot
#' @param graphTit (character) depreciated custom title to plot, please use 'tit'
#' @param wex (integer) relative expansion factor of the violin-plot (will be passed to \code{\link[wrGraph]{vioplotW}})
#' @param specPref (character or list) define characteristic text for recognizing (main) groups of species (1st for comtaminants - will be marked as 'conta', 2nd for main species- marked as 'mainSpe',
#'  and optional following ones for supplemental tags/species - maked as 'species2','species3',...);
#'  if list and list-element has multiple values they will be used for exact matching of accessions (ie 2nd of argument \code{annotCol})
#' @param plotGraph (logical or integer) optional plot of type vioplot of initial and normalized data (using \code{normalizeMeth}); if integer, it will be passed to \code{layout} when plotting
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a list with \code{$raw} (initial/raw abundance values), \code{$quant} with final normalized quantitations, \code{$annot}, \code{$counts} an array with number of peptides, \code{$quantNotes}
#'  and \code{$notes}; or if \code{separateAnnot=FALSE} the function returns a data.frame with annotation and quantitation only
#' @seealso \code{\link[utils]{read.table}}, \code{\link[wrMisc]{normalizeThis}}) , \code{\link{readMaxQuantFile}}, \code{\link{readProlineFile}}
#' @examples
#' path1 <- system.file("extdata", package="wrProteo")
#'
#' @export
readProtDiscovPeptides <- function(fileName, path=NULL, normalizeMeth="median", sampleNames=NULL, suplAnnotFile=TRUE, sdrf=NULL,read0asNA=TRUE, quantCol="^Abundances*",
  annotCol=NULL, contamCol="Contaminant", refLi=NULL, separateAnnot=TRUE, FDRCol=list(c("^Protein.FDR.Confidence","High"), c("^Found.in.Sample.","High")), plotGraph=TRUE, tit="Proteome Discoverer", graphTit=NULL, wex=1.6,
  specPref=c(conta="CON_|LYSC_CHICK", mainSpecies="OS=Homo sapiens"), silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## read ProteomeDiscoverer exported txt
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readProtDiscovPeptides")
  oparMar <- if(plotGraph) graphics::par("mar") else NULL       # only if figure might be drawn

  reqPa <- c("utils","wrMisc")
  chPa <- sapply(reqPa, requireNamespace, quietly=TRUE)
  if(any(!chPa)) stop("package(s) '",paste(reqPa[which(!chPa)], collapse="','"),"' not found ! Please install first from CRAN")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  excluCol <- "^Abundances.Count"   # exclude this from quantifications columns  # needed ?
  cleanDescription <- TRUE          # clean 'Description' for artifacts of truncated text (tailing ';' etc)
  modifSensible <- TRUE             # separate modified from unmodified peptides (by attaching modif to seq)
  setupSd <- quant <-NULL           # initialize

  ## check if path & file exist
  msg <- "Invalid entry for 'fileName'"
  if(length(fileName) >1) { fileName <- fileName[1]
    if(!silent) message(fxNa," 'fileName' shoud be of length=1, using 1st value")
  } else { if(length(fileName) <1) stop(msg) else if(nchar(fileName) <0) stop(msg)}
  paFi <- fileName                      # presume (& correct if path is given)
  if(debug) {message(fxNa,"rPDP1"); rPDP1 <- list(fileName=fileName,path=path,normalizeMeth=normalizeMeth,sampleNames=sampleNames,suplAnnotFile=suplAnnotFile,read0asNA=read0asNA,quantCol=quantCol,cleanDescription=cleanDescription)  }

  if(length(path) >0) if(!dir.exists(path[1])) { path <- NULL
    if(!silent) message(fxNa,"Invalid path '",path[1],"'  (not existing), ignoring...") }
  if(length(path) >0) { chFi <- file.exists(file.path(path[1], fileName))
    if(chFi) paFi <- file.path(path[1], fileName) else {
      if(file.exists(fileName)) {paFi <- fileName
        if(!silent) message(fxNa,"NOTE : Unable to find file '",fileName,"' in path '",path,"' but found in current path !")
      } else chFi <- FALSE                      # if path+fileName not found, check without path
    }
  } else chFi <- file.exists(fileName)         #
  if(!chFi) stop(" File ",fileName," was NOT found ",if(length(path) >0) paste(" in path ",path)," !")
  if(!grepl("\\.txt$|\\.txt\\.gz$", fileName)) message(fxNa,"Trouble ahead, expecting tabulated text file ('",fileName,"' might not be in right format) !!")
  if(debug) message(fxNa,"rPDP2a ..")
  ## prepare for reading files
  if(debug) { message(fxNa,"rPDP3 .. Ready to read", if(length(path) >0) c(" from path ",path[1])," the file  ",fileName[1])
    }

  ## read (main) file
  ## future: look for fast reading of files
  tmp <- try(utils::read.delim(file.path(paFi), stringsAsFactors=FALSE, header=TRUE), silent=TRUE)
  if(inherits(tmp, "try-error")) stop("Unable to read input file ('",paFi,"')!")
  if(debug) { message(fxNa,"rPDP3b .. dims of initial data : ", nrow(tmp)," li and ",ncol(tmp)," col ")}
  ## extract peptide sequence
  pepSe <- sub("\\.\\[[[:upper:]]\\]$","", sub("^\\[[[:upper:]]\\]\\.","", tmp[,"Annotated.Sequence"]))
  precAA <- postAA <- rep("",nrow(tmp))
  ch1 <- grep("\\.\\[[[:upper:]]\\]$", tmp[,"Annotated.Sequence"])
  if(length(ch1) >0) precAA[ch1] <- substr(tmp[ch1,"Annotated.Sequence"], 2, 2)
  ch1 <- grep("^\\[[[:upper:]]\\]\\.", tmp[,"Annotated.Sequence"])
  if(length(ch1) >0) postAA[ch1] <- substr(sub(".*\\.\\[","",tmp[ch1,"Annotated.Sequence"]), 1, 1)
  ## other peptide/protein info
  #txtCol <- c("Modifications", "Master.Protein.Accessions","Positions.in.Master.Proteins","Master.Protein.Descriptions")
  seqCol <- c("Sequence","Annotated.Sequence","Modifications", "Qvality.PEP","Qvality.q.value", "Number.of.Protein.Groups","Number.of.Proteins","Number.of.PSMs",  # no1-8
    "Master.Protein.Accessions","Positions.in.Master.Proteins","Modifications.in.Master.Proteins","Master.Protein.Descriptions",  #no 9-12;  last (ie 12th) missing in data from LN
    "Number.of.Missed.Cleavages","Theo.MHplus.in.Da","Contaminant",                                # no 13-14
    "Charge.by.Search.Engine.A5.Sequest.HT","XCorr.by.Search.Engine.A10.Sequest.HT","XCorr.by.Search.Engine.A5.Sequest.HT", "Top.Apex.RT.in.min" ) # no 15-18; 15 & 16 are currently not used, but use grep for 'Charge'
  maSeCo1 <- match(seqCol, colnames(tmp))
  maSeCo2 <- match(gsub("",".",seqCol), colnames(tmp))
  maSeCo <- if(sum(is.na(maSeCo1)) > sum(is.na(maSeCo2))) maSeCo2 else maSeCo1   # switch betw R-friendly and std
     #quanCo <- "Abundance.F62.Sample.na"
  AbundCol <- "^Abundance"              # use as pattern
  IdentTyCol <- "Found.in.Sample"       # use as pattern
  ## need other example for extracting quantifications ?
    #"Confidence.by.Search.Engine.Sequest.HT","Percolator.q.Value.by.Search.Engine.Sequest.HT","Percolator.PEP.by.Search.Engine.Sequest.HT", "XCorr.by.Search.Engine.Sequest.HT","Channel.Occupancy.in.Percent")
  if(debug) {message(fxNa,"rPDP4 .. Ready to read", if(length(path) >0) c(" from path ",path[1])," the file  ",fileName[1])
     rPDP4 <- list(fileName=fileName,path=path,chFi=chFi,paFi=paFi,normalizeMeth=normalizeMeth,sampleNames=sampleNames,suplAnnotFile=suplAnnotFile,read0asNA=read0asNA,quantCol=quantCol,cleanDescription=cleanDescription,tmp=tmp,seqCol=seqCol,maSeCo=maSeCo,modifSensible=modifSensible)}
  .chColNa <- function(x, mat, renameTo=NULL, silent=FALSE, fxNa=NULL){
    ## check in 'matr' for column-name 'x', if required rename best hit (if no direct hit look using grep, then grep wo case); return corrected mat
    chX <- x %in% colnames(mat)
    if(all(chX)) {
      if(is.character(renameTo) & length(renameTo) ==1) colnames(mat)[match(x, colnames(mat))] <- renameTo   # juste simple rename
    } else {     # try to localize column to use
      chX <- grep(x, colnames(mat))
      if(length(chX) >0) {
        if(is.character(renameTo) & length(renameTo) ==1) colnames(mat)[chX[1]] <- renameTo else x
        if(!silent & length(chX) >1) message(fxNa,"Found multiple columns containing '",x,"' : ",wrMisc::pasteC(colnames(mat)[chX], quoteC="'"),", using 1st")
      } else {
        chX <- grep(tolower(x), tolower(colnames(mat)))
        if(length(chX) >0) {
          if(is.character(renameTo) & length(renameTo) ==1) colnames(mat)[chX[1]] <- renameTo else x
          if(!silent & length(chX) >1) message(fxNa,"Found multiple columns containing '",tolower(x),"' : ",wrMisc::pasteC(colnames(mat)[chX], quoteC="'"),", using 1st")
        } else stop("Could NOT find column '",x,"' !!\n  (available columns ",wrMisc::pasteC(colnames(mat), quoteC="'"),")") }
    }
  mat }
  ## EXTRACT PEPTIDE SEQUENCES
  ## extract peptide sequences
  if(debug) {message(fxNa,"rPDP4a .. Ready to start extracting pep seq ")
     rPDP4a <- list(fileName=fileName,path=path,chFi=chFi,paFi=paFi,normalizeMeth=normalizeMeth,sampleNames=sampleNames,suplAnnotFile=suplAnnotFile,read0asNA=read0asNA,quantCol=quantCol,cleanDescription=cleanDescription,tmp=tmp,seqCol=seqCol,maSeCo=maSeCo,modifSensible=modifSensible)}
  if(is.na(maSeCo[1])) { if(is.na(maSeCo[2])) {if(!silent) message(fxNa,"Invalid type of data"); pepSeq <- NULL
    } else {
      pepSeq <- tmp[,maSeCo[2]] } #sub("\\.\\[A-Z\\]$", "", sub("^\\[A-Z\\]\\.", "", tmp[,maSeCo[2]])) }
  } else pepSeq <- tmp[,maSeCo[1]]
  fxPrecAA <- function(x) {   ## separate/extract note of preceeding & following AA; take char vector, returns 3-column matrix
    chPre <- grep("^\\[([[:upper:]]|\\-)\\]\\.", x)         # has note of preceeding AA
    chFoll <- grep("\\.\\[([[:upper:]]|\\-)\\]($|_)", x)    # has note of following AA
    out <- cbind(pep=sub("\\.\\[([[:upper:]]|\\-)\\]","", sub("^\\[([[:upper:]]|\\-)\\]\\.","", x)), prec=NA, foll=NA)
    if(length(chPre) >0) out[chPre,2] <- sub(".*\\[","", sub("\\]\\..+","", x[chPre]))   # the preceeding AA
    if(length(chFoll) >0) out[chFoll,3] <- sub("\\].*","", sub(".+\\.\\[","", x[chFoll]))
    out }
  annot1 <- fxPrecAA(pepSeq)
  pepSeq <- annot1[,1]

  if(modifSensible) { hasMod <- nchar(tmp[,maSeCo[3]]) >0
    if(any(hasMod, na.rm=TRUE)) pepSeq[which(hasMod)] <- paste(pepSeq[which(hasMod)],tmp[which(hasMod),maSeCo[3]],sep="_") }     # modificaton-separator
  if(debug) {message(fxNa,"rPDP4b .. Done extracting pep seq ") }
  ## ANNOATION (protein oriented)
   usColAnn <- maSeCo[c(3,6:7,9:14)]

  if(any(is.na(usColAnn), na.rm=TRUE)) {
    if(!silent) message(fxNa,"Note : ",sum(is.na(usColAnn))," protein-annotation columns (typically exported) were NOT found in this data-set !")
    usColAnn <- wrMisc::naOmit(usColAnn) }
  if(length(usColAnn) >0) { annot <- if(sum(!is.na(usColAnn)) >1) tmp[, wrMisc::naOmit(usColAnn)] else as.matrix(tmp[, wrMisc::naOmit(usColAnn)])
  } else annot <- NULL
  chPrecAA <- !is.na(annot1[,2])
  chFollAA <- !is.na(annot1[,3])
  if(any(chPrecAA)) if("precAA" %in% colnames(annot)) annot[,"precAA"] <- annot1[,2] else annot <- cbind(annot, prec.AA=annot1[,2])
  if(any(chFollAA)) if("follAA" %in% colnames(annot)) annot[,"follAA"] <- annot1[,3] else annot <- cbind(annot, prec.AA=annot1[,3])

  usColCha <- grep("^charge",tolower(colnames(tmp)))           # include charge
  if(length(usColCha) >0) { char <- tmp[,usColCha]
    if(length(usColCha) >1) {     ## more than 1 cols, need to find best col : choose with fewest NAs
      usColCha <- usColCha[which.min(colSums(is.na(char)))] }
    if(debug) message(fxNa," Column for Charge found & added")
    annot <- cbind(annot, Charge=tmp[,usColCha])
  }
  rm(annot1)
  if(debug) {message(fxNa,"rPDP4c .. Done extracting peptide annotation ")
     rPDP4c <- list(fileName=fileName,path=path,chFi=chFi,paFi=paFi,normalizeMeth=normalizeMeth,sampleNames=sampleNames,suplAnnotFile=suplAnnotFile,read0asNA=read0asNA,quantCol=quantCol,cleanDescription=cleanDescription,tmp=tmp,seqCol=seqCol,maSeCo=maSeCo,modifSensible=modifSensible, pepSeq=pepSeq,hasMod=hasMod, annot=annot,AbundCol=AbundCol)}

  ## ABUNDANCE
  abundCols <- grep(AbundCol, colnames(tmp))
  if(length(abundCols) >0) { abund <- if(length(abundCols) >1)  tmp[,abundCols] else {
      matrix(tmp[,abundCols], ncol=1, dimnames=list(rownames(tmp),NULL))}   # how to know column-name if single sample ?
    rownames(abund) <- wrMisc::correctToUnique(pepSeq, silent=silent, callFrom=fxNa)
    ## exculde using excluCol ??
    exclCo <- grepl(excluCol, colnames(abund))
    if(any(exclCo, na.rm=TRUE)) abund <- abund[,which(!exclCo)]
    if(length(dim(abund)) !=2) abund <- matrix(abund, ncol=1, dimnames=list(names(abund),NULL))   # how to know column-name if single sample ?
  } else abund <- NULL

  ## take log2 & normalize
  if(length(abund) >0) { 
    quant <- if(utils::packageVersion("wrMisc") > "1.10") { 
      try(wrMisc::normalizeThis(log2(abund), method=normalizeMeth, mode="additive", refLines=refLi, silent=silent, callFrom=fxNa), silent=TRUE)
    } else try(wrMisc::normalizeThis(log2(abund), method=normalizeMeth, refLines=refLi, silent=silent, callFrom=fxNa), silent=TRUE)       #
    if(debug) { message(fxNa,"rPDP4d .. dim quant: ", nrow(quant)," li and  ",ncol(quant)," cols; colnames : ",wrMisc::pasteC(colnames(quant))," ")} }

  ## PD colnames are typically very cryptic, replace ..
  if(length(sampleNames)==ncol(abund) & all(!is.na(sampleNames)) ) {   # custom sample names given
    colnames(abund) <- colnames(abund) <- sampleNames
    if(length(counts) >0) colnames(counts) <- sampleNames }
    ### GROUPING OF REPLICATES AND SAMPLE META-DATA
    ## META-DATA : read additional annotation & documentation files produced by PD
  if(length(suplAnnotFile) >0 | length(sdrf) >0) {
    setupSd <- readSampleMetaData(sdrf=sdrf, suplAnnotFile=suplAnnotFile, quantMeth="PD", path=path, abund=utils::head(abund), silent=silent, debug=debug, callFrom=fxNa)
  }
  ## finish sample-names: use file-names from meta-data if no custom 'sampleNames' furnished
  if(any(length(sampleNames) <1, length(sampleNames) != ncol(quant), na.rm=TRUE)) {
    if(length(setupSd) >0) {
      tmp <- if(length(setupSd$sdrf) >0) setupSd$sdrf$simpleNa else setupSd$sampleSummary[,"File.Name"]
      colnames(abund) <- wrMisc::trimRedundText(gsub("\\\\","/",tmp))      # no possibility to match colnames
      if(length(counts) >0) colnames(counts) <- colnames(abund) }
  }

  if(debug) {message(fxNa,"rPDP5")
    }
  ## add count-data ?
  counts=NULL
  ## meta-data
  notes <- c(inpFile=paFi, qmethod="ProteomeDiscoverer", qMethVersion=if(length(setupSd) >0) unique(setupSd$Software.Revision) else NA,
        rawFilePath= if(length(setupSd) >0) setupSd$File.Name[1] else NA, normalizeMeth=normalizeMeth, call=match.call(),
    created=as.character(Sys.time()), wrProteo.version=utils::packageVersion("wrProteo"), machine=Sys.info()["nodename"])
  if(debug) {message(fxNa,"rPDP6")}
  ## final output
  if(isTRUE(separateAnnot)) list(raw=NA, quant=quant, pepSeq=pepSeq, annot=annot, counts=counts, notes=notes) else data.frame(quant, annot)
}
  
