#' Read Tabulated Files Exported by ProteomeDiscoverer At Peptide Level, Deprecated
#'
#' Depreciated old version of Peptide identification and quantification results from \href{https://www.thermofisher.com/order/catalog/product/OPTON-30812}{Thermo ProteomeDiscoverer}
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
#' @param gr (character or factor) custom defined pattern of replicate association, will override final grouping of replicates from \code{sdrf} and/or \code{suplAnnotFile} (if provided)   \code{}
#' @param sdrf (character, list or data.frame) optional extraction and adding of experimenal meta-data: if character, this may be the ID at ProteomeExchange,
#'   the second element may give futher indicatations for automatic organization of groups of replicates.
#'   Besides, the output from \code{readSdrf} or a list from \code{defineSamples} may be provided; if \code{gr} is provided, \code{gr} gets priority for grouping of replicates
#' @param suplAnnotFile (logical or character) optional reading of supplemental files produced by ProteomeDiscoverer; however, if \code{gr} is provided, \code{gr} gets priority for grouping of replicates;
#'  if \code{TRUE} defaults to file '*InputFiles.txt' (needed to match information of \code{sdrf}) which can be exported next to main quantitation results;
#'  if \code{character} the respective file-name (relative or absolute path)
#' @param read0asNA (logical) decide if initial quntifications at 0 should be transformed to NA
#' @param quantCol (character or integer) exact col-names, or if length=1 content of \code{quantCol} will be used as pattern to search among column-names for $quant using \code{grep}
#' @param contamCol (character or integer, length=1) which columns should be used for contaminants marked by ProteomeDiscoverer.
#'  If a column named \code{contamCol} is found, the data will be lateron filtered to remove all contaminants, set to \code{NULL} for keeping all contaminants
#' @param refLi (character or integer) custom specify which line of data is main species, if character (eg 'mainSpe'), the column 'SpecType' in $annot will be searched for exact match of the (single) term given
#' @param separateAnnot (logical) if \code{TRUE} output will be organized as list with \code{$annot}, \code{$abund} for initial/raw abundance values and \code{$quant} with final normalized quantitations
#' @param annotCol (character) column names to be read/extracted for the annotation section (default  c("Accession","Description","Gene","Contaminant","Sum.PEP.Score","Coverage....","X..Peptides","X..PSMs","X..Unique.Peptides", "X..AAs","MW..kDa.") )
#' @param FDRCol (list) optional indication to search for protein FDR information
#' @param titGraph (character) custom title to plot
#' @param titGraph (character) depreciated custom title to plot, please use 'tit'
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
#' @seealso \code{\link[utils]{read.table}}, \code{\link[wrMisc]{normalizeThis}}) , \code{\link{readMaxQuantFile}}, \code{\link{readProteomeDiscovererFile}}
#' @examples
#' path1 <- system.file("extdata", package="wrProteo")
#'
#' @export
readProtDiscovPeptides <- function(fileName, path=NULL, normalizeMeth="median", sampleNames=NULL, suplAnnotFile=TRUE, gr=NULL, sdrf=NULL, read0asNA=TRUE, quantCol="^Abundances*",
  annotCol=NULL, contamCol="Contaminant", refLi=NULL, separateAnnot=TRUE, FDRCol=list(c("^Protein.FDR.Confidence","High"), c("^Found.in.Sample.","High")), plotGraph=TRUE,
  titGraph="Proteome Discoverer", wex=1.6, specPref=c(conta="CON_|LYSC_CHICK", mainSpecies="OS=Homo sapiens"), silent=FALSE, debug=FALSE, callFrom=NULL) {  
  message("+++ This function is depreciated, it has been replaced by readProtDiscovererPeptides() from the same package ! \n +++ Synthax and structure of output remain the same ! \n")
  ## read ProteomeDiscoverer exported txt
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readProtDiscovPeptides")

  .Deprecated(new="readProtDiscovererPeptides", package="wrProteo", msg="The function readProtDiscovPeptides() has been deprecated and replaced by readProtDiscovererPeptides() from the same package
     \n +++ Synthax and structure of output remain the same ! \n")

  oparMar <- if(plotGraph) graphics::par("mar") else NULL       # only if figure might be drawn
  reqPa <- c("utils","wrMisc")
  chPa <- sapply(reqPa, requireNamespace, quietly=TRUE)
  if(any(!chPa)) stop("package(s) '",paste(reqPa[which(!chPa)], collapse="','"),"' not found ! Please install first from CRAN")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  excluCol <- c("^Abundance\\.Count","^Abundances\\.Count","^Abundance\\.Ratio","^Abundances\\.Ratio","^Abundance\\.Grouped","^Abundances\\.Grouped")   # exclude this from quantifications columns
  cleanDescription <- TRUE          # clean 'Description' for artifacts of truncated text (tailing ';' etc)
  infoDat <- infoFi <- setupSd <- parametersD <- quant <- counts <- NULL        # initialize
  modifSensible <- TRUE             # separate modified from unmodified peptides (by attaching modif to seq)
 .corPathW <- function(x) gsub("\\\\", "/", x)


  ## check if path & (tsv) file exist
  if(!grepl("\\.txt$|\\.txt\\.gz$", fileName)) message(fxNa,"Trouble ahead, expecting tabulated text file (the file'",fileName,"' might not be right format) !!")
  paFi <- wrMisc::checkFilePath(fileName, path, expectExt="txt", compressedOption=TRUE, stopIfNothing=TRUE, callFrom=fxNa, silent=silent,debug=debug)
  if(debug) message(fxNa,"rPDP2a ..")

  ## prepare for reading files
  if(debug) { message(fxNa,"rPDP3 .. Ready to read", if(length(path) >0) c(" from path ",path[1])," the file  ",fileName[1])}

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
  if(debug) {message(fxNa,"rPDP3z  length(seqCol) ",length(seqCol))
     rPDP3z <- list(tmp=tmp,fileName=fileName,path=path, paFi=paFi,normalizeMeth=normalizeMeth,sampleNames=sampleNames,suplAnnotFile=suplAnnotFile,read0asNA=read0asNA,quantCol=quantCol,seqCol=seqCol,cleanDescription=cleanDescription,tmp=tmp,seqCol=seqCol,modifSensible=modifSensible)}
  maSeCo1 <- match(seqCol, colnames(tmp))
  maSeCo2 <- match(gsub("",".",seqCol), colnames(tmp))
  maSeCo <- if(sum(is.na(maSeCo1)) > sum(is.na(maSeCo2))) maSeCo2 else maSeCo1   # switch betw R-friendly and std
     #quanCo <- "Abundance.F62.Sample.na"
  quantCol <- "^Abundance"              # use as pattern
  IdentTyCol <- "Found.in.Sample"       # use as pattern
  ## need other example for extracting quantifications ?
    #"Confidence.by.Search.Engine.Sequest.HT","Percolator.q.Value.by.Search.Engine.Sequest.HT","Percolator.PEP.by.Search.Engine.Sequest.HT", "XCorr.by.Search.Engine.Sequest.HT","Channel.Occupancy.in.Percent")
  if(debug) {message(fxNa,"rPDP4 ")
     rPDP4 <- list(fileName=fileName,path=path, paFi=paFi,normalizeMeth=normalizeMeth,sampleNames=sampleNames,suplAnnotFile=suplAnnotFile,read0asNA=read0asNA,quantCol=quantCol,seqCol=seqCol,cleanDescription=cleanDescription,tmp=tmp,seqCol=seqCol,maSeCo=maSeCo,modifSensible=modifSensible)}
  .chColNa <- function(x, mat, renameTo=NULL, silent=FALSE, fxNa=NULL){
    ## check in 'matr' for column-name 'x', if required rename best hit (if no direct hit look using grep, then grep wo case); return corrected mat
    chX <- x %in% colnames(mat)
    if(all(chX)) {
      if(is.character(renameTo) && length(renameTo) ==1) colnames(mat)[match(x, colnames(mat))] <- renameTo   # juste simple rename
    } else {     # try to localize column to use
      chX <- grep(x, colnames(mat))
      if(length(chX) >0) {
        if(is.character(renameTo) && length(renameTo) ==1) colnames(mat)[chX[1]] <- renameTo else x
        if(!silent && length(chX) >1) message(fxNa,"Found multiple columns containing '",x,"' : ",wrMisc::pasteC(colnames(mat)[chX], quoteC="'"),", using 1st")
      } else {
        chX <- grep(tolower(x), tolower(colnames(mat)))
        if(length(chX) >0) {
          if(is.character(renameTo) && length(renameTo) ==1) colnames(mat)[chX[1]] <- renameTo else x
          if(!silent && length(chX) >1) message(fxNa,"Found multiple columns containing '",tolower(x),"' : ",wrMisc::pasteC(colnames(mat)[chX], quoteC="'"),", using 1st")
        } else stop("Could NOT find column '",x,"' !!\n  (available columns ",wrMisc::pasteC(colnames(mat), quoteC="'"),")") }
    }
  mat }

  ## EXTRACT PEPTIDE SEQUENCES
  ## extract peptide sequences
  if(debug) {message(fxNa,"rPDP4a .. Ready to start extracting pep seq ")
     rPDP4a <- list(fileName=fileName,path=path, paFi=paFi,normalizeMeth=normalizeMeth,sampleNames=sampleNames,suplAnnotFile=suplAnnotFile,read0asNA=read0asNA,quantCol=quantCol,seqCol=seqCol,cleanDescription=cleanDescription,tmp=tmp,seqCol=seqCol,maSeCo=maSeCo,modifSensible=modifSensible)}
  if(is.na(maSeCo[1])) { if(is.na(maSeCo[2])) {if(!silent) message(fxNa,"Invalid type of data")#; pepSeq <- NULL
    } else pepSeq <- tmp[,maSeCo[2]]
    #else { pepSeq <- tmp[,maSeCo[2]] } #sub("\\.\\[A-Z\\]$", "", sub("^\\[A-Z\\]\\.", "", tmp[,maSeCo[2]])) }
  } else pepSeq <- tmp[,maSeCo[1]]
  fxPrecAA <- function(x) {   ## separate/extract note of preceeding & following AA; take char vector, returns 3-column matrix
    chPre <- grep("^\\[([[:upper:]]|\\-)\\]\\.", x)         # has note of preceeding AA
    chFoll <- grep("\\.\\[([[:upper:]]|\\-)\\]($|_)", x)    # has note of following AA
    out <- cbind(pep=sub("\\.\\[([[:upper:]]|\\-)\\]","", sub("^\\[([[:upper:]]|\\-)\\]\\.","", x)), prec=NA, foll=NA, modifSeq=NA)
    if(length(chPre) >0) out[chPre,2] <- sub(".*\\[","", sub("\\]\\..+","", x[chPre]))   # the preceeding AA
    if(length(chFoll) >0) out[chFoll,3] <- sub("\\].*","", sub(".+\\.\\[","", x[chFoll]))
    out }
  annot1 <- fxPrecAA(pepSeq)      # split
  pepSeq <- annot1[,4] <- annot1[,1]            # also used lateron for rownames of quant

  if(modifSensible) { hasMod <- nchar(tmp[,maSeCo[3]]) >0
    if(any(hasMod, na.rm=TRUE)) annot1[which(hasMod),4] <- gsub(" ","", paste(annot1[which(hasMod),1], tmp[which(hasMod),maSeCo[3]], sep="_"))         # add separator & modification
  }
    #old#if(any(hasMod, na.rm=TRUE)) pepSeq[which(hasMod)] <- paste(pepSeq[which(hasMod)],tmp[which(hasMod),maSeCo[3]],sep="_") }     # modification-separator
  if(debug) {message(fxNa,"Done extracting pep seq    rPDP4b"); rPDP4b <- list(fileName=fileName,path=path, paFi=paFi,normalizeMeth=normalizeMeth,sampleNames=sampleNames,suplAnnotFile=suplAnnotFile,read0asNA=read0asNA,quantCol=quantCol,seqCol=seqCol,pepSeq=pepSeq,annot1=annot1,cleanDescription=cleanDescription,tmp=tmp,seqCol=seqCol,maSeCo=maSeCo,modifSensible=modifSensible) }
  ## ANNOATION (peptide/protein oriented)
   usColAnn <- maSeCo[c(3,6:7,9:14)]

  if(any(is.na(usColAnn), na.rm=TRUE)) {
    if(!silent) message(fxNa,"Note : ",sum(is.na(usColAnn))," protein-annotation columns (typically exported) were NOT found in this data-set !")
    usColAnn <- wrMisc::naOmit(usColAnn) }
  if(length(usColAnn) >0) { annot <- if(sum(!is.na(usColAnn)) >1) tmp[, wrMisc::naOmit(usColAnn)] else as.matrix(tmp[, wrMisc::naOmit(usColAnn)])
  } else annot <- NULL
  chPrecAA <- !is.na(annot1[,2])
  chFollAA <- !is.na(annot1[,3])
  if(any(chPrecAA)) if("precAA" %in% colnames(annot)) annot[,"precAA"] <- annot1[,2] else annot <- cbind(annot, prec.AA=annot1[,2])
  if(any(chFollAA)) if("follAA" %in% colnames(annot)) annot[,"follAA"] <- annot1[,3] else annot <- cbind(annot, foll.AA=annot1[,3])
  annot <- if(ncol(annot1) >3) cbind(annot, seq=annot1[,1], modifSeq=annot1[,4]) else cbind(annot, seq=annot1[,1])
  chDuNa <- duplicated(annot1[,4])
  if(any(chDuNa)) { if(!silent) message(fxNa,"Note : Some 'modifSeq' appear duplicated !!")
      rownames(annot) <-  wrMisc::correctToUnique(annot1[,4], silent=silent, callFrom=fxNa)    # "modifSeq"
  } else rownames(annot) <- annot1[,4]   # "modifSeq"

  usColCha <- grep("^charge", tolower(colnames(tmp)))           # include charge
  if(length(usColCha) >0) { char <- tmp[,usColCha]
    if(length(usColCha) >1) {     ## more than 1 cols, need to find best col : choose with fewest NAs
      usColCha <- usColCha[which.min(colSums(is.na(char)))] }
    if(debug) message(fxNa,"Column for Charge found & added", if(debug) "    rPDP4c")
    annot <- cbind(annot, Charge=tmp[,usColCha])
  }
  rm(annot1)
  if(debug) {message(fxNa,"rPDP4c .. Done extracting peptide annotation ")
     rPDP4c <- list(fileName=fileName,path=path, paFi=paFi,normalizeMeth=normalizeMeth,sampleNames=sampleNames,suplAnnotFile=suplAnnotFile,read0asNA=read0asNA,quantCol=quantCol,cleanDescription=cleanDescription,tmp=tmp,seqCol=seqCol,pepSeq=pepSeq,annot=annot,maSeCo=maSeCo,modifSensible=modifSensible, pepSeq=pepSeq,hasMod=hasMod, annot=annot,quantCol=quantCol)}

  ## ABUNDANCE
    ## locate & extract abundance/quantitation data
    msg <- " CANNOT find ANY quantification columns"
    if(length(quantCol) >1) {
      ## explicit columns (for abundance/quantitation data)
      ## problem : extract '^Abundances*' but NOT 'Abundances.Count.*'
      quantColIni <- quantCol <- grep(quantCol[1], colnames(tmp))
      if(length(quantCol) <1) stop(msg,"  ('",quantCol,"')")
    } else {
      ## pattern search (for abundance/quantitation data)
      if(length(quantCol) <1) { quantCol <- "^Abundance"
        if(!silent) message(fxNa,"Setting argument 'quantCol' to '^Abundance'")}
      quantCol <- grep(quantCol, colnames(tmp))
      if(length(quantCol) <1) quantCol <- grep("^abundance", tolower(colnames(tmp)))
      if(length(quantCol) <1) quantCol <- grep("Intensity$", colnames(tmp))
      if(length(quantCol) <1) quantCol <- grep("intensity$", tolower(colnames(tmp)))
      quantColIni <- quantCol
      if(length(quantCol) <1) stop(msg," specified in argument 'quantCol' !") }
    ## check for columns to exclude (like 'Abundances.Count.')
    if(length(excluCol) >1) {
      excCo <- unique(unlist(lapply(excluCol, grep, colnames(tmp))))
      if(length(excCo) >0) {
        quantCol <- quantCol[-wrMisc::naOmit(match(excCo, quantCol))]
        if(length(quantCol) <1) stop(msg," (all match to 'excluCol')") else {
          if(!silent) message(fxNa,"Removed ",length(quantColIni) -length(quantCol)," columns")}
      }
    }

  if(length(quantCol) >0) { abund <- if(length(quantCol) >1)  tmp[,quantCol] else {
      matrix(tmp[,quantCol], ncol=1, dimnames=list(rownames(tmp),NULL))}   # how to know column-name if single sample ?
    rownames(abund) <-  rownames(annot) #wrMisc::correctToUnique(pepSeq, silent=silent, callFrom=fxNa)
    ## check for columns to exclude (like 'Abundances.Count.')
    if(length(excluCol)==1) {
      excCo <- grep(excluCol, colnames(tmp))
      if(any(duplicated(excCo, quantCol), na.rm=TRUE)) {
        quantCol <- quantCol[-match(excCo, quantCol)]
        if(length(quantCol) <1) stop(msg," (all match to 'excluCol')") else {
          if(!silent) message(fxNa,"Removed ",length(quantColIni) -length(quantCol)," columns")}
      }
    }
    abund <- as.matrix(tmp[,quantCol])                      # abundance val
    rownames(abund) <- rownames(annot)
    if(debug) {message(fxNa,"rPDP8 ..  "); rPDP8 <- list(tmp=tmp,annot=annot,specPref=specPref,abund=abund,quantCol=quantCol)}

    ## check & clean abudances
    chNorm <- grep("\\.Normalized\\.", colnames(abund))
    if(length(chNorm)*2 == ncol(abund)) {              # in case Normalized makes 1/2 of columns use non-normalized
      abund <- abund[,-chNorm]
    }
    colnames(abund) <- sub("^Abundances\\.Normalized\\._{0,1}|^abundances\\.Normalized\\._{0,1}|^Abundances{0,1}_{0,1}|^abundances{0,1}_{0,1}","",colnames(abund))
    chNum <- is.numeric(abund)
    if(!chNum) {abund <- apply(tmp[,quantCol], 2, wrMisc::convToNum, convert="allChar", silent=silent, callFrom=fxNa)}


    ## remove heading 'X..' from headers (only if header won't get duplicated
    ### why here ??? 24mar23
    chXCol <- grep("^X\\.\\.",colnames(annot))
    if(length(chXCol) >0) {
      newNa <- sub("^X\\.\\.","",colnames(annot)[chXCol])
      chDu <- duplicated(c(newNa, colnames(annot)), fromLast=TRUE)
      if(any(chDu, na.rm=TRUE)) newNa[which(chDu)] <- colnames(annot)[chXCol][which(chDu)]
      colnames(annot)[chXCol] <- newNa }
    ## remove heading/tailing spaces (first look which columns might be subject to this treatment)
    ch1 <- list(A=grep("^ +",annot[1,]), B=grep("^ +",annot[2,]), C=grep("^ +",annot[floor(mean(nrow(annot))),]), D=grep("^ +",annot[nrow(annot),]) )
    chCo <- unique(unlist(ch1))
    annot[,chCo] <- sub("^ +","",sub(" +$","",annot[,chCo]))   # remove heading/tailing spaces
    if(debug) { message(fxNa,"rPDP9 .. dim annot ",nrow(annot)," and ",ncol(annot)); rPDP9 <- list(annot=annot,tmp=tmp,abund=abund,sampleNames=sampleNames,specPref=specPref,annotCol=annotCol,contamCol=contamCol,infoDat=infoDat) }

    ## add custom sample names (if provided)
    if(length(sampleNames) ==ncol(abund) && ncol(abund) >0) {
      if(debug) { message(fxNa,"rPDP9b") }
      if(length(unique(sampleNames)) < length(sampleNames)) {
        if(!silent) message(fxNa,"Custom sample names not unique, correcting to unique")
        sampleNames <- wrMisc::correctToUnique(sampleNames, callFrom=fxNa) }
      colnames(abund) <- sampleNames
      if(debug) { message(fxNa,"rPDP9c") }
    } else {
      colnames(abund) <- sub("Abundance\\.F[[:digit:]]+\\.Sample\\.|Abundances\\.F[[:digit:]]+\\.Sample\\.","Sample.", colnames(abund))
    }
  } else abund <- NULL


  ## take log2 & normalize
  if(length(abund) >0) {
    quant <- if(utils::packageVersion("wrMisc") > "1.10") {
      try(wrMisc::normalizeThis(log2(abund), method=normalizeMeth, mode="additive", refLines=refLi, silent=silent, callFrom=fxNa), silent=TRUE)
    } else try(wrMisc::normalizeThis(log2(abund), method=normalizeMeth, refLines=refLi, silent=silent, callFrom=fxNa), silent=TRUE)       #
    if(debug) { message(fxNa,"rPDP9d .. dim quant: ", nrow(quant)," li and  ",ncol(quant)," cols; colnames : ",wrMisc::pasteC(colnames(quant))," ")} }

  ## PD colnames are typically very cryptic, replace ..
  if(length(sampleNames)==ncol(abund) && all(!is.na(sampleNames)) ) {   # custom sample names given
    colnames(abund) <- colnames(abund) <- sampleNames
    if(length(counts) >0) colnames(counts) <- sampleNames }

  ### GROUPING OF REPLICATES AND SAMPLE META-DATA
  ## META-DATA : read additional annotation & documentation files produced by PD
  if(length(suplAnnotFile) >0 || length(sdrf) >0) {
    setupSd <- readSampleMetaData(sdrf=sdrf, suplAnnotFile=suplAnnotFile, quantMeth="PD", path=path, abund=utils::head(abund), silent=silent, debug=debug, callFrom=fxNa)
  }

    ## finish groups of replicates & annotation setupSd
    setupSd <- .checkSetupGroups(abund=abund, setupSd=setupSd, gr=gr, sampleNames=sampleNames, quantMeth="PD", silent=silent, debug=debug, callFrom=fxNa)
    colNa <- if(length(setupSd$sampleNames)==ncol(abund)) setupSd$sampleNames else setupSd$groups
    chGr <- grepl("^X[[:digit:]]", colNa)                                                # check & remove heading 'X' from initial column-names starting with digits
    if(any(chGr)) colNa[which(chGr)] <- sub("^X","", colNa[which(chGr)])                 # 
    colnames(quant) <- colnames(abund) <- colNa
    if(length(setupSd$sampleNames)==ncol(abund)) setupSd$sampleNames <- colNa else setupSd$groups <- colNa
    if(length(dim(counts)) >1 && length(counts) >0) colnames(counts) <- colNa

    if(debug) {message(fxNa,"Read sample-meta data, rPDP14"); rPDP14 <- list(sdrf=sdrf,suplAnnotFile=suplAnnotFile,abund=abund, quant=quant,refLi=refLi,annot=annot,setupSd=setupSd,sampleNames=sampleNames)}

    ## main plotting of distribution of intensities
    custLay <- NULL
    if(is.numeric(plotGraph) && length(plotGraph) >0) {custLay <- as.integer(plotGraph); plotGraph <- TRUE} else {
        if(!isTRUE(plotGraph)) plotGraph <- FALSE}
    if(plotGraph) .plotQuantDistr(abund=abund, quant=quant, custLay=custLay, normalizeMeth=normalizeMeth, softNa="Proteome Discoverer",
      refLi=refLi, refLiIni=nrow(abund), tit=titGraph, silent=silent, callFrom=fxNa, debug=debug)

    ## meta-data
    notes <- c(inpFile=paFi, qmethod="ProteomeDiscoverer", qMethVersion=if(length(infoDat) >0) unique(infoDat$Software.Revision) else NA,
    	rawFilePath= if(length(infoDat) >0) infoDat$File.Name[1] else NA, normalizeMeth=normalizeMeth, call=match.call(),
      created=as.character(Sys.time()), wrProteo.version=utils::packageVersion("wrProteo"), machine=Sys.info()["nodename"])
    ## final output
    if(isTRUE(separateAnnot)) list(raw=abund, quant=quant, annot=annot, counts=counts, sampleSetup=setupSd, quantNotes=parametersD, notes=notes) else data.frame(quant,annot)
}

