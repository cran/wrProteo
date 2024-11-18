#' Read Tabulated Files Exported by Ionbot At Peptide Level
#'
#' This function allows importing initial petide identification and quantification results from  \href{https://ionbot.cloud/about}{Ionbot} 
#' which were exported as tabulated tsv can be imported and relevant information extracted.
#' The final output is a list containing 3 main elements: \code{$annot}, \code{$raw} and optional \code{$quant}, or returns data.frame with entire content of file if \code{separateAnnot=FALSE}.
#'
#' @details
#' Using the argument \code{suplAnnotFile} it is possible to specify a specific file (or search for default file) to read for extracting file-names as sample-names 
#' and other experiment realted information.
#'
#' 
#'
#'
#' @param fileName (character) name of file to be read
#' @param path (character) path of file to be read
#' @param normalizeMeth (character) normalization method, defaults to \code{median}, for more details see \code{/link[wrMisc]{normalizeThis}}
#' @param sampleNames (character) new column-names for quantification data (ProteomeDiscoverer does not automatically use file-names from spectra); this argument has priority over \code{suplAnnotFile}
#' @param gr (character or factor) custom defined pattern of replicate association, will override final grouping of replicates from \code{sdrf} and/or \code{suplAnnotFile} (if provided) 
#' @param sdrf (character, list or data.frame) optional extraction and adding of experimenal meta-data: if character, this may be the ID at ProteomeExchange,
#'   the second & third elements may give futher indicatations for automatic organization of groups of replicates.
#'   Besides, the output from \code{readSdrf} or a list from \code{defineSamples} may be provided;
#'   if \code{gr} is provided, \code{gr} gets priority for grouping of replicates;
#'   if \code{sdrfOrder=TRUE} the output will be put in order of sdrf
#' @param read0asNA (logical) decide if initial quntifications at 0 should be transformed to NA
#' @param quantCol (character or integer) exact col-names, or if length=1 content of \code{quantCol} will be used as pattern to search among column-names for $quant using \code{grep}
#' @param contamCol (character or integer, length=1) which columns should be used for contaminants marked by ProteomeDiscoverer.
#'  If a column named \code{contamCol} is found, the data will be lateron filtered to remove all contaminants, set to \code{NULL} for keeping all contaminants
#' @param refLi (character or integer) custom specify which line of data is main species, if character (eg 'mainSpe'), the column 'SpecType' in $annot will be searched for exact match of the (single) term given
#' @param separateAnnot (logical) if \code{TRUE} output will be organized as list with \code{$annot}, \code{$abund} for initial/raw abundance values and \code{$quant} with final normalized quantitations
#' @param annotCol (character) column names to be read/extracted for the annotation section (default 
#'   c("Accession","Description","Gene","Contaminant","Sum.PEP.Score","Coverage....","X..Peptides","X..PSMs","X..Unique.Peptides", "X..AAs","MW..kDa.") )
#' @param FDRCol (list) optional indication to search for protein FDR information
#' @param groupPref (list) additional parameters for interpreting meta-data to identify structure of groups (replicates), will be passed to \code{readSampleMetaData}.
#'   May contain \code{lowNumberOfGroups=FALSE} for automatically choosing a rather elevated number of groups if possible (defaults to low number of groups, ie higher number of samples per group)
#'   May contain \code{chUnit} (logical or character) to be passed to \code{readSampleMetaData()} for (optional) adjustig of unit-prefixes in meta-data group labels, in case multiple different unit-prefixes 
#'   are used (eg '100pMol' and '1nMol').
#' @param titGraph (character) custom title to plot
#' @param titGraph (character) depreciated custom title to plot, please use 'tit'
#' @param wex (integer) relative expansion factor of the violin-plot (will be passed to \code{/link[wrGraph]{vioplotW}})
#' @param suplAnnotFile (logical or character) optional reading of supplemental files produced by ProteomeDiscoverer; however, if \code{gr} is provided, \code{gr} gets priority for grouping of replicates;
#'  if \code{TRUE} defaults to file '*InputFiles.txt' (needed to match information of \code{sdrf}) which can be exported next to main quantitation results;
#'  if \code{character} the respective file-name (relative or absolute path)
#' @param specPref (character or list) define characteristic text for recognizing (main) groups of species (1st for comtaminants - will be marked as 'conta', 2nd for main species- marked as 'mainSpe',
#'   and optional following ones for supplemental tags/species - maked as 'species2','species3',...);
#'   if list and list-element has multiple values they will be used for exact matching of accessions (ie 2nd of argument \code{annotCol})
#' @param plotGraph (logical or integer) optional plot of type vioplot of initial and normalized data (using \code{normalizeMeth}); if integer, it will be passed to \code{layout} when plotting
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This function returns a list with \code{$raw} (initial/raw abundance values), \code{$quant} with final normalized quantitations, \code{$annot}, \code{$counts} an array with number of peptides, \code{$quantNotes}
#'   and \code{$notes}; or if \code{separateAnnot=FALSE} the function returns a data.frame with annotation and quantitation only
#' @seealso \code{/link[utils]{read.table}}, \code{/link{readMaxQuantFile}}, \code{/link{readProteomeDiscovererFile}}, \code{/link[wrMisc]{normalizeThis}}) 
#' @examples
#' path1 <- system.file("extdata", package="wrProteo")
#' fiIonbot <- "tinyIonbotFile1.tsv.gz"
#' datIobPep <- readIonbotPeptides(fiIonbot, path=path1) 
#'
#' @export
readIonbotPeptides <- function(fileName, path=NULL, normalizeMeth="median", sampleNames=NULL, gr=NULL, sdrf=NULL, read0asNA=TRUE, quantCol="^Abundances*",
  annotCol=NULL, contamCol="Contaminant", refLi=NULL, separateAnnot=TRUE, FDRCol=list(c("^Protein.FDR.Confidence","High"), c("^Found.in.Sample.","High")), plotGraph=TRUE, 
  suplAnnotFile=TRUE, groupPref=list(lowNumberOfGroups=TRUE, chUnit=TRUE), titGraph="Ionbot", wex=1.6, specPref=c(conta="CON_|LYSC_CHICK", mainSpecies="OS=Homo sapiens"), silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## read ProteomeDiscoverer exported txt
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readIonbotPeptides")    # Ionbot
  oparMar <- if(plotGraph) graphics::par("mar") else NULL       # only if figure might be drawn

  reqPa <- c("utils","wrMisc")
  chPa <- sapply(reqPa, requireNamespace, quietly=TRUE)
  if(any(!chPa)) stop("package(s) '",paste(reqPa[which(!chPa)], collapse="','"),"' not found ! Please install first from CRAN")
  if(isTRUE(debug)) silent <- FALSE else { debug <- FALSE
    if(!isTRUE(silent)) silent <- FALSE }
  excluCol <- c()
  #excluCol <- c("^Abundance//.Count","^Abundances//.Count","^Abundance//.Ratio","^Abundances//.Ratio","^Abundance//.Grouped","^Abundances//.Grouped")   # exclude this from quantifications columns
  cleanDescription <- TRUE          # clean 'Description' for artifacts of truncated text (tailing ';' etc)
  infoDat <- infoFi <- setupSd <- parametersD <- quant <- counts <- NULL        # initialize
  modifSensible <- TRUE             # separate modified from unmodified peptides (by attaching modif to seq)
 .corPathW <- function(x) gsub("\\", "/", x)

  contamFilter <- TRUE             # filter contaminants away
  
  ## check if path & (tsv) file exist
  #if(!grepl("(\\.tsv$)|(\\.tsv\\.gz$))", fileName)) message(fxNa,"Trouble ahead, expecting '.tsv' tabulated text file (the file'",fileName,"' might not be right format) !!")
  if(!grepl("\\.tsv(\\.gz){0,1}$", fileName)) message(fxNa,"Trouble ahead, expecting '.tsv' tabulated text file (the file'",fileName,"' might not be right format) !!")
  paFi <- wrMisc::checkFilePath(fileName, path, expectExt="tsv", mode=c("compressed"), callFrom=fxNa, silent=silent,debug=debug)  # compressedOption=TRUE, stopIfNothing=TRUE
  if(debug) { message(fxNa,"rIBP2a .."); rIBP2a <- list(paFi=paFi,fileName=fileName, path=path, normalizeMeth=normalizeMeth, sampleNames=sampleNames,suplAnnotFile=suplAnnotFile, gr=gr, 
    sdrf=sdrf,quantCol=quantCol,annotCol=annotCol,contamCol=contamCol,refLi=refLi,separateAnnot=separateAnnot,FDRCol=FDRCol )}    

  ## prepare for reading files
  if(debug) { message(fxNa,"rIBP3 .. Ready to read", if(length(path) >0) c(" from path '",path[1]),"' the file  '",fileName[1],"'")}

  ## read (main) file
  ## future: look for fast reading of files
  tmp <- try(utils::read.delim(file.path(paFi), stringsAsFactors=FALSE, header=TRUE), silent=TRUE)
  if(inherits(tmp, "try-error")) stop("Unable to read input file ('",paFi,"')!")
  if(debug) { message(fxNa,"rIBP3b .. dims of initial data : ", nrow(tmp)," li and ",ncol(tmp)," col ")}
  if(debug) { message(fxNa,"rIBP3b .. Ready to start extracting pep seq"); rIBP3b <- list(tmp=tmp,paFi=paFi,fileName=fileName, path=path, normalizeMeth=normalizeMeth, sampleNames=sampleNames,suplAnnotFile=suplAnnotFile, gr=gr, 
    sdrf=sdrf,quantCol=quantCol,annotCol=annotCol,contamCol=contamCol,refLi=refLi,separateAnnot=separateAnnot,FDRCol=FDRCol )}
  
  ## EXTRACT PEPTIDE SEQUENCES
  ## extract peptide sequence - not needed, stripped seq already in col 'Base.Sequence'
  #pepSe <- gsub("\\[[[:digit:]].+","", tmp[,"Sequence"])   # ionbot strip modif name & position (eg 'AAAGSVLLEDCK[4]Carbamidomethyl[C]')
  
  ## no obj precAA or postAA !!
  if(debug) { message(fxNa,"rIBP3c ..")}

  ## other peptide/protein info
  #txtCol <- c("Modifications", "Master.Protein.Accessions","Positions.in.Master.Proteins","Master.Protein.Descriptions")
  seqCol <- c("Sequence",	"Base.Sequence","Protein.Groups","Gene.Names","Organism")
  if(debug) {message(fxNa,"rIBP3z  length(seqCol) ",length(seqCol))
     rIBP3z <- list(tmp=tmp,fileName=fileName,path=path, paFi=paFi,normalizeMeth=normalizeMeth,sampleNames=sampleNames,suplAnnotFile=suplAnnotFile,read0asNA=read0asNA,quantCol=quantCol,seqCol=seqCol,cleanDescription=cleanDescription,tmp=tmp,seqCol=seqCol,modifSensible=modifSensible)}
  maSeCo1 <- match(seqCol, colnames(tmp))
  maSeCo2 <- match(gsub("",".",seqCol), colnames(tmp))
  maSeCo <- if(sum(is.na(maSeCo1)) > sum(is.na(maSeCo2))) maSeCo2 else maSeCo1   # switch betw R-friendly and std
  names(maSeCo) <- seqCol
     #quanCo <- "Abundance.F62.Sample.na"
  quantCol <- "Intensity_"                  # use as pattern
  identifTyCol <- "Detection\\.Type_"       # use as pattern (Ionbot specific, eg 'Detection.Type_QEKAC160601_34' ) ==>  STILL NEED TO INTEGRATE  !!!

  ## need other example for extracting quantifications ?
  if(debug) {message(fxNa,"rIBP4 ")
     rIBP4 <- list(fileName=fileName,path=path, paFi=paFi,normalizeMeth=normalizeMeth,sampleNames=sampleNames,suplAnnotFile=suplAnnotFile,read0asNA=read0asNA,quantCol=quantCol,seqCol=seqCol,cleanDescription=cleanDescription,tmp=tmp,seqCol=seqCol,maSeCo=maSeCo,modifSensible=modifSensible)}

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
        } else stop("Could NOT find column '",x,"' !!/n  (available columns ",wrMisc::pasteC(colnames(mat), quoteC="'"),")") }
    }
  mat }

  if(debug) {message(fxNa,"rIBP4a  ")
     rIBP4a <- list(fileName=fileName,path=path, paFi=paFi,normalizeMeth=normalizeMeth,sampleNames=sampleNames,suplAnnotFile=suplAnnotFile,read0asNA=read0asNA,quantCol=quantCol,seqCol=seqCol,cleanDescription=cleanDescription,tmp=tmp,seqCol=seqCol,maSeCo=maSeCo,modifSensible=modifSensible)}
     
  ## extract peptide sequences ????
  if(is.na(maSeCo[1])) { if(is.na(maSeCo[2])) {if(!silent) message(fxNa,"Invalid type of data")
    } else pepSeq <- tmp[,maSeCo[2]]
  } else pepSeq <- tmp[,maSeCo[1]]
  ## fxPrecAA not needed with ionbot ...
    
  annot1 <- cbind(tmp[,wrMisc::naOmit(maSeCo[2:1])], modif=NA, tmp[,wrMisc::naOmit(maSeCo[3:length(maSeCo)])]) 
  hasMod <- grep("\\[[[:digit:]]", tmp[,maSeCo[1]])
  if(length(hasMod) >0) annot1[hasMod,3] <- substring(tmp[hasMod, maSeCo[1]], nchar(tmp[hasMod, maSeCo[2]]) +1)

  ## modifs : tmp[,1] eg 'AAFTECCQAADKAACLLPK[4]Carbamidomethyl[C]|[4]Carbamidomethyl[C]|[4]Carbamidomethyl[C]'
  if(modifSensible) { hasMod <- grepl("\\[[[:digit:]]+\\]", tmp[,maSeCo[1]])   
    ## should one keep the modifications for sep col ??
  }
  if(debug) {message(fxNa,"Done extracting pep seq    rIBP4b"); rIBP4b <- list(fileName=fileName,path=path, paFi=paFi,normalizeMeth=normalizeMeth,sampleNames=sampleNames,suplAnnotFile=suplAnnotFile,read0asNA=read0asNA,quantCol=quantCol,seqCol=seqCol,pepSeq=pepSeq,annot1=annot1,cleanDescription=cleanDescription,tmp=tmp,seqCol=seqCol,maSeCo=maSeCo,modifSensible=modifSensible,specPref=specPref) }

    ## Ionbot : Organism & Gene.Nameq are all NA !!
    ## annot$Protein.Groups) :  "sp|P00925|ENO2_YEAST" , "sp|P0C0T4|RS25B_YEAST__sp|Q3E792|RS25A_YEAST",  
    ## UniProt : db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion 
  annot <- cbind(annot1, uniqueIdentifier=NA, entryName=NA, proteinName=NA, Species=NA, Contam=NA) 
  ## extract UniqueIdentifier & EntryName  (note ProteinName only avail from fasta)
  chAnn <- nchar(annot1[,"Protein.Groups"])
  if(any(chAnn) >0) {
    strictSpecPattern <- TRUE
    paUI <- "^[[:upper:]]+[[:digit:]]+([[:upper:]]|[[:digit:]])*"
    paEN <- if(isTRUE(strictSpecPattern)) "^[[:upper:]]+[[:digit:]]([[:upper:]]|[[:digit:]])*_[[:upper:]]+$" else "^[[:upper:]]+[[:digit:]]([[:upper:]]|[[:digit:]])*(_[[:upper:]]+){0,1}$"
    annX <- sub("([[:lower:]]+)|([[:upper:]]+[[:lower:]]+)\\|","",annot1[,"Protein.Groups"])    # remove heading db-sign
    
    annot[,"entryName"] <- sub("\\|.+","", sub(paste0(paUI,"\\|"),"", annX))   # first remove heading UniqueIdent and than reomve info from mult prot
    annot[,"uniqueIdentifier"] <- sub(paste0(paEN,"\\|*.+"),"", annX)         # remove EntryName (2nd pos) and all thereafter (info from mult prot)
    chUI <- grepl(paUI, annot[,"uniqueIdentifier"])
    if(any(!chUI)) annot[which(!chUI),"uniqueIdentifier"] <- NA   # double-check  (& set NA if not fitting)
    chEN <- grepl(paEN, annot[,"entryName"])
    if(any(!chEN)) annot[which(!chEN),"entryName"] <- NA          # double-check  (& set NA if not fitting
    rm(annX)
    ## add info from fasta
        
  }
  ## extract Species from col Protein.Groups
  .annSpecies <- function(spe=c("_HUMAN","Homo sapiens"), anno, exCoNa) {
    ## extract species tags out of annot[,"Majority.protein.IDs"], place as convert to regular name in anno, return matrix anno
    ch1 <- grep(spe[1], anno[,exCoNa[1]])
    if(length(ch1) >0) anno[ch1,exCoNa[2]] <- spe[2]  # eg "Homo sapiens"
    anno }
  commonSpec <- .commonSpecies()
  for(i in 1:nrow(commonSpec)) annot[,c("Protein.Groups","Species")] <- .annSpecies(commonSpec[i,], anno=annot[,c("Protein.Groups","Species")], exCoNa=1:2) 
  
  ## look for tags from  specPref & create col SpecType
  if(length(specPref) >0) {
    ## look if available, for specif tags 
    annot <- .extrSpecPref(specPref, annot, useColumn=c("Protein.Groups","Species"), silent=silent, debug=debug, callFrom=fxNa)   # set annot[,"specPref"] according to specPref
  } else if(debug) message(fxNa,"Note: Argument 'specPref' not specifed (empty)")
  
  ## mine info to add 'EntryName' and 'Accession' 
  if("Protein.Groups" %in% colnames(annot)) {
    ## problem : may contain multiple ids from mult prot (separated by '__') - use last 
    entNa <- sub(".+\\|","", annot[,"Protein.Groups"])    # get last (assuming as 'EntryName')
    annot <- cbind(Accession=sub(".+\\|","", substr(annot[,"Protein.Groups"], 1, nchar(annot[,"Protein.Groups"]) - nchar(entNa) -1)), EntryName=entNa, annot )  # get entry before last, assuming as 'Accession'
  }
  ## 
    
  ##  include charge - not avail in ionbot.peptides
  rm(annot1)
  if(debug) {message(fxNa,"rIBP4c .. Done extracting peptide annotation ")
     rIBP4c <- list(fileName=fileName,path=path, paFi=paFi,normalizeMeth=normalizeMeth,sampleNames=sampleNames,suplAnnotFile=suplAnnotFile,read0asNA=read0asNA,quantCol=quantCol,cleanDescription=cleanDescription,tmp=tmp,seqCol=seqCol,pepSeq=pepSeq,annot=annot,maSeCo=maSeCo,modifSensible=modifSensible, pepSeq=pepSeq,hasMod=hasMod, annot=annot,specPref=specPref)}

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
    if(debug) {message(fxNa,"rIBP8 ..  "); rIBP8 <- list(tmp=tmp,annot=annot,specPref=specPref,abund=abund,quantCol=quantCol)}

    ## check & clean abudances
    chNorm <- grep("//.Normalized//.", colnames(abund))
    if(length(chNorm)*2 == ncol(abund)) {              # in case Normalized makes 1/2 of columns use non-normalized
      abund <- abund[,-chNorm]
    }
    colnames(abund) <- sub("^Abundances//.Normalized//._{0,1}|^abundances//.Normalized//._{0,1}|^Abundances{0,1}_{0,1}|^abundances{0,1}_{0,1}","",colnames(abund))
    chNum <- is.numeric(abund)
    if(!chNum) {abund <- apply(tmp[,quantCol], 2, wrMisc::convToNum, convert="allChar", silent=silent, callFrom=fxNa)}


    ## remove heading 'X..' from headers (only if header won't get duplicated
    ### why here ??? 24mar23
    chXCol <- grep("^X//.//.",colnames(annot))
    if(length(chXCol) >0) {
      newNa <- sub("^X//.//.","",colnames(annot)[chXCol])
      chDu <- duplicated(c(newNa, colnames(annot)), fromLast=TRUE)
      if(any(chDu, na.rm=TRUE)) newNa[which(chDu)] <- colnames(annot)[chXCol][which(chDu)]
      colnames(annot)[chXCol] <- newNa }
    ## remove heading/tailing spaces (first look which columns might be subject to this treatment)
    ch1 <- list(A=grep("^ +",annot[1,]), B=grep("^ +",annot[2,]), C=grep("^ +",annot[floor(mean(nrow(annot))),]), D=grep("^ +",annot[nrow(annot),]) )
    chCo <- unique(unlist(ch1))
    annot[,chCo] <- sub("^ +","",sub(" +$","",annot[,chCo]))   # remove heading/tailing spaces
    if(debug) { message(fxNa,"rIBP9 .. dim annot ",nrow(annot)," and ",ncol(annot)); rIBP9 <- list(annot=annot,tmp=tmp,abund=abund,sampleNames=sampleNames,specPref=specPref,annotCol=annotCol,contamCol=contamCol,infoDat=infoDat) }

    ## add custom sample names (if provided)
    if(length(sampleNames) ==ncol(abund) && ncol(abund) >0) {
      if(debug) { message(fxNa,"rIBP9b") }
      if(length(unique(sampleNames)) < length(sampleNames)) {
        if(!silent) message(fxNa,"Custom sample names not unique, correcting to unique")
        sampleNames <- wrMisc::correctToUnique(sampleNames, callFrom=fxNa) }
      colnames(abund) <- sampleNames
      if(debug) { message(fxNa,"rIBP9c") }
    } else {
      colnames(abund) <- sub("Abundance//.F[[:digit:]]+//.Sample//.|Abundances//.F[[:digit:]]+//.Sample//.","Sample.", colnames(abund))
    }
  } else abund <- NULL

  ## check for reference for normalization (need annotation)
  ##### ADJUST !!
  refLiIni <- refLi
  if(is.character(refLi) && length(refLi)==1) { refLi <- which(annot[,"SpecType"]==refLi)
    if(length(refLi) <1) message(fxNa,"Could not find any peptide matching argument 'refLi', ignoring ...") else {
      if(!silent) message(fxNa,"Normalize using subset of ",length(refLi)) } }           # may be "mainSpe"
  if(length(refLi) <1) refLi <- NULL

  ## take log2 & normalize
  if(length(abund) >0) {
    quant <- if(utils::packageVersion("wrMisc") > "1.10") {
      try(wrMisc::normalizeThis(log2(abund), method=normalizeMeth, mode="additive", refLines=refLi, silent=silent, callFrom=fxNa), silent=TRUE)
    } else try(wrMisc::normalizeThis(log2(abund), method=normalizeMeth, refLines=refLi, silent=silent, callFrom=fxNa), silent=TRUE)       #
    if(debug) { message(fxNa,"rIBP9d .. dim quant: ", nrow(quant)," li and  ",ncol(quant)," cols; colnames : ",wrMisc::pasteC(colnames(quant))," ")} }

  ##  colnames may be cryptic, replace .. (needed with IB ??)
  if(length(sampleNames)==ncol(abund) && all(!is.na(sampleNames)) ) {   # custom sample names given
    colnames(abund) <- colnames(abund) <- sampleNames
    if(length(counts) >0) colnames(counts) <- sampleNames }
  if(debug) {message(fxNa,"  rIBP12b"); rIBP12b <- list(sdrf=sdrf,suplAnnotFile=suplAnnotFile,abund=abund, quant=quant,refLi=refLi,annot=annot,sampleNames=sampleNames)}

  ### GROUPING OF REPLICATES AND SAMPLE META-DATA
  ## META-DATA : read additional annotation & documentation files produced by IB
    #sdrfFi <- wrMisc::checkFilePath(fileName, path, expectExt="tsv", compressedOption=TRUE, stopIfNothing=TRUE, callFrom=fxNa, silent=silent,debug=debug)

  if(length(suplAnnotFile) >0 || length(sdrf) >0) {
    if(length(sampleNames) %in% c(1, ncol(abund))) groupPref$sampleNames <- sampleNames
    if(length(gr) %in% c(1, ncol(abund))) groupPref$gr <- gr
    setupSd <- readSampleMetaData(sdrf=sdrf, suplAnnotFile=suplAnnotFile, quantMeth="IB", path=path, abund=utils::head(quant), chUnit=isTRUE(groupPref$chUnit), groupPref=groupPref, silent=silent, debug=debug, callFrom=fxNa)
  }
  if(debug) {message(fxNa,"  rIBP13"); rIBP13 <- list(sdrf=sdrf,suplAnnotFile=suplAnnotFile,abund=abund, quant=quant,refLi=refLi,annot=annot,sampleNames=sampleNames,setupSd=setupSd)}

    ## finish groups of replicates & annotation setupSd
    setupSd <- .checkSetupGroups(abund=utils::head(abund), setupSd=setupSd, gr=gr, sampleNames=sampleNames, quantMeth="IB", silent=silent, debug=debug, callFrom=fxNa)
    if(debug) {message(fxNa,"  rIBP13b"); rIBP13b <- list(sdrf=sdrf,suplAnnotFile=suplAnnotFile,abund=abund, quant=quant,refLi=refLi,annot=annot,sampleNames=sampleNames,setupSd=setupSd)}

    ## harmonize sample-names/1
    colNa <- if(length(setupSd$sampleNames)==ncol(abund)) setupSd$sampleNames else setupSd$groups


    ## option : choose (re-)naming of levels & columns on most redundance
    if(length(setupSd$sampleNames)==ncol(quant) && length(setupSd$sampleNaSdrf)==ncol(quant)) {
      ## check if setupSd$sampleNaSdrf  or   setupSd$sampleNames contain better info (some info on replicates)
      chRed1 <- sum(duplicated(sub("(_|\\.| )[[:digit:]]+.*","", setupSd$sampleNaSdrf)), na.rm=TRUE)
      chRed2 <- sum(duplicated(sub("(_|\\.| )[[:digit:]]+.*","", setupSd$sampleNames)), na.rm=TRUE)
      if(chRed2 < chRed1) {                                          # use info for levels depending on where more 
        colNa <- colnames(abund) <- setupSd$sampleNames <- setupSd$sampleNaSdrf                     ## take sample names from sdrf via  setupSd$sampleNaSdrf
      } else {
        colNa <- colnames(abund) <- setupSd$sampleNames }          ## take sample names from sdrf via  setupSd$sampleNaSdrf
      #setupSd$level  
    }    
    ## option : set order of samples as sdrf
    if("sdrfOrder" %in% names(sdrf) && isTRUE(as.logical(sdrf["sdrfOrder"])) && length(setupSd$iniSdrfOrder)==ncol(abund) && ncol(abund) >1) {  # set order according to sdrf (only if >1 samples)
      nOrd <- order(setupSd$iniSdrfOrder)
      abund <- abund[,nOrd]
      if(length(quant) >0) quant <- quant[,nOrd]
      #setupSd$level <- setupSd$level[nOrd]
      ## rename columns according to sdrf and set order of quant and abund ..
      ## now adapt order of setupSd, incl init Sdrf
      if(length(setupSd) >0) { 
        is2dim <- sapply(setupSd, function(x,le) length(dim(x))==2 && nrow(x)==le, le=length(nOrd))    # look for matr or df to change order of lines
        if(any(is2dim) >0) for(i in which(is2dim)) setupSd[[i]] <- setupSd[[i]][nOrd,]
        isVe <- sapply(setupSd, function(x,le) length(x)==le && length(dim(x)) <1, le=length(nOrd))    # look for vector to change order in setupSd
        if(any(isVe) >0) for(i in which(isVe)) setupSd[[i]] <- setupSd[[i]][nOrd] }
      gr <- gr[nOrd]

      if(length(counts) >0 && length(dim(counts))==3) counts <- array(counts[,nOrd,], dim=c(nrow(counts), length(nOrd), dim(counts)[3]), 
        dimnames=list(rownames(counts), colnames(counts)[nOrd], dimnames(counts)[[3]]))
      if(debug) {message(fxNa,"rIBP13c .."); rIBP13c <- list(path=path,chPa=chPa,tmp=tmp,groupPref=groupPref,quantCol=quantCol,abund=abund,chNum=chNum,
        sdrf=sdrf,quant=quant,annot=annot,setupSd=setupSd)}
        
      ## try re-adjusting levels
      tm1 <- sub("^[[:alpha:]]+( |_|-|\\.)+[[:alpha:]]+","", colnames(abund))  # remove heading text
      if(all(grepl("^[[:digit:]]", tm1))) {
        tm1 <- try(as.numeric(sub("( |_|-|\\.)*[[:alpha:]].*","", tm1)), silent=TRUE)   # remove tailing text and try converting to numeric
        if(!inherits(tm1, "try-error")) {
          setupSd$level <- match(tm1, sort(unique(tm1)))
          names(setupSd$level) <- tm1
          if(!silent) message(fxNa,"Sucessfully re-adjusted levels after bringing in order of Sdrf")}
      }     
    } else {     # no sdrf-info

      ## harmonize sample-names/2
      colNa <- colnames(abund)
      chGr <- grepl("^X[[:digit:]]", colNa)                                                # check & remove heading 'X' from initial column-names starting with digits
      if(any(chGr)) colNa[which(chGr)] <- sub("^X","", colNa[which(chGr)])                 #
      colnames(quant) <- colNa
      if(length(abund) >0) colnames(abund) <- colNa  
    }  
    if(length(setupSd$sampleNames)==ncol(abund)) setupSd$sampleNames <- colNa #no#else setupSd$groups <- colNa
    if(length(dim(counts)) >1 && length(counts) >0) colnames(counts) <- colNa

    if(debug) {message(fxNa,"Read for plotting, rIBP14"); rIBP14 <- list(sdrf=sdrf,suplAnnotFile=suplAnnotFile,abund=abund, quant=quant,refLi=refLi,annot=annot,setupSd=setupSd,sampleNames=sampleNames, plotGraph=plotGraph,normalizeMeth=normalizeMeth,titGraph=titGraph)}

    ## main plotting of distribution of intensities
    custLay <- NULL
    if(is.numeric(plotGraph) && length(plotGraph) >0) {custLay <- as.integer(plotGraph); plotGraph <- TRUE} else {
        if(!isTRUE(plotGraph)) plotGraph <- FALSE}
    if(plotGraph) .plotQuantDistr(abund=abund, quant=quant, custLay=custLay, normalizeMeth=normalizeMeth, softNa="Ionbot peptides",
      refLi=refLi, refLiIni=nrow(abund), tit=titGraph, silent=silent, callFrom=fxNa, debug=debug)

    ## meta-data
    notes <- c(inpFile=paFi, qmethod="Ionbot", qMethVersion=if(length(infoDat) >0) unique(infoDat$Software.Revision) else NA,
    	identType="peptide", rawFilePath= if(length(infoDat) >0) infoDat$File.Name[1] else NA, normalizeMeth=normalizeMeth, call=deparse(match.call()),
      created=as.character(Sys.time()), wrProteo.version=utils::packageVersion("wrProteo"), machine=Sys.info()["nodename"])

    ## final output
    if(isTRUE(separateAnnot)) list(raw=abund, quant=quant, annot=annot, counts=counts, sampleSetup=setupSd, quantNotes=parametersD, notes=notes) else data.frame(quant,annot)
}
