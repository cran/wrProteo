#' Read (Normalized) Quantitation Data Files Produced By AlphaPept 
#'
#' Protein quantification results from \href{https://github.com/MannLabs/alphapept}{AlphaPept} can be read using this function.
#' Input files compressed as .gz can be read as well.
#' The protein abundance values (XIC) get extracted. Since protein annotation is not very extensive with this format of data, the function allows reading the
#' initial fasta files (from the directory above the quantitation-results) allowing to extract more protein-annotation (like species).
#' Sample-annotation (if available) can be extracted from  sdrf files, too.
#' The protein abundance values may be normalized using multiple methods (median normalization as default), the determination of normalization factors can be restricted to specific proteins
#' (normalization to bait protein(s), or to invariable matrix of spike-in experiments).
#' The protein annotation data gets parsed to extract specific fields (ID, name, description, species ...).
#' Besides, a graphical display of the distribution of protein abundance values may be generated before and after normalization.
#'
#' @details
#' 
#' Meta-data describing the samples and experimental setup may be available from a sdrf-file (from the directory above the analysis/quantiication results)
#' If available, the meta-data will be examined for determining groups of replicates and
#' the results thereof can be found in $sampleSetup$levels.
#' Alternatively, a dataframe formatted like sdrf-files (ie for each sample a separate line, see also function \code{readSdrf}) may be given, too.
#'
#' This import-function has been developed using AlphaPept version x.x.
#' The final output is a list containing these elements: \code{$raw}, \code{$quant}, \code{$annot}, \code{$counts}, \code{$sampleSetup}, \code{$quantNotes}, \code{$notes}, or (if \code{separateAnnot=FALSE}) data.frame
#'   with annotation- and main quantification-content. If \code{sdrf} information has been found, an add-tional list-element \code{setup}
#' will be added containg the entire meta-data as \code{setup$meta} and the suggested organization as \code{setup$lev}.
#'
#'
#' @param fileName (character) name of file to be read (default 'results_proteins.csv'). Gz-compressed files can be read, too.
#' @param path (character) path of file to be read
#' @param fasta (logical or character) if \code{TRUE} the (first) fasta from one direcory higher than \code{fileName} will be read as fasta-file to extract further protein annotation;
#'    if \code{character} a fasta-file at this location will be read/used/
#' @param isLog2 (logical) typically data read from AlphaPept are expected NOT to be \code{isLog2=TRUE}
#' @param normalizeMeth (character) normalization method, defaults to \code{median}, for more details see \code{\link[wrMisc]{normalizeThis}})
#' @param extrColNames (character or \code{NULL}) custom definition of col-names to extract
#' 
#' 
#' @param quantCol (character or integer) exact col-names, or if length=1 content of \code{quantCol} will be used as pattern to search among column-names for $quant using \code{grep}
#' @param contamCol (character or integer, length=1) which columns should be used for contaminants
#' @param read0asNA (logical) decide if initial quntifications at 0 should be transformed to NA (thus avoid -Inf in log2 results)
#' @param sampleNames (character) custom column-names for quantification data; this argument has priority over \code{suplAnnotFile}
#' @param specPref (character) prefix to identifiers allowing to separate i) recognize contamination database, ii) species of main identifications and iii) spike-in species
#' @param refLi (character or integer) custom specify which line of data should be used for normalization, ie which line is main species; if character (eg 'mainSpe'), the column 'SpecType' in $annot will be searched for exact match of the (single) term given
#' @param remRev (logical) option to remove all protein-identifications based on reverse-peptides
#' @param remConta (logical) option to remove all proteins identified as contaminants
#' @param separateAnnot (logical) if \code{TRUE} output will be organized as list with \code{$annot}, \code{$abund} for initial/raw abundance values and \code{$quant} with final normalized quantitations
#' @param gr (character or factor) custom defined pattern of replicate association, will override final grouping of replicates from \code{sdrf} and/or \code{suplAnnotFile} (if provided)   \code{}
#' @param sdrf (character, list or data.frame) optional extraction and adding of experimenal meta-data: if character, this may be the ID at ProteomeExchange,
#'   the second & third elements may give futher indicatations for automatic organization of groups of replicates.
#'   Besides, the output from \code{readSdrf} or a list from \code{defineSamples} may be provided;
#'   if \code{gr} is provided, \code{gr} gets priority for grouping of replicates;
#'   if \code{sdrfOrder=TRUE} the output will be put in order of sdrf
#' @param suplAnnotFile (logical or character) optional reading of supplemental files produced by Compomics; if \code{gr} is provided, it gets priority for grouping of replicates
#'  if \code{TRUE} default to files 'summary.txt' (needed to match information of \code{sdrf}) and 'parameters.txt' which can be found in the same folder as the main quantitation results;
#'  if \code{character} the respective file-names (relative ro absolute path), 1st is expected to correspond to 'summary.txt' (tabulated text, the samples as given to Compomics) and 2nd to 'parameters.txt' (tabulated text, all parameters given to Compomics)  
#' @param groupPref (list) additional parameters for interpreting meta-data to identify structure of groups (replicates), will be passed to \code{readSampleMetaData}.
#'   May contain \code{lowNumberOfGroups=FALSE} for automatically choosing a rather elevated number of groups if possible (defaults to low number of groups, ie higher number of samples per group)
#'   May contain \code{chUnit} (logical or character) to be passed to \code{readSampleMetaData()} for (optional) adjustig of unit-prefixes in meta-data group labels, in case multiple different unit-prefixes 
#'   are used (eg '100pMol' and '1nMol').
#' @param plotGraph (logical) optional plot vioplot of initial and normalized data (using \code{normalizeMeth}); alternatively the argument may contain numeric details that will be passed to \code{layout} when plotting
#' @param titGraph (character) custom title to plot of distribution of quantitation values
#' @param wex (numeric)  relative expansion factor of the violin in plot
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a list with  \code{$raw} (initial/raw abundance values), \code{$quant} with final normalized quantitations, \code{$annot} (columns ), \code{$counts} an array with 'PSM' and 'NoOfRazorPeptides',
#'   \code{$quantNotes}, \code{$notes} and optional \code{setup} for meta-data from \code{sdrf}; or a data.frame with quantitation and annotation if \code{separateAnnot=FALSE}
#' @seealso \code{\link[utils]{read.table}}, \code{\link[wrMisc]{normalizeThis}}) , \code{\link{readProteomeDiscovererFile}}; \code{\link{readProlineFile}} (and other import-functions), \code{\link{matrixNAinspect}}
#' @examples
#' path1 <- system.file("extdata", package="wrProteo")
#' # Here we'll load a short/trimmed example file
#' fiNaAP <- "tinyAlpaPeptide.csv.gz"
#' dataAP <- readAlphaPeptFile(file=fiNaAP, path=path1, tit="tiny AlphaPaptide ")
#' summary(dataAP$quant)
#' @export
readAlphaPeptFile <- function(fileName="results_proteins.csv", path=NULL, fasta=NULL, isLog2=FALSE, normalizeMeth="none", quantCol="_LFQ$", contamCol=NULL,
  read0asNA=TRUE, refLi=NULL, sampleNames=NULL,                                       # pepCountCol=c("number_of_peptides"), extrColNames=c("protein_group"),
  specPref=NULL, extrColNames=NULL,
  remRev=TRUE, remConta=FALSE, separateAnnot=TRUE, gr=NULL, sdrf=NULL, suplAnnotFile=NULL, groupPref=list(lowNumberOfGroups=TRUE, chUnit=TRUE),
  titGraph=NULL, wex=1.6, plotGraph=TRUE, silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## prepare
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readAlphaPeptFile")
  oparMar <- graphics::par("mar")                     # old margins, for rest after figure
  oparLayout <- graphics::par("mfcol")                # old layout, for rest after figure
  on.exit(graphics::par(mar=oparMar, mfcol=oparLayout))            # restore old mar settings
  remStrainNo <- TRUE                   # if TRUE extract Species in very stringent pattern
  cleanDescription <- TRUE              # clean 'Description' for artifacts of truncated text (tailing ';' etc)
  fixSpeciesNames <- TRUE  
  trimColNames <- TRUE      ## further trim quantitation colnames
  chCol <- NULL

  ## functions
  .cleanMQann <- function(x, sep="\\|", silent=FALSE, debug=FALSE, callFrom=NULL) {
    ## split multiple protein entries as with 1st column of MaxQuant data
    ## return matrix with
    ## example ann1 <- read.delim(file.path(system.file("extdata", package="wrProteo"), "tinyWombCompo1.csv.gz"), sep=",", stringsAsFactors=FALSE)[,1]
    ##   .cleanMQann(ann1)
    #  x=rAP4a$tmp[c(5,31:32,81:82,111:114),1]
    xIni <- x       # keep backup for recupera-ting bizzare nonparsed
    isCont <- grepl("CON__", x)
    mult <- nchar(x) - nchar(gsub(";", "", x))
    chMult <- mult >0
    if(any(chMult)) {
      spl1 <- strsplit(x[which(chMult)], ";")
      ## use entry with most separators (when multiple entries, eg 'sp|P00761|CON__TRYP_PIG;CON__P00761')
      spl1 <- sapply(spl1, function(y) { nSep <- nchar(y) - nchar(gsub("|","",y)); y[which.max(nSep)] })
      x[which(chMult)] <- spl1 }
    ## split separators
    chSpl <- function(y) {chID <- grepl("^[[:upper:]]{1,3}[[:digit:]]{2,}|^[[:upper:]]{1,3}[[:digit:]]+[[:upper:]]+[[:digit:]]*", y); chName <- grepl("[A-Z0-9]_[[:upper:]]",y);   # extract db, ID & prot-name
      c(dbIni= if((length(y) >1 && grepl("^[[:lower:]]{1,8}$", y[1])) || length(y) >2 && grepl("^[[:lower:]]{2}|[[:lower:]]{2}$",
        y[1])) y[1] else NA, IDini=if(any(chID)) y[which(chID)[1]] else NA, nameIni=if(any(chName)) y[which(chName)[1]] else NA) }
    x <- t(sapply(strsplit(x, sep), chSpl))
    nColIni <- ncol(x)
    cleanID <- function(y, useCol=c(db=1, ID=2, name=3)) {
      ext <- grepl("[[:lower:]]+$", y[,useCol[2]])    # look for extension like 'P08758ups'
      extNoDb <- which(ext & is.na(y[,useCol[1]]))
      if(any(ext)) { cleanID <- sub("[[:lower:]]+$","", y[which(ext), useCol[2]])
        if(length(extNoDb) >0) y[which(ext), useCol[1]] <- substring(y[which(ext), useCol[2]], nchar(cleanID) +1 )
        y[which(ext), useCol[2]] <- cleanID }
      prefi <- grepl("^[[:upper:]]+__[[:upper:]]", y[,useCol[3]])       # look for prefix like 'CON__FA5_BOVIN'
      if(any(prefi)) { ch2 <- grepl("[A-Z0-9]_[[:upper:]]",  y[which(prefi), useCol[3]]); if(any(ch2)) {
        y[which(prefi)[which(ch2)], useCol[1]] <- tolower(sub("__[[:upper:]].+","", y[which(prefi)[which(ch2)], useCol[3]]))
        y[which(prefi)[which(ch2)], useCol[3]] <- sub("^[[:upper:]]+__","", y[which(prefi)[which(ch2)], useCol[3]])}}
      colnames(y) <- c("db","ID","name")
      y }
    x <- cbind(x, cleanID(x, useCol=c(db=1, ID=2, name=3)))
    x <- cbind(x, conta=grepl("^con|^REV_", x[,"db"]) | grepl("__CON__",xIni))
    ## recuperate all (bizarre) non-parsed into ID
    isNa <- rowSums(is.na(x)) > nColIni -2
    if(any(isNa)) x[which(isNa),c(2+nColIni)] <- xIni[which(isNa)]
    cbind(x[,c((nColIni+1):ncol(x), 1:nColIni)], iniSoftAnn=xIni) }

  ## end functions


  ## init check
  reqPa <- c("utils","wrMisc")
  chPa <- sapply(reqPa, requireNamespace, quietly=TRUE)
  if(any(!chPa)) stop("Package(s) '",paste(reqPa[which(!chPa)], collapse="','"),"' not found ! Please install first from CRAN")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
         excluCol <- "^Abundances.Count"   # exclude this from quantifications columns
  cleanDescription <- TRUE                 # clean 'Description' for artifacts of truncated text (tailing ';' etc)
  infoDat <- infoFi <- setupSd <- parametersD <- annot <- annotMQ <- annotMQ <- NULL        # initialize

  ## check if path & file exist
  paFi <- wrMisc::checkFilePath(fileName, path, expectExt="csv", compressedOption=TRUE, stopIfNothing=TRUE, callFrom=fxNa, silent=silent,debug=debug)
  ## read (main) file
  ## future: look for fast reading of files
    #  read.csv("C:\\E\\projects\\MassSpec\\smallProj\\testAlphaPept\\demoData_aug23\\testHSE_A1-3\\results_proteins.csv", header=TRUE)
  tmp <- try(utils::read.csv(paFi, header=TRUE), silent=TRUE)

  if(length(tmp) <1 || inherits(tmp, "try-error") || length(dim(tmp)) <2) {
    if(inherits(tmp, "try-error")) warning("Unable to read input file ('",paFi,"')!  (check format or if rights to read)") else {
      if(!silent) message(fxNa,"Content of  file '",paFi,"' seeps empty or non-conform !  Returning NULL; check if this is really a Compomics-file") }
    tmp <- NULL
    return(NULL)
  } else {
    ## start checking format of initial read 
    if(debug) { message(fxNa,"rAP1 .. dims of initial data : ", nrow(tmp)," li and ",ncol(tmp)," col "); rAP1 <- list(fileName=fileName,path=path,paFi=paFi,tmp=tmp,normalizeMeth=normalizeMeth,read0asNA=read0asNA,quantCol=quantCol,
      refLi=refLi,separateAnnot=separateAnnot   )}  # annotCol=annotCol,FDRCol=

    ## check which columns can be extracted (for annotation)
    isSummaryCsv <- any(grepl("^n_sequence\\.", colnames(tmp)))    # check if "results_protein_summary.csv (has cols heading '^n_sequence\\.' and cols '$LFQ.intensity.' and' $intensity. 
    ## note that 'results_proteins.csv' has 'XXX_LFQ' and 'XXX'
    if(length(extrColNames) <1) { extrColNames <-            # default as list for patterns
      if(isSummaryCsv) list(lfq="^LFQ\\.intensity\\.",nSeq="^n_sequence\\.",annot=1) else list(LFQ="_LFQ$", annot=1)
    } else {
      ## check custom entry (& format as list)
      ch1 <- extrColNames %in% colnames(tmp)
      if(!all(ch1) && "LFQ" %in% names(extrColNames)) {       # seems to be pattern ?

      } 
      if(!"annot" %in% names(extrColNames)) if(is.list(extrColNames)) extrColNames$annot <- 1 else c(extrColNames, annot=1)
    
    }
    if(debug) { message(fxNa,"rAP1b"); rAP1b <- list(tmp=tmp,paFi=paFi,fasta=fasta,quantCol=quantCol,extrColNames=extrColNames,isSummaryCsv=isSummaryCsv)}

  if(length(tmp) >0) {
    ## check for lines with absent IDs  => eliminate
    chNa <- is.na(tmp[,extrColNames$annot]) | nchar(tmp[,extrColNames$annot]) <2
    if(any(chNa)) {
      if(!silent) message(fxNa,"Removing ",sum(chNa)," lines since absent ID or all NA  (won't be able to do anything lateron withour ID ..)")
      tmp <- tmp[-which(chNa),]
    }
    if(debug) {message(fxNa,"rAP1d"); rAP1d <- list(tmp=tmp,paFi=paFi,fasta=fasta,quantCol=quantCol)}

    ## further extracting : quantitation
    useDCol <- grep(extrColNames$LFQ, colnames(tmp))
    if(length(useDCol) <1) stop("NO columns matching term ",wrMisc::pasteC(quantCol, quoteC="'")," from argument 'quantCol' found !")
    abund <- as.matrix(tmp[,useDCol])               #  abundances  (not normalized, not log2)
    if(debug) {message(fxNa,"rAP1d2"); rAP1d2 <- list(tmp=tmp,paFi=paFi,fasta=fasta,abund=abund,quantCol=quantCol,extrColNames=extrColNames)}

    chNum <- try(is.numeric(abund), silent=TRUE)
    if(!chNum) {abund <- try(apply(tmp[,useDCol], 2, wrMisc::convToNum, convert="allChar", silent=debug, callFrom=fxNa), silent=TRUE)
      if(inherits(abund, "try-error")) {datOK <- FALSE; warning(fxNa,"CANNOT transform 'abund' to numeric data !'")} }
    if(length(dim(abund)) <2 && !is.numeric(abund)) abund <- matrix(as.numeric(abund), ncol=ncol(abund), dimnames=dimnames(abund))
    if(debug) {message(fxNa,"rAP1e"); rAP1e <- list(tmprAP1e=tmp,abund=abund,paFi=paFi,fasta=fasta,quantCol=quantCol,extrColNames=extrColNames)}
    
    ## adjust colnames of abund
    colnames(abund) <- sub(extrColNames$LFQ,"", sub("^LFQ.intensity\\.","", colnames(abund)))

    if(trimColNames) {  ## further trim
      colnames(abund) <- wrMisc::.trimLeft(wrMisc::.trimRight( sub(paste0("^",quantCol),"", colnames(abund))))
      ## no trim needed for AlphaPept ?
    }
    if(debug) {message(fxNa,"rAP3"); rAP3 <- list(abund=abund,paFi=paFi,path=path,chPa=chPa,tmp=tmp,remConta=remConta,fasta=fasta)}
    ## check for neg values
    chNeg <- which(abund  <0)
    if(length(chNeg) ==prod(dim(abund))) { read0asNA <- FALSE
      message(fxNa,"NOTE : Bizzare, ALL values are NEGATIVE !!  omit transforming neg values to NA)")}
    ## convert 0 to NA
    if(!isFALSE(read0asNA)) {
      ch1 <- abund <= 0
      if(any(ch1, na.rm=TRUE)) { abund[which(ch1)] <- NA
        if(!silent) message(fxNa,"Transform ",sum(ch1),"(",100*round(sum(ch1)/length(ch1),3),"%) initial '0' values to 'NA'")}}

    ## further extracting : prepare for countig data
    if("nSeq" %in% names(extrColNames)) {
      counts <- as.matrix(tmp[,extrColNames$nSeq])
      if(is.numeric(counts)) counts <- try(matrix(as.numeric(counts, ncol=ncol(counts)), dimnames=dimnames(counts)))
      if(inherits(counts, "try-error")) { counts <- NULL
        message(fxNa,"Unable to extract counts (not numeric ?)")}
    } else { counts <- NULL}

    if(debug) {message(fxNa,"rAP4"); rAP4 <- list(abund=abund,counts=counts,annot=annot,paFi=paFi,path=path,chPa=chPa,tmp=tmp,remConta=remConta,fasta=fasta,annotMQ=annotMQ)}


    ## Annotation
    annot <- .cleanMQann(tmp[,extrColNames$annot])    # typically 1st column entitled 'X'
    colnames(annot)[2] <- "protein_group"
   
    if(debug) {message(fxNa,"rAP4aa"); rAP4aa <- list(abund=abund,counts=counts,annot=annot,paFi=paFi,path=path,chPa=chPa,tmp=tmp,remConta=remConta,fasta=fasta,annot=annot)}

    ## read fasta from same dir (AlphaPept)
    if(length(fasta) >0) {fasta <- fasta[1]; if(isFALSE(fasta) || is.na(fasta)) fasta <- NULL}
    if(isTRUE(fasta)) {
      chFa <- grep("\\.fasta$", dirname(paFi))
      faFi <- if(length(chFa) >0) dir(dirname(paFi), pattern="\\.fasta$")[1]
    } else faFi <- fasta
    if(length(faFi) >0) {     # has fasta for recuperating annotation
      fasta <- try(readFasta2(filename=faFi, tableOut=TRUE, silent=silent,debug=debug,callFrom=fxNa), silent=TRUE)
      ## Potential problem with inconsistent format of fasta
      if(inherits(fasta, "try-error")) { fasta <- NULL
        if(!silent) message(fxNa,"Unable to read/open fasta file '",faFi,"'  (check rights to read ?)")
      } else {
        tmpAnn <- if(length(annot) >0) annot[,2] else {if(length(annot) >0) annot[,2] else tmp[,1]}      # 'P02768' still missing
        tm2 <- wrMisc::concatMatch(tmpAnn, fasta[,2], sepPattern=NULL, globalPat="digitExtension", silent=silent, debug=debug, callFrom=fxNa)   # clean protein-names (eg digit extensions, concateneated IDs) & match to data
        iniAnn <- if(length(annot) >0) annot else {if(length(annot) >0) annot else cbind(iniSoftAnn=tmp[,1])}
        #colnames(iniAnn) <- c("iniSoftAnn", if(ncol(iniAnn) >1) paste0(colnames(iniAnn)[-1],".",quantSoft))

        useFaCol <- match(c("uniqueIdentifier","entryName","proteinName","OS","OX","GN","database"), colnames(fasta))             # do not export full 'sequence'
        annot <- cbind(trimIdentifier=names(tm2), fasta[tm2, useFaCol], iniAnn=tmpAnn)
        if(debug) {message(fxNa,"rAP4ab"); rAP4ab <- list(abund=abund,counts=counts,annot=annot,fasta=fasta,tmp=tmp,chNa=chNa,annot=annot) }

        foundFastaCol <- !is.na(useFaCol)
        if(any(foundFastaCol)) colnames(annot)[1:sum(foundFastaCol)] <- c("Accession","AccessionFull","Description","EntryName","Species","OX","GeneName","Database")[which(foundFastaCol)]
        ##  strip species details
        if("Species" %in% annot) annot[,"Species"] <- sub(" \\(.+", "", annot[,"Species"])
      }
    } else {
      #annot <- if(length(annot) >0) annot 
      if(debug) message(fxNa,"NOTE : No fasta-file found in main directory ...")
    }
    if(debug) {message(fxNa,"dim annot  ",nrow(annot)," ",ncol(annot),"  rAP4b"); rAP4b <- list(annot=annot,faFi=faFi,abund=abund,tmp=tmp,fasta=fasta,annot=annot)}

    ## check ID col of annot
    #chID <- match(c("Accession","protein_group","uniqueIdentifier"), colnames(annot))
    #if(all(is.na(chID))) { chID <- wrMisc::naOmit(match(c("protein_group","ID"), colnames(annot)))
    #  if(length(chID) < 1) warning("PROBLEM : UNEXPECTED colnames in annot") #else colnames(annot)[chID][1] <- "Accession"
    #}
    annot <- as.matrix(tmp[,1])

    ## remove lines wo IDs
    chNa <- is.na(annot[,1])
    if(any(chNa)) {
      if(!silent) message(fxNa,"Removing ",sum(chNa)," out of ",nrow(abund)," lines wo ID")
      rmLi <- which(chNa)
      tmp <- tmp[-rmLi,]
      annot <- annot[-rmLi,]
      if(length(dim(annot)) <2) annot <- matrix(annot, ncol=1, dimnames=list(NULL,colnames(tmp)[1]))
      abund <- abund[-rmLi,]
      if(length(counts) >0) counts <- if(length(dim(counts))==3) counts[-rmLi,,] else counts[-rmLi,]
    }
    if(debug) {message(fxNa,"dim annot",nrow(annot)," ",ncol(annot),"  rAP4d"); rAP4d <- list(annot=annot,faFi=faFi,abund=abund,tmp=tmp)}

    ## unique ID
    chD <- duplicated(annot[,1])
    uniqueID <- if(any(chD, na.rm=TRUE)) wrMisc::correctToUnique(annot[,1], silent=silent, callFrom=fxNa) else annot[,1]  # extrColNames[1]
    rownames(annot) <- rownames(abund) <- uniqueID
    if(length(counts) >0) rownames(counts) <- uniqueID
    if(debug) {message(fxNa,"rAP4e"); rAP4e <- list(paFi=paFi,path=path,chPa=chPa,tmp=tmp,counts=counts,
      quantCol=quantCol,abund=abund,chNum=chNum,annot=annot,remConta=remConta,specPref=specPref)}

    ## remove Wombat contaminants
   #useColumn <- wrMisc::naOmit(match(c("Accession","protein_group"), colnames(annot)))
    conLi <- grep("CON__[[:alnum:]]", annot[, if(ncol(annot) >1) wrMisc::naOmit(match(c("Accession","protein_group"), colnames(annot)))[1] else 1])
    if(remConta) {
      if(length(conLi) >0) {
        iniLi <- nrow(annot)
        annot <- annot[-conLi,]
        abund <- abund[-conLi,]
        counts <- if(length(dim(counts))==3) counts[-conLi,,] else counts[-conLi,,]
        if(debug) message(fxNa,"Removing ",length(conLi)," instances of contaminants to final ",nrow(annot)," lines/IDs")}
    }

    ## split Annotation
    if(debug) {message(fxNa,"rAP4f"); rAP4f <- list(path=path,chPa=chPa,tmp=tmp,counts=counts,
      quantCol=quantCol,abund=abund,chNum=chNum,annot=annot,remConta=remConta,specPref=specPref)}

    ## finalize annotation
    chCols <- c("EntryName","GeneName","Species","Contam","Description")
    chCol2 <- chCols %in% colnames(annot)
    if(any(!chCol2)) annot <- cbind(annot, matrix(NA, nrow=nrow(annot), ncol=sum(!chCol2), dimnames=list(NULL, chCols[which(!chCol2)]))) # add columns so far not present
    if(!remConta && length(conLi) >0) annot[conLi, "Contam"] <- "TRUE"

    if(debug) {message(fxNa,"rAP5"); rAP5 <- list(path=path,chPa=chPa,tmp=tmp,chCol=chCol,counts=counts,
      quantCol=quantCol,abund=abund,chNum=chNum,annot=annot,remConta=remConta,remStrainNo=remStrainNo, specPref=specPref)}

    ## extract species according to custom search parameters 'specPref'
    if(remStrainNo && any(!is.na(annot[,"Species"]))) {
      annot[,"Species"] <- sub(" \\(strain [[:alnum:]].+","", annot[,"Species"])
    }
    ## complete species annot   by info extracted from fasta : ' OS='
    .completeSpeciesAnnot <- function(spe=c("Homo sapiens", "_HUMAN"), anno=annot, exCoNa=c("Species", "EntryName","name","proteinName")) {    # re-written 12jun23
      ## complete species if missing in anno[,exCoNa[2]] but found in anno[,exCoNa[1]]; return corrected anno
      chNa <- is.na(anno[,exCoNa[1]]) | nchar(anno[,exCoNa[1]]) <1             # missing (species) annotation
      if(any(chNa, na.rm=TRUE)) {        # suppose that all 'exCoNa' are present as colnames in 'annot'
        useColumn <- if(all(is.na(anno[,exCoNa[2]]))) wrMisc::naOmit(match(exCoNa[3:length(exCoNa)], colnames(anno))) else exCoNa[2]
        if(length(useColumn) >1) useColumn <- useColumn[1]
        chS <- grep(spe[1], anno[,useColumn])
        if(length(chS) >0) anno[chS, exCoNa[1]] <- spe[2]
      }
      anno }
    if(isTRUE(fixSpeciesNames)) {          # try to recuperate/fix non-given/bad formatted species
      chNa <- is.na(annot[,"Species"])
      if(any(chNa)) {
        commonSpec <- .commonSpecies()
        for(i in 1:nrow(commonSpec)) annot[which(chNa),] <- .completeSpeciesAnnot(commonSpec[i,], annot[which(chNa),], exCoNa=c("Species","EntryName","name","proteinName")) }
      if(debug) {message(fxNa,"rAP6"); rAP6 <- list(path=path,chPa=chPa,tmp=tmp,chCol=chCol,counts=counts,
        quantCol=quantCol,abund=abund,chNum=chNum,annot=annot,remConta=remConta,remStrainNo=remStrainNo, specPref=specPref)}

      ## check/complete for truncated species names (ie names found inside other ones)
      chSpe <- which(!is.na(annot[,"Species"]) & nchar(annot[,"Species"]) >0)
      if(length(chSpe) >0) {
        OS <- gsub(";{1,5}$", "", annot[chSpe,"Species"])              # remove tailing separators
        OSna <- unique(OS)
        ch1 <- nchar(OSna) <1
        if(debug) {message(fxNa,"rAP6b")}
        if(any(ch1, na.rm=TRUE)) OSna <- OSna[which(nchar(OSna) >0)]     # (just in case) remove empty tags
        ch2 <- lapply(OSna, grep, OSna)
        chTr <- sapply(ch2, length) >1
        if(any(chTr, na.rm=TRUE)) { if(!silent) message(fxNa,"Found ",sum(chTr)," species name(s) appearing inside other ones, assume as truncated (eg  ",OSna[which(chTr)[1]],")")
          for(i in which(chTr)) OS[which(OS==OSna[i])] <- OSna[ch2[[i]][1]]
        }
        annot[chSpe,"Species"] <- OS }
    }
    ## in case "Accession" is avail not "EntryName" is not

    if(debug) {message(fxNa,"rAP7"); rAP7 <- list(path=path,chPa=chPa,tmp=tmp,chCol=chCol,quantCol=quantCol,remStrainNo=remStrainNo,
      abund=abund,chNum=chNum,specPref=specPref, annot=annot,remConta=remConta,counts=counts)}

    ## look for tags from  specPref
    if(length(specPref) >0) {
      ## set annot[,"specPref"] according to specPref
      annot <- .extrSpecPref(specPref, annot, useColumn=c("Description","Species","EntryName","GeneName"), silent=silent, debug=debug, callFrom=fxNa)
    } else if(debug) message(fxNa,"Note: Argument 'specPref' not specifed (empty)")
    if(debug) {message(fxNa,"rAP7b") }

    if(!silent) { chSp <- sum(is.na(annot[,"Species"]))
    if(chSp >0) message(fxNa,"Note: ",chSp," proteins with unknown species")

    tab <- table(annot[,"Species"])
    if(length(tab) >0) {
      tab <- rbind(names(tab), paste0(": ",tab,",  "))
      if(!silent) message("     data by species : ", apply(tab, 2, paste)) } }              # all lines assigned

    if(debug) {message(fxNa,"rAP8"); rAP8 <- list(path=path,chPa=chPa,tmp=tmp,chCol=chCol,quantCol=quantCol,remStrainNo=remStrainNo,
      abund=abund,chNum=chNum, annot=annot,remConta=remConta,counts=counts) }

    ## look for unique col from $annot to use as rownames
    if(nrow(annot) <1) warning("annot is empty (NO lines)")
     ## maybe annot is empty ?

    chAn <- colSums(apply(annot[,c(1:min(ncol(annot),7))], 2, duplicated), na.rm=TRUE)          # look at first 6 cols : how many elements per column duplicated
    if(!silent) message(fxNa,"Use column '",colnames(annot)[which.min(chAn)],"' as identifyer (has fewest, ie ",chAn[which.min(chAn)]," duplicated entries) as rownames")
    rownames(abund) <- rownames(annot) <- if(any(chAn==0)) annot[,which(chAn==0)[1]] else wrMisc::correctToUnique(annot[,which.min(chAn)], callFrom=fxNa)
    if(length(counts) >0) rownames(counts) <- rownames(annot)
    if(debug) {message(fxNa,"rAP9"); rAP9 <- list(path=path,chPa=chPa,tmp=tmp,chCol=chCol,quantCol=quantCol,abund=abund,chNum=chNum,
      annot=annot,refLi=refLi,remConta=remConta,normalizeMeth=normalizeMeth)}

    ## check for reference for normalization
    refLiIni <- refLi
    if(is.character(refLi) && length(refLi)==1) {
       refLi <- which(annot[,"SpecType"]==refLi)
      if(length(refLi) <1 ) { refLi <- 1:nrow(abund)
        if(!silent) message(fxNa,"Could not find any proteins matching argument 'refLi=",refLiIni,"', ignoring ...")
      } else {
        if(!silent) message(fxNa,"Normalize using (custom) subset of ",length(refLi)," lines specified as '",refLiIni,"'")}}    # may be "mainSpe"

    ## take log2 & normalize
    quant <- try(wrMisc::normalizeThis(if(isLog2) abund else log2(abund), method=normalizeMeth, mode="additive", refLines=refLi, silent=silent, debug=debug, callFrom=fxNa), silent=TRUE)
    if(inherits(quant, "try-error")) { warning(fxNa,"PROBLEMS ahead : Unable to normalize as log2-data !!") }

    if(debug) {message(fxNa,"rAP10"); rAP10 <- list(path=path,chPa=chPa,tmp=tmp,chCol=chCol,quantCol=quantCol,abund=abund,chNum=chNum,
      quant=quant,annot=annot,remConta=remConta,groupPref=groupPref,suplAnnotFile=suplAnnotFile,normalizeMeth=normalizeMeth, sdrf=sdrf,paFi=paFi )}

    ### GROUPING OF REPLICATES AND SAMPLE META-DATA
    ## prepare for sdrf (search in directory above)
    if(isTRUE(sdrf)) {
      hiDir <- dir(file.path(dirname(paFi),".."))
      chFa <- grep("^sdrf.+\\.tsv$", hiDir)
      if(length(chFa) >0) sdrf <- file.path(dirname(paFi),"..",hiDir[chFa[1]]) else {sdrf <- NULL
        if(!silent) message(fxNa,"NO sdrf file found in directory above main data !")}
    }

    if(length(suplAnnotFile) >0 || length(sdrf) >0) {
      headAbund <- utils::head(quant)
      chX <- grepl("^X[[:digit:]]",colnames(quant))                     #check for heading X in all colnames
      if(any(chX)) colnames(headAbund)[which(chX)] <- sub("^X", "", colnames(headAbund)[which(chX)])
      if(length(sampleNames) %in% c(1, ncol(abund))) groupPref$sampleNames <- sampleNames
      if(length(gr) %in% c(1, ncol(abund))) groupPref$gr <- gr
    ## check for matching : (as done within readSampleMetaData) - can't , sdrf not read yet ...
      setupSd <- readSampleMetaData(sdrf=sdrf, suplAnnotFile=suplAnnotFile, quantMeth="AP", path=path, abund=utils::head(quant), chUnit=isTRUE(groupPref$chUnit), groupPref=groupPref, silent=silent, debug=debug, callFrom=fxNa)
    }
    if(debug) {message(fxNa,"rAP13 .."); rAP13 <- list(sdrf=sdrf,gr=gr,suplAnnotFile=suplAnnotFile,abund=abund, quant=quant,refLi=refLi,annot=annot,setupSd=setupSd,sampleNames=sampleNames)}

    ## finish groups of replicates & annotation setupSd
    setupSd <- .checkSetupGroups(abund=abund, setupSd=setupSd, gr=gr, sampleNames=sampleNames, quantMeth="AP", silent=silent, debug=debug, callFrom=fxNa)
    if(debug) {message(fxNa,"rAP13a .."); rAP13a <- list(sdrf=sdrf,gr=gr,suplAnnotFile=suplAnnotFile,abund=abund, quant=quant,refLi=refLi,annot=annot,setupSd=setupSd,sampleNames=sampleNames)}

    ## harmonize sample-names/1
    colNa <- if(length(setupSd$sampleNames)==ncol(abund)) setupSd$sampleNames else setupSd$groups

    ## option : set order of samples as sdrf
    if("sdrfOrder" %in% names(sdrf) && isTRUE(as.logical(sdrf["sdrfOrder"])) && length(setupSd$iniSdrfOrder)==ncol(abund) && ncol(abund) >1) {  # set order according to sdrf (only if >1 samples)
      nOrd <- order(setupSd$iniSdrfOrder)
      ## rename columns according to sdrf and set order of quant and abund ..
      abund <- abund[,nOrd]
      if(length(quant) >0) quant <- quant[,nOrd]
      if(length(setupSd$sampleNames)==ncol(quant)) {
        colNa <- colnames(abund) <- setupSd$sampleNames <- setupSd$sampleNaSdrf[nOrd]  #old# setupSd$sampleNames[nOrd]         ## take sample names from sdrf via  setupSd$sampleNaSdrf
        if(length(quant) >0) colnames(quant) <- setupSd$sampleNaSdrf[nOrd]  #old# setupSd$sampleNames[nOrd]
      } else colNa <- colnames(abund)
      ## now adapt order of setupSd, incl init Sdrf
      if(length(setupSd) >0) { 
        is2dim <- sapply(setupSd, function(x,le) length(dim(x))==2 && nrow(x)==le, le=length(nOrd))    # look for matr or df to change order of lines
        if(any(is2dim) >0) for(i in which(is2dim)) setupSd[[i]] <- setupSd[[i]][nOrd,]
        isVe <- sapply(setupSd, function(x,le) length(x)==le && length(dim(x)) <1, le=length(nOrd))    # look for vector to change order in setupSd
        if(any(isVe) >0) for(i in which(isVe)) setupSd[[i]] <- setupSd[[i]][nOrd] }
      gr <- gr[nOrd]

      if(length(counts) >0 && length(dim(counts))==3) counts <- array(counts[,nOrd,], dim=c(nrow(counts), length(nOrd), dim(counts)[3]), 
        dimnames=list(rownames(counts), colnames(counts)[nOrd], dimnames(counts)[[3]]))
      ## try re-adjusting levels
      tm1 <- sub("^[[:alpha:]]+( |_|-|\\.)+[[:alpha:]]+","", colnames(abund))  # remove heading text
      if(all(grepl("^[[:digit:]]", tm1))) {
        tm1 <- try(as.numeric(sub("( |_|-|\\.)*[[:alpha:]].*","", tm1)), silent=TRUE)   # remove tailing text and try converting to numeric
        if(!inherits(tm1, "try-error")) {
          setupSd$level <- match(tm1, sort(unique(tm1)))
          names(setupSd$level) <- tm1
          if(!silent) message(fxNa,"Sucessfully re-adjusted levels after bringing in order of Sdrf")}
      }     
    } else {

      ## harmonize sample-names/2
      colNa <- colnames(abund)
      chGr <- grepl("^X[[:digit:]]", colNa)                                                # check & remove heading 'X' from initial column-names starting with digits
      if(any(chGr)) colNa[which(chGr)] <- sub("^X","", colNa[which(chGr)])                 #
      if(length(dim(abund)) >1) colnames(abund) <- colNa else if(debug) message(fxNa,"Note : No abundance data found rAP13b")
      if(length(dim(quant)) >1) colnames(quant) <- colNa
    }  
    if(length(setupSd$sampleNames)==ncol(abund)) setupSd$sampleNames <- colNa #no#else setupSd$groups <- colNa
    if(length(dim(counts)) >1 && length(counts) >0) colnames(counts) <- colNa

    if(debug) {message(fxNa,"Read sample-meta data, rAP14"); rAP14 <- list(sdrf=sdrf,suplAnnotFile=suplAnnotFile,abund=abund, quant=quant,refLi=refLi,annot=annot,setupSd=setupSd,plotGraph=plotGraph,normalizeMeth=normalizeMeth,isLog2=isLog2)}

    ## main plotting of distribution of intensities
    custLay <- NULL
    if(is.numeric(plotGraph) && length(plotGraph) >0) {custLay <- as.integer(plotGraph); plotGraph <- TRUE} else {
        if(!isTRUE(plotGraph)) plotGraph <- FALSE}
    if(debug) message(fxNa," rAP15   normalizeMeth= ",normalizeMeth," ;  plotGraph ", plotGraph)
    ## need to plot 2 diustribution3s ?
    ## a) data are same ?
    chSame <- (identical(abund, quant) || identical(log2(abund), quant)) && "none" %in% normalizeMeth
    if(plotGraph) .plotQuantDistr(abund=abund, quant=if(chSame) NULL else quant, custLay=custLay, normalizeMeth=normalizeMeth, notLogAbund=TRUE, softNa=paste("AlphaPept"),
      refLi=refLi, refLiIni=refLiIni, tit=titGraph, silent=debug, callFrom=fxNa, debug=debug)
    
    ## meta-data
    notes <- c(inpFile=paFi, qmethod=paste("AlphaPept"), qMethVersion=if(length(infoDat) >0) unique(infoDat$Software.Revision) else NA,
    	identType="protein", rawFilePath= if(length(infoDat) >0) infoDat$File.Name[1] else NA, normalizeMeth=normalizeMeth, call=deparse(match.call()),
      created=as.character(Sys.time()), wrProteo.version=paste(utils::packageVersion("wrProteo"), collapse="."), machine=Sys.info()["nodename"])
    ## final output
    if(isTRUE(separateAnnot)) list(raw=abund, quant=quant, annot=annot, counts=counts, sampleSetup=setupSd, quantNotes=parametersD, notes=notes) else data.frame(quant,annot) }
  }  
}
    
