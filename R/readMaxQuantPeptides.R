#' Read Peptide Identification and Quantitation Data-Files (peptides.txt) Produced By MaxQuant
#'
#' Peptide level identification and quantification data produced by \href{https://www.maxquant.org/}{MaxQuant} can be read using
#' this function and relevant information extracted.
#' Input files compressed as .gz can be read as well.
#' The peptide abundance values (XIC), peptide counting information and sample-annotation (if available) can be extracted, too.
#'
#'
#' The peptide annotation data gets parsed to extract specific fields (ID, name, description, species ...).
#' Besides, a graphical display of the distribution of peptide abundance values may be generated before and after normalization.
#'
#' @details
#' \href{https://www.maxquant.org/}{MaxQuant} is proteomics quantification software provided by the MaxPlanck institute.
#' By default MaxQuant write the results of each run to the path \code{combined/txt}, from there (only) the files
#'  'peptides.txt' (main quantitation at peptide level), 'summary.txt' and 'parameters.txt' will be used for this function.
#'
#' Meta-data describing the samples and experimental setup may be available from two sources :
#' a) The file \code{summary.txt} which gets produced by MaxQuant in the same folder as the main quantification data.
#' b) Furthermore, meta-data deposited as \code{sdrf} at Pride can be retreived (via the respective github page) when giving
#' the accession number in argument \code{sdrf}.
#' Then, the meta-data will be examined for determining groups of replicates and
#' the results thereof can be found in $sampleSetup$levels.
#' Alternatively, a dataframe formatted like sdrf-files (ie for each sample a separate line, see also function \code{readSdrf}) may be given.
#' In tricky cases it is also possible to precise the column-name to use for defining the groups of replicates or the method for automatically choosing
#'  the most suited column via the 2nd value of the argument \code{sdrf}, see also the function \code{defineSamples} (which gets used internally).
#' Please note, that sdrf is still experimental and only a small fraction of proteomics-data on Pride have been annotated accordingly.
#' If a valid sdrf is furnished, it's information has priority over the information extracted from the MaxQuant produced file summary.txt.
#'
#' This function has been developed using MaxQuant versions 1.6.10.x to 2.0.x, the format of the resulting file 'peptides.txt'
#' is typically well conserved between versions.
#' The final output is a list containing these elements: \code{$raw}, \code{$quant}, \code{$annot}, \code{$counts}, \code{$sampleSetup},
#' \code{$quantNotes}, \code{$notes}, or (if \code{separateAnnot=FALSE}) data.frame
#'   with annotation- and main quantification-content. If \code{sdrf} information has been found, an add-tional list-element \code{setup}
#' will be added containg the entire meta-data as \code{setup$meta} and the suggested organization as \code{setup$lev}.
#'
#'
#' @param path (character) path of file to be read
#' @param fileName (character) name of file to be read (default 'peptides.txt' as typically generated by MaxQuant in txt folder). Gz-compressed files can be read, too.
#' @param normalizeMeth (character) normalization method (for details see \code{\link[wrMisc]{normalizeThis}})
#' @param quantCol (character or integer) exact col-names, or if length=1 content of \code{quantCol} will be used as pattern to search among column-names for $quant using \code{grep}
#' @param contamCol (character or integer, length=1) which columns should be used for contaminants
#' @param pepCountCol (character) pattern to search among column-names for count data (defaults to 'Experiment')
#' @param sampleNames (character) custom column-names for quantification data; this argument has priority over \code{suplAnnotFile}
#' @param extrColNames (character) column names to be read (1st position: prefix for quantitation, default 'intensity';
#'   2nd: column name for peptide-IDs, default )
#' @param specPref (character) prefix to identifiers allowing to separate i) recognize contamination database,
#'   ii) species of main identifications and iii) spike-in species
#' @param refLi (character or integer) custom specify which line of data should be used for normalization, ie which line is main species; if character (eg 'mainSpe'), the column 'SpecType' in $annot will be searched for exact match of the (single) term given
#' @param remRev (logical) option to remove all peptide-identifications based on reverse-peptides
#' @param remConta (logical) option to remove all peptides identified as contaminants
#' @param separateAnnot (logical) if \code{TRUE} output will be organized as list with \code{$annot}, \code{$abund}
#'   for initial/raw abundance values and \code{$quant} with final normalized quantitations
#' @param gr (character or factor) custom defined pattern of replicate association, will override final grouping of
#'   replicates from \code{sdrf} and/or \code{suplAnnotFile} (if provided)   \code{}
#' @param sdrf (character, list or data.frame) optional extraction and adding of experimenal meta-data: if character, this may be the ID at ProteomeExchange,
#'   the second & third elements may give futher indicatations for automatic organization of groups of replicates.
#'   Besides, the output from \code{readSdrf} or a list from \code{defineSamples} may be provided;
#'   if \code{gr} is provided, \code{gr} gets priority for grouping of replicates;
#'   if \code{sdrfOrder=TRUE} the output will be put in order of sdrf
#' @param suplAnnotFile (logical or character) optional reading of supplemental files produced by MaxQuant; if \code{gr} is provided, it gets priority for grouping of replicates
#'   if \code{TRUE} default to files 'summary.txt' (needed to match information of \code{sdrf}) and 'parameters.txt' which can be found in the same folder as the main quantitation results;
#'   if \code{character} the respective file-names (relative ro absolute path), 1st is expected to correspond to 'summary.txt' (tabulated text, the samples as given to MaxQuant) and 2nd to 'parameters.txt' (tabulated text, all parameters given to MaxQuant)
#' @param groupPref (list) additional parameters for interpreting meta-data to identify structure of groups (replicates), will be passed to \code{readSampleMetaData}.
#'   May contain \code{lowNumberOfGroups=FALSE} for automatically choosing a rather elevated number of groups if possible (defaults to low number of groups, ie higher number of samples per group)
#'   May contain \code{chUnit} (logical or character) to be passed to \code{readSampleMetaData()} for (optional) adjustig of unit-prefixes in meta-data group labels, in case multiple different unit-prefixes 
#'   are used (eg '100pMol' and '1nMol').
#' @param plotGraph (logical) optional plot vioplot of initial and normalized data (using \code{normalizeMeth}); alternatively the argument may contain numeric details that will be passed to \code{layout} when plotting
#' @param titGraph (character) custom title to plot
#' @param wex (numeric)  relative expansion factor of the violin in plot
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This function returns a list with  \code{$raw} (initial/raw abundance values), \code{$quant} with final normalized quantitations, \code{$annot} (columns ), \code{$counts} an array with 'PSM' and 'NoOfRazorPeptides',
#'   \code{$quantNotes}, \code{$notes} and optional \code{setup} for meta-data from \code{sdrf}; or a data.frame with quantitation and annotation if \code{separateAnnot=FALSE}
#' @seealso \code{\link[utils]{read.table}}, \code{\link[wrMisc]{normalizeThis}}), for reading protein level \code{\link{readMaxQuantFile}}, \code{\link{readProlineFile}}
#' @examples
#' # Here we'll load a short/trimmed example file (thus not the MaxQuant default name)
#' MQpepFi1 <- "peptides_tinyMQ.txt.gz"
#' path1 <- system.file("extdata", package="wrProteo")
#' specPref1 <- c(conta="conta|CON_|LYSC_CHICK", mainSpecies="YEAST", spec2="HUMAN")
#' dataMQpep <- readMaxQuantPeptides(path1, file=MQpepFi1, specPref=specPref1,
#'   tit="Tiny MaxQuant Peptides")
#' summary(dataMQpep$quant)
#' @export
readMaxQuantPeptides <- function(path, fileName="peptides.txt", normalizeMeth="median", quantCol="Intensity", contamCol="Potential.contaminant",
  pepCountCol="Experiment", refLi=NULL,  sampleNames=NULL,
  extrColNames=c("Sequence","Proteins","Leading.razor.protein","Start.position","End.position","Mass","Missed.cleavages","Unique..Groups.","Unique..Proteins.","Charges"),
  specPref=c(conta="conta|CON_|LYSC_CHICK", mainSpecies="HUMAN"),
  remRev=TRUE, remConta=FALSE, separateAnnot=TRUE, gr=NULL, sdrf=NULL, suplAnnotFile=NULL, groupPref=list(lowNumberOfGroups=TRUE, chUnit=TRUE),
  titGraph=NULL, wex=1.6, plotGraph=TRUE, silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## prepare
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readMaxQuantPeptides")
  opar <- graphics::par(no.readonly=TRUE)
  remStrainNo <- TRUE                   # if TRUE extract Species in very stringent pattern
  cleanDescription <- TRUE              # clean 'Description' for artifacts of truncated text (tailing ';' etc)
  setupSd <- NULL                   # initialize

  oparMar <- graphics::par("mar")                     # old margins, for rest after figure
  oparLayout <- graphics::par("mfcol")                # old layout, for rest after figure
  on.exit(graphics::par(mar=oparMar, mfcol=oparLayout))            # restore old mar settings
  reqPa <- c("utils","wrMisc")
  chPa <- sapply(reqPa, requireNamespace, quietly=TRUE)
  if(any(!chPa)) stop("package(s) '",paste(reqPa[which(!chPa)], collapse="','"),"' not found ! Please install first from CRAN")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  excluCol <- "^Abundances.Count"   # exclude this from quantifications columns
  setupSd <- setupSdmq <- summaryD <- parametersD <- NULL       # initialize (meta-data)
  if(debug) {message("rMQP1"); rMQP1 <- list(path=path,chPa=chPa,setupSd=setupSd,excluCol=excluCol,remStrainNo=remStrainNo,cleanDescription=cleanDescription)}

  ## check if path & file exist
  if(!grepl("\\.txt$|\\.txt\\.gz$",fileName)) message(fxNa,"Suspicious filename, this function was designed for reading tabulated text files produced by MaxQuant")
  paFi <- wrMisc::checkFilePath(fileName, path, expectExt="txt", compressedOption=TRUE, stopIfNothing=TRUE, callFrom=fxNa, silent=silent,debug=debug)  
  if(debug) {message("rMQP1d")}

  chPa <- try(find.package("utils"), silent=TRUE)
  if(inherits(chPa, "try-error")) stop("package 'utils' not found ! Please install first from CRAN")

  ## initial read MaxQuant
  tmp <- try(utils::read.delim(file.path(path,fileName), stringsAsFactors=FALSE), silent=TRUE)
  if(debug) {message(fxNa,"rMQP2"); rMQP2 <- list()}

  if(length(tmp) <1 | inherits(tmp, "try-error")) {
    message(fxNa,"Unable to read file 'fileName' !  Returning NULL; check if this is original MaxQuant-file and/or rights to read")
    tmp <- NULL
    return(NULL)
  } else {
    if(any(extrColNames[2:3] %in% colnames(tmp), na.rm=TRUE)) {
      ## check which columns can be extracted (for annotation)
      if(is.integer(contamCol)) {contamCol <- colnames(tmp)[contamCol]
        if(debug) message(fxNa," Custom 'contamCol' points to ",contamCol)}
      extrColNames <- wrMisc::naOmit(union(extrColNames, contamCol))                 # add contamCol if not included in extrColNames
      chCol <- extrColNames %in% colnames(tmp)
      if(any(!chCol, na.rm=TRUE)) { extrColNames[which(!chCol)] <- gsub("\\."," ",extrColNames[which(!chCol)])
        chCol <- extrColNames %in% colnames(tmp) }
      if(all(!chCol)) stop("Problem locating annotation columns (",wrMisc::pasteC(extrColNames, quoteC="''"),")")
      if(any(!chCol, na.rm=TRUE) ) {
        if(!silent) message(fxNa,"Note: Can't find columns ",wrMisc::pasteC(extrColNames[!chCol], quoteC="'")," !")
      }
      if(debug) {message(fxNa,"rMQP3"); rMQP3 <- list(path=path,chPa=chPa,tmp=tmp,extrColNames=extrColNames,chCol=chCol,remConta=remConta,quantCol=quantCol)}

      ## 'REVERSE' peptides
      chRazProCol <- extrColNames[3] %in% colnames(tmp)
      chRev <- grep("REV__", tmp[,extrColNames[if(chRazProCol) 3 else 2]])
        if(debug) message(fxNa,"rMQP3b    chRev : ",utils::head(paste(chRev, collapse=" ")),"\n")
      if(length(chRev) >0) {
        if(!silent) message(fxNa,"Note: Found ",length(chRev)," out of ",nrow(tmp)," peptides with proteins marked as 'REV_' (reverse peptide identification)", if(isTRUE(remRev)) " - Removing")
        if(isTRUE(remRev)) tmp <- if(length(chRev) < nrow(tmp) -1)  tmp[-1*chRev,] else matrix(tmp[-1*chRev,], nrow=nrow(tmp)-length(remRev), dimnames=list(rownames(tmp)[-1*chRev], colnames(tmp)))
      }
      ## remove MaxQuant internal contaminants CON__
      if(isTRUE(remConta) && nrow(tmp) >0) { isConta <- grep("CON__{0,1}[[:alpha:]]+", tmp[,extrColNames[2]])
        if(length(isConta) >0) {
          if(!silent) message(fxNa,"Note: Found ",length(isConta)," out of ",nrow(tmp)," proteins marked as 'CON_' (contaminants) - Removing")
          tmp <- if(length(isConta) < nrow(tmp) -1) tmp[-1*isConta,] else matrix(tmp[-1*isConta,], nrow=nrow(tmp)-length(isConta), dimnames=list(rownames(tmp)[-1*isConta], colnames(tmp)))
        } }
      } else if(!silent) message(fxNa,"BIZZARE, trouble ahead : Unable to find columns ",wrMisc::pasteC(extrColNames[2:3],quoteC="'")," (from argument 'extrColNames')")
    if(debug) {message(fxNa,"rMQP4"); rMQP4 <- list(path=path,chPa=chPa,tmp=tmp,extrColNames=extrColNames,chCol=chCol,chRazProCol=chRazProCol,chRev=chRev,remConta=remConta,quantCol=quantCol,pepCountCol=pepCountCol,specPref=specPref)}
  }
  if(length(tmp) >0) {
    ## further extracting : quantitation
    grepX <- function(x) grep(paste0(x,"\\."), colnames(tmp))
    useDCol <- if(length(quantCol)==1) grepX(quantCol) else unique(as.integer(sapply(quantCol, grepX)))
    if(length(useDCol) <1) { useDCol <- if(length(quantCol)==1) grepX(tolower(quantCol)) else unique(as.integer(sapply(tolower(quantCol), grepX)))
      if(length(useDCol) >0) {  # need to re-adjust 'quantCol'
        ch1 <- regexpr(tolower(quantCol), tolower(colnames(tmp)[useDCol])) 
        quantCol <- unique(substr(colnames(tmp)[useDCol], ch1, ch1+ nchar(quantCol)-1) )
        if(length(quantCol) >1) {if(!silent) message(fxNa,"Need to adjust 'quantCol' , but found ",length(quantCol)," variants, using '",quantCol[1],"' ")}
      } 
    }
   if(debug) {message(fxNa,"rMQP4a"); rMQP4a <- list()}
    if(length(useDCol) <1) { abund <- NULL
      warning(fxNa, "NO columns matching terms ",wrMisc::pasteC(quantCol, quoteC="'")," from argument 'quantCol' found !  rMQP4a")
    } else {
      quantColP <- NULL                           # initialize
      if(length(quantCol) <1) stop(" 'quantCol' must be provided !")
      if(length(quantCol) >1) { abund <- as.matrix(wrMisc::extrColsDeX(tmp, extrCol=quantCol, doExtractCols=TRUE, silent=silent, callFrom=fxNa))
      } else { chP <- substr(quantCol, nchar(quantCol), nchar(quantCol)) != "."
        #quantColP <- quantCol
        quantCol <- if(chP) grep(paste0(quantCol,"\\."), colnames(tmp)) else grep(quantCol, colnames(tmp))
        chNa <- is.na(quantCol)
        if(all(chNa) | length(quantCol) <1) stop("Could not find any of the columns specified in argument 'quantCol' !")
        if(any(chNa)) {
          if(!silent) message(fxNa,"Could not find columns ",wrMisc::pasteC(quantCol[which(chNa)],quote="'")," .. omit")
          quantCol <- wrMisc::naOmit(quantCol)}
        abund <- as.matrix(tmp[,quantCol]) }           # abundance val
        if(debug) {message(fxNa,"rMQP4b"); rMQP4b <- list(path=path,chPa=chPa,tmp=tmp,abund=abund,extrColNames=extrColNames,chCol=chCol,chRazProCol=chRazProCol,chRev=chRev,remConta=remConta,quantCol=quantCol,pepCountCol=pepCountCol,specPref=specPref)}
      chNum <- is.numeric(abund)
      if(!chNum) {abund <- apply(tmp[,quantCol], 2, wrMisc::convToNum, convert="allChar", silent=silent, callFrom=fxNa)}
      if(length(dim(abund)) <2) abund <- matrix(as.numeric(abund), ncol=ncol(abund), dimnames=dimnames(abund))
      colnames(abund) <- if(length(quantColP)==1) sub(paste0(quantColP,"\\."),"", colnames(abund)) else wrMisc::.trimLeft(wrMisc::.trimRight(colnames(abund)))
      if(debug) {message(fxNa,"rMQP5")}

      ## convert 0 to NA
      ch1 <- abund <= 0
      if(any(ch1, na.rm=TRUE)) { abund[which(ch1)] <- NA
        if(!silent) message(fxNa,"Transform ",sum(ch1),"(",100*round(sum(ch1)/length(ch1),3),"%) initial '0' values to 'NA'")}
    }

    ## further extracting : prepare for countig data
    ch1 <- !grepl("\\^",pepCountCol)
    if(any(ch1, na.rm=TRUE)) {pepCountCol[which(ch1)] <- paste0("^",pepCountCol[which(ch1)])}  # add heading '^' (if not yet present)
    ch1 <- !grepl(" $",pepCountCol)
    if(any(ch1, na.rm=TRUE)) {pepCountCol[which(ch1)] <- paste0(pepCountCol[which(ch1)]," ")}  # add tailing ' ' (if not yet present)
    if(length(grep("\\\\",pepCountCol)) <1) pepCountCol <- gsub("\\.","\\\\.",pepCountCol)       # protect '.' (if not yet protected)
    ## prepare for column-name style with '.' or '...'
    tm2 <- lapply(as.list(pepCountCol), function(x) c(x, gsub(" ",".", sub(" \\+ ","...",x))) )
    names(tm2) <- pepCountCol
    usePCol <- lapply(tm2, function(x) {ch1 <- lapply(x, grep, colnames(tmp)); if(length(ch1) >1) ch1[[which.max(sapply(ch1,length))]] else ch1[[1]]})
    usePCol <- lapply(usePCol, wrMisc::naOmit)
    ch2 <- sapply(usePCol, length)  -ncol(abund)                # take abund as ref for not extracting more cols
    if(any(ch2 >0, na.rm=TRUE)) usePCol[which(ch2 >0)] <- lapply(usePCol[which(ch2 >0)], function(x) x[-1])
    ch2 <- sapply(usePCol, length) ==ncol(abund)      # update
    if(!silent && any(!ch2, na.rm=TRUE)) message(fxNa,"Could not find peptide counts columns (argument 'pepCountCol') matching to '",pepCountCol[which(!ch2)],"'")
    if(debug) {message(fxNa,"rMQP6"); rMQP6 <- list(path=path,chPa=chPa,tmp=tmp,extrColNames=extrColNames,chCol=chCol,chRazProCol=chRazProCol,chRev=chRev,quantCol=quantCol,abund=abund,chNum=chNum,ch2=ch2,usePCol=usePCol,pepCountCol=pepCountCol,specPref=specPref,remConta=remConta)}
    ## make array of PSM counts etc
    if(any(ch2, na.rm=TRUE)) {
      counts <- array(dim=c(nrow(tmp),ncol(abund),sum(ch2)), dimnames=list(NULL, colnames(abund), pepCountCol[which(ch2)]))
      for(i in 1:sum(ch2)) counts[,,i] <- as.numeric(as.matrix(tmp[,usePCol[[which(ch2)[i]]] ]))
    } else counts <- NULL
    if(debug) {message(fxNa,"rMQP6a")}

    ## Annotation
    useACol <- list(annC=match(extrColNames, colnames(tmp)) )
    annot <- as.matrix(tmp[,useACol$annC])
    specMQ <- rep(NA, nrow(abund))         # initialize
    if(debug) {message(fxNa,"rMQP6b")}

    useProCo <- 2 + extrColNames[3] %in% colnames(tmp)         # specific to eptide reading
    .MultGrep <- function(pat, y) if(length(pat)==1) grep(pat, y) else unlist(sapply(pat, grep, y))  # (multiple) grep() when length of pattern 'pat' >0

    ## MaxQuant internal contaminants specific : remove non-protein DB parts - if possible, eg "CON__ENSEMBL:ENSBTAP00000007350;CON__P01030" -> "CON__P01030"
    conID <- paste0("CON__",c("ENSEMBL","REFSEQ","H-INV"),":")
    ch2 <- sapply(sapply(conID, grep, annot[,useProCo]), length) >0
    if(any(ch2, na.rm=TRUE)) {
      conID <- conID[which(ch2)]
      conID <- paste0(conID, c("[[:upper:]]+[[:digit:]]*;{0,1}", "[[:upper:]]+_[[:digit:]]+;{0,1}", "[[:upper:]]+[[:digit:]]+;{0,1}"))
      acc1 <- annot[,useProCo]
      for(i in 1:length(conID)) {
        acc2 <- acc1                     # need previous 'status' to compare if all text disappears
        acc1 <- sub(conID[i], "", acc1)
        chLe <- nchar(acc1)  <2
        if(any(chLe, na.rm=TRUE)) acc1[which(chLe)] <- sub("CON__","", acc2[which(chLe)]) }  # remove entire entry only if something (else) remains
      ## remove first of CON_ entries (only if min 3 characters 3 remain)
      ch2 <- grep("CON__{0,1}[A-Z0-9]+;", acc1)
      if(length(ch2) >0) { acc2 <- acc1
        acc1 <- sub("CON__{0,1}[A-Z0-9]+;", "", acc1)
        chLe <- nchar(acc1) <2
        if(any(chLe, na.rm=TRUE)) acc1[which(chLe)] <- sub("CON__","", acc2[which(chLe)]) }  # remove entire entry only if something (else) remains
      ## remove first of "CON_" marks
      ch2 <- grep("CON_", acc1)
      if(length(ch2) >0) acc1[ch2] <- sub("CON__{0,1}","", acc1[ch2])
      annot[,useProCo] <- acc1 }

    if(length(specPref) >0) {
      ## look if available, for specif tags (otherwise look in 'Proteins')
      specMQ0 <- lapply(specPref, .MultGrep, annot[,useProCo])    # in 'Proteins'
      for(i in 1:length(specMQ0)) {if(length(specMQ0[[i]]) >0) specMQ[as.integer(specMQ0[[i]])] <- names(specMQ0)[i]}
    }
    if(debug) {message(fxNa,"rMQP6c")}
    annot <- cbind(SpecType=specMQ, annot)                                       # better to name column 'species' ??
    if(debug) {message(fxNa,"rMQP7"); rMQP7 <- list(path=path,chPa=chPa,tmp=tmp,extrColNames=extrColNames,chCol=chCol,chRazProCol=chRazProCol,counts=counts,specPref=specPref,
      chRev=chRev,quantCol=quantCol,abund=abund,chNum=chNum,ch2=ch2,annot=annot,remConta=remConta,specPref=specPref)}

    ## remove MQ-internal contaminants
    if(remConta & extrColNames[useProCo] %in% colnames(annot)) {
      conLi <- grep("CON__[[:alnum:]]", annot[,extrColNames[2]])
      if(length(conLi) >0) {
        iniLi <- nrow(annot)
        annot <- annot[-conLi,]
        abund <- abund[-conLi,]
        specMQ <- specMQ[-conLi]       # needed ??
        #specMQ0 <- specMQ0[-conLi]
        counts <- if(length(dim(counts))==3) counts[-conLi,,] else counts[-conLi,,]
        if(debug) message(fxNa,"Removing ",length(conLi)," instances of MaxQuant-contaminants to final ",nrow(annot)," lines/IDs")} }
    if(debug) {message(fxNa,"rMQP7b"); rMQP7b <- list()}

    ## split Annotation
    remHeader <- c("^conta\\|","^sp\\|")
    MQan2 <- strsplit(sub(remHeader[1], "", sub(remHeader[2], "", if(useProCo==2) sub(";.+", "", annot[,useProCo+1]) else annot[,useProCo+1])), "\\|") # separate AccessionNumber (eg P02768) and EntryName (eg ALBU_HUMAN)       
    MQan2 <- t(sapply(MQan2, function(x) if(length(x)==1) { c(NA,x)} else x[1:2]) )
    colnames(MQan2) <- c("Accession","EntryName")
    MQan2 <- cbind(MQan2, Species=NA)
    hasSpe <- grep("._", MQan2[,2])
    if(length(hasSpe) >0) MQan2[hasSpe,3] <- gsub(".*_", "", MQan2[hasSpe,2])       # separate 'YEAST' from "MTNA_YEAST
    
    ## Problem : "CHICK" stays in annot[,"Species"]
    if("Species" %in% colnames(annot)) { specNa <- unique(wrMisc::naOmit(annot[,"Species"]))
      ch2 <- grepl("$[[:upper:]]+$", specNa)
      if(any(ch2, na.rm=TRUE)) for(i in specNa[which(ch2)]) {annot[which(annot[,"Species"] %in% i),"Species"] <- inspectSpeciesIndic(i)}
    }

    #gsub(".*_", "", MQan2[,2]))   # separate AccessionNumber (eg P02768) and EntryName (eg ALBU_HUMAN)
    #MQan2 <- cbind(MQan2, Species=sub("_.+|[[:punct:]].+","", sub("[[:upper:]]+[[:digit:]]*_", "", MQan2[,2]))) # separate AccessionNumber (eg P02768) and EntryName (eg ALBU_HUMAN)
    if(debug) {message(fxNa,"rMQP7c"); rMQP7c <- list(specPref=specPref,annot=annot,path=path,MQan2=MQan2,chPa=chPa,tmp=tmp,extrColNames=extrColNames,chCol=chCol,chRazProCol=chRazProCol,counts=counts,remStrainNo=remStrainNo,specPref=specPref)}



    
    ##  prepare for reading fasta
    ## look for tags from  specPref
    if(length(specPref) >0) {
      ## set annot[,"specPref"] according to specPref
      annot <- .extrSpecPref(specPref, annot, useColumn=c("Leading.razor.protein","Proteins","Majority.protein.IDs","Fasta.headers"), silent=silent, debug=debug, callFrom=fxNa)  # useful ,"Majority.protein.IDs","Fasta.headers"
    } else if(debug) message(fxNa,"Note: Argument 'specPref' not specifed (empty)")
    if(debug) { message(fxNa,"rMQP7c2"); rMQP7c2 <- list(path=path,chPa=chPa,tmp=tmp,abund=abund,specPref=specPref,MQan2=MQan2,annot=annot,specMQ0=specMQ0,remStrainNo=remStrainNo,remConta=remConta,extrColNames=extrColNames,remHeader=remHeader)}
    
    
    


    ## extract species according to custom search parameters 'specPref'
    .annSpecies <- function(spe=c("_HUMAN","Homo sapiens"), anno=annot, exCoNa=extrColNames) {
      ## extract species tags out of annot[,exCoNa[2]], place as convert to regular name in anno, return matrix anno
      ch1 <- grep(spe[1], anno[,exCoNa[2]])
      if(length(ch1) >0) anno[ch1,"Species"] <- spe[2]  #"Homo sapiens"
      anno }
    if(remStrainNo) {
      commonSpec <- .commonSpecies()
      commonSpec[,1] <- sub("^_","", commonSpec[,1])   # '_' has already been stripped off during strsplit
      convSpe <- which(commonSpec[,1] %in% unique(MQan2[hasSpe,"Species"]))
      if(length(convSpe) >0) for(i in convSpe) MQan2[hasSpe,"Species"] <- sub(commonSpec[i,1], commonSpec[i,2], MQan2[hasSpe,"Species"]) }

    annot <- cbind(annot, MQan2)
    if(debug) { message(fxNa,"rMQP7d"); rMQP7d <- list(path=path,chPa=chPa,tmp=tmp,abund=abund,specPref=specPref,MQan2=MQan2,annot=annot,specMQ0=specMQ0,remStrainNo=remStrainNo,remConta=remConta,extrColNames=extrColNames,remHeader=remHeader)}

    ## contaminants (fuse from column 'Potential.contaminant' and those found via specPref[1])
    contam <- rep(FALSE, nrow(annot))
    if(!all(is.na(specMQ0))) if(length(specMQ0$conta) >0) contam[specMQ0$conta] <- TRUE         ## from 'specPref' search
    if("Potential.contaminant" %in% colnames(annot)) { chCo <- grepl("+",annot[,"Potential.contaminant"])
      if(any(chCo, na.rm=TRUE)) contam[which(chCo)[1]] <- TRUE
      annot[,"Potential.contaminant"] <- contam
    } else annot <- cbind(annot, Potential.contaminant=contam)
    if(debug) {message(fxNa,"rMQP9");  rMQP9 <- list(annot=annot,tmp=tmp,abund=abund,MQan2=MQan2,specMQ0=specMQ0,remStrainNo=remStrainNo,remConta=remConta,extrColNames=extrColNames,remHeader=remHeader)}
    
    
    ## mine info to add 'EntryName' and 'Accession' 
    if("Leading.razor.protein" %in% colnames(annot)) {       
      ## rm heading db
      ann3 <- ann2 <- sub("^[[:lower:]]+\\|","", annot[,"Leading.razor.protein"])
      hasSep <- grep("\\|", ann2)      
      if(length(hasSep) >0) { ann3[hasSep] <- sub("[[:upper:]]+([[:digit:]]|[[:upper:]])+\\|", "", ann2[hasSep])
        ann2[hasSep] <- sub("\\|.+","", ann2[hasSep]) } 
      annot <- cbind(Accession=ann2, EntryName=ann3, annot )  
    }
    if(debug) {message(fxNa,"rMQP9b");  rMQP9b <- list(path=path,chPa=chPa,tmp=tmp,extrColNames=extrColNames,chCol=chCol,chRazProCol=chRazProCol,counts=counts,contam=contam,
      quantCol=quantCol,abund=abund,chNum=chNum,ch2=ch2,annot=annot,specMQ=specMQ,remConta=remConta,ann3=ann3,ann2=ann2)}

    ## look for unique col from $annot to use as rownames
    rowNa <- annot[,2]
    chAn <- sum(duplicated(rowNa))
    if(chAn >0) {
      chMod <- grep("^Oxidation", colnames(tmp))     # look for column  "Oxidation..M..site.IDs"
      if(length(chMod) >0) {
        hasMod <- which(nchar(tmp[,chMod[1]]) >0)
        if(length(hasMod) >0) {rowNa[hasMod] <- paste0(rowNa[hasMod],".ox")
          chAn <- sum(duplicated(rowNa))
          if(chAn >0 && !silent) message(fxNa," Still ",chAn, " duplicated rownames")}
    } }
    if(chAn >0) {
      chAn <- duplicated(rowNa, fromLast=FALSE)
      if(any(chAn, na.rm=TRUE)) {
        chAn <- which(chAn | duplicated(rowNa, fromLast=TRUE))
        ## note : this is peptide-based - so it is normal that many proteins appear multiple times (as diff peptides)
        rowNa[chAn] <- wrMisc::correctToUnique(rowNa[chAn]) 
        if(!silent) message(fxNa,"Note : Some peptide sequences appear duplicated (despite considering for oxidation)")
    } }

    rownames(abund) <- rownames(annot) <- rowNa
    if(length(counts) >0) rownames(counts) <- rownames(annot)
    if(debug) {message(fxNa,"rMQP9c");  rMQP9c <- list(path=path,chPa=chPa,tmp=tmp,extrColNames=extrColNames,chCol=chCol,chRazProCol=chRazProCol,counts=counts,contam=contam,
      rowNa=rowNa,refLi=refLi,chRev=chRev,quantCol=quantCol,abund=abund,chNum=chNum,ch2=ch2,annot=annot,specMQ=specMQ,remConta=remConta)}


    ## check for reference for normalization
    refLiIni <- refLi
    if(is.character(refLi) && length(refLi)==1) { refLi <- which(annot[,"SpecType"]==refLi)
      if(length(refLi) <1) message(fxNa,"Could not find any peptide matching argument 'refLi', ignoring ...") else {
        if(!silent) message(fxNa,"Normalize using subset of ",length(refLi)) } }           # may be "mainSpe"
    if(length(refLi) <1) refLi <- NULL


    ## take log2 & normalize
    quant <- try(wrMisc::normalizeThis(log2(abund), method=normalizeMeth, mode="additive", refLines=refLi, silent=silent, debug=debug, callFrom=fxNa), silent=TRUE)
    if(inherits(quant, "try-error")) { warning(fxNa,"PROBLEMS ahead : Unable to normalize as log2-data !!") }

    if(debug) {message(fxNa,"rMQP10"); rMQP10 <- list(path=path,chPa=chPa,tmp=tmp,extrColNames=extrColNames,groupPref=groupPref,chCol=chCol,chRev=chRev,quantCol=quantCol,abund=abund,chNum=chNum,ch2=ch2,
      quant=quant,annot=annot,MQan2=MQan2,contam=contam,remConta=remConta)}

    ### GROUPING OF REPLICATES AND SAMPLE META-DATA
    if(length(suplAnnotFile) >0 || length(sdrf) >0) {
      if(length(sampleNames) %in% c(1, ncol(abund))) groupPref$sampleNames <- sampleNames
      if(length(gr) %in% c(1, ncol(abund))) groupPref$gr <- gr
      if(debug) {message(fxNa,"rMQP11"); rMQP11 <- list()}
      setupSd <- readSampleMetaData(sdrf=sdrf, suplAnnotFile=suplAnnotFile, quantMeth="MQ", path=path, abund=utils::head(quant), chUnit=isTRUE(groupPref$chUnit), groupPref=groupPref, silent=silent, debug=debug, callFrom=fxNa)
    }
    if(debug) {message(fxNa,"rMQP13 .."); rMQP13 <- list(path=path,chPa=chPa,tmp=tmp,extrColNames=extrColNames,groupPref=groupPref,chCol=chCol,chRev=chRev,quantCol=quantCol,abund=abund,chNum=chNum,ch2=ch2,
      quant=quant,annot=annot,MQan2=MQan2,contam=contam,remConta=remConta,setupSd=setupSd)}

    ## finish groups of replicates & annotation setupSd
    setupSd <- .checkSetupGroups(abund=abund, setupSd=setupSd, gr=gr, sampleNames=sampleNames, quantMeth="MQ", silent=silent, debug=debug, callFrom=fxNa)
    if(debug) {message(fxNa,"rMQP13b .."); rMQP13b <- list(path=path,chPa=chPa,tmp=tmp,extrColNames=extrColNames,groupPref=groupPref,chCol=chCol,chRev=chRev,quantCol=quantCol,abund=abund,chNum=chNum,ch2=ch2,
      sdrf=sdrf,quant=quant,annot=annot,MQan2=MQan2,contam=contam,remConta=remConta,setupSd=setupSd)}

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
      if(debug) {message(fxNa,"rMQP13c .."); rMQP13c <- list(path=path,chPa=chPa,tmp=tmp,extrColNames=extrColNames,groupPref=groupPref,chCol=chCol,chRev=chRev,quantCol=quantCol,abund=abund,chNum=chNum,ch2=ch2,
        sdrf=sdrf,quant=quant,annot=annot,MQan2=MQan2,contam=contam,remConta=remConta,setupSd=setupSd)}
        
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

    if(debug) {message(fxNa,"Read sample-meta data, rMQP14"); rMQP14 <- list(sdrf=sdrf,suplAnnotFile=suplAnnotFile,abund=abund, quant=quant,refLi=refLi,annot=annot,setupSd=setupSd)}

    ## main plotting of distribution of intensities
    custLay <- NULL
    if(is.numeric(plotGraph) && length(plotGraph) >0) {custLay <- as.integer(plotGraph); plotGraph <- TRUE} else {
        if(!isTRUE(plotGraph)) plotGraph <- FALSE}
    if(plotGraph) .plotQuantDistr(abund=abund, quant=quant, custLay=custLay, normalizeMeth=normalizeMeth, softNa="MaxQuant Peptides",
      refLi=refLi, refLiIni=refLiIni, tit=titGraph, silent=debug, callFrom=fxNa, debug=debug)


    ## meta-data
    notes <- c(inpFile=file.path(path,fileName), qmethod="MaxQuant", qMethVersion=if(length(parametersD) >0) "xx" else NA,
      identType="peptide", rawFilePath=if(length(parametersD) >0) "xx" else NA, normalizeMeth=normalizeMeth, call=deparse(match.call()), created=as.character(Sys.time()),
      wrProteo.version=paste(utils::packageVersion("wrProteo"), collapse="."), machine=Sys.info()["nodename"])
    ## prepare for final output
    if(isTRUE(separateAnnot)) list(raw=abund, quant=quant, annot=annot, counts=counts, sampleSetup=setupSd, quantNotes=parametersD, notes=notes) else data.frame(abund, annot) }
}
    
                                                                                                       