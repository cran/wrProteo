#' Read Tabulated Files Exported By ProteomeDiscoverer At Protein Level
#'
#' Depreciated old version of Protein identification and quantification results from
#' \href{https://www.thermofisher.com/order/catalog/product/OPTON-30812}{Thermo ProteomeDiscoverer}
#' which were exported as tabulated text can be imported and relevant information extracted.
#' The final output is a list containing 3 elements: \code{$annot}, \code{$raw} and optional \code{$quant},
#' or returns data.frame with entire content of file if \code{separateAnnot=FALSE}.
#' Please use readProteomeDiscovererFile() from the same package instead !
#'
#' @details
#' This function has been replaced by \code{readProteomeDiscovererFile} (from the same package) !
#' The syntax  and strcuture of output has remained the same, you can simply replace the name of the function called.
#' 
#' This function has been developed using Thermo ProteomeDiscoverer versions 2.2 to 2.5.
#' The format of resulting files at export also depends which columns are chosen as visible inside ProteomeDiscoverer and subsequently get chosen for export.
#' Using the argument \code{suplAnnotFile} it is possible to specify a specific file (or search for default file) to read for extracting file-names as sample-names and other experiment realted information.
#' If a column named \code{contamCol} is found, the data will be lateron filtered to remove all contaminants, set to \code{NULL} for keeping all contaminants
#' This function replaces the depreciated function \code{readPDExport}.
#'
#' @param fileName (character) name of file to be read
#' @param path (character) path of file to be read
#' @param normalizeMeth (character) normalization method, defaults to \code{median}, for more details see \code{\link[wrMisc]{normalizeThis}})
#' @param sampleNames (character) custom column-names for quantification data (ProteomeDiscoverer does not automatically use file-names from spectra); this argument has priority over \code{suplAnnotFile}
#' @param read0asNA (logical) decide if initial quntifications at 0 should be transformed to NA
#' @param quantCol (character or integer) exact col-names, or if length=1 content of \code{quantCol} will be used as pattern to search among column-names for $quant using \code{grep}
#' @param contamCol (character or integer, length=1) which columns should be used for contaminants marked by ProteomeDiscoverer.
#'  If a column named \code{contamCol} is found, the data will be lateron filtered to remove all contaminants, set to \code{NULL} for keeping all contaminants
#' @param refLi (character or integer) custom specify which line of data is main species, if character (eg 'mainSpe'), the column 'SpecType' in $annot will be searched for exact match of the (single) term given
#' @param separateAnnot (logical) if \code{TRUE} output will be organized as list with \code{$annot}, \code{$abund} for initial/raw abundance values and \code{$quant} with final log2 (normalized) quantitations
#' @param annotCol (character) column names to be read/extracted for the annotation section (default  c("Accession","Description","Gene","Contaminant","Sum.PEP.Score","Coverage....","X..Peptides","X..PSMs","X..Unique.Peptides", "X..AAs","MW..kDa.") )
#' @param FDRCol (list) optional indication to search for protein FDR information
#' @param specPref (character or list) define characteristic text for recognizing (main) groups of species (1st for comtaminants - will be marked as 'conta', 2nd for main species- marked as 'mainSpe',
#'  and optional following ones for supplemental tags/species - maked as 'species2','species3',...);
#'  if list and list-element has multiple values they will be used for exact matching of accessions (ie 2nd of argument \code{annotCol})
#' @param gr (character or factor) custom defined pattern of replicate association, will override final grouping of replicates from \code{sdrf} and/or \code{suplAnnotFile} (if provided)   \code{}
#' @param sdrf (character, list or data.frame) optional extraction and adding of experimenal meta-data: if character, this may be the ID at ProteomeExchange,
#'   the second element may give futher indicatations for automatic organization of groups of replicates.
#'   Besides, the output from \code{readSdrf} or a list from \code{defineSamples} may be provided; if \code{gr} is provided, \code{gr} gets priority for grouping of replicates
#' @param suplAnnotFile (logical or character) optional reading of supplemental files produced by ProteomeDiscoverer; however, if \code{gr} is provided, \code{gr} gets priority for grouping of replicates;
#'  if \code{TRUE} defaults to file '*InputFiles.txt' (needed to match information of \code{sdrf}) which can be exported next to main quantitation results;
#'  if \code{character} the respective file-name (relative or absolute path)
#' @param groupPref (list) additional parameters for interpreting meta-data to identify structure of groups (replicates), will be passed to \code{readSampleMetaData}.
#'   May contain \code{lowNumberOfGroups=FALSE} for automatically choosing a rather elevated number of groups if possible (defaults to low number of groups, ie higher number of samples per group)
#' @param plotGraph (logical or integer) optional plot of type vioplot of initial and normalized data (using \code{normalizeMeth}); if integer, it will be passed to \code{layout} when plotting
#' @param titGraph (character) custom title to plot of distribution of quantitation values
#' @param wex (integer) relative expansion factor of the violin-plot (will be passed to \code{\link[wrGraph]{vioplotW}})
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a list with \code{$raw} (initial/raw abundance values), \code{$quant} with final normalized quantitations, \code{$annot}, \code{$counts} an array with number of peptides, \code{$quantNotes}
#'  and \code{$notes}; or if \code{separateAnnot=FALSE} the function returns a data.frame with annotation and quantitation only
#' @seealso \code{\link[utils]{read.table}}, \code{\link[wrMisc]{normalizeThis}}) , \code{\link{readMaxQuantFile}}, \code{\link{readProlineFile}}, \code{\link{readFragpipeFile}}
#' @examples
#' path1 <- system.file("extdata", package="wrProteo")
#' fiNa <- "tinyPD_allProteins.txt.gz"
#' ## Please use the function readProteinDiscovererFile(), as shown below (same syntax)
#' dataPD <- readProteomeDiscovererFile(file=fiNa, path=path1, suplAnnotFile=FALSE)
#' summary(dataPD$quant)
#'
#' @export
readProtDiscovFile <- function(fileName, path=NULL, normalizeMeth="median", sampleNames=NULL, read0asNA=TRUE, quantCol="^Abundances*", annotCol=NULL, contamCol="Contaminant",
  refLi=NULL, separateAnnot=TRUE, FDRCol=list(c("^Protein.FDR.Confidence","High"), c("^Found.in.Sample.","High")), gr=NULL, sdrf=NULL, suplAnnotFile=TRUE,
  groupPref=list(lowNumberOfGroups=TRUE), specPref=c(conta="CON_|LYSC_CHICK", mainSpecies="OS=Homo sapiens"), plotGraph=TRUE, wex=1.6, titGraph="Proteome Discoverer",
  silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## read ProteomeDiscoverer exported txt

  message("+++ This function is depreciated, it has been replaced by readProteomeDiscovererFile() from the same package ! \n +++ Synthax and structure of output remain the same ! \n")
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readProtDiscovFile")
  reqPa <- c("utils","wrMisc")
  chPa <- sapply(reqPa, requireNamespace, quietly=TRUE)
  if(any(!chPa)) stop("package(s) '",paste(reqPa[which(!chPa)], collapse="','"),"' not found ! Please install first from CRAN")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  excluCol <- "^Abundances.Count"   # exclude this from quantifications columns
  cleanDescription <- TRUE          # clean 'Description' for artifacts of truncated text (tailing ';' etc)
  infoDat <- infoFi <- setupSd <- parametersD <- NULL        # initialize
  .corPathW <- function(x) gsub("\\\\", "/", x)

  ## check if path & file exist
  msg <- "Invalid entry for 'fileName'"
  if(length(fileName) >1) { fileName <- fileName[1]
    if(!silent) message(fxNa," 'fileName' shoud be of length=1, using 1st value")    # nolint # nolint
  } else { if(length(fileName) <1) stop(msg) else if(nchar(fileName) <0) stop(msg)}
  paFi <- fileName                      # presume (& correct if path is given)
  chFi <- file.exists(fileName)         # presume (& correct otherwise)
  if(length(path) >0) if(!dir.exists(path[1])) { path <- NULL
    if(!silent) message(fxNa,"Invalid path '",path[1],"'  (not existing), ignoring...") }
  if(length(path) >0) { chFi <- file.exists(file.path(path[1], fileName))
    if(chFi) paFi <- file.path(path[1], fileName) else {
      if(file.exists(fileName)) {paFi <- fileName
        if(!silent) message(fxNa,"Note : Unable to find file '",fileName,"' in path '",path,"' but found without specified path !")
      } else chFi <- FALSE                      # if path+fileName not found, check without path
  } }
  if(!chFi) stop(" File ",fileName," was NOT found ",if(length(path) >0) paste(" in path ",path)," !") # nolint: line_length_linter.
  chFi <- file.info(file.path(path[1], fileName))$size > 1
  if(!chFi) stop(" File ",fileName," was found BUT TOO SMALL (size ",file.info(file.path(path[1], fileName))$size," bytes) !")
  if(!grepl("\\.txt$|\\.txt\\.gz$", fileName)) message(fxNa,"Trouble ahead, expecting tabulated text file (the file'",fileName,"' might not be right format) !!")
  if(debug) message(fxNa,"rpd0a ..")

  ## note : reading sample-setup from 'suplAnnotFile' at this place won't allow comparing if number of  samples/columns corresponds to data; do after reading main data
  if(debug) message(fxNa,"rpd0 .. Ready to read", if(length(path) >0) c(" from path ",path[1])," the file  ",fileName[1])


  ## read (main) file
  ## future: look for fast reading of files
  tmp <- try(utils::read.delim(file.path(paFi), stringsAsFactors=FALSE), silent=TRUE)

  if(length(tmp) <1 || inherits(tmp, "try-error") || length(dim(tmp)) <2) {
    if(inherits(tmp, "try-error")) warning("Unable to read input file ('",paFi,"')!  (check if rights to read)") else {
      if(!silent) message(fxNa,"Content of  file '",paFi,"' seeps empty or non-conform !  Returning NULL; check if this is really a ProteomeDiscoverer-file") }
    NULL
  } else {
    if(debug) { message(fxNa,"rpd1 ... dims of initial data : ", nrow(tmp)," li and ",ncol(tmp)," col ")
      rpd1 <- list(tmp=tmp,paFi=paFi,annotCol=annotCol,fileName=fileName) }

    ## locate & extract annotation
    if(length(annotCol) <1) annotCol <- c("Protein.ID","Description","Gene","Contaminant","Sum.PEP.Score","Coverage....","X..Peptides","X..PSMs","X..Unique.Peptides", "X..AAs","MW..kDa.")
    ## option for future: also extract column "MarkedAs"
    PSMCol <- "^Number.of.PSMs.by.Search.Engine"             # pattern searching tag for PSM-data
    PepCol <- "^Number.of.Peptides.by.Search.Engine"         # pattern searching tag for Number of peptides
    ## future option : lateron rename columns called as "Description" to annotCol[2]
    ## below use explicit colnames "Accession","Description", rename if tolower() fits

    .chColNa <- function(x, mat, altern=NULL, renameTo=NULL, silent=FALSE, fxNa=NULL){
      ## check in 'matr' for column-name 'x', if required rename best hit (if no direct hit look using grep, then grep wo case); return corrected mat
      chX <- x %in% colnames(mat)
      if(all(chX)) {
        if(is.character(renameTo) && length(renameTo) ==1) colnames(mat)[match(x, colnames(mat))] <- renameTo   # juste simple rename # nolint # nolint: line_length_linter.
      } else {     # try to localize column to use
        chX <- grep(x, colnames(mat))
        if(length(chX) >0) {
          if(is.character(renameTo) && length(renameTo) ==1) colnames(mat)[chX[1]] <- renameTo else x
          if(!silent & length(chX) >1) message(fxNa,"Found multiple columns containing '",x,"' : ",wrMisc::pasteC(colnames(mat)[chX], quoteC="'"),", using 1st")
        } else {
          if(length(altern==1)) chX <- grep(altern, colnames(mat))
          if(length(chX) > 0) {
            if(is.character(renameTo) && length(renameTo) ==1) colnames(mat)[chX[1]] <- renameTo else x
            if(!silent && length(chX) >1) message(fxNa,"Found multiple columns containing '",x,"' : ",wrMisc::pasteC(colnames(mat)[chX], quoteC="'"),", using 1st")
          } else {
            chX <- grep(tolower(x), tolower(colnames(mat)))
            if(length(chX) > 0) {
              if(is.character(renameTo) && length(renameTo) ==1) colnames(mat)[chX[1]] <- renameTo else x
              if(!silent && length(chX) >1) message(fxNa,"Found multiple columns containing '",tolower(x),"' : ",wrMisc::pasteC(colnames(mat)[chX], quoteC="'"),", using 1st")
            } else stop("Could NOT find column '",x,"' !!\n  (available columns ",wrMisc::pasteC(colnames(mat), quoteC="'"),")") } }
      }
    mat }

    ## check for essential colnames !
    tmp <- .chColNa(annotCol[1], tmp, altern="Accession", rename="Accession", silent=silent, fxNa=fxNa)
    tmp <- .chColNa(annotCol[2], tmp, rename="Description", silent=silent, fxNa=fxNa)
    annotCol[1:2] <- c("Accession","Description")                     # update (just in case..)
    if(is.character(annotCol)) annotColNo <- match(annotCol, colnames(tmp))
    if(length(contamCol) >0) {
      contamCol <- if(is.character(contamCol)) which(colnames(tmp)==contamCol[1]) else as.integer(contamCol[1])
      if(length(contamCol) >0) contamFilter <- TRUE
      annotColNo <- union(annotColNo, contamCol)
    }
    if(debug) { message(fxNa,"rpd2 .. annotColNo : ", wrMisc::pasteC(annotColNo),"     contamCol : ",wrMisc::pasteC(contamCol)," ")
      rpd2 <- list(tmp=tmp,annotCol=annotCol,PSMCol=PSMCol,PepCol=PepCol,fileName=fileName)}

    ## check for R-friendly export
    Rfriendly <- FALSE
    specRepl <- cbind(ini=c("Coverage...."), new=c("Coverage.in.Percent"))
    annotCol2 <- unique(c(sub("X\\.\\.","Number.of.",annotCol), apply(specRepl, 1, function(x) sub(x[1], x[2], annotCol)) ))
    annotColN2 <- match(annotCol2, colnames(tmp))
    if(sum(!is.na(annotColN2)) > sum(!is.na(annotColNo))) { Rfriendly <- TRUE
      annotColNo <- annotColN2
      annotCol <- annotCol2
      if(!silent) message(fxNa,"Setting 'annotCol' to export of 'R-friendly' colnames")}
    if(all(is.na(annotColNo))) stop("Problem with 'annotCol' : Could NOT find any annotation-column")           # nolint: line_length_linter.
    if(any(is.na(annotColNo), na.rm=TRUE)) { if(!identical(annotCol, annotCol2)) message(fxNa,"Can't find column(s) ",wrMisc::pasteC(annotCol[is.na(annotColNo)],quote="'"))
      annotCol <- annotCol[!is.na(annotColNo)] }
    annot <- as.matrix(tmp[,wrMisc::naOmit(annotColNo)])
    if(debug) { message(fxNa,"rpd3 .. Rfriendly: ",Rfriendly,"   ncol annot ",ncol(annot)," cols; colnames : ",wrMisc::pasteC(colnames(annot))," ")
      rpd3 <- list(tmp=tmp,annotCol=annotCol,PSMCol=PSMCol,PepCol=PepCol,fileName=fileName)}

    ## clean 'Description' entries: remove tailing punctuation or open brackets (ie not closed) at end of (truncated) fasta header
    if(cleanDescription) {
      if(debug) { message(fxNa,"rpd3a") }
      annot[,"Description"] <- sub("[[:punct:]]+$","", sub("\\ +$", "", annot[,"Description"]))   # tailing ';' and/or tailing space
      annot[,"Description"] <- sub(" \\([[:alpha:]]*$", "", annot[,"Description"])                # tailing (ie truncated) open '(xxx'
    }

    tmp <- .chColNa("Accession", tmp, silent=silent, fxNa=fxNa)
    annot <- cbind(Accession=annot[,"Accession"], EntryName=NA, GeneName=NA, Species=NA, Contam=NA, SpecType=NA, annot[,-1])   # may be better to name column 'species'

    if(debug) {message(fxNa,"rpd4 .. dim annot: ", nrow(annot)," li and  ",ncol(annot)," cols; colnames : ",wrMisc::pasteC(colnames(annot))," "); rpd4 <- list(annot=annot,tmp=tmp,specPref=specPref,annotCol=annotCol,Rfriendly=Rfriendly,contamCol=contamCol,PSMCol=PSMCol,PepCol=PepCol)}

    if("Contaminant" %in% colnames(annot)) annot[,"Contam"] <- toupper(gsub(" ","",annot[,colnames(tmp)[contamCol]]))

    ## try extract GeneNames from 'Description'
    chPrNa <- is.na(annot[,"GeneName"])
    if(all(chPrNa)) { grLi <- grep("\\ GN=[[:upper:]]{2,}[[:digit:]]", annot[which(chPrNa),"Description"])
      if(length(grLi) >0) { zz <- sub("[[:print:]]+\\ GN=", "", annot[which(chPrNa)[grLi],"Description"])    # remove surplus to left
        annot[which(chPrNa)[grLi],"GeneName"] <- sub("\\ [[:print:]]+","",zz)                                                                    # remove surplus to right
      } }
    if(debug) { message(fxNa,"rpd6 ..  "); rpd6 <- list(annot=annot,specPref=specPref)}
    ## try extract species from 'Description'
    DescrIni <- annot[,"Description"]
    chSpe <- grep("OS=[[:upper:]][[:lower:]]+\\ [[:lower:]]+", DescrIni)
    if(length(chSpe) >0) {         # term OS= exists,
      annot[chSpe,"Description"] <- sub("OS=[[:upper:]][[:lower:]]+\\ [[:lower:]][[:print:]]+", "", DescrIni[chSpe])     # everything left of OS=
      annot[chSpe,"Species"] <- sub("\\ {0,1}[[:upper:]]{2}=[[:print:]]+", "", substr(DescrIni[chSpe], nchar(annot[chSpe,"Description"]) +4, nchar(DescrIni[chSpe])) ) # all right of OS= until next tag
      ## remmove ' (strain ..) ' specification
      annot[chSpe,"Species"] <- sub("\\ \\(strain\\ [[:print:]]+\\)\\ {0,1}$","", annot[chSpe,"Species"])
      annot[chSpe,"Description"] <- sub("\\ $","",annot[chSpe,"Description"])   # remove taining space
    }
    if(debug) {message(fxNa,"rpd6b ..  "); rpdb6 <- list(annot=annot,specPref=specPref)}

    ## separate multi-species (create columns 'Accession','GeneName','Species','SpecType')
    if(!silent) { chSp <- is.na(annot[,"Species"])
      if(any(chSp, na.rm=TRUE) && !all(chSp)) message(fxNa,"Note: ",sum(chSp)," (out of ",nrow(tmp),") lines with unrecognized species")
      if(!all(chSp)) { tab <- table(annot[,"Species"])
        tab <- rbind(names(tab),": ",tab," ;  ")
        if(!silent) message(fxNa,"Count by 'specPref' : ",apply(tab,2,paste)) }}             # all lines assigned

    if(length(specPref) >0) {
      annot <- .extrSpecPref(specPref, annot, useColumn=c("Species","EntryName","GeneName","Accession","Majority.protein.IDs","Fasta.headers"), suplInp=tmp, silent=silent, debug=debug, callFrom=fxNa) }


    if(debug) {message(fxNa,"rpd7 ..  "); rpd7 <- list(annot=annot,specPref=specPref,chSp=chSp,tmp=tmp,quantCol=quantCol )}

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
    if(length(excluCol)==1) {
      excCo <- grep(excluCol, colnames(tmp))
      if(any(duplicated(excCo, quantCol), na.rm=TRUE)) {
        quantCol <- quantCol[-match(excCo, quantCol)]
        if(length(quantCol) <1) stop(msg," (all match to 'excluCol')") else {
          if(!silent) message(fxNa,"Removed ",length(quantColIni) -length(quantCol)," columns")}
      }
    }
    abund <- as.matrix(tmp[,quantCol])                      # abundance val
    if(debug) {message(fxNa,"rpd8 ..  "); rpd8 <- list(annot=annot,specPref=specPref,abund=abund,quantCol=quantCol)}

    ## check & clean abudances
    chNorm <- grep("\\.Normalized\\.", colnames(abund))
    if(length(chNorm)*2 == ncol(abund)) {              # in case Normalized makes 1/2 of columns use non-normalized
      abund <- abund[,-chNorm]
    }
    colnames(abund) <- sub("^Abundances\\.Normalized\\._{0,1}|^abundances\\.Normalized\\._{0,1}|^Abundances{0,1}_{0,1}|^abundances{0,1}_{0,1}","",colnames(abund))
    chNum <- is.numeric(abund)
    if(!chNum) {abund <- apply(tmp[,quantCol], 2, wrMisc::convToNum, convert="allChar", silent=silent, callFrom=fxNa)}

    ## remove heading 'X..' from headers (only if header won't get duplicated
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
    if(debug) { message(fxNa,"rpd9 .. dim annot ",nrow(annot)," and ",ncol(annot)); rpd9 <- list(annot=annot,tmp=tmp,abund=abund,sampleNames=sampleNames,specPref=specPref,annotCol=annotCol,Rfriendly=Rfriendly,contamCol=contamCol,PSMCol=PSMCol,PepCol=PepCol,infoDat=infoDat) }

    ## add custom sample names (if provided)
    if(length(sampleNames) ==ncol(abund) && ncol(abund) >0) {
      if(debug) { message(fxNa,"rpd9b") }
      if(length(unique(sampleNames)) < length(sampleNames)) {
        if(!silent) message(fxNa,"Custom sample names not unique, correcting to unique")
        sampleNames <- wrMisc::correctToUnique(sampleNames, callFrom=fxNa) }
      colnames(abund) <- sampleNames
      if(debug) { message(fxNa,"rpd9c") }
    }

    ## (optional) filter by FDR  (so far use 1st of list where matches are found from argument FDRCol)
    if(length(FDRCol) >0) {
      ## stand : "Found.in.Sample...S32..F32..Sample"   Rfriendly : "Found.in.Sample.in.S33.F33.Sample"
      ## stand : "X..Protein.Groups"                    Rfriendly : "Number.of.Protein.Groups"
      chFDR <- lapply(FDRCol, function(x) {z <- grep(x[1], colnames(tmp)); if(length(z) ==ncol(abund)) z else NULL})
      names(chFDR) <- sapply(FDRCol, function(x) x[1])
      chFDR <- chFDR[which(sapply(chFDR, length) >0)]
      if(length(chFDR) >0) {
        i <- 1    # so far just use 1st instance matching
        searchFor <- FDRCol[[which(sapply(FDRCol, function(x) x[1]) %in% names(chFDR)[i])]]
        filtFdrHi <- tmp[,chFDR[[i]]] == searchFor[2]  # find occurances of best tag 'High'
        roSu <- rowSums(filtFdrHi) <1
        if(all(roSu) && !silent) message(fxNa,"NONE of the lines/proteins had any '",searchFor[1],"' in column(s) '",searchFor[2],"' !!  This is probably not a good filtering-parameter, ignoring")
        if(any(roSu, na.rm=TRUE) && !all(roSu)) { if(!silent) message(fxNa,"Removing ",sum(roSu)," lines/proteins without ANY '",searchFor[2],"' in columns '",searchFor[1],"'")
          rmLi <- -1*which(roSu)
          annot <- annot[rmLi,]
          abund <- abund[rmLi,]
          filtFdrHi <- filtFdrHi[rmLi,]   # useful lateron ?
          tmp <- tmp[rmLi,] }
      }
    }
    if(debug) { message(fxNa,"rpd11 .. length(FDRCol) ",length(FDRCol),"   dim annot ",nrow(annot)," and ",ncol(annot))}

    ## rownames : check if Accession is unique
    chAc <- duplicated(annot[,"Accession"], fromLast=FALSE)
    if(any(chAc, na.rm=TRUE)) {
       getLiToRemove <- function(x,useCol=c("rowNo","Contaminant","SpecType")) {  # return index for all lines to remove from matrix ...
          if(is.data.frame(x)) x <- as.matrix(x)
          spe <- grep("^species", x[,useCol[3]])
          if(length(spe) >0) {
            rmLi <- x[which(1:nrow(x) != spe[1]), useCol[1]]
          } else {                        ## look for any lines marked as Contaminant="true", then mark other(s) for remove
            rmLi <- if(any(tolower(x[,useCol[2]])=="true", na.rm=TRUE)) x[which(tolower(x[,useCol[2]]) !="true") ,useCol[1]]  }
          as.integer(rmLi) }

      ## check if one of duplicated lines is marked as Contaminant -> remove non-contaminant, BUT NOT 'speciesX' ?
      if(contamFilter) {                # ready to correct (if possible) duplicated 'Accession' entries
        ## elaborate procedure for removing duplicate Accession lines : 'fuse' annot where no NA & use quantification-line with fewest NAs
        ## need to separate all groups of repeated IDs & treat separately
        annot <- cbind(annot, rowNo=1:nrow(tmp))
        duplAc <- unique(annot[which(chAc), "Accession"])
        ## need to remove duplicated lines which are not marked as Contaminant="True"
        chAc2 <- duplicated(annot[,"Accession"], fromLast=TRUE)
        rmLi <- chAc | chAc2
        ## find lines where is not Contaminant="True" (and keep contaminant)
        annot <- cbind(annot, iniIndex=1:nrow(annot), nNA=rowSums(is.na(abund)))
        useCol2 <- c("Accession","GeneName","Species","Contam","SpecType","Description","Contaminant", "iniIndex","nNA")  # the last 2 are added within function
        useCol2 <- wrMisc::naOmit(match(useCol2,colnames(annot)))
        abund <- cbind(abund,iniIndex=1:nrow(abund))
        rmAbund <- as.integer(unlist(by(abund[which(rmLi),], annot[which(rmLi),"Accession"], function(x) x[(1:nrow(x))[-which.min(rowSums(is.na(x)))],ncol(x)])))
        rmAnnot2 <- as.integer(unlist(by(annot[which(rmLi),], annot[which(rmLi),"Accession"], function(x) x[2:nrow(x),ncol(x) -1])))
        rmAnnot <- which(chAc)
        for(j in unique(annot[which(rmLi),useCol2][1])) {    # need loop for 'fusing' columns with fewes NAs and recording which lines should be removed
          x <- annot[which(annot[,"Accession"] %in% j),]
          useLi <- apply(0+is.na(x),2,which.min)
          if(any(useLi >1, na.rm=TRUE)) for(i in 2:max(useLi)) annot[as.integer(x[1,"iniIndex"]),which(useLi==i)] <- annot[as.integer(x[i,"iniIndex"]),which(useLi==i)]
        }
        if(length(rmAnnot) >0) {annot <- annot[-rmAnnot,]; tmp <- tmp[-rmAnnot,]
          abund <- abund[-rmAbund,]
          if(!silent) message(fxNa,"Removing ",length(rmAnnot)," lines due to duplicated Accessions (typically due to contaminants)")
        }
        annot <- annot[,-ncol(annot) +(1:0)]                  # remove extra columns (ie "iniIndex","nNA")
        abund <- abund[,-ncol(abund)]                         # remove extra column (ie "iniIndex")
        chAc <- duplicated(annot[,"Accession"], fromLast=FALSE)
        if(debug) { message(fxNa,"rpd11b .. dim abund ",nrow(abund)," and ",ncol(abund)); rpd9 <- list(annot=annot,tmp=tmp,abund=abund,sampleNames=sampleNames,specPref=specPref,annotCol=annotCol,Rfriendly=Rfriendly,contamCol=contamCol,PSMCol=PSMCol,PepCol=PepCol,infoDat=infoDat)}
        }}
    ## Now we are ready to add unique rownames
    if(any(chAc, na.rm=TRUE)) {
      if(!silent) message(fxNa,sum(chAc)," (out of ",length(chAc),") cases of duplicated 'Accession' exist, adding extensions for use as rownames")
      rownames(tmp) <- rownames(annot) <- wrMisc::correctToUnique(annot[,"Accession"], sep="_", atEnd=TRUE, callFrom=fxNa)
    } else rownames(abund) <- rownames(annot) <- annot[,"Accession"]

    ## optional/additional counting results (PSM, no of peptides)
    PSMCol <- if(length(PSMCol) ==1) grep(PSMCol,colnames(tmp)) else NULL
    PepCol <- if(length(PepCol) ==1) grep(PepCol,colnames(tmp)) else NULL
    usTy <- c("PSM","NoOfPeptides")[which(c(length(PSMCol),length(PepCol)) ==ncol(abund))]
    if(length(usTy) >0) {
      counts <- array(NA,dim=c(nrow(abund),ncol(abund),length(usTy)), dimnames=list(rownames(abund),colnames(abund),usTy))
      if("PSM" %in% usTy) counts[,,"PSM"] <- as.matrix(tmp[,PSMCol])
      if("NoOfPeptides" %in% usTy) counts[,,"NoOfPeptides"] <- as.matrix(tmp[,PepCol])
    } else counts <- NULL
    if(debug) {message(fxNa,"rpd12 .. "); rpd12 <- list(annot=annot,tmp=tmp,abund=abund,sampleNames=sampleNames,specPref=specPref,annotCol=annotCol,refLi=refLi,Rfriendly=Rfriendly,contamCol=contamCol,PSMCol=PSMCol,PepCol=PepCol,infoDat=infoDat)}

    ## check for reference for normalization
    refLiIni <- refLi
    if(is.character(refLi) && length(refLi)==1) {
      refLi <- which(annot[,"SpecType"]==refLi)
      if(length(refLi) <1 && identical(refLiIni, "mainSpe")) refLi <- which(annot[,"SpecType"] =="mainSpecies")              # fix compatibility problem  'mainSpe' to 'mainSpecies'
      if(length(refLi) <1 ) { refLi <- 1:nrow(abund)
        if(!silent) message(fxNa,"Could not find any proteins matching argument 'refLi=",refLiIni,"', ignoring ...")
      } else {
        if(!silent) message(fxNa,"Normalize using (custom) subset of ",length(refLi)," lines specified as '",refLiIni,"'")}}    # may be "mainSpe"

    ## take log2 & normalize
    quant <- try(wrMisc::normalizeThis(log2(abund), method=normalizeMeth, mode="additive", refLines=refLi, silent=silent, debug=debug, callFrom=fxNa), silent=TRUE)
    if(debug) { message(fxNa,"rpd13 .. dim quant: ", nrow(quant)," li and  ",ncol(quant)," cols; colnames : ",wrMisc::pasteC(colnames(quant))," "); rpd13 <- list(annot=annot,tmp=tmp,abund=abund,quant=quant,sampleNames=sampleNames,specPref=specPref,annotCol=annotCol,Rfriendly=Rfriendly,contamCol=contamCol,PSMCol=PSMCol,PepCol=PepCol,infoDat=infoDat)}

    ### GROUPING OF REPLICATES AND SAMPLE META-DATA
    if(length(suplAnnotFile) >0 || length(sdrf) >0) {
      setupSd <- readSampleMetaData(sdrf=sdrf, suplAnnotFile=suplAnnotFile, quantMeth="PD", path=path, abund=utils::head(quant), groupPref=groupPref, silent=silent, debug=debug, callFrom=fxNa)
    }
    if(debug) {message(fxNa,"rpd13b .."); rpd13b <- list()}

    ## finish groups of replicates & annotation setupSd
    setupSd <- .checkSetupGroups(abund=abund, setupSd=setupSd, gr=gr, sampleNames=sampleNames, quantMeth="PD", silent=silent, debug=debug, callFrom=fxNa)
    colnames(quant) <- colnames(abund) <- if(length(setupSd$sampleNames)==ncol(abund)) setupSd$sampleNames else setupSd$groups
    if(length(dim(counts)) >1 && length(counts) >0) colnames(counts) <- setupSd$sampleNames

    if(debug) {message(fxNa,"Read sample-meta data, rpd14"); rpd14 <- list(sdrf=sdrf,suplAnnotFile=suplAnnotFile,abund=abund, quant=quant,refLi=refLi,annot=annot,setupSd=setupSd,sampleNames=sampleNames)}

    ## main plotting of distribution of intensities
    custLay <- NULL
    if(is.numeric(plotGraph) && length(plotGraph) >0) {custLay <- as.integer(plotGraph); plotGraph <- TRUE} else {
        if(!isTRUE(plotGraph)) plotGraph <- FALSE}
    if(plotGraph) .plotQuantDistr(abund=abund, quant=quant, custLay=custLay, normalizeMeth=normalizeMeth, softNa="Proteome Discoverer",
      refLi=refLi, refLiIni=refLiIni, tit=titGraph, silent=silent, callFrom=fxNa, debug=debug)

    ## meta-data
    notes <- c(inpFile=paFi, qmethod="ProteomeDiscoverer", qMethVersion=if(length(infoDat) >0) unique(infoDat$Software.Revision) else NA,
    	rawFilePath= if(length(infoDat) >0) infoDat$File.Name[1] else NA, normalizeMeth=normalizeMeth, call=match.call(),
      created=as.character(Sys.time()), wrProteo.version=utils::packageVersion("wrProteo"), machine=Sys.info()["nodename"])
    ## final output
    if(isTRUE(separateAnnot)) list(raw=abund, quant=quant, annot=annot, counts=counts, sampleSetup=setupSd, quantNotes=parametersD, notes=notes) else data.frame(quant,annot)
  }
}
  
