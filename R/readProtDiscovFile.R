#' Read tabulated files exported by ProteomeDiscoverer
#'
#' Protein quantification results form \href{https://www.thermofisher.com/order/catalog/product/OPTON-30812}{Thermo ProteomeDiscoverer} 
#' which were exported as tabulated text can be imported and relevant information extracted. 
#' The final output is a list containing 3 elements: \code{$annot}, \code{$raw} and optional \code{$quant}, or returns data.frame with entire content of file if \code{separateAnnot=FALSE}.
#' @details
#' This function has been developed using Thermo ProteomeDiscoverer versions 2.2 to 2.5.
#' The format of resulting files at export also depends which columns are chosen as visible inside ProteomeDiscoverer and subsequently get chosen for export.
#' Using the argument \code{infoFile} it is possible to specify a specific file (or search for default file) to read for extracting file-names as sample-names and other experiment realted information.  
#' If a column named \code{contamCol} is found, the data will be lateron filtered to remove all contaminants, set to \code{NULL} for keeping all contaminants 
#' This function replaces the depreciated function \code{readPDExport}.
#' 
#' @param fileName (character) name of file to be read  
#' @param path (character) path of file to be read
#' @param normalizeMeth (character) normalization method (will be sent to  \code{\link[wrMisc]{normalizeThis}}) 
#' @param sampleNames (character) new column-names for quantification data (ProteomeDiscoverer does not automatically use file-names from spectra); this argument has priority over \code{infoFile}
#' @param infoFile (character or logical) filename containing additional information about MS-samples (produced by ProteomeDiscoverer default '*.InputFiles.txt'), 
#'  if \code{TRUE} the first file in \code{path} containing the default name will be used. If no specific \code{sampleNames} given, the filenames will will be trimmed to remove redundant text and used as \code{sampleNames}. 
#'  Besides, ProteomeDiscoverer version number and full raw-file path will be extracted for $notes in fial output.
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
#' fiNa <- "tinyPD_allProteins.txt.gz"
#' dataPD <- readProtDiscovFile(file=fiNa, path=path1)
#' summary(dataPD$quant)
#' 
#' @export
readProtDiscovFile <- function(fileName, path=NULL, normalizeMeth="median", sampleNames=NULL, infoFile=TRUE, read0asNA=TRUE, quantCol="^Abundances*", 
  annotCol=NULL, contamCol="Contaminant", refLi=NULL, separateAnnot=TRUE, FDRCol=list(c("^Protein.FDR.Confidence","High"), c("^Found.in.Sample.","High")), plotGraph=TRUE, tit="Proteome Discoverer", graphTit=NULL, wex=1.6,
  specPref=c(conta="CON_|LYSC_CHICK", mainSpecies="OS=Homo sapiens"), silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## read ProteomeDiscoverer exported txt
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readProtDiscovFile")
  oparMar <- graphics::par("mar")
  reqPa <- c("utils","wrMisc")
  chPa <- sapply(reqPa, requireNamespace, quietly=TRUE)
  if(any(!chPa)) stop("package(s) '",paste(reqPa[which(!chPa)], collapse="','"),"' not found ! Please install first from CRAN")  
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  excluCol <- "^Abundances.Count"   # exclude this from quantifications columns
  cleanDescription <- TRUE          # clean 'Description' for artifacts of truncated text (tailing ';' etc)         
  
  ## check if path & file exist
  msg <- "invalid entry for 'fileName'"
  if(length(fileName) >1) { fileName <- fileName[1]
    if(!silent) message(fxNa," 'fileName' shoud be of length=1, using 1st value")
  } else { if(length(fileName) <1) stop(msg) else if(nchar(fileName) <0) stop(msg)}
  paFi <- fileName                      # presume (& correct if path is given)
  chFi <- file.exists(fileName)         # presume (& correct otherwise)
  if(length(path) >0) if(!file.exists(path[1])) { path <- NULL
    if(!silent) message(fxNa,"Invalid path '",path[1],"'  (not existing), ignoring...") }
  if(length(path) >0) { chFi <- file.exists(file.path(path[1], fileName))
    if(chFi) paFi <- file.path(path[1], fileName) else {
      if(grepl(paste0("^",path[1]), fileName)) {chFi <- file.exists(fileName); path <- NULL} else chFi <- FALSE       # if path+fileName not found, check if 'path' should be omitted if already contained in fileName
  } } 
  if(!chFi) stop(" File ",fileName," was NOT found ",if(length(path) >0) paste(" in path ",path)," !")
  if(!grepl("\\.txt$|\\.txt\\.gz$", fileName)) message(fxNa,"Trouble ahead, expecting tabulated text file (this file might not be right format) !!")

  ## prepare for reading files
  ## check for separate file 'infoFile' and store as infoFile
  infoDat <- infoFi <- NULL                     # initialize
  if(length(infoFile)==1) {
    infoFile <- if(isTRUE(infoFile)) "InputFiles.txt$" else NULL
    infoFi <- if(length(path) >0) file.path(path[1],infoFile) else infoFile    # with path (if applicable)
    if(!file.exists(infoFi)) { 
      infoFi <- if(length(path) >0) dir(path, pattern=infoFile) else dir(pattern=infoFile) 
      if(length(infoFi) >0) message(fxNa,"Setting 'infoFile' to '",infoFi[1],"'")}
    if(debug) message(fxNa,"rpd0a ..")
    if(length(infoFi) >0) { if(file.exists(infoFi)) {   ## try reading infoFile 
      infoDat <- try(utils::read.delim(file=if(length(path) >0)file.path(path[1],infoFi[1]) else infoFi[1], stringsAsFactors=FALSE), silent=TRUE)
      if(inherits(infoDat, "try-error")) {if(!silent) message(fxNa,"Failed to open file '",infoFi[1],"' for getting additional information about experiment")
      } else {
        useColInf <- c("File.ID","File.Name","Creation.Date","RT.Range.in.min","RT.Range..min.","Instrument.Name","Software.Revision","Ref.File.ID","Ref..File.ID")
        chCol <- match(useColInf, colnames(infoDat))
        if(all(is.na(chCol))) { 
          if(!silent & is.character(infoFile)) message(fxNa," Can't find any of columns expected in file '",infoFi[1],"' for additional information about experiment; ignoring 'infoFile'")
          infoDat <- NULL
        } else infoDat <- infoDat[,wrMisc::naOmit(chCol)]  
      }
    } else if(!silent) message(fxNa,"Can't find file '",file.path(path,infoFi[1]),"' for additional information about experiment; ignoring 'infoFile'")}
  }
  if(debug) message(fxNa,"rpd1 .. Ready to read", if(length(path) >0) c(" from path ",path[1])," the file  ",fileName[1])
  
  ## read (main) file
  ## future: look for fast reading of files
  tmp <- try(utils::read.delim(file.path(paFi), stringsAsFactors=FALSE), silent=TRUE)
  if(inherits(tmp, "try-error")) stop("Unable to read input file ('",paFi,"')!") 
  if(debug) { message(fxNa,"rpd2 .. dims of initial data : ", nrow(tmp)," li and ",ncol(tmp)," col ")}
  
  ## locate & extract annotation 
  if(length(annotCol) <1) annotCol <- c("Accession","Description","Gene","Contaminant","Sum.PEP.Score","Coverage....","X..Peptides","X..PSMs","X..Unique.Peptides", "X..AAs","MW..kDa.")
  ## option for future: also extract column "MarkedAs"
  PSMCol <- "^Number.of.PSMs.by.Search.Engine"             # pattern searching tag for PSM-data
  PepCol <- "^Number.of.Peptides.by.Search.Engine"         # pattern searching tag for Number of peptides  
  ## future option : lateron rename columns called as "Description" to annotCol[2]    
  ## below use explicit colnames "Accession","Description", rename if tolower() fits

  .chColNa <- function(x, mat, renameTo=NULL, silent=FALSE, fxNa=NULL){
    ## check for column-name 'x' in 'matr', if required rename best hit (if no direct hit look using grep, then grep wo case); return corrected mat
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
    
  ## check for essential colnames !
  tmp <- .chColNa(annotCol[1], tmp, rename="Accession", silent=silent, fxNa=fxNa) 
  tmp <- .chColNa(annotCol[2], tmp, rename="Description", silent=silent, fxNa=fxNa) 
  if(is.character(annotCol)) annotColNo <- match(annotCol, colnames(tmp))
  if(length(contamCol) >0) {
    contamCol <- if(is.character(contamCol)) which(colnames(tmp)==contamCol[1]) else as.integer(contamCol[1])
    if(length(contamCol) >0) contamFilter <- TRUE
    annotColNo <- union(annotColNo, contamCol)
  }
  if(debug) { message(fxNa,"rpd3 .. annotColNo : ", wrMisc::pasteC(annotColNo),"     contamCol : ",wrMisc::pasteC(contamCol)," ")}
  
  ## check for R-friendly export
  Rfriendly <- FALSE
  specRepl <- cbind(ini=c("Coverage...."), new=c("Coverage.in.Percent"))
  annotCol2 <- unique(c(sub("X\\.\\.","Number.of.",annotCol), apply(specRepl, 1, function(x) sub(x[1], x[2], annotCol)) ))
  annotColN2 <- match(annotCol2, colnames(tmp))
  if(sum(!is.na(annotColN2)) > sum(!is.na(annotColNo))) { Rfriendly <- TRUE
    annotColNo <- annotColN2
    annotCol <- annotCol2
    if(!silent) message(fxNa," setting 'annotCol' to export of 'R-friendly' colnames")}  
  if(all(is.na(annotColNo))) stop(" Problem with 'annotCol' : Could NOT find any annotation-column")
  if(any(is.na(annotColNo))) { if(!identical(annotCol, annotCol2)) message(fxNa,"Can't find column(s) ",wrMisc::pasteC(annotCol[is.na(annotColNo)],quote="'"))
    annotCol <- annotCol[!is.na(annotColNo)] }
  annot <- as.matrix(tmp[,wrMisc::naOmit(annotColNo)])
  if(debug) { message(fxNa,"rpd4 .. Rfriendly: ",Rfriendly,"   ncol annot ",ncol(annot)," cols; colnames : ",wrMisc::pasteC(colnames(annot))," ") }

  ## clean 'Description' entries: remove tailing punctuation or open backets (ie not closed) at end of fasta header
  if(cleanDescription) {
    if(debug) { message(fxNa,"rpd4b") }   
    annot[,"Description"] <- sub("[[:punct:]]+$","", sub("\\ +$", "", annot[,"Description"]))                 # tailing ';' of tailing space 
    annot[,"Description"] <- sub(" \\([[:alpha:]]*$", "", annot[,"Description"])           # tailing (ie truncated) open '(xxx'
  }

if(debug) { message(fxNa,"rpd4c") }
  tmp <- .chColNa("Accession",tmp,silent=silent,fxNa=fxNa) 
  annot <- cbind(Accession=annot[,"Accession"], EntryName=NA, GeneName=NA, Species=NA, Contam=NA, SpecType=NA, annot[,-1])   # may be better to name column 'species'
  if(length(specPref) >0) for(i in 1:length(specPref)) {         # locate specPref
    naSp <- if(i==1) "conta" else {if(i==2) "mainSpe" else paste0("species",i-1)}
    chSp <- unlist(specPref[i]) 
    chSp <- if(length(chSp) >1) match(chSp, annot[,annotCol[1]]) else grep(unlist(specPref[i]), annot[,annotCol[2]])    # line-no to set tag a 'SpecType'
    if(length(chSp) >0) {
      chNa <- is.na(annot[chSp,"SpecType"])
      if(any(!chNa) & !silent) message(fxNa,"Beware, ",sum(!chNa)," 'SpecType' will be overwritten by ",naSp)
      annot[chSp,"SpecType"] <- naSp }
  }
  if("Contaminant" %in% colnames(annot)) annot[,"Contam"] <- toupper(gsub(" ","",annot[,colnames(tmp)[contamCol]]))
  if(debug) {message(fxNa,"rpd5 .. dim annot: ", nrow(annot)," li and  ",ncol(annot)," cols; colnames : ",wrMisc::pasteC(colnames(annot))," ")}

  ## separate multi-species (create columns 'Accession','GeneName','Species','SpecType')
  #if(is.character(contamCol)) contamCol <- which(colnames(annot)==contamCol)
  ## had already message if column for contaminations not soecified/found
  if(length(specPref) <1) message(fxNa,"Note: argument 'specPref' not specifed (empty)")

  if(debug) {message(fxNa,"rpd6 ..  ")}
  
  ## try extract GeneNames from 'Descripton'
  chPrNa <- is.na(annot[,"GeneName"])
  if(all(chPrNa)) { grLi <- grep("\\ GN=[[:upper:]]{2,}[[:digit:]]", annot[which(chPrNa),"Description"])
    if(length(grLi) >0) { zz <- sub("[[:print:]]+\\ GN=", "", annot[which(chPrNa)[grLi],"Description"])    # remove surplus to left
      annot[which(chPrNa)[grLi],"GeneName"] <- sub("\\ [[:print:]]+","",zz)                                                                    # remove surplus to right
    } }
  ## try extract species from 'Descripton'
  DescrIni <- annot[,"Description"]
  chSpe <- grep("OS=[[:upper:]][[:lower:]]+\\ [[:lower:]]+", DescrIni)
  if(length(chSpe) >0) {  # term OS= exists, 
    annot[chSpe,"Description"] <- sub("OS=[[:upper:]][[:lower:]]+\\ [[:lower:]][[:print:]]+", "", DescrIni[chSpe])     # everything left of OS=
    annot[chSpe,"Species"] <- sub("\\ {0,1}[[:upper:]]{2}=[[:print:]]+", "", substr(DescrIni[chSpe], nchar(annot[chSpe,"Description"]) +4, nchar(DescrIni[chSpe])) ) # all right of OS= until next tag
    if(TRUE) annot[chSpe,"Species"] <- sub("\\ \\(strain\\ [[:print:]]+\\)\\ {0,1}$","", annot[chSpe,"Species"])
    annot[chSpe,"Description"] <- sub("\\ $","",annot[chSpe,"Description"])
  }
  
  if(!silent) { chSp <- is.na(annot[,"Species"])
    if(any(chSp) & !all(chSp)) message(fxNa,"Note: ",sum(chSp)," (out of ",nrow(tmp),") lines with unrecognized species")
    if(!all(chSp)) { tab <- table(annot[,"Species"])
      tab <- rbind(names(tab),": ",tab," ;  ")
      if(!silent) message(fxNa,"Count by 'specPref' : ",apply(tab,2,paste)) }}             # all lines assigned   
  if(debug) { message(fxNa,"rpd7 ..  ") }
   
  ## locate & extract abundance/quantitation data
  msg <- " CANNOT find ANY quantification columns" 
  if(length(quantCol) >1) { 
    ## explicit columns (for abundance/quantitation data)
    ## problem : extract '^Abundances*' but NOT 'Abundances.Count.*'
    quantColIni <- quantCol <- grep(quantCol[1],colnames(tmp))
    if(length(quantCol) <1) stop(msg,"  ('",quantCol,"')")
  } else {
    ## pattern search (for abundance/quantitation data)
    if(length(quantCol) <1) { quantCol <- "^Abundance" 
      if(!silent) message(fxNa,"Setting argument 'quantCol' to '^Abundance'")}
    quantCol <- grep(quantCol, colnames(tmp))
    quantColIni <- quantCol
    if(length(quantCol) <1) { quantCol <- grep(tolower(quantColIni), tolower(colnames(tmp)))
      if(!silent) message(fxNa,msg,", trying all as lower caps : now ",length(quantCol), "columns") }    
    if(length(quantCol) <1) stop(msg," specified in argument 'quantCol' !") }
  ## check for columns to exclude (like 'Abundances.Count.')                                       
  if(length(excluCol)==1) {
    excCo <- grep(excluCol, colnames(tmp))
    if(any(duplicated(excCo, quantCol))) {
      quantCol <- quantCol[-match(excCo, quantCol)]
      if(length(quantCol) <1) stop(msg," (all match to 'excluCol')") else {
        if(!silent) message(fxNa,"Removed ",length(quantColIni) -length(quantCol)," columns")}
    }    
  }
  abund <- as.matrix(tmp[,quantCol])                      # abundance val
  if(debug) {message(fxNa,"rpd8 ..  ")}

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
    if(any(chDu)) newNa[which(chDu)] <- colnames(annot)[chXCol][which(chDu)]
    colnames(annot)[chXCol] <- newNa }
  ## remove heading/tailing spaces (first look which columns might be subject to this treatment)
  ch1 <- list(A=grep("^ +",annot[1,]), B=grep("^ +",annot[2,]), C=grep("^ +",annot[floor(mean(nrow(annot))),]), D=grep("^ +",annot[nrow(annot),]) )
  chCo <- unique(unlist(ch1))
  annot[,chCo] <- sub("^ +","",sub(" +$","",annot[,chCo]))   # remove heading/tailing spaces
  if(debug) { message(fxNa,"rpd9 .. dim annot ",nrow(annot)," and ",ncol(annot)) }

  ## add custom sample names
  if(length(sampleNames) ==ncol(abund) & ncol(abund) >0) {
    if(length(unique(sampleNames)) < length(sampleNames)) {
      if(!silent) message(fxNa,"Custom sample names not unique, correcting to unique")
      sampleNames <- wrMisc::correctToUnique(sampleNames, callFrom=fxNa) }
    colnames(abund) <- sampleNames 
  } else if(length(infoDat) >0) if(nrow(infoDat)==ncol(abund) & ncol(abund) >0) {  # otherwise try using original raw-file names form file  *InputFiles.txt
    sampleNames <- wrMisc::correctToUnique(wrMisc::.trimFromEnd(wrMisc::.trimFromStart(basename(infoDat$File.Name))), callFrom=fxNa) 
    if(!silent) message(fxNa,"Using sample names based on .raw filenames from file '",infoFi[1],"'")
    colnames(abund) <- sampleNames }   
  if(isTRUE(read0asNA)) { is0 <- which(abund <= 0)
    if(length(is0) >0) { if(!silent) message(fxNa,"Replacing ",length(is0)," (",round(100*length(is0)/prod(dim(abund)),3),"%) instances of '0' by NA")
      abund[which(is0)] <- NA }  }    
  if(debug) { message(fxNa,"rpd10 .. sampleNames",wrMisc::pasteC(sampleNames,quoteC="'"))}

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
      if(all(roSu) & !silent) message(fxNa,"NONE of the lines/proteins had any '",searchFor[1],"' in column(s) '",searchFor[2],"' !!  This is probably not a good filtering-parameter, ignoring")
      if(any(roSu) & !all(roSu)) { if(!silent) message(fxNa,"Removing ",sum(roSu)," lines/proteins without ANY '",searchFor[2],"' in columns '",searchFor[1],"'")
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
  if(any(chAc)) {
     getLiToRemove <- function(x,useCol=c("rowNo","Contaminant","SpecType")) {  # return index for all lines to remove from matrix ...
        if(is.data.frame(x)) x <- as.matrix(x)
        spe <- grep("^species", x[,useCol[3]])
        if(length(spe) >0) { 
          rmLi <- x[which(1:nrow(x) != spe[1]), useCol[1]]
        } else {                        ## look for any lines marked as Contaminant="true", then mark other(s) for remove
          rmLi <- if(any(tolower(x[,useCol[2]])=="true")) x[which(tolower(x[,useCol[2]]) !="true") ,useCol[1]]  }
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
        if(any(useLi >1)) for(i in 2:max(useLi)) annot[as.integer(x[1,"iniIndex"]),which(useLi==i)] <- annot[as.integer(x[i,"iniIndex"]),which(useLi==i)]
      } 
      if(length(rmAnnot) >0) {annot <- annot[-rmAnnot,]; tmp <- tmp[-rmAnnot,]
        abund <- abund[-rmAbund,]
        if(!silent) message(fxNa,"Removing ",length(rmAnnot)," lines due to duplicated Accessions (typically due to contaminants)")
      }
      annot <- annot[,-ncol(annot) +(1:0)]                  # remove extra columns (ie "iniIndex","nNA")
      abund <- abund[,-ncol(abund)]                         # remove extra column (ie "iniIndex")
      chAc <- duplicated(annot[,"Accession"], fromLast=FALSE)
      if(debug) { message(fxNa,"rpd11b .. dim abund ",nrow(abund)," and ",ncol(abund))}
      }}
  ## Now we are ready to add unique rownames
  if(any(chAc)) {
    if(!silent) message(fxNa,sum(chAc)," (out of ",length(chAc),") cases of duplicated 'Accession' exist, adding extensions for use as rownames")
    rownames(tmp) <- rownames(annot) <- wrMisc::correctToUnique(annot[,"Accession"], sep="_", atEnd=TRUE, callFrom=fxNa)    
  } else rownames(abund) <- rownames(annot) <- annot[,"Accession"]

  ## optional counting results (PSM, no of peptides)
  PSMCol <- if(length(PSMCol) ==1) grep(PSMCol,colnames(tmp)) else NULL
  PepCol <- if(length(PepCol) ==1) grep(PepCol,colnames(tmp)) else NULL
  usTy <- c("PSM","NoOfPeptides")[which(c(length(PSMCol),length(PepCol)) ==ncol(abund))] 
  if(length(usTy) >0) {
    counts <- array(NA,dim=c(nrow(abund),ncol(abund),length(usTy)), dimnames=list(rownames(abund),colnames(abund),usTy))
    if("PSM" %in% usTy) counts[,,"PSM"] <- as.matrix(tmp[,PSMCol])
    if("NoOfPeptides" %in% usTy) counts[,,"NoOfPeptides"] <- as.matrix(tmp[,PepCol])
  } else counts <- NULL 
  if(debug) {message(fxNa,"rpd12 .. ")}
        
  ## check for reference for normalization
  refLiIni <- refLi
  if(is.character(refLi) & length(refLi)==1) { refLi <- which(annot[,"SpecType"]==refLi)
    if(length(refLi) <1) message(fxNa,"Could not find any protein matching argument 'refLi', ignoring ...") else {
      if(!silent) message(fxNa,"Normalize using (custom) subset of ",length(refLi)," lines")}}    # may be "mainSpe"
   #makes no sense#if(length(refLi) <1) refLi <- NULL
  ## take log2 & normalize
  quant <- wrMisc::normalizeThis(log2(abund), method=normalizeMeth, refLines=refLi, silent=silent, callFrom=fxNa) 
  if(debug) { message(fxNa,"rpd13 .. dim quant: ", nrow(quant)," li and  ",ncol(quant)," cols; colnames : ",wrMisc::pasteC(colnames(quant))," ")}

  ## plot distribution of intensities
  custLay <- NULL
  if(is.numeric(plotGraph) & length(plotGraph) >0) {custLay <- as.integer(plotGraph); plotGraph <- TRUE} else {
    if(!isTRUE(plotGraph)) plotGraph <- FALSE}
  if(plotGraph) {
    if(debug) {message(fxNa,"rpd14 .. length custLay ", length(custLay) )}
    if(length(custLay) >0) graphics::layout(custLay) else {if(!identical(normalizeMeth,"none") & length(quant) >0) graphics::layout(1:2)}
    graphics::par(mar=c(3, 3, 3, 1))                          # mar: bot,le,top,ri
    if(length(graphTit) >0) message(fxNa,"Argument 'graphTit' is depreciated, please rather use 'tit'")
    if(is.null(tit) & !is.null(graphTit)) tit <- graphTit     # for derpreciated argument
    if(is.null(tit)) tit <- "ProteomeDiscoverer quantification "
    chGr <- try(find.package("wrGraph"), silent=TRUE)
    chSm <- try(find.package("sm"), silent=TRUE)
    misPa <- c(inherits(chGr, "try-error"), inherits(chSm, "try-error"))
    titSu <- if(length(refLi) >0) paste0(c(" by ",if(length(refLiIni) >1) c(length(refLi)," selected lines") else c("'",refLiIni,"'")),collapse="")  else NULL
    if(debug) { message(fxNa,"rpd15 .. misPa ", wrMisc::pasteC(misPa,quoteC="'") )}
    if(any(misPa)) { 
      if(!silent) message(fxNa,"Missing package ",wrMisc::pasteC(c("wrGraph","sm")[which(misPa)],quoteC="'")," for drawing vioplots")
      ## wrGraph not available : simple boxplot  
      graphics::boxplot(log2(abund), main=paste(tit,"(initial)",sep=" "), las=1, outline=FALSE) 
      graphics::abline(h=round(log2(stats::median(abund,na.rm=TRUE))) +c(-2:2), lty=2, col=grDevices::grey(0.6)) 
      ## plot normalized
      if(identical(normalizeMeth,"none") | length(quant) <0) {
        graphics::boxplot(quant, main=paste(tit," (",normalizeMeth,"-normalized",titSu,")"), las=1, outline=FALSE) 
        if(debug) {message(fxNa,"rpd16 .. dim quant: ", nrow(quant)," li and  ",ncol(quant)," cols; colnames : ",wrMisc::pasteC(colnames(quant))," ")}
        graphics::abline(h=round(stats::median(quant, na.rm=TRUE)) +c(-2:2), lty=2, col=grDevices::grey(0.6)) }
    } else {                                            # wrGraph and sm are available
      if(debug) {message(fxNa,"rpd17  print vioplotW "  )}
      wrGraph::vioplotW(log2(abund), tit=paste(tit,"(initial)",sep=" "), wex=wex, silent=silent, callFrom=fxNa) 
      graphics::abline(h=round(stats::median(log2(abund), na.rm=TRUE)) +c(-2:2), lty=2, col=grDevices::grey(0.6)) 
      ## now normalized
      if(debug) {message(fxNa,"rpd18  print norm vioplotW() ",identical(normalizeMeth,"none")," ou ", length(quant) <0)}
      if(!identical(normalizeMeth,"none") | length(quant) >0) {      
        if(debug) {message(fxNa,"rpd19  print vioplotW() for normalized")}
          wrGraph::vioplotW(quant, tit=paste(tit,", ",normalizeMeth,"-normalized",titSu), wex=wex, silent=silent, callFrom=fxNa) 
          graphics::abline(h=round(stats::median(quant, na.rm=TRUE)) +c(-2:2), lty=2, col=grDevices::grey(0.6)) }    
    }
    on.exit(graphics::par(oparMar)) }   # restaure old settings
  ## meta-data
  notes <- c(inpFile=paFi, qmethod="ProteomeDiscoverer", qMethVersion=if(length(infoDat) >0) unique(infoDat$Software.Revision) else NA,
  	rawFilePath= if(length(infoDat) >0) infoDat$File.Name[1] else NA, normalizeMeth=normalizeMeth, call=match.call(), 
    created=as.character(Sys.time()), wrProteo.version=utils::packageVersion("wrProteo"), machine=Sys.info()["nodename"])
  ## final output
  if(isTRUE(separateAnnot)) list(raw=abund, quant=quant, annot=annot, counts=counts, notes=notes) else data.frame(quant,annot)  
}
    
