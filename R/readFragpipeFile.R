#' Read Tabulated Files Exported by FragPipe At Protein Level
#'
#' This function allos importing protein identification and quantification results form \href{https://fragpipe.nesvilab.org/}{Fragpipe}
#' which were previously exported as tabulated text (tsv). Quantification data and other relevant information will be extracted similar like the other import-functions from this package.
#' The final output is a list containing the elements: \code{$annot}, \code{$raw} and \code{$quant}, or a data.frame with the quantication data and a part of the annotation if argument \code{separateAnnot=FALSE}.
#'
#' @details
#' This function has been developed using Fragpipe versions 18.0 and 19.0.
#'
#' Using the argument \code{suplAnnotFile} it is possible to specify a specific file (or search for default file) to read for extracting file-names as sample-names and other experiment related information.
#'
#' @param fileName (character) name of file to be read
#' @param path (character) path of file to be read
#' @param normalizeMeth (character) normalization method, defaults to \code{median}, for more details see \code{\link[wrMisc]{normalizeThis}})
#' @param sampleNames (character) new column-names for quantification data (ProteomeDiscoverer does not automatically use file-names from spectra); this argument has priority over \code{suplAnnotFile}
#' @param read0asNA (logical) decide if initial quntifications at 0 should be transformed to NA (thus avoid -Inf in log2 results)
#' @param quantCol (character or integer) exact col-names, or if length=1 content of \code{quantCol} will be used as pattern to search among column-names for $quant using \code{grep}
#' @param refLi (character or integer) custom specify which line of data is main species, if character (eg 'mainSpe'), the column 'SpecType' in $annot will be searched for exact match of the (single) term given
#' @param separateAnnot (logical) if \code{TRUE} output will be organized as list with \code{$annot}, \code{$abund} for initial/raw abundance values and \code{$quant} with final log2 (normalized) quantitations
#' @param annotCol (character) column names to be read/extracted for the annotation section (default  c("Accession","Description","Gene","Contaminant","Sum.PEP.Score","Coverage....","X..Peptides","X..PSMs","X..Unique.Peptides", "X..AAs","MW..kDa.") )
#' @param FDRCol (list) optional indication to search for protein FDR information
#' @param tit (character) custom title to plot
#' @param graphTit (character) depreciated custom title to plot, please use 'tit'
#' @param wex (integer) relative expansion factor of the violin-plot (will be passed to \code{\link[wrGraph]{vioplotW}})
#' @param specPref (character or list) define characteristic text for recognizing (main) groups of species (1st for comtaminants - will be marked as 'conta', 2nd for main species- marked as 'mainSpe',
#'  and optional following ones for supplemental tags/species - maked as 'species2','species3',...);
#'  if list and list-element has multiple values they will be used for exact matching of accessions (ie 2nd of argument \code{annotCol})
#' @param gr (character or factor) custom defined pattern of replicate association, will override final grouping of replicates from \code{sdrf} and/or \code{suplAnnotFile} (if provided)   \code{}
#' @param sdrf (character, list or data.frame) optional extraction and adding of experimenal meta-data: if character, this may be the ID at ProteomeExchange. Besides, the output from \code{readSdrf} or a list from \code{defineSamples} may be provided; if \code{gr} is provided, it gets priority for grouping of replicates
#' @param suplAnnotFile (logical or character) optional reading of supplemental files; however, if \code{gr} is provided, \code{gr} gets priority for grouping of replicates;
#'  if \code{character} the respective file-name (relative or absolute path)
#' @param groupPref (list) additional parameters for interpreting meta-data to identify structure of groups (replicates), will be passed to \code{readSampleMetaData}.
#'   May contain \code{lowNumberOfGroups=FALSE} for automatically choosing a rather elevated number of groups if possible (defaults to low number of groups, ie higher number of samples per group)
#' @param plotGraph (logical or integer) optional plot of type vioplot of initial and normalized data (using \code{normalizeMeth}); if integer, it will be passed to \code{layout} when plotting
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a list with \code{$raw} (initial/raw abundance values), \code{$quant} with final normalized quantitations, \code{$annot}, \code{$counts} an array with number of peptides, \code{$quantNotes}
#'  and \code{$notes}; or if \code{separateAnnot=FALSE} the function returns a data.frame with annotation and quantitation only
#' @seealso \code{\link[utils]{read.table}}, \code{\link[wrMisc]{normalizeThis}}) , \code{\link{readMaxQuantFile}}, \code{\link{readProlineFile}}
#' @examples
#' FPproFi1 <- "tinyFragpipe1.tsv.gz"
#' path1 <- system.file("extdata", package="wrProteo")
#' ## let's define the main species and allow tagging some contaminants
#' specPref1 <- c(conta="conta|CON_|LYSC_CHICK", mainSpecies="MOUSE")
#' dataFP <- readFragpipeFile(path1, file=FPproFi1, specPref=specPref1, tit="Tiny Fragpipe Data")
#' summary(dataFP$quant)
#'
#' @export
readFragpipeFile <- function(fileName, path=NULL, normalizeMeth="median", sampleNames=NULL, read0asNA=TRUE, quantCol="Intensity$",
  annotCol=NULL,  refLi=NULL, separateAnnot=TRUE, FDRCol=list("Protein.Probability", lim=0.99),    # contamCol="Contaminant",
  groupPref=list(lowNumberOfGroups=TRUE), plotGraph=TRUE, tit="FragPipe", graphTit=NULL, wex=1.6, specPref=c(conta="CON_|LYSC_CHICK", mainSpecies="OS=Homo sapiens"),
  gr=NULL, sdrf=NULL, suplAnnotFile=FALSE, silent=FALSE, debug=FALSE, callFrom=NULL) {

  ## read Fragpipe exported txt
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readFragpipeFile")
  oparMar <- if(plotGraph) graphics::par("mar") else NULL       # only if figure might be drawn

  reqPa <- c("utils","wrMisc")
  chPa <- sapply(reqPa, requireNamespace, quietly=TRUE)
  if(any(!chPa)) stop("package(s) '",paste(reqPa[which(!chPa)], collapse="','"),"' not found ! Please install first from CRAN")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
         excluCol <- "^Abundances.Count"   # exclude this from quantifications columns
  cleanDescription <- TRUE          # clean 'Description' for artifacts of truncated text (tailing ';' etc)
  infoDat <- infoFi <- setupSd <- parametersD <- NULL        # initialize

  ## check if path & file exist
  msg <- "Invalid entry for 'fileName'"
  if(length(fileName) >1) { fileName <- fileName[1]
    if(!silent) message(fxNa," 'fileName' shoud be of length=1, using 1st value")
  } else { if(length(fileName) <1) stop(msg) else if(is.na(fileName) | nchar(fileName) <1) stop(msg)}
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
  if(!chFi) stop(" File ",fileName," was NOT found ",if(length(path) >0) paste(" in path ",path)," !")
  if(!grepl("\\.tsv$|\\.tsv\\.gz$", fileName)) message(fxNa,"Trouble ahead, expecting tabulated text file (the file'",fileName,"' might not be right format) !!")
  if(debug) message(fxNa,"rfp0a ..")

  ## note : reading sample-setup from 'suplAnnotFile' at this place won't allow comparing if number of  samples/columns corresponds to data; do after reading main data
  if(debug) message(fxNa,"rfp0 .. Ready to read", if(length(path) >0) c(" from path ",path[1])," the file  ",fileName[1])


  ## read (main) file
  ## future: look for fast reading of files
  tmp <- try(utils::read.delim(file.path(paFi), stringsAsFactors=FALSE), silent=TRUE)

  if(length(tmp) <1 | inherits(tmp, "try-error") | length(dim(tmp)) <2) {
    if(inherits(tmp, "try-error")) warning("Unable to read input file ('",paFi,"')!  (check if rights to read)") else {
      if(!silent) message(fxNa,"Content of  file '",paFi,"' seeps empty or non-conform !  Returning NULL; check if this is really a Fragpipe-file") }
    NULL
  } else {
    if(debug) { message(fxNa,"rfp1 .. dims of initial data : ", nrow(tmp)," li and ",ncol(tmp)," col "); rfp1 <- list(fileName=fileName,path=path,paFi=paFi,tmp=tmp,normalizeMeth=normalizeMeth,sampleNames=sampleNames,read0asNA=read0asNA,quantCol=quantCol,
      annotCol=annotCol,refLi=refLi,separateAnnot=separateAnnot,FDRCol=FDRCol   )}

    ## locate & extract annotation
    ## note : space (' ') in orig colnames are transformed to '.'
    if(length(annotCol) <1) annotCol <- c("Protein","Protein.ID","Entry.Name","Description","Gene","Organism", "Protein.Length","Protein.Existence","Protein.Probability","Top.Peptide.Probability")
    ## note cols 2-6 are part to common format wrProteo
    PSMCol <- "\\.Spectral\\.Count$"                   # pattern searching tag for PSM-data
    PepCol <- "Unique\\.Spectral\\.Count$"             # pattern searching tag for Number of peptides
    ## future option : lateron rename columns called as "Description" to annotCol[2]
    ## below use explicit colnames "Accession","Description", rename if tolower() fits

    .chColNa <- function(x, mat, renameTo=NULL, silent=FALSE, fxNa=NULL){
      ## check in 'matr' for column-name 'x', if required rename best hit (if no direct hit look using grep, then grep wo case); return corrected mat
      chX <- x %in% colnames(mat)
      if(all(chX)) {
        if(is.character(renameTo) & length(renameTo) ==1) colnames(mat)[match(x, colnames(mat))] <- renameTo   # juste simple rename (single col only)
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
    if(is.character(annotCol)) annotColNo <- match(annotCol, colnames(tmp))
    chNa <- is.na(annotColNo)
    if(any(chNa) & silent) message(fxNa,"Missing ",sum(chNa)," annotation columns:  ",wrMisc::pasteC(annotCol[chNa], quoteC="'"))
    ## rename to wrProteo format
    tmp <- .chColNa(annotCol[2], tmp, renameTo="Accession", silent=silent, fxNa=fxNa)           # rename 'Protein ID' to 'Accession'  (Uniprot ID)
    tmp <- .chColNa(annotCol[3], tmp, renameTo="EntryName", silent=silent, fxNa=fxNa)           # like  THOC2_MOUSE
    tmp <- .chColNa(annotCol[4], tmp, renameTo="Description", silent=silent, fxNa=fxNa)         # full (long) name

    annot <- cbind(Accession=tmp[,"Accession"], EntryName=tmp[,"EntryName"], GeneName=NA, Species=NA, Contam=NA, SpecType=NA,
      Description=tmp[,"Description"], tmp[,wrMisc::naOmit(annotColNo[-(1:6)])])   # may be better to name column 'species'
    if(debug) { message(fxNa,"rfp2 .. annotColNo : ", wrMisc::pasteC(annotColNo)); rfp2 <- list(annot=annot,annotCol=annotCol,tmp=tmp,specPref=specPref )}

    ## Species  (need to run before reparsing nadly parsed)
    if(!is.na(annotColNo[6])) { spec <- tmp[,annotColNo[6]]
      spec <- sub("^\ +|\ +$","",spec)          # remove heading or tailing (white) space
      chOX <- grep(" OX=", spec)
      if(length(chOX) >0) { OX <- sub(" OX=", "", spec[chOX])
        spec[chOX] <- sub(" OX=[[:digit:]]+[[:print:]]*","", spec[chOX])
        chO2 <- nchar(spec[chOX]) <3 & nchar(OX) >1
        if(any(chO2)) spec[chOX[which(chO2)]] <- OX[which(chO2)]    # use OX=.. in case no other information available
      }
      if(TRUE) spec <- sub(" \\([[:alpha:]][[:print:]]+\\).*", "", spec)   # remove ' (..)'
      annot[,"Species"] <- spec
    }

    ## look for not well parsed (use separator '|' as indicator)
    chPa <- grep("\\|", annot[,"Accession"])
    if(length(chPa) >0) {
      chSp <- grep(" ", annot[chPa,"Accession"])
      if(length(chSp) >0) {
        # extract species
        chOS <- grep("[[:print:]]+ OS=[[:alpha:]]", annot[chPa[chSp],"Accession"])
        if(length(chOS) >0) annot[chPa[chSp[chOS]],"Species"] <- sub(" [[:upper:]]{2}=.+","", sub("[[:print:]]+ OS=","", annot[chPa[chSp[chOS]],"Accession"]))   # extract species
        ## extract GeneName
        chGn <- grep("[[:print:]]+ GN=", annot[chPa[chSp],"Accession"])
        if(length(chGn) >0) annot[chPa[chSp[chGn]],"GeneName"] <- sub(" [[:upper:]]{2}=.+","", sub("[[:print:]]+ GN=","", annot[chPa[chSp[chGn]],"Accession"]))
        ## extract Description
        annot[chPa[chSp],"Description"] <-  sub(".*? ", "", sub(" [[:upper:]]{2}=.+","", annot[chPa[chSp],"Accession"]))
        ## extract EntryName (option 1)
        annot[chPa[chSp],"EntryName"] <- gsub(".*\\|","", sub(" .+","", annot[chPa,"Accession"]))
      } else {
        annot[chPa,"EntryName"] <- gsub(".*\\|","", annot[chPa,"Accession"])     ## extract EntryName (option 2)
      }
      ## extract Accession
      annot[chPa,"Accession"] <- sapply(strsplit(annot[chPa,"Accession"], "\\|"), function(x) if(length(x) >1) x[2] else NA)
    }


    ## clean 'Description' entries: remove tailing punctuation or open brackets (ie not closed) at end of (truncated) fasta header
    if(cleanDescription) {
      if(debug) { message(fxNa,"rfp3a") }
      annot[,"Description"] <- sub("[[:punct:]]+$","", sub("\\ +$", "", annot[,"Description"]))     # tailing ';' and/or tailing space
      annot[,"Description"] <- sub(" \\([[:alpha:]]*$", "", annot[,"Description"])                  # tailing (ie truncated) open '(xxx'
    }

    if(debug) { message(fxNa,"rfp3b"); rfp3b <- list(annot=annot,tmp=tmp,specPref=specPref,annotCol=annotCol,PSMCol=PSMCol,PepCol=PepCol) }

    if(debug) {message(fxNa,"rfp4 .. dim annot: ", nrow(annot)," li and  ",ncol(annot)," cols; colnames : ",wrMisc::pasteC(colnames(annot))," ")}
    .MultGrep <- function(pat, y) if(length(pat)==1) grep(pat, y) else unlist(sapply(pat, grep, y))  # (multiple) grep() when length of pattern 'pat' >0

    ## Contam
    if("Contaminant" %in% colnames(annot)) {                         # just in case there is a column called 'Contaminant' (so far not seen)
      useLi <- which[nchar(annot[,"Contaminant"]) >0 & !is.na(annot[,"Contaminant"])]
      if(length(useLi) >0) annot[useLi,"Contam"] <- toupper(gsub(" ","",annot[useLi,"Contaminant"]))}
    chConta <- grep("^contam", tmp[,annotCol[1]])                    # specific to  Fragpipe
    if(length(chConta) >0) annot[chConta,"Contam"] <- TRUE

    ## get more species annot;  separate multi-species (create columns 'Accession','GeneName','Species','SpecType')
    chSp <- is.na(annot[,"Species"]) | nchar(annot[,"Species"]) <2
    if(any(chSp)) { chSep <- grep("_", annot[which(chSp),"EntryName"])                       # look for eg 'TRY1_BOVIN'
      if(length(chSep) >0) { chSep <- which(chSp)[chSep]
        spe2 <- sub("[[:alnum:]]+_", "", annot[chSep,"EntryName"])
        if(debug) message(fxNa,"Recover Species name for ",length(chSep)," entries based on 'EntryName'")
        commonSpec <- .commonSpecies()
        chSp3 <- which(sub("^_","",commonSpec[,1]) %in% spe2)
        if(length(chSp3) >0) for(i in chSp3) annot[chSep,"Species"] <- commonSpec[i,2]
      }
      chSp <- is.na(annot[,"Species"]) | nchar(annot[,"Species"]) <2 }     # update
    if(debug) {message(fxNa,"rfp6d ..  "); rfp6d <- list(annot=annot,tmp=tmp,chSp=chSp,specPref=specPref,annotCol=annotCol,PSMCol=PSMCol,PepCol=PepCol)}

    ## look for tags from  specPref
    if(length(specPref) >0) {
      ## set annot[,"specPref"] according to specPref
      annot <- .extrSpecPref(specPref, annot, silent=silent, debug=debug, callFrom=fxNa)
    } else if(debug) message(fxNa,"Note: Argument 'specPref' not specifed (empty)")
    if(debug) {message(fxNa,"rfp6b ..  ")}

    if(!silent) {
      if(any(chSp, na.rm=TRUE) & !all(chSp)) message(fxNa,"Note: ",sum(chSp)," (out of ",nrow(tmp),") lines with unrecognized species")
      if(!all(chSp)) { tab <- table(annot[,"Species"])
        tab <- rbind(names(tab), paste0(": ",tab," ;  "))
        if(!silent) message(fxNa,"Count by 'specPref' : ",apply(tab, 2, paste)) }}             # all lines assigned
    if(debug) {message(fxNa,"rfp6e ..  ")}

    ## check for unique annot[,"Accession"]
    chDu <- duplicated(annot[,"Accession"], fromLast=FALSE)
    if(any(chDu)) { warning(fxNa," NOTE : ",sum(chDu)," entries have same '",annotCol[2],"' (ie Accession) - correcting to UNIQUE !")
      rownames(tmp) <- rownames(annot) <- wrMisc::correctToUnique(annot[,"Accession"], sep="_", atEnd=TRUE, callFrom=fxNa)
    } else { rownames(annot) <- rownames(tmp) <- annot[,"Accession"] }

    if(debug) { message(fxNa,"rfp7 .. dim annot ",nrow(annot)," and ",ncol(annot)); rfp7 <- list(annot=annot,tmp=tmp,annot=annot,specPref=specPref) }


    ## locate & extract abundance/quantitation data
    msg <- " CANNOT find ANY quantification columns"
    if(length(quantCol) >1) {
      ## explicit columns (for abundance/quantitation data)
      if(is.character(quantCol)) quantCol <- match(quantCol, colnames(tmp))
    } else {
      ## pattern search (for abundance/quantitation data)
      ## problem : extract 'xx1.Intensity' but NOT 'xx.MaxLFQ.Intensity'
      useMaxLFQItens <- FALSE
      quantColIni <- quantCol <- grep(quantCol, colnames(tmp))
      chLFQ <- grep("MaxLFQ\\.", colnames(tmp)[quantCol])
      if(length(chLFQ) >0) { if(!silent & length(chLFQ)==length(quantCol)) message(fxNa,"All quantification columns are MaxLFQ !")
        if(length(chLFQ) < length(quantCol)) quantCol <- quantCol[(if(useMaxLFQItens) 1 else -1) *chLFQ] else warning("No non-MaxLFQ data available, using MaxLFQ.Intensity instead !") }
    }
    if(length(quantCol) <1) stop(msg,"  ('",quantCol,"')")
    abund <- as.matrix(tmp[, quantCol])
    rownames(abund) <- annot[,"Accession"]
    if(debug) { message(fxNa,"rfp8 .. dim abund ",nrow(abund)," and ",ncol(abund)) }

    ## check & clean abundances

    ## add custom sample names (if provided)
    if(length(sampleNames) ==ncol(abund) & ncol(abund) >0) {
      if(debug) { message(fxNa,"rfp8b") }
      if(length(unique(sampleNames)) < length(sampleNames)) {
        if(!silent) message(fxNa,"Custom sample names not unique, correcting to unique")
        sampleNames <- wrMisc::correctToUnique(sampleNames, callFrom=fxNa) }
      colnames(abund) <- sampleNames
    }
    if(debug) { message(fxNa,"rfp8c")}

    ## (optional) filter by FDR  (so far use 1st of list where matches are found from argument FDRCol)
    if(length(FDRCol) >0) {
      if(FDRCol[[1]] %in% colnames(tmp)) {
        if(length(FDRCol[[2]]) >0 & is.numeric(FDRCol[[2]])) FdrLim <- FDRCol[[2]][1] else {
          if(!silent) message(fxNa,"No valid FDR limit found, using default 0.95 (ie 5% filter)")
          FdrLim <- 0.95 }
        rmLi <- which(as.numeric(tmp[,FDRCol[[1]]]) < FdrLim)    # default 5% 'FDR' filter
        if(length(rmLi) == nrow(abund)) warning(fxNa,"Omit FDR-filter; otherwise NO MORE LINES/proteins remaining !!!")  else {
          if(!silent) message(fxNa,"Removing ",length(rmLi)," lines/proteins removed as NOT passing protein identification filter at ",FdrLim, if(debug) "   rfp8d")
          abund <- abund[-rmLi,]
          if(length(dim(abund)) <2) abund <- matrix(abund, nrow=1, dimnames=list(rownames(annot)[-rmLi], names(abund)))
          annot <- if(nrow(abund) ==1) matrix(annot[-rmLi,], nrow=1, dimnames=list(rownames(abund), colnames(annot))) else annot[-rmLi,]
          tmp <- if(nrow(abund) ==1) matrix(tmp[-rmLi,], nrow=1, dimnames=list(rownames(abund), colnames(tmp))) else tmp[-rmLi,]
        }
      }
    }
    if(debug) { message(fxNa,"rfp11 .. length(FDRCol) ",length(FDRCol),"   dim annot ",nrow(annot)," and ",ncol(annot)); rfp11 <- list(annot=annot,tmp=tmp,abund=abund)}

    PSMCol <- "\\.Spectral\\.Count$"                   # pattern searching tag for PSM-data
    PepCol <- "Unique\\.Spectral\\.Count$"             # pattern searching tag for Number of peptides
    PSMColExcl <- "Total\\.Spectral\\.Count$"          # exclude this pattern searching tag for PSM
    usTy <- c("PSM", "UniquePeptides")

    ## optional/additional counting results (PSM, no of peptides)
    PSMExl <- grep(paste0("Combined",PSMCol), colnames(tmp))
    PepExl <- grep(paste0("Combined\\.",PepCol), colnames(tmp))
    PSMCol <- if(length(PSMCol) ==1) grep(PSMCol, colnames(tmp)) else NULL
    PepCol <- if(length(PepCol) ==1) grep(PepCol, colnames(tmp)) else NULL
    if(any(c(length(PSMExl), length(PSMColExcl)) >0)) PSMCol <- PSMCol[-which(PSMCol %in% c(PepCol, PSMExl, grep(PSMColExcl, colnames(tmp))))]    # remove unwanted columns
    if(length(PepExl) >0) PepCol <- PepCol[-which(PepCol %in% PepExl)]
    if(any(c(length(PSMCol), length(PepCol)) >0)) {
      counts <- array(NA, dim=c(nrow(abund), ncol(abund), length(usTy)), dimnames=list(rownames(abund),colnames(abund), usTy))
      if(length(PSMCol) >0) counts[,,"PSM"] <- as.matrix(tmp[,PSMCol])
      if(length(PepCol) >0) counts[,,"UniquePeptides"] <- as.matrix(tmp[,PepCol])
    } else counts <- NULL
    if(debug) {message(fxNa,"rfp12 .. ");
      rfp12 <- list(tmp=tmp,abund=abund,annot=annot,sdrf=sdrf, fileName=fileName,path=path,paFi=paFi,normalizeMeth=normalizeMeth,sampleNames=sampleNames,
            refLi=refLi,specPref=specPref,read0asNA=read0asNA,quantCol=quantCol,annotCol=annotCol,refLi=refLi,separateAnnot=separateAnnot,FDRCol=FDRCol,gr=gr) }


    ## check for reference for normalization
    refLiIni <- refLi
    if(is.character(refLi) & length(refLi)==1) { refLi <- which(annot[,"SpecType"]==refLi)
      if(length(refLi) <1 & !silent) message(fxNa,"Could not find any protein matching argument 'refLi', ignoring ...") else {
        if(!silent) message(fxNa,"Normalize using (custom) subset of ",length(refLi)," lines specified as '",refLiIni,"'")}}    # may be "mainSpe"
    if(length(refLi) <1) refLi <- 1:nrow(abund)

    ## set 0 values to NA (avoid -Inf at log2)
    if(!isFALSE(read0asNA)) { ch0 <- abund ==0
      if(any(ch0, na.rm=TRUE)) abund[which(ch0)] <- NA }

    ## take log2 & normalize
    quant <- if(utils::packageVersion("wrMisc") > "1.10") {
        try(wrMisc::normalizeThis(log2(abund), method=normalizeMeth, mode="additive", refLines=refLi, silent=silent, callFrom=fxNa), silent=TRUE)
      } else try(wrMisc::normalizeThis(log2(abund), method=normalizeMeth, refLines=refLi, silent=silent, callFrom=fxNa), silent=TRUE)       #
    if(debug) { message(fxNa,"rfp13 .. dim quant: ", nrow(quant)," li and  ",ncol(quant)," cols; colnames : ",wrMisc::pasteC(colnames(quant))," ")
      rfp13 <- list(tmp=tmp,quant=quant,abund=abund,annot=annot,sdrf=sdrf, fileName=fileName,path=path,paFi=paFi,normalizeMeth=normalizeMeth,sampleNames=sampleNames,
            refLi=refLi,specPref=specPref,read0asNA=read0asNA,quantCol=quantCol,annotCol=annotCol,refLi=refLi,separateAnnot=separateAnnot,FDRCol=FDRCol,gr=gr) }

    ### GROUPING OF REPLICATES AND SAMPLE META-DATA
    if(any((!isFALSE(suplAnnotFile) & length(suplAnnotFile) >0), length(sdrf) >0)) {
      setupSd <- readSampleMetaData(sdrf=sdrf, suplAnnotFile=suplAnnotFile, quantMeth="FP", path=path, abund=utils::head(quant), groupPref=groupPref, silent=silent, debug=debug, callFrom=fxNa)
    }
    if(debug) {message(fxNa,"rfp13b"); rfp13b <- list(tmp=tmp,quant=quant,abund=abund,setupSd=setupSd,annot=annot,sdrf=sdrf, fileName=fileName,path=path,paFi=paFi,normalizeMeth=normalizeMeth,sampleNames=sampleNames,
            specPref=specPref,read0asNA=read0asNA,quantCol=quantCol,annotCol=annotCol,refLi=refLi,separateAnnot=separateAnnot,FDRCol=FDRCol, gr=gr )}

    ## finish groups of replicates & annotation setupSd
    if(length(setupSd) >0 & length(setupSd$groups) <1) {                       # if nothing found/matching from sdrf & file, try getting sample-setup from colnames (use if at least 1 replicate found)
      if(length(gr) ==ncol(abund)) setupSd$groups <- gr else {
        if(debug) {message(fxNa,"rfp13e  Note: setupSd$groups is still empty ! ")}
        if(length(setupSd$lev) ==ncol(abund)) setupSd$groups <- setupSd$lev else {
          if(length(setupSd$lev) >0 & ncol(abund) >0) warning(fxNa," sdrf sample-metadata dot not fit (sdrf suggests ",length(setupSd$lev)," but real data consist of ",ncol(abund)," samples) !")
          ## try defining groups based on colnames
          if(debug) {message(fxNa,"rfp13f  Note: setupSd$lev is still empty/not matching !  .. try getting sample-setup from colnames");
            rfp13f <- list(tmp=tmp,quant=quant,abund=abund,setupSd=setupSd,annot=annot,sdrf=sdrf, ileName=fileName,path=path,paFi=paFi,tmp=tmp,normalizeMeth=normalizeMeth,sampleNames=sampleNames,
            read0asNA=read0asNA,quantCol=quantCol,annotCol=annotCol,refLi=refLi,separateAnnot=separateAnnot,FDRCol=FDRCol )}
          delPat <- "_[[:digit:]]+$|\\ [[:digit:]]+$|\\-[[:digit:]]+$"       # remove enumerators, ie trailing numbers after separator
          grou <- sub(delPat,"", colnames(abund))
          if(length(unique(grou)) >1 & length(unique(grou)) < ncol(abund)) setupSd$groups <- grou else {
            grou <- sub("[[:digit:]]+$","", colnames(abund))
            if(length(unique(grou)) >1 & length(unique(grou)) < ncol(abund)) setupSd$groups <- grou }
        } }
    }
    if(debug) {message(fxNa,"rfp13g"); rfp13g <- list(tmp=tmp,quant=quant,abund=abund,setupSd=setupSd,annot=annot,sdrf=sdrf, ileName=fileName,path=path,paFi=paFi,tmp=tmp,normalizeMeth=normalizeMeth,sampleNames=sampleNames,
            read0asNA=read0asNA,quantCol=quantCol,annotCol=annotCol,refLi=refLi,separateAnnot=separateAnnot,FDRCol=FDRCol )}

    ## finish sample-names: use file-names from meta-data if no custom 'sampleNames' furnished
    ## One more check for colnames & sampleNames
    if(length(sampleNames) <1) sampleNames <- wrMisc::trimRedundText(
       if(length(setupSd$sdrfDat$comment.data.file.) ==ncol(abund)) setupSd$sdrfDat$comment.data.file. else colnames(abund), minNchar=2, spaceElim=TRUE, silent=silent, callFrom=fxNa, debug=debug)

    if(debug) {message(fxNa,"rfp13h"); rfp13h <- list(tmp=tmp,quant=quant,abund=abund,setupSd=setupSd,annot=annot,sdrf=sdrf, fileName=fileName,path=path,paFi=paFi,tmp=tmp,normalizeMeth=normalizeMeth,sampleNames=sampleNames,
            read0asNA=read0asNA,quantCol=quantCol,annotCol=annotCol,refLi=refLi,separateAnnot=separateAnnot,FDRCol=FDRCol )}

    if(length(sampleNames) == ncol(quant)) {
      colnames(quant) <- colnames(abund) <- sampleNames     # use custom provided sampleNames
      if(length(dim(counts)) >1 & length(counts) >0) colnames(counts) <- sampleNames
    } else {
      if(debug) { message(fxNa,"rfp13i") }
      if("annotBySoft$File.Name" %in% names(setupSd)) {
        sampleNames <- if(length(setupSd$annotBySoft$File.Name)==ncol(quant)) {    # try taking sampleNames from sdrf
          basename(sub("\\.raw$|\\.Raw$|\\.RAW$","", setupSd$annotBySoft$File.Name)) } else colnames(abund)     # otherwise remain with colnames from main file
        sampleNames <- wrMisc::trimRedundText(sampleNames, minNchar=2, spaceElim=TRUE, silent=silent, callFrom=fxNa, debug=debug)
        colnames(quant) <- colnames(abund) <- sampleNames
      } else  colnames(quant) <- sampleNames <- colnames(abund)
      if(length(dim(counts)) >1 & length(counts) >0) colnames(counts) <- sampleNames

    }
    if(debug) {message(fxNa,"Read sample-meta data, rfp14"); rfp14 <- list(sdrf=sdrf,suplAnnotFile=suplAnnotFile,abund=abund, quant=quant,refLi=refLi,annot=annot,setupSd=setupSd,sampleNames=sampleNames)}

    ## plot distribution of intensities
    custLay <- NULL
    if(is.numeric(plotGraph) & length(plotGraph) >0) {custLay <- as.integer(plotGraph); plotGraph <- TRUE} else {
      if(!isTRUE(plotGraph)) plotGraph <- FALSE}
    if(plotGraph) {
      if(debug) {message(fxNa,"rfp15 .. length custLay ", length(custLay) )}
      if(length(custLay) >0) graphics::layout(custLay) else {if(!identical(normalizeMeth,"none") & length(quant) >0) graphics::layout(1:2)}
      ch1 <- try(graphics::par(mar=c(3, 3, 3, 1)), silent=TRUE)                          # mar: bot,le,top,ri
      if(inherits(ch1, "try-error")) {plotGraph <- FALSE; message(fxNa,"UNABLE to graphical paramters !! Abandoning plot ! check current graphics device ...")} }      
    if(plotGraph) {
      if(length(graphTit) >0) message(fxNa,"Argument 'graphTit' is depreciated, please rather use 'tit'")
      if(is.null(tit) & !is.null(graphTit)) tit <- graphTit     # for derpreciated argument
      if(is.null(tit)) tit <- "Fragpipe quantification "
      chGr <- try(find.package("wrGraph"), silent=TRUE)
      chSm <- try(find.package("sm"), silent=TRUE)
      misPa <- c(inherits(chGr, "try-error"), inherits(chSm, "try-error"))
      titSu <- if(length(refLi) >0) paste0(c(" by ",if(length(refLiIni) >1) c(length(refLi)," selected lines") else c("'",refLiIni,"'")),collapse="")  else NULL
      if(debug) { message(fxNa,"rfp16 .. misPa ", wrMisc::pasteC(misPa,quoteC="'") )}
      if(any(misPa, na.rm=TRUE)) {
        if(!silent) message(fxNa,"Missing package ",wrMisc::pasteC(c("wrGraph","sm")[which(misPa)],quoteC="'")," for drawing vioplots")
        ## wrGraph not available : simple boxplot
        graphics::boxplot(log2(abund), main=paste(tit,"(initial)",sep=" "), las=1, outline=FALSE)
        graphics::abline(h=round(log2(stats::median(abund,na.rm=TRUE))) +c(-2:2)*3, lty=2, col=grDevices::grey(0.6))
        ## plot normalized
        if(identical(normalizeMeth, "none") | length(quant) <0) {
          graphics::boxplot(quant, main=paste(tit," (",normalizeMeth,"-normalized",titSu,")"), las=1, outline=FALSE)
          if(debug) {message(fxNa,"rfp17 .. dim quant: ", nrow(quant)," li and  ",ncol(quant)," cols; colnames : ",wrMisc::pasteC(colnames(quant))," ")}
          graphics::abline(h=round(stats::median(quant, na.rm=TRUE)) +c(-2:2)*3, lty=2, col=grDevices::grey(0.6)) }
      } else {                                            # wrGraph and sm are available
        if(debug) {message(fxNa,"rfp18  print vioplotW "  )}
        ch1 <- try(wrGraph::vioplotW(log2(abund), tit=paste(tit,"(initial)",sep=" "), wex=wex, silent=silent, callFrom=fxNa), silent=TRUE)
        if(inherits(ch1, "try-error")) {plotGraph <- FALSE; message(fxNa,"UNABLE to plot vioplotW !!")
        } else graphics::abline(h=round(stats::median(log2(abund), na.rm=TRUE)) +c(-2:2)*3, lty=2, col=grDevices::grey(0.6))
        ## now normalized
        if(debug) {message(fxNa,"rfp19  print norm vioplotW() ",identical(normalizeMeth,"none")," ou ", length(quant) <0)}
        if(!identical(normalizeMeth,"none") | length(quant) >0) {
          if(debug) {message(fxNa,"rfp20  print vioplotW() for normalized")}
            ch1 <- try(wrGraph::vioplotW(quant, tit=paste(tit,", ",normalizeMeth,"-normalized",titSu), wex=wex, silent=silent, callFrom=fxNa), silent=TRUE)
            if(inherits(ch1, "try-error")) {plotGraph <- FALSE; message(fxNa,"UNABLE to plot vioplotW !!")
            } else graphics::abline(h=round(stats::median(quant, na.rm=TRUE)) +c(-2:2)*3, lty=2, col=grDevices::grey(0.6)) }
        on.exit(graphics::par(mar=oparMar)) }            # restaure old settings
    }

    ## meta-data
    notes <- c(inpFile=paFi, qmethod="Fragpipe", qMethVersion=if(length(infoDat) >0) unique(infoDat$Software.Revision) else NA,
    	rawFilePath= if(length(infoDat) >0) infoDat$File.Name[1] else NA, normalizeMeth=normalizeMeth, call=match.call(),
      created=as.character(Sys.time()), wrProteo.version=utils::packageVersion("wrProteo"), machine=Sys.info()["nodename"])
    ## final output
    if(isTRUE(separateAnnot)) list(raw=abund, quant=quant, annot=annot, counts=counts, sampleSetup=setupSd, quantNotes=parametersD, notes=notes) else data.frame(quant,annot) }
}

