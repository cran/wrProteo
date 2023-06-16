#' Read Tabulated Files Exported by DiaNN At Peptide Level
#'
#' This function allows importing peptide identification and quantification results from \href{https://github.com/vdemichev/DiaNN}{DiaNN}, see also \doi{10.1038/s41592-019-0638-x}{Demichev et al, 2020}.
#' Data should be exported as tabulated text (tsv) to allow import by thus function.
#' Quantification data and other relevant information will be extracted similar like the other import-functions from this package.
#' The final output is a list containing the elements: \code{$annot}, \code{$raw} and \code{$quant}, or a data.frame with the quantication data and a part of the annotation if argument \code{separateAnnot=FALSE}.
#'
#' @details
#' This function has been developed using DiaNN version 1.8.x.
#'
#' Using the argument \code{suplAnnotFile} it is possible to specify a specific file (or search for default file) to read for extracting file-names as sample-names and other experiment related information.
#'
#' @param fileName (character) name of file to be read
#' @param path (character) path of file to be read
#' @param normalizeMeth (character) normalization method, defaults to \code{median}, for more details see \code{\link[wrMisc]{normalizeThis}})
#' @param sampleNames (character) custom column-names for quantification data; this argument has priority over \code{suplAnnotFile}
#' @param read0asNA (logical) decide if initial quntifications at 0 should be transformed to NA (thus avoid -Inf in log2 results)
#' @param quantCol (character or integer) exact col-names, or if length=1 content of \code{quantCol} will be used as pattern to search among column-names for $quant using \code{grep}
#' @param refLi (character or integer) custom specify which line of data is main species, if character (eg 'mainSpe'), the column 'SpecType' in $annot will be searched for exact match of the (single) term given
#' @param separateAnnot (logical) if \code{TRUE} output will be organized as list with \code{$annot}, \code{$abund} for initial/raw abundance values and \code{$quant} with final log2 (normalized) quantitations
#' @param annotCol (character) column names to be read/extracted for the annotation section (default  c("Accession","Description","Gene","Contaminant","Sum.PEP.Score","Coverage....","X..Peptides","X..PSMs","X..Unique.Peptides", "X..AAs","MW..kDa.") )
#' @param FDRCol (list) - not used
#' @param wex (integer) relative expansion factor of the violin-plot (will be passed to \code{\link[wrGraph]{vioplotW}})
#' @param specPref (character or list) define characteristic text for recognizing (main) groups of species (1st for comtaminants - will be marked as 'conta', 2nd for main species- marked as 'mainSpe',
#'  and optional following ones for supplemental tags/species - maked as 'species2','species3',...);
#'  if list and list-element has multiple values they will be used for exact matching of accessions (ie 2nd of argument \code{annotCol})
#' @param gr (character or factor) custom defined pattern of replicate association, will override final grouping of replicates from \code{sdrf} and/or \code{suplAnnotFile} (if provided)   \code{}
#' @param sdrf (character, list or data.frame) optional extraction and adding of experimenal meta-data: if character, this may be the ID at ProteomeExchange,
#'   the second element may give futher indicatations for automatic organization of groups of replicates.
#'   Besides, the output from \code{readSdrf} or a list from \code{defineSamples} may be provided; if \code{gr} is provided, \code{gr} gets priority for grouping of replicates
#' @param suplAnnotFile (logical or character) optional reading of supplemental files; however, if \code{gr} is provided, \code{gr} gets priority for grouping of replicates;
#'  if \code{character} the respective file-name (relative or absolute path)
#' @param groupPref (list) additional parameters for interpreting meta-data to identify structure of groups (replicates), will be passed to \code{readSampleMetaData}.
#'   May contain \code{lowNumberOfGroups=FALSE} for automatically choosing a rather elevated number of groups if possible (defaults to low number of groups, ie higher number of samples per group)
#' @param plotGraph (logical or integer) optional plot of type vioplot of initial and normalized data (using \code{normalizeMeth}); if integer, it will be passed to \code{layout} when plotting
#' @param titGraph (character) custom title to plot of distribution of quantitation values
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a list with \code{$raw} (initial/raw abundance values), \code{$quant} with final normalized quantitations, \code{$annot}, \code{$counts} an array with number of peptides, \code{$quantNotes}
#'  and \code{$notes}; or if \code{separateAnnot=FALSE} the function returns a data.frame with annotation and quantitation only
#' @seealso \code{\link[utils]{read.table}}, \code{\link[wrMisc]{normalizeThis}}) , \code{\link{readMaxQuantFile}}, \code{\link{readProtDiscovFile}}, \code{\link{readProlineFile}}
#' @examples
#' diaNNFi1 <- "tinyDiaNN1.tsv.gz"
#' ## This file contains much less identifications than one may usually obtain
#' path1 <- system.file("extdata", package="wrProteo")
#' ## let's define the main species and allow tagging some contaminants
#' specPref1 <- c(conta="conta|CON_|LYSC_CHICK", mainSpecies="HUMAN")
#' dataNN <- readDiaNNFile(path1, file=diaNNFi1, specPref=specPref1, tit="Tiny DIA-NN Data")
#' summary(dataNN$quant)
#'
#' @export
readDiaNNPeptides <- function(fileName, path=NULL, normalizeMeth="median", sampleNames=NULL, read0asNA=TRUE, quantCol="\\.raw$",
  annotCol=NULL,  refLi=NULL, separateAnnot=TRUE, FDRCol=NULL,    # contamCol="Contaminant",
  groupPref=list(lowNumberOfGroups=TRUE), plotGraph=TRUE, titGraph="DiaNN", wex=1.6, specPref=c(conta="CON_|LYSC_CHICK", mainSpecies="OS=Homo sapiens"),
  gr=NULL, sdrf=NULL, suplAnnotFile=FALSE, silent=FALSE, debug=FALSE, callFrom=NULL) {

  ## read DiaNN exported txt
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readDiaNNPeptides")
  oparMar <- if(plotGraph) graphics::par("mar") else NULL       # only if figure might be drawn

  reqPa <- c("utils","wrMisc")
  chPa <- sapply(reqPa, requireNamespace, quietly=TRUE)
  if(any(!chPa)) stop("package(s) '",paste(reqPa[which(!chPa)], collapse="','"),"' not found ! Please install first from CRAN")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
         excluCol <- "^Abundances.Count"   # exclude this from quantifications columns
  cleanDescription <- TRUE          # clean 'Description' for artifacts of truncated text (tailing ';' etc)
  infoDat <- infoFi <- setupSd <- parametersD <- NULL        # initialize

  ## check if path & (tsv) file exist
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
  if(debug) message(fxNa,"rdnp0a ..")

  ## note : reading sample-setup from 'suplAnnotFile' at this place won't allow comparing if number of  samples/columns corresponds to data; do after reading main data
  if(debug) message(fxNa,"rdnp0 .. Ready to read", if(length(path) >0) c(" from path ",path[1])," the file  ",fileName[1])


  ## read (main) file
  ## future: look for fast reading of files
  tmp <- try(utils::read.delim(file.path(paFi), stringsAsFactors=FALSE), silent=TRUE)

  if(length(tmp) <1 || inherits(tmp, "try-error") || length(dim(tmp)) <2) {
    if(inherits(tmp, "try-error")) warning("Unable to read input file ('",paFi,"')!  (check if rights to read)") else {
      if(!silent) message(fxNa,"Content of  file '",paFi,"' seeps empty or non-conform !  Returning NULL; check if this is really a Fragpipe-file") }
    NULL
  } else {
    if(debug) { message(fxNa,"rdnp1 .. dims of initial data : ", nrow(tmp)," li and ",ncol(tmp)," col "); rdnp1 <- list(fileName=fileName,path=path,paFi=paFi,tmp=tmp,normalizeMeth=normalizeMeth,sampleNames=sampleNames,read0asNA=read0asNA,quantCol=quantCol,
      annotCol=annotCol,refLi=refLi,separateAnnot=separateAnnot,FDRCol=FDRCol   )}

    ## locate & extract annotation
    ## note : space (' ') in orig colnames are transformed to '.'
    if(length(annotCol) <1) annotCol <- c("Protein.Group","Protein.Ids","Protein.Names","Genes","First.Protein.Description","Proteotypic","Stripped.Sequence","Precursor.Id","Precursor.Charge")
    ## note : 'Precursor.Id' contains modifications (Unimod number) and terminal charge number

    ## check for essential colnames !
      ## 'Accesion' (eg "P00498")  .. missing
      ## 'Description' (eg "Cyclin-dependent kinase 1")   .. missing
      ## no PSM o spectral counts data in file
      ##
    if(is.character(annotCol)) annotColNo <- match(annotCol, colnames(tmp))
    chNa <- is.na(annotColNo)
    if(any(chNa) && silent) message(fxNa,"Missing ",sum(chNa)," annotation columns:  ",wrMisc::pasteC(annotCol[chNa], quoteC="'"))
    ## rename columns to wrProteo format
    annot <- cbind(Accession=NA, EntryName=tmp[,annotCol[3]], GeneName=tmp[,annotCol[4]], Species=NA, Contam=NA, SpecType=NA,
      Description=NA,  UniProtID=tmp[,annotCol[2]], EntryNamesAll=tmp[,annotCol[3]],  GeneNameAll=tmp[,annotCol[4]], tmp[,wrMisc::naOmit(annotColNo[c(5:length(annotCol))])])   # may be better to name column 'species'
    if(debug) { message(fxNa,"rdnp2 .. annotColNo : ", wrMisc::pasteC(annotColNo)); rdnp2 <- list(annot=annot,annotCol=annotCol,tmp=tmp,specPref=specPref )}

    ## 'EntryName' & 'GeneName' may contain multiple proteins, pick 1st
    chMult <- grep(";",annot[,2])
    if(length(chMult) >0) annot[,2] <- sub(";.+","", annot[,2])
    chMult <- grep(";",annot[,3])
    if(length(chMult) >0) annot[,3] <- sub(";.+","", annot[,3])

    ## Species  (need to run before reparsing badly parsed)
      ## is there anaything that can be done for  annotColNo[3], ie 'First.Protein.Description' ?? was all NA in data provided
    if(!is.na(annotColNo[3])) {
      chSp <- grep("^[[:alnum:]]+_[[:upper:]]", tmp[,annotColNo[3]])
      if(length(chSp) >0) { commonSpec <- .commonSpecies()
        spe2 <- sub("^[[:alnum:]]+_", "", tmp[chSp,annotColNo[3]])
        chSp3 <- which(sub("^_","", commonSpec[,1]) %in% spe2)
      if(length(chSp3) >0) for(i in chSp3) annot[chSp,"Species"] <- commonSpec[i,2]}
    }
    ## clean 'Description' entries: remove tailing punctuation or open brackets (ie not closed) at end of (truncated) fasta header - applicable ?
    if(debug) {message(fxNa,"rdnp6d ..  "); rdnp6d <- list(annot=annot,tmp=tmp,chSp=chSp,specPref=specPref,annotCol=annotCol)}


    ## look for tags from  specPref
    if(length(specPref) >0) {
      ## set annot[,"specPref"] according to specPref
      annot <- .extrSpecPref(specPref, annot, silent=silent, debug=debug, callFrom=fxNa)
    } else if(debug) message(fxNa,"Note: Argument 'specPref' not specifed (empty)")
    if(debug) {message(fxNa,"rdnp6b ..  ")}

    if(!silent) {
      if(any(chSp, na.rm=TRUE) && !all(chSp)) message(fxNa,"Note: ",sum(chSp)," (out of ",nrow(tmp),") lines with unrecognized species")
      if(!all(chSp)) { tab <- table(annot[,"Species"])
        tab <- rbind(names(tab), paste0(": ",tab," ;  "))
        if(!silent) message(fxNa,"Count by 'specPref' : ",apply(tab, 2, paste)) }}             # all lines assigned
    if(debug) {message(fxNa,"rdnp6e ..  ")}

    ## check for unique annot[,"Accession"] - not applicable
    if(debug) { message(fxNa,"rdnp7 .. dim annot ",nrow(annot)," and ",ncol(annot)); rdnp7 <- list(annot=annot,tmp=tmp,annot=annot,specPref=specPref) }


    ## locate & extract abundance/quantitation data
    msg <- " CANNOT find ANY quantification columns"
    if(length(quantCol) <1) quantCol <- "\\.raw$"

    if(length(quantCol) ==1) {
      ## pattern search (for abundance/quantitation data)
      if(is.character(quantCol)) quantCol <- grep(quantCol, tolower(colnames(tmp)))
    }
    if(length(quantCol) <1) stop(msg,"  ('",quantCol,"')")
    abund <- as.matrix(tmp[, quantCol])
    rownames(abund) <- annot[,"EntryName"]
    if(debug) { message(fxNa,"rdnp8 .. dim abund ",nrow(abund)," and ",ncol(abund)) }

    ## check & clean abundances
    ## add custom sample names (if provided)
    if(length(sampleNames) ==ncol(abund) && ncol(abund) >0) {
      if(debug) { message(fxNa,"Valid 'sampleNames' were provided   rdnp8b") }
      if(length(unique(sampleNames)) < length(sampleNames)) {
        if(!silent) message(fxNa,"Custom sample names not unique, correcting to unique")
        sampleNames <- wrMisc::correctToUnique(sampleNames, callFrom=fxNa) }
      colnames(abund) <- sampleNames
      rownames(abund) <- annot[,"Precursor.Id"]
    }
    counts <- NULL                                 # not available
    if(debug) { message(fxNa,"rdnp8c")}

    ## (optional) filter by FDR  (so far use 1st of list where matches are found from argument FDRCol) - not applicable ?
    if(debug) { message(fxNa,"rdnp11"); rdnp11 <- list(annot=annot,tmp=tmp,abund=abund)}
    if(debug) {message(fxNa,"rdnp12 .. ");
      rdnp12 <- list(tmp=tmp,abund=abund,annot=annot,sdrf=sdrf, fileName=fileName,path=path,paFi=paFi,normalizeMeth=normalizeMeth,sampleNames=sampleNames,
            refLi=refLi,specPref=specPref,read0asNA=read0asNA,quantCol=quantCol,annotCol=annotCol,refLi=refLi,separateAnnot=separateAnnot,FDRCol=FDRCol,gr=gr) }

    ## correct colnames from 'Pathabc.raw' to 'abc'
    colnames(abund) <- wrMisc::.trimLeft(sub("\\.raw$|\\.RAW$","", colnames(abund)), silent=silent, debug=debug, callFrom=fxNa)

    ## check for reference for normalization
    refLiIni <- refLi
    if(is.character(refLi) && length(refLi)==1) {
      refLi <- which(annot[,"SpecType"]==refLi)
      if(length(refLi) <1 ) { refLi <- 1:nrow(abund)
        if(!silent) message(fxNa,"Could not find any peptides matching argument 'refLi=",refLiIni,"', ignoring ...")
      } else {
        if(!silent) message(fxNa,"Normalize using (custom) subset of ",length(refLi)," lines specified as '",refLiIni,"'")}}    # may be "mainSpe"

    ## set 0 values to NA (avoid -Inf at log2)
    if(!isFALSE(read0asNA)) { ch0 <- abund ==0
      if(any(ch0, na.rm=TRUE)) abund[which(ch0)] <- NA }

    ## take log2 & normalize
    quant <- try(wrMisc::normalizeThis(log2(abund), method=normalizeMeth, mode="additive", refLines=refLi, silent=silent, callFrom=fxNa), silent=TRUE)
    if(debug) { message(fxNa,"rdnp13 .. dim quant: ", nrow(quant)," li and  ",ncol(quant)," cols; colnames : ",wrMisc::pasteC(colnames(quant))," ")
      rdnp13 <- list(tmp=tmp,quant=quant,abund=abund,annot=annot,sdrf=sdrf, fileName=fileName,path=path,paFi=paFi,normalizeMeth=normalizeMeth,sampleNames=sampleNames,groupPref=groupPref,
            refLi=refLi,refLiIni=refLiIni,specPref=specPref,read0asNA=read0asNA,quantCol=quantCol,annotCol=annotCol,separateAnnot=separateAnnot,FDRCol=FDRCol,gr=gr,silent=silent,debug=debug) }

    ### GROUPING OF REPLICATES AND SAMPLE META-DATA
    if(length(suplAnnotFile) >0 || length(sdrf) >0) {
      #setupSd <- readSampleMetaData(sdrf=rdnp13$sdrf, suplAnnotFile=T, quantMeth="DN", path=rdnp13$path, abund=utils::head(rdnp13$quant), groupPref=rdnp13$groupPref, silent=rdnp13$silent, debug=rdnp13$debug, callFrom="readDN")
      setupSd <- readSampleMetaData(sdrf=sdrf, suplAnnotFile=separateAnnot, quantMeth="DN", path=path, abund=utils::head(quant), groupPref=groupPref, silent=silent, debug=debug, callFrom=fxNa)
    }
    if(debug) {message(fxNa,"rdnp13b .."); rdnp13b <- list()}

    ## finish groups of replicates & annotation setupSd
    setupSd <- .checkSetupGroups(abund=abund, setupSd=setupSd, gr=gr, sampleNames=sampleNames, quantMeth="DN", silent=silent, debug=debug, callFrom=fxNa)
    colnames(quant) <- colnames(abund) <- if(length(setupSd$sampleNames)==ncol(abund)) setupSd$sampleNames else setupSd$groups
    if(length(dim(counts)) >1 && length(counts) >0) colnames(counts) <- setupSd$sampleNames

    if(debug) {message(fxNa,"Read sample-meta data, rdnp14"); rdnp14 <- list(setupSd=setupSd, sdrf=sdrf, suplAnnotFile=suplAnnotFile,quant=quant,abund=abund,plotGraph=plotGraph)}

    ## main plotting of distribution of intensities
    custLay <- NULL
    if(is.numeric(plotGraph) && length(plotGraph) >0) {custLay <- as.integer(plotGraph); plotGraph <- TRUE} else {
        if(!isTRUE(plotGraph)) plotGraph <- FALSE}
    if(plotGraph) .plotQuantDistr(abund=abund, quant=quant, custLay=custLay, normalizeMeth=normalizeMeth, softNa="DiaNN",
      refLi=refLi, refLiIni=refLiIni, tit=titGraph, silent=silent, callFrom=fxNa, debug=debug)
    if(debug) {message(fxNa,"Read sample-meta data, rdnp15"); rdnp15 <- list()}


    ## meta-data
    notes <- c(inpFile=paFi, qmethod="DiaNN", qMethVersion=if(length(infoDat) >0) unique(infoDat$Software.Revision) else NA,
    	rawFilePath= if(length(infoDat) >0) infoDat$File.Name[1] else NA, normalizeMeth=normalizeMeth, call=match.call(),
      created=as.character(Sys.time()), wrProteo.version=utils::packageVersion("wrProteo"), machine=Sys.info()["nodename"])
    ## final output
    if(isTRUE(separateAnnot)) list(raw=abund, quant=quant, annot=annot, counts=counts, sampleSetup=setupSd, quantNotes=parametersD, notes=notes) else data.frame(quant,annot) }
}

