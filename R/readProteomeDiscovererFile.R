#' Read Tabulated Files Exported By ProteomeDiscoverer At Protein Level
#'
#' Protein identification and quantification results from
#' \href{https://www.thermofisher.com/order/catalog/product/OPTON-30812}{Thermo ProteomeDiscoverer}
#' which were exported as tabulated text can be imported and relevant information extracted.
#' 
#' @details
#' This function has been developed using Thermo ProteomeDiscoverer versions 2.2 to 2.5.
#' The format of resulting files at export also depends which columns are chosen as visible inside ProteomeDiscoverer and subsequently get chosen for export.
#' Using the argument \code{suplAnnotFile} it is possible to specify a specific file (or search for default file) to read for extracting file-names as sample-names and other experiment realted information.
#' If a column named \code{contamCol} is found, the data will be lateron filtered to remove all contaminants, set to \code{NULL} for keeping all contaminants.
#' 
#' The final output is a list containing as (main) elements: \code{$annot}, \code{$raw} and optional \code{$quant},
#' or returns data.frame with entire content of file if \code{separateAnnot=FALSE}.
#'
#' This function replaces the depreciated function \code{readProtDiscovFile} which will soon be retracted from this package.
#'
#' @param fileName (character) name of file to be read
#' @param path (character) path of file to be read
#' @param normalizeMeth (character) normalization method, defaults to \code{median}, for more details see \code{\link[wrMisc]{normalizeThis}})
#' @param sampleNames (character) custom column-names for quantification data (ProteomeDiscoverer does not automatically use file-names from spectra); this argument has priority over \code{suplAnnotFile}
#' @param read0asNA (logical) decide if initial quntifications at 0 should be transformed to NA
#' @param quantCol (character or integer) define ywhich columns should be extracted as quantitation data : The argument may be the exact column-names to be used, or if length=1 
#'   content of \code{quantCol} will be used as pattern to search among column-names for $quant using \code{grep};
#'   if \code{quantCol='allAfter_calc.pI'} all columns to the right of the column 'calc.pI' will be interpreted as quantitation data 
#'   (may be useful with files that have been manually edited before passing to wrProteo)
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
#'   the second & third elements may give futher indicatations for automatic organization of groups of replicates.
#'   Besides, the output from \code{readSdrf} or a list from \code{defineSamples} may be provided;
#'   if \code{gr} is provided, \code{gr} gets priority for grouping of replicates;
#'   if \code{sdrfOrder=TRUE} the output will be put in order of sdrf
#' @param suplAnnotFile (logical or character) optional reading of supplemental files produced by ProteomeDiscoverer; however, if \code{gr} is provided, \code{gr} gets priority for grouping of replicates;
#'  if \code{TRUE} defaults to file '*InputFiles.txt' (needed to match information of \code{sdrf}) which can be exported next to main quantitation results;
#'  if \code{character} the respective file-name (relative or absolute path)
#' @param groupPref (list) additional parameters for interpreting meta-data to identify structure of groups (replicates), will be passed to \code{readSampleMetaData}.
#'   May contain \code{lowNumberOfGroups=FALSE} for automatically choosing a rather elevated number of groups if possible (defaults to low number of groups, ie higher number of samples per group)
#'   May contain \code{chUnit} (logical or character) to be passed to \code{readSampleMetaData()} for (optional) adjustig of unit-prefixes in meta-data group labels, in case multiple different unit-prefixes 
#'   are used (eg '100pMol' and '1nMol').
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
#' dataPD <- readProteomeDiscovererFile(file=fiNa, path=path1, suplAnnotFile=FALSE)
#' summary(dataPD$quant)
#'
#' @export
readProteomeDiscovererFile <- function(fileName, path=NULL, normalizeMeth="median", sampleNames=NULL, read0asNA=TRUE, quantCol="^Abundance", annotCol=NULL, contamCol="Contaminant",
  refLi=NULL, separateAnnot=TRUE, FDRCol=list(c("^Protein.FDR.Confidence","High"), c("^Found.in.Sample.","High")), gr=NULL, sdrf=NULL, suplAnnotFile=TRUE,
  groupPref=list(lowNumberOfGroups=TRUE, chUnit=TRUE), specPref=c(conta="CON_|LYSC_CHICK", mainSpecies="OS=Homo sapiens"), plotGraph=TRUE, wex=1.6, titGraph="Proteome Discoverer",
  silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## read ProteomeDiscoverer exported txt
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readProteomeDiscovererFile")
  oparMar <- graphics::par("mar")                     # old margins, for rest after figure
  oparLayout <- graphics::par("mfcol")                # old layout, for rest after figure
  on.exit(graphics::par(mar=oparMar, mfcol=oparLayout))            # restore old mar settings

  reqPa <- c("utils","wrMisc")
  chPa <- sapply(reqPa, requireNamespace, quietly=TRUE)
  if(any(!chPa)) stop("package(s) '",paste(reqPa[which(!chPa)], collapse="','"),"' not found ! Please install first from CRAN")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  contamFilter <- TRUE             # filter contaminants away
  #excluCol <- c("^Abundance\\.Count","^Abundances\\.Count","^Abundance\\.Ratio","^Abundances\\.Ratio","^Abundance\\.Grouped","^Abundances\\.Grouped")   # exclude this from quantifications columns
  cleanDescription <- TRUE          # clean 'Description' for artifacts of truncated text (tailing ';' etc)
  infoDat <- infoFi <- setupSd <- parametersD <- NULL        # initialize

  .corPathW <- function(x) gsub("\\\\", "/", x)

  ## check if path & (tsv) file exist
  if(!grepl("\\.txt$|\\.txt\\.gz$", fileName)) message(fxNa,"Trouble ahead, expecting tabulated text file (the file'",fileName,"' might not be right format) !!")
  paFi <- wrMisc::checkFilePath(fileName, path, expectExt="txt", compressedOption=TRUE, stopIfNothing=TRUE, callFrom=fxNa, silent=silent, debug=debug)

  ## note : reading sample-setup from 'suplAnnotFile' at this place won't allow comparing if number of samples/columns corresponds to data; do after reading main data
  if(debug) {message(fxNa,"rpd0 .. Ready to read", if(length(path) >0) c(" from path ",path[1])," the file  ",fileName[1]); 
    rpd0 <- list(fileName=fileName,path=path,paFi=paFi,chPa=chPa) }

  ## read (main) file
  ## future: look for fast reading of files
  tmp <- try(utils::read.delim(paFi, stringsAsFactors=FALSE), silent=TRUE)

  if(length(tmp) <1 || inherits(tmp, "try-error") || length(dim(tmp)) <2) {
    if(inherits(tmp, "try-error")) warning("Unable to read input file ('",paFi,"')!  (check if rights to read)") else {
      if(!silent) message(fxNa,"Content of  file '",paFi,"' seeps empty or non-conform !  Returning NULL; check if this is really a ProteomeDiscoverer-file") }
    NULL
  } else {
    if(debug) { message(fxNa,"rpd1 ... dims of initial data : ", nrow(tmp)," li and ",ncol(tmp)," col ")
      rpd1 <- list(tmp=tmp,paFi=paFi,annotCol=annotCol,fileName=fileName) }

    ## locate & extract annotation
    ## default as R-friendly (convert standard cols later to this format)
    if(length(annotCol) <1) annotCol <- c("Accession","Description","Gene","GeneName","Gene.Name","Marked.as", "Number.of.Peptides","Number.of.PSMs","Number.of.Unique.Peptides","Number.of.AAs","Coverage.in.Percent")
    ## also  ??    "Exp.q.value.Combined","Sum.PEP.Score"
    if(debug) message(fxNa,"rpd1a")

    ## option for future: also extract column "MarkedAs"
    PSMCol <- "^Number\\.of\\.PSMs."             # pattern searching tag for PSM-data  (was previously .. .by.Search.Engine)
    PepCol <- "^Number\\.of\\.Peptides."         # pattern searching tag for Number of peptides (was preiously .. .by.Search.Engine)
    ## future option : lateron rename columns called as "Description" to annotCol[2]
    ## below use explicit colnames "Accession","Description", rename if tolower() fits

    .extrCol <- function(x, mat, asIndex=FALSE, notFound=NA, silent=FALSE, fxNa=NULL) {
      ## integrate to wrMisc ?
      ## look for column 'x' in matrix, extract index or content of column
      if(length(x) >1) x <- x[1]
      isNum <- grepl("^[[:digit:]]+$", x)
      if(isTRUE(isNum)) {
        x <- as.integer(x)
        if(x > ncol(mat)) stop(fxNa,"Invalid column number chosen (",x," but only ",ncol(mat)," available)")
        out <- if(isTRUE(asIndex)) x else matrix(mat[,x], ncol=1, dimnames=list(rownames(mat),colnames(mat)[x]))
      } else {
        ch1 <- colnames(mat) %in% x        # direct match
        if(any(ch1)) { if(sum(ch1) >1 && !silent) message(fxNa,"Note: ",sum(ch1)," columns at direct match, using 1st")
          out <- if(isTRUE(asIndex)) which(ch1)[1] else matrix(mat[,which(ch1)[1]], ncol=1, dimnames=list(rownames(mat),colnames(mat)[which(ch1)[1]]))
        } else {
          ch1 <- grep(x, mat)              # partial match (grep)
          if(length(ch1) >0) {  if(length(ch1) >1  && !silent)  message(fxNa,"Note: ",length(ch1)," columns found by grep, using 1st")
            out <- if(isTRUE(asIndex)) which(ch1)[1] else matrix(mat[,ch1[1]], ncol=1, dimnames=list(rownames(mat),colnames(mat)[ch1[1]]))
          } else {
            ch1 <- grep(tolower(x), tolower(mat))   # case-tolerant grep
            if(length(ch1) >0) {  if(length(ch1) >1  && !silent)  message(fxNa,"Note: ",length(ch1)," columns found by grep (case-independent), using 1st")
              out <- if(isTRUE(asIndex)) which(ch1)[1] else  matrix(mat[,ch1[1]], ncol=1, dimnames=list(rownames(mat),colnames(mat)[ch1[1]]))
            } else {
              if(grepl("^error|^stop", notFound)) {
                msg <- c(fxNa,"Could NOT find column '",x,"' !!\n  (available columns ",wrMisc::pasteC(colnames(mat), quoteC="'"),")")
                stop(msg)
              } else { if(is.null(notFound)) out <- NULL else {
                out <- if(isTRUE(asIndex)) NA else matrix(NA, ncol=1, dimnames=list(rownames(mat),x))
                if(!silent) if(grepl("^warning", notFound)) message(msg) else warning(msg)}}
            }
          }
        }
      }
      out }

     
    ## check for essential colnames !
    annot <- cbind(Accession=.extrCol(annotCol[1], tmp, silent=silent, fxNa=fxNa),    # 1st typically 'Accession'
      EntryName=NA, GeneName=NA, Description=.extrCol(annotCol[2], tmp, silent=silent, fxNa=fxNa), Species=NA, Contam=NA)  #, iniDescription=NA)
    
    ## address different naming/writing of GeneName
    GNco <- colnames(tmp) %in% c("Gene.Name","GeneName")
    if(any(GNco)) { annot[,"GeneName"] <- sub(" +$","",tmp[,GNco])       # also remove tailing space
      annotCol <- annotCol[-which(annotCol %in% c("Gene.Name","GeneName"))]
    }
    ## add other cols form annotCol
    if(debug) { message(fxNa,"rpd1b "); rpd1b <- list(tmp=tmp,annot=annot,paFi=paFi,annotCol=annotCol,fileName=fileName) }
    #notAnnCol <- c( "Accession","Description","Gene","Number.of.Peptides","Number.of.PSMs","Number.of.Unique.Peptides")
    if(length(annotCol) >2) { chSupCol <- match(annotCol[-(1:2)], colnames(tmp))
      if(sum(is.na(chSupCol)) < length(annotCol) -2) annot <- cbind(annot, tmp[,wrMisc::naOmit(chSupCol)])
      if(!silent) message(fxNa,"Adding supl annotation-columns ",wrMisc::pasteC(annotCol[wrMisc::naOmit(chSupCol)], quoteC="'")) }
    ## note 'EntryName' (eg 'UBE2C_HUMAN' not avail in PD)
    ## note 'GeneName' (eg 'UBE2C' can be extracted out of 'Description' after GN=
    ## extract GN
    if(debug) { message(fxNa,"rpd1c "); rpd1c <- list()  }

    ##  add more annot cols
    ## check for R-friendly
    Rfriendly <- !(length(grep("^X\\.\\.", colnames(tmp))) >0 || length(grep("Abundance\\.\\.F[[:digit:]]+\\.\\.", colnames(tmp)) >0))
    annotColNo <- NULL               # initialize for message rpd2
    if(length(annotCol) >2) {
      if(is.numeric(annotCol)) { message(fxNa,"Extraction by column index not yet finished")
      } else {
        annotColNo <- match(annotCol[3:length(annotCol)], colnames(tmp))
        if(!Rfriendly) {
          ## assume that annotation is in first 19 columns (example until 16) !!
          maxNCol <- min(ncol(tmp), 19)
          ## covert standard output to R-friendly
          colnames(tmp[1:maxNCol]) <-  sub("^X\\.\\.","Number.of.", colnames(tmp[1:maxNCol]))  #   sub("\\.\\.\\.$",".in.Percent",
          ## define table of specific terms to substitute ..
          subst=cbind(std=c("Exp.q.value.Combined","Coverage....","MW..kDa.","calc..pI"),
            rfriendly=c("Exp..q.value..Combined","Coverage.in.Percent","MW.in.kDa","calc.pI"))
          chSubst <- match(subst[,1], colnames(tmp)[1:maxNCol])
          if(any(!is.na(chSubst))) colnames(tmp)[wrMisc::naOmit(chSubst)] <- subst[wrMisc::naOmit(match(colnames(tmp)[1:maxNCol], subst[,1])), 2]
        }
        ## extract contam separately  & remove from annotCol
        if(contamCol %in% colnames(tmp)) annot[,"Contam"] <- as.logical(gsub(" ","",tmp[,contamCol]))
        annotCol <- annotCol[-which(annotCol %in% contamCol)]
      }
    }
    if(length(annotCol) >2) {
      annotColNo <- wrMisc::naOmit(match(annotCol[-(1:2)], colnames(tmp)))
      if(length(annotColNo) >0) annot <- cbind(annot, tmp[,annotColNo])
    }
    if(debug) { message(fxNa,"rpd2 .. annotColNo : ", wrMisc::pasteC(annotColNo))  #,"     contamCol : ",wrMisc::pasteC(contamCol)," ")
      rpd2 <- list(tmp=tmp,annot=annot,annotCol=annotCol,fileName=fileName,Rfriendly=Rfriendly)}  # PSMCol=PSMCol,PepCol=PepCol,

    ## extract GN from  'Description'  (like ...GN=LIPK )
    chGN <- grep(" GN=[[:alpha:]]", annot[,"Description"])
    if(length(chGN) >0)  annot[chGN,"GeneName"] <- sub(" [[:upper:]]{2}=.*","", sub(".* GN=","", annot[chGN,"Description"]))

    ## extract Species from  'Description'
    chSpe <- grep(" OS=[[:alpha:]]", annot[,"Description"])
    if(length(chSpe) >0)  { annot[chSpe,"Species"] <- sub(" [[:upper:]]{2}=.*","", sub(".* OS=","", annot[chSpe,"Description"]))
      chStr <- grep("^[[:upper:]][[:lower:]]+ [[:lower:]]+ \\(.", annot[chSpe,"Species"])                           # check for add'l strain name
      if(length(chStr) >0) annot[chSpe[chStr],"Species"]  <- substr(annot[chSpe[chStr],"Species"], 1, nchar(annot[chSpe[chStr],"Species"]) -2 -nchar(sub("^[[:upper:]][[:lower:]]+ [[:lower:]]+ \\(","", annot[chSpe[chStr],"Species"] )))       # remove strain no
    }

    ## clean 'Description' entries: remove tailing punctuation or open brackets (ie not closed) at end of (truncated) fasta header
    if(cleanDescription) {
      if(debug) { message(fxNa,"rpd3a"); rpd3a <- list() }
      chD <- grep(" [[:upper:]]{1,2}=", annot[,"Description"])
      if(length(chD) >0) annot[chD,"Description"] <- sub(" [[:upper:]]{1,2}=.*","", annot[chD,"Description"])   # remove all add'l fields (eg OS=... OX=... GN=...)
    }
    if(debug) {message(fxNa,"rpd4 .. dim annot: ", nrow(annot)," li and  ",ncol(annot)," cols; colnames : ",wrMisc::pasteC(colnames(annot))," ");
      rpd4 <- list(annot=annot,tmp=tmp,specPref=specPref,annotCol=annotCol,Rfriendly=Rfriendly,contamCol=contamCol)}

    ## identify user defined (main) groups of proteins (eg species, function, etc), create column 'SpecType'
    if(length(specPref) >0) {
      annot <- .extrSpecPref(specPref=specPref, annot=annot, useColumn=c("Species","EntryName","GeneName","Accession","Marked.as"), suplInp=tmp, silent=silent, debug=debug, callFrom=fxNa)   #"Majority.protein.IDs","Fasta.headers",
      chCont <- grep("^conta$|^contaminant", tolower(annot[,"SpecType"]))
      if(length(chCont) >0) {
        if(length(chCont)==nrow(tmp)) {warning(fxNa,"All proteins/peptides were marked as contaminants based on 'specPref' and will be removed !  Trouble ahead ?")
        } else if(!silent) message(fxNa,"Note : ",length(chCont)," proteins/peptide(s) were marked (in addition) as contaminants based on 'specPref'")
        annot[chCont,"Contam"] <- TRUE
      }
    }
   if(debug) {message(fxNa,"rpd5"); rpd5 <- list(annot=annot,tmp=tmp,specPref=specPref,annotCol=annotCol,Rfriendly=Rfriendly,contamCol=contamCol)}

    ## filter/remove contaminants unless in SpecTy
    if(any(as.logical(annot[,"Contam"]), na.rm=TRUE))  {    # filter contam
      filt1 <- if("SpecType" %in% colnames(annot)) {
        which( !as.logical(annot[,"Contam"]) | !grepl("^conta", tolower(annot[,"SpecType"])))     ## do not filter away proteins/peptides anotated by SpecTy (exept maked as 'conta')
      } else which(!as.logical(annot[,"Contam"]))
      if(length(filt1)==0) {
        warning(fxNa,"All lines will be removed based on contaminant-filter (column '",contamCol,"') - TROUBLE AHEAD ?")
      } else if(!silent) message(fxNa,"Note : ",nrow(tmp) - length(filt1)," contaminant protein(s)/peptide(s) will be removed, ",length(filt1)," remain")
      tmp <- tmp[filt1,]
      annot <- annot[filt1,]
    }
    if(debug) {message(fxNa,"rpd6b ..  "); rpd6b <- list()}

    ## report by species
    if(!silent) { chSp <- is.na(annot[,"Species"])
      if(!all(chSp)) { tab <- table(annot[,"Species"])
        if(any(chSp)) tab <- c(tab, na=sum(chSp))
        tab <- rbind(names(tab),": ",tab," ;  ")
        tab[length(tab)] <- ""       # no separator at last position
        message(fxNa,"Count by 'Species' : ",apply(tab, 2, paste)) }}             # all lines assigned

    if(debug) {message(fxNa,"rpd7 ..  "); rpd7 <- list(annot=annot,specPref=specPref,chSp=chSp,tmp=tmp,quantCol=quantCol )}

    extraQuantColCheck <- TRUE      # avoid  'Abundances.Count.' and 'Abundances.Normalized.' when pattern searching for quant-columns
    ## locate & extract abundance/quantitation data
    msg <- " CANNOT find ANY quantification columns"
    if(length(quantCol) <1) {
      ## default pattern search (for abundance/quantitation data)
      quantCol <- "^Abundance"
      if(debug) message(fxNa,"Setting argument 'quantCol' to '^Abundance'")
    }
    ## construct  quantColInd for index to use
    if(length(quantCol) >1) {
      ## explicit columns (for abundance/quantitation data)
      ch1 <- match(quantCol, colnames(tmp))          # look for direct match (rarley ok)
      chNA <- is.na(ch1)
          if(debug) {message(fxNa,"rpd7a ..  "); rpd7a <- list(annot=annot,specPref=specPref,chSp=chSp,tmp=tmp,quantCol=quantCol,ch1=ch1,chNA=chNA )}
      if(any(chNA)) {
        ch2 <- chNA & grepl("^[[:digit:]]", quantCol)
        if(any(ch2)) { quantCol[which(ch2)] <- paste0("X", quantCol[which(ch2)])
          ch1 <- match(quantCol, colnames(tmp))          # (update) look for direct match (rarley ok) x
          chNA <- is.na(ch1)
          if(debug) message(fxNa,"Trying to recuperate ",sum(ch2)," quantCol-names due to start with digits")}
      }     
      if(all(chNA)) {                          # no direct match found
        ch1 <- match(paste0(quantCol,".F.Sample"), sub("\\.F[[:digit:]]+\\.Sample",".F.Sample",colnames(tmp)))    # match to simplified "Cancer01.F.Sample" instead of "Cancer01.F1.Sample"
        chNA <- is.na(ch1)
        if(all(chNA)) {                     # no composed match found
          ## run grep on each indiv quantCol
          message(fxNa,"grep on each indiv quantCol - Not yet developed")
          ## develop further ?
          if(debug) {message(fxNa,"rpd7b ..  "); rpd7b <- list(annot=annot,specPref=specPref,chSp=chSp,tmp=tmp,quantCol=quantCol,ch1=ch1,chNA=chNA )}
          stop(fxNa,"None of samples found !!")
        }
      }
      if(any(chNA)) { message(fxNa,"Note : ",sum(chNA)," out of ",length(chNA)," samples NOT found !")
        ch1 <- ch1[which(!chNA)]
      }
      quantColInd <- ch1
    } else {    ## single value => use as search pattern
      if(debug) {message(fxNa,"rpd7c"); rpd7c <- list(annot=annot,specPref=specPref,chSp=chSp,tmp=tmp,quantCol=quantCol) }
      if(identical("allAfter_calc.pI", quantCol)) {
        calcPIcol <- c("calc.pI","calc..pI","calc pI","calc_pI")    # possible variations of 'calc.pI'
        chCol <- calcPIcol %in% colnames(tmp)
        endCol <- c("Number.of.Protein.Groups","Modifications")
        chEnd <- colnames(tmp) %in% endCol
        lastCol <- if(any(chEnd)) min(which(chEnd)) -1 else ncol(tmp)
        if(any(chCol, na.rm=TRUE)) {
          quantColInd <- (which(colnames(tmp)==calcPIcol[which(chCol)[1]]) +1) : lastCol
        } else {quantColInd <- NULL; stop(fxNa,"Cannot find column called 'calc_pI' (or 'cal.pI); don't know which columns to choose as quantitation data !")}
      } else {
        ## need to avoid "Abundances.Normalized.F1.Sample" or 'Abundances.Count.*'
        ch1 <- grep(paste0(if(substr(quantCol,1,1) !="^") "^",quantCol,"\\.F[[:digit:]]+\\.Sample"), colnames(tmp))      # construct pattern code to  "Abundance.F1.Sample"
        if(length(ch1) <1) {
          ch1 <- grep(paste0(if(substr(quantCol,1,1) !="^") "^",quantCol,"s\\.F[[:digit:]]+\\.Sample"), colnames(tmp))   # construct pattern code to  "Abundances.F1.Sample"
          if(length(ch1) <1) {
            ch1 <- grep(paste0(if(substr(quantCol,1,1) !="^") "^",quantCol,"s{0,1}\\.F[[:digit:]]"), colnames(tmp))      # construct pattern code to  "Abundances.F1"
            ## challange : one may pick too many columns (type "Abundances.Normalized.F1.Sample" or 'Abundances.Count.*')
            if(length(ch1) <1) {
              qC1 <- paste0(sub("\\\\.","", sub("S$","", quantCol)),"\\.S")
              ch1 <- grep(qC1, colnames(tmp))      # construct pattern code to  "Abundance.S"
              excluPat <- "\\.Normalized|\\.Count|\\.Ratio|\\.Grouped"
              if(length(ch1) >0) { ch2 <- grepl(excluPat, colnames(tmp)[ch1])    # avoid  '.Normalized' (if not concerning all)
                if(length(ch2) >0 && !all(ch2)) ch1 <- ch1[which(!ch2)]
              } else {
                ch1 <- grep(quantCol, colnames(tmp))      # use directly as pattern (may find too many)
                if(length(ch1) >0) { ch2 <- grepl(excluPat, colnames(tmp)[ch1])    # avoid  '.Normalized' (if not concerning all)
                  if(length(ch2) >0 && !all(ch2)) ch1 <- ch1[which(!ch2)]
                } else {
                  warning(fxNa,"Unable to find any matches to '",quantCol,"' !") }
              }
            }
          }
        if(debug) message(fxNa,"Found ",length(ch1)," quantitation-columns", if(length(ch1) >0) c(" (eg ",wrMisc::pasteC(colnames(tmp)[utils::head(ch1)], quoteC="'"),")"))
        }
        quantColInd <- ch1
      }
    }
    
    if(length(quantColInd) <1) stop(msg,"  ('",quantCol,"') NOT FOUND !")
        ## look for columns endig with 'intensity'  ?
    quantCol <- quantColInd

    ## extract quantitation
    if(length(quantCol) >0) { abund <- if(length(quantCol) >1)  tmp[,quantCol] else {
      matrix(tmp[,quantCol], ncol=1, dimnames=list(rownames(tmp),NULL))}}   # how to know column-name if single sample ?
    #rownames(abund) <- wrMisc::correctToUnique(pepSeq, silent=silent, callFrom=fxNa)
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
    if(identical("sdrf",sampleNames)) {specPref$sampleNames <- sampleNames; sampleNames <- NULL}
    if(length(sampleNames) ==ncol(abund) && ncol(abund) >0) {
      if(debug) { message(fxNa,"rpd9b") }
      if(length(unique(sampleNames)) < length(sampleNames)) {
        if(!silent) message(fxNa,"Custom sample names not unique, correcting to unique")
        sampleNames <- wrMisc::correctToUnique(sampleNames, callFrom=fxNa) }
      colnames(abund) <- sampleNames
      if(debug) { message(fxNa,"rpd9c") }
    } else {
      colnames(abund) <- sub("Abundance\\.F[[:digit:]]+\\.Sample\\.|Abundances\\.F[[:digit:]]+\\.Sample\\.","Sample.", colnames(abund))
    }

    ## (optional) filter by FDR  (so far use 1st of list where matches are found from argument FDRCol)
    if(length(FDRCol) >0) {
      ## stand : "Found.in.Sample...S32..F32..Sample"   Rfriendly : "Found.in.Sample.in.S33.F33.Sample"
      ## stand : "X..Protein.Groups"                    Rfriendly : "Number.of.Protein.Groups"
      chFDR <- lapply(FDRCol, function(x) {z <- grep(x[1], colnames(tmp)); if(length(z) ==ncol(abund)) z else NULL})
      names(chFDR) <- sapply(FDRCol, function(x) x[1])
      chFDR <- chFDR[which(sapply(chFDR, length) >0)]
      if(length(chFDR) >0) {
        i <- 1              # so far just use 1st instance matching
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
    } }
    if(debug) { message(fxNa,"rpd11b .. dim abund ",nrow(abund)," and ",ncol(abund)); rpd11b <- list()}

    ## Now we are ready to add unique rownames
    if(any(chAc, na.rm=TRUE)) {
      if(!silent) message(fxNa,sum(chAc)," (out of ",length(chAc),") cases of duplicated 'Accession' exist, adding extensions for use as rownames")
      rownames(tmp) <- rownames(annot) <- wrMisc::correctToUnique(annot[,"Accession"], sep="_", atEnd=TRUE, callFrom=fxNa)
    } else rownames(abund) <- rownames(annot) <- annot[,"Accession"]

    ## optional/additional counting results (PSM, no of peptides)
    PSMCol <- if(length(PSMCol) ==1) grep(PSMCol, colnames(tmp)) else NULL
    PepCol <- if(length(PepCol) ==1) grep(PepCol, colnames(tmp)) else NULL
    usTy <- c("PSM","NoOfPeptides")[which(c(length(PSMCol), length(PepCol)) ==ncol(abund))]
    if(length(usTy) >0) {
      counts <- array(NA,dim=c(nrow(abund), ncol(abund), min(1, length(usTy))), dimnames=list(rownames(abund), colnames(abund),usTy))
      if("PSM" %in% usTy) counts[,,"PSM"] <- as.matrix(tmp[,PSMCol])
      if("NoOfPeptides" %in% usTy) counts[,,"NoOfPeptides"] <- as.matrix(tmp[,PepCol])
    } else counts <- NULL
    if(debug) {message(fxNa,"rPD12 .. "); rPD12 <- list(annot=annot,tmp=tmp,abund=abund,sampleNames=sampleNames,specPref=specPref,annotCol=annotCol,refLi=refLi,Rfriendly=Rfriendly,contamCol=contamCol,PSMCol=PSMCol,PepCol=PepCol,counts=counts,infoDat=infoDat)}

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
    if(inherits(quant, "try-error")) quant <- NULL 
    if(debug) { message(fxNa,"rPD13 .. dim quant: ", nrow(quant)," li and  ",ncol(quant)," cols; colnames : ",wrMisc::pasteC(colnames(quant))," ");
      rPD13 <- list(annot=annot,tmp=tmp,abund=abund,quant=quant,sampleNames=sampleNames,normalizeMeth=normalizeMeth,specPref=specPref,groupPref=groupPref,annotCol=annotCol,Rfriendly=Rfriendly,contamCol=contamCol,PSMCol=PSMCol,PepCol=PepCol,counts=counts,infoDat=infoDat, refLi=refLi,suplAnnotFile=suplAnnotFile,sdrf=sdrf, gr=gr)}

    ### GROUPING OF REPLICATES AND SAMPLE META-DATA
    if(length(specPref) >0 && "sampleNames" %in% names(specPref)) groupPref$sampleNames <- specPref$sampleNames
    if(length(suplAnnotFile) >0 || length(sdrf) >0) {
      if(length(sampleNames) %in% c(1, ncol(abund))) groupPref$sampleNames <- sampleNames
      if(length(gr) %in% c(1, ncol(abund))) groupPref$gr <- gr
      setupSd <- readSampleMetaData(sdrf=sdrf, suplAnnotFile=suplAnnotFile, quantMeth="PD", path=path, abund=utils::head(quant), chUnit=isTRUE(groupPref$chUnit), groupPref=groupPref, silent=silent, debug=debug, callFrom=fxNa)
    }
    if(debug) {message(fxNa,"rPD13b ..");  rPD13b <- list(annot=annot,tmp=tmp,abund=abund,suplAnnotFile=suplAnnotFile,quant=quant,sampleNames=sampleNames,specPref=specPref,groupPref=groupPref,annotCol=annotCol,Rfriendly=Rfriendly,contamCol=contamCol,PSMCol=PSMCol,PepCol=PepCol,counts=counts,infoDat=infoDat, refLi=refLi,setupSd=setupSd,gr=gr)}

    ## finish groups of replicates & annotation setupSd
    setupSd <- .checkSetupGroups(abund=abund, setupSd=setupSd, gr=gr, sampleNames=sampleNames, quantMeth="PD", silent=silent, debug=debug, callFrom=fxNa)
    if(debug) {message(fxNa,"rPD13c .."); rPD13c <- list(annot=annot,tmp=tmp,abund=abund,quant=quant,gr=gr,sdrf=sdrf,setupSd=setupSd,sampleNames=sampleNames,specPref=specPref,annotCol=annotCol,Rfriendly=Rfriendly,contamCol=contamCol,PSMCol=PSMCol,PepCol=PepCol,counts=counts,infoDat=infoDat, refLi=refLi)}

    ## harmonize sample-names/1
    colNa <- if(length(setupSd$sampleNames)==ncol(abund)) setupSd$sampleNames else setupSd$groups     # initialize

    ## option : choose (re-)naming of levels & columns on most redundance
    if(length(setupSd$sampleNames)==ncol(quant) && length(setupSd$sdrfSampleNames)==ncol(quant) && any(duplicated(setupSd$level))) {
      ## check if setupSd$sdrfSampleNames  or   setupSd$sampleNames contain better info (some info on replicates)
      chRed1 <- sum(duplicated(sub("(_|\\.| )[[:digit:]]+.*","", setupSd$sdrfSampleNames)), na.rm=TRUE)
      chRed2 <- sum(duplicated(sub("(_|\\.| )[[:digit:]]+.*","", setupSd$sampleNames)), na.rm=TRUE)
      if(chRed1 <1) chRed1 <- sum(duplicated(sub("[[:digit:]]+$","", setupSd$sdrfSampleNames)), na.rm=TRUE)   # remove terminal digits
      if(chRed2 <1) chRed2 <- sum(duplicated(sub("[[:digit:]]+$","", setupSd$sampleNames)), na.rm=TRUE)    # remove terminal digits
      if(any(c(chRed1,chRed2) >0)) {
        if(chRed2 < chRed1) {                                          # use info for levels depending on where more replicates found 
          colNa <- colnames(abund) <- setupSd$sampleNames                                  ## take sample names from sdrf via  setupSd$sampleNames
        } else {
          colNa <- colnames(abund) <- setupSd$sdrfSampleNames <- setupSd$sampleNames }        ## take sample names from sdrf via  setupSd$sdrfSampleNames
        if(length(quant) >0) colnames(quant) <- colNa  }
      #setupSd$level  
    }    
    if(debug) {message(fxNa,"rPD13d .."); rPD13d <- list(annot=annot,tmp=tmp,abund=abund,quant=quant,gr=gr,sdrf=sdrf,setupSd=setupSd,sampleNames=sampleNames,specPref=specPref,annotCol=annotCol,Rfriendly=Rfriendly,contamCol=contamCol,PSMCol=PSMCol,PepCol=PepCol,counts=counts,infoDat=infoDat, refLi=refLi, colNa=colNa)}
    ## option : set order of samples as (init) sdrf
    if(length(quant) >0 && "sdrfOrder" %in% names(sdrf) && isTRUE(as.logical(sdrf["sdrfOrder"])) && length(setupSd$iniSdrfOrder)==ncol(abund) && ncol(abund) >1) {  # set order according to sdrf (only if >1 samples)
      nOrd <- (setupSd$iniSdrfOrder)
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
      colNa <- colNa[nOrd]

      if(length(counts) >0 && length(dim(counts))==3) counts <- array(counts[,nOrd,], dim=c(nrow(counts), length(nOrd), dim(counts)[3]), 
        dimnames=list(rownames(counts), colnames(counts)[nOrd], dimnames(counts)[[3]]))
      if(debug) {message(fxNa,"rPD13e .."); rPD13e <- list(path=path,chPa=chPa,tmp=tmp,groupPref=groupPref,quantCol=quantCol,abund=abund,chNum=chNum,
        sdrf=sdrf,quant=quant,annot=annot,setupSd=setupSd)}  #
        
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
      ### begin old
      ## harmonize sample-names/2
      colNa <- colnames(abund)
      chGr <- grepl("^X[[:digit:]]", colNa)                                                # check & remove heading 'X' from initial column-names starting with digits
      if(any(chGr)) colNa[which(chGr)] <- sub("^X","", colNa[which(chGr)])                 #
    }  
    if(length(colNa)==ncol(abund)) { colnames(abund) <- colNa   
      if(length(setupSd$sampleNames)==length(colNa)) setupSd$sampleNames <- colNa
      if(length(quant) >0 && ncol(quant)==length(colNa)) colnames(quant) <- colNa }
    if(length(dim(counts)) >1 && length(counts) >0) colnames(counts) <- colNa

    if(debug) {message(fxNa,"Read sample-meta data, rPD14"); rPD14 <- list(sdrf=sdrf,suplAnnotFile=suplAnnotFile,abund=abund, quant=quant,refLi=refLi,annot=annot,setupSd=setupSd,counts=counts,sampleNames=sampleNames)}

    ## main plotting of distribution of intensities
    custLay <- NULL
    if(is.numeric(plotGraph) && length(plotGraph) >0) {custLay <- as.integer(plotGraph); plotGraph <- TRUE} else {
        if(!isTRUE(plotGraph)) plotGraph <- FALSE}
    if(plotGraph) .plotQuantDistr(abund=abund, quant=quant, custLay=custLay, normalizeMeth=normalizeMeth, softNa="Proteome Discoverer",
      refLi=refLi, refLiIni=refLiIni, notLogAbund=TRUE, tit=titGraph, las=NULL, silent=silent, callFrom=fxNa, debug=debug)

    ## meta-data
    notes <- c(inpFile=paFi, qmethod="ProteomeDiscoverer", qMethVersion=if(length(infoDat) >0) unique(infoDat$Software.Revision) else NA,
    	identType="protein", rawFilePath= if(length(infoDat) >0) infoDat$File.Name[1] else NA, normalizeMeth=normalizeMeth, call=deparse(match.call()),
      created=as.character(Sys.time()), wrProteo.version=utils::packageVersion("wrProteo"), machine=Sys.info()["nodename"])

    ## final output
    if(isTRUE(separateAnnot)) list(raw=abund, quant=quant, annot=annot, counts=counts, sampleSetup=setupSd, quantNotes=parametersD, notes=notes) else data.frame(quant,annot)
  }
}




#' Additional/final Check And Adjustments To Sample-order After readSampleMetaData()
#'
#' This (low-level) function performs an additional/final check & adjustments to sample-names after readSampleMetaData()
#'
#' @param abund (matrix or data.frame) abundance data, only the colnames will be used
#' @param setupSd (list) describing sammple-setup, typically produced by \code{readSampleMetaData()} (from this package)
#' @param gr (factor) optional custom information about replicate-layout, has priority over setupSd
#' @param sampleNames (character) custom sample-names, has priority over abund and setuoSd
#' @param quantMeth (character) 2-letter abbreviation of name of quantitation-software (eg 'MQ')
#' @param silent (logical) suppress messages
#' @param debug (logical) display additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns an enlaged/updated list 'setupSd' (set setupSd$sampleNames,  setupSd$groups)
#' @seealso used in \code{readProtDiscovererFile},  \code{\link{readMaxQuantFile}}, \code{\link{readProlineFile}}, \code{\link{readFragpipeFile}}
#' @examples
#' abun1 <- matrix(1:16, ncol=8, dimnames=list(NULL,paste("samp", LETTERS[8:1], sep="_")))
#' sdrf1 <- data.frame(source.name=paste(rep(LETTERS[1:4],each=2), 1:2, sep="_"), 
#'   assay.name=paste0("run", 1:8), comment.data.file.=paste0("MSrun", 8:1))
#' setU1 <- list(level=gl(4,2), meth="lowest", sampleNames=paste("samp", LETTERS[1:8], sep="_"), 
#'   sdrfDat=sdrf1, annotBySoft=NULL)
#' .checkSetupGroups(abun1, setU1)
#' @export
.checkSetupGroups <- function(abund, setupSd, gr=NULL, sampleNames=NULL, quantMeth=NULL, silent=FALSE, callFrom=NULL, debug=FALSE) {
  ## additional/final check & adjustments to sample-names after readSampleMetaData()
  ## returns enlaged/updated list 'setupSd' (set setupSd$sampleNames,  setupSd$groups)
  ## 
  ## NOTE : 
  ## $level should be ordered integer (based on $groups)
  ## $groups should be group-names with $level as names
  ## assume 'abund' to be valid matrix of data.frame !!
  ## 'setupSd' as list produced by readSampleMetaData()
  ## 'gr' .. (char vector or factor) designating who should be considered as replicate/same level, as optional custom entry (has prior over automatic groups/levels)
  ## 'sampleNames' .. optional custom entry for sample names (has prior over automatic names)
  ## 'quantMeth' .. (char vect, length=1) design method via abbreviation like 'PD' for ProteomeDiscoverer, 'MQ' for MaxQuant, 'PL' for Proline, etc (used for automatic trimming of default column/sample-names)
  fxNa <- wrMisc::.composeCallName(callFrom, newNa=".checkSetupGroups")
  delPat <- "_[[:digit:]]+$|\\ [[:digit:]]+$|\\-[[:digit:]]+$"             # remove enumerators, ie tailing numbers after separator
  rawExt <- "(\\.RAW$)|(\\.Raw$)|(\\.raw$)"
  
  .corPathW <- function(x) gsub("\\\\", "/", x)
  .adjPat <- function(x) { out <- match(x, unique(x)); names(out) <- names(x); out}  # used
    
  ## start main function
  if(!is.list(setupSd)) { if(length(setupSd) >0) warning(fxNa,"BIZZARE format of 'setupSd', it's content will be lost")
    setupSd <- list()}
  ## finish groups of replicates & annotation setupSd
  if(debug) { message(fxNa,"cSG0"); cSG0 <- list(gr=gr,abund=abund,sampleNames=sampleNames, setupSd=setupSd, quantMeth=quantMeth) }          # sampleNames=sampleNames

  ## grX .. (new) pattern
  ## $level should be ordered integer (based on $groups)
  ## $groups should be group-names with $level as names
  ## special case : fractionated samples : need to change assignment of sample-names here ???

  colNa <- grou <- grou2 <- saNa <- grX <- NULL
  iniSaNa <- iniGr <- FALSE
  iniSetupSd <- length(setupSd) >0   # check if init 'setupSd' given
  if(!is.list(setupSd)) { if(iniSetupSd) as.list(setupSd) else setupSd <- list()}
    

  ## OPTIONS : 
  ## 1) use custom -values (if avail)
  ## 2) extract from 'setupSd' (if avail)
  ## 3) from prev mined from setupSd (groups & levels, not for sampleNames)
  ## 4) pick/mine from colnames of 'abund'

  ## sampleNames/colnames  : if valid user-defied sampleNames are given => use (and define groups based on stripping enumerators)
  if(length(gr)==ncol(abund)) setupSd$groups <- gr
  if(length(sampleNames) == ncol(abund)) {
    ##  1) use externally given sampleName ...
    iniSaNa <- TRUE                                                                                   
    setupSd$sampleNames <- sampleNames
    if(length(gr) !=ncol(abund)) { 
      grpNa <- wrMisc::rmEnumeratorName(sampleNames, nameEnum=c("","Number","No","no","number", "#", "Replicate","Rep","Re","R","replicate","rep","re", "Sample","Samp","Sa","S"), sepEnum=c(" ","-","_","/"), incl="rmEnum")  # remove enum
      setupSd$groups <- grpNa }
  } else {
    ## check if colnames(abund) or sdrf more useful
    ## basic strategy : try picking series of names where some group-names appear when removing enumerators; avoid single specimen groups (no duplication of group-names)
    ## check colnames of abund
    colNaA <- wrMisc::rmSharedWords(colnames(abund), sep=c("_", " ", ".", "-", ";"), callFrom=fxNa, silent=!debug)  # simplify
    ## check for enumerators    
    colNaA2 <- wrMisc::rmEnumeratorName(colNaA, nameEnum=c("","Number","No","no","number", "#", "Replicate","Rep","Re","R","replicate","rep","re", "Sample","Samp","Sa","S"), sepEnum=c(" ","-","_","/"), incl="rmEnum")  # remove enum
    chDuA2 <- duplicated(colNaA2)

    ## initialize
    if(length(setupSd$sampleNames) !=ncol(abund)) setupSd$sampleNames <- colNaA else {
      ## may choose from prev sampleNames or based on colnames
      if(sum(chDuA2, na.rm=TRUE) >= sum(duplicated(setupSd$level), na.rm=TRUE)) {setupSd$sampleNames <- colNaA} # use colnames(abund) only if at least same number of duplicates (for groups) 
    }
    if(length(setupSd$groups) !=ncol(abund)) setupSd$groups <- colNaA2 else {
      if(sum(duplicated(colNaA2), na.rm=TRUE) >= sum(duplicated(setupSd$groups), na.rm=TRUE)) {setupSd$groups <- colNaA2} # use colnames(abund) only if at least same number of duplicates (for groups) 
    }    
    if(debug) { message(fxNa,"cSG1a"); cSG1a <- list(gr=gr,abund=abund,sampleNames=sampleNames, setupSd=setupSd, quantMeth=quantMeth) }          # sampleNames=sampleNames
    
    ## check sdrf
    colNaS <- repNaS <- NULL
    if(length(setupSd$sdrfDat) >0) {   
      ## rather try to find sampleNames in sdrf
      chSdrfCol <- c("source.name", "characteristics.organism.part.","characteristics.cell.type.","characteristics.disease.","assay.name","comment.label.","comment.file.uri.","comment.data.file.","characteristics.spiked.compound." )  
      sdrfT <- setupSd$sdrfDat[,wrMisc::naOmit(match(chSdrfCol, colnames(setupSd$sdrfDat))), drop=FALSE]
      ## check for all different names (potential sample-names)
      chSd1 <- colSums(apply(sdrfT, 2, duplicated))
      if(any(chSd1 ==0)) { 
        sdrfT <- sdrfT[,which(chSd1==0), drop=FALSE]}    # retain only cols with some duplicated names (wo enumerators)
      if(length(sdrfT) >0) {  
        ## now check if containing enumerators
        sdrfT <- apply(sdrfT, 2, wrMisc::rmSharedWords, sep=c("_", " ", ".", "-", ";"), callFrom=fxNa, silent=!debug)
        sdrfT2 <- apply(sdrfT, 2, wrMisc::rmEnumeratorName, nameEnum=c("","Number","No","N","no","number", "#", "Replicate","Rep","Re","R","replicate","rep","re", "Sample","Samp","Sa","S"), sepEnum=c(" ","-","_","/"), incl="rmEnum")  # remove enum
        ## choose likely columns
        if(length(dim(sdrfT2)) >1 && ncol(sdrfT2) >1) {
          chDu <- colSums(apply(sdrfT2, 2, duplicated))
          ## check for containing text
          chTe <- colSums(apply(sdrfT2, 2, function(x) grepl("[[:alpha:]]", x))) ==nrow(setupSd$sdrfDat)
          if(any(!chTe)) chDu[which(!chTe)] <- 0
          if(any(chDu >0)) { 
            useColS <- chDu[which(chDu >0 & chDu < nrow(sdrfT) -1)]
            useColS <- which(chDu==  if(length(useColS) >1) sort(useColS)[ceiling(length(useColS)/2)] else useColS)[1]   # so far : choose median number of replicates
            colNaS <- sdrfT[,useColS]
            repNaS <- sdrfT2[,useColS]
            ## integrate sampleNames & replicateNames from sdrf to setupSd
            setupSd$sdrfSampleNames <- colNaS
            setupSd$sdrfReplicateNames <- repNaS
          }       
        } else { colNaS <- sdrfT; repNaS <- sdrfT2 }         
      }                                                   
      ## now compare to colnames of abund & choose between colnames(abund) or sdrf : choose from sdrf if at least same number of grouops/levels
      if(length(repNaS) >0 && sum(duplicated(repNaS)) < length(repNaS) -1 && length(unique(repNaS)) >= length(unique(colNaA2)) ) { 
        if(!silent) message(fxNa,"Setting setupSd$sampleNames based on sdrf    cSG1a2")
        setupSd$sampleNames <- colNaS 
        if(length(gr) !=ncol(abund)) { setupSd$groups <- repNaS }         
      }
    }
    if(debug) { message(fxNa,"cSG1b"); cSG1b <- list(gr=gr,abund=abund,sampleNames=sampleNames, setupSd=setupSd, quantMeth=quantMeth) }          # sampleNames=sampleNames

    ## PD : mine suplData
    if("PD" %in% quantMeth && length(setupSd$annotBySoft$File.Name) ==ncol(abund)) {
      chOr <- try(match(basename(setupSd$annotBySoft$File.Name), basename(if("sdrfDat" %in% names(setupSd)) setupSd$sdrfDat$comment.data.file. else setupSd$annotBySoft$filePath)), silent=TRUE)
      if(inherits(chOr, "try-error") || !any(is.na(chOr))) {
        if(all(chOr ==1:ncol(abund), na.rm=TRUE)) {
          colNaP <-  sub("\\.raw$|\\.RAW$","", setupSd$annotBySoft$File.Name)
          colNaP <- wrMisc::rmSharedWords(colNaP, sep=c("_", " ", ".", "-", ";"), callFrom=fxNa, silent=!debug)  # simplify
          repNaP <- wrMisc::rmEnumeratorName(colNaP, nameEnum=c("","Number","No","no","number", "#", "Replicate","Rep","Re","R","replicate","rep","re", "Sample","Samp","Sa","S"), sepEnum=c(" ","-","_","/"), incl="rmEnum")  # remove enum
          ## now compare to colnames of abund & choose between colnames(abund) or sdrf : choose where more levels/groups
          chDu <- sum(duplicated(repNaP))
          if(chDu < length(colNaP) -1 && chDu > sum(duplicated(setupSd$groups))) { 
            setupSd$sampleNames <- colNaP 
            if(length(gr) !=ncol(abund)) { setupSd$groups <- repNaP }
          }
        }
      }
    }   ## finish PD specific check

    ## potentially need to adjust units for groups (only groups !)
    chUnit <- wrMisc::checkUnitPrefix(setupSd$groups, stringentSearch=TRUE, callFrom=fxNa, silent=!debug)
    if(length(chUnit) ==1) {
      setupSd$groups <- wrMisc::adjustUnitPrefix(setupSd$groups, unit=chUnit)
    }
     
  }     
  setupSd$level <- match(setupSd$groups, unique(setupSd$groups))
  names(setupSd$level) <- setupSd$groups
  names(setupSd$groups) <- setupSd$level       
  ## done setting sample-names and groups of replicates    
    

  if(debug) {message(fxNa,"cSG2"); cSG2 <- list(abund=abund,setupSd=setupSd,gr=gr,quantMeth=quantMeth,saNa=saNa)}

  #if(debug) {message(fxNa,"cSG2"); cSG2 <- list(sampleNames=sampleNames,gr=gr,grX=grX,abund=abund,iniSaNa=iniSaNa,saNa=saNa,iniGr=iniGr,setupSd=setupSd)}

  #if(debug) {message(fxNa,"cSG3"); cSG3 <- list(sampleNames=sampleNames,gr=gr,grX=grX,abund=abund,iniSaNa=iniSaNa,saNa=saNa,iniGr=iniGr,setupSd=setupSd)}

  ## simplify terms ?
  #if(TRUE) {
  #  setupSd$sampleNames <- wrMisc::rmSharedWords(setupSd$sampleNames, sep=c("_"," ","-","/"))
  #  setupSd$groups <-  names(setupSd$level) <- wrMisc::rmSharedWords(setupSd$groups, sep=c("_"," ","-","/"))
  #}
  if(debug) {message(fxNa,"cSG4"); cSG4 <- list(abund=abund,setupSd=setupSd,gr=gr,quantMeth=quantMeth,saNa=saNa)}
  ## return results-object (finish .checkSetupGroups() )
  setupSd
  }



#' Generic Plotting Of Density Distribution For Quantitation Import-functions
#'
#' This (low-level) function allows (generic) plotting of density distribution for quantitation import-functions
#'
#' @param abund (matrix or data.frame) abundance data, will be plottes as distribution
#' @param quant (matrix or data.frame) optional additional abundance data, to plot 2nd distribution, eg of normalized data
#' @param custLay (matrix) describing sammple-setup, typically produced by
#' @param normalizeMeth (character, length=1) name of normalization method (will be displayed in title of figure)
#' @param softNa (character, length=1) name of quantitation-software (typically 2-letter abbreviation, eg 'MQ')
#' @param refLi (integer) to display number reference lines
#' @param refLiIni (integer) to display initial number reference lines
#' @param notLogAbund (logical) set to \code{TRUE} if \code{abund} is linear but should be plotted as log2 
#' @param figMarg (numeric, length=4) custom figure margins (will be passed to \code{\link[graphics]{par}}), defaults to c(3.5, 3.5, 3, 1)
#' @param tit (character) custom title
#' @param las (integer) indicate orientation of text in axes
#' @param cexAxis (numeric) size of numeric axis labels as cex-expansion factor (see also \code{\link[graphics]{par}})
#' @param nameSer (character) custom label for data-sets or columns (length must match number of data-sets)
#' @param cexNameSer (numeric) size of individual data-series labels as cex-expansion factor (see also \code{\link[graphics]{par}})
#' @param silent (logical) suppress messages
#' @param debug (logical) display additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns logical value (if data were valid for plotting) and produces a density dustribution figure (if data were found valid)
#' @seealso used in \code{readProtDiscovererFile},  \code{\link{readMaxQuantFile}}, \code{\link{readProlineFile}}, \code{\link{readFragpipeFile}}
#' @examples
#' set.seed(2018);  datT8 <- matrix(round(rnorm(800) +3,1), nc=8, dimnames=list(paste(
#'   "li",1:100,sep=""), paste(rep(LETTERS[1:3],c(3,3,2)),letters[18:25],sep="")))
#' .plotQuantDistr(datT8, quant=NULL, refLi=NULL, tit="Synthetic Data Distribution")                                
#' @export
.plotQuantDistr <- function(abund, quant, custLay=NULL, normalizeMeth=NULL, softNa=NULL, refLi=NULL, refLiIni=NULL, notLogAbund=NA, figMarg=c(3.5, 3.5, 3, 1), tit=NULL, las=NULL, cexAxis=0.8, nameSer=NULL, cexNameSer=NULL, silent=FALSE, callFrom=NULL, debug=FALSE) {
  ## generic plotting of densirt distribution for quantitation import-functions
  ## assume 'abund' (raw, non-normalized) and 'quant' (final normalized) to be valid matrix of data.frame !!
  ## 'custLay' ..(matrix) for layout()
  ## 'normalizeMeth' (character)
  ## 'softNa' .. (char vect, length=1) design method, used in display only
  ## 'refLi' .. (integer) reference line
  ## 'refLiIni' .. (integer) initial reference line(s)

  fxNa <- wrMisc::.composeCallName(callFrom, newNa=".plotQuantDistr")
  oparMar <- graphics::par("mar")                     # old margins, for rest after figure
  oparLayout <- graphics::par("mfcol")                # old layout, for rest after figure
  on.exit(graphics::par(mar=oparMar, mfcol=oparLayout))            # restore old mar settings
  
  if(debug) {message(fxNa,"pQD0 .. length custLay ", length(custLay)); pQD0 <- list(abund=abund,quant=quant,custLay=custLay,normalizeMeth=normalizeMeth,softNa=softNa,refLi=refLi)}
  plotGraph <- TRUE                                         # useful ? (to report if plot can be drawn as output ?)
  ## check if abund & quant are redundant (while normalizeMeth="none")
  abundQuantRed <- FALSE       # see if abund and quant are redundant (ie no need to plot twice)
  if("none" %in% normalizeMeth && length(abund) >0 && length(quant) >0) {
    if(identical(quant, abund) || identical(quant, log2(abund))) { abundQuantRed <- TRUE  
      if(debug) message(fxNa,"No need for 2nd plot, 'abund' and 'quant' are identical ..") }  
  }
      
  if(length(custLay) >0) graphics::layout(custLay) else {
    if(!identical(normalizeMeth,"none") && (length(abund) >0 && length(quant) >0) && !abundQuantRed) { ch1 <- try(graphics::layout(1:2), silent=TRUE)
      if(inherits(ch1, "try-error")) message(fxNa,"Problem with figure, Need to restore layout ..") else if(debug) message(fxNa,"Setting layout to 1:2") }}
  if(debug) { message(fxNa,"pQD1"); pQD1 <- list(abund=abund,quant=quant,custLay=custLay,normalizeMeth=normalizeMeth,softNa=softNa,refLi=refLi)}

  if(length(las) != 1 || !is.integer(las)) {
     ch1 <- c(if(length(abund) >0) ncol(abund) >7 || stats::median(nchar(colnames(abund)), na.rm=TRUE) >8 else FALSE,
       if(length(quant) >0) ncol(quant) >7 || stats::median(nchar(colnames(quant)), na.rm=TRUE) >8 else FALSE)
    las <- if(any(ch1)) 2 else 1 }
  if(length(figMarg) != 4 || !is.numeric(figMarg)) { figMarg <- c(3.5, 3.8, 3, 1)     # mar: bot,le,top,ri
    if(debug) message(fxNa,"Invalid entry for argument 'figArg' (must be umeric of length=4), setting to default")}
  ch1 <- try(graphics::par(mar=figMarg), silent=TRUE)                     
  if(inherits(ch1, "try-error")) message(fxNa,"Problem with figure, Need to restore mar ..")
  if(is.null(tit)) tit <- paste(softNa," Quantification ")
  titSu <- if(length(refLi) >0) paste0(c(if(length(refLiIni) ==1) c("'",refLiIni,"'") else c(" by ",length(refLi)," selected lines")),collapse="")  else NULL #
  usePa <- c("wrGraph","sm")
  misPa <- !sapply(usePa, function(pkg) requireNamespace(pkg, quietly=TRUE))
  titQ <- if(length(abund) >0) paste(tit, "(initial)",sep=" ") else tit
  .findSuplLi <- function(x, n=3) {
    if(length(x) >2) pretty(stats::quantile(x, c(0.15, 0.85), na.rm=TRUE), n) else NULL  
  }
  ## check log-status, ie notLogAbund 
  if(length(abund) >0 && length(quant) >0) {
   if(length(notLogAbund) <1) { notLogAbund <- NA
     if(!silent) message(fxNa,"Invalid entry for 'notLogAbund', setting to default =NA") }
   if(is.na(notLogAbund)) {
    dAb <- diff(range(abund, na.rm=TRUE))
    dQa <- diff(range(quant, na.rm=TRUE))
    sugLogAbun <- dAb > 2^(signif(dQa,3) -2) && dAb > 3*dQa
    if(debug) message(fxNa,"Trying to figure out if log2 should be taken for 'abund': ",dAb > 2^(signif(dQa,3) -2)," and ", dAb > 3*dQa)
    if(sugLogAbun) {
      notLogAbund <- TRUE }   # abund is very likely linear scale, while quant is likely log-scale
  } }
  msg <- c("Invalid entry for 'notLogAbund',","Unknown log-status,"," assuming data is log2, ie notLogAbund=FALSE")
  if(length(notLogAbund) <1 && !silent) { notLogAbund <- FALSE
    if(!silent) message(fxNa, msg[-2])
  } else if(any(is.na(notLogAbund))) { notLogAbund <- notLogAbund[which(is.na(notLogAbund))] <- FALSE
    if(!silent) message(fxNa, msg[-1]) }      
  
  colNa <- if(length(abund) >0) colnames(abund) else colnames(quant)
  if(length(colNa) >0) {
    ncharAx <- stats::quantile(nchar(colNa), c(0.5,0.9), na.rm=TRUE)
    if(length(cexNameSer) !=1 || any(is.na(cexNameSer))) cexNameSer <- sort(round( c(7/ncharAx[2], 1.9/log(length(colNa)), 1, 0.48), 2))[2]     # adjust cex series-names to length and ncol
    suplLineA <- if(length(abund) >0) {.findSuplLi(if(isTRUE(notLogAbund)) log2(abund) else abund)} else NA
    suplLineQ <- if(length(quant) >0)  .findSuplLi(quant) else NA

  } else { ncharAx <- cexNameSer <- suplLineA <- suplLineQ <- NA
    if(debug) message(fxNa," Note : BOTH 'abund' and 'quant' are EMPTY - NOTHING TO PLOT")}  
  
  ## check axis cex
  if(debug) { message(fxNa,"pQD2 .. misPa ", wrMisc::pasteC(misPa,quoteC="'"),";  abund is linear ",notLogAbund); pQD2 <- list(abund=abund,quant=quant,custLay=custLay,normalizeMeth=normalizeMeth,softNa=softNa,refLi=refLi,tit=tit,titQ=titQ,suplLineA=suplLineA,suplLineQ=suplLineQ,ncharAx=ncharAx) }
  
  chNeg <- sum(abund <0, na.rm=TRUE)
  if(chNeg >0 && notLogAbund) { notLogAbund <- FALSE
    if(!silent) message(fxNa,"Data suggest taking log2 might be useful, but due to presence of ",chNeg," negative values this is not possible")  
  } 
  if(any(misPa, na.rm=TRUE)) {
    if(!silent) message(fxNa,"Please install package(s) ", wrMisc::pasteC(usePa[which(misPa)])," from CRAN for drawing vioplots")
    if(length(abund) >0) {
      ## wrGraph not available : simple boxplot
      ch1 <- try(graphics::boxplot(if(isTRUE(notLogAbund)) suppressWarnings(log2(abund)) else abund, main=titQ, las=1, outline=FALSE), silent=TRUE)
      if(inherits(ch1, "try-error")) {plotGraph <- FALSE; if(!silent) message(fxNa,"UNABLE to draw boxplot of distribution !!  ", sub("^Error in",":",ch1))} else {
        graphics::abline(h=suplLineA, lty=2, col=grDevices::grey(0.6))} }

    ## plot normalized
    if(length(quant) >0 && !abundQuantRed) {      
      if(identical(normalizeMeth,"none") || length(quant) <0) {
        if(debug) {message(fxNa,"pQD3 .. dim quant: ", nrow(quant)," li and  ",ncol(quant)," cols; colnames : ",wrMisc::pasteC(colnames(quant))," ")}
        ch1 <- try(graphics::boxplot(quant, main=paste(tit," (",normalizeMeth,"-normalized",titSu,")"), las=1, outline=FALSE), silent=TRUE)
        if(inherits(ch1, "try-error")) if(!silent) message(fxNa,"UNABLE to draw boxplot of normalized distribution !!  ", sub("^Error in",":",ch1)) else {
          graphics::abline(h=suplLineQ, lty=2, col=grDevices::grey(0.6)) } } }

  } else {                                            # wrGraph and sm are available
    if(debug) { message(fxNa,"pQD4  draw vioplotW ,  abund is linear ",notLogAbund); pQD4 <- list(abund=abund,quant=quant,tit=tit,normalizeMeth=normalizeMeth,softNa=softNa,refLi=refLi,titQ=titQ,notLogAbund=notLogAbund,titSu=titSu,las=las, cexAxis=cexAxis, nameSer=nameSer, cexNameSer=cexNameSer,notLogAbund=notLogAbund,abundQuantRed=abundQuantRed ) }
    if(length(abund) >0) {       
      if(length(chNeg) <1) abund <- suppressWarnings(log2(abund))           #  presume as true 'raw', ie  NOT log2
      if(debug) message(fxNa," Try 1st Vioplot ,  is linear abundance data ",notLogAbund,",  cexNameSer ",cexNameSer)
      ch1 <- try(wrGraph::vioplotW(if(isTRUE(notLogAbund)) log2(abund) else abund, tit=titQ, wex=NULL, las=las, cexAxis=cexAxis, nameSer=nameSer, cexNameSer=cexNameSer, horizontal=FALSE, silent=debug, debug=debug, callFrom=fxNa), silent=TRUE)
      if(inherits(ch1, "try-error")) {plotGraph <- FALSE; if(!silent) message(fxNa,"UNABLE to plot vioplotW !!  ", sub("^Error in",":",ch1))
      } else graphics::abline(h=suplLineA, lty=2, col=grDevices::grey(0.6))
    }
    ## now normalized (and/or log-scale)
    if(length(quant) >0  && !abundQuantRed) {
      if(debug) {message(fxNa,"pQD5  draw norm/quant as vioplotW() ", length(quant) >0)}
      tit1 <- if(length(normalizeMeth) ==1) paste(tit,", ",normalizeMeth,"-normalized",titSu) else paste(tit,", ",titSu)
      if(debug) message(fxNa," Try 2nd Vioplot ")
      ch1 <- try(wrGraph::vioplotW(quant, tit=tit1, wex=NULL, las=las, cexAxis=cexAxis, nameSer=nameSer, cexNameSer=cexNameSer, horizontal=FALSE, silent=debug, debug=debug,callFrom=fxNa), silent=TRUE)
      if(inherits(ch1, "try-error")) { if(!silent) message(fxNa,"UNABLE to plot vioplotW for normalized data !!  ", sub("^Error in",":",ch1))
      } else graphics::abline(h=suplLineQ, lty=2, col=grDevices::grey(0.6))
    }
  }
  plotGraph }
   


#' Extract Additional Information To Construct The Colum 'SpecType', Allows Adding Information From Fasta
#'
#' This (low-level) function creates the column annot[,'SpecType'] which may help distinguishing different lines/proteins.
#' This information may, for example, be used to normalize only to all proteins of a common backgroud matrix (species).
#' In order to compare \code{specPref} a species-column will be added to the annotation (\code{annot}) - if not already present
#' If $mainSpecies or $conta: match to annot[,"Species"], annot[,"EntryName"], annot[,"GeneName"], if length==1 grep in  annot[,"Species"]
#'
#' @details
#' Different to \code{readSampleMetaData} this function also considers the main annotation as axtracted with main quantification data.
#' For example, this function can complement protein annotation data if columns 'Accession','EntryName' or 'SpecType' are missing
#' 
#' 
#' @param specPref (list) may contain $mainSpecies, $conta ...
#' @param annot (matrix) main protein annotation
#' @param useColumn (factor) columns from annot to use/mine
#' @param suplInp (matrix) additional custom annotation
#' @param soft (character, length=1) additional info which software was initially used (so far only special treatmentr for IB)
#' @param silent (logical) suppress messages
#' @param debug (logical) display additional messages for debugging (starting with 'mainSpecies','conta' and others - later may overwrite prev settings)
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a matrix with additional column 'SpecType'
#' @seealso used in \code{readProtDiscovererFile},  \code{\link{readMaxQuantFile}}, \code{\link{readProlineFile}}, \code{\link{readFragpipeFile}}
#' @examples
#' annot1 <- cbind( Leading.razor.protein=c("sp|P00925|ENO2_YEAST",
#'   "sp|Q3E792|RS25A_YEAST", "sp|P09938|RIR2_YEAST", "sp|P09938|RIR2_YEAST",
#'   "sp|Q99186|AP2M_YEAST", "sp|P00915|CAH1_HUMAN"), 
#'   Species= rep(c("Saccharomyces cerevisiae","Homo sapiens"), c(5,1)))
#' specPref1 <- list(conta="CON_|LYSC_CHICK", 
#'   mainSpecies="OS=Saccharomyces cerevisiae", spike="P00915")   # MQ type
#' .extrSpecPref(specPref1, annot1, useColumn=c("Species","Leading.razor.protein"))  
#' @export
.extrSpecPref <- function(specPref, annot, useColumn=c("Species","EntryName","GeneName","Accession"), suplInp=NULL, soft=NA, silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## create column annot[,'SpecType']
  ## if $mainSpecies or $conta: match to annot[,"Species"], annot[,"EntryName"], annot[,"GeneName"], if length==1 grep in  annot[,"Species"]
  ## if other : match to annot[,"Species"], annot[,"Accession"], annot[,"EntryName"], annot[,"GeneName"], if length==1 grep in   annot[,"EntryName"], annot[,"GeneName"], annot[,"Species"]
  ## 'suplInp' add'l matrix of annot (really needed ?)
  ## return results in column annot[,"SpecType"] (starting with 'mainSpecies','conta' and others - later may overwrite prev settings)
  ## special for PD : optional useColumn[5:6] : look by grep for specPref tags in cols "Majority.protein.IDs" and "Fasta.headers"
  ## 'specPref' ..(list) may contain $mainSpecies, $conta ...
  ## 'annot' ..(matrix) main protein annotation
  ## 'useColumn' ..(character) columns from annot to use/mine
  ## 'suplInp' ..(matrix) additional custom annotation
  ## aim : add new column 'SpecType' 
  fxNa <- wrMisc::.composeCallName(callFrom, newNa=".extrSpecPref")
  if(debug) {message(fxNa," eSP0"); eSP0 <- list(specPref=specPref,annot=annot,useColumn=useColumn,suplInp=suplInp)}

  if(length(annot) <1 || length(dim(annot)) !=2) stop("invalid 'annot' (must be matrix or data.frame)")
  ## check suplInp & match to useColumn, add to useColumn
  if(length(useColumn) >4 && length(suplInp) >0 && length(dim(suplInp))==2) {
    chAnn <- useColumn[5:length(useColumn)] %in% colnames(suplInp)
    if(any(!chAnn)) useColumn <- c(useColumn[1:4], useColumn[(5:length(useColumn))[which(chAnn)]])
    if(length(useColumn) <5) suplInp <- NULL
  } else suplInp <- NULL

  ## check useColumn
  chAnn <- useColumn[1:min(length(useColumn), 4)] %in% colnames(annot)   # check for length useCol=0  ?? # nolint
  if(all(!chAnn)) stop("Unknown/Non-standard 'annot' (missing colnames ",wrMisc::pasteC(useColumn[which(!chAnn)], quoteC="'"),")")
  if(debug) {message(fxNa,"eSP1"); eSP1 <- list(specPref=specPref,annot=annot,useColumn=useColumn,suplInp=suplInp,chAnn=chAnn, soft=soft)}

  if(!"SpecType" %in% colnames(annot)) {annot <- cbind(annot, SpecType=rep(NA, nrow(annot))); if(debug) message(fxNa,"Adding column 'SpecType' to 'annot'")}
  if(length(specPref) > 0) specPref <- specPref[which(sapply(specPref, length) >0)]     # remove empty .
  
  ## split Leading.razor.protein  if "Accession","EntryName" NOT in annot
  if(!all(c("Accession","EntryName") %in% colnames(annot))) {
    chSplCol <- c("Leading.razor.protein","Proteins","Protein.Groups")    # columns for choice to use ..
    chSplCol <- which(colnames(annot) %in% chSplCol)
    if(length(chSplCol) >0) {
      chSplCol <- chSplCol[1]
      if(isTRUE("IB" %in% soft) && "EntryName" %in% colnames(annot)) {
        ## IB : rather use last protein
        nc1 <- nchar(annot[,"EntryName"])
        tmpA <- sub(".+_{2,4}","", sub(".+\\|","", annot[,chSplCol]))    #  last EntryName from eg  "sp|P0C0T4|RS25B_YEAST__sp|Q3E792|RS25A_YEAST", "sp|P09938|RIR2_YEAST"
        Accession <- sub(".+_{1,3}","", gsub("[[:alpha:]]+\\|","", sub("\\|CON_$","", substr(annot[,chSplCol], 1,nchar(annot[,chSplCol]) - nc1 -1))  ))
        newEntryNa <- annot[,"EntryName"]
      } else {
        ## note : init column may contain multiple annot
        tmpA <- sub("^[[:lower:]]+\\|","", annot[,chSplCol])                              # remove (init) DB ID
        paMult <- "_{2,}[[:lower:]]{2,}\\|[A-Z,0-9]+\\|[[:upper:]].+"
        chMult <- grep(paMult, tmpA)
        if(length(chMult) >0) tmpA[chMult] <- sub(paMult,"", tmpA[chMult])   # use 1st of combined multiple entries: contains Accession + EntryName
        Accession <- sub("(\\|.+)|(;.+)","", tmpA) 
        newEntryNa <- sub("(\\|.+)|(;.+)","", substring(tmpA, nchar(sub("\\|.+","", Accession)) +2))
      }
      chCol <- match("entryName", colnames(annot))
      if(!is.na(chCol)) colnames(annot)[chCol] <- "EntryName" 

      chCol <- colnames(annot) %in% c("EntryName")
      if(any(chCol)) { useLi <- which(is.na(annot[,which(chCol)]) | nchar(annot[,which(chCol)]) <1)
        if(length(useLi) >0) annot[useLi,which(chCol)] <- newEntryNa[useLi]     # add new to incomplete EntryName
        newEntryNa <- annot[,which(chCol)]                                      # update
      } else annot <- cbind(annot, EntryName =newEntryNa)

      #why#chCol <- which(colnames(annot) %in% c("EntryName","entryName"))
      #why#annot <- annot[,-chCol, drop=FALSE]      
      annotTail <- if(ncol(annot) <3) NULL else annot[,3:ncol(annot)]
      if(length(annotTail) >0) {
        chCol <- which(colnames(annotTail) %in% c("Accession","EntryName"))
        if(length(chCol) > 0) annotTail <- annotTail[,-chCol, drop=FALSE]   # avoid duplicating columns (with old content)
      }      
      annot <- cbind(annot[,1:min(2, ncol(annot))], Accession=Accession, EntryName=annot[,"EntryName"])
      if(length(annotTail) >0) annot <- cbind(annot, annotTail)
      
      rm(Accession, tmpA, newEntryNa) }
  }   
  if(!"Species" %in% colnames(annot)) annot <- cbind(annot, Species=NA)
  
  if(debug) {message(fxNa,"eSP1a"); eSP1a <- list(specPref=specPref,annot=annot,useColumn=useColumn,suplInp=suplInp,chAnn=chAnn)}

  ## inspect specPref
  matrixNa <- c("matrix","Matrix", "mainSpecies","MainSpecies")                           # equivalent names for treating as $matrix
  if(is.character(specPref)) specPref <- as.list(specPref)
  if(any(matrixNa %in% names(specPref))) {
    chOthNames <- matrixNa[-1] %in% names(specPref)
    if(any(chOthNames))  specPref$matrix <- specPref[which(chOthNames)]     # copy content of similar names into specPref$matrix  
    if(length(specPref$matrix) ==1) {           # if length=1 this must be a name for a species or group/collection of proteins
      specPref$IniMatrix <- specPref$matrix
      ## check for known collections (so far only UPSx) to extract specific IDs
      if(grepl("^UPS", specPref$matrix[1])) {         # realistic that UPS1 used as matrix ??
        specPref$SpeciesMatrix <- "Homo sapiens"   
        specPref$EntryNameMatrix <- getUPS1acc()$EntryName
        specPref$matrix <- getUPS1acc()$acOld
        names(specPref$matrix) <- names(specPref$EntryNumberMatrix) <- getUPS1acc()$EntryName      # eg 'CAH1'
      }                                                                     # add similar to UPS1 if available ...
    } else if(length(specPref$matrix) <6) specPref$matrix <- sapply(specPref$matrix, inspectSpeciesIndic, silent=TRUE, debug=debug)
  }
  if("spike" %in% names(specPref)) {
    if(length(specPref$spike) ==1) {
      specPref$IniSpike <- specPref$spike
      if(grepl("^UPS", specPref$spike[1])) { 
        if(debug) {message(fxNa,"Found instance of UPS as spike ..")}
        specPref$SpeciesSpike <- "Homo sapiens"       # add more if available ...
        specPref$EntryNameSpike <- getUPS1acc()$EntryName   # eg 'CAH1_HUMAN'
        specPref$EntryNumberSpike <- getUPS1acc()$acNew     # eg 'CAH1_HUMAN'
        specPref$spike <- getUPS1acc()$acOld
        names(specPref$spike) <- names(specPref$EntryNumberSpike) <- getUPS1acc()$EntryName      # eg 'CAH1'
        ## note: single Gene name  may correspond to multiple UniProt EntryName or accession number
      }
    } else if(length(specPref$spike) <6) specPref$spike <- sapply(specPref$spike, inspectSpeciesIndic, silent=TRUE, debug=debug)
  }
  if(debug) {message(fxNa,"eSP1c"); eSP1c <- list(specPref=specPref,annot=annot,useColumn=useColumn,suplInp=suplInp,chAnn=chAnn)}


  ## add fasta & integrate to annot as Species
  chFasta <- length(specPref) > 0 && "fasta" %in% names(specPref)
  fasta <- NULL     # initialize (for debuging)
  eqiCol <- matrix(c("Accession","EntryName","Description","Species","Modif","GN",  "uniqueIdentifier","entryName","proteinName","OS","Modif","GN"), ncol=2, dimnames=list(NULL, c("annot","fasta")))

  if(chFasta) {
    if(is.list(specPref) && length(dim(specPref$fasta)==2)) fasta <- specPref$fasta
    if(file.exists(specPref$fasta)) fasta <- readFasta2(specPref$fasta, delim="|", databaseSign=c("sp","tr","generic","gi","syn"), removeEntries=NULL, tableOut=TRUE, 
      UniprSep=c("OS=","OX=","GN=","PE=","SV="), callFrom=fxNa,silent=silent,debug=debug)  
    if(length(fasta) >0 && !inherits(fasta,"try-error")) {
      if(debug) {message(fxNa,"eSP2"); eSP2 <- list(specPref=specPref,annot=annot,useColumn=useColumn,suplInp=suplInp,chAnn=chAnn,fasta=fasta, debug=debug)}
        #fasta0= fasta
       ## "O76070" is in annot$Accession, thus it should be ok to mark as spike
       ##
       ## (option) NEED TO ADD spike-proteins like UPS1 to fasta ("O76070" not part of fasta !!)
       ## check if specPref$spike contains Accessions not in fasta AND neither in ; how to know we're in UPS1 ?
      if(length(fasta) >1 && all(c("spike", "IniSpike") %in% names(specPref))) {
        #chSp <- markLi[which(!annot[markLi,"Accession"] %in% getUPS1acc()$acFull)] 


        faNew <- matrix(NA, nrow=length(specPref$spike), ncol=ncol(fasta), dimnames=list(NULL,colnames(fasta)))
        faNew[,"database"] <- rep("sp", nrow(faNew))             
        faNew[,"OS"] <- rep(specPref$SpeciesSpike, nrow(faNew))             

        if(all(c("EntryNameSpike","EntryNumberSpike") %in% names(specPref))) {
          faNew <- cbind(database=rep("sp",length(specPref$spike)), uniqueIdentifier=as.character(specPref$EntryNumberSpike), entryName=specPref$EntryNameSpike,
            proteinName=rep(NA,length(specPref$spike)), sequence=rep(NA,length(specPref$spike)), OS=rep(specPref$SpeciesSpike, length(specPref$spike)), OX=NA,GN=NA,PE=NA,SV=NA)
        } else {
          if(length(names(specPref$spike)) >0) {
            faNew["uniqueIdentifier"] <- as.character(specPref$spike)
            faNew["entryName"] <- names(specPref$spike)
          } else {
            if(debug) message(fxNa,"Note : Very poorly documented spike-in !  (Trouble ahead ?)")
            faNew["uniqueIdentifier"] <- specPref$spike
          }
        }
        ## check if new lines for spike not already present in fasta
        chLi <- !faNew[,"uniqueIdentifier"] %in% fasta[,"uniqueIdentifier"]      ## need to correct faNew[,"uniqueIdentifier"] to old writing ?      
        if(any(chLi)) fasta <- rbind(fasta, faNew[which(chLi),])
        if(debug) {message(fxNa,length(specPref$spike)," Additional ",specPref$IniSpike," proteins were added to Fasta information","  eSP2aa");
          eSP2aa <- list(specPref=specPref,annot=annot,useColumn=useColumn,suplInp=suplInp,chAnn=chAnn,fasta=fasta, debug=debug)}

        ## complete annot / EntryName (if absent, like with PD)
        chLi <- is.na(annot[,"EntryName"]) | nchar(annot[,"EntryName"]) <1
        if(any(chLi)) {
          chID <- annot[which(chLi),"Accession"]
          chLi2 <- match(annot[which(chLi),"Accession"], fasta[,"uniqueIdentifier"])
          if(length(wrMisc::naOmit(chLi2)) >0) annot[which(chLi),"EntryName"] <- fasta[chLi2,"entryName"] 
        }

        ## complete annot / Description/proteinName (if absent...)
        if(!"Description" %in% colnames(annot)) { annot <- cbind(annot, Description=NA)
          chLi <- 1:nrow(annot)
        } else { chLi <- which(is.na(annot[,"Description"]) | nchar(annot[,"Description"]) <1)}
        if(length(chLi) >0) {
          chID <- annot[chLi,"Description"]
          chLi2 <- match(annot[chLi,"Accession"], fasta[,"uniqueIdentifier"])
          if(length(wrMisc::naOmit(chLi2)) >0) annot[chLi,"Description"] <- fasta[chLi2,"proteinName"] 
        }

        ## complete annot / GeneName/GN (if absent...)
        if(!"GeneName" %in% colnames(annot)) { annot <- cbind(annot, GeneName=NA)
          chLi <- 1:nrow(annot)
        } else { chLi <- which(is.na(annot[,"GeneName"]) | nchar(annot[,"GeneName"]) <1)}
        if(length(chLi) >0) {
          chID <- annot[chLi,"GeneName"]
          chLi2 <- match(annot[chLi,"Accession"], fasta[,"uniqueIdentifier"])
          if(length(wrMisc::naOmit(chLi2)) >0) annot[chLi,"GeneName"] <- fasta[chLi2,"GN"] 
        }

        if(debug) {message(fxNa,"eSP2ab"); eSP2ab <- list(specPref=specPref,annot=annot,useColumn=useColumn,suplInp=suplInp,chAnn=chAnn,fasta=fasta,eqiCol=eqiCol, debug=debug)}


        ## complete annot / species (if absent, like with PD)
        ## add species to annot (if not already present)
        chSpNA <- is.na(annot[,"Species"]) 
        if(any(chSpNA) && "OS" %in% colnames(fasta) && !all(is.na(fasta[,"OS"]))) {
          ## complete Species annot from fasta  
          spec <- spc2 <- rep(NA, nrow(annot))
          if(eqiCol[2,1] %in% colnames(annot)) { 
            hasSep <- grep("_[[:alpha:]]",annot[,eqiCol[2,1]])      # "EntryName"
            if(length(hasSep) >0) {    # extract species from col "EntryName" eg "ENO2_YEAST"
              spc2[hasSep] <- sub(".+_","", annot[hasSep,"EntryName"])
              spe2 <- unique(sub("^_","", wrMisc::naOmit(spc2)))
              spe3 <- sapply(spe2, inspectSpeciesIndic, silent=TRUE)
              for(i in 1:length(spe2)) spec[which(spc2 %in% spe2[i])] <- spe3[i]
              if("Species" %in% colnames(annot)) {isNA <- is.na(spec) | length(nchar(spec)) <1 ; if(any(!isNA)) annot[which(!isNA),"Species"] <- spec[which(!isNA)] 
                } else annot <-  cbind(annot, Species=spec) }}
        }

        ## adjust fasta to annot for special case of UPS1
        if("P10636-8" %in% annot[,"Accession"] && "P10636" %in% fasta[,"uniqueIdentifier"]) fasta[grep("P10636",fasta[,"uniqueIdentifier"]),"uniqueIdentifier"] <- "P10636-8"
      }     
      if(debug) {message(fxNa,"eSP2b"); eSP2b <- list(specPref=specPref,annot=annot,useColumn=useColumn,suplInp=suplInp,chAnn=chAnn,fasta=fasta)}

    }     
  } else {if(debug) fasta <- NULL}         ## finish adding fasta  
 
  ## main use of specPref to assign info to column 'SpecType'
  if(length(specPref) > 0) { if(is.list(specPref)) {    # remove NA from specPref
    chNA <- sapply(specPref, is.na)
    if(any(unlist(chNA))) specPref <- sapply(specPref, wrMisc::naOmit)
  } else { chNA  <- is.na(specPref)
    if(all(chNA)) specPref <- NULL else { spNames <- names(specPref[which(!chNA)]); specPref <- as.list(specPref[which(!chNA)]); names(specPref) <- spNames}}}

  if(length(specPref) > 0) {
    if(debug) {message(fxNa,"eSP3"); eSP3 <- list(specPref=specPref,annot=annot,useColumn=useColumn)}
    ## convert  "OS=Saccharomyces cerevisiae" to  "Saccharomyces cerevisiae"
    chPr <- sapply(specPref, function(x) grep("^OS=[[:upper:]][[:lower:]]+", x))          #   look for species indication of type 'OS=Ho...' , need to remove 'OS='
    chPr2 <- sapply(chPr, length) > 0
    if(any(chPr2)) {for(i in which(chPr2)) specPref[[i]] <- sub("^OS=", "", specPref[[i]]) }
    
    ## if specPref$mainSpecies not existing : copy/convert content specPref$matrix to specPref$mainSpecies (redundant)
    if("matrix" %in% names(specPref) && !"mainSpecies" %in% names(specPref)) specPref$mainSpecies <- specPref$matrix
    chNa <- c("mainSpecies","conta") %in% names(specPref)
    if(any(!chNa) && !silent) message(fxNa," ",wrMisc::pasteC(c("mainSpecies","conta")[which(!chNa)], quoteC="'")," Seem absent from 'specPref' !")
    if(debug) {message(fxNa,"eSP3a"); eSP3a <- list(specPref=specPref,annot=annot,eqiCol=eqiCol,useColumn=useColumn,suplInp=suplInp,fasta=fasta)}
    
    ## add species to annot (if not already present)
    if(!"Species" %in% colnames(annot)) annot <- cbind(annot, Species=NA)
    chSpNA <- is.na(annot[,"Species"]) 
    
        
    ## need to distinguish spike-protein collections from full species mix
    spikeTy <- NA
    if("IniSpike" %in% names(specPref)) {
      if(length(specPref$IniSpike)==1 && grepl("^UPS", specPref$IniSpike)) spikeTy <- "collection" else {
        if(length(specPref$IniSpike) >=1 && all(grepl("^[[:upper:]][[:lower:]]+ [[:upper:]]+$", specPref$IniSpike))) spikeTy <- "entireSpecies"
      }
    } else if("spike" %in% names(specPref)) {
      if(length(specPref$spike)==1 && grepl("^UPS", specPref$spike)) spikeTy <- "collection" else {
        if(length(specPref$spike) >=1 && all(grepl("^[[:upper:]][[:lower:]]+ [[:upper:]]+$", specPref$spike))) spikeTy <- "entireSpecies"
      }
    }

    ## final assigning content of annot[,"SpecType"] 
    #old# multGrepFields <- c("mainSpecies","spike","conta")   # fields from specPref to consider (now has redundant & other ...)
    rmSpecPrefFi <- c("matrix","sampleNames","IniMatrix","IniSpike","SpeciesSpike","EntryNameSpike")    # rather define fields which are technical copies (to exclude some terms)    
     
    multGrepFields <- -1*which(names(specPref) %in% rmSpecPrefFi)
    if(length(multGrepFields) <1) multGrepFields <- 1:length(specPref)
    if(debug) {message(fxNa,"eSP3aa"); eSP3aa <- list(specPref=specPref,annot=annot,eqiCol=eqiCol,useColumn=useColumn,suplInp=suplInp,fasta=fasta,rmSpecPrefFi=rmSpecPrefFi,multGrepFields=multGrepFields)}
        ##
    .MultGrep2 <- function(pat, y) { y <- as.matrix(y); z <- if(length(pat)==1) grepl(pat, y) else rowSums(sapply(pat, grepl, y)) >0
      if(length(dim(y)) >1) rowSums(matrix(z, ncol=ncol(y))) >0 else z }  # (multiple) grepl() when length of pattern 'pat' >0
    mulP <- lapply(specPref[multGrepFields], .MultGrep2, if("collection" %in% spikeTy) annot[,which(colnames(annot) %in% c("Accession","Master.Protein.Accessions")), drop=FALSE] else annot)

    combMulP <- function(rmWord, mul) {   # combine elements of 'mul' with 'spike'
      if(all(c("spike",rmWord) %in% names(mul))) {
        mul[["spike"]] <- mul[["spike"]] | mul[[rmWord]]   # combine
        mul <- mul[-which(names(mulP)==rmWord)]             # remove this part
      }
      mul }
    mulP <- combMulP("EntryNameSpike", mulP)  
    mulP <- combMulP("EntryNumberSpike", mulP)  

   chLe <- sapply(mulP, length)
    ## assign to annot   
    if(debug) {message(fxNa,"eSP3ab"); eSP3ab <- list(specPref=specPref,annot=annot,eqiCol=eqiCol,useColumn=useColumn,suplInp=suplInp,fasta=fasta,mulP=mulP,chLe=chLe,multGrepFields=multGrepFields)}
    ## check  which(mulP$spike) if correct    18may25
    
    if(!"SpecType" %in% colnames(annot)) annot <- cbind(annot, SpecType=NA)      # check & ajdust
    
    ## finally fill column 'SpecType'
    for(i in which(chLe >0)) { 
      markLi <- which(as.logical(mulP[[i]]))   # which lines should get new/corrected label
      if(length(markLi) >0) { 
        if("collection" %in% spikeTy) {
          if(names(mulP)[i] !="mainSpecies") annot[markLi,"SpecType"] <- names(mulP)[i] 
        } else if("collection" %in% spikeTy) {  # entire Species collections
          annot[markLi,"SpecType"] <- names(mulP)[i] 
        }        
        ## old
        #if("spike" %in% names(chLe)[i] && grepl("^UPS",specPref$IniSpike) && "Master.Protein.Accessions" %in% colnames(annot)) {   # all are so far NA
        #  ##  Accession doesn't fit to value from getUPS1acc : split annot[,"Master.Protein.Accessions"] and replace Accession by one fitting
        #  chSp <- markLi[which(!annot[markLi,"Accession"] %in% getUPS1acc()$acFull)]   #
        #  if(length(chSp) >0) for(j in chSp) annot[j,"Accession"] <- names(wrMisc::naOmit(sapply(gsub(" ","", unlist(strsplit(annot[j,"Master.Protein.Accessions"],";"))), match, getUPS1acc()$acFull )))
        #}
        #annot[markLi,"SpecType"] <- names(mulP)[i] 
      }
    }  
    ## fix problem with specPref$mainSpecies not found    
    if(length(specPref$mainSpecies) >0 && !any("mainSpecies" %in% annot[,"SpecType"], na.rm=TRUE)) {
      isMain <- which(annot[,"Species"] %in% specPref$mainSpecies)
      if(length(isMain) >0)  annot[isMain,"SpecType"] <- "mainSpecies"
    }  
    if(debug) {message(fxNa,"eSP4"); eSP4 <- list(specPref=specPref,annot=annot,useColumn=useColumn,suplInp=suplInp,mulP=mulP, fasta=fasta)}
    rm(mulP)
  }
  annot
}



#' Get Matrix With UniProt Abbreviations For Selected Species As Well As Simple Names
#'
#' This (low-level) function allows accessing matrix with UniProt abbreviations for species frequently used in research.
#' This information may be used to harmonize species descriptions or extract species information out of protein-names.
#'
#' @return This function returns a 2-column matrix with species names
#' @seealso used eg in \code{readProtDiscovererFile},  \code{\link{readMaxQuantFile}}, \code{\link{readProlineFile}}, \code{\link{readFragpipeFile}}
#' @examples
#' .commonSpecies()
#' @export
.commonSpecies <- function() {
  ## matrix with UniProt abbreviations for common species
  cbind(ext=c("_HUMAN","_MOUSE","_RAT","_CAVPO","_PIG","_BOVIN","_RABIT", "_SHEEP",
      "_CHICK", "_CAEEL", "_DROME","_DANRE",
      "_XENLA","_AMBME", "_ARATH","_SOLTU","_BETVV",  "_YEAST", "_ECOLI","_MYCTU"),
    name=c("Homo sapiens","Mus muscullus","Rattus norvegicus","Cavia porcellus","Sus scrofa","Bos taurus","Oryctolagus cuniculus","Ovis aries",
      "Gallus gallus", "Caenorhabditis elegans","Drosophila melanogaster","Danio rerio", 
      "Xenopus laevis","Ambystoma mexicanum",
      "Arabidopsis thaliana","Solanum tuberosum","Beta vulgaris",  "Saccharomyces cerevisiae",  "Escherichia coli", "Mycobacterium tuberculosis"),
    simple=c("Human","Mouse","Rat","Pig","Guinea pigs", "Cow","Rabit", "Sheep",
      "Chicken", "Celegans","Droso","Zebrafish",
      "Frog","Axolotl",  "Arabidopsis","Potato","Sugar beet",   "Yeast","Ecoli","Mtuberculosis")  )
}
    


