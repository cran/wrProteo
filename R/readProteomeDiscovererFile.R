#' Read Tabulated Files Exported By ProteomeDiscoverer At Protein LevelrSM2
#'
#' Protein identification and quantification results from
#' \href{https://www.thermofisher.com/order/catalog/product/OPTON-30812}{Thermo ProteomeDiscoverer}
#' which were exported as tabulated text can be imported and relevant information extracted.
#' The final output is a list containing 3 elements: \code{$annot}, \code{$raw} and optional \code{$quant},
#' or returns data.frame with entire content of file if \code{separateAnnot=FALSE}.
#'
#' @details
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
#' dataPD <- readProteomeDiscovererFile(file=fiNa, path=path1, suplAnnotFile=FALSE)
#' summary(dataPD$quant)
#'
#' @export
readProteomeDiscovererFile <- function(fileName, path=NULL, normalizeMeth="median", sampleNames=NULL, read0asNA=TRUE, quantCol="^Abundance", annotCol=NULL, contamCol="Contaminant",
  refLi=NULL, separateAnnot=TRUE, FDRCol=list(c("^Protein.FDR.Confidence","High"), c("^Found.in.Sample.","High")), gr=NULL, sdrf=NULL, suplAnnotFile=TRUE,
  groupPref=list(lowNumberOfGroups=TRUE), specPref=c(conta="CON_|LYSC_CHICK", mainSpecies="OS=Homo sapiens"), plotGraph=TRUE, wex=1.6, titGraph="Proteome Discoverer",
  silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## read ProteomeDiscoverer exported txt

  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readProteomeDiscovererFile")
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
  if(!grepl("\\.txt$|\\.txt\\.gz$", fileName)) message(fxNa,"Trouble ahead, expecting tabulated text file (the file'",fileName,"' might not be right format) !!")

  ## note : reading sample-setup from 'suplAnnotFile' at this place won't allow comparing if number of  samples/columns corresponds to data; do after reading main data
  if(debug) message(fxNa,"rdn0 .. Ready to read", if(length(path) >0) c(" from path ",path[1])," the file  ",fileName[1])


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
    #if(length(annotCol) <1) annotCol <- c("Protein.ID","Description","Gene","Gene..name","Contaminant","Sum.PEP.Score","Coverage....", "X..Peptides","X..PSMs","X..Unique.Peptides", "X..AAs","MW..kDa.","Marked..as","Marked as")
    ## default as R-friendly (convert standard cols later to this format)
    if(length(annotCol) <1) annotCol <- c("Accession","Description","Gene","Gene.Name","Marked.as", "Number.of.Peptides","Number.of.PSMs","Number.of.Unique.Peptides","Number.of.AAs","Coverage.in.Percent")
    # alse  ??    "Exp.q.value.Combined","Sum.PEP.Score"
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
    #annot[,"iniDescription"] <- annot[,"Description"]
    ## add other cols form annotCol

    if(debug) { message(fxNa,"rpd1b "); rpd1b <- list() }
    #notAnnCol <- c( "Accession","Description","Gene","Number.of.Peptides","Number.of.PSMs","Number.of.Unique.Peptides")
    if(length(annotCol) >2) { chSupCol <- match(annotCol[-(1:2)], colnames(tmp))
      if(sum(is.na(chSupCol)) < length(annotCol) -2) annot <- cbind(annot, tmp[,wrMisc::naOmit(chSupCol)])
      if(!silent) message(fxNa,"Adding supl annotation-columns ",wrMisc::pasteC(annotCol[wrMisc::naOmit(chSupCol)], quoteC="'")) }
    ## note 'EntryName' (eg 'UBE2C_HUMAN' not avail in PD)
    ## note 'GeneName' (eg 'UBE2C' can be extracted out of 'Description' after GN=
    ## extract GN
    if(debug) { message(fxNa,"rpd1c "); rpd1c <- list() }

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
          ## covert standard output to R-friendly
          colnames(tmp[1:19]) <-  sub("^X\\.\\.","Number.of.", colnames(tmp[1:19]))  #   sub("\\.\\.\\.$",".in.Percent",
          ## define table of specific terms to substitute ..
          subst=cbind(std=c("Exp.q.value.Combined","Coverage....","MW..kDa.","calc..pI"),
            rfriendly=c("Exp..q.value..Combined","Coverage.in.Percent","MW.in.kDa","calc.pI"))
          chSubst <- match(subst[,1], colnames(tmp)[1:19])
          if(any(!is.na(chSubst))) colnames(tmp)[wrMisc::naOmit(chSubst)] <- subst[wrMisc::naOmit(match(colnames(tmp)[1:19], subst[,1])), 2]
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
    if(length(chSpe) >0)  annot[chSpe,"Species"] <- sub(" [[:upper:]]{2}=.*","", sub(".* OS=","", annot[chSpe,"Description"]))

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
      if(all(chNA)) {                          # no direct match found

        ch1 <- match(paste0(quantCol,".F.Sample"), sub("\\.F[[:digit:]]+\\.Sample",".F.Sample",colnames(tmp)))    # match to simplified "Cancer01.F.Sample" instead of "Cancer01.F1.Sample"
        chNA <- is.na(ch1)
        if(all(chNA)) {                     # no composed match found
          ## run grep on each indiv quantCol
          message(fxNa,"grep on each indiv quantCol - Not yet developed")
          ## develop further ?

          stop("None of samples found !!")
        }
      }
      if(any(chNA)) { message(fxNa,"Note : ",sum(chNA)," out of ",length(chNA)," samples NOT found !")
        ch1 <- ch1[which(!chNA)]
      }
    } else {    ## single value => use as search pattern
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
              if(length(ch2) >0 & !all(ch2)) ch1 <- ch1[which(!ch2)]
            } else {
              ch1 <- grep(quantCol, colnames(tmp))      # use directly as pattern (may find too many)
              if(length(ch1) >0) { ch2 <- grepl(excluPat, colnames(tmp)[ch1])    # avoid  '.Normalized' (if not concerning all)
                if(length(ch2) >0 & !all(ch2)) ch1 <- ch1[which(!ch2)]
              } else {
                warning(fxNa,"Unable to find any matches to '",quantCol,"' !") }
            }
          }
        }
      if(debug) message(fxNa,"Found ",length(ch1)," quantitation-columns", if(length(ch1) >0) c(" (eg ",wrMisc::pasteC(colnames(tmp)[utils::head(ch1)], quoteC="'"),")"))
      }

    }
    quantColInd <- ch1
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
    } }
    if(debug) { message(fxNa,"rpd11b .. dim abund ",nrow(abund)," and ",ncol(abund)); rpd11b <- list(annot=annot,tmp=tmp,abund=abund,sampleNames=sampleNames,specPref=specPref,annotCol=annotCol,Rfriendly=Rfriendly,contamCol=contamCol,PSMCol=PSMCol,PepCol=PepCol,infoDat=infoDat)}

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
    if(debug) { message(fxNa,"rpd13 .. dim quant: ", nrow(quant)," li and  ",ncol(quant)," cols; colnames : ",wrMisc::pasteC(colnames(quant))," "); 
      rpd13 <- list(annot=annot,tmp=tmp,abund=abund,quant=quant,sampleNames=sampleNames,specPref=specPref,annotCol=annotCol,Rfriendly=Rfriendly,contamCol=contamCol,PSMCol=PSMCol,PepCol=PepCol,infoDat=infoDat, refLi=refLi)}

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
      refLi=refLi, refLiIni=refLiIni, tit=titGraph, las=NULL, silent=silent, callFrom=fxNa, debug=debug)

    ## meta-data
    notes <- c(inpFile=paFi, qmethod="ProteomeDiscoverer", qMethVersion=if(length(infoDat) >0) unique(infoDat$Software.Revision) else NA,
    	rawFilePath= if(length(infoDat) >0) infoDat$File.Name[1] else NA, normalizeMeth=normalizeMeth, call=match.call(),
      created=as.character(Sys.time()), wrProteo.version=utils::packageVersion("wrProteo"), machine=Sys.info()["nodename"])

    ## final output
    if(isTRUE(separateAnnot)) list(raw=abund, quant=quant, annot=annot, counts=counts, sampleSetup=setupSd, quantNotes=parametersD, notes=notes) else data.frame(quant,annot)
  }
}




#' Additional/final chek and adjustments to sample-order after readSampleMetaData()
#'
#' This (low-level) function performs an additional/final chek & adjustments to sample-names after readSampleMetaData()
#'
#' @param abund (matrix or data.frame) abundance data, only the colnames will be used
#' @param setupSd (list) describing sammple-setup, typically produced by   from package wrMisc
#' @param gr (factor) optional custom information about replicate-layout, has priority over setuoSd
#' @param sampleNames (character) custom sample-names, has priority over abund and setuoSd
#' @param quantMeth (character) 2-letter abbreviation of name of quantitation-software (eg 'MQ')
#' @param silent (logical) suppress messages
#' @param debug (logical) display additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns an enlaged/updated list 'setupSd' (set setupSd$sampleNames,  setupSd$groups)
#' @seealso used in \code{readProtDiscovererFile},  \code{\link{readMaxQuantFile}}, \code{\link{readProlineFile}}, \code{\link{readFragpipeFile}}
#' @examples
#'set.seed(2021)
#' @export
.checkSetupGroups <- function(abund, setupSd, gr=NULL, sampleNames=NULL, quantMeth=NULL, silent=FALSE, callFrom=NULL, debug=FALSE) {
  ## additional/final chek & adjustments to sample-names after readSampleMetaData()
  ## examine
  ## returns enlaged/updated list 'setupSd' (set setupSd$sampleNames,  setupSd$groups)

  ## assume 'abund' to be valid matrix of data.frame !!
  ## 'setupSd' as list produced by readSampleMetaData()
  ## 'gr' .. (char vector or factor) designating who should be considered as replicate/same level, as optional custom entry (has prior over automatic groups/levels)
  ## 'sampleNames' .. optional custom entry for sample names (has prior over automatic names)
  ## 'quantMeth' .. (char vect, length=1) design method via abbreviation like 'PD' for ProteomeDiscoverer, 'MQ' for MaxQuant, 'PL' for Proline, etc (used for automatic trimming of default column/sample-names)
  fxNa <- wrMisc::.composeCallName(callFrom, newNa=".checkSetupGroups")
  delPat <- "_[[:digit:]]+$|\\ [[:digit:]]+$|\\-[[:digit:]]+$"             # remove enumerators, ie tailing numbers after separator
  rawExt <- "\\.raw$|\\.Raw$|\\.RAW$"     # paste(paste0("\\.",c("Raw","raw","RAW"),"$"), collapse="|")
  .corPathW <- function(x) gsub("\\\\", "/", x)

  .extrColNames <- function(abun, meth, silent=FALSE, callFrom=NULL, debug=FALSE) {   ## get sampleNames  from colnames of abund & clean
    fxNa <- wrMisc::.composeCallName(callFrom, newNa=".extrColNames")
      # abun=abund; meth=quantMeth
    grou <- grou2 <- NULL
    colNa2 <- sub(rawExt,"", .corPathW(colnames(abun)))
    colNa <- basename(colNa2)                   # trim to filename
    chDu <- duplicated(colNa)
    if(debug) {message(fxNa,"eCN1"); eCN1 <- list(abun=abun,meth=meth,colNa2=colNa2,colNa=colNa,chDu=chDu)}
    if(any(chDu)) {
     if(debug) message(fxNa,"Some stripped filenames appear duplicated, need to keep with path ...")
       colNa <- colNa2
       chDu <- duplicated(colNa) }
    if(any(chDu) && !silent) message(fxNa,"Some filenames appear duplicated !!   (eg  ",wrMisc::pasteC(utils::head(unique(colNa[which(chDu)]), 3))," )")
    if(debug) {message(fxNa,"eCN2"); eCN2 <- list(abun=abun,meth=meth,colNa2=colNa2,colNa=colNa,chDu=chDu)}

    if(is.null(colNa) || sum(duplicated(colNa)) ==length(colNa) -1) {  sampleNames <- NULL
      if(debug) message(fxNa,"All colnames of abund are empty or identical ! (can't use to figure out pattern of replicates/levels)")
    } else {
      ## colNa seem usable
      if(any(c("FP") %in% meth)) colNa <- sub("MaxLFQ Intensity$|Intensity$","", colNa)
      if(any(c("MQ") %in% meth)) colNa <- sub("^LFQ Intensity","", colNa)
      if("PL" %in% meth) colNa <- sub("^abundance|^Abundance","", colNa)
      colNa <- sub(" +$|\\.+$|_+$|\\-+$","", colNa)         # remove tailing separators (' ','.','_','-')
      colNa <- sub("^ +|^\\.+|^_+|^\\-+","", colNa)         # heading tailing separators
      chDu <- duplicated(colNa)
      if(!silent && any(chDu)) message(fxNa,"NOTE : ",sum(chDu)," DUPLICATED colnames for abund !!   (eg  ",wrMisc::pasteC(utils::head(unique(colNa[which(chDu)]), 3))," )")
      if(debug) {message(fxNa,"eCN3"); eCN3 <- list()}

      ## now address group names/levels
      if("PD" %in% meth) {
        grou2 <- sub("rep[[:digit:]]+$","", colnames(abun))
        grou2 <- sub(" +$|\\.+$|_+$|\\-+$","", grou2)       # remove tailing separators (' ','.','_','-')
        if(all(grepl("^\\.F[[:digit:]]+\\.Sample$", colnames(abun)))) grou2 <- NULL        #  PD default names like   '.F1.Sample', '.F2.Sample' etc
      } else {
        if("PL" %in% meth) colNa <- sub("\\.mzDB\\.t\\.xml$","", colNa)
        colNa <- sub("\\.raw$|\\.RAW$","", colNa)
        sep <- c("_","\\-","\\.")
         #sub("\\.mzDB|\\.t\\.xml","",colNa)
        rmTx <- paste0(c("Sample","Samp","Replicate","Rep"),"$")
        rmTx <- paste(paste0(rep(sep, each=length(rmTx)), rep(rmTx, length(sep))), collapse="|")
        grou2 <- sub(rmTx,"",sub(delPat,"", colNa))                               # remove tailing enumerators..
      }
      if(debug) message(fxNa,"Based on colnames(abund) : ",length(unique(grou2))," levels for ",ncol(abun)," samples")
    }
    if(debug) {message(fxNa,"eCN4 done"); eCN4 <- list()}
    list(sampleNames=colNa, grou=grou2) }

  extrSamNaSetup <- function(setS, meth) {
    ## extract (use-given) sampleNames out of setupSd$annotBySoft (colnames may depend on quant-method)
    if(!any(c("PD","MQ","PL","FP") %in% meth, na.rm=TRUE)) meth <- "other"
    switch(meth, PD = NULL,
      MQ = setupSd$annotBySoft$Experiment , PL = setupSd$annotBySoft$Experiment, FP = setupSd$annotBySoft$Experiment, other=NULL)}

  defColNa <- function(colN, meth) {                             # check if colN may represent default colnames (ie not useful since wo any indication about samples)
    ## note MQ : requires setExperiment to be defined by user, set each sample as different for getting quant by sample !
    ## note FP : requires setExperiment to be defined by user (defining different bioreplictates sufficient for getting quant by sample)
    if(!any(c("PD","MQ","PL","FP") %in% meth, na.rm=TRUE)) meth <- "other"
    switch(meth, PD = all(grepl("F[[:digit:]]+\\.Sample$", colN), na.rm=TRUE),
      MQ = FALSE , PL = FALSE, FP = FALSE, other=FALSE)}

  ## finish groups of replicates & annotation setupSd
  if(debug) { message(fxNa,"cSG0"); cSG0 <- list(gr=gr,abund=abund,sampleNames=sampleNames, setupSd=setupSd) }          # sampleNames=sampleNames

  ## grou2 ... colnames modified to pattern ; grou .. pattern as index, + names (level names)

  colNa <- grou <- grou2 <- NULL
  iniSaNa <- iniGr <- FALSE
  ## if valid user-defied sampleNames is given => use
  if(length(sampleNames) != ncol(abund)) {
    if(debug && length(sampleNames) >0) message(fxNa,"Invalid entry of 'sampleNames' (length= ",length(sampleNames),"  but  ",ncol(abund)," expected)  ...ignoring")
    sampleNames <- NULL
  } else { iniSaNa <- TRUE
    setupSd$sampleNames <- sampleNames }
  if(debug) {message(fxNa,"cSG1"); cSG1 <- list()}

  ## if valid user-defied grouping is given => use
  if(length(gr) != ncol(abund)) {
    if(debug && length(gr) >0) message(fxNa,"Invalid entry of 'gr' (length= ",length(gr),"  but  ",ncol(abund)," expected)  ...ignoring")
    gr <- NULL
  } else {
    ## check if setupSd has 'prioritized' grouping
    if("groups" %in% names(setupSd)) {
      setupSd$level <- setupSd$groups
    } else {
      iniGr <- TRUE
      setupSd$level <- match(gr, unique(gr))
      setupSd$groups <- names(setupSd$level) <- gr } }
  if(debug) {message(fxNa,"cSG2"); cSG2 <- list()}

  defaultColNa <- defColNa(colN=colnames(abund), meth=quantMeth)
  saNa <- .extrColNames(abund, meth=quantMeth, silent=silent,callFrom=fxNa,debug=debug)      # sampleNames  from colnames of abund & clean
  if(debug) {message(fxNa,"cSG3"); cSG3 <- list(sampleNames=sampleNames,gr=gr,abund=abund,iniSaNa=iniSaNa,iniGr=iniGr,setupSd=setupSd,defaultColNa=defaultColNa,saNa=saNa)}

  ## sampleNames/colnames : use orig colnames if avail  (priority to colnames)
  #colNa <- if(defaultColNa) NULL else sub(rawExt,"", colnames(abund))  # redundant to saNa$sampleNames

  if(length(setupSd$level) ==ncol(abund)) gr <- setupSd$lev <- setupSd$level else {
    if(length(setupSd$groups) ==ncol(abund)) gr <- setupSd$lev <- setupSd$groups }
  if(debug) { message(fxNa,"cSG3b"); cSG3b <- list()}

  #if("lev" %in% names(setupSd)) {
  if(length(setupSd$lev) ==ncol(abund)) {
    ## thus, we do have setupSd
    ## sampleNames/colnames : use orig colnames if avail  (priority to colnames)
    if(!iniSaNa && !defaultColNa && length(saNa$sampleNames) ==ncol(abund)) sampleNames <- saNa$sampleNames
    if(!iniSaNa && length(setupSd$sampleNames) ==ncol(abund)) sampleNames <- setupSd$sampleNames           # setupSd$sampleNames has prioroty (if defined)
      ## compare grouping of orig colnames to sdrf ?

    if(length(sampleNames) !=ncol(abund)) {
      ## get sampleNames from setupSd (as far as possible)
      saNa2 <- extrSamNaSetup(setupSd, quantMeth)           # from setupSd$annotBySoft (by quant method)
      sampleNames <- if(length(saNa2) ==ncol(abund)) saNa2 else if(length(setupSd$sdrfDat$comment.data.file.)==ncol(abund)) sub(rawExt,"", setupSd$sdrfDat$comment.data.file.)
    }

    ## now for gr
    if(length(gr) != ncol(abund)){
      if(!iniSaNa && !defaultColNa && length(saNa$grou) ==ncol(abund)) gr <- saNa$grou
      if(!iniSaNa && length(setupSd$groups) ==ncol(abund)) gr <- setupSd$groups           # setupSd$groups has prioroty (if defined)
    }
    if(debug) {message(fxNa,"cSG4a"); cSG4a <- list()}

  } else {
    ## (no setupSd)
    ## get sampleNames from abund (as far as possible)
    if(iniSaNa && length(sampleNames) != ncol(abund)) {
      if(defaultColNa) { if(length(gr)==ncol(abund)) sampleNames <- wrMisc::correctToUnique(gr)   # case of PD : use gr if suitable
      } else sampleNames <- saNa$sampleNames                                                      # other use colnames of abund
    }
    ## now for gr  (no setupSd)
    if(!iniGr) {
      if(!defaultColNa) {         ## standard case (eg MQ)
        ## guess gr from colnames
        if(debug) message(fxNa,"Guess 'gr' from colnames,  ",quantMeth,"")
        gr <- saNa$grou
      } else {     ##  (PD:) no way to guess groups
        if(!silent) message(fxNa,"Difficulty to identify groups of replicates (no setupSd) in case of absence of metadata by method ",quantMeth,"")
      }
    }
    if(debug) {message(fxNa,"cSG4b"); cSG4b <- list()}

  }
  ## case of PD : check if  setupSd$annotBySoft$File.Name  usable
  if(defaultColNa && length(sampleNames) != ncol(abund)) {
    chOr <- NA                           # initialize
    if(length(setupSd$annotBySoft$File.Name) ==ncol(abund)) chOr <- match(setupSd$annotBySoft$File.Name, setupSd$sdrfDat$comment.data.file.)
    if(!any(is.na(chOr))) {
      if(!all(chOr ==1:ncol(abund), na.rm=TRUE)) {
        sampleNames <-  sub("\\.raw$|\\.RAW$","", setupSd$annotBySoft$File.Name)
        ## try extracting pattern of replicates
        colNa <- sub("\\.raw$|\\.RAW$","", basename(sub(rawExt,"", .corPathW(colnames(sampleNames)))))
        sep <- c("_","\\-","\\.")
        rmTx <- paste0(c("Sample","Samp","Replicate","Rep"),"$")
        rmTx <- paste(paste0(rep(sep, each=length(rmTx)), rep(rmTx, length(sep))), collapse="|")
        gr <- sub(rmTx,"", sub(delPat,"", colNa))                               # remove tailing enumerators..
        if(debug) message(fxNa,"cSG4c   Method ",quantMeth," : Extracted ",length(unique(gr)), " groups of replicates based on meta-data")
      }
  } }


  if(debug) {message(fxNa,"cSG5"); cSG5 <- list(sampleNames=sampleNames,gr=gr,abund=abund,iniSaNa=iniSaNa,iniGr=iniGr,setupSd=setupSd,defaultColNa=defaultColNa,saNa=saNa)}
  if(length(sampleNames) != ncol(abund) && defaultColNa & !silent) message(fxNa,"Still UNABLE to find suitable colnames")
  if(length(gr) != ncol(abund) && defaultColNa && debug) message(fxNa,"Still UNABLE to find suitable groups")
  ## set result to object
  if(!is.list(setupSd)) { if(length(setupSd) >0) warning(fxNa,"BIZZARE format of 'setupSd', it's content will be lost")
    setupSd <- list()}
  setupSd$sampleNames <- sampleNames
  setupSd$groups <- gr
  setupSd
  }




#' Generic plotting of density distribution for quantitation import-functions
#'
#' This (low-level) function allows (generic) plotting of density distribution for quantitation import-functions
#'
#' @param abund (matrix or data.frame) abundance data, will be plottes as distribution
#' @param quant (matrix or data.frame) optional additional abundance data, to plot 2nd distribution, eg of normalized data
#' @param custLay (matrix) describing sammple-setup, typically produced by
#' @param normalizeMeth (charactert, length=1) name of normalization method (will be displayed in title of figure)
#' @param softNa (charactert, length=1) 2-letter abbreviation of name of quantitation-software (eg 'MQ')
#' @param refLi (integer) to display number reference lines
#' @param refLiIni (integer) to display initial number reference lines
#' @param tit (character) custom title
#' @param las (integer) indicate orientation of text in axes
#' @param silent (logical) suppress messages
#' @param debug (logical) display additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns logical value (if data were valid for plotting) and produces a density dustribution figure (if data were found valid)
#' @seealso used in \code{readProtDiscovererFile},  \code{\link{readMaxQuantFile}}, \code{\link{readProlineFile}}, \code{\link{readFragpipeFile}}
#' @examples
#' set.seed(2018);  datT8 <- matrix(round(rnorm(800)+3,1), nc=8, dimnames=list(paste(
#'   "li",1:100,sep=""), paste(rep(LETTERS[1:3],c(3,3,2)),letters[18:25],sep="")))
#' .plotQuantDistr(datT8, quant=NULL, refLi=NULL, tit="Synthetic Data")
#' @export
.plotQuantDistr <- function(abund, quant, custLay=NULL, normalizeMeth=NULL, softNa=NULL, refLi=NULL, refLiIni=NULL, tit=NULL, las=NULL, silent=FALSE, callFrom=NULL, debug=FALSE) {
  ## generic plotting of densirt distribution for quantitation import-functions
  ## assume 'abund' (raw, non-normalized) and 'quant' (final normalized) to be valid matrix of data.frame !!
  ## 'custLay' ..(matrix) for layout()
  ## 'normalizeMeth' (character)
  ## 'softNa' .. (char vect, length=1) design method, used in display only
  ## 'refLi' .. (integer) reference line
  ## 'refLiIni' .. (integer) initial reference line(s)

  fxNa <- wrMisc::.composeCallName(callFrom, newNa=".plotQuantDistr")
  oparMar <- graphics::par("mar")                # only if figure might be drawn
  on.exit(graphics::par(mar=oparMar))            # restore old mar settings
  if(debug) {message(fxNa,"pQD1 .. length custLay ", length(custLay) )}
  ## check if abund & quant are redundant (while normalizeMeth="none")
  if("none" %in% normalizeMeth && length(abund) >0 && (identical(quant,abund) || identical(quant,log2(abund)))) {
    quant <- NULL
    if(debug) message(fxNa,"No need for 2nd plot, 'abund' and 'quant' are identical ..")}

  if(length(custLay) >0) graphics::layout(custLay) else {
    if(!identical(normalizeMeth,"none") && length(quant) >0) { ch1 <- try(graphics::layout(1:2), silent=TRUE)
      if(inherits(ch1, "try-error")) message(fxNa,"Problem with figure, Need to restore layout ..")}}
  plotGraph <- TRUE                                         # useful ? (to report if plot can be drawn as output ?)
  if(length(las) != 1 || !is.integer(las)) las <- if(ncol(abund) >7 || stats::median(nchar(colnames(abund)), na.rm=TRUE) >8 ) 2 else 1
  ch1 <- try(graphics::par(mar=c(3, 3, 3, 1)), silent=TRUE)                          # mar: bot,le,top,ri
  if(inherits(ch1, "try-error")) message(fxNa,"Problem with figure, Need to restore mar ..")
  if(is.null(tit)) tit <- paste(softNa," quantification ")
  titSu <- if(length(refLi) >0) paste0(c(if(length(refLiIni) ==1) c("'",refLiIni,"'") else c(" by ",length(refLi)," selected lines")),collapse="")  else NULL #
  usePa <- c("wrGraph","sm")
  misPa <- !sapply(usePa, function(pkg) requireNamespace(pkg, quietly=TRUE))
  if(debug) { message(fxNa,"pQD2 .. misPa ", wrMisc::pasteC(misPa,quoteC="'") )}

  if(any(misPa, na.rm=TRUE)) {
    if(!silent) message(fxNa,"Please install package(s) ", wrMisc::pasteC(usePa[which(misPa)])," from CRAN for drawing vioplots")
    ## wrGraph not available : simple boxplot
    ch1 <- try(graphics::boxplot(log2(abund), main=paste(tit,"(initial)", sep=" "), las=1, outline=FALSE), silent=TRUE)
    if(inherits(ch1, "try-error")) {plotGraph <- FALSE; message(fxNa,"UNABLE to draw boxplot of distribution !!  ", sub("^Error in",":",ch1))} else {
      graphics::abline(h=round(log2(stats::median(abund,na.rm=TRUE))) +c(-2:2), lty=2, col=grDevices::grey(0.6))}
    ## plot normalized
    if(identical(normalizeMeth,"none") | length(quant) <0) {
      if(debug) {message(fxNa,"pQD3 .. dim quant: ", nrow(quant)," li and  ",ncol(quant)," cols; colnames : ",wrMisc::pasteC(colnames(quant))," ")}
      ch1 <- try(graphics::boxplot(quant, main=paste(tit," (",normalizeMeth,"-normalized",titSu,")"), las=1, outline=FALSE), silent=TRUE)
      if(inherits(ch1, "try-error")) message(fxNa,"UNABLE to draw boxplot of normalized distribution !!  ", sub("^Error in",":",ch1)) else {
        graphics::abline(h=round(stats::median(quant, na.rm=TRUE)) +c(-2:2), lty=2, col=grDevices::grey(0.6)) }}

  } else {                                            # wrGraph and sm are available
    if(debug) {message(fxNa,"pQD4  draw vioplotW "  )}
    tit1 <- if(length(quant) >0) paste(tit, "(initial)",sep=" ") else tit
    ch1 <- try(wrGraph::vioplotW(log2(abund), tit=tit1, wex=NULL, silent=silent, callFrom=fxNa), silent=TRUE)
    if(inherits(ch1, "try-error")) {plotGraph <- FALSE; message(fxNa,"UNABLE to plot vioplotW !!  ", sub("^Error in",":",ch1))
    } else graphics::abline(h=round(stats::median(log2(abund), na.rm=TRUE)) +c(-2:2), lty=2, col=grDevices::grey(0.6))
    ## now normalized
    if(debug) {message(fxNa,"pQD5  draw norm vioplotW() ", length(quant) <0)}
    if(length(quant) >0) if(identical(quant, abund)) {quant <- NULL;
      if(!silent) message(fxNa,"Note : 'quant' and 'abund' were found identical, omit 2nd (redundant) graph ..") }
    if(length(quant) >0) {
      if(debug) {message(fxNa,"pQD6  draw vioplotW() for normalized")}            ## now normalized
      tit1 <- if(length(normalizeMeth) ==1) paste(tit,", ",normalizeMeth,"-normalized",titSu) else paste(tit,", ",titSu)
      ch1 <- try(wrGraph::vioplotW(quant, tit=tit1, wex=NULL, silent=silent, callFrom=fxNa), silent=TRUE)
      if(inherits(ch1, "try-error")) { message(fxNa,"UNABLE to plot vioplotW for normalized data !!  ", sub("^Error in",":",ch1))
      } else graphics::abline(h=round(stats::median(quant, na.rm=TRUE)) +c(-2:2), lty=2, col=grDevices::grey(0.6)) }
  }
  plotGraph }



#' Extract additional information to construct colum SpecType
#'
#' This (low-level) function creates the column annot[,'SpecType'] which may help distinguishing different lines/proteins.
#' This information may, for example, be used to normalize only to all proteins of a common backgroud matrix (species).
#' If $mainSpecies or $conta: match to annot[,"Species"], annot[,"EntryName"], annot[,"GeneName"], if length==1 grep in  annot[,"Species"]
#'
#' @param specPref (list) may contain $mainSpecies, $conta ...
#' @param annot (matrix) main protein annotation
#' @param useColumn (factor) columns from annot to use/mine
#' @param suplInp (matrix) additional custom annotation
#' @param silent (logical) suppress messages
#' @param debug (logical) display additional messages for debugging (starting with 'mainSpecies','conta' and others - later may overwrite prev settings)
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a matrix with additional column 'SpecType'
#' @seealso used in \code{readProtDiscovererFile},  \code{\link{readMaxQuantFile}}, \code{\link{readProlineFile}}, \code{\link{readFragpipeFile}}
#' @examples
#'.checkKnitrProt()
#' @export
.extrSpecPref <- function(specPref, annot, useColumn=c("Species","EntryName","GeneName","Accession"), suplInp=NULL, silent=FALSE, debug=FALSE, callFrom=NULL) {
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
  fxNa <- wrMisc::.composeCallName(callFrom, newNa=".extrSpecPref")
  if(debug){ message(fxNa," eSP0"); eSP0 <- list(specPref=specPref,annot=annot,useColumn=useColumn,suplInp=suplInp)}

  if(length(annot) <1 || length(dim(annot)) !=2) stop("invalid 'annot' (must be matrix or data.frame)")
  ## check suplInp & match to useColumn, add to useColumn
  if(length(useColumn) >4 && length(suplInp) >0 && length(dim(suplInp))==2) {
    chAnn <- useColumn[5:length(useColumn)] %in% colnames(suplInp)
    if(any(!chAnn)) useColumn <- c(useColumn[1:4], useColumn[(5:length(useColumn))[which(chAnn)]])
    if(length(useColumn) <5) suplInp <- NULL
  } else suplInp <- NULL
  ## check useColumn
  chAnn <- useColumn[1:min(length(useColumn), 4)] %in% colnames(annot)   # check for length useCol=0  ?? # nolint
  if(any(!chAnn)) stop("Unknown/Non-standard 'annot' (missing colnames ",wrMisc::pasteC(useColumn[which(!chAnn)], quoteC="'"),")")

  if(!"SpecType" %in% colnames(annot)) {annot <- cbind(annot, SpecType=rep(NA, nrow(annot))); if(debug) message(fxNa,"Adding column 'SpecType' to 'annot'")}
  if(length(specPref) > 0) specPref <- specPref[which(sapply(specPref, length) >0)]     # remove empty ..
  if(length(specPref) > 0) if(is.list(specPref)) {    # remove NA from specPref
    chNA <- sapply(specPref, is.na)
    if(any(unlist(chNA))) specPref <- sapply(specPref, wrMisc::naOmit)
  } else { chNA  <- is.na(specPref)
    if(all(chNA)) specPref <- NULL else { spNames <- names(specPref[which(!chNA)]); specPref <- as.list(specPref[which(!chNA)]); names(specPref) <- spNames}}
  if(length(specPref) > 0) {
    if(debug) {message(fxNa,"eSP1"); eSP1 <- list()}
    chNa <- c("mainSpecies","conta") %in% names(specPref)
    if(any(!chNa) && ! silent) message(fxNa," ",wrMisc::pasteC(c("mainSpecies","conta")[which(!chNa)], quoteC="'")," Seem absent from 'specPref' !")

    #old# .MultGrepM <- function(patt, anno) if(length(dim(anno)) >1) apply(anno, 2, .MultGrep2, patt) else  .MultGrep2(anno, patt)
    .MultGrep2 <- function(pat, y) {y <- as.matrix(y); z <- if(length(pat)==1) grepl(pat, y) else rowSums(sapply(pat, grepl, y)) >0
      if(length(dim(y)) >1) rowSums(matrix(z, ncol=ncol(y))) >0 else z }  # (multiple) grepl() when length of pattern 'pat' >0

    mulP <- lapply(specPref, .MultGrep2, annot)
    chLe <- sapply(mulP, length)
    ## assign to annot
    for(i in which(chLe >0)) { markLi <- which(as.logical(mulP[[i]]))
      if(length(markLi) >0)  annot[markLi,"SpecType"] <- names(mulP)[i] }
    if(debug) {message(fxNa,"eSP2"); eSP2 <- list(specPref=specPref,annot=annot,useColumn=useColumn,suplInp=suplInp,mulP=mulP)}
    rm(mulP)
  }
  annot
}



#' Get matrix with UniProt abbreviations for common species
#'
#' This (low-level) function allows accessing matrix with UniProt abbreviations for common species
#' This information maybe used to harmonize species descriptions.
#'
#' @return This function returns a 2-column matrix with species names
#' @seealso used eg in \code{readProtDiscovererFile},  \code{\link{readMaxQuantFile}}, \code{\link{readProlineFile}}, \code{\link{readFragpipeFile}}
#' @examples
#'.commonSpecies()
#' @export
.commonSpecies <- function() {
  ## matrix with UniProt abbreviations for common species
  cbind(c("_HUMAN","_MOUSE","_RAT","_PIG","_BOVIN","_SHEEP","_CAEEL", "_DROME","_YEAST","_ARATH","_ECOLI","_MYCTU"),
      c("Homo sapiens","Mus muscullus","Rattus norvegicus","Sus scrofa","Bos taurus","Ovis aries","Caenorhabditis elegans",
        "Drosophila melanogaster","Saccharomyces cerevisiae", "Arabidopsis thaliana", "Escherichia coli", "Mycobacterium tuberculosis"))
}

