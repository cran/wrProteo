#' Read xlsx, csv or tsv files exported from Proline and MS-Angel
#'
#' Quantification results from Proline \href{http://www.profiproteomics.fr/proline/}{Proline} and MS-Angel exported as xlsx format can be read directly using this function.
#' Besides, files in tsv, csv (European and US format) or tabulated txt can be read, too.
#' Then relevant information gets extracted, the data can optionally normalized and displayed as boxplot or vioplot.
#' The final output is a list containing 6 elements: \code{$raw}, \code{$quant},  \code{$annot}, \code{$counts}, \code{$quantNotes} and \code{$notes}.
#' Alternatively, a data.frame with annotation and quantitation data may be returned if \code{separateAnnot=FALSE}.
#' Note: There is no normalization by default since quite frequently data produced by Proline are already sufficiently normalized.
#' The figure produced using the argument \code{plotGraph=TRUE} may help judging if the data appear sufficiently normalized (distribtions should align).
#'
#' @details
#' This function has been developed using Proline version 1.6.1 coupled with MS-Angel 1.6.1.
#' The classical way of using ths function consists in exporting results produced by Proline and MS-Angel as xlsx file.
#' Besides, other formats may be read, too. This includes csv (eg the main sheet/table of ths xlsx exported file saved as csv).
#' \href{https://github.com/wombat-p}{WOMBAT} represents an effort to automatize quantitative proteomics experiments, using this route
#' data get exported as txt files which can be read, too.
#'
#' @param fileName (character) name of file to read; .xlsx-, .csv-, .txt- and .tsv can be read (csv, txt and tsv may be gz-compressed). Reading xlsx requires package 'readxl'.
#' @param path (character) optional path (note: Windows backslash sould be protected or written as '/')
#' @param normalizeMeth (character) normalization method (for details and options see \code{\link[wrMisc]{normalizeThis}})
#' @param logConvert (logical) convert numeric data as log2, will be placed in $quant
#' @param sampleNames (character) custom column-names for quantification data; this argument has priority over \code{suplAnnotFile}
#' @param quantCol (character or integer) colums with main quantitation-data : precise colnames to extract, or if length=1 content of \code{quantCol} will be used as pattern to search among column-names for $quant using \code{grep}
#' @param annotCol (character) precise colnames or if length=1 pattern to search among column-names for $annot
#' @param remStrainNo (logical) if \code{TRUE}, the organism annotation will be trimmed to uppercaseWord+space+lowercaseWord (eg Homo sapiens)
#' @param pepCountCol (character) pattern to search among column-names for count data of PSM and NoOfPeptides
#' @param trimColnames (logical) optional trimming of column-names of any redundant characters from beginning and end
#' @param refLi (integer) custom decide which line of data is main species, if single character entry it will be used to choose a group of species (eg 'mainSpe')
#' @param separateAnnot (logical) separate annotation form numeric data (quantCol and annotCol must be defined)
#' @param plotGraph (logical or matrix of integer) optional plot vioplot of initial data; if integer, it will be passed to \code{layout} when plotting
#' @param titGraph (character) custom title to plot of distribution of quantitation values
#' @param wex (integer) relative expansion factor of the violin-plot (will be passed to \code{\link[wrGraph]{vioplotW}})
#' @param specPref (character or list) define characteristic text for recognizing (main) groups of species (1st for comtaminants - will be marked as 'conta', 2nd for main species- marked as 'mainSpe',
#'  and optional following ones for supplemental tags/species - maked as 'species2','species3',...);
#'  if list and list-element has multiple values they will be used for exact matching of accessions (ie 2nd of argument \code{annotCol})
#' @param gr (character or factor) custom defined pattern of replicate association, will override final grouping of replicates from \code{sdrf} and/or \code{suplAnnotFile} (if provided)   \code{}
#' @param sdrf (character, list or data.frame) optional extraction and adding of experimenal meta-data: if character, this may be the ID at ProteomeExchange, 
#'   the second element may give futher indicatations for automatic organization of groups of replicates. 
#'   Besides, the output from \code{readSdrf} or a list from \code{defineSamples} may be provided; if \code{gr} is provided, \code{gr} gets priority for grouping of replicates
#' @param suplAnnotFile (logical or character) optional reading of supplemental files produced by quantification software; however, if \code{gr} is provided, \code{gr} gets priority for grouping of replicates;
#'  if \code{TRUE} defaults to file '*InputFiles.txt' (needed to match information of \code{sdrf}) which can be exported next to main quantitation results;
#'  if \code{character} the respective file-name (relative or absolute path)
#' @param groupPref (list) additional parameters for interpreting meta-data to identify structure of groups (replicates), will be passed to \code{readSampleMetaData}.
#'   May contain \code{lowNumberOfGroups=FALSE} for automatically choosing a rather elevated number of groups if possible (defaults to low number of groups, ie higher number of samples per group)
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of messages produced
#' @param debug (logical) display additional messages for debugging
#' @return This function returns a list with \code{$raw} (initial/raw abundance values), \code{$quant} with final normalized quantitations, \code{$annot} (columns ), \code{$counts} an array with 'PSM' and 'NoOfPeptides', \code{$quantNotes} and \code{$notes}; or a data.frame with quantitation and annotation if \code{separateAnnot=FALSE}
#' @seealso \code{\link[utils]{read.table}}
#' @examples
#' path1 <- system.file("extdata", package="wrProteo")
#' fiNa <- "exampleProlineABC.csv.gz"
#' dataABC <- readProlineFile(path=path1, file=fiNa)
#' summary(dataABC$quant)
#' @export
readProlineFile <- function(fileName, path=NULL, normalizeMeth="median", logConvert=TRUE, sampleNames=NULL, quantCol="^abundance_",
  annotCol=c("accession","description","is_validated","protein_set_score","X.peptides","X.specific_peptides"), remStrainNo=TRUE,
  pepCountCol=c("^psm_count_","^peptides_count_"), trimColnames=FALSE, refLi=NULL, separateAnnot=TRUE, plotGraph=TRUE, titGraph=NULL,
  wex=2, specPref=c(conta="_conta\\|", mainSpecies="OS=Homo sapiens"), gr=NULL, sdrf=NULL, suplAnnotFile=TRUE, groupPref=list(lowNumberOfGroups=TRUE), silent=FALSE, callFrom=NULL, debug=FALSE) {
  ## 'quantCol', 'annotCol' (character) exact col-names or if length=1 pattern to search among col-names for $quant or $annot
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readProlineFile")
  oparMar <- if(plotGraph) graphics::par("mar") else NULL       # only if figure might be drawn
  reqPa <- c("utils","wrMisc")
  chPa <- sapply(reqPa, requireNamespace, quietly=TRUE)
  if(any(!chPa)) stop("package(s) '",paste(reqPa[which(!chPa)], collapse="','"),"' not found ! Please install first from CRAN")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  if(debug) {message(fxNa," rpf1")}

  .adjustDecUnit <- function(txt, abbr=NULL, unit=NULL, silent=TRUE, callFrom=NULL) {
    ## (move to wrMisc?)  Replace differing pairs of units (appearing at different spots of txt)
    ## note: won't work when spanning 3 types of decimal units;  won't translate use of comma !!
    ## note: 'micro' written as 'u'
    ##  .adjustDecUnit(c("5kg","500g","2kg","100g"), silent=FALSE)
    fxNa <- wrMisc::.composeCallName(callFrom, newNa="adjustDecUnit")
    if(length(abbr) <1) abbr <- c("k","","m","u","n","p","f","a","z")
    if(length(unit) <1) unit <- c("mol","Mol","l","liter","g","V","A","s","K")
    chX <- matrix(rep(1:length(abbr), each=2)[-1*c(1,length(abbr)*2)], nrow=2)
    chPa <- sapply(unit, function(y) {apply(chX, 2,  function(x) { length(c(grep(paste0("[[:digit:]]",abbr[x[1]],y), txt), grep(paste0("[[:digit:]]",abbr[x[2]],y), txt)))==length(txt)})})
    if(any(chPa)) {
      whUnit <- which.max(colSums(chPa))
      sear <- chX[1, which.max(chPa[,whUnit])]    # index of abbrev
      repl <- paste0("000",abbr[sear +1], names(whUnit))
       #message("xx2"); xx2 <- list(txt=txt,abbr=abbr,chX=chX,chPa=chPa,unit=unit,whUnit=whUnit,sear=sear,repl=repl)
      if(!silent) message("replace '",paste0(abbr[sear], names(whUnit)),"'  by  '", repl,"'")
      sear <- paste0(abbr[sear], names(whUnit))
      txt <- sub(sear, repl, txt)
    }
  txt }


  ## check if path & file exist
  if(!grepl("\\.csv$|\\.xlsx$|\\.tsv$|\\.csv\\.gz$|\\.tsv\\.gz$", fileName)) message(fxNa,"Trouble ahead ? Expecting .xlsx, .csv or .tsv file (the file'",fileName,"' might not be right format) !!")
  paFi <- wrMisc::checkFilePath(fileName, path, expectExt=NULL, compressedOption=TRUE, stopIfNothing=TRUE, callFrom=fxNa, silent=silent,debug=debug)
  if(debug) {message(fxNa," rpf2")}

  ## read file
  out <- counts <- fileTy <- NULL                    # initialize default
  if(length(grep("\\.xlsx$", paFi)) >0) {
    ## Extract out of Excel
    reqPa <- c("readxl")
    chPa <- sapply(reqPa, requireNamespace, quietly=TRUE)
    if(any(!chPa)) stop("package( '",paste(reqPa[which(!chPa)], collapse="','"),"' not found ! Please install first from CRAN")
    sheets <- if(debug) try(readxl::excel_sheets(paFi[1]), silent=TRUE) else suppressMessages(try(readxl::excel_sheets(paFi[1]), silent=TRUE))
    if(debug) {message(fxNa," rpf3")}
    if(inherits(sheets, "try-error")) { message(fxNa,"Unable to read file '",paFi[1],"' ! Returning NULL; check format & rights to read")
    } else {
      if(!silent) message(fxNa,"Found sheets ",wrMisc::pasteC(sheets,quoteC="'"))
      prSh <- which(sheets == "Protein sets")
      if(length(prSh) <1) prSh <- grep("Protein", sheets)
      if(length(prSh) >1) {
        if(!silent) message(fxNa,"Multipe sheets containing 'Protein' found, using 1st :",sheets[prSh])
        prSh <- prSh[1]
      } else if(length(prSh) <1) {
        prSh <- length(sheets)
        if(!silent) message(fxNa,"No sheets containing 'Protein' in name found, try using : '",sheets[prSh],"'")
      }
      if(debug) message(fxNa,"Ready to use sheet ",prSh,"   rpf4")
      out <- as.data.frame(if(debug) readxl::read_xlsx(paFi[1], sheet=prSh) else suppressMessages(readxl::read_xlsx(paFi[1], sheet=prSh)))
      qCoSh <- which(sheets=="Quant config")
      quantConf <- if(length(qCoSh) >0) { if(debug) readxl::read_xlsx(paFi[1], sheet=qCoSh[1]) else suppressMessages(readxl::read_xlsx(paFi[1], sheet=qCoSh[1]))} else NULL
      if(length(quantConf) >0) {
        tmp <- quantConf[[1]]
        quantConf <- quantConf[[2]]
        names(quantConf) <- sub("^#","X.",tmp)
      }
      if(debug) message(fxNa,"Initial xlsx read as ",nrow(tmp)," lines & ",ncol(tmp)," rows   rpf5")
      fileTy <- "xlsx"
      ## tibble colnames may include/start with '#'
      ## adopt annotCol
      annotCol <- sub("^X.","#",annotCol)                                      # adopt to real reading of colnames
      bestT <- 0 }                                                               # initialize
  } else {
    ## anythng else but Excel xlsx ...
    ## try to find out which input format; read according to extension multiple times and choose the one with most columns
    tmp <- list()
    if(grepl("\\.txt$|\\.txt\\.gz$",paFi[1])) tmp[[1]] <- try(utils::read.delim(paFi[1], stringsAsFactors=FALSE), silent=TRUE)             # read tabulated text-file (from Proline/parse_Proline.R)
    if(grepl("\\.csv$|\\.csv\\.gz$",paFi[1])) {
      tmp[[2]] <- try(utils::read.csv(file=paFi[1], stringsAsFactors=FALSE), silent=TRUE)               # read US csv-file
      tmp[[3]] <- try(utils::read.csv2(file=paFi[1], stringsAsFactors=FALSE), silent=TRUE)}             # read Euro csv-file
    if(grepl("\\.tsv$|\\.tsv\\.gz$",paFi[1])) {
      tmp[[4]] <- try(utils::read.csv(file=paFi[1], stringsAsFactors=FALSE, sep='\t', header=TRUE))    # read US comma tsv-file
      tmp[[5]] <- try(utils::read.csv2(file=paFi[1], stringsAsFactors=FALSE, sep='\t', header=TRUE))}  # read Euro tsv-file

    if(debug) {message(fxNa," rpf6"); rpf6 <- list(fileName=fileName,paFi=paFi,chPa=chPa,path=path,tmp=tmp)}
    chCl <- sapply(tmp, inherits, "try-error") | sapply(tmp, length) <2
    if(all(chCl)) stop(" Failed to extract data from '",fileName,"'  (check format & rights to read)")
    nCol <- sapply(tmp, function(x) if(length(x) >0) {if(! inherits(x, "try-error")) ncol(x) else NA} else NA)
    bestT <- which.max(nCol)
    if(length(bestT) <1) stop("Problem when reading flat file : insufficient columns, check type of input !")
    if(debug) message(fxNa,"Reading flat file, best as type no ",bestT," , ie as ",c("txt","US-csv","Euro-csv","tsv (US)","tsv (Euro)")[bestT],"  rpf7")

    out <- tmp[[bestT]]
    fileTy <- c("txt","US.csv","Euro.csv")[bestT]
    if(bestT==1) {
      if(length(annotCol) >1) { if(debug) message(fxNa,"Set 1st letter of annotCol[1:2] as upercase")
        substr(annotCol[1:2], 1, 1) <- toupper(substr(annotCol[1:2], 1, 1))}    # capitalize 1st letter to get  c("Accession","Description",...
      chQuCol <- if(length(quantCol)==1) grep(quantCol, colnames(out)) else match(quantCol, colnames(out))
      if(length(chQuCol) <1 & length(quantCol)==1) {  quantCol <- "^Intensity"                      # adjust
        chQuCol <- grep(quantCol, colnames(out))
        if(!silent) message(fxNa,"Data read as txt, adjusting argument 'quantCol' to '^Intensity' ",c("NOT ")[length(chQuCol) <1]," succeful")
      }
    }
    if(debug) {message(fxNa," rpf8"); rpf8 <- list(fileName=fileName,chPa=chPa,path=path,tmp=tmp,out=out,bestT=bestT)}

    ## quant info settings from separate sheet:
    quanConfFile <- "Quant config.tsv"
    if(length(fileName) >1) if(all(!is.na(fileName[2]) & nchar(fileName[2]) >0)) quanConfFile <- fileName[2]
    paFi2 <- if(all(chPa)) {  # extract & use path out of fileName of csv
      if(is.null(grep(dirname(fileName[1]),quanConfFile))) file.path(path[1],quanConfFile) else file.path(dirname(fileName),quanConfFile) # use path of fileName -if present, otherwise add path
    } else file.path(path[1], quanConfFile)
    if(file.exists(paFi2)) {
      quantConf <- try(utils::read.csv(file=paFi2, stringsAsFactors=FALSE, sep='\t', header=TRUE))
      if(inherits(quantConf, "try-error")) quantConf <- NULL else {
        tmp <- colnames(quantConf)[-1]
        quantConf <- as.character(quantConf)[-1]
        names(quantConf) <- tmp
      }
    } else {quantConf <- NULL
      if(debug) message(fxNa," file '",paFi2,"' not found for quantifaction specific information   rpf8c")}
    ## adopt annotCol
    annotCol <- sub("^#","X.",annotCol)
  }
  if(debug) {message(fxNa," rpf9")}

  if(length(tmp) >0) {
    ## look for  annotation columns
    if(!isTRUE(separateAnnot)) separateAnnot <- FALSE
    if(any(c(length(quantCol), length(annotCol)) <1)) separateAnnot <- FALSE else {
      if(any(all(is.na(quantCol)), all(is.na(annotCol)))) separateAnnot <- FALSE
    }
    if(debug) {message(fxNa," rpf10")}

    ## extract meta-data and abundances from main table
    metaCo <- wrMisc::naOmit(if(length(annotCol) >1) wrMisc::extrColsDeX(out, extrCol=annotCol, doExtractCols=FALSE, silent=silent, callFrom=fxNa) else grep(annotCol,colnames(out)))
    chMissCo <- !(annotCol %in% colnames(out)[metaCo])
    if(any(chMissCo)) {
      if(any(chMissCo[1:2])) {
        metaCo2 <- match(tolower(annotCol[1:2]), tolower(colnames(out)))
        if(length(metaCo2) <1) stop("Can't find column ",wrMisc::pasteC(annotCol[which(chMissCo[1:2])])," in data")
        if(!silent) message(fxNa,"Found only partial match of 'annotCol' (upper/lower-case issues), adjusting names ")
        colnames(out)[metaCo2] <- annotCol[1:2]
        metaCo <- union(metaCo2,metaCo) }
    }
    if(debug) {message(fxNa," rpf11")}

    ## further refine : separate Accession (ie annotCol[1]) (eg P00359) from ProteinName (eg TDH3)
    annot <- as.matrix(out[,metaCo])
    chSep <- nchar(utils::head(annot[,annotCol[1]])) - nchar(gsub("\\|","",utils::head(annot[,annotCol[1]])))
    tmp <- if(any(chSep >1)) sub("^[[:alpha:]]+\\||^[[:alpha:]]+_[[:alpha:]]+\\|","", annot[,annotCol[1]])        # presume presence of database origin-tag -> remove  (eg sp|...)
    tmp3 <- sub("^[[:alpha:]]+_{0,1}[[:alpha:]]*\\|[[:upper:]][[:digit:]]+-{0,1}[[:digit:]]*\\|[[:alnum:]]+_{0,1}[[:alnum:]]*_{0,1}[[:alnum:]]*\\ ",
      "",annot[,annotCol[2]])        # supposes annotCol[2], without db|accession|EntryName
    tmp2 <- cbind(Accession=sub("\\|[[:print:]]+$","",tmp), EntryName=sub("^[[:alnum:]]+-{0,1}[[:digit:]]{0,2}\\|","",tmp),
      ProteinName=sub("\\ [[:upper:]]{2}=[[:print:]]+","",tmp3) , GN=NA, Species=NA, Contam=NA, SpecType=NA) #
    if(debug) {message(fxNa," rpf12")}

    ## recover GN
    GNLi <- grep("\\ GN=[[:upper:]]{2,}[[:digit:]]*", tmp3)
      if(length(GNLi) >0) { zz <- sub("[[:print:]]+\\ GN=", "",tmp3[GNLi])             # remove surplus to left
        tmp2[GNLi,"GN"] <- sub("\\ +$","",sub("\\ [[:print:]]+","",zz)) }              # remove surplus to right (and right trailing space)
    ## recover OS
    OSLi <- grep("\\ OS=[[:upper:]][[:lower:]]+", tmp3)
      if(length(OSLi) >0) { zz <- sub("[[:print:]]+\\ OS=", "",tmp3[OSLi])             # remove surplus to left
        tmp2[OSLi,"Species"] <- sub("\\ [[:upper:]]{2}=[[:print:]]+","",zz) }          # remove surplus to right (next tag)
    if(debug) {message(fxNa," rpf13")}

    ## extract/check on Description
    annot <- cbind(tmp2,annot[,-1*match(annotCol[1], colnames(annot))])
    chDesc <- colnames(annot) %in% "description"
    if(any(chDesc)) colnames(annot)[which(chDesc)] <- "Description"                # make uniform to other
    chUniNa <- duplicated(annot[,1])
    if(any(chUniNa) & !silent) message(fxNa,"Caution: ",sum(chUniNa)," Accession entries appear repeatedly !  Making unique rownames by adding counter ...")
    rowNa <- if(any(chUniNa)) wrMisc::correctToUnique(annot[,1], callFrom=fxNa) else annot[,1]
    rownames(annot) <- rowNa

    if(isTRUE(remStrainNo)) {
      chSp <- grep("^[[:upper:]][[:lower:]]*\\ [[:lower:]]+\\ [[:print:]]*", annot[,"Species"])
      if(length(chSp) >0) { nch1 <- nchar(sub("^[[:upper:]][[:lower:]]+\\ [[:lower:]]+", "", annot[chSp,"Species"]))
        annot[chSp,"Species"] <- substr(annot[chSp,"Species"], 1, nchar(annot[chSp,"Species"]) - nch1) }
    }
    if(debug) {message(fxNa," rpf14"); rpf14 <- list(out=out,annot=annot,specPref=specPref,quantCol=quantCol,rowNa=rowNa,tmp=tmp)}

    ## set protein/gene annotation to common format
    chNa <- match(c("GN","ProteinName"), colnames(annot))
    colnames(annot)[chNa] <- c("GeneName","Description")

    ## locate special groups, column "SpecType"
    if(length(specPref) >0) {
      annot <- .extrSpecPref(specPref, annot, useColumn=c("Species","EntryName","GeneName","Accession"), silent=silent, debug=debug, callFrom=fxNa) }

    ## locate & extract abundance/quantitation data
    if(length(quantCol) >1) { abund <- as.matrix(wrMisc::extrColsDeX(out, extrCol=quantCol, doExtractCols=TRUE, silent=silent, callFrom=fxNa))
    } else {
      quantCol2 <- grep(quantCol, colnames(out))
      chNa <- is.na(quantCol2)
      if(all(chNa)) stop("Could not find any of of the columns specified in argument 'quantCol' !")
      if(any(chNa)) {
        if(!silent) message(fxNa,"Could not find columns ",wrMisc::pasteC(quantCol2[which(chNa)],quote="'")," .. omit")
        quantCol2 <- wrMisc::naOmit(quantCol2)}
      abund <- as.matrix(out[,quantCol2]) }                    # abundance val
    ## peptide counts
    if(length(pepCountCol) >1) { supCol <- lapply(pepCountCol, grep, colnames(out))
      chLe <- sapply(supCol,length)
      chLe <- chLe==ncol(abund)
      if(any(chLe)) { pepCount <- array(dim=c(nrow(out), ncol(abund), sum(chLe)),
        dimnames=list(rowNa, colnames(abund), sub("\\^peptides_count","NoOfPeptides",sub("\\^psm_count","PSM",sub("_$","",pepCountCol)))[which(chLe)]))
        for(i in 1:sum(chLe)) pepCount[,,i] <- as.matrix(out[,supCol[[which(chLe)[i]]]])
      } else pepCount <- NULL
    }
    if(debug) {message(fxNa," rpf15"); rpf15 <- list(out=out,annot=annot,specPref=specPref,quantCol=quantCol,rowNa=rowNa,tmp=tmp)}

    ## check abundance/quantitation data
    chNum <- is.numeric(abund)
    if(!chNum) {abund <- apply(out[,quantCol2], 2, wrMisc::convToNum, convert="allChar", silent=silent, callFrom=fxNa)}
    rownames(abund) <- rowNa
    ##
    ## Custom (alternative) colnames
    if(length(sampleNames) ==ncol(abund)) {
      colnames(abund) <- colnames(pepCount) <- sampleNames
    } else {
      colnames(abund) <- sub(if(length(quantCol) >1) "^abundance" else quantCol, "", colnames(abund))
      if(trimColnames) colnames(abund) <- wrMisc::.trimFromStart(wrMisc::.trimFromEnd(colnames(abund)))
    }

    ## treat colnames, eg problem with pxd001819 : some colnames writen as amol, other as fmol
    adjDecUnits <- FALSE
    if(adjDecUnits) {
      colnames(abund) <- sub("^ +","", .adjustDecUnit(wrMisc::trimRedundText(colnames(abund), silent=silent, debug=debug, callFrom=fxNa)))
    }
    if(debug) { message(fxNa," rpf16"); rpf16 <- list(abund=abund,annot=annot,sampleNames=sampleNames,normalizeMeth=normalizeMeth,refLi=refLi,sdrf=sdrf,suplAnnotFile=suplAnnotFile,path=path,fxNa=fxNa,silent=silent)}

    ## check for reference for normalization
    refLiIni <- refLi
    if(is.character(refLi) && length(refLi) ==1) {
      refLi <- which(annot[,"SpecType"] ==refLi)
      if(length(refLi) <1 && identical(refLiIni, "mainSpe")) refLi <- which(annot[,"SpecType"] =="mainSpecies")              # fix compatibility problem  'mainSpe' to 'mainSpecies'
      if(length(refLi) <1 ) { refLi <- 1:nrow(abund)
        if(!silent) message(fxNa,"Could not find any proteins matching argument 'refLi=",refLiIni,"', ignoring ...")
      } else {
        if(!silent) message(fxNa,"Normalize using (custom) subset of ",length(refLi)," lines specified as '",refLiIni,"'")}}    # may be "mainSpe"
    if(debug) { message(fxNa," rpf16b"); rpf16b <- list()}

    ## take log2 & normalize
    if(length(normalizeMeth) <1) normalizeMeth <- "median"
    quant <- try(wrMisc::normalizeThis(log2(abund), method=normalizeMeth, mode="additive", refLines=refLi, silent=silent, debug=debug, callFrom=fxNa), silent=!debug)      
    if(debug) { message(fxNa,"rpf16c .. dim quant: ", nrow(quant)," li and  ",ncol(quant)," cols; colnames : ",wrMisc::pasteC(colnames(quant))," ")}

    ### GROUPING OF REPLICATES AND SAMPLE META-DATA
    if(length(suplAnnotFile) >0 || length(sdrf) >0) {
      if(isTRUE(suplAnnotFile[1])) suplAnnotFile <- paFi
      setupSd <- readSampleMetaData(sdrf=sdrf, suplAnnotFile=suplAnnotFile, quantMeth="PL", path=path, abund=utils::head(quant), groupPref=groupPref, silent=silent, debug=debug, callFrom=fxNa)
    }
    if(debug) {message(fxNa,"rpf16d .."); rpf16d <- list(sdrf=sdrf,gr=gr,suplAnnotFile=suplAnnotFile, quant=quant,refLi=refLi,annot=annot,setupSd=setupSd,sampleNames=sampleNames)}

    ## finish groups of replicates & annotation setupSd
    ## need to transform fmol in amol for proper understanding & matching with metadata
    setupSd <- .checkSetupGroups(abund=abund, setupSd=setupSd, gr=gr, sampleNames=sampleNames, quantMeth="PL", silent=silent, debug=debug, callFrom=fxNa)
    colNa <- if(length(setupSd$sampleNames)==ncol(abund)) setupSd$sampleNames else setupSd$groups
    chGr <- grepl("^X[[:digit:]]", colNa)                                                # check & remove heading 'X' from initial column-names starting with digits
    if(any(chGr)) colNa[which(chGr)] <- sub("^X","", colNa[which(chGr)])                 # 
    colnames(quant) <- colnames(abund) <- colNa
    if(length(setupSd$sampleNames)==ncol(abund)) setupSd$sampleNames <- colNa else setupSd$groups <- colNa
    if(length(dim(counts)) >1 && length(counts) >0) colnames(counts) <- colNa

    if(debug) {message(fxNa,"Read sample-meta data, rpf17"); rpf17 <- list(sdrf=sdrf,suplAnnotFile=suplAnnotFile,abund=abund, quant=quant,refLi=refLi,annot=annot,setupSd=setupSd,sampleNames=sampleNames)}

    ## main plotting of distribution of intensities
    custLay <- NULL
    if(is.numeric(plotGraph) && length(plotGraph) >0) {custLay <- as.integer(plotGraph); plotGraph <- TRUE} else {
        if(!isTRUE(plotGraph)) plotGraph <- FALSE}
    if(plotGraph) .plotQuantDistr(abund=abund, quant=quant, custLay=custLay, normalizeMeth=normalizeMeth, softNa="Proline",
      refLi=refLi, refLiIni=refLiIni, tit=titGraph, silent=silent, callFrom=fxNa, debug=debug)


    ## meta-data to export
    notes <- c(inpFile=paFi, qmethod="Proline", normalizeMeth="none", call=match.call(),
      created=as.character(Sys.time()), wrProteo.version=utils::packageVersion("wrProteo"), machine=Sys.info()["nodename"])
    ##
    if(separateAnnot) {
      if(!is.numeric(quant) && logConvert) { message(fxNa,"Problem: Abundance data seem not numeric, can't transform log2 !")}
      out <- list(raw=abund, quant=quant, annot=annot, counts=pepCount, sampleSetup=setupSd, quantNotes=quantConf, notes=notes)
      if(!logConvert) out$quant <- 2^out$quant }
  }
  ## final Proline-import
  out }
  
