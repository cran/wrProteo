#' Read csv or txt files exported from MS-Angel and Proline
#'
#' Quantification results form MS-Angel and Proline \href{http://www.profiproteomics.fr/proline/}{Proline} exported as xlsx format can be read directly.
#' Besides, files in csv (European and US format), tsv or tabulated txt can be read, too.
#' Then relevant information gets extracted, the data can optionally normalized and displayed as boxplot or vioplot. 
#' The final output is a list containing 6 elements: \code{$raw}, \code{$quant},  \code{$annot}, \code{$counts}, \code{$quantNotes} and \code{$notes}. 
#' Alternatively, a data.frame with annotation and quantitation data may be returned if \code{separateAnnot=FALSE}.
#' Note: There is no normalization by default since quite frequently data produced by Proline are already sufficiently normalized. 
#' The figure produced using the argument \code{plotGraph=TRUE} may help judging if the data appear sufficiently normalized (distribtions should align).
#' 
#' @details
#' This function has been developed using Proline version 1.6.1 coupled with MS-Angel 1.6.1. 
#' The format of the exported file depends on the columns chosen for export, default settings from Proline and MS-Angel work fine. 
#' 
#' @param fileName (character) name of file to read
#' @param path (character) optional path (note: Windows backslash sould be protected or written as '/')
#' @param normalizeMeth (character) normalization method (will be sent to  \code{\link[wrMisc]{normalizeThis}}) 
#' @param logConvert (logical) convert numeric data as log2, will be placed in $quant
#' @param sampleNames (character) new column-names for quantification data (ProteomeDiscoverer does not automatically use file-names from spectra); Please use with care since order of samples might be different as you expect
#' @param quantCol (character or integer) colums with main quantitation-data : precise colnames to extract, or if length=1 content of \code{quantCol} will be used as pattern to search among column-names for $quant using \code{grep} 
#' @param annotCol (character) precise colnames or if length=1 pattern to search among column-names for $annot
#' @param pepCountCol (character) pattern to search among column-names for count data of PSM and NoOfPeptides
#' @param trimColnames (logical) optional trimming of column-names of any redundant characters from beginning and end 
#' @param refLi (integer) custom decide which line of data is main species, if single character entry it will be used to choose a group of species (eg 'mainSpe')
#' @param separateAnnot (logical) separate annotation form numeric data (quantCol and annotCol must be defined)
#' @param plotGraph (logical or matrix of integer) optional plot vioplot of initial data; if integer, it will be passed to \code{layout} when plotting
#' @param tit (character) custom title to plot
#' @param wex (integer) relative expansion factor of the violin-plot (will be passed to \code{\link[wrGraph]{vioplotW}})
#' @param graphTit (character) (depreciated custom title to plot), please use 'tit'
#' @param specPref (character or list) define characteristic text for recognizing (main) groups of species (1st for comtaminants - will be marked as 'conta', 2nd for main species- marked as 'mainSpe', 
#'  and optional following ones for supplemental tags/species - maked as 'species2','species3',...); 
#'  if list and list-element has multiple values they will be used for exact matching of accessions (ie 2nd of argument \code{annotCol})
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @return list with \code{$raw} (initial/raw abundance values), \code{$quant} with final normalized quantitations, \code{$annot} (columns ), \code{$counts} an array with 'PSM' and 'NoOfPeptides', \code{$quantNotes} and \code{$notes}; or a data.frame with quantitation and annotation if \code{separateAnnot=FALSE}
#' @seealso \code{\link[utils]{read.table}} 
#' @examples
#' path1 <- system.file("extdata", package="wrProteo")
#' fiNa <- "exampleProlineABC.csv"
#' dataABC <- readProlineFile(file.path(path1,fiNa))
#' summary(dataABC$quant)
#' matrixNAinspect(dataABC$quant, gr=as.factor(substr(colnames(dataABC$quant),1,1))) 
#' @export
readProlineFile <- function(fileName, path=NULL, normalizeMeth=NULL, logConvert=TRUE, sampleNames=NULL, quantCol="^abundance_", 
  annotCol=c("accession","description","is_validated","protein_set_score","X.peptides","X.specific_peptides"), pepCountCol=c("^psm_count_","^peptides_count_"), 
  trimColnames=FALSE, refLi=NULL,separateAnnot=TRUE,plotGraph=TRUE,tit=NULL,graphTit=NULL,wex=2, specPref=c(conta="_conta\\|", mainSpecies="OS=Homo sapiens"), silent=FALSE,callFrom=NULL){
  ## 'quantCol', 'annotCol' (character) exact col-names or if length=1 pattern to search among col-names for $quant or $annot
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readProlineFile")
  opar <- graphics::par(no.readonly=TRUE)
  ## check & read file
  chPa <- length(grep("/",fileName)) >0 | length(grep("\\\\",fileName)) >0      # check for path already in fileName "
  if(length(path) <1) path <- "."
  paFi <- if(!chPa) file.path(path[1],fileName[1]) else fileName[1]             # use path only when no path joined to fileName
  chFi <- file.exists(paFi)
  if(!chFi) stop(" file ",fileName," was NOT found ",if(length(path) >0) paste(" in path ",path)," !")
  counts <- NULL
  if(length(grep("\\.xlsx$",fileName)) >0) {
    ## Extract out of Excel
    chPa <- try(find.package("readxl"), silent=TRUE)
    if("try-error" %in% class(chPa)) stop(fxNa,"package 'readxl' not found ! Please install first") 
    sheets <- readxl::excel_sheets(paFi)
    if(!silent) message(fxNa,"Found sheets ",wrMisc::pasteC(sheets,quoteC="'")) 
    prSh <- which(sheets == "Protein sets")
    if(length(prSh) <1) prSh <- grep("Protein", sheets)
    if(length(prSh) >1) { if(!silent) message(fxNa,"No sheet named 'Protein sets' found, multiple sheets contain term 'Sheet' :",paste(sheets[prSh])," using first")
      prSh <- prSh[1] }
    if(length(prSh) <1) {prSh <- length(sheets)
      if(!silent) message(fxNa,"No sheets with 'Protein' in name found, using : '",paste(sheets[prSh]),"'")
    }
    out <- readxl::read_xlsx(paFi, sheet=prSh)
    qCoSh <- which(sheets=="Quant config")
    quantConf <- if(length(qCoSh) >0) readxl::read_xlsx(paFi, sheet=qCoSh[1]) else NULL
    if(length(quantConf) >0) {
      tmp <- quantConf[[1]]
      quantConf <- quantConf[[2]]
      names(quantConf) <- sub("^#","X.",tmp)
    }
    ## adopt annotCol
    annotCol <- sub("^X.","#",annotCol)                                         # adopt to real reading of colnames
  } else {
    ## try to find out which input format; read according to extension multiple times and choose the one with most columns 
    tmp <- list()
    if(length(grep("\\.txt$",fileName)) >0) tmp[[1]] <- try(utils::read.delim(paFi, stringsAsFactors=FALSE), silent=TRUE)             # read tabulated text-file
    if(length(grep("\\.csv$",fileName)) >0) tmp[[2]] <- try(utils::read.csv(paFi, stringsAsFactors=FALSE), silent=TRUE)               # read US csv-file
    if(length(grep("\\.csv$",fileName)) >0) tmp[[3]] <- try(utils::read.csv2(paFi, stringsAsFactors=FALSE), silent=TRUE)              # read Euro csv-file
    if(length(grep("\\.tsv$",fileName)) >0) tmp[[4]] <- try(utils::read.csv(file=paFi, stringsAsFactors=FALSE, sep='\t', header=TRUE)) # read US comma tsv-file
    chCl <- sapply(tmp,class) =="try-error" 
    if(length(chCl) <1) stop("Failed to recognize file extensions of input data (unknown format)")
    if(any(chCl)) {if(all(chCl)) stop(" Failed to extract data (unknown format) from ",fileName)}
    nCol <- sapply(tmp, function(x) if(length(x) >0) {if(class(x) != "try-error") ncol(x) else NA} else NA)
    bestT <- which.max(nCol)
    out <- tmp[[bestT]]
    ## tibble colnames may include/start with '#' ... adopt to rest
    corColNa <- grep("^#",colnames(out))
    if(length(corColNa) >0) colnames(out)[which(corColNa)] <- sub("^#","X.",colnames(out)[which(corColNa)])    # make colnames alike

    ## quant info settings from separate sheet:
    quanConfFile <- if(length(fileName) >1 & !is.na(fileName[2])) fileName[2] else "Quant config.tsv"
    paFi <- if(chPa) {  # extract & use path out of fileName
      if(is.null(grep(dirname(fileName[1]),quanConfFile))) file.path(path[1],quanConfFile) else file.path(dirname(fileName),quanConfFile) # use path of fileName -if present, otherwise add path    
    } else file.path(path[1],quanConfFile)
    if(file.exists(paFi)) {
      quantConf <- try(utils::read.csv(file=paFi, stringsAsFactors=FALSE, sep='\t', header=TRUE))
      if("try-error" %in% class(quantConf)) quantConf <- NULL else {
        tmp <- colnames(quantConf)[-1]
        quantConf <- as.character(quantConf)[-1]
        names(quantConf) <- tmp 
      }
    } else quantConf <- NULL
    ## adopt annotCol
    annotCol <- sub("^#","X.",annotCol) 
  }  
  ## look for separate annotation columns ?
  if(any(c(length(quantCol), length(annotCol)) <1)) separateAnnot <- FALSE else {
    if(any(all(is.na(quantCol)), all(is.na(annotCol)))) separateAnnot <- FALSE
  }  
  ## extract meta-data and abundances from main table
  metaCo <- wrMisc::naOmit(if(length(annotCol) >1) wrMisc::extrColsDeX(out, extrCol=annotCol, doExtractCols=FALSE, callFrom=fxNa) else grep(annotCol,colnames(out)))
  ## further refine : separate Accesion (ie annotCol[1]) (eg P00359) from ProteinName (eg TDH3)
  annot <- as.matrix(out[,metaCo])
  chSep <- nchar(utils::head(annot[,annotCol[1]])) - nchar(gsub("\\|","",utils::head(annot[,annotCol[1]])))
  tmp <- if(any(chSep >1)) sub("^[[:alpha:]]+\\||^[[:alpha:]]+_[[:alpha:]]+\\|","", annot[,annotCol[1]])        # presume presence of database origin-tag -> remove  (eg sp|...)
  tmp3 <- sub("^[[:alpha:]]+\\|[[:alnum:]]+\\|{0,1}[[:alnum:]]*\\_{0,1}[[:upper:]]*\\ ","",annot[,annotCol[2]]) # without db|accession|EntryName
  tmp2 <- cbind(Accession=sub("\\|[[:print:]]+$","",tmp), EntryName=sub("^[[:alnum:]]+\\|","",tmp), 
    ProteinName=sub("\\ [[:upper:]]{2}=[[:print:]]+","",tmp3) , GN=NA, Species=NA, Contam=NA, SpecType=NA) #
  ## recover GN
  GNLi <- grep("\\ GN=[[:upper:]]{2,}[[:digit:]]", tmp3)
    if(length(GNLi) >0) { zz <- sub("[[:print:]]+\\ GN=", "",tmp3[GNLi])             # remove surplus to left
      tmp2[GNLi,"GN"] <- sub("\\ +$","",sub("\\ [[:print:]]+","",zz)) }              # remove surplus to right (and right trailing space)
  ## recover OS
  OSLi <- grep("\\ OS=[[:upper:]][[:lower:]]+", tmp3)
    if(length(OSLi) >0) { zz <- sub("[[:print:]]+\\ OS=", "",tmp3[OSLi])             # remove surplus to left
      tmp2[OSLi,"Species"] <- sub("\\ [[:upper:]]{2}=[[:print:]]+","",zz) }          # remove surplus to right  
  ## locate special groups, column "SpecType"
  if(length(specPref) >0) for(i in 1:length(specPref)) {         # locate specPref
    naSp <- if(i==1) "conta" else {if(i==2) "mainSpe" else paste0("species",i-1)}
    chSp <- unlist(specPref[i]) 
    chSp <- if(length(chSp) >1) match(chSp,annot[,annotCol[1]]) else grep(unlist(specPref[i]), annot[,annotCol[2]])    # line-no to set tag a 'SpecType'
    if(length(chSp) >0) {
      chNa <- is.na(tmp2[chSp,"SpecType"])
      if(any(!chNa) & !silent) message(fxNa," Beware, ",sum(!chNa)," 'SpecType' will be overwritten by ",naSp)
      tmp2[chSp,"SpecType"] <- naSp }
  }  
  
  annot <- cbind(tmp2,annot[,-1*match(annotCol[1], colnames(annot))])
  chDesc <- colnames(annot) %in% "description" 
  if(any(chDesc)) colnames(annot)[which(chDesc)] <- "Description"                # make uniform to other
  chUniNa <- duplicated(annot[,1])
  if(any(chUniNa) & !silent) message(fxNa," Caution: ",sum(chUniNa)," Accession entries appear repeatedly !  Making unique rownames by adding counter ...") 
  rowNa <- if(any(chUniNa)) wrMisc::correctToUnique(annot[,1],callFrom=fxNa) else annot[,1] 
  rownames(annot) <- rowNa
    
  ## locate & extract abundance/quantitation data
  if(length(quantCol) >1) { abund <- as.matrix(wrMisc::extrColsDeX(out, extrCol=quantCol, doExtractCols=TRUE, callFrom=fxNa))
  } else {
    quantCol2 <- grep(quantCol, colnames(out)) 
    chNa <- is.na(quantCol2)
    if(all(chNa)) stop("Could not find any of of the columns specified in argument 'quantCol' !")
    if(any(chNa)) { 
      if(!silent) message(fxNa," Could not find columns ",wrMisc::pasteC(quantCol2[which(chNa)],quote="'")," .. omit")
      quantCol2 <- wrMisc::naOmit(quantCol2)} 
    abund <- as.matrix(out[,quantCol2]) }                    # abundance val
  ## peptide counts
  if(length(pepCountCol) >1) { supCol <- lapply(pepCountCol, grep, colnames(out))
    chLe <- sapply(supCol,length) ==ncol(abund)
    if(any(chLe)) {pepCount <- array(dim=c(nrow(out), ncol(abund), sum(chLe)), dimnames=list(rowNa,colnames(abund),c("PSM","NoOfPeptides")[which(chLe)]))
      for(i in which(chLe)) pepCount[,,i] <- as.matrix(out[,supCol[[i]]]) }                                                 
  }    
  ## check abundance/quantitation data
  chNum <- is.numeric(abund)
  if(!chNum) {abund <- apply(out[,quantCol2], 2, wrMisc::convToNum, convert="allChar", callFrom=fxNa)} 
  rownames(abund) <- rowNa
  ##
  ## Custom (alternative) colnames
  if(length(sampleNames) ==ncol(abund)) {   
    colnames(abund) <- colnames(pepCount) <- sampleNames
  } else {
    colnames(abund) <- sub(if(length(quantCol) >1) "^abundance" else quantCol, "", colnames(abund))
    if(trimColnames) colnames(abund) <- wrMisc::.trimFromStart(wrMisc::.trimFromEnd(colnames(abund)))
    }
  
  ## normalize (if chosen)
  if(length(normalizeMeth)==1){ 
    ## check for reference for normalization
    refLiIni <- refLi
    if(is.character(refLi) & length(refLi)==1) refLi <- which(annot[,"SpecType"]==refLi)
    if(length(refLi) <1) { normalizeMeth <- NULL
      message(fxNa," could not find any protein matching argument 'refLi', ignoring ...") } else {
      if(!silent) message(fxNa," normalize using subset of ",length(refLi)) } }   # may be "mainSpe"
  if(length(normalizeMeth)==1) { rawD <- abund
    abund <- wrMisc::normalizeThis(log2(abund), method=normalizeMeth, refLines=refLi, callFrom=fxNa)
  } else {rawD <- abund; abund <- log2(abund)}
  ##  
  ## main plotting of distribution of intensities
  custLay <- NULL
  if(length(plotGraph) >0) {if(is.numeric(plotGraph)) {custLay <- plotGraph; plotGraph <- TRUE
    } else { plotGraph <- as.logical(plotGraph[1])}}
  if(plotGraph){
    if(length(custLay) >0) graphics::layout(custLay) else if(length(normalizeMeth) >0) graphics::layout(1:2) 
    graphics::par(mar=c(3, 3, 3, 1))                          # mar: bot,le,top,ri
    if(length(graphTit) >0) message(fxNa,"argument 'graphTit' is depreciated, please rather use 'tit' ")
    if(is.null(tit) & !is.null(graphTit)) tit <- graphTit
    if(is.null(tit)) tit <- "Distribution of quantification values"
    chGr <- try(find.package("wrGraph"), silent=TRUE)
    chSm <- try(find.package("sm"), silent=TRUE)
    misPa <- c("try-error" %in% class(chGr),"try-error" %in% class(chSm))
    if(any(misPa)) if(!silent) message(fxNa," missing package ",wrMisc::pasteC(c("wrGraph","sm")[which(misPa)],quoteC="'")," for drawing vioplots")
    if(length(normalizeMeth) >0) {
      if(any(misPa)) { 
        ## wrGraph not available : simple boxplot  
        graphics::boxplot(log2(rawD), main=paste(tit," (initial)"), las=1, outline=FALSE) 
        graphics::abline(h=round(log2(stats::median(rawD, na.rm=TRUE))) +c(-1:1), lty=2, col=grDevices::grey(0.6))      
      } else {                                          # wrGraph & sm are available
        wrGraph::vioplotW(log2(rawD), tit=paste(tit," (initial)"), wex=wex) 
        graphics::abline(h=round(stats::median(rawD, na.rm=TRUE)) +(-1:1), lty=2, col=grDevices::grey(0.6))           
      }
    }
    if(length(normalizeMeth) >0) tit <- paste(tit," (normalized)")
    if(any(misPa)) { 
      ## wrGraph not available : simple boxplot  
      graphics::boxplot(abund, main=tit, las=1, outline=FALSE) 
      graphics::abline(h=round(stats::median(abund, na.rm=TRUE)) +c(-1:1), lty=2, col=grDevices::grey(0.6))      
    } else {                                          # wrGraph & sm are available
      wrGraph::vioplotW(abund, tit=paste(tit, sep=" "), wex=wex) 
      graphics::abline(h=round(stats::median(abund, na.rm=TRUE)) +(-1:1), lty=2, col=grDevices::grey(0.6))           
    }        
    on.exit(graphics::par(opar)) }                         #
  ## meta-data to export
  notes <- c(qmethod="Proline", normalizeMeth="none", call=match.call(), created=as.character(Sys.time()), 
    wrProteo.version=utils::packageVersion("wrProteo"), machine=Sys.info()["nodename"])
  ##
  if(separateAnnot) {
    if(!is.numeric(abund) & logConvert) {message(fxNa," Problem: Abundance data seem not numeric, can't transform log2 !")}
    out <- list(raw=rawD, quant=abund, annot=annot, counts=pepCount, quantNotes=quantConf, notes=notes)
    if(!logConvert) out$quant <- 2^out$quant
    }
  ## final  
  out }
    
