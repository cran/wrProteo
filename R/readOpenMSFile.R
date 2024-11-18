#' Read csv files exported by OpenMS
#'
#' Protein quantification results form \href{https://openms.de/}{OpenMS} 
#' which were exported as \code{.csv} can be imported and relevant information extracted. 
#' Peptide data get summarized by protein by top3 or sum methods.
#' The final output is a list containing the elements: \code{$annot}, \code{$raw}, \code{$quant} ie normaized final quantifications, or returns data.frame with entire content of file if \code{separateAnnot=FALSE}.
#' 
#' @details
#' This function has been developed based on the OpenMS peptide-identification and label-free-quantification module.
#' Csv input files may also be compresses as .gz.
#' 
#' Note: With this version the information about protein-modifications (PTMs) may not yet get exploited fully.
#' 
#' @param fileName (character) name of file to be read 
#' @param path (character) path of file to be read
#' @param normalizeMeth (character) normalization method (will be sent to  \code{\link[wrMisc]{normalizeThis}}) 
#' @param refLi (character or integer) custom specify which line of data is main species, if character (eg 'mainSpe'), the column 'SpecType' in $annot will be searched for exact match of the (single) term given
#' @param sampleNames (character) new column-names for quantification data (by default the names from files with spectra will be used)
#' @param quantCol (character or integer) exact col-names, or if length=1 content of \code{quantCol} will be used as pattern to search among column-names for $quant using \code{grep} 
#' @param sumMeth (character) method for summarizing peptide data (so far 'top3' and 'sum' available)
#' @param minPepNo (integer) minumun number of peptides to be used for retruning quantification
#' @param protNaCol (character) column name to be read/extracted for the annotation section (default "ProteinName")
#' @param separateAnnot (logical) if \code{TRUE} output will be organized as list with \code{$annot}, \code{$abund} for initial/raw abundance values and \code{$quant} with final normalized quantitations 
#' @param tit (character) custom title to plot
#' @param wex (integer) relative expansion factor of the violin-plot (will be passed to \code{\link[wrGraph]{vioplotW}})
#' @param specPref (character or list) define characteristic text for recognizing (main) groups of species (1st for comtaminants - will be marked as 'conta', 2nd for main species- marked as 'mainSpe', 
#'  and optional following ones for supplemental tags/species - maked as 'species2','species3',...); 
#'  if list and list-element has multiple values they will be used for exact matching of accessions (ie 2nd of argument \code{annotCol})
#' @param separateAnnot (logical) if \code{TRUE} output will be organized as list with \code{$annot}, \code{$abund} for initial/raw abundance values and \code{$quant} with final normalized quantitations
#' @param plotGraph (logical) optional plot of type vioplot of initial and normalized data (using \code{normalizeMeth}); if integer, it will be passed to \code{layout} when plotting
#' @param silent (logical) suppress messages
#' @param debug (logical) display additional messages for debugging
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @return This function returns a list with \code{$raw} (initial/raw abundance values), \code{$quant} with final normalized quantitations, \code{$annot}, \code{$counts} an array with number of peptides, \code{$quantNotes},\code{$expSetup} and \code{$notes}; or if \code{separateAnnot=FALSE} the function returns a data.frame with annotation and quantitation only 
#' @seealso \code{\link[utils]{read.table}}, \code{\link[wrMisc]{normalizeThis}}) , \code{\link{readMaxQuantFile}}, \code{\link{readProlineFile}}, \code{\link{readProtDiscovFile}} 
#' @examples
#' path1 <- system.file("extdata", package="wrProteo")
#' fiNa <- "OpenMS_tiny.csv.gz"
#' dataOM <- readOpenMSFile(file=fiNa, path=path1, tit="tiny OpenMS example")
#' summary(dataOM$quant)
#' 
#' @export
readOpenMSFile <- function(fileName=NULL, path=NULL, normalizeMeth="median", refLi=NULL, sampleNames=NULL, quantCol="Intensity", 
  sumMeth="top3", minPepNo=1, protNaCol="ProteinName", separateAnnot=TRUE, plotGraph=TRUE, tit="OpenMS", wex=1.6,
  specPref=c(conta="LYSC_CHICK", mainSpecies="OS=Homo sapiens"), silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## read OpenMS exported csv
  #fileName <- "C:\\E\\projects\\MassSpec\\smallProj\\2021\\OpenMS\\beforeMSstats\\proteomics_lfq\\out.csv"; 
  #specPref <- list(conta="CON_|LYSC_CHICK", mainSpecies="OS=Saccharomyces cerevisiae", spike="UPS")

  ## initialize
  pepNaCol <- "PeptideSequence"; preCol <- "PrecursorCharge"; condCol <- "Condition"; runCol <- c("BioReplicate","Run"); mzMLCol <- "Reference"
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readOpenMSFile")
  opar <- graphics::par(no.readonly=TRUE)
  chPa <- try(find.package("utils"), silent=TRUE)
  if(inherits(chPa, "try-error")) stop("package 'utils' not found ! Please install first")   
  if(isTRUE(debug)) silent <- FALSE
  if(!isTRUE(silent)) silent <- FALSE
  infoDat <- NULL
  if(debug) message("rmso1")
  ## check & read file
  paFi <- wrMisc::checkFilePath(fileName, path, expectExt="csv", compressedOption=TRUE, stopIfNothing=TRUE, callFrom=fxNa, silent=silent,debug=debug)
  
  ##
  tmp <- list()
  tmp[[1]] <- try(utils::read.csv(paFi, stringsAsFactors=FALSE), silent=TRUE)                # read US csv-file
  tmp[[2]] <- try(utils::read.csv2(paFi, stringsAsFactors=FALSE), silent=TRUE)              # read Euro csv-file
  chCl <- sapply(tmp, inherits, "try-error")
  if(all(chCl)) stop(" Failed to extract data from '",fileName,"'  check format (should be .csv or .csv.gz) & rights to read")
  nCol <- sapply(tmp, function(x) if(length(x) >0) {if(! inherits(x, "try-error")) ncol(x) else NA} else NA)
  bestT <- which.max(nCol)
  datA <- tmp[[bestT]]
  tmp <- NULL                         # reset

  ## locating of columns to treat
  chCoS <- match(runCol, colnames(datA))
  if(all(is.na(chCoS))) stop("Cannot find column indicating replicates/runs")
  chCoS <- if(sum(!is.na(chCoS)) >1) chCoS[1] else wrMisc::naOmit(chCoS)
  chCoQ <- match(quantCol, colnames(datA))
  if(is.na(chCoQ)) stop("Cannot find column indicating intensity values")
  chCoPr <- match(protNaCol, colnames(datA))
  if(is.na(chCoPr)) stop("Cannot find column indicating protein names")
  chCoPe <- match(pepNaCol, colnames(datA))
  if(is.na(chCoPe)) message(fxNa,"Cannot find column indicating peptide names")
  chCoPc <- match(preCol, colnames(datA))
  if(is.na(chCoPc)) message(fxNa,"Cannot find column indicating precursor charge")

  ## extract mzML filenames to 'Run'
  chCoMz <- match(mzMLCol, colnames(datA))
  if(is.na(chCoMz)) {
    message(fxNa,"Cannot find column indicating mzML-file names")
    expSetup <- datA[match(unique(datA[,chCoS]), datA[,chCoS]), match(c(condCol,runCol), colnames(datA))]
  } else {
    expSetup <- datA[match(unique(datA[,chCoMz]), datA[,chCoMz]), match(c(condCol,runCol,mzMLCol), colnames(datA))]
    datA <- datA[,-1*chCoMz]                                            # remove column "Reference" (redundant info, captured in expSetup)
  }
  
  ## main quantitation
  if(length(chCoPe >0)) {                   # "PeptideSequence" found
    ## need to check if 1 prot has mult peptides
      #datA[which(datA[,1]=="P02768ups|ALBU_HUMAN_UPS;sp|ALBU_HUMAN|" & datA[,chCoS]==1),]
    ch1 <- duplicated(paste(datA[,chCoPr],datA[,chCoS]), fromLast=TRUE)
    if(any(ch1)) {
      ## re-sort for easier/faster treating
      nPep <- nrow(datA)
      datA <- datA[order(paste(datA[,chCoPr],datA[,chCoS])),]
      
      useLi <- which(ch1 | duplicated(paste(datA[,chCoPr],datA[,chCoS]), fromLast=FALSE))
        ## summarize peptides ..
        ## ? how to count same peptide if (in same sample) at diff charge states (at diff quant) ?
        ##
        ## NOTE : not yet developed for distinguishing PTM variants !!
        if(!silent) message(fxNa,"summarizing ",nrow(datA)," peptides using ",sumMeth)
        IDs <- unique(datA[,chCoPr]) 
        ## much faster with separate quant summarization and attaching of 'constant' info (compared to 'integrated' function)
        firOf <- match(unique(paste(datA[,chCoPr],datA[,chCoS])), paste(datA[,chCoPr],datA[,chCoS]))
        sumPro <- unlist(by(datA[useLi,c(chCoQ)], paste(datA[useLi,chCoPr],datA[useLi,chCoS]), 
          function(x) if(sumMeth=="top3") {v <- sort(x,decreasing=TRUE); c(mean(v[1:min(3,length(v))]),length(x))} else c(sum(x),length(x)) ) )
        newOrd <- match(paste0(paste(datA[firOf,chCoPr],datA[firOf,chCoS]),"2"), names(sumPro)[2*(1:(length(sumPro)/2))]) 
        byProtSum <- cbind(datA[firOf, c(chCoPr,chCoS)], matrix(as.numeric(sumPro), ncol=2, byrow=TRUE)[newOrd,])
        colnames(byProtSum) <- c(colnames(datA[,c(chCoPr,chCoS,chCoQ)]),"PSM")
        useLr <- if(length(useLi)==nrow(datA)) NULL else (1:nrow(datA))[-1*useLi]           # lines with just one pep
        datA <- cbind(datA[firOf,-ncol(datA)],  byProtSum[,c(-1:0) +ncol(byProtSum)])
        ## optinal filter for min no of peptides (using PSM)
        if(minPepNo[1] > 1) datA <- datA[which(datA[,ncol(datA)] >=minPepNo[1]),]
      } }    
          
  ##
  ## extract quantity (ie unstack)
  samp <- unique(datA[,chCoS])
  IDs <- unique(datA[,chCoPr])
  abund <- array(dim=c(length(IDs),nrow(expSetup),2), dimnames=list(IDs, sub("\\.mzML","",expSetup[,ncol(expSetup)]) ,c("XIC","PSM")))
  for(i in 1:length(samp)) {tmp <- datA[which(datA[,chCoS]==samp[i]), c(chCoPr, c(-1,0) +ncol(datA))]; abund[match(tmp[,1],IDs),i,] <- as.matrix(tmp[,-1]) }

  ## need to recuperate info "BioReplicate" vs "Condition" to grp
  grp <- datA[match(samp,datA[,chCoS]),"Condition"]
  names(grp) <- samp
    
  ## sparate ID from 'ProteinName' entries: remove tailing punctuation or open brackets (ie not closed) at end of fasta header
  annot <- matrix(NA,nrow=nrow(abund), ncol=6, dimnames=list(rownames(abund),c("Accession","EntryName","ProteinName","Species","Contam","SpecType")))
  annot[,1] <- sub("\\|[[:print:]]*$","", sub("^[[:lower:]]{2}\\|","", rownames(abund)))      # extract only 1st ID
  annot[,"ProteinName"] <- sub("[[:upper:]][[:alnum:]]+\\|","", sub("^[[:lower:]]{2}\\|","", rownames(abund)))  # extract ProteinName

  chMult <- grep(";[[:alpha:]]", annot[,chCoPr])                                              # check for multiple concatenated protein-names protNa
  if(length(chMult) >0) {
    annot[chMult,"ProteinName"] <- sub(";[[:alpha:]]+\\|[[:print:]]*","",annot[chMult,"ProteinName"]) }           # extract 1st ProteinName

  ## extract common organism species names (unfortunately OS= term can't be found in csv')
  if(TRUE) {
    commonSpec <- .commonSpecies()                 # allow customizing via argument ?
    for(i in 1:nrow(commonSpec)) {useLi <- grep(commonSpec[i,1], annot[,"ProteinName"])
      if(length(useLi) >0) annot[useLi,"Species"] <- commonSpec[i,2]} }
  
  ## look for special annotation terms
  if(length(specPref) >0) {
    for(i in 1:min(length(specPref),6)) { useLi <- grep(specPref[i], rownames(abund))
      if(length(useLi) >0) annot[useLi,"SpecType"] <- if(nchar(names(specPref)[i]) >0) names(specPref)[i] else c("conta","mainSpe","species2","species3","species4")[i]} 
  }  
          
  ## check for reference for normalization
  refLiIni <- refLi
  if(is.character(refLi) & length(refLi)==1) { refLi <- which(annot[,"SpecType"]==refLi)
    if(length(refLi) <1) message(fxNa," could not find any protein matching argument 'refLi', ignoring ...") else {
      if(!silent) message(fxNa," normalize using subset of ",length(refLi))}}    # may be "mainSpe"
  if(length(refLi) <1) refLi <- NULL
  ## take log2 & normalize
  quant <- wrMisc::normalizeThis(log2(abund[,,1]), method=normalizeMeth, refLines=refLi, callFrom=fxNa) 

  ## plot distribution of intensities
  custLay <- NULL
  if(length(plotGraph) >0) { if(is.numeric(plotGraph)) { custLay <- plotGraph; plotGraph <- TRUE
    } else  {plotGraph <- as.logical(plotGraph[1])}}
  if(plotGraph) {
    if(length(custLay) >0) graphics::layout(custLay) else graphics::layout(1:2)
    graphics::par(mar=c(3, 3, 3, 1))                          # mar: bot,le,top,ri
    if(is.null(tit)) tit <- "OpenMS quantification"
    
    
    
    chGr <- try(find.package("wrGraph"), silent=TRUE)
    chSm <- try(find.package("sm"), silent=TRUE)
    misPa <- c(inherits(chGr, "try-error"), inherits(chSm, "try-error"))
    titSu <- if(length(refLi) >0) paste0(c(" by ",if(length(refLiIni) >1) c(length(refLi)," selected lines") else c("'",refLiIni,"'")),collapse="")  else NULL
    if(any(misPa)) { 
      if(!silent) message(fxNa," missing package ",wrMisc::pasteC(c("wrGraph","sm")[which(misPa)],quoteC="'")," for drawing vioplots")
      ## wrGraph not available : simple boxplot  
      graphics::boxplot(log2(abund[,,1]), main=paste(tit,"(initial)",sep=" "), las=1, outline=FALSE) 
      graphics::abline(h=round(log2(stats::median(abund[,,1],na.rm=TRUE))) +c(-2:2), lty=2, col=grDevices::grey(0.6)) 
      ## plot normalized
      graphics::boxplot(quant, main=paste0(tit," (",normalizeMeth,"-normalized",titSu,")"), las=1, outline=FALSE) 
      graphics::abline(h=round(stats::median(quant, na.rm=TRUE)) +c(-2:2), lty=2, col=grDevices::grey(0.6)) 
    } else {                                            # wrGraph and sm are available
      wrGraph::vioplotW(log2(abund[,,1]), tit=paste(tit,"(initial)",sep=" "), wex=wex, callFrom=fxNa) 
      graphics::abline(h=round(stats::median(log2(abund[,,1]), na.rm=TRUE)) +c(-2:2), lty=2, col=grDevices::grey(0.6)) 
      ## now normalized
      wrGraph::vioplotW(quant, tit=paste0(tit," (",normalizeMeth,"-normalized",titSu,")"), wex=wex, callFrom=fxNa) 
      graphics::abline(h=round(stats::median(quant, na.rm=TRUE)) +c(-2:2), lty=2, col=grDevices::grey(0.6))    
    }
    on.exit(graphics::par(opar)) }   #
  ## meta-data
  notes <- c(inpFile=paFi, qmethod="OpenMS", qMethVersion=if(length(infoDat) >0) unique(infoDat$Software.Revision) else NA, 
    identType="protein", normalizeMeth=normalizeMeth, pepSumMeth=sumMeth, nIniPep=nPep,
    call=deparse(match.call()), created=as.character(Sys.time()), wrProteo.version=paste(utils::packageVersion("wrProteo"), collapse="."), machine=Sys.info()["nodename"])
  ## final output
  if(separateAnnot) list(raw=abund[,,1], quant=quant, annot=annot, counts=abund[,,2], expSetup=expSetup, notes=notes) else data.frame(quant,annot)  
}
  
