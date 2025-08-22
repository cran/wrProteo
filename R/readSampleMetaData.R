#' Read Sample Meta-data from Quantification-Software And/Or Sdrf And Align To Experimental Data
#'
#' Sample/experimental annotation meta-data form \href{https://www.maxquant.org/}{MaxQuant}, ProteomeDiscoverer, FragPipe, Proline or similar, can be read using this function and relevant information extracted.
#' Furthermore, annotation in \href{https://github.com/bigbio/proteomics-sample-metadata}{sdrf-format} can be added (the order of sdrf will be adjated automatically, if possible).
#' This functions returns a list with grouping of samples into replicates and additional information gathered.
#' Input files compressed as .gz can be read as well.
#'
#' @details
#'
#' When initally reading/importing quantitation data, typically very little is known about the setup of different samples in the underlying experiment.
#' The overall aim is to read and mine the corresponding sample-annotation documeneted by the quantitation-software and/or from n sdrf repository and to attach it to the experimental data.
#' This way, in subsequent steps of analysis (eg PCA, statictical tests) the user does not have to bother stuying the experimental setup to figure out which
#' samples should be considered as relicate of whom.
#'
#' Sample annotation meta-data can be obtained from two sources :
#'  a) form additional files produced (and exported) by the initial quantitation software (so far MaxQuant and ProteomeDiscoverer have een implemeneted) or
#'  b) from the universal sdrf-format (from Pride or user-supplied).
#' Both types can be imported and checked in the same run, if valid sdrf-information is found this will be given priority.
#' For more information about the sdrf format please see \href{https://github.com/bigbio/proteomics-sample-metadata}{sdrf on github}.
#'
#'
#' @param quantMeth (character, length=1) quantification method used; 2-letter abbreviations like 'MQ','PD','PL','FP' etc may be used
#' @param sdrf (character, list or data.frame) optional extraction and adding of experimenal meta-data:
#'  This may be a matrix or data.frame with information respective to the experimental setup (to understand which lines=samples should be evaluated as replicates).
#'  If _\code{sdrf} is a character vector, the first entry will be interpreted as path to a local file or the ID/name at ProteomeExchange. 
#'  Here it possible to indicate if specific columns of the sdrf should be skipped (_sdrf=c(myFile, skipCol="abc")_)
#'  For mining the experimental setup/structure \code{sdrf} will get priority over \code{suplAnnotFile}.
#' @param suplAnnotFile (logical or character) optional reading of supplemental files produced by MaxQuant; if \code{gr} is provided, it gets priority for grouping of replicates
#'  if \code{TRUE} in case of \code{method=='MQ'} (MaxQuant) default to files 'summary.txt' (needed to match information of \code{sdrf}) and 'parameters.txt' which can be found in the same folder as the main quantitation results;
#'  if \code{character} the respective file-names (relative ro absolute path), 1st is expected to correspond to 'summary.txt' (tabulated text, the samples as given to MaxQuant) and 2nd to 'parameters.txt' (tabulated text, all parameters given to MaxQuant)
#'  in case of \code{method=='PL'} (Proline), this argument should contain the initial file-name (for the identification and quantification data) in the first position
#' @param path (character) optional path of file(s) to be read
#' @param abund (matrix or data.frame) experimental quantitation data; only column-names will be used for aligning order of annotated samples
#' @param groupPref (list) additional parameters for interpreting meta-data to identify structure of groups (replicates);
#'   May contain \code{lowNumberOfGroups=FALSE} for automatically choosing a rather elevated number of groups if possible (defaults to low number of groups, ie higher number of samples per group).
#'   A vector of custom sample-names may be provided via \code{sampleNames=...} (must be of correct length);
#'   if contains \code{sampleNames="sdrf"} sample-names will be used from trimmed file-names.
#' @param chUnit (logical or character) optional adjustig of group-labels from sample meta-data in case multipl different unit-prefixes are used to single common prefix 
#'   (eg adjust '100pMol' and '1nMol' to '100pMol' and '1000pMol') for better downstream analysis. This option will call \code{\link[wrMisc]{adjustUnitPrefix}} and \code{\link[wrMisc]{checkUnitPrefix}} from package \code{wrMisc}
#'   If \code{character} exatecly this/these unit-names will be searched in sample-names and checked if multiple different decimal prefixes are used; 
#'   if \code{TRUE} the default set of unit-names ('Mol','mol', 'days','day','m','sec','s','h') will be checked in the sample-names for different decimal prefixes
#' @param silent (logical) suppress messages if \code{TRUE}
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This function returns a list with \code{$groups} and \code{$level} (grouping of samples given as integer), and \code{$meth} (method by which grouping as determined).
#'  If valid \code{sdrf} was given, the resultant list contains in addition \code{$sdrfDat} (data.frame of annotation).
#'  Alternatively it may contain a \code{$sdrfExport} if sufficient information has been gathered (so far only for MaxQuant) for a draft sdrf for export (that should be revised and completed by the user).
#'  If software annotation has been found it will be shown in \code{$annotBySoft}.
#'  If all entries are invalid or entries do not pass the tests, this functions returns an empty \code{list}.
#' @seealso This function is used internally by \code{\link{readMaxQuantFile}},\code{/link{readProteomeDiscovererFile}} etc; uses \code{\link{readSdrf}} for reading sdrf-files, \code{\link[wrMisc]{replicateStructure}} for mining annotation columns
#' @examples
#' sdrf001819Setup <- readSampleMetaData(quantMeth=NA, sdrf="PXD001819")
#' str(sdrf001819Setup)
#'
#' @export
readSampleMetaData <- function(quantMeth, sdrf=NULL, suplAnnotFile=NULL, path=".", abund=NULL, groupPref=list(lowNumberOfGroups=TRUE, sampleNames=NULL, gr=NULL), chUnit=TRUE, silent=FALSE, debug=FALSE, callFrom=NULL)  {
  ##  sdrf..()
  ##  suplAnnotFile..(character or logical)
  ##  quantMeth..(character)
  ##  abund..(matrix or data.frame) column-names will be used to comapre & align sample meta-data)

  #### CONTENT
  ## 1.1 SOFTWARE specific META-DATA .. read additional annotation & documentation files
  ##   Aim : extract/build 'summaryD' (& parametersD) allowing to match colnames of 'abund' to suplAnnotFile and/or sdrf
  ## 1.2 basic check of summaryD to quant data, extract supl info for sdrf
  ##   evaluate summaryD to consistent format
  ## 1.3 TRY CHECKING/ADJUSTING ORDER of summaryD
  ## 1.4   replicateStructure
  ### 2  READ SDRF annotation & pick groups of replicates; has priority over grouping based on summary.txt
  ## 2.1 basic check (distinguish full $sampleSetup) form custom data.frame
  ## 2.2 need to match lines (samples) of sdrf (setupDat) to summaryD and/or colnames of abund
  ## 2.3 ready to make setupSd

  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readSampleMetaData")
  if(isTRUE(debug)) silent <- FALSE
  if(!isTRUE(silent)) silent <- FALSE
  summaryD <- parametersD <- setupSdSoft <- setupSd <- sdrfInf <- annSh <- parametersSd <- useSdrfCol <- sdrfDat <- NULL         # initialize    (setupSd needed ?)
  nonUsefulColNa <- c("\\.technical\\.replicate\\.")   # (sdfr-) colnames that might give counter-productive interpretation of experimental structure
  ## checks
  if(length(suplAnnotFile) >1) if(is.na(suplAnnotFile[1])) suplAnnotFile <- NULL
  datOK <- length(sdrf) >0 || length(suplAnnotFile) >0
  if(length(quantMeth) <1) quantMeth <- NA
  if(length(abund) >0 && any(length(dim(abund)) !=2, dim(abund) < 1, na.rm=TRUE)) { datOK <- FALSE
    warning("Invalid argument 'abund'; must be matrix or data.frame with min 1 line and 1 col")}
  if(debug) {message(fxNa,"Ready search & extract sample meta-data   rSM0"); rSM0 <- list(sdrf=sdrf,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,path=path,abund=abund)}


  .corPathW <- function(x) gsub("\\\\", "/", x)
  .adjPat <- function(x) { out <- match(x, unique(x)); names(out) <- if(length(names(x)) >0) names(x) else x; out}   # reduce to integer-pattern (with names)

  .adjTrimPat <- function(x) { x <- wrMisc::rmEnumeratorName(x, incl=c("anyCase","trim0","rmEnum"), sepEnum=c(" ","-","_"), nameEnum=c("Number","No","#","","Replicate","Sample"))
    out <- match(x, unique(x)); names(out) <- names(x); out}   # used

  .redLstToDf <- function(lst) {    # transform lst to data.frame; in case some list-entries have different length, choose the entries with most feq common length
    leL <- sapply(lst, length)
    if(any(duplicated(leL))) {      # need to reduce : find most frequent
      leL2 <- tabulate(leL)
      lst <- lst[which(leL==which.max(leL2))] }
    as.data.frame(lst) }
    
  .replSingleWord <- function(ii, tx, se, replBy="") {   ## for single word in all tx; moove to wrMisc?
    ## NOT USED ANY MORE (when using  wrMisc::rmSharedWords);  moove to wrMisc? 
    ## ii .. word to remove
    ## tx .. ini char-vector
    ## se .. possible separators                           
    if(length(se) >1) se <- paste0("(",paste(se,collapse="|"),")")
    i2 <- grepl(paste0("^",ii), tx)      # heading
    out <- rep(NA, length(tx))
    ## need to protect special characters in ii
    ii <- wrMisc::protectSpecChar(ii)
    if(any(i2)) out[which(i2)] <- sub(paste0("^",ii,se), replBy, tx[which(i2)])       # heading : remove with following sep (if avail)
    if(any(!i2)) out[which(!i2)] <- sub(paste0("(",se,ii,")|(^",ii,")"), replBy, tx[which(!i2)])   # not heading (may be heding now..): remove preceeding sep
    out }

 .trimRedWord <- function(txt, sep=c(" ","_","-","/"), minLe=3, strict=TRUE, silent=TRUE, callFrom=NULL, debug=FALSE) {
    ## replaced by wrMisc::rmSharedWords  (also used inside .checkSetupGroups () )
    ## NOT USED ANY MORE (when using  wrMisc::rmSharedWords) 
    ## function to trim redundant words (@separator) similar to wrMisc::trimRedundText()
    ## strict .. (logical) requires separator to occur in each single character-string to be considered
    ## minLe .. min length for words to be considered (otherwise frequently problem with '1')
    ##
    datOK <- length(txt) >0
    if(datOK) { chNA <- is.na(txt)
      if(all(chNA)) datOK <- FALSE else tx1 <- txt[which(!chNA)] 
    }     
    if(datOK) {
      #strict=TRUE
      chSe <- sapply(sep, function(x) nchar(tx1) > nchar(gsub(x,"",tx1)))
      chS2 <- if(strict) colSums(chSe) ==length(tx1) else colSums(chSe) >0       # if strict require at least instace of 'sep' in each element      
      if(debug) {message(fxNa,"tRW1"); tRW1 <- list(txt=txt,sep=sep,chSe=chSe,chS2=chS2,strit=strict,tx1=tx1)}
      if(any(chS2)) sep <- sep[which(chS2)] else datOK <- FALSE           # reduce to sep found
    }
    if(datOK) {
      allW <- unique(unlist(strsplit(tx1, paste(sep, collapse="|")), use.names=FALSE))
      ## keep only >2 char words
      chLe <- nchar(allW) >= minLe
      if(any(!chLe)) allW <- allW[which(chLe)]       
      ## check all 'words' for recurring in each char-string
      if(length(allW) >0) { 
        chW <- colSums(sapply(allW, grepl, tx1)) ==length(tx1)
        if(any(chW)) {          
          rmWo <- names(chW[which(chW)])
          #chLe <- nchar(rmWo) >0
          if(debug) {message(fxNa,"tRW2"); tRW2 <- list()}
          if(length(rmWo) >0) {
            for(wo in rmWo) tx1 <- .replSingleWord(wo, tx1, sep) 
            txt[which(!chNA)] <- tx1 
          }
          
          if(any(chLe)) {
            rmWo <- rmWo[which(chLe)]
            for(wo in rmWo) tx1 <- .replSingleWord(wo, tx1, sep) 
            txt[which(!chNA)] <- tx1 }
      } }
    }
    txt }

  if( utils::packageVersion("wrMisc") > "1.15.1.1")  .trimRedWord <- wrMisc::rmSharedWords

  .chColOrder <- function(sdr1, sdr2, colNa=c("comment.file.uri.","comment.data.file.")) { out <- NULL
        ## use sdr1 as old/inital, sdr2 as new; return vector for re-establishing init order
        ## NOT USED ANY MORE !!
        for(i in colNa) {
          if(i %in% colnames(sdr1) && i %in% colnames(sdr2)  && sum(duplicated(sdr1[,i])) <1 && length(out) <1) {out <- match(sdr1[,i], sdr2[,i]); break }
        }
        out }
  ## end suppl fx
  

  path <- if(length(path) <1) "." else path[1]
  nSamp0 <- if(length(dim(abund)) >1) ncol(abund) else 0
  chSoft <- c("MQ", "PD", "PL", "FP","MC","AP","IB","NN","SA")
  defUnits <- c("Mol","mol", "days","day","m","sec","s","h")   # for unit-conversion of sample/column-names
  syncColumns <- c(sdrfDat=NA, annotBySoft=NA)
  if(datOK) {
    if("maxquant" %in% tolower(quantMeth)) quantMeth <- "MQ"              # MaxQuant
    if("proteomediscoverer" %in% tolower(quantMeth)) quantMeth <- "PD"    # ProteomeDiscoverer
    if("proline" %in% tolower(quantMeth)) quantMeth <- "PL"               # Proline
    if("fragpipe" %in% tolower(quantMeth)) quantMeth <- "FP"              # FragPipe
    if("masschroq" %in% tolower(quantMeth)) quantMeth <- "MC"             # MassChroq
    if("alphapept" %in% tolower(quantMeth)) quantMeth <- "AP"             # AlphaPeptide
    if("ionbot" %in% tolower(quantMeth)) quantMeth <- "IB"                # IonBot
    if("dia-nn" %in% tolower(quantMeth) || "diann" %in% tolower(quantMeth)) quantMeth <- "NN"   # DiaNN
    if("sage" %in% tolower(quantMeth)) quantMeth <- "SA"                  # Sage
  }

  if(datOK) {  if(length(abund) >0) if(is.null(colnames(abund))) { abund <- NULL
      if(!silent) message(fxNa,"Invalid 'abund' : has NO colnames !") }
    refNSamp <- if(length(dim(abund)) >1) ncol(abund) else NULL  

    ### IMPORT SAMPLE META-DATA, if possible GROUPING OF REPLICATES
    if(length(suplAnnotFile) ==1) {
      if(isFALSE(suplAnnotFile)) suplAnnotFile <- NULL else if(is.na(suplAnnotFile)) suplAnnotFile <- NULL }

   ## 1.1 SOFTWARE specific META-DATA : read additional annotation & documentation files produced by var software as  summaryD & parametersD
    if(length(suplAnnotFile) >0)  {     # read quant software-generated sample annotation
      chFiNa <- NULL                    # initialize
      if(debug) {message(fxNa,"rSM1"); rSM1 <- list(sdrf=sdrf,abund=abund,path=path,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,syncColumns=syncColumns,groupPref=groupPref) }
           ## option 1 : suplAnnotFile has path (do not use default 'path'), use same path for default suplAnnotFile (if applicable)
           ## option 2 : suplAnnotFile has no path, use 'path' for sdrf & suplAnnotFile

      ## Aim : extract/build 'summaryD' (& parametersD) allowing to match colnames of 'abund' to suplAnnotFile and/or sdrf
      

      ## MaxQuant :    (summary.txt & parameters.txt)
      if("MQ" %in% quantMeth && length(suplAnnotFile) >0) {
        isDir <- if(is.character(suplAnnotFile)) utils::file_test("-d",suplAnnotFile[1]) else FALSE
        if(isDir) { path <- suplAnnotFile[1]; suplAnnotFile <- TRUE}
        if(isTRUE(suplAnnotFile)) {      # automatic search for standard file-names ('summary.txt','parameters.txt') in same dir as main MaxQuant data
          chFiNa <- c("summary.txt","summary.txt.gz","parameters.txt","parameters.txt.gz")
          chFi <- file.exists(file.path(path, chFiNa))
          if(debug) {message(fxNa,"rSM0mq\n"); rSM0mq <- list(path=path,sdrf=sdrf,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,chFi=chFi,chFiNa=chFiNa )}
          if(any(chFi, na.rm=TRUE)) { suplAnnotFile <- c(summary=chFiNa[1:2][which(chFi[1:2])[1]], parameters=chFiNa[3:4][which(chFi[3:4])[1]] )
            if(all(names(suplAnnotFile)=="parameters")) suplAnnotFile <- c(NA, parameters=suplAnnotFile$parameters)   # make length=2
            chFi <- c(chFi[1] || chFi[2], chFi[3] || chFi[4])    #needed ?
          } else suplAnnotFile <- NULL
        } else {      # specific/non-default file given
          if(length(suplAnnotFile) >2) suplAnnotFile <- suplAnnotFile[1:2]   # use max length=2
            chFi <- rep(FALSE, 2)
    	    if(!is.na(suplAnnotFile[1])) chFi[1] <- file.exists(file.path(path, suplAnnotFile[1]))
    	    if(!is.na(suplAnnotFile[2])) chFi[2] <- file.exists(file.path(path, suplAnnotFile[2]))
        }
        if(debug) {message(fxNa,"rSM1mq"); rSM1mq <- list(path=path,sdrf=sdrf,summaryD=summaryD,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,path=path,nSamp0=nSamp0,chFiNa=chFiNa,chFi=chFi )}

        ## main reading of MQ sample meta-data
        if(chFi[1]) summaryD <- try(utils::read.delim(file.path(path, suplAnnotFile[1]), stringsAsFactors=FALSE), silent=TRUE)              # 'summary.txt'
        if(chFi[2]) parametersD <- try(utils::read.delim(file.path(path, suplAnnotFile[2]), stringsAsFactors=FALSE), silent=TRUE)           # 'parameters.txt'
        if(inherits(summaryD, "try-error")) {summaryD <- NULL; if(!silent) message(fxNa,"Meta-data: Failed to read '",suplAnnotFile[1],"'  for getting additional information about experiment !")} else {
          summaryD <- if(nrow(summaryD) >2) summaryD[-nrow(summaryD),] else matrix(summaryD[-nrow(summaryD),], nrow=1,dimnames=list(NULL,colnames(summaryD)))  # need to remove last summary-line
          if(debug) message(fxNa,"Successfully read sample annotation from '",suplAnnotFile[1],"'") }
        if(inherits(parametersD, "try-error")) {if(!silent) message(fxNa,"Meta-data: Failed to read '",suplAnnotFile[2],"' !")} else {
          if(debug && chFi[2]) message(fxNa,"Successfully read ",quantMeth," parameters from '",suplAnnotFile[2],"'") }
        syncColumns["annotBySoft"] <- FALSE
        if(debug) { message(fxNa,"rSM1mq2"); rSM1mq2 <- list()}
      }

      ## ProteomeDiscoverer
      ## uses suplAnnotFile as path for '.InputFiles\\.txt'
      if("PD" %in% quantMeth && length(suplAnnotFile) >0) {
        if(debug) {message(fxNa,"rSM1pd"); rSM1pd <- list(sdrf=sdrf,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth)}
        if(length(suplAnnotFile) >1) { if(!silent) message(fxNa,"Only 1st value of argument 'suplAnnotFile' can be used with quantMeth=PD")
          suplAnnotFile <- suplAnnotFile[1] }
        if(isTRUE(suplAnnotFile)) {      # automatic search for standard file-name ('InputFiles.txt') in same dir as main MaxQuant data
          suplAnnotFile <- list.files(path=path, pattern=".InputFiles\\.txt$|.InputFiles\\.txt\\.gz$")
          if(length(suplAnnotFile) >1) { if(!silent) message(fxNa,"Found ",length(suplAnnotFile)," files matching general patter, using ONLY 1st, ie ",suplAnnotFile[1])
            suplAnnotFile <- suplAnnotFile[1] }
          chFi <- length(suplAnnotFile) >0
          if(!chFi && !silent) message(fxNa,"Note: Unable to (automatically) find sample-annotation file. Maybe it was not exported from ProteomeDiscoverer ?")
        } else chFi <- try(file.exists(file.path(path, suplAnnotFile)), silent=TRUE)
        if(inherits(chFi, "try-error") && silent) {chFi <- FALSE; message(fxNa,"Meta-data: Failed to see file '",suplAnnotFile[1]," ! (check if file exists or rights to read directory ?)")}
        if(debug) {message(fxNa,"rSM1pd2");  rSM1pd2 <- list(sdrf=sdrf,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth, chFi=chFi)}

        ## main reading of PD sample meta-data
        if(chFi) summaryD <- try(utils::read.delim(file.path(path, suplAnnotFile[1]), stringsAsFactors=FALSE), silent=TRUE)
        if(inherits(summaryD, "try-error")) {summaryD <- NULL; if(!silent) message(fxNa,"Meta-data: Failed to read '",suplAnnotFile[1],"' !")
        } else {
          syncColumns["annotBySoft"] <- FALSE
          if("File.Name" %in% colnames(summaryD)) {
            chRa <- grep("\\.raw$", tolower(summaryD[,"File.Name"]))
            if(length(chRa) <1) chRa <- grep("\\.raw", tolower(summaryD[,"File.Name"]))
            if(length(chRa) < nrow(summaryD) && length(chRa) >0) {
              if(debug) message(fxNa,"Filter summaryD to '.raw' from ",nrow(summaryD)," to ",length(chRa))
              summaryD <- if(length(chRa) > 1) summaryD[chRa,] else matrix(summaryD[chRa,], nrow=length(chRa), dimnames=list(rownames(summaryD)[chRa], colnames(summaryD))) } 
          } else if("Input.Files.Workflow.ID" %in% colnames(summaryD)) {
            chNeg <- try(as.integer(summaryD[,"Input.Files.Workflow.ID"]), silent=TRUE)
            if(!inherits(chNeg, "try-error")) { chNeg <- chNeg <0
              if(any(chNeg)) summaryD <- if(sum(chNeg) > nrow(summaryD) -2) matrix(summaryD[which(!chNeg),], nrow=sum(!chNeg), dimnames=list(rownames(summaryD)[which(!chNeg)], colnames(summaryD))) else summaryD[which(!chNeg),] }
          }  
          if(debug)  message(fxNa,"ProteomeDiscoverer Meta-data successfully read '",suplAnnotFile[1])}
        if(debug) {message(fxNa,"rSM1pd3"); rSM1pd3 <- list(summaryD=summaryD,parametersD=parametersD,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth, sdrf=sdrf,path=path,nSamp0=nSamp0,chSoft=chSoft,syncColumns=syncColumns)}
      }

      ## Proline
      ## so far only for reading out of xslx
      if("PL" %in% quantMeth && length(suplAnnotFile) >0) {
        if(debug) {message(fxNa,"rSM0pl"); rSM0pl <- list(sdrf=sdrf,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth)}
        summaryD <- NULL
        ## need init filename given via suplAnnotFile
        if(length(grep("\\.xlsx$", suplAnnotFile[1])) >0) {           # won't enter here if suplAnnotFile==NULL
          ## Extract out of Excel
          reqPa <- c("readxl")
          chPa <- sapply(reqPa, requireNamespace, quietly=TRUE)
          if(any(!chPa)) message(fxNa,"package( '",paste(reqPa[which(!chPa)], collapse="','"),"' not found ! Please install first from CRAN") else {
            sheets <- if(debug) try(readxl::excel_sheets(suplAnnotFile[1]), silent=TRUE) else suppressMessages(try(readxl::excel_sheets(suplAnnotFile[1]), silent=TRUE))
            if(debug) {message(fxNa,"rSM2pl"); rSM2pl <- list()}
            if(inherits(sheets, "try-error")) { message(fxNa,"Unable to read file '",suplAnnotFile,"' ! Returning NULL; check format & rights to read")
            } else {
              annShe <- c("Import and filters", "Search settings and infos")     # sheets from xslx to try reading for sample/meta-information
              annSh <- wrMisc::naOmit(match(annShe, sheets))
              annSh <- grep("Import", if(length(annSh) <1) sheets else sheets[annSh])
              if(length(annSh) >1) {
                if(!silent) message(fxNa,"Multipe sheets containing 'Import' found, using 1st :",sheets[annSh[1]])
                annSh <- annSh[1]
              } else if(length(annSh) <1 && !silent) {
                message(fxNa,"Note: NONE of ANNOTATION SHEETS (",wrMisc::pasteC(annShe),") in '",suplAnnotFile,"' FOUND !  Can't check Matching order of samples to sdrf-anotation !")
              }
              summaryD <- as.matrix(as.data.frame(if(debug) readxl::read_xlsx(suplAnnotFile[1], sheet=annSh, col_names=FALSE) else suppressMessages(readxl::read_xlsx(suplAnnotFile[1], sheet=annSh, col_names=FALSE))))
              rownames(summaryD) <- summaryD[,1]
              summaryD <- t(summaryD[,-1])
              rownames(summaryD) <- 1:nrow(summaryD)
              #syncColumns["annotBySoft"] <- FALSE
            }
          }
        } else if(debug) message(fxNa,"Unknown type of sample/experiment annotation file ('",suplAnnotFile[1],"') for Proline, ignoring !!")
      }                 # finish PL

      ## FragPipe
      ##
      if("FP" %in% quantMeth && length(suplAnnotFile) >0) {
        if(debug) { message(fxNa,"rSM1fp1"); rSM1fp1 <- list()}
        ## option 1 : suplAnnotFile has path (do not use default 'path'), use same path for default suplAnnotFile (if applicable)
        ## option 2 : sdrf has no path, use 'path' for sdrf & suplAnnotFile
        ## Aim : extract/build 'summaryD' allowing to match colnames of 'abund' to suplAnnotFile and/or sdrf
        ##  filelist_ionquant.txt & fragpipe-files.fp-manifest
        isDir <- if(is.character(suplAnnotFile)) utils::file_test("-d", suplAnnotFile[1]) else FALSE
        if(isDir) { path <- suplAnnotFile[1]; suplAnnotFile <- TRUE}
        if(isTRUE(suplAnnotFile)) {      # automatic search for standard file-names ('summary.txt','parameters.txt') in same dir as main MaxQuant data
          chFiNa <- c("doNotUseDoNotUse","doNotUseDoNotUse", "fragpipe-files.fp-manifest","fragpipe-files.fp-manifest.gz", "fragpipe.workflow","fragpipe.workflow.gz")
          chFi <- file.exists(file.path(path, chFiNa))
          if(debug) {message(fxNa,"rSM1fp2"); rSM1fp2 <- list(path=path,sdrf=sdrf,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,chFi=chFi,chFiNa=chFiNa )}
          if(any(chFi, na.rm=TRUE)) { suplAnnotFile <- c(summary=chFiNa[1:4][which(chFi[1:4])[1]], parameters=chFiNa[5:6][which(chFi[5:6])[1]] )
            if(all(names(suplAnnotFile)=="parameters")) suplAnnotFile <- c(NA, parameters=suplAnnotFile$parameters)   # make length=2
            chFi <- c(chFi[1] || chFi[2] || chFi[3] || chFi[4], chFi[5] || chFi[6])    # reduce to length=2 (1st for summary, 2nd for parameters)
          } else suplAnnotFile <- NULL
        } else {             # specific/non-default file given (1st for summary, 2nd for parameters)
          if(length(suplAnnotFile) >2) suplAnnotFile <- suplAnnotFile[1:2]   # use max length=2
          chFi <- rep(FALSE, 2)
    	    if(!is.na(suplAnnotFile[1])) chFi[1] <- file.exists(file.path(path, suplAnnotFile[1]))
    	    if(!is.na(suplAnnotFile[2])) chFi[2] <- file.exists(file.path(path, suplAnnotFile[2]))
        }
        if(debug) {message(fxNa,"rSM1fp3"); rSM1fp3 <- list()}

        ## main reading of FP sample meta-data
        if(chFi[1]) summaryD <- try(utils::read.delim(file.path(path, suplAnnotFile[1]), header=FALSE, stringsAsFactors=FALSE), silent=TRUE)
        if(chFi[2]) parametersD <- try(utils::read.delim(file.path(path, suplAnnotFile[2]), header=FALSE, stringsAsFactors=FALSE), silent=TRUE)
        if(inherits(summaryD, "try-error")) { summaryD <- NULL; if(!silent) message(fxNa,"Meta-data: Failed to read '",suplAnnotFile[1],"'  for getting additional information about experiment !")
        } else if(!is.null(summaryD)) {
          msg <- c("File '",suplAnnotFile[1],"' is NOT good annotation file !  Ignoring")
          if(identical(summaryD[1,], c("flag","value"))) { warning(fxNa, msg); summaryD <- NULL}
          if(sum(dim(summaryD) >1) <2) { warning(fxNa, msg); summaryD <- NULL}
          if(length(summaryD) >0) {
            colnames(summaryD) <- c("file","experiment","bioreplicate","dataType")[1:min(ncol(summaryD), 4)]
            summaryD <- as.matrix(summaryD)
            summaryD[,1] <- .corPathW(summaryD[,1])
          }
          #syncColumns["annotBySoft"] <- FALSE
          if(debug) message(fxNa,"Successfully read sample annotation from '",suplAnnotFile[1],"'") }
        if(inherits(parametersD, "try-error")) {if(!silent) message(fxNa,"Meta-data: Failed to read '",suplAnnotFile[2],"' !")
        } else if(!is.null(parametersD))  {
          parametersD <- sub("\\\\:",":", gsub("\\\\\\\\","/", as.character(as.matrix(parametersD))[-(2:3)]))
          if(debug && chFi[2]) message(fxNa,"Successfully read ",quantMeth," parameters from '",suplAnnotFile[2],"'") }
        if(debug) { message(fxNa,"rSM1fp4")}
      }

      ## MassChroq
      if("MC" %in% quantMeth && length(suplAnnotFile) >0) {
        warning(fxNa,"Reading supplemental meta-data from MassChroq is currently not implemented") }

      ## FragPipe
      if("FP" %in% quantMeth && length(suplAnnotFile) >0) {
        warning(fxNa,"Reading supplemental meta-data from FragPipe is currently not implemented") }

      ## Dia-NN
      if("NN" %in% quantMeth && length(suplAnnotFile) >0) {
        warning(fxNa,"Reading supplemental meta-data from Dia-NN is currently not implemented") }

      ## OTHER software ? ..
      if(!any(quantMeth %in% chSoft, na.rm=TRUE) && !silent) message(fxNa,"Note: No specific procedure has been implemented so far for gathering meta-data by the analysis-software/method '",quantMeth,"'")
    }            ## finished main reading of suplAnnotFile into summaryD
    if(debug) { message(fxNa,"rSM2"); rSM2 <- list(sdrf=sdrf,abund=abund,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,summaryD=summaryD,parametersD=parametersD,suplAnnotFile=suplAnnotFile,syncColumns=syncColumns) }
#}}  #  exit after 1.x



    ## 1.2 basic check of summaryD to quant data, extract supl info for sdrf
    if(length(summaryD) >0) {      ##  more checks
      if(length(dim(summaryD)) !=2) summaryD <- matrix(summaryD, ncol=1, dimnames=list(names(summaryD),NULL))
      if(length(abund) <1) { refNSamp <- nrow(summaryD)
          message(fxNa,"Can't verify/correct names of annotation since content of 'abund' has was not given (ie NULL) or has no colnames")
        } else {        
          if(!identical(ncol(abund), nrow(summaryD))) { summaryD <- NULL
		        if(!silent) message(fxNa,"Note : Number of columns of 'abund' does NOT FIT to number of samples in annotation-data !")	}  }
	    }
	  }
    if(debug) { message(fxNa,"rSM3"); rSM3 <- list(sdrf=sdrf,abund=abund,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,summaryD=summaryD,parametersD=parametersD,suplAnnotFile=suplAnnotFile,syncColumns=syncColumns) }

    ## continue evaluating summaryD to consistent format
    if(length(summaryD) >0) {      ## define  setupSdSoft
      ## need to match colnames(abund) to (MQ:) $Raw.file or $Experiment  .. need to find best partial match
      MStype <- "FTMS"                  # used for extracting (more) sdrf info out of parametersSd

      if("MQ" %in% quantMeth) {         ## NOT IN SAME ORDER !!
        useMQSuCol <- c("Raw.file","Experiment","Enzyme","Variable.modifications", "Fixed.modifications","Multi.modifications")
        summaryD <- summaryD[,wrMisc::naOmit(match(useMQSuCol, colnames(summaryD)))]           # cor 21oct22, more cols 7jun23
        chSd <- length(abund) >0 && nrow(summaryD) == ncol(abund)
        ## normally  colnames(abund) and summaryD should alread be in correct order
        if(isTRUE(chSd)) {
          if(!silent && length(abund) >0) if(nrow(summaryD) != ncol(abund)) { message(fxNa,"PROBLEM : Meta-data and abundance data do not match !  ",
            "Number of samples from ",suplAnnotFile[1]," (",nrow(summaryD),") and from main data (",ncol(abund),") do NOT match !! .. ignoring") }
          #if(debug) save(sdrf,abund,suplAnnotFile,quantMeth,summaryD,quantMeth,syncColumns, file="C:\\E\\projects\\TCAmethods\\wrProteoRamus\\rSM4mq.RData")
        }
        if(length(parametersD) >0) {   ## create 'parametersSd' for sdrf
          parametersCol <- paste0(c("MS/MS tol.","MS/MS deisotoping tolerance","MS/MS deisotoping tolerance unit")," (",MStype,")")    # also "Top MS/MS peaks per Da interval." ?
          parametersCol <- c("Modifications included in protein quantification","Match between runs","Fasta file", parametersCol)
          parametersSd <- if(parametersCol[4] %in% parametersD[,1]) parametersD[match(parametersCol[4],parametersD[,1]) ,2] else NA     # eg '20 ppm'
          if(!is.na(parametersSd)) if(grepl("ppm$", parametersSd)) parametersSd <- paste0(1/as.numeric(sub(" ppm$","",parametersSd))," Da")
          fragMassT <- if(all(parametersCol[5:6] %in% parametersD[,1])) paste0( parametersD[match(parametersCol[5:6],parametersD[,1]) ,2], collapse=" ") else NA
          supPar <- parametersD[match(c("Modifications included in protein quantification","Match between runs"), parametersD[,1]), 2]
          parametersSd <- c(precMassTol=parametersSd, fragMassTol=fragMassT, modifs=supPar[1], matchBetwRun=toupper(supPar[2]) )
        } else parametersSd <- c(precMassTol=NA, fragMassTol=NA)
        parametersSd <- c(assayName="run1", label="NT=label free sample (check if correct)", instrum=NA, parametersSd, cleavAgent=paste0("NT=",summaryD[2,"Enzyme"]) )
        ## add PTM modifs ...


        if(debug) { message(fxNa," .. rSM4mq"); rSM4mq <- list(sdrf=sdrf,abund=abund,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,summaryD=summaryD,parametersD=parametersD,syncColumns=syncColumns,parametersSd=parametersSd,MStype=MStype)}
      }

      if("PD" %in% quantMeth) { useCo <- c("Input.Files.","File.ID","File.Name","Instrument.Name")    # no suitable 2nd column ...
        useCo <- wrMisc::naOmit(match(useCo, colnames(summaryD)))
        summaryD <- if(length(useCo) >1) summaryD[,useCo] else matrix(summaryD, ncol=1, dimnames=list(rownames(summaryD), colnames(summaryD)[useCo]))
        if(debug) { message(fxNa,"rSM4pd"); rSM4pd <- list(sdrf=sdrf,useCo=useCo,abund=abund,uplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,summaryD=summaryD,chFiNa=chFiNa) }
        colNa <- wrMisc::trimRedundText(gsub("\\\\","/",as.character(summaryD[,"File.Name"])), silent=debug, debug=debug, callFrom=fxNa)
        if(length(colNa) < ncol(abund)) warning(fxNa,"Trouble ahead : Sample annotation data from ProteomeDiscoverer has FEWER samples than data read !") else {
          if(length(colNa) > ncol(abund)) { message(fxNa,"Note : Sample annotation data from ProteomeDiscoverer has MORE samples than data read, using only first (might be incorrect)")
            colNa <- colNa[1:ncol(abund)]
            summaryD <- summaryD[1:ncol(abund),]
            } }
        ## presume that filenames (from summaryD) are in same order as abund, then trim to file-names (if all in same path)
        ## potential check of order via 'File.ID' to colnames(abund)
        coNa1 <- sub("\\.Sample$","", sub("^Abun[[:lower:]]+\\.","", colnames(abund)))
        sumDOrd <- match(summaryD$File.ID, coNa1)
        chNA <- is.na(sumDOrd)
        if(any(chNA)) { 
          if(!silent) message(fxNa,"NOTE : Unable to match colnames of 'abund' to 'summaryD$File.ID' (NAs at attempt to match) !!")
        } else if(!identical(sumDOrd, 1:ncol(abund))) {
          summaryD <- summaryD[sumDOrd,]
          colNa <- colNa[sumDOrd] }

        colnames(abund) <- colNa                   #
    	  summaryD <- cbind(summaryD, filePath= summaryD[,"File.Name"])               # copy filename+path first to new column
        summaryD[,"File.Name"] <- basename(.corPathW(summaryD[,"File.Name"]))                  # correct to filename only
        syncColumns["annotBySoft"] <- TRUE
        if(debug) { message(fxNa," .. rSM4pd")}
      }

      if("PL" %in% quantMeth) {   ## order OK ?
        chSd <- length(abund) >0 && nrow(summaryD) == ncol(abund)
        ## normally  colnames(abund) and summaryD should alread be in correct order
        if(chSd) {
          # still need to develope extra verification ?
          chCol <- match(c("result_file_name" ,"quant_channel_name","import_params"), colnames(summaryD))
          if(debug) { message(fxNa,"rSM4pl"); rSM4pl <- list(sdrf=sdrf,abund=abund,uplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,summaryD=summaryD,parametersD=parametersD)
          if(all(is.na(chCol))) summaryD <- NULL else {
            parametersD <- summaryD[1, 3:ncol(summaryD)]                                        # how to integrate this later ??
            summaryD <- summaryD[, chCol]
            summaryD[,1] <- sub("\\.mzDB\\.t\\.xml", "", summaryD[,1] )                  # remove Proline spefic file-format extensons
            chFiNa <- colnames(summaryD) %in% "result_file_name"
            if(any(chFiNa)) colnames(summaryD)[which(chFiNa)] <- "File.Name"             # this column should be called 'File.Name'
            summaryD <- as.data.frame(summaryD)
          }                                                                              # adjust to original raw names
        } else {
          if(!silent && nrow(summaryD) != ncol(abund)) message(fxNa,"PROBLEM : Invalid meta-data !  ", "Number of samples from ",
            suplAnnotFile[1]," (",nrow(summaryD),") and from main data (",ncol(abund),") do NOT match !! .. ignoring") }
        }
        syncColumns["annotBySoft"] <- TRUE
      }

      if("FP" %in% quantMeth) {         ## NOT IN SAME ORDER !!
        mat1 <- match(c("file","experiment"), colnames(summaryD))
        if(all(is.na(mat1))) { message(fxNa,"UNABLE to interpret content of ",suplAnnotFile[1]); summaryD <- NULL
        } else {
          summaryD <- cbind(path=dirname(summaryD[,mat1[1]]), Raw.file= basename(summaryD[,mat1[2]]), Experiment=summaryD[,mat1[2]], trimExp=NA)
          summaryD[,4] <- gsub("_+$|-+$|\\.+$| +$|","", sub("[[:digit:]]+$","", wrMisc::trimRedundText(summaryD[,3], side="right", callFrom=fxNa, silent=debug, debug=debug)))   # remove tailing numbers (and tailing redundant text to get to numbers)
          chSd <- length(abund) >0 && nrow(summaryD) == ncol(abund)
          if(length(chSd) <1) chSd <- FALSE
          ## normally  colnames(abund) and summaryD should alread be in same/correct order
          if(!chSd) {
            if(!silent && length(abund) >0) if(nrow(summaryD) == ncol(abund)) { message(fxNa,"PROBLEM : meta-data and abundance data do not match !  ",
              "Number of samples from ",suplAnnotFile[1]," (",nrow(summaryD),") and from main data (",ncol(abund),") do NOT match !! .. ignoring") }
            #if(debug) save(sdrf,abund,suplAnnotFile,quantMeth,summaryD,quantMeth, file="C:\\E\\projects\\TCAmethods\\wrProteoRamus\\rSM4mq.RData")
          }
        }
        syncColumns["annotBySoft"] <- TRUE
        if(debug) { message(fxNa," .. rSM4mq"); rSM4fp <- list()}
      }

      ## other software ? ...

      if(debug) { message(fxNa,"rSM4d"); rSM4d <- list(sdrf=sdrf,abund=abund,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,summaryD=summaryD,parametersD=parametersD,groupPref=groupPref,syncColumns=syncColumns,setupSdSoft=setupSdSoft,sdrf=sdrf) }
#}}}  #  exit after 1.2


      ## 1.3 TRY CHECKING/ADJUSTING ORDER of summaryD
      if(length(abund) >0 && length(summaryD) >0) {
        ## some software specific options otherwise check if filenames can be matched to colnames ?
        ## PD not much possible since colnames  ".F1.Sample",".F2.Sample",".F3.Sample",...
        ## most other software has summaryD in same order as abund
        if("MQ" %in% quantMeth) {   # colnames of abund not necessarly found in summaryD
          summaryD <- wrMisc::matchMatrixLinesToRef(mat=summaryD, ref=colnames(abund), inclInfo=TRUE, silent=TRUE, debug=FALSE, callFrom=fxNa)
          syncColumns["annotBySoft"] <- length(summaryD$newOrder)  >0
          summaryD <- summaryD$mat  }
        if(any(c("PL","FP") %in% quantMeth, na.rm=TRUE)) {
          summaryD <- wrMisc::matchMatrixLinesToRef(mat=summaryD, ref=colnames(abund), inclInfo=TRUE, silent=TRUE, debug=FALSE, callFrom=fxNa)
          syncColumns["annotBySoft"] <- length(summaryD$newOrder)  >0
          summaryD <- summaryD$mat  }
      }
      if(debug) { message(fxNa,"rSM4e"); rSM4e <- list(sdrf=sdrf,abund=abund,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,summaryD=summaryD,parametersD=parametersD,groupPref=groupPref,syncColumns=syncColumns,setupSdSoft=setupSdSoft,sdrf=sdrf) }
#}}}  #  exit after 1.3

      ## 1.4   replicateStructure  (produce setupSdSoft$level)
      grp <- NULL

      ## 1.4.1  (special case) replicate structure based on custom gr (if given as groupPref$gr)
      if(length(groupPref$gr) <1) { setupSdSoft$level <- grp <- .adjPat(groupPref$gr)
      } else { 
        if(length(groupPref$gr)==1) {
          ## now length(groupPref$gr)==1
          if(groupPref$gr=="colnames") grp <- .adjPat(wrMisc::rmEnumeratorName(colnames(abund), incl=c("anyCase","trim0","rmEnum"), sepEnum=c(" ","-","_"), nameEnum=c("Number","No","#","","Replicate","Sample"), silent=silent, debug=debug, callFrom=fxNa))
  
          if(grepl("^sdrf", groupPref$gr)) setupSdSoft$level <- grp <- rep(groupPref$gr, ncol(abund))  # temporal fill (sdrf not read yet)
   
          #if(grepl("^groupPref\\$[[:alpha:]]", groupPref$gr)) {
          #} else if(groupPref$gr=="sdrf")
  
      } }
#}}}  #  exit after 1.4.1      
      
      ## 1.4.2  (special case) replicate structure based on custom colNames (if 'sampleNames' given and 'gr' not given)
      if(length(grp) <1 && length(summaryD) >0 && length(groupPref$sampleNames)==nrow(summaryD) && length(groupPref$gr) <1) {
        grp <- .adjPat(wrMisc::rmEnumeratorName(groupPref$sampleNames, incl=c("anyCase","trim0","rmEnum"), sepEnum=c(" ","-","_"), nameEnum=c("Number","No","#","","Replicate","Sample"), silent=silent, debug=debug, callFrom=fxNa))
        if(sum(duplicated(grp), na.rm=TRUE) <1) {
          grp <- wrMisc::trimRedundText(txt=groupPref$sampleNames, spaceElim=TRUE, silent=debug, debug=debug, callFrom=fxNa)  
          grp <- .adjPat(wrMisc::rmEnumeratorName(groupPref$sampleNames, incl=c("anyCase","trim0","rmEnum"), sepEnum=c(" ","-","_"), nameEnum=c("Number","No","#","","Replicate","Sample"), silent=silent, debug=debug, callFrom=fxNa))
        }
      }
      if(debug) { message(fxNa,"rSM4f"); rSM4f <- list(grp=grp,sdrf=sdrf,abund=abund,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,summaryD=summaryD,parametersD=parametersD,groupPref=groupPref,syncColumns=syncColumns,setupSdSoft=setupSdSoft,sdrf=sdrf) }

      ## 1.4.3  replicateStructure   (of summaryD)
      if(length(grp) <1) {
        setupSdSoft <- wrMisc::replicateStructure(summaryD, silent=silent, debug=debug, callFrom=fxNa)
        chLe <- names(setupSdSoft) %in% "lev"
        if(any(chLe)) names(setupSdSoft)[which(chLe)] <- "level"        # rename 'lev' to 'level' ..
        if(debug) { message(fxNa,"rSM4f"); rSM4f <- list() }

        ## so far no direct information about groups (all filenames are different), need to try to find out (remove enumerators)
        if(all(!duplicated(setupSdSoft$level)) && length(abund) >0) {
          setupSdSoft$sampleNames <- grpA <- wrMisc::trimRedundText(txt=colnames(abund), spaceElim=TRUE, silent=debug, debug=debug, callFrom=fxNa)                 # 26oct22
          colNaGrp <- wrMisc::rmEnumeratorName(grpA, incl=c("anyCase","trim0","rmEnum"), sepEnum=c(" ","-","_"), nameEnum=c("Number","No","#","","Replicate","Sample"), silent=silent, debug=debug, callFrom=fxNa)
          colNaGrPref <- TRUE                                  ## preferential to colnames for searching groups
          if(any(duplicated(colNaGrp)) && colNaGrPref) {       # colnames may be used for designing groups
            setupSdSoft$level <- grp <- .adjPat(colNaGrp)
          } else {
            if(all(setupSdSoft$level ==1:ncol(abund))) {
              ## note : .adjTrimPat() does NOT allow keeping names of levels
              grp2 <- if(ncol(summaryD) >1) apply(summaryD, 2, .adjTrimPat) else as.matrix(.adjTrimPat(summaryD))
              if(ncol(summaryD) >1) { grp3 <- apply(grp2, 2, function(x) length(unique(x)))
                if(any(grp3 < ncol(abund))) {
                  if(length(grp3) >0) { useCol <- if(isTRUE(groupPref$lowNumberOfGroups)) which.min(grp3) else which(grp3 ==stats::median(grp3))[1]
                    setupSdSoft$level <- grp2[,useCol]
                    names(setupSdSoft$level) <- wrMisc::rmEnumeratorName(wrMisc::trimRedundText(txt=summaryD[,useCol], spaceElim=TRUE, silent=debug, debug=debug, callFrom=fxNa),
                      incl=c("anyCase","trim0","rmEnum"), sepEnum=c(" ","-","_"), nameEnum=c("Number","No","#","","Replicate","Sample"), silent=silent, debug=debug, callFrom=fxNa)
                } }
              } else {
                names(grp2) <- wrMisc::rmEnumeratorName(wrMisc::trimRedundText(txt=as.character(summaryD), spaceElim=TRUE, silent=debug, debug=debug, callFrom=fxNa),
                      incl=c("anyCase","trim0","rmEnum"), sepEnum=c(" ","-","_"), nameEnum=c("Number","No","#","","Replicate","Sample"), silent=silent, debug=debug, callFrom=fxNa)
                grp <- setupSdSoft$level <- grp2
              }
            }
          }
          if(debug) { message(fxNa,"rSM4h"); rSM4h <- list() }

        } else { if(!silent) message(fxNa,"Note : Abundance data are ABSENT, CANNOT adjust order of annotation to abundance data")}
      } 

      if(length(grp) >0) { if(length(names(grp)) ==0) names(grp) <- grp
        summaryD <- as.data.frame(cbind(summaryD, grp=names(grp), grpInd=grp))}                 # add presumed grouping to summaryD
    } else grp <- NULL
    if(debug) { message(fxNa,"rSM5"); rSM5 <- list(grp=grp,sdrf=sdrf,grp=grp,abund=abund,groupPref=groupPref,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,summaryD=summaryD,parametersD=parametersD,setupSdSoft=setupSdSoft) }
#}}  #  exit after 1.4.3 
   


    ### 2  READ SDRF annotation & pick groups of replicates; has priority over grouping based on summary.txt
    ###
    if(length(sdrf) >0) {
      useSdrfCol <- NULL          # reset
      ## priority column for groups from sdrf (+ define default colnames for priority)
      chGr <- c("sdrfColumn","sdrfCol")
      chGrPref <- chGr %in% names(groupPref)
      if(any(chGrPref)) {
        groupPref$sdrfColumn <- if(length(groupPref[[chGr[which(chGrPref)[1]]]]) >0) groupPref[[chGr[which(chGrPref)[1]]]] else c("factor.value.disease.","characteristics.disease.", "factor.value.treatment.","characteristics.treatment.","comment.technical.replicate.")
      }

      ## 2.0  Check if 'functional' sdrf (ie list) is provided -> use as is
      iniSdrfOrder <- sdrfDaIni <- sdrfDat <- NULL
      if(is.list(sdrf) && all(c("sdrfDat","col","level") %in% names(sdrf), na.rm=TRUE)) {
        if(debug) message(fxNa,"Custom setupSd provided as sdrf")
        if(all(dim(sdrf$sdrfDat)) >0) {
          sdrfDat <- sdrf$sdrfDat
        } else { sdrf <- NULL
          if(!silent) message(fxNa,"PROBLEM : Invalid custom-sdrf  (should be list containing sdrf$sdrfDat with matrix or data.frame)")          
        }  
        ## ? keep initial sdrf ?# sdrf <- "user provided custom object"
      } else {
        ## 'sdrf' may be  character vector (length <3) => assume path or sdrf accession, 2nd as sdrf-column to use
        ## 'sdrf' may be  (matrix or) data.frame => to use as table to exploit
        if(is.character(sdrf) && length(sdrf) <5 && length(dim(sdrf)) <2) {
          ## read sdrf from file or github
          sdrfDat <- readSdrf(sdrf, silent=silent, debug=debug, callFrom=fxNa)
          if(length(sdrf) >1 &&  !is.list(sdrf)) sdrf <- as.list(sdrf)
          ## skip specified columns
          if(length(sdrf) >1 && "skipCol" %in% names(sdrf)) {
            chCol <- wrMisc::naOmit(if(grepl("^[[:digit:]]+$", sdrf$skipCol)) as.integer(sdrf$skipCol) else which(colnames(sdrfDat) ==sdrf$skipCol) )
            if(length(chCol) >0) { 
              if(!silent) message(fxNa,"Skipping columns ",wrMisc::pasteC(colnames(sdrfDat)[chCol],quoteC="'")," of initial sdrfDat")
              sdrfDat <- sdrfDat[,-chCol, drop=FALSE] }
          }

          ## check for priority columns, retain 1st of them as sdrf[2]
          if(length(groupPref$sdrfColumn) ==1 && length(sdrf) <2) { ch1 <- groupPref$sdrfColumn %in% colnames(sdrfDat)
            if(any(ch1)) sdrf[2] <- which(ch1)[1]
          }
        } else {
          ## user provided custom sample annotation object
          if(length(dim(sdrf)) <2 && !silent) message(fxNa,"Note: 'sdrf' looks bizarre (trouble ahead ?), expecting either file, data.frame or complete list")
          sdrfDat <- sdrf
          sdrf[1] <- "User provided custom object"}
      }
      if(debug) { message(fxNa,"rSM6  dim sdrfDat ",nrow(sdrfDat)," ",ncol(sdrfDat)); rSM6 <- list(sdrf=sdrf,sdrfDat=sdrfDat,setupSd=setupSd,abund=abund,suplAnnotFile=suplAnnotFile,groupPref=groupPref,quantMeth=quantMeth,abund=abund,summaryD=summaryD,parametersD=parametersD,syncColumns=syncColumns) }
#}}}  #  exit after 2.0   ##

      ## 2.1 basic check (distinguish full $sampleSetup) from custom data.frame
      if(length(sdrfDat) >0) {
        syncColumns["sdrfDat"] <- FALSE    # initialize
        if(is.list(sdrfDat) && "sdrfDat" %in% names(sdrfDat)) { 
          if("groups" %in% names(sdrfDat)) groupPref$groups <- sdrfDat$groups
          sdrfDat <- sdrfDat$sdrfDat
          if(debug) message(fxNa,"It seems a full $sampleSetup has been given") }
        if(length(dim(sdrfDat)) <2) sdrfDat <- as.matrix(sdrfDat)
        if(length(abund) >0 && nrow(sdrfDat) != ncol(abund)) {
          if(!silent) message(fxNa,"Note : Ignoring 'sdrf'  : it does NOT have the expected number or rows (",nrow(sdrfDat)," given but ",ncol(abund)," expected !)")
          sdrf <- sdrfDat <- NULL }}
      if(debug) {message(fxNa,"rSM6a"); rSM6a <- list(sdrf=sdrf,sdrfDat=sdrfDat,setupSd=setupSd,abund=abund,suplAnnotFile=suplAnnotFile,groupPref=groupPref,quantMeth=quantMeth,abund=abund,summaryD=summaryD,parametersD=parametersD,syncColumns=syncColumns) }

      ## 2.2  need to match lines (samples) of sdrf (setupDat) to summaryD and/or colnames of abund  (ie bring sdrf in order of abund)
      if(length(sdrfDat) >0) {
        if(length(summaryD) >0) {                             
          ## 2.2.1  summaryD exists,  try matching by file-names  (ie column 'File.Name' and/or col 'filePath')
          chFiNames <- c("File.Name","File","FileName",MQ="Raw.file",PL="raw_file_name","Raw.File","rawfile")     # search in summaryD
          chFiNa <- chFiNames %in% colnames(summaryD)
          if(debug) {message(fxNa,"rSM6a1") }
          if(any(chFiNa, na.rm=TRUE) && "comment.file.uri." %in% colnames(sdrfDat)) {
            ## align sdrfDat by filenames in summaryD
            chFi <- match(sub("\\.zip$|\\.gz$","", basename(.corPathW(summaryD[,chFiNames[which(chFiNa)[1]]]))),
              sub("\\.zip$|\\.gz$","", basename(.corPathW(sdrfDat[,"comment.file.uri."]))))  # new order
            if(any(is.na(chFi)) && any(grepl("\\.raw",sdrfDat[,"comment.file.uri."]), na.rm=TRUE)) {
              ## try to align fileNames to adjust order of sdrf :
              sumDaFiNa <- sub("\\.raw","",sub("\\.zip$|\\.gz$","", basename(.corPathW(summaryD[,chFiNames[which(chFiNa)[1]]]))))
              sdrfFiNa <- sub("\\.raw","",sub("\\.zip$|\\.gz$","", basename(.corPathW(sdrfDat[,"comment.file.uri."]))))
              chFi <- match(sumDaFiNa, sdrfFiNa)  # new order
              if(any(is.na(chFi)) && "FP" %in% quantMeth) {
                sumDaFiNa <- wrMisc::rmEnumeratorName(wrMisc::trimRedundText(txt=sumDaFiNa, spaceElim=TRUE, silent=debug, debug=debug, callFrom=fxNa), newSep="_", incl=c("anyCase","trim0"), silent=silent, debug=debug, callFrom=fxNa)
                sdrfFiNa <- wrMisc::rmEnumeratorName(wrMisc::trimRedundText(txt=sdrfFiNa, spaceElim=TRUE, silent=debug, debug=debug, callFrom=fxNa), newSep="_", incl=c("anyCase","trim0"), silent=silent, debug=debug, callFrom=fxNa)
                chFi <- match(sumDaFiNa, sdrfFiNa)     # new order
                if(debug) { message(fxNa,"rSM6aa  dim sdrfDat ",nrow(sdrfDat)," ",ncol(sdrfDat)); rSM6aa <- list(sdrfDat=sdrfDat,summaryD=summaryD,setupSd=setupSd,quantMeth=quantMeth, sumDaFiNa=sumDaFiNa,sdrfFiNa=sdrfFiNa,chFi=chFi)}
              }
              rmRaw <- TRUE
            } else rmRaw <- FALSE
            if(sum(is.na(chFi)) >0) { warning(fxNa,"UNABLE to match all filenames from sdrf and ",basename(.corPathW(suplAnnotFile)),
               " ! \n  ++ BEWARE : Grouping of replicates may be incorrect !! \n") 
            } else {
              ## Adjust order of sdrf to data (summaryD)
              if(!silent && rmRaw) message(fxNa,"Note : Some filenames contain '.raw', others do NOT; solved inconsistency ..")
              iniSdrfOrder <- order(chFi)   #(1:nrow(sdrfDat))[chFi]
              #later#iniSdrfOrder <- sdrfDat$iniSdrfOrd
              #sdrfDat$iniSdrfOrd <- iniSdrfOrder 
              sdrfDat <- sdrfDat[chFi,]
              if(!silent) message(fxNa,"Successfully adjusted order of sdrf to content of ",basename(.corPathW(suplAnnotFile)))
            }
            syncColumns["sdrfDat"] <- TRUE
          } else if(!silent) message(fxNa, if(debug) "rSM6a  "," summaryD exists, but unable to find file-names")
          if(debug) { message(fxNa,"rSM6a1  dim sdrfDat ",nrow(sdrfDat)," ",ncol(sdrfDat)); rSM6a1 <- list(sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,suplAnnotFile=suplAnnotFile,groupPref=groupPref,quantMeth=quantMeth,abund=abund,summaryD=summaryD,parametersD=parametersD,syncColumns=syncColumns,iniSdrfOrder=iniSdrfOrder) }

        } else {             
          ## 2.2.2  no summaryD, try colnames of abund to check order
          if(length(abund) >0 && length(dim(abund)) >1 && ncol(abund)==nrow(sdrfDat)) {   ## valid abund, dimensions matching sdrf
            if(debug) { message(fxNa,"rSM6b  dim sdrfDat ",nrow(sdrfDat)," ",ncol(sdrfDat)); rSM6b <- list(sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,uplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,summaryD=summaryD,parametersD=parametersD,syncColumns=syncColumns) }
            ## requires utils::packageVersion("wrMisc") > "1.11.1"
            sdrfDaIni <- sdrfDat
            chIB <- "IB" %in% quantMeth || all(grepl("(^Intensity_)", colnames(abund)) & grepl("\\.raw",sdrfDaIni$comment.data.file.))   # facilitate aligning for case of IB 

            if(chIB) {
              sdrfDat <- cbind(sdrfDat, iniSdrfOrd=1:nrow(sdrfDat), matched=match(sub("\\.raw","",sdrfDaIni$comment.data.file.),sub("^Intensity_","",colnames(abund))))  # step 1
              if(all(is.na(sdrfDat[,"matched"]))) warning(fxNa,"Unable to match sample-names to sdrf !!  rSM6b2")
              if(sum(duplicated(sdrfDat[,"matched"])) ==0) { sdrfDat <- sdrfDat[order(sdrfDat[,"matched"]), 1:(ncol(sdrfDaIni) +1)]  # step 2, in order of abund
              } else sdrfDat <- wrMisc::matchMatrixLinesToRef(mat=sdrfDat, ref=colnames(abund), addRef=TRUE, silent=silent, debug=debug, callFrom=fxNa)  # 2way-grep
            } else sdrfDat <- wrMisc::matchMatrixLinesToRef(mat=cbind(sdrfDat, iniSdrfOrd=1:nrow(sdrfDat)), ref=colnames(abund), addRef=TRUE, silent=silent, debug=debug, callFrom=fxNa)  # 2way-grep
            if(debug) { message(fxNa,"rSM6b2  dim sdrfDat ",nrow(sdrfDat)," ",ncol(sdrfDat)); rSM6b2 <- list(sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,uplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,summaryD=summaryD,parametersD=parametersD,syncColumns=syncColumns) }
          
          } else {  
            ## check matching ?
            if(length(sdrfDat) <1 && length(abund) >0) {                                       ## so far failed to align - further trim names
              if(debug) { message(fxNa,"Failed to align to sdrf - further trim names from abund    rSM6b3   ")}
              ## now look for bad separator '.' before text and remove
              colNaAbund <- colnames(abund)
              ch1 <- grep("[[:digit:]]\\.[[:alpha:]]", colNaAbund)
              if(any(ch1)) {
                selLoc <- sapply(gregexpr("[[:digit:]]\\.[[:alpha:]]", colNaAbund[ch1]), function(x) x[[1]])
                colNaAbund[ch1] <- paste0(substr(colNaAbund[ch1],1,selLoc), substring(colNaAbund[ch1], selLoc+2)) }
              sdrfDat <- wrMisc::matchMatrixLinesToRef(mat=cbind(sdrfDaIni, 1:nrow(sdrfDat)), ref=colNaAbund, exclCol=ncol(sdrfDat)+1, addRef=TRUE, silent=silent, debug=debug, callFrom=fxNa)  # 2way-grep
              if(length(sdrfDat) <1) {
                colNaEnum <- all(grepl("_[[:digit:]]+$", colNaAbund))
                if(colNaEnum) { tm1 <- sub("_[[:digit:]]+$","", colNaAbund)
                  colNaAbund2 <- sub("\\..+","", substr(colNaAbund, 1, nchar(tm1)))
                  colNaAbund3 <- paste0(colNaAbund2,substring(colNaAbund, nchar(colNaAbund) -1),"$")    # without repeated text after 1st '.'
                  ## Adjust order of sdrf to data ()
                  sdrfDat <- wrMisc::matchMatrixLinesToRef(mat=cbind(sdrfDaIni, 1:nrow(sdrfDat)), ref=colNaAbund3, addRef=TRUE, exclCol=ncol(sdrfDat)+1, silent=silent, debug=debug, callFrom=fxNa)  # 2way-grep
                }
              }
              if(length(sdrfDat) <1 && !silent) message(fxNa,"PROBLEM : FAILED to align sdrf to actual colnames of data !!!  rSM6b4")
            }
            if(debug) { message(fxNa,"rSM6b3  dim sdrfDat ",nrow(sdrfDat)," ",ncol(sdrfDat)); rSM6b3 <- list(sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,uplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,summaryD=summaryD,parametersD=parametersD,syncColumns=syncColumns) }
             
            #iniSdrfOrder <- sdrfDat$iniSdrfOrd
            sdrfDat <- sdrfDat[,-ncol(sdrfDat), drop=FALSE]   # why remove col 'ref' ?
            if(debug) { message(fxNa,"rSM6c1  dim sdrfDat ",nrow(sdrfDat)," ",ncol(sdrfDat)); rSM6c1 <- list() }            
            rm(sdrfDaIni)                                                                                                                                                                             
            syncColumns["sdrfDat"] <- TRUE                         # really sure that synchronization successful ?
          } #else  if(!silent) message(fxNa,"Note : NO Additional information on filenames-order found, can't correct/adjust sdrf (ie sdrfDat) !!", if(debug) "   rSM6a3")
          refNSamp <- nrow(sdrfDat)
          if(debug) { message(fxNa,"rSM6c2  dim sdrfDat ",nrow(sdrfDat)," ",ncol(sdrfDat)); rSM6c2 <- list(sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,suplAnnotFile=suplAnnotFile,groupPref=groupPref,quantMeth=quantMeth,abund=abund,summaryD=summaryD,parametersD=parametersD,syncColumns=syncColumns,iniSdrfOrder=iniSdrfOrder,chUnit=chUnit) }
        }
      }
      if(debug) { message(fxNa,"rSM6d  dim sdrfDat ",nrow(sdrfDat)," ",ncol(sdrfDat)); rSM6d <- list(sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,groupPref=groupPref, nonUsefulColNa=nonUsefulColNa, summaryD=summaryD,parametersD=parametersD,syncColumns=syncColumns,iniSdrfOrder=iniSdrfOrder,setupSd=setupSd,defUnits=defUnits,chUnit=chUnit) }
#}}}  #  exit after 2.3   ##


      ## 2.3  ready to make setupSd
      chLe <- names(setupSd) %in% "lev"
      if(any(chLe)) setupSd <- setupSd[-which(chLe)]       # remove $lev - if existing
      
      if(length(sdrfDat) >0 && all(dim(sdrfDat) >0)) {     
        ## 2.3.1  sdrfDat is matrix or data.frame

        ## check for custom-provided sampleNames  (priority)
        if(length(groupPref$sampleNames) ==nrow(sdrfDat)) { 
          ## 2.3.1.1  custom sampleNames => use later (instead of colnames/file-names)
          setupSd$sampleNames <- groupPref$sampleNames          # 
        } else if((length(groupPref$sampleNames) <1 || "sdrf" %in% groupPref$sampleNames) && length(dim(sdrfDat))==2) {
          ## 2.3.1.2  use sample--names from sdrf (groupPref$sampleNames :  empty or "sdrf" specified)
          if("source.name" %in% colnames(sdrfDat)) useSdrfCol <- which(colnames(sdrfDat)=="source.name") else {
            ## try finding potential column for sample names
            ch1 <- colSums(apply(sdrfDat, 2, duplicated)) ==0
            if(any(ch1)) useSdrfCol <- which(ch1)[1] else {useSdrfCol <- which.min(ch1)[1]; 
              if(!silent) message(fxNa,"Bizzare, can't find any column with all different (unique) content to use for Sample Names !!??!!    (..using column '",colnames(sdrfDat)[useSdrfCol],"' with least repetitions)")}
          }  
          setupSd$sampleNames <- sdrfDat[,useSdrfCol[1]]
          ## further check if useSdrfCol just contains 'Sample 1' or 'run 1' etc
          ch1 <- grepl("(^sample)|(^run) {0,2}[[:digit:]]+$", tolower(setupSd$sampleNames))
          avoidCol <- c( unlist(lapply(tolower(nonUsefulColNa), grep, tolower(colnames(sdrfDat)))))
          if(all(ch1)) {
            ch2 <- colSums(apply(sdrfDat[,-avoidCol, drop=FALSE], 2, duplicated))
            if(any(ch2 ==0)) {
              ## possible new sample-names
              setupSd$sampleNames <- sdrfDat[,((1:ncol(sdrfDat))[-avoidCol])[which(ch2 == 0)[1]]]
            } else {
              if(!silent) message(fxNa,"Having trouble identifying possibly meaningful sample-names (all look like 'sample 1', etc)")
            } 
            ### mine for groups while info is available - but mines only single column !!
            #if(any(ch2 > 0 & ch2 < nrow(sdrfDat) -1)) {  # use only if some replicates found
            #  ## use new sampleNames from here
            #  setupSd$groups <- setupSd$iniGroups <- sdrfDat[,((1:ncol(sdrfDat))[-avoidCol])[which(ch2 > 0 & ch2 < nrow(sdrfDat) -1)[1]]]
            #} 
          }          
          setupSd$iniSampleNames <- setupSd$sampleNames          
          setupSd$sampleNames <- wrMisc::rmSharedWords(setupSd$sampleNames, sep=c("_","-",",",".","=",";"), callFrom=fxNa, silent=!debug)
          if(!identical(sdrfDat[,useSdrfCol], setupSd$sampleNames)) setupSd$fullSampleNames <- sdrfDat[,useSdrfCol]

          #if(length(setupSd$groups) ==nrow(sdrfDat)) {
          #  setupSd$level <- names(setupSd$groups) <- .adjPat(setupSd$groups)
          #  names(setupSd$level) <- setupSd$groups
          #}
          
        } else avoidCol <- unlist(lapply(tolower(nonUsefulColNa), grep, tolower(colnames(sdrfDat))))
        if(debug) { message(fxNa,"rSM6d0  dim sdrfDat ",nrow(sdrfDat)," ",ncol(sdrfDat)); rSM6d0 <- list(sdrf=sdrf,sdrfDat=sdrfDat,useSdrfCol=useSdrfCol,abund=abund,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,groupPref=groupPref, nonUsefulColNa=nonUsefulColNa, summaryD=summaryD,parametersD=parametersD,syncColumns=syncColumns,iniSdrfOrder=iniSdrfOrder,setupSd=setupSd,defUnits=defUnits,chUnit=chUnit,avoidCol=avoidCol) }
        

        ## 2.3.2  groups of replicates 
        ##  when to consider/integrate info from summaryD ?
        newNa <- grp <- chColNa <- useSdrfCol <- NULL               ## initialize
        ## groupPref$gr could be text for colname(s) or integer for col-index; if nothing found use all cols

        if(isTRUE(length(groupPref$gr)==refNSamp)) {
          ## groupPref$gr is given; check if index
          setupSd$gr <- setupSd$iniGr <- groupPref$gr
          if(is.character(groupPref$gr)  && all(grepl("^[[:digit:]]+$", groupPref$gr))) { 
            setupSd$gr <- try(as.integer(groupPref$gr))
            if(inherits(setupSd$gr, "try-error") || any(setupSd$gr <0 | setupSd$gr > refNSamp, na.rm=TRUE)) { setupSd$gr <- NULL
              if(!silent) message("Having difficulty understanding argument groupPref$gr (all digits but not valid as index), ignoring")
            } else if(length(colnames(abund)) >0) setupSd$gr <- colnames(abund)[setupSd$gr]
          } 
          setupSd$gr <- wrMisc::rmSharedWords(setupSd$gr, sep=c("_","-",",",".","=",";"), callFrom=fxNa, silent=!debug)
        } else {
          ## 2.3.2.2 mining sdrf for groups of replicates 
          ## groupPref$gr not explicitely given (but may be colname(s) from sdrfDat to use), need to define which cols of sdrfDat to mine
          if(debug) { message(fxNa,"rSM6d1  dim sdrfDat ",nrow(sdrfDat)," ",ncol(sdrfDat)); rSM6d1 <- list(sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,groupPref=groupPref,avoidCol=avoidCol,summaryD=summaryD,parametersD=parametersD,syncColumns=syncColumns,iniSdrfOrder=iniSdrfOrder,setupSd=setupSd,chColNa=chColNa,chUnit=chUnit) }

          ## automatic mining (from sdrf) (eg groupPref=list(gr="sdrf") or groupPref=list(gr="sdrf$source.name") or groupPref=list(gr=c("sdrf$source.name","sdrf$xy"),combMeth="median") .. last will be used for choosing
          
          ## extact columns names of sdrf to examine, return index
          if(length(groupPref$gr) >0 && all(grepl("^sdrf\\$",groupPref$gr))) {
            chColNa <- sub("^sdrf\\$","", groupPref$gr)
            if(all(grepl("^[[:digit:]]+$"))) { chColNa <- try(as.integer(groupPref$gr))
              if(inherits(chColNa, "try-error") || any(chColNa <0 | chColNa > refNSamp, na.rm=TRUE)) { chColNa <- NULL
                if(!silent) message("Having difficulty understanding argument groupPref$gr (invalid index for sdrfDat), ignoring") }            

            } else chColNa <- chColNa[which(colnames(sdrfDat) %in% chColNa)]    # transform to index
          } else chColNa <- NULL      # custom given like : groupPref=list(gr="sdrf$source.name")
          #chColNa <- if(length(groupPref$gr) >0 && all(grepl("^sdrf\\$",groupPref$gr))) sub("^sdrf\\$","", groupPref$gr) else NULL      # custom given like : groupPref=list(gr="sdrf$source.name")
          ## define which cols of sdrf to check further
          nonUsefulColNa <- c("\\.technical\\.replicate\\.")
          if(length(nonUsefulColNa) >0) nonUsefulColNa <- grep(nonUsefulColNa, colnames(sdrfDat))
          if(length(avoidCol) >0) {avoidCol <- if(is.character(avoidCol) && !all(grepl("^[[:digit:]]+$", avoidCol))) wrMisc::naOmit(match(sub("^sdrf\\$","",avoidCol), colnames(sdrfDat))) else as.integer(avoidCol)}   # combine avoidCol (txt or index) & nonUsefulColNa (txt)
          avoidCol <- wrMisc::naOmit(unique(nonUsefulColNa, avoidCol))
          #if(length(nonUsefulColNa) >0) avoidCol <- wrMisc::naOmit(unique(avoidCol, unlist(lapply(tolower(nonUsefulColNa), grep, tolower(colnames(sdrfDat))))))
          ## final assignment of columns to check for defining groups
          if(length(chColNa) ==0) {
            chColNa <- if(length(avoidCol) >0) (1:ncol(sdrfDat))[-avoidCol] else 1:ncol(sdrfDat)
            if(!silent) message(fxNa,"Unable to find initially designed colnames for mining of sdrf, now using all")
            groupPref$meth <-paste("Using all columns (except avoidCol,nonUsefulColNa) from sdrf for defining groups", if(length(setupSd$iniGr) >0) "  (had trouble finding columns indicated -> default uninsg all)")
          } else {
            if(length(avoidCol) >0) chColNa <- chColNa[which(!chColNa %in% avoidCol)]
            if(length(chColNa) >0) groupPref$meth <- paste("Using specified ",length(chColNa)," column(s) from sdrf for defining groups") else {
              chColNa <- if(length(avoidCol) >0) (1:ncol(sdrfDat))[-avoidCol] else 1:ncol(sdrfDat)
              if(!silent) message(fxNa,"Unable to find initially designed colnames for mining of sdrf, now using all (except avoidCol,nonUsefulColNa)")
            }
          }
          groupPref$chCol <- chColNa    # (verified) integer realtive to column of sdrfDat

          #if(length(chColNa) ==1) {
          #  ## if single is meant to be column index
          #  if(chColNa %in% names(sdrfDat)) groupPref$chCol <- which(colnames(sdrfDat) ==chColNa)[1] else groupPref$chCol <- chColNa
          #} else {
          #  chColNa <- NULL
          #  if(!silent && !identical(groupPref$gr, "sdrf")) message(fxNa,"NOTE : '",groupPref$gr[1],"' does NOT match column of sdrfDat !  (ignoring)")
          #  groupPref$chCol <-  if(length(avoidCol) >0) 1:nrow(sdrfDat)[,-avoidCol] else 1:nrow(sdrfDat) }  # set to default  
          
          #if(length(groupPref$gr) ==0 || length(chColNa)==0) {                      
          #  ## 2.3.2.2.1  use all columnds of sdrf
          #  groupPref$chCol <- if(length(avoidCol) >0) 1:nrow(sdrfDat)[,-avoidCol] else 1:nrow(sdrfDat)
          #  groupPref$meth <- "using all columns from sdrf"
          #} else {
          #  ## 2.3.2.2.2  use specific column of sdrf
          #  groupPref$meth <- paste0("groups as specified from sdrf-column '",setupSd$gr,"'")
          #  groupPref$chCol <- sdrfDat[,chColNa]
          #}
          #if(length(groupPref$gr) <1) {
          #  groupPref$chCol <- if(length(avoidCol) >0) 1:nrow(sdrfDat)[,-avoidCol] else 1:nrow(sdrfDat)
          #  groupPref$meth <- "(had trouble finding columns indicated -> default) using all columns from sdrf"
          #}

          if(length(groupPref$combMeth) ==0) groupPref$combMeth <- "def"
          if(debug) {message(fxNa,"length setupSd ", length(setupSd),"  rSM6g"); rSM6g <- list(setupSd=setupSd,sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,suplAnnotFile=suplAnnotFile,groupPref=groupPref,summaryD=summaryD,parametersD=parametersD,setupSdSoft=setupSdSoft,quantMeth=quantMeth,useSdrfCol=useSdrfCol,chUnit=chUnit) }  #useSdrfCol=useSdrfCol,
          replStrOpt <- c("default","def","highest","high","lowest","min","max","median","combAll","combineAll","combNonOrth")  # possible options for groupPref$combMeth
          chMeth <- groupPref$combMeth %in% replStrOpt
          if(length(chMeth) <1 || !isTRUE(all(chMeth))) { 
            if(length(chMeth) >0) message(fxNa,"Unable to understand argument groupPref$combMeth, setting to default")
            groupPref$combMeth <- "def"
          } else groupPref$combMeth <- groupPref$combMeth[which(chMeth)]      # default in doubt

          ## which columns of sdfr to consider
          if(length(groupPref$chCol) >1) {    
            ## 2.3.2.3   setupSd$gr : consider multiple columns of sdrf
            ## use/mine MULTIPLE/all columns of sdrf :  groupPref=list(gr=c("sdrf$source.name","sdrf$xy"))
            useSdrfCol <- if("useCol" %in% names(groupPref)) match(groupPref$useCol, colnames(sdrfDat)) else {
              if(is.numeric(groupPref$chCol)) groupPref$chCol else match(sub("^sdrf\\$","",groupPref$chCol), colnames(sdrfDat))}  ## check & tranform to index
            useSdrfCol <- wrMisc::naOmit(useSdrfCol)
            if(length(useSdrfCol) <1) {useSdrfCol <- 1:ncol(sdrfDat); groupPref$meth <- "(nothing matching initially) using all columns from sdrf"
               if(!silent) message(fxNa,"none of groupPref$gr match colnames of sdrfDat, using all columns")} 
            #useSdrf <- sdrfDat[,useSdrfCol, drop=FALSE]
            
            nonUsefulColNa <- c("\\.technical\\.replicate\\.")     # coluln name of sdrf to exclude ( to beginning ?)

            ## exclude 'technicalReplicates' from argument groupPref
            if(!("excludeReplicates" %in% names(groupPref) && isFALSE(groupPref$excludeTechnReplicates))) {       # exclude columns with replicate-numbers unless   groupPref=list(excludeTechReplicates=FALSE)
              ## exclude columns labeled 'replicate'
              chCo2 <- unlist(lapply(nonUsefulColNa, grep, colnames(sdrfDat)[useSdrfCol]))
              if(length(chCo2) >0) useSdrfCol <- useSdrfCol[-1*(chCo2)]             
            }             
            
            ## check groupPref$combMeth
            if(length(groupPref$combMeth) ==0) groupPref$combMeth <- "def"
            if(debug) {message(fxNa,"length setupSd ", length(setupSd),"  rSM6g"); rSM6g <- list(setupSd=setupSd,sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,suplAnnotFile=suplAnnotFile,groupPref=groupPref,summaryD=summaryD,parametersD=parametersD,setupSdSoft=setupSdSoft,quantMeth=quantMeth,useSdrfCol=useSdrfCol,chUnit=chUnit) }  #useSdrfCol=useSdrfCol,
            replStrOpt <- c("default","def","highest","lowest","min","max","median","combAll","combNonOrth")  # possible options for groupPref$combMeth
            chMeth <- groupPref$combMeth %in% replStrOpt
            groupPref$combMeth <- if(length(chMeth) <1 || !isTRUE(all(chMeth))) "def" else groupPref$combMeth[which(chMeth)]      # default in doubt

            setupSdIni <- setupSd
            replMeth <- "failed"             # initialize
            ## check if 'combMeth' is 'default' or 'def'
            if(isTRUE(grepl("^def", groupPref$combMeth))) groupPref$combMeth <- c("combNonOrth","lowest")
            
            if(length(useSdrfCol) >1) {
              if(length(groupPref$combMeth) >1) {
                ## 'default : test 2 methods ('lowest' & 'combNonOrth'), then 
                tmp <- list(combNonOrth=try(wrMisc::replicateStructure(sdrfDat[,useSdrfCol, drop=FALSE], method=groupPref$combMeth[1], silent=silent, callFrom=fxNa, debug=debug), silent=TRUE),
                  lowest=try(wrMisc::replicateStructure(sdrfDat[,useSdrfCol , drop=FALSE], method=groupPref$combMeth[2], silent=silent, callFrom=fxNa, debug=debug), silent=TRUE))
              } else {
                tmp <- list(a=try(wrMisc::replicateStructure((sdrfDat[,useSdrfCol, drop=FALSE]), method=groupPref$combMeth[1], silent=silent, callFrom=fxNa, debug=debug), silent=TRUE))}
              if(debug) {message(fxNa," rSM6g2"); rSM6g2 <- list(setupSd=setupSd,sdrf=sdrf,sdrfDat=sdrfDat,tmp=tmp,abund=abund,suplAnnotFile=suplAnnotFile,groupPref=groupPref,summaryD=summaryD,parametersD=parametersD,setupSdSoft=setupSdSoft,quantMeth=quantMeth,useSdrfCol=useSdrfCol,chUnit=chUnit) }  #useSdrfCol=useSdrfCol,
              ## evaluate ..
              ch1 <- sapply(tmp, inherits, "try-error")
              if(all(ch1)) { message(fxNa,"UNABLE to understand replicate-structure from sdrf !!"); setupSd <- NULL; syncColumns["sdrfDat"] <- FALSE
              } else {
                ## 2.3.2.3.2  init mining was successful, now continue extracting
                if(any(ch1)) tmp <- tmp[which(!ch1)]
                if(length(tmp) >1) {
                  ## choose among multiple options for grouping (number of groups) based on groupPref$lowNumberOfGroups 
                  ch1 <- sapply(tmp, function(x) length(x$lev[which(!duplicated(x$lev))]))
                  lowNumberOfGroups <- if(length(groupPref$lowNumberOfGroups)==1) isTRUE(groupPref$lowNumberOfGroups) else TRUE     # lowNumberOfGroups defaults to TRUE
                  useSe <- if(any(ch1 ==1, na.rm=TRUE)) which(ch1 !=1) else if(lowNumberOfGroups) which.min(ch1) else which.max(ch1)
                  replMeth <- useSe <- useSe[1]
                  if(!silent) message(fxNa,"Choosing model '",names(useSe),"' for evaluating replicate-structure (ie ",ch1[useSe[1]]," groups of samples)", if(debug) "rSM6g2")         
                  tmp <- tmp[[useSe[1]]]
                }
                if(debug) {message(fxNa,"length setupSd ", length(setupSd),"  rSM6h"); rSM6h <- list(setupSd=setupSd,sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,suplAnnotFile=suplAnnotFile,groupPref=groupPref,summaryD=summaryD,parametersD=parametersD,setupSdSoft=setupSdSoft,quantMeth=quantMeth,tmp=tmp,replMeth=replMeth) }  #useSdrfCol=useSdrfCol,
  
                ## 2.3.2.2.2.2  now need to extract/evaluate 'tmp' into setupSd
                ## NEED TO CHECK below  17aug25
                setupSd$level <- tmp$lev
                setupSd$methGr <- paste0("groups based on mining '",groupPref$gr,"', ie '",tmp$meth,"'")
                ## need to construct/extract likely names
                setupSd$groups <- names(tmp$lev)
                if(debug) message(fxNa," rSM6h.3")
              }

            }               
          } else {     ## single column
            ## 2.3.2.4  
            useSdrfCol <- groupPref$chCol
            setupSd$level <- .adjPat(sdrfDat[,useSdrfCol])
            setupSd$methGr <- paste0("groups based on mining '",groupPref$gr,"', ie '",tmp$meth,"'")
            setupSd$groups <- wrMisc::rmSharedWords(sdrfDat[,useSdrfCol], sep=c("_","-",",",".","=",";"), callFrom=fxNa, silent=!debug)    # may still contain enumerators
            ## try to find names of groups
            if(length(!duplicated(setupSd$groups)) != length(!duplicated(setupSd$level))) { grNew <- sub("[[:digit:]]+$","", setupSd$groups)       # use rather more elaborate fx ??
              if(length(!duplicated(grNew)) != length(!duplicated(setupSd$level))) setupSd$groups <- grNew else if(!silent) message("Problem: names of groups don't match pettern of levels") }          
            if(debug) message(fxNa," rSM6h.4")
          }
          ## 2.3.2.5  simplify result of mining sdrf
          setupSd$iniGr <- setupSd$groups
          setupSd$groups <- wrMisc::rmSharedWords(setupSd$groups, sep=c("_","-",",",".","=",";"), callFrom=fxNa, silent=!debug)
          ## levels
          setupSd$level <- names(setupSd$groups) <- .adjPat(setupSd$groups)
          names(setupSd$level) <- setupSd$groups
          setupSd$useCol <- useSdrfCol          
          if(debug) message(fxNa," rSM6h.5")

        }      ## finished 2.3.2, determining groups  based on groupPref$gr
        if(debug) {message(fxNa,"length setupSd ", length(setupSd),"  rSM6i"); rSM6i <- list(setupSd=setupSd,sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,suplAnnotFile=suplAnnotFile,groupPref=groupPref,summaryD=summaryD,parametersD=parametersD,setupSdSoft=setupSdSoft,quantMeth=quantMeth) }  #useSdrfCol=useSdrfCol,
 

        ## 2.3.3  consider/integrate info from summaryD if prev options failed
        if(length(setupSd)==0 && length(summaryD) >0 && isTRUE("gr" %in% colnames(summaryD))) {
          setupSd <- list()
          sampNa <- wrMisc::rmSharedWords(summaryD[,"gr"], sep=c("_","-",",",".","=",";"), callFrom=fxNa, silent=!debug)
          if(sum(duplicated(sampNa)) ==0) {
            setupSd$iniSampleName <- summaryD[,"gr"] 
            setupSd$sampleNames <- sampNa 
            setupSd$groups <- wrMisc::rmEnumeratorName(sampNa, nameEnum=c("","Number","No","N","no","number", "#", "Replicate","Rep","Re","R","replicate","rep","re", "Sample","Samp","Sa","S"),
            sepEnum=c(""," ","-","_","/"), incl="rmEnum")  # remove enum
            ## levels
            setupSd$level <- names(setupSd$groups) <- .adjPat(setupSd$groups)
            names(setupSd$level) <- setupSd$groups
          }  
        }
        
        if(length(setupSd$sampleNames) >0) {
          ## check & adjust sampleNames as unique
          chDu <- duplicated(setupSd$sampleNames)
          if(sum(chDu) >0) { setupSd$iniSampleName <- setupSd$sampleNames
            setupSd$sampleNames <- wrMisc::correctToUnique(setupSd$sampleNames, callFrom=fxNa, silent=!debug)}
        }
        if(debug) {message(fxNa,"length setupSd ", length(setupSd),"  rSM6j"); rSM6j <- list(setupSd=setupSd,sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,suplAnnotFile=suplAnnotFile,groupPref=groupPref,summaryD=summaryD,parametersD=parametersD,setupSdSoft=setupSdSoft,quantMeth=quantMeth) }  #useSdrfCol=useSdrfCol,

        ## 2.3.4  check for units & re-adjus setupSd
        grp <- setupSd$groups
        chUnit <- if(isTRUE(chUnit) && all(grepl("[[:punct:]]", grp), na.rm=TRUE)) wrMisc::checkUnitPrefix(grp, if(isTRUE(chUnit)) defUnits else as.character(chUnit), stringentSearch=TRUE, callFrom=fxNa) else NULL
        newNa <- if(length(chUnit) ==1) try(wrMisc::adjustUnitPrefix(grp, unit=chUnit[1], returnType=c("NAifInvalid"), silent=TRUE, callFrom=fxNa), silent=TRUE) else grp
        if(inherits(newNa, "try-error")) {if(!silent) message(fxNa,"Failed to adjust unit-prefixes")} else grp <- newNa
        setupSd$groups <- grp
        setupSd$level <- names(setupSd$groups) <- .adjPat(setupSd$groups)
        names(setupSd$level) <- setupSd$groups
        setupSd$sdrfDat <- sdrfDat
        if(length(iniSdrfOrder) >0) setupSd$iniSdrfOrder <- iniSdrfOrder
      }   ## finished making setupSd
  
#}}}  #  exit after 2.3   ## h

      ## 2.4  ...

    }    ## reading sdrf done

    ## 3.0  allow export of sdrf-draft (based on summaryD, so far only from MQ)
    if(length(parametersSd) >0 && length(setupSd$sdrfDat) <1) {
      setupSd$sdrfExport <- parametersSd
      setupSd$summaryD <- summaryD
    }
    if(debug) { message(fxNa,"rSM8  head of setupSd$level : ",wrMisc::pasteC(utils::head(setupSd$level))); rSM8 <- list(setupSd=setupSd)}
  ## finished readSampleMetaData
  setupSd }


      
       
#
#              } else {   # use specific single method    
#                  tmp <- try(wrMisc::replicateStructure(sdrfDat[,useCol, drop=FALSE], method=groupPref$combMeth[1], silent=silent, callFrom=fxNa, debug=debug), silent=TRUE)
#                  if(inherits(tmp, "try-error")) {
#                    tmp <- try(wrMisc::replicateStructure(sdrfDat[,useCol, drop=FALSE], method="median", silent=silent, callFrom=fxNa, debug=debug), silent=TRUE)
#                    if(!silent) message(fxNa,"Failed determining replicates structure using '",,"' , trying as (default) 'median")
#                  }
#                  if(inherits(tmp, "try-error")) { 
#                    if(!silent) message(fxNa,"Failed determining replicates structure, abandon")
#                    setupSd <- NULL
#                  }
#                }
#              } else {
#                ## simple case
#                
#              }   
#
# 
#                  ## choose among multiple options for grouping (number of groups)
#                  ch1 <- sapply(tmp, function(x) length(x$lev[which(!duplicated(x$lev))]))
#                  lowNumberOfGroups <- if(length(groupPref$lowNumberOfGroups)==1) isTRUE(groupPref$lowNumberOfGroups) else TRUE
#                  useSe <- if(any(ch1 ==1, na.rm=TRUE)) which(ch1 !=1) else if(lowNumberOfGroups) which.min(ch1) else which.max(ch1)
#                  replMeth <- useSe <- useSe[1]
#                  if(!silent) message(fxNa,"Choosing model '",names(useSe),"' for evaluating replicate-structure (ie ",ch1[useSe[1]]," groups of samples)" )         
#                  tmp <- tmp[[useSe]]
#                }   # {message(fxNa,"REMOVING one attempt of understanding replicate-structure") }
#              } else {  
#                ## either just one column to mine or single specific method (eg groupPref$combMeth=="lowest"))
#                if(length(groupPref$combMeth)==1 && nchar(groupPref$combMeth) <1) groupPref$combMeth <- NULL
#                replMeth <- if(isTRUE(groupPref$lowNumberOfGroups) || isTRUE(nchar(groupPref$combMeth) ==1)) "lowest" else groupPref$combMeth
#                tmp <- try(wrMisc::replicateStructure(sdrfDat[,useCol, drop=FALSE], method=replMeth, silent=TRUE, callFrom=fxNa), silent=TRUE)
#                if(inherits(tmp, "try-error")) {syncColumns["sdrfDat"] <- FALSE; if(!silent) message(fxNa,"UNABLE to understand replicate-structure from sdrf (based on method '",replMeth,"')")}
#              }
#
#
#
#
#
#
#            } else {
#              ## no specific options recognized, default use of single col 
#              setupSd$meth <- paste0("groups derived from sdrf-column '",chColNa,"'")
#              setupSd$gr <- sdrfDat[,chColNa] }}
#
#          ## transform sample-names to groups : remove terminal enumerators
#          setupSd$gr <- wrMisc::rmSharedWords(setupSd$gr, sep=c("_","-"," ",".","=",";"), callFrom=fxNa, silent=!debug)   # try making more compact  
#          #setupSd$gr <- wrMisc::rmEnumeratorName(setupSd$gr, sepEnum=c(""," ","-","_"), newSep="_", incl=c("anyCase","trim1"), callFrom="fxNa", silent=silent) 
#          grp <- setupSd$gr <- wrMisc::rmEnumeratorName(setupSd$gr, nameEnum=c("","Number","No","N","no","number", "#", "Replicate","Rep","Re","R","replicate","rep","re", "Sample","Samp","Sa","S"),
#            sepEnum=c(""," ","-","_","/"), incl="rmEnum"))  # remove enum
#
#          ## check for units
#          chUnit <- if(isTRUE(chUnit) && all(grepl("[[:punct:]]", grp), na.rm=TRUE)) wrMisc::checkUnitPrefix(grp, if(isTRUE(chUnit)) defUnits else as.character(chUnit), stringentSearch=TRUE, callFrom=fxNa) else NULL
#          newNa <- if(length(chUnit) ==1) try(wrMisc::adjustUnitPrefix(grp, unit=chUnit[1], returnType=c("NAifInvalid"), silent=TRUE, callFrom=fxNa), silent=TRUE) else grp
#          if(inherits(newNa, "try-error")) {if(!silent) message(fxNa,"Failed to adjust unit-prefixes")} else grp <- newNa
#          setupSd$gr <- grp
#        }
#        
#        ## as groups are set, define levels (sorted integer) & give names
#        setupSd$level <- names(setupSd$gr) <- .adjPat(setupSd$gr) 
#        names(setupSd$level) <- setupSd$gr
#        if(debug) { message(fxNa,"rSM6d2  dim sdrfDat ",nrow(sdrfDat)," ",ncol(sdrfDat)); rSM6d2 <- list(sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,groupPref=groupPref,summaryD=summaryD,parametersD=parametersD,syncColumns=syncColumns,iniSdrfOrder=iniSdrfOrder,setupSd=setupSd, grp=grp,chUnit=chUnit,newNa=newNa,chColNa=chColNa) }
#        
#       
#      }    ## finish setting sample-names and groups from sdrf (summaryD not yet considered)
#
#
#
#
#        if(length(groupPref$gr)==0 || (length(groupPref$gr)==1 && grepl("^sdrf", groupPref$gr[1]))) {     
#          ## check for custom-provided gr  (priority)   wr modif 25sep24
#          ## check for picking specific column of sdrf (eg 'sdrf$comment.disease')
#          if(grepl("^sdrf\\$[[:alpha:]]", groupPref$gr[1])) {                                         
#            chColNa <- sub("^sdrf\\$","", groupPref$gr[1])
#            if(chColNa %in% colnames(sdrfDat)) { 
#              ## groupPref$gr points to specific column of sdrf
#              
#              grp <- wrMisc::rmSharedWords(sdrfDat[,chColNa], sep=c("_","-"," ",".","=",";"), callFrom=fxNa, silent=!debug)   # try making more compact
#              chUnit <- if(isTRUE(chUnit) && all(grepl("[[:punct:]]", grp), na.rm=TRUE)) wrMisc::checkUnitPrefix(grp, if(isTRUE(chUnit)) defUnits else as.character(chUnit), stringentSearch=TRUE) else NULL
#              newNa <- if(length(chUnit) ==1) try(wrMisc::adjustUnitPrefix(grp, unit=chUnit[1], returnType=c("NAifInvalid"), silent=TRUE, callFrom=fxNa), silent=TRUE) else grp
#              if(inherits(newNa, "try-error")) {if(!silent) message(fxNa,"Failed to adjust unit-prefixes")} else grp <- newNa
#              setupSd$gr <- grp
#              setupSd$level <- .adjPat(grp) 
#              setupSd$meth <- paste0("custom groups from sdrf-column '",chColNa,"'")
#              if(debug) { message(fxNa,"rSM6d1  dim sdrfDat ",nrow(sdrfDat)," ",ncol(sdrfDat)); rSM6d1 <- list(sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,groupPref=groupPref,summaryD=summaryD,parametersD=parametersD,syncColumns=syncColumns,iniSdrfOrder=iniSdrfOrder,setupSd=setupSd, grp=grp,chUnit=chUnit,newNa=newNa,chColNa=chColNa) }
#            }          
#          } else {
#            ## regular minig of sdrf
#            groupPref$gr <- "sdrf"
#          }
#        } else  {    # no custom provided  groupPref$gr
#          if(length(abund) <1) message(fxNa,"Can't check if length of $gr is OK since 'abund' empty or not given ..")
#          setupSd$gr <- groupPref$gr
#          setupSd$level <- grp <- .adjPat(setupSd$gr) 
#          setupSd$meth <- paste0("custom from groupPref$gr ")
#          ## Note :  groupPref$gr may ny absent =>  setupSd empty      
#        }  
#        if(debug) { message(fxNa,"rSM6d2  dim sdrfDat ",nrow(sdrfDat)," ",ncol(sdrfDat)); rSM6d2 <- list(sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,groupPref=groupPref,summaryD=summaryD,parametersD=parametersD,syncColumns=syncColumns,iniSdrfOrder=iniSdrfOrder,setupSd=setupSd,chUnit=chUnit) }
#
#        ## DEFAULT mining (if $levels and/or $sampleNames not yet set)
#        #if(length(setupSd$level) != ncol(abund) || length(setupSd$sampleNames) != ncol(abund)) {
#        if(length(setupSd$sampleNames) <1 || length(setupSd$level) <1) {
#          ## check for column 'comment.technical.replicate.' (to exclude from using)
#          ## look for and/or 'characteristics.spiked.compound.' ?
#          useSdrfCol <- 1:ncol(sdrfDat)  
#          chTechnReplCol <- colnames(sdrfDat) %in% "comment.technical.replicate."
#          if(any(chTechnReplCol, na.rm=TRUE)) useSdrfCol <- useSdrfCol[-which(chTechnReplCol)]    # exclude designation of technical replicates
#  
#          ## option : use sampleNames from sdrf-fileNames
#          if("sdrf" %in% groupPref$sampleNames && length(groupPref$sampleNames)==1) {
#            ## use 1st of 'comment.data.file.' and 'comment.file.uri.'
#            colNa <- c("comment.data.file.", "comment.file.uri.")
#            chCol <- colNa %in% colnames(sdrfDat)
#            if(any(chCol)) { 
#              sampleNames <- sdrfDat[,which(chCol[1])]        # as in order of abund
#              ## try to simplify (remove redundant words, try adjusting varying units)
#              iniNa <- wrMisc::rmSharedWords(sampleNames, sep=c(" ","_","-","/"), silent=silent, debug=debug, callFrom=fxNa)
#              if(!isFALSE(chUnit)) chUnit <- wrMisc::checkUnitPrefix(iniNa, unit=if(isTRUE(chUnit)) defUnits else as.character(chUnit))
#              if(length(chUnit) >0) {
#                newNa <- try(wrMisc::adjustUnitPrefix(iniNa, unit=chUnit[1], returnType=c("NAifInvalid"), silent=silent, debug=debug, callFrom=fxNa), silent=TRUE)
#                if(inherits(newNa, "try-error")) {if(!silent) message(fxNa,"Failed to adjust unit-prefixes")} else setupSd$sampleNames <- groupPref$sampleNames <- newNa           
#              }
#              ## need adjust to summaryD ??
#            } else groupPref <- groupPref[-1*which(names(groupPref) ==sampleNames)] 
#          } else {
#            ## sampleNames not yet set, need to pick sampleNames from other ressources !      
#          }         
#        }  
#        if(debug) { message(fxNa,"rSM6e  dim sdrfDat ",nrow(sdrfDat)," ",ncol(sdrfDat)); rSM6e <- list(sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,groupPref=groupPref,summaryD=summaryD,parametersD=parametersD,syncColumns=syncColumns,iniSdrfOrder=iniSdrfOrder,setupSd=setupSd) }
#      
#        setupSd$sdrfDat <- sdrfDat   # useful here ?
#        
#          
#        ## now determine gr  (if not yet present)
#        if(length(setupSd$gr) < (if(length(abund) >1) ncol(abund) else 1)) {
#        #old# if(length(setupSd$gr) <1 || identical(groupPref$gr,"sdrf")) {               #  && length(dim(sdrf)) >1
#          replStrOpt <- c("highest","lowest","min","max","median","combAll","combNonOrth")
#          if(debug) {message(fxNa,if(debug)"rSM6g  ","length setupSd ", length(setupSd)); rSM6g <- list(setupSd=setupSd,sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,suplAnnotFile=suplAnnotFile,groupPref=groupPref,summaryD=summaryD,parametersD=parametersD,setupSdSoft=setupSdSoft,quantMeth=quantMeth) }  #useSdrfCol=useSdrfCol,
#          ## check for custom provided method for sdrf-mining : (is it risky to search in 2nd value of sdrf ?)
#          ## 'groupPref$useCol' design custom choice of multiple columns 
#          ## 'groupPref$combMeth' designs method for choosing from automatic mining            
#          if("useCol" %in% names(groupPref) && all(nchar(groupPref$useCol) >0)) {   ## custom choice for sdrf-column
#            chCol <- which(colnames(sdrfDat) %in% groupPref$useCol)
#          } else chCol <- which(!colnames(sdrfDat) %in% c("comment.technical.replicate."))            
#          if(length(chCol) >0) {  # chCol <- sdrfDat[, useSdrfCol[which(chCol)]]
#            setupSdIni <- setupSd
#            replMeth <- "failed"   # initialize
#            ## check content of 'combMeth' ?
#            #old#if(!isTRUE(groupPref$lowNumberOfGroups) &&  isTRUE(groupPref$combMeth != "lowest")) {
#            if(length(chCol) >1 && isTRUE(groupPref$combMeth != "lowest")) {
#              ## test 2 methods
#              tmp <- list(combNonOrth=try(wrMisc::replicateStructure(sdrfDat[,chCol], method=groupPref$combMeth, silent=silent, callFrom=fxNa, debug=debug), silent=TRUE),
#                lowest=try(wrMisc::replicateStructure(sdrfDat[,chCol], method="lowest", silent=silent, callFrom=fxNa, debug=debug), silent=TRUE))  
#              ch1 <- sapply(tmp, inherits, "try-error")
#              if(all(ch1)) { message(fxNa,"UNABLE to understand replicate-structure from sdrf !!"); setupSd <- NULL; syncColumns["sdrfDat"] <- FALSE
#              } else {
#                ## choose among multiple options for grouping (number of groups)
#                ch1 <- sapply(tmp, function(x) length(x$lev[which(!duplicated(x$lev))]))
#                lowNumberOfGroups <- if(length(groupPref$lowNumberOfGroups)==1) isTRUE(groupPref$lowNumberOfGroups) else TRUE
#                useSe <- if(any(ch1 ==1, na.rm=TRUE)) which(ch1 !=1) else if(lowNumberOfGroups) which.min(ch1) else which.max(ch1)
#                replMeth <- useSe <- useSe[1]
#                if(!silent) message(fxNa,"Choosing model '",names(useSe),"' for evaluating replicate-structure (ie ",ch1[useSe[1]]," groups of samples)" )         
#                tmp <- tmp[[useSe]]
#              }   # {message(fxNa,"REMOVING one attempt of understanding replicate-structure") }
#            } else {  
#              ## either just one column to mine or single specific method (eg groupPref$combMeth=="lowest"))
#              if(length(groupPref$combMeth)==1 && nchar(groupPref$combMeth) <1) groupPref$combMeth <- NULL
#              replMeth <- if(isTRUE(groupPref$lowNumberOfGroups) || isTRUE(nchar(groupPref$combMeth) ==1)) "lowest" else groupPref$combMeth
#              tmp <- try(wrMisc::replicateStructure(sdrfDat[,chCol], method=replMeth, silent=TRUE, callFrom=fxNa), silent=TRUE)
#              if(inherits(tmp, "try-error")) {syncColumns["sdrfDat"] <- FALSE; if(!silent) message(fxNa,"UNABLE to understand replicate-structure from sdrf (based on method '",replMeth,"')")}
#            }
#            if(debug) {message(fxNa,if(debug)"rSM6g2  ","length setupSd ", length(setupSd)); rSM6g2 <- list(setupSd=setupSd,setupSdIni=setupSdIni,sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,suplAnnotFile=suplAnnotFile,groupPref=groupPref,summaryD=summaryD,parametersD=parametersD,setupSdSoft=setupSdSoft,quantMeth=quantMeth,tmp=tmp) }  #useSdrfCol=useSdrfCol,
#               
#            if(!inherits(tmp, "try-error") && length(tmp) >0) {   # valid mining
#              if(isTRUE(tmp$col==1) && all(grepl("^Sample {0,1}[[:digit:]]+$", names(tmp$lev))) ) {  # the column with sample-numbers was picked, not very informative; try to find better one
#                ## Try finding better content/text for names of levels, ie tmp$lev
#                chPa <- apply(sdrfDat[,chCol], 2, .adjPat)
#                ch1 <- apply(chPa[,-1], 2, function(x) all(chPa[,1]==x, na.rm=TRUE))
#                if(any(ch1)) names(tmp$lev) <- wrMisc::rmSharedWords(sdrfDat[,chCol[which(ch1)[1] +1]], sep=c("_"," ",".","=",";"), silent=silent,callFrom=fxNa) else {
#                  ## otherwise try to find column/pattern fitting to prev hit after removing redundant text & removing enumerators
#                  ref2 <- apply(sdrfDat[,chCol[-1]], 2, wrMisc::rmSharedWords, sep=c("_"," ",".","=",";"), silent=silent,callFrom=fxNa)
#                  ref2 <- apply(ref2, 2, wrMisc::rmEnumeratorName, silent=TRUE)   
#                  chPa2 <- apply(ref2, 2, .adjPat)
#                  ch1 <- apply(chPa2, 2, function(x) all(chPa[,1]==x, na.rm=TRUE))
#                  if(any(ch1)) {names(tmp$lev) <- chPa2[,which(ch1)[1]] }  # also document column used ?
#                }                 
#              }
#
#              if(length(abund) >0 && length(setupSdIni$sampleNames) == ncol(abund)) setupSd$sampleNames <- setupSdIni$sampleNames     # priority to specified $sampleNames
#              setupSd$level <- tmp$lev                # grouping correct,   names may be not meaningful (if column picked contains 'Sample 1' etc)
#              if(debug) {message(fxNa,if(debug)"rSM6g3  ","length tmp$lev ", length(tmp$lev)); rSM6g3 <- list(setupSd=setupSd,sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,suplAnnotFile=suplAnnotFile,
#                groupPref=groupPref,summaryD=summaryD,parametersD=parametersD,setupSdSoft=setupSdSoft,quantMeth=quantMeth,tmp=tmp,chUnit=chUnit) }
#
#              ## optional adjusting of units to readjust levels
#              grp <- wrMisc::rmSharedWords(tmp$lev, sep=c("_"," ",".","=",";"))            # try making more compact
#              if(!isFALSE(chUnit[1])) { 
#                if(all(grepl("[[:punct:]]", grp), na.rm=TRUE) && (length(chUnit) <1 || isTRUE(chUnit))) chUnit <- wrMisc::checkUnitPrefix(grp, if(isTRUE(chUnit)) defUnits else as.character(chUnit), stringentSearch=TRUE)
#                if(length(chUnit) ==1 && nchar(chUnit) >0) {newNa <- try(wrMisc::adjustUnitPrefix(grp, unit=chUnit[1], returnType=c("NAifInvalid"), silent=TRUE, callFrom=fxNa), silent=TRUE)
#                  if(inherits(newNa, "try-error")) {if(!silent) message(fxNa,"Failed to adjust unit-prefixes")} else setupSd$level <- grp <- .adjPat(newNa) }} 
#              setupSd$col <- tmp$col  
#              setupSd$meth <- replMeth
#            }
#          } 
#          if(debug) {message(fxNa,if(debug)"rSM6h  ","length setupSd ", length(setupSd)); rSM6h <- list(setupSd=setupSd,sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,suplAnnotFile=suplAnnotFile,
#            groupPref=groupPref,summaryD=summaryD,parametersD=parametersD,setupSdSoft=setupSdSoft,quantMeth=quantMeth,tmp=tmp,chUnit=chUnit) }
#          
#
#          ## remove/rename $lev (normally not expected any more)          
#          chLe <- names(setupSd) %in% "lev"
#          if(any(chLe)) { 
#            if(debug) message(fxNa,"Still found setupSd$lev !! Used to replace setupSd$level")
#            names(setupSd$lev) <- gsub("^[[:space:]]*","", names(setupSd$lev))
#            setupSd$level <- setupSd$lev
#            setupSd <- setupSd[-which(chLe)]
#          }
#          if(debug) {message(fxNa,if(debug) "rSM6i  ","length setupSd ", length(setupSd)); rSM6i <- list() }
#        }   # finish determining setupSd$gr 
#                 
#        if(debug) {message(fxNa,"rSM6j"); rSM6j <- list(setupSd=setupSd,sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,suplAnnotFile=suplAnnotFile,groupPref=groupPref,summaryD=summaryD,parametersD=parametersD,setupSdSoft=setupSdSoft,quantMeth=quantMeth)}
#        if("setupSd" %in% names(setupSd)) { setupSd <- wrMisc::partUnlist(setupSd, callFrom=fxNa,debug=debug);
#          if(debug) message(fxNa,"rSM6j2  - not expecting list of list(s) for setupSd ! .. correcting")}
#          
#      
#        if(debug) {message(fxNa,"rSM6i   names setupSd : ", wrMisc::pasteC(names(setupSd))); rSM6g <- list() }
#        if(!is.list(setupSd)) {setupSd <- as.list(setupSd); if(debug) message(fxNa,"rSM6i  'setupSd' should be list, but was NOT !!")}
#        if(!"sdrfDat" %in% names(setupSd)) setupSd$sdrfDat <- sdrfDat
#        if(debug) {message(fxNa, "rSM6k  ")
#          rSM6k <- list(sdrf=sdrf,setupSd=setupSd,sdrfDat=sdrfDat,abund=abund,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,groupPref=groupPref,newNa=newNa,summaryD=summaryD,parametersD=parametersD,syncColumns=syncColumns,iniSdrfOrder=iniSdrfOrder)
#        }
#
#        setupSd$annotBySoft <- as.data.frame(summaryD)
#        setupSd$syncColumns <- syncColumns
#
#      } else {      ## sdrf was given - but NOT conform : (no soft-generated sample annot available) try to match colnames of abund
#        if(debug) message(fxNa, if(debug) "rSM6l ","NO valid sdrf found")
#        ## ie single source of info
#        if(length(summaryD) <1) {                                              ## ie no  summaryD
#          if(debug) message(fxNa, if(debug) "rSM6m ","NO valid sdrf and NO valid information (summaryD) from quant-software found")
#        } else {                         # ie summaryD is available
#          setupSd <- setupSdSoft
#          if(!silent) message(fxNa, if(debug) "rSM6n ","Reading of sdrf was NOT successful and no summaryD available => nothing can be done to mine experimental setup...")
#        }
#      }
#            
#      if(length(iniSdrfOrder) >0) setupSd$iniSdrfOrder <- iniSdrfOrder
#      if(debug) { message(fxNa,"rSM7  head of setupSd$level : ",wrMisc::pasteC(utils::head(setupSd$level))); rSM7 <- list(setupSd=setupSd,sdrf=sdrf,sdrfDat=sdrfDat,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,summaryD=summaryD,nSamp0=nSamp0,iniSdrfOrder=iniSdrfOrder)}
#      if(length(setupSd) >0) if(length(setupSd$level) != nSamp0 && length(abund) >0) {               ## keep this ? - redundant !
#        if(!silent) warning(fxNa, if(debug) "rSM7  ","Invalid information from sample meta-data or wrong experiment ! Number of samples from sdrf ",
#          " (",length(setupSd$level),") and from experimental data (",ncol(abund),") don't match !")
#        setupSd <- NULL } else {
#          if(length(abund) <1 && !silent) message(fxNa,"Note: Order of lines in sdrf not ajusted since no valid 'abund' given...")
#      }
#
#    } else { setupSd <- setupSdSoft; setupSd$annotBySoft <- summaryD }
#    
#    ## allow export of sdrf-draft (so far only from MQ)
#    if(length(parametersSd) >0 && length(setupSd$sdrfDat) <1) {
#      setupSd$sdrfExport <- parametersSd
#      setupSd$summaryD <- summaryD
#      }
#
#    if(debug) { message(fxNa,"rSM8  head of setupSd$level : ",wrMisc::pasteC(utils::head(setupSd$level))); rSM8 <- list(setupSd=setupSd)}
#
#  }
#  ## finished readSampleMetaData
#  setupSd }
