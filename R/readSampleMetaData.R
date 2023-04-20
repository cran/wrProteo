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
#'  if character, this may be the ID at ProteomeExchange or a similarly formatted local file. \code{sdrf} will get priority over \code{suplAnnotFile}, if provided.
#' @param suplAnnotFile (logical or character) optional reading of supplemental files produced by MaxQuant; if \code{gr} is provided, it gets priority for grouping of replicates
#'  if \code{TRUE} in case of \code{method=="MQ"} (MaxQuant) default to files 'summary.txt' (needed to match information of \code{sdrf}) and 'parameters.txt' which can be found in the same folder as the main quantitation results;
#'  if \code{character} the respective file-names (relative ro absolute path), 1st is expected to correspond to 'summary.txt' (tabulated text, the samples as given to MaxQuant) and 2nd to 'parameters.txt' (tabulated text, all parameters given to MaxQuant)
#'  in case of \code{method=="PL"} (Proline), this argument should contain the initial file-name (for the identification and quantification data) in the first position
#' @param path (character) optional path of file(s) to be read
#' @param abund (matrix or data.frame) experimental quantitation data; only column-names will be used for aligning order of annotated samples
#' @param groupPref (list) additional parameters for interpreting meta-data to identify structure of groups (replicates);
#'   May contain \code{lowNumberOfGroups=FALSE} for automatically choosing a rather elevated number of groups if possible (defaults to low number of groups, ie higher number of samples per group)
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a list with \code{$lev} and \code{$level} (grouping of samples given as integer), and \code{$meth} (method by which grouping as determined).
#'  If valid \code{sdrf} was given, the resultant list contains in addition \code{$sdrfDat} (data.frame of annotation). If software annotation has been found it will be shown in \code{$annotBySoft}.
#'  If all entries are invalid or entries do not pass the tests, this functions returns an empty \code{list}.
#' @seealso this function is used internally by \code{\link{readMaxQuantFile}}, \code{\link{readProtDiscovFile}} etc; use \code{\link{readSdrf}} for reading sdrf-files, \code{\link[wrMisc]{replicateStructure}} for mining annotation columns
#' @examples
#' sdrf001819Setup <- readSampleMetaData(quantMeth=NA, sdrf="PXD001819")
#' str(sdrf001819Setup)
#'
#' @export
readSampleMetaData <- function(quantMeth, sdrf=NULL, suplAnnotFile=NULL, path=".", abund=NULL, groupPref=list(lowNumberOfGroups=TRUE), silent=FALSE, debug=FALSE, callFrom=NULL)  {
  ##  sdrf..()
  ##  suplAnnotFile..(character or logical)
  ##  quantMeth..(character)
  ##  abund..(matrix or data.frame) column-names will be used to comapre & align sample meta-data)

  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readSampleMetaData")
  if(isTRUE(debug)) silent <- FALSE
  if(!isTRUE(silent)) silent <- FALSE
  summaryD <- parametersD <- setupSdSoft <- setupSd <- sdrfInf <- annSh <- NULL         # initialize    (setupSd needed ?)
  ## checks
  if(length(suplAnnotFile) >1) if(is.na(suplAnnotFile[1])) suplAnnotFile <- NULL
  datOK <- length(sdrf) >0 || length(suplAnnotFile) >0
  if(length(quantMeth) <1) quantMeth <- NA
  if(length(abund) >0 && any(length(dim(abund)) !=2, dim(abund) < 1, na.rm=TRUE)) { datOK <- FALSE
    warning("Invalid argument 'abund'; must be matrix or data.frame with min 1 line and 1 col")}
  if(debug) {message(fxNa,"Ready search & extract sample meta-data   rSM0"); rSM0 <- list(sdrf=sdrf,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,path=path,abund=abund)}


  .corPathW <- function(x) gsub("\\\\", "/", x)
  .adjPat <- function(x) { out <- match(x, unique(x)); names(out) <- names(x); out}  # needed ??
  .adjTrimPat <- function(x) { x <- wrMisc::rmEnumeratorName(x, incl=c("anyCase","trim0","rmEnum"), sepEnum=c(" ","-","_"), nameEnum=c("Number","No","#","","Replicate","Sample"))
    out <- match(x, unique(x)); names(out) <- names(x); out}   # used

  .redLstToDf <- function(lst) {    # transform lst to data.frame; in case some list-entries have different length, choose the entries with most feq common length
    leL <- sapply(lst, length)
    if(any(duplicated(leL))) {      # need to reduce : find most frequent
      leL2 <- tabulate(leL)
      lst <- lst[which(leL==which.max(leL2))] }
    as.data.frame(lst) }
  ## end suppl fx


  path <- if(length(path) <1) "." else path[1]
  nSamp0 <- if(length(dim(abund)) >1) ncol(abund) else 0
  chSoft <- c("MQ", "PD", "PL", "FP", "MC")
  syncColumns <- c(sdrfDat=NA, annotBySoft=NA)
  if(datOK) {
    if("maxquant" %in% tolower(quantMeth)) quantMeth <- "MQ"
    if("proteomediscoverer" %in% tolower(quantMeth)) quantMeth <- "PD"
    if("proline" %in% tolower(quantMeth)) quantMeth <- "PL"
    if("fragpipe" %in% tolower(quantMeth)) quantMeth <- "FP"
    if("masschroq" %in% tolower(quantMeth)) quantMeth <- "MC"
    }

  if(datOK) {
    if(length(abund) >0) if(is.null(colnames(abund))) { abund <- NULL
      if(!silent) message(fxNa,"Invalid 'abund' : has NO colnames !") }

    ## priority column for groups from sdrf (+ define default colnames for priority)
    chGr <- c("sdrfColumn","sdrfCol")
    chGrPref <- chGr %in% names(groupPref)
    if(any(chGrPref)) {
      groupPref$sdrfColumn <- if(length(groupPref[[chGr[which(chGrPref)[1]]]]) >0) groupPref[[chGr[which(chGrPref)[1]]]] else c("factor.value.disease.","characteristics.disease.", "factor.value.treatment.","characteristics.treatment.")
    }

    ### IMPORT SAMPLE META-DATA, if possible GROUPING OF REPLICATES
    if(length(suplAnnotFile) ==1) {
      if(isFALSE(suplAnnotFile)) suplAnnotFile <- NULL else if(is.na(suplAnnotFile)) suplAnnotFile <- NULL }

    ## 1.1 SOFTWARE specific META-DATA : read additional annotation & documentation files produced by var software as  summaryD & parametersD
    if(length(suplAnnotFile) >0)  {     # read quant software-generated sample annotation
      chFiNa <- NULL                    # initialize
      if(debug) {message(fxNa,"rSM1"); rSM1 <- list(sdrf=sdrf,abund=abund,path=path,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth) }
           ## option 1 : suplAnnotFile has path (do not use default 'path'), use same path for default suplAnnotFile (if applicable)
           ## option 2 : suplAnnotFile has no path, use 'path' for sdrf & suplAnnotFile

      ## Aim : extract/build 'summaryD' allowing to match colnames of 'abund' to suplAnnotFile and/or sdrf

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
            chFi <- c(chFi[1] | chFi[2], chFi[3] | chFi[4])    #needed ?
          } else suplAnnotFile <- NULL
        } else {      # specific/non-default file given
          if(length(suplAnnotFile) >2) suplAnnotFile <- suplAnnotFile[1:2]   # use max length=2
            chFi <- rep(FALSE, 2)
    	    if(!is.na(suplAnnotFile[1])) chFi[1] <- file.exists(file.path(path, suplAnnotFile[1]))
    	    if(!is.na(suplAnnotFile[2])) chFi[2] <- file.exists(file.path(path, suplAnnotFile[2]))
        }
        if(debug) {message(fxNa,"rSM1mq"); rSM1mq <- list(path=path,sdrf=sdrf,summaryD=summaryD,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,path=path,nSamp0=nSamp0,chFiNa=chFiNa,chFi=chFi )}

        ## main reading of MQ sample meta-data
        if(chFi[1]) summaryD <- try(utils::read.delim(file.path(path, suplAnnotFile[1]), stringsAsFactors=FALSE), silent=TRUE)
        if(chFi[2]) parametersD <- try(utils::read.delim(file.path(path, suplAnnotFile[2]), stringsAsFactors=FALSE), silent=TRUE)
        if(inherits(summaryD, "try-error")) {summaryD <- NULL; if(!silent) message(fxNa,"Meta-data: Failed to read '",suplAnnotFile[1],"'  for getting additional information about experiment !")} else {
          summaryD <- if(nrow(summaryD) >2) summaryD[-nrow(summaryD),] else matrix(summaryD[-nrow(summaryD),], nrow=1,dimnames=list(NULL,colnames(summaryD)))  # need to remove last summary-line
          if(debug) message(fxNa,"Successfully read sample annotation from '",suplAnnotFile[1],"'") }
        if(inherits(parametersD, "try-error")) {if(!silent) message(fxNa,"Meta-data: Failed to read '",suplAnnotFile[2],"' !")} else {
          if(debug && chFi[2]) message(fxNa,"Successfully read ",quantMeth," parameters from '",suplAnnotFile[2],"'") }
        syncColumns["annotBySoft"] <- FALSE
        if(debug) { message(fxNa,"rSM1mq2")}
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
        if(inherits(chFi, "try-error") & silent) {chFi <- FALSE; message(fxNa,"Meta-data: Failed to see file '",suplAnnotFile[1]," ! (check if file exists or rights to read directory ?)")}
        if(debug) {message(fxNa,"rSM1pd2") }

        ## main reading of PD sample meta-data
        if(chFi) summaryD <- try(utils::read.delim(file.path(path, suplAnnotFile[1]), stringsAsFactors=FALSE), silent=TRUE)
        if(inherits(summaryD, "try-error")) {summaryD <- NULL; if(!silent) message(fxNa,"Meta-data: Failed to read '",suplAnnotFile[1],"' !")
        } else {
          syncColumns["annotBySoft"] <- FALSE
          if(debug)  message(fxNa,"ProteomeDiscoverer Meta-data successfully read '",suplAnnotFile[1])}
        if(debug) {message(fxNa,"rSM1pd3")}
      }

      ## Proline
      ## so far only for reading out of xslx
      if("PL" %in% quantMeth && length(suplAnnotFile) >0) {
        if(debug) {message(fxNa,"rSM0pl"); rSM0pl <- list(sdrf=sdrf,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth)}
        summaryD <- NULL
        ## need init filename given via suplAnnotFile
        #if(isTRUE(suplAnnotFile)) { tmp <- grep("\\.xlsx$", dir(path))
        #  suplAnnotFile <- if(length(tmp) >0) dir(path, full.names=TRUE)[grep("\\.xlsx$", dir(path))[1]] else NULL }
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

        ## FragPipe :    (filelist_ionquant.txt & fragpipe-files.fp-manifest)
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

      if("MC" %in% quantMeth && length(suplAnnotFile) >0) {
        warning(fxNa,"Reading supplemental meta-data from MassChroq is currently not implemented") }

      ## OTHER software ? ..
      if(!any(quantMeth %in% chSoft, !silent, na.rm=TRUE)) message(fxNa,"Note: No specific procedure has been implemented so far for gathering meta-data by the analysis-software/method '",quantMeth,"'")
    }            ## finished main reading of suplAnnotFile into summaryD
    if(debug) { message(fxNa,"rSM2"); rSM2 <- list(sdrf=sdrf,abund=abund,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,summaryD=summaryD,suplAnnotFile=suplAnnotFile,syncColumns=syncColumns) }


    ## 1.2 basic check of summaryD to quant data
    if(length(summaryD) >0) {      ##  more checks
      if(length(abund) <1) message(fxNa,"Can't verify/correct names of annotation since content of 'abund' has was not given (ie NULL) or has no colnames") else {
        if(!identical(ncol(abund), nrow(summaryD))) { summaryD <- NULL
		      if(!silent) message(fxNa,"Note : Number of columns of 'abund' does NOT FIT to number of samples in annotation-data !")	}
	    }
      if(length(dim(summaryD)) !=2) summaryD <- matrix(summaryD, ncol=1, dimnames=list(names(summaryD),NULL))
	  }
    if(debug) { message(fxNa,"rSM3"); rSM3 <- list(sdrf=sdrf,abund=abund,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,summaryD=summaryD,suplAnnotFile=suplAnnotFile,syncColumns=syncColumns) }

    ## continue evaluating summaryD to consistent format
    if(length(summaryD) >0) {      ## define  setupSdSoft
      ## need to match colnames(abund) to (MQ:) $Raw.file or $Experiment  .. need to find best partial match
      if("MQ" %in% quantMeth) {         ## NOT IN SAME ORDER !!
        summaryD <- summaryD[,wrMisc::naOmit(match(c("Raw.file","Experiment"), colnames(summaryD)))]           # cor 21oct22
        chSd <- length(abund) >0 && nrow(summaryD) == ncol(abund)
        if(length(chSd) <1) chSd <- FALSE
        ## normally  colnames(abund) and summaryD should alread be in correct order
        if(!chSd) {
          if(!silent && length(abund) >0) if(nrow(summaryD) == ncol(abund)) { message(fxNa,"PROBLEM : meta-data and abundance data do not match !  ",
             "Number of samples from ",suplAnnotFile[1]," (",nrow(summaryD),") and from main data (",ncol(abund),") do NOT match !! .. ignoring") }
            if(debug) save(sdrf,abund,suplAnnotFile,quantMeth,summaryD,quantMeth,syncColumns, file="C:\\E\\projects\\TCAmethods\\wrProteoRamus\\rSM4mq.RData") }
        #syncColumns["annotBySoft"] <- TRUE

        if(debug) { message(fxNa," .. rSM4mq"); rSM4mq <- list()}
      }

      if("PD" %in% quantMeth) { useCo <- c("Input.Files.","File.ID","File.Name","Instrument.Name")    # no suitable 2nd column ...
        useCo <- wrMisc::naOmit(match(useCo, colnames(summaryD)))
        summaryD <- if(length(useCo) >1) summaryD[,useCo] else matrix(summaryD, ncol=1, dimnames=list(rownames(summaryD), colnames(summaryD)[useCo]))
        ## presume that filenames (from summaryD) are in same order as abund, then trim to file-names (if all in same path)
        if(debug) { message(fxNa,"rSM4pd"); rSM4pd <- list(sdrf=sdrf,useCo=useCo,abund=abund,uplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,summaryD=summaryD,chFiNa=chFiNa) }
        colNa <- wrMisc::trimRedundText(gsub("\\\\","/",as.character(summaryD[,"File.Name"])), silent=silent, debug=debug, callFrom=fxNa)
        if(length(colNa) < ncol(abund)) warning(fxNa,"Trouble ahead : Sample annotation data from ProteomeDiscoverer has FEWER samples than data read !") else {
          if(length(colNa) > ncol(abund)) { message(fxNa,"note : Sample annotation data from ProteomeDiscoverer has MORE samples than data read, using only first (might be incorrect)")
            colNa <- colNa[1:ncol(abund)]
            summaryD <- summaryD[1:ncol(abund),]
            } }
        colnames(abund) <- colNa                    # no possibility to match colnames at this point
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
          if(debug) { message(fxNa,"rSM4pl"); rSM4pl <- list(sdrf=sdrf,abund=abund,uplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,summaryD=summaryD)
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
          summaryD[,4] <- gsub("_+$|-+$|\\.+$| +$|","", sub("[[:digit:]]+$","", wrMisc::trimRedundText(summaryD[,3], side="right", callFrom=fxNa, silent=silent)))   # remove tailing numbers (and tailing redundant text to get to numbers)
          chSd <- length(abund) >0 && nrow(summaryD) == ncol(abund)
          if(length(chSd) <1) chSd <- FALSE
          ## normally  colnames(abund) and summaryD should alread be in same/correct order
          if(!chSd) {
            if(!silent && length(abund) >0) if(nrow(summaryD) == ncol(abund)) { message(fxNa,"PROBLEM : meta-data and abundance data do not match !  ",
               "Number of samples from ",suplAnnotFile[1]," (",nrow(summaryD),") and from main data (",ncol(abund),") do NOT match !! .. ignoring") }
              if(debug) save(sdrf,abund,suplAnnotFile,quantMeth,summaryD,quantMeth, file="C:\\E\\projects\\TCAmethods\\wrProteoRamus\\rSM4mq.RData") }
        }
        syncColumns["annotBySoft"] <- TRUE
        if(debug) { message(fxNa," .. rSM4mq"); rSM4fp <- list()}
      }

      ## other software ? ...

      if(debug) { message(fxNa,"rSM4d"); rSM4d <- list(sdrf=sdrf,abund=abund,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,summaryD=summaryD,syncColumns=syncColumns) }


      ## 1.3 TRY CHECKING/ADJUSTING ORDER of summaryD
      if(length(abund) >0 && length(summaryD) >0) {
        ## some software specific options otherwise check if filenames can be matched to colnames ?
        ## PD not much possible since colnames  ".F1.Sample",".F2.Sample",".F3.Sample",...
        ## most other software has summaryD in same order as abund
        if("MQ" %in% quantMeth) {   # colnames of abund not necessarly found in summaryD
          summaryD <- wrMisc::matchMatrixLinesToRef(mat=summaryD, ref=colnames(abund), inclInfo=TRUE, silent=silent, debug=debug, callFrom=fxNa)
          syncColumns["annotBySoft"] <- length(summaryD$newOrder)  >0
          summaryD <- summaryD$mat  }
        if(any(c("PL","FP") %in% quantMeth, na.rm=TRUE)) {
          summaryD <- wrMisc::matchMatrixLinesToRef(mat=summaryD, ref=colnames(abund), inclInfo=TRUE, silent=silent, debug=debug, callFrom=fxNa)
          syncColumns["annotBySoft"] <- length(summaryD$newOrder)  >0
          summaryD <- summaryD$mat  }
      }



      ## 1.4   replicateStructure
      setupSdSoft <- wrMisc::replicateStructure(summaryD, silent=silent, debug=debug, callFrom=fxNa)
      if(debug) { message(fxNa,"rSM4e"); rSM4e <- list() }

      ## so far no direct information about groups (all filenames are different), need to try to find out (remove enumerators)
      if(length(abund) >0) {
        grpA <- wrMisc::trimRedundText(txt=colnames(abund), spaceElim=TRUE, silent=silent, debug=debug, callFrom=fxNa)                 # 26oct22
        colNaGrp <- wrMisc::rmEnumeratorName(grpA, incl=c("anyCase","trim0","rmEnum"), sepEnum=c(" ","-","_"), nameEnum=c("Number","No","#","","Replicate","Sample"), silent=silent, debug=debug, callFrom=fxNa)
        colNaGrPref <- TRUE                                  ## preferential to colnames for searching groups
        if(any(duplicated(colNaGrp)) && colNaGrPref) {       # colnames may be used for designing groups
           grp <- .adjPat(colNaGrp)
           names(grp) <- colNaGrp
           setupSdSoft$lev <- grp
        } else {
          if(all(setupSdSoft$lev ==1:ncol(abund))) {
            ## note : .adjTrimPat() does NOT allow keeping names of levels
            grp2 <- if(ncol(summaryD) >1) apply(summaryD, 2, .adjTrimPat) else as.matrix(.adjTrimPat(summaryD))
            if(ncol(summaryD) >1) { grp3 <- apply(grp2, 2, function(x) length(unique(x)))
              if(any(grp3 < ncol(abund))) {
                if(length(grp3) >0) { useCol <- if(isTRUE(groupPref$lowNumberOfGroups)) which.min(grp3) else which(grp3 ==stats::median(grp3))[1]
                  setupSdSoft$lev <- grp2[,useCol]
                  names(setupSdSoft$lev) <- wrMisc::rmEnumeratorName(wrMisc::trimRedundText(txt=summaryD[,useCol], spaceElim=TRUE, silent=silent, debug=debug, callFrom=fxNa),
                    incl=c("anyCase","trim0","rmEnum"), sepEnum=c(" ","-","_"), nameEnum=c("Number","No","#","","Replicate","Sample"), silent=silent, debug=debug, callFrom=fxNa)
              } }
            } else {
              names(grp2) <- wrMisc::rmEnumeratorName(wrMisc::trimRedundText(txt=as.character(summaryD), spaceElim=TRUE, silent=silent, debug=debug, callFrom=fxNa),
                    incl=c("anyCase","trim0","rmEnum"), sepEnum=c(" ","-","_"), nameEnum=c("Number","No","#","","Replicate","Sample"), silent=silent, debug=debug, callFrom=fxNa)
              grp <- setupSdSoft$lev <- grp2
            }
          }

        }
         if(debug) { message(fxNa,"rSM4f"); rSM4f <- list() }

        summaryD <- as.data.frame(cbind(summaryD, grp=names(setupSdSoft$lev), grpInd=setupSdSoft$lev))                 # add presumed grouping to summaryD
      } else { if(!silent) message(fxNa,"Note : Abundance data are absent, adjust order of annotation to abundance data")}
    }
    if(debug) { message(fxNa,"rSM5"); rSM5 <- list(sdrf=sdrf,abund=abund,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,summaryD=summaryD,setupSdSoft=setupSdSoft) }


    ### 2  READ SDRF annotation & pick groups of replicates; has priority over grouping based on summary.txt
    if(length(sdrf) >0) {
      ## check if 'functional' sdrf (ie list) is provided -> use as is
      if(is.list(sdrf) && all(c("sdrfDat","col","lev") %in% names(sdrf), na.rm=TRUE)) {
        if(debug) message(fxNa,"Custom setupSd povdided as sdrf")
        sdrfDat <- sdrf$sdrfDat
      } else {
        ## may be : character vector (length <3) => assume path or sdrf accession, 2nd as sdrf-column to use
        ## may be : (matrix or) data.frame to use as table to exploit
        if(is.character(sdrf) && length(sdrf) <3 && length(dim(sdrf)) <2) {
          ## read sdrf from file or github
          sdrfDat <- readSdrf(sdrf, silent=silent, debug=debug, callFrom=fxNa)

          ## check for priority columns
          if(length(groupPref$sdrfColumn) ==1 && length(sdrf) <2) { ch1 <- groupPref$sdrfColumn %in% colnames(sdrfDat)
            if(any(ch1)) sdrf[2] <- which(ch1)[1]
          }
        } else {
          ## user provided custom sample annotation object
          if(length(dim(sdrf)) <2  && !silent) message(fxNa,"Note: 'sdrf' looks bizarre (trouble ahead ?), expecting either file, data.frame or complete list")
          sdrfDat <- sdrf
          sdrf <- "user provided custom object"}
      }
      if(debug) { message(fxNa,"rSM6  dim sdrfDat ",nrow(sdrfDat)," ",ncol(sdrfDat)); rSM6 <- list(sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,uplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,summaryD=summaryD,syncColumns=syncColumns) }


      ## 2.1 basic check (distinguish full $sampleSetup) form custom data.frame
      if(length(sdrfDat) >0) {
        syncColumns["sdrfDat"] <- FALSE    # initialize
        if(is.list(sdrfDat) && "sdrfDat" %in% names(sdrfDat)) { sdrfDat <- sdrfDat$sdrfDat
          if("groups" %in% names(sdrfDat)) groupPref$groups <- sdrfDat$groups
          if(debug) message(fxNa,"It seems a full $sampleSetup may have been given") }
        if(length(dim(sdrfDat)) <2) sdrfDat <- as.matrix(sdrfDat)
        if(length(abund) >0 && nrow(sdrfDat) != ncol(abund)) {
          if(!silent) message(fxNa,"NOTE : Ignoring 'sdrf'  : it does NOT have the expected number or rows (",nrow(sdrfDat)," given but ",ncol(abund)," expected !)")
          sdrf <- sdrfDat <- NULL }}

      ## 2.2  need to match lines (samples) of sdrf (setupDat) to summaryD and/or colnames of abund
      if(length(sdrfDat) >0) {
        if(length(summaryD) >0) {                   ## summaryD exist try matching by file-names
          chFiNames <- c("File.Name","File","FileName","Raw.file")     # search in summaryD
          chFiNa <- chFiNames %in% colnames(summaryD)
          if(debug) {message(fxNa,"rSM6a1") }
          if(any(chFiNa, na.rm=TRUE) && "comment.file.uri." %in% colnames(sdrfDat)) {
            ## align by filenames
            chFi <- match(sub("\\.zip$|\\.gz$","", basename(.corPathW(summaryD[,chFiNames[which(chFiNa)[1]]]))),
              sub("\\.zip$|\\.gz$","", basename(.corPathW(sdrfDat[,"comment.file.uri."]))))  # new order
            if(any(is.na(chFi)) && any(grepl("\\.raw",sdrfDat[,"comment.file.uri."]), na.rm=TRUE)) {
              sumDaFiNa <- sub("\\.raw","",sub("\\.zip$|\\.gz$","", basename(.corPathW(summaryD[,chFiNames[which(chFiNa)[1]]]))))
              sdrfFiNa <- sub("\\.raw","",sub("\\.zip$|\\.gz$","", basename(.corPathW(sdrfDat[,"comment.file.uri."]))))
              chFi <- match(sumDaFiNa, sdrfFiNa)  # new order
              if(any(is.na(chFi)) && "FP" %in% quantMeth) {
                sumDaFiNa <- wrMisc::rmEnumeratorName(wrMisc::trimRedundText(txt=sumDaFiNa, spaceElim=TRUE, silent=silent, debug=debug, callFrom=fxNa), newSep="_", incl=c("anyCase","trim0"), silent=silent, debug=debug, callFrom=fxNa)
                sdrfFiNa <- wrMisc::rmEnumeratorName(wrMisc::trimRedundText(txt=sdrfFiNa, spaceElim=TRUE, silent=silent, debug=debug, callFrom=fxNa), newSep="_", incl=c("anyCase","trim0"), silent=silent, debug=debug, callFrom=fxNa)
                chFi <- match(sumDaFiNa, sdrfFiNa)     # new order
                if(debug) { message(fxNa,"rSM6aa  dim sdrfDat ",nrow(sdrfDat)," ",ncol(sdrfDat))}
              }
              rmRaw <- TRUE
            } else rmRaw <- FALSE
            if(sum(is.na(chFi)) >0) warning(fxNa,"UNABLE to match all filenames from sdrf and ",basename(.corPathW(suplAnnotFile)),
               " ! \n  ++ BEWARE : Grouping of replicates may be incorrect !! \n") else {
              if(!silent && rmRaw) message(fxNa,"Note : Some filenames contain '.raw', others do NOT; solved inconsistency ..")
               # sdrfDat[chFi, c(10,22:24)]
              sdrfDat <- sdrfDat[chFi,]
              if(!silent) message(fxNa,"Successfully adjusted order of sdrf to content of ",basename(.corPathW(suplAnnotFile)))
            }
            syncColumns["sdrfDat"] <- TRUE
          } else if(!silent) message(fxNa, if(debug) "rSM6a  "," summaryD exists, but unable to find file-names")
          if(debug) {message(fxNa,"rSM6a"); rSM6a <- list() }

          #message("==save Image  rSM6a.Rdata ==");  save(sdrf,sdrfDat,abund,suplAnnotFile,quantMeth,abund,summaryD,setupSdSoft, file="C:\\E\\projects\\TCAmethods\\wrProteoRamus\\rSM6a.Rdata")
        } else {             # no summaryD, try clanames of abund
          if(length(abund) >0 && length(dim(abund)) >1) {
            ## requires utils::packageVersion("wrMisc") > "1.11.1"
            sdrfDat <- wrMisc::matchMatrixLinesToRef(mat=sdrfDat, ref=colnames(abund), addRef=TRUE, silent=silent, debug=debug, callFrom=fxNa)  # 2way-grep
            syncColumns["sdrfDat"] <- TRUE                         # really sure that synchronization successful ?
            if(debug) message(fxNa, "rSM6b "," Custom experimental setup was given")
          } else  message(fxNa, if(debug) "rSM6b ","Note : NO Additional information on filenames-order found, can't correct/adjust sdrf (ie sdrfDat) !!   rSM6x")
        }
      }
      if(debug) { message(fxNa,"rSM6c  dim sdrfDat ",nrow(sdrfDat)," ",ncol(sdrfDat)); rSM6c <- list() }

      ## 2.3  ready to make setupSd
      if(length(sdrfDat) >0) {
        if(TRUE) {
          if(length(sdrf) >1 ) {     # sdrf as data.frame (or list)           #  && length(dim(sdrf)) >1
            replStrOpt <- c("highest","lowest","min","max","median","combAll","combNonOrth")
            ## is it risky to search in 2nd value of sdrf ?
            chCol <- grep(paste0(sdrf[2],"$"), sub("\\.$","", colnames(sdrfDat)))
            if(length(chCol) >0) {
              if(debug) message(fxNa,"Using column '",colnames(sdrfDat)[chCol[1]],"' for setupSd")
              sdrf[2] <- chCol[1]
            } else if(!any(replStrOpt %in% sdrf[2])) {
              if(!silent) message(fxNa,"Invalid entry for 'sdrf[2]', resetting to default ('median')")
              sdrf[2] <- "median"
            }
            setupSd <- try(wrMisc::replicateStructure(sdrfDat, method=if(length(sdrf) >1) sdrf[2], silent=silent, callFrom=fxNa, debug=debug), silent=TRUE)
            if(inherits(setupSd, "try-error")) {message(fxNa,"Unable to figure out replicate Structure"); setupSd <- NULL; syncColumns["sdrfDat"] <- FALSE}

          } else {                  # sdrf given as file-name (eg git)
            combMeth <- if(length(groupPref$combMeth) ==1) groupPref$combMeth else "combNonOrth"
            setupSd <- list(combNonOrth=try(wrMisc::replicateStructure(sdrfDat, method=combMeth, silent=silent, callFrom=fxNa, debug=debug)),
              lowest=try(wrMisc::replicateStructure(sdrfDat, method="lowest", silent=silent, callFrom=fxNa, debug=debug)))
            ch1 <- sapply(setupSd, inherits, "try-error")
            if(all(ch1)) {message(fxNa,"UNABLE to understand replicate-structure from sdrf !!"); setupSd <- NULL; syncColumns["sdrfDat"] <- FALSE
            } else if(any(ch1)) {setupSd <- setupSd[which(!ch1)]; message(fxNa,"REMOVING one attempt of understanding replicate-structure")}
            if(debug) {message(fxNa,if(debug) "rSM6d  ","length setupSd ", length(setupSd)); rSM6d <- list() }
            ## choose among multiple options for grouping (number of groups)
            ch1 <- sapply(setupSd, function(x) length(x$lev[which(!duplicated(x$lev))]))
            lowNumberOfGroups <- if(length(groupPref) >0 && is.list(groupPref)) isTRUE(groupPref$lowNumberOfGroups) && length(groupPref$combMeth) !=1 else TRUE
            useSe <- if(any(ch1 ==1, na.rm=TRUE)) which(ch1 !=1) else if(isTRUE(lowNumberOfGroups)) which.min(ch1) else which.max(ch1)
            if(!silent) message(fxNa,"Choosing model '",names(useSe),"' for evaluating replicate-structure (ie ",ch1[useSe[1]]," groups of samples)" )
            setupSd <- setupSd[[useSe[1]]]             # select appropriate model
            if(debug) {message(fxNa,if(debug)"rSM6e  ","length setupSd ", length(setupSd)); rSM6e <- list() }
          }
          if("setupSd" %in% names(setupSd)) {setupSd <- wrMisc::partUnlist(setupSd, callFrom=fxNa,debug=debug);
            if(debug) message(fxNa,"rSM6e  - not expecting list of list for setupSd ! .. correcting"); rSM6e <- list(setupSd=setupSd,sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,uplAnnotFile=suplAnnotFile,summaryD=summaryD,setupSdSoft=setupSdSoft,quantMeth=quantMeth)}
          ## try to find usable names for setupSd$lev
          ## still need usable names for factor-levels
          newLe <- as.character(match(setupSd$lev, unique(setupSd$lev)))
          names(newLe) <- gsub("[[:space:]]|[[:punct:]]", "_", sdrfDat[,setupSd$col[1]])      # recuperate names out of sdrfDat
          if(debug) message(fxNa, if(debug)"rSM6f  ","Using as names for groups of replicates/levels :  ",wrMisc::pasteC(utils::head(unique(names(newLe)), 5)))
          setupSd[["lev"]] <- newLe
          if("setupSd" %in% names(setupSd)) { setupSd <- wrMisc::partUnlist(setupSd, callFrom=fxNa,debug=debug); if(debug) message(fxNa," rSM6g - not expecting list of list for setupSd ! .. correcting")}
          ## NOTE : names of levels may not be very meaningful/optimal
          ## option (future) : search in file-names for similar pattern
          if(debug) {message(fxNa,"rSM6h   names setupSd : ", wrMisc::pasteC(names(setupSd))); rSM6h <- list() }
        }
        if(debug) {message(fxNa,"rSM6i   names setupSd : ", wrMisc::pasteC(names(setupSd))); rSM6i <- list() }
        if(!is.list(setupSd)) {setupSd <- as.list(setupSd); if(debug) message(fxNa,"rSM6f  'setupSd' should be list, but was NOT !!")}
        if(!"sdrfDat" %in% names(setupSd)) setupSd$sdrfDat <- sdrfDat

        ## re-adjust numbers of levels
        iniNa <- names(setupSd$lev)
        newLev <- match(setupSd$lev, unique(setupSd$lev))
        chNa1 <- paste0("Sample_", sort(unique(sub("^Sample_","", iniNa))))
        chNa2 <- paste0("Group_", sort(unique(sub("^Group_","", iniNa))))
        if(all(chNa1 %in% iniNa) || all(chNa2 %in% iniNa)) { names(newLev) <- paste0("Group_", newLev)
        } else names(newLev) <- iniNa
        setupSd$lev <- newLev
        if(length(abund) >0) if(length(groupPref$groups) == ncol(abund)) {  setupSd$lev <- setupSd$groups <- groupPref$groups
          if(debug) message(fxNa,"Using custom groups (for extracted from full $sampleSetup) for  setupSd$lev ")
        }
        setupSd$annotBySoft <- as.data.frame(summaryD)
        setupSd$syncColumns <- syncColumns

      } else {      ## sdrf was given - but NOT conform : (no soft-generated sample annot available) try to match colnames of abund
        if(debug) message(fxNa, if(debug) "rSM6j ","NO valid sdrf found")
        ## ie single source of info
        if(length(summaryD) <1) {                                              ## ie no  summaryD
          if(debug) message(fxNa, if(debug) "rSM6k ","NO valid sdrf and NO valid information (summaryD) from quant-software found")

        } else {                         # ie summaryD is available
          setupSd <- setupSdSoft
          #message(fxNa,"Reading of sdrf was NOT successful and no summaryD => nothing can be done to mine experimental setup...")
        }
      }
      if(debug) { message(fxNa,"rSM7  head of setupSd$lev : ",wrMisc::pasteC(utils::head(setupSd$lev))); rSM7 <- list(setupSd=setupSd,sdrf=sdrf,sdrfDat=sdrfDat,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,summaryD=summaryD,nSamp0=nSamp0)}

      if(length(setupSd) >0) if(length(setupSd$lev) != nSamp0 && length(abund) >0) {               ## keep this ? - redundant !
        if(!silent) warning(fxNa, if(debug) "rSM7  ","Invalid information from sample meta-data or wrong experiment ! Number of samples from sdrf ",
          " (",length(setupSd$lev),") and from experimental data (",ncol(abund),") don't match !")
        setupSd <- NULL } else {
          if(length(abund) <1 && !silent) message(fxNa,"Note: Order of lines in sdrf not ajusted since no valid 'abund' given...")
          setupSd$level <- setupSd$lev
      }

    } else { setupSd <- setupSdSoft; setupSd$annotBySoft <- summaryD; setupSd$lev <- setupSd$lev }
    if(debug) { message(fxNa,"rSM8  head of setupSd$lev : ",wrMisc::pasteC(utils::head(setupSd$lev))); rSM_ <- list()}

  }
  ## finished readSampleMetaData
  setupSd }
  
