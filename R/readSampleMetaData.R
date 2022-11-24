#' Read Sample Meta-data from Quantification-Software And/Or Sdrf And Align To Experimental Data
#'
#' Sample annotation meta-data form \href{https://www.maxquant.org/}{MaxQuant}, ProteomeDiscoverer or similar, can be read using this function and relevant information extracted.
#' Furthermore, annotation in \href{https://github.com/bigbio/proteomics-metadata-standard}{sdrf-format} can be added (the order of sdrf will be adjated automatically, if possible).
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
#' For more information about the sdrf format please see \href{https://github.com/bigbio/proteomics-metadata-standard}{sdrf on github}.
#'
#'
#' @param sdrf (character, list or data.frame) optional extraction and adding of experimenal meta-data:
#'  if character, this may be the ID at ProteomeExchange or a similarly formatted local file. \code{sdrf} will get priority over \code{suplAnnotFile}, if provided.
#' @param suplAnnotFile (logical or character) optional reading of supplemental files produced by MaxQuant; if \code{gr} is provided, it gets priority for grouping of replicates
#'  if \code{TRUE} default to files 'summary.txt' (needed to match information of \code{sdrf}) and 'parameters.txt' which can be found in the same folder as the main quantitation results;
#'  if \code{character} the respective file-names (relative ro absolute path), 1st is expected to correspond to 'summary.txt' (tabulated text, the samples as given to MaxQuant) and 2nd to 'parameters.txt' (tabulated text, all parameters given to MaxQuant)
#' @param quantMeth (character) quantification method used
#' @param path (character) optional path of file(s) to be read
#' @param abund (matrix or data.frame) experimental quantitation data; only column-names will be used for aligning order of annotated samples
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a list with \code{$lev} and \code{$level} (grouping of samples given as integer), and \code{$meth} (method by which grouping as determined).
#'  If valid \code{sdrf} was given, the resultant list contains in addition \code{$sdrfDat} (data.frame of annotation). If software annotation has been found it will be shown in \code{$annotBySoft}.
#'  If all entries are invalid or entries do not pass the tests, this functions returns an empty \code{list}.
#' @seealso this function is used internally by \code{\link{readMaxQuantFile}} or \code{\link{readProtDiscovFile}}; \code{\link{readSdrf}} for reading sdrf-files, \code{\link[wrMisc]{replicateStructure}} for mining annotation columns
#' @examples
#' sdrf001819Setup <- readSampleMetaData("PXD001819")
#' str(sdrf001819Setup)
#'
#' @export
readSampleMetaData <- function(sdrf=NULL, suplAnnotFile=NULL, quantMeth="MQ", path=".", abund=NULL, silent=FALSE, debug=FALSE, callFrom=NULL)  {
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
  datOK <- all(c(length(sdrf) >0 | length(suplAnnotFile) >0, length(quantMeth)==1))
  if(length(abund) >0 & any(length(dim(abund)) !=2, dim(abund) < 1, na.rm=TRUE)) { datOK <- FALSE
    warning("Invalid argument 'abund'; must be matrix or data.frame with min 1 line and 1 col")}
  if(debug) {message(fxNa,"Ready search & extract sample meta-data   rSM0"); rSM0 <- list(sdrf=sdrf,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,path=path,abund=abund)}

  .corPathW <- function(x) gsub("\\\\", "/", x)
  .adjPat <- function(x) { out <- match(x, unique(x)); names(out) <- names(x); out}  # needed ??

  .compToRef <- function(mat, ref, addRef=TRUE, silent=FALSE, debug=FALSE, callFrom=NULL) {
    ## find best column for matching mat (matrix) to ref (char vector) via two way grep
    ## return list $grep (matched matrix), $col best colum
    fxNa <- wrMisc::.composeCallName(callFrom, newNa=".compToRef")
    out <- NULL              # initialize
    leM <- apply(nchar(as.matrix(mat)), 2, stats::median, na.rm=TRUE)
    leR <- stats::median(nchar(ref), na.rm=TRUE)     # median length
    ## grep ref in each col of mat
    gM <- apply(as.matrix(mat), 2, function(x) sapply(ref, grep, x))
    lM <- sapply(gM, sapply, length)          # matrix of counts
    ## 'reverse' grep in each col of mat in ref
    gR <- apply(as.matrix(mat), 2, sapply, grep, ref)
    lR <- if(is.list(gR)) {if(is.list(gR[[1]])) sapply(gR, sapply, length) else sapply(gR, length) } else sapply(gR, length)
    #if(is.list(lR)) lT <- sapply(gR, sapply, length) else if(length(dim(lR)) <2) lR <- matrix(lR, nrow=1,dimnames=list(NULL, names(lR)))
    if(is.list(lR)) warning(fxNa, " Trouble ahead, can't make matrix/vector of counts from grep in each col of mat in ref")
    if(debug) {message(fxNa,"..cTR1\n"); cTR1 <- list(mat=mat,ref=ref,leM=leM,leR=leR,gM=gM,lM=lM,gR=gR,lR=lR)}
    ##                         leM=leM,leR=leR,gM=gM,lM=lM,lR=lR,gR=gR,
    if(any(lM >0, lR >0)) {           # some (perfect or partial) solutions
      chM1 <- apply(lM, 2, function(x) all(x==1))
      chR1 <- apply(lR, 2, function(x) all(x==1))
      if(sum(chM1, chR1) >0) {         # ideal solution exists
         if(any(chM1)) { out <- list(by="mat", colNo=which(chM1), le=if(length(dim(lM)) >1) lM[,which(chM1)] else lM[which(chM1)], grep=gM[[which(chM1)]])} else {
           out <- list(by="ref", colNo=which(chR1), le=if(length(dim(lR)) >1) lR[,which(chR1)] else lR[which(chR1)], grep=gR[[which(chR1)]]) }
      } else {                        # less ideal (multiple hits and/or empty)
        chM1 <- colSums(lM >0)
        chR1 <- colSums(lR >0)
        if(any(chM1==nrow(mat), chR1==nrow(mat))) {   # multiple hits, no empty, need further refinement
           if(any(chM1==nrow(mat))) { out <- list(by="mat", colNo=which(chM1==nrow(mat))[1],
             le=if(length(dim(lM)) >1) lM[,which(chM1==nrow(mat))[1]] else lM[which(chM1==nrow(mat))[1]], grep=gM[[which(chM1==nrow(mat))[1]]])
           } else if(any(chR1==nrow(mat))) { out <- list(by="mat", colNo=which(chR1==nrow(mat))[1],
             le=if(length(dim(lR)) >1) lR[,which(chR1==nrow(mat))[1]] else lR[which(chR1==nrow(mat))[1]], grep=gR[[which(chR1==nrow(mat))[1]]])}
          ## resolve multiple hits in $grep by picking 1st not yet used
          newGr <- sapply(out$grep, function(x) x[1])
          multX <- unique(names(newGr[which(newGr >1)]))                 # 5jul22
          for(i in 1:length(multX)) newGr[which(names(newGr)==multX[i])] <- out$grep[[multX[i]]]
          out$grep <- newGr
        } else warning(fxNa,"Impossible to find complete solution, returning NULL")
      }
      if(out$by=="ref") { message(fxNa,"Best matching by direct matching")
         out <- mat[out$grep,]
        if(!isFALSE(addRef)) out <- if("ref" %in% colnames(out)) cbind(out, ref) else cbind(out, ref=ref)
      } else {  message(fxNa,"Best matching by reverse matching")
        ref2 <- ref[out$grep]
          #still need prev resuts ?# outIni=out
        out <- mat
        if(!isFALSE(addRef)) out <- if("ref" %in% colnames(out)) cbind(out, ref2) else cbind(out, ref=ref2)
        out <- mat[match(ref,ref2),]
      }
    } else warning(fxNa,"Impossible to find any matches, returning NULL")
    out }
    ## example
    # mat1 <- matrix(paste0("__",letters[rep(c(1,1,2,2,3),3) +rep(0:2,each=5)], rep(1:5)), ncol=3)
    # .compToRef(mat1, paste0(letters[c(3,4,5,3,4)],c(1,3,5,2,4)))
    # mat2 <- matrix(paste0("__",letters[rep(c(1,1,2,2,3),3) +rep(0:2,each=5)], c(rep(1:5,2),1,1,3:5 )), ncol=3)
    # .compToRef(mat2, paste0(letters[c(3,4,5,3,4)],c(1,3,5,1,4)))
    # mat3 <- matrix(paste0(letters[rep(c(1,1,2,2,3),3) +rep(0:2,each=5)], c(rep(1:5,2),1,1,3,3,5 )), ncol=3)
    # .compToRef(mat3, paste0("__",letters[c(3,4,5,3,4)],c(1,3,5,1,3)))

  ## end suppl fx


  path <- if(length(path) <1) "." else path[1]
  nSamp0 <- if(length(dim(abund)) >1) ncol(abund) else 0

  if(datOK) {
    ### GROUPING OF REPLICATES AND SAMPLE META-DATA

    ## SOFTWARE specific META-DATA : read additional annotation & documentation files produced by MaxQuant (summary.txt & parameters.txt) or PD (".InputFiles.txt")
    if(length(suplAnnotFile) >0)  {     # read quant software-generated sample annotation
      chFiNa <- NULL                    # initialize
      if(debug) {message(fxNa,"rSM1"); rSM1 <- list(sdrf=sdrf,abund=abund,path=path,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth) }

      ## MaxQuant
      if("MQ" %in% quantMeth) {
        isDir <- if(is.character(suplAnnotFile)) utils::file_test("-d",suplAnnotFile[1]) else FALSE
        if(isDir) { path <- suplAnnotFile[1]; suplAnnotFile <- TRUE}
        if(isTRUE(suplAnnotFile)) {      # automatic search for standard file-names ('summary.txt','parameters.txt') in same dir as main MaxQuant data
          chFiNa <- c("summary.txt","summary.txt.gz","parameters.txt","parameters.txt.gz")
          chFi <- file.exists(file.path(path, chFiNa))
          if(debug) {message(fxNa,"rSM0a\n"); rSM0a <- list(path=path,sdrf=sdrf,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,chFi=chFi,chFiNa=chFiNa )}
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
          if(debug & chFi[2]) message(fxNa,"Successfully read ",quantMeth," parameters from '",suplAnnotFile[2],"'") }
        if(debug) { message(fxNa,"rSM1mq2")}
      }

      ## ProteomeDiscoverer
      ## uses suplAnnotFile as path for '.InputFiles\\.txt'
      if("PD" %in% quantMeth) {
        if(debug) {message(fxNa,"rSM1pd"); rSM1pd <- list(sdrf=sdrf,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth)}
        if(length(suplAnnotFile) >1) { if(!silent) message(fxNa,"Only 1st value of argument 'suplAnnotFile' can be used with quantMeth=PD")
          suplAnnotFile <- suplAnnotFile[1] }
        if(isTRUE(suplAnnotFile)) {      # automatic search for standard file-name ('InputFiles.txt') in same dir as main MaxQuant data
          suplAnnotFile <- list.files(path=path, pattern=".InputFiles\\.txt$|.InputFiles\\.txt\\.gz$")
          if(length(suplAnnotFile) >1) { if(!silent) message(fxNa,"Found ",length(suplAnnotFile)," files matching general patter, using ONLY 1st, ie ",suplAnnotFile[1])
            suplAnnotFile <- suplAnnotFile[1] }
          chFi <- length(suplAnnotFile) >0
          if(!chFi & !silent) message(fxNa,"Note: Unable to (automatically) find sample-annotation file. Maybe it was not exported from ProteomeDiscoverer ?")
        } else chFi <- try(file.exists(file.path(path, suplAnnotFile)), silent=TRUE)
        if(inherits(chFi, "try-error") & silent) {chFi <- FALSE; message(fxNa,"Meta-data: Failed to see file '",suplAnnotFile[1]," ! (check if file exists or rights to read directory ?)")}
        if(debug) {message(fxNa,"rSM1pd2") }

        ## main reading of PD sample meta-data
        if(chFi) summaryD <- try(utils::read.delim(file.path(path, suplAnnotFile[1]), stringsAsFactors=FALSE), silent=TRUE)
        if(inherits(summaryD, "try-error")) {summaryD <- NULL; if(!silent) message(fxNa,"Meta-data: Failed to read '",suplAnnotFile[1],"' !")
        } else {
          if(debug)  message(fxNa,"ProteomeDiscoverer Meta-data sucessfully read '",suplAnnotFile[1])}
        if(debug) {message(fxNa,"rSM1pd3")}
      }

      ## Proline
      ## so far only for reading out of xslx
      if("PL" %in% quantMeth) {
        if(debug) {message(fxNa,"rSM0pl"); rSM0pl <- list(sdrf=sdrf,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth)}
        summaryD <- NULL
        ## need suplAnnotFile
        if(length(grep("\\.xlsx$", suplAnnotFile)) >0) {
          ## Extract out of Excel
          reqPa <- c("readxl")
          chPa <- sapply(reqPa, requireNamespace, quietly=TRUE)
          if(any(!chPa)) message("package( '",paste(reqPa[which(!chPa)], collapse="','"),"' not found ! Please install first from CRAN") else {

            sheets <- if(debug) try(readxl::excel_sheets(suplAnnotFile), silent=TRUE) else suppressMessages(try(readxl::excel_sheets(suplAnnotFile), silent=TRUE))
            if(debug) {message(fxNa," ww3")}
            if(inherits(sheets, "try-error")) { message(fxNa,"Unable to read file '",suplAnnotFile,"' ! Returning NULL; check format & rights to read")
            } else {
              annShe <- c("Import and filters", "Search settings and infos")     # sheets from xslx to try reading for sample/meta-information
              annSh <- wrMisc::naOmit(match(annShe, sheets))
              if(length(annSh) <1) annSh <- grep("Import", sheets)
              if(length(annSh) >1) {
                if(!silent) message(fxNa,"Multipe sheets containing 'Import' found, using 1st :",sheets[annSh[1]])
                annSh <- annSh[1]
              } else if(length(annSh) <1 & !silent) {
                message(fxNa,"NOTE: NONE of ANNOTATION SHEETS (",wrMisc::pasteC(annShe),") in '",suplAnnotFile,"' FOUND !  Can't check Matching order of samples to sdrf-anotation !")
              }
              summaryD <- as.matrix(as.data.frame(if(debug) readxl::read_xlsx(suplAnnotFile, sheet=annSh, col_names=FALSE) else suppressMessages(readxl::read_xlsx(suplAnnotFile, sheet=annSh, col_names=FALSE))))
              rownames(summaryD) <- summaryD[,1]
              summaryD <- t(summaryD[,-(1:2)])
              rownames(summaryD) <- 1:nrow(summaryD)
            }
          }
        }
        if(debug) {message(fxNa,"rSM1pl")}
      }                 # finish PL

      ## other software ? ..
      chSoft <- c("MQ","PD","PL")
      if(!any(quantMeth %in% chSoft, !silent, na.rm=TRUE)) message(fxNa,"Note: No specific procedure has been implemented so far for gathering meta-data by the analysis-software/method '",quantMeth,"'")
    }
    if(debug) { message(fxNa,"rSM2"); rSM2 <- list(sdrf=sdrf,abund=abund,uplAnnotFile=suplAnnotFile,quantMeth=quantMeth,summaryD=summaryD,suplAnnotFile=suplAnnotFile) }

    ## adjust summaryD to quant data
    if(length(summaryD) >0) {      ##  more checks
      if(length(abund) <0) warning(fxNa,"Can't verify/correct names of anotation since content of 'abund' has no colnames or is NULL") else {
        if(!identical(ncol(abund), nrow(summaryD))) { summaryD <- NULL
		  if(!silent) message(fxNa,"NOTE : Number of columns of 'abund' does NOT FIT to number of samples in annotation-data !")	}
	  }
      if(length(dim(summaryD)) !=2) summaryD <- matrix(summaryD, ncol=1, dimnames=list(names(summaryD),NULL))
	}

    if(length(summaryD) >0) {      ## define sdrf as setupSd
      ## need to match colnames(abund) to (MQ:) $Raw.file or $Experiment  .. need to find best partial match
      if("MQ" %in% quantMeth) {         ## NOT IN SAME ORDER !!
        summaryD <- summaryD[,wrMisc::naOmit(match(c("Raw.file","Experiment"), colnames(summaryD)))]           # cor 21oct22
        chSd <- length(abund) >0 & nrow(summaryD) == ncol(abund)
        if(length(chSd) <1) chSd <- FALSE
        ## normally  colnames(abund) and summaryD should alread be in correct order
        if(!chSd) {
          if(!silent & length(abund) >0) if(nrow(summaryD) == ncol(abund)) message(fxNa,"PROBLEM : meta-data and abundance data do not match !  ",
             "Number of samples from ",suplAnnotFile[1]," (",nrow(summaryD),") and from main data (",ncol(abund),") do NOT match !! .. ignoring") }
        if(debug) { message(fxNa," .. rSM4mq"); rSM4mq <- list(sdrf=sdrf,abund=abund,uplAnnotFile=suplAnnotFile,quantMeth=quantMeth,summaryD=summaryD)}
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
      }

      if("PL" %in% quantMeth) {   ## order OK ?
        chSd <- length(abund) >0 & nrow(summaryD) == ncol(abund)
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
          if(!silent & nrow(summaryD) == ncol(abund)) message(fxNa,"PROBLEM : Invalid meta-data !  ", "Number of samples from ",
          suplAnnotFile[1]," (",nrow(summaryD),") and from main data (",ncol(abund),") do NOT match !! .. ignoring") }
        }
      }
      ## other software ? ...

      if(debug) { message(fxNa,"rSM4d"); rSM4d <- list(sdrf=sdrf,abund=abund,uplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,summaryD=summaryD) }

      ## finish adjusting
      #setupSdSoft <- replicateStructure(summaryD, silent=silent, debug=debug, callFrom=fxNa)
      setupSdSoft <- wrMisc::replicateStructure(summaryD, silent=silent, debug=debug, callFrom=fxNa)

      if(debug) { message(fxNa,"rSM4e"); rSM4e <- list(sdrf=sdrf,setupSdSoft=setupSdSoft,abund=abund,uplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,summaryD=summaryD) }
      ## so far no direct information about groups (all filenames are different), need to try to find out (remove enumerators)
      if(length(abund) >0) {
        #grpA <- .trimToRedund(colNa=colnames(abund), silent=silent, debug=debug, callFrom=fxNa)
        grpA <- wrMisc::trimRedundText(txt=colnames(abund), spaceElim=TRUE, silent=silent, debug=debug, callFrom=fxNa)                 # 26oct22
        grpS <- setupSdSoft$lev
        grpLe <- c(abu=length(unique(grpA)), sum=length(unique(setupSdSoft$lev)))
        if(any(grpLe > 1 & grpLe < ncol(abund), na.rm=TRUE)) {
          grp <- list(grpA,grpS)[[which(grpLe > 1 & grpLe < ncol(abund))[1]]]
          setupSdSoft$lev <- grp
        } else grp <- setupSdSoft$lev
        ## re-adjust numbers of levels
        grpNa <- if(length(names(grp)) >0) names(grp) else grp
        grp <- match(grp, unique(grp))
        names(grp) <- grpNa
        setupSdSoft$lev <- grp
        summaryD <- as.data.frame(cbind(summaryD, grp=grp))                 # add presumed grouping to summaryD
      } else { if(!silent) message(fxNa,"Note : abundance data are absent, adjust order of annotation to abundance data")}
    }
    if(debug) { message(fxNa,"rSM5"); rSM5 <- list(sdrf=sdrf,abund=abund,uplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,summaryD=summaryD,setupSdSoft=setupSdSoft) }

    ### READ SDRF annotation & pick groups of replicates; has priority over grouping based on summary.txt
    if(length(sdrf) >0) {
      ## check if 'functional' sdrf (ie list) is provided -> use as is
      if(is.list(sdrf) & all(c("sdrfDat","col","lev") %in% names(sdrf), na.rm=TRUE)) sdrfDat <- sdrf$sdrfDat else {
        ## may be : character vector (length <3) => assume path or sdrf accession, 2nd as sdrf-column to use
        ## may be : (matrix or) data.frame to use as table to exploit
        if(all(is.character(sdrf) & length(sdrf) <3, length(dim(sdrf)) <2, na.rm=TRUE)) {  # length(dim(sdrf)) ==2 & dim(sdrf) > 1
          ## read sdrf from file or github
          sdrfDat <- readSdrf(sdrf, silent=silent, debug=debug, callFrom=fxNa)
        } else {sdrfDat <- sdrf}
      }
      if(debug) { message(fxNa,"rSM6"); rSM6 <- list(sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,uplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,summaryD=summaryD) }

      ## need to match lines (samples) of sdrf (setupDat) to colnames of abund or summaryD
      if(length(summaryD) >0 & length(sdrfDat) >0) {                   ## summaryD exist try matching by file-names
        chFiNames <- c("File.Name","File","FileName")     # search in summaryD
        chFiNa <- chFiNames %in% colnames(summaryD)
        if(any(chFiNa, na.rm=TRUE) & "comment.file.uri." %in% colnames(sdrfDat)) {
          ## align by filenames
          chFi <- match(sub("\\.zip$|\\.gz$","", basename(.corPathW(summaryD[,chFiNames[which(chFiNa)[1]]]))), 
            sub("\\.zip$|\\.gz$","", basename(.corPathW(sdrfDat[,"comment.file.uri."]))))  # new order
          if(sum(is.na(chFiNa)) >0) warning(fxNa,"Unable to match all filenames from sdrf and ",basename(.corPathW(suplAnnotFile)),
             " ! \n  Beware : Grouping of replicates may be incorrect !!") else {
             # sdrfDat[chFi, c(10,22:24)]
            sdrfDat <- sdrfDat[chFi,]
            if(!silent) message(fxNa,"Sucessfully adjusted order of sdrf to content of ",basename(.corPathW(suplAnnotFile)))
          }
        } else message(fxNa," summaryD exists, but unable to find file-names")
        if(debug) {message(fxNa,"rSM6b   length setupSd ", length(setupSd))}

        if(length(setupSd) <1) setupSd <- wrMisc::replicateStructure(sdrfDat, method=if(length(sdrf) >1) sdrf[2] else "combNonOrth", silent=silent, callFrom=fxNa, debug=debug)
        setupSd$sdrfDat <- sdrfDat
        if(debug) {message(fxNa,"rSM6c   length setupSd ", length(setupSd)); rSM6c <- list(setupSd=setupSd,sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,uplAnnotFile=suplAnnotFile,summaryD=summaryD,setupSdSoft=setupSdSoft,quantMeth=quantMeth) }
        ## try to find usable names for setupSd$lev
        ## check with cols of sdrfDat have same pattern as setupSd$lev
        chSdrfCol <- c("comment.data.file.","factor.value.spiked.compound.")
        rePa <- apply(sdrfDat, 2, function(x) all(match(x, unique(x)) == setupSd$lev, na.rm=TRUE))
        if(sum(rePa) >1) {            # multiple columns for choice
          rePa2 <- apply(sdrfDat[,which(rePa)], 2, wrMisc::trimRedundText, spaceElim=TRUE, silent=silent, callFrom=fxNa, debug=debug)           # 26 oct22
          #rePa2 <- apply(sdrfDat[,which(rePa)], 2, .trimToRedund, highSep=TRUE, silent=silent, callFrom=fxNa, debug=debug)

          if(any(chSdrfCol %in% colnames(rePa2), na.rm=TRUE)) rePa <- rePa2[,which(colnames(rePa2) %in% chSdrfCol)[1]] else {
            rePa3 <- apply(sdrfDat[,which(rePa)], 2, function(x) stats::median(nchar(x), na.rm=TRUE))
            rePa <- rePa2[,order(rePa3)[ceiling(length(rePa3)/2)]] }  # pick column with ceiling of median length of names
        }
        if(debug) {message(fxNa,"rSM6e   length setupSd ", length(setupSd)); rSM6e <- list(setupSd=setupSd,sdrf=sdrf,sdrfDat=sdrfDat,abund=abund,uplAnnotFile=suplAnnotFile,summaryD=summaryD,setupSdSoft=setupSdSoft,quantMeth=quantMeth) }
        ## re-adjust numbers of levels
        setupSd$lev <- match(setupSd$lev, unique(setupSd$lev))
        names(setupSd$lev) <- rePa
        setupSd$annotBySoft <- as.data.frame(summaryD)

      } else {      ## (no soft-generated sample annot available) try to match colnames of abund
        if(debug) message(fxNa,"NO valid sdrf and/or summaryD found")
        ## ie no  summaryD
        if(length(abund) >0 & length(dim(abund)) >1 & length(sdrfDat) >0) sdrfDat <- .compToRef(mat=sdrfDat, ref=colnames(abund), addRef=TRUE, silent=silent, debug=debug, callFrom=fxNa) else {
          if(!silent) message(fxNa,"note : No (valid) experimental data given for adjusting order of samples from Sdrf annotation")
          if(length(setupSdSoft) >0) setupSd <- setupSdSoft else {
            if(length(sdrfDat) >0) {
              setupSd <- wrMisc::replicateStructure(sdrfDat, method=if(length(sdrf) >1) sdrf[2] else "combNonOrth", silent=silent, callFrom=fxNa, debug=debug)
              setupSd$sdrfDat <- sdrfDat }
          }
        }
      }
      if(debug) { message(fxNa,"rSM7"); rSM7 <- list(setupSd=setupSd,sdrf=sdrf,sdrfDat=sdrfDat,suplAnnotFile=suplAnnotFile,quantMeth=quantMeth,abund=abund,summaryD=summaryD,nSamp0=nSamp0)}

      if(length(setupSd) >0) if(length(setupSd$lev) != nSamp0 & length(abund) >0) {               ## keep this ? - redundant !
        if(!silent) warning(fxNa,"Invalid information from sample meta-data or wrong experiment ! Number of samples from sdrf ",
          " (",length(setupSd$lev),") and from experimental data (",ncol(abund),") don't match !")
        setupSd <- NULL } else {
          if(length(abund) <1 & !silent) message(fxNa,"Note: Order of lines in sdrf not ajusted since no valid 'abund' given...")
          setupSd$level <- setupSd$lev }

    } else { setupSd <- setupSdSoft; setupSd$annotBySoft <- summaryD; setupSd$level <- setupSd$lev }
    if(debug) { message(fxNa,"rSM8")}
  }
  ## finished readSampleMetaData
  setupSd }
  
