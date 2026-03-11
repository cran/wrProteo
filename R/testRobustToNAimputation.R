#' Pairwise Testing Robust To NA-Imputation
#'
#' This function replaces \code{NA} values based on group neighbours (based on grouping of columns in argument \code{gr}), following overall assumption of close to Gaussian distribution.
#' Furthermore, it is assumed that \code{NA}-values originate from experimental settings where measurements at or below detection limit are recoreded as \code{NA}.
#' In such cases (eg in proteomics) it is current practice to replace \code{NA}-values by very low (random) values in order to be able to perform t-tests.
#' However, random normal values used for replacing may in rare cases deviate from the average (the 'assumed' value) and in particular, if multiple \code{NA} replacements are above the average, 
#' may look like induced biological data and be misinterpreted as so.   
#'    
#' @details
#' The statistical testing uses \code{eBayes} from Bioconductor package \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{limma} for robust testing in the context of small numbers of replicates. 
#' By repeating multiple times the process of replacing \code{NA}-values and subsequent testing the results can be sumarized afterwards by median over all repeated runs to remmove the stochastic effect of individual NA-imputation.
#' Thus, one may gain stability towards random-character of \code{NA} imputations by repeating imputation & test 'nLoop' times and summarize p-values by median (results stabilized at 50-100 rounds).
#' It is necessary to define all groups of replicates in \code{gr} to obtain all possible pair-wise testing (multiple columns in $BH, $lfdr etc). 
#' The modified testing-procedure of Bioconductor package \href{https://bioconductor.org/packages/release/bioc/html/ROTS.html}{ROTS} may optionaly be included, if desired.
#' This function returns a \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{limma}-like S3 list-object further enriched by additional fields/elements.
#' 
#' The argument \code{multCorMeth} allows to choose which multiple correction algorimths will be used and included to the final results.
#' Possible options are 'lfdr','BH','BY','tValTab', ROTSn='100' (name to element necessary) or 'noLimma' (to add initial p.values and BH to limma-results). By default 'lfdr' (local false discovery rate from package 'fdrtools') and 'BH' (Benjamini-Hochberg FDR) are chosen.
#' The option 'BY' referrs to Benjamini-Yakuteli FDR, 'tValTab' allows exporting all individual t-values from the repeated NA-substitution and subsequent testing.
#' 
#' This function is compatible with automatic extraction of experimental setup based on sdrf or other quantitation-specific sample annotation.
#' In this case, the results of automated importing and mining of sample annotation should be stored as \code{$sampleSetup$groups} or \code{$sampleSetup$lev}  
#' 
#' It is possible to limit the pairwise combinations to a custom designed set using the argument \code{useComparison}.
#' This may be a matrix where very line designs a new pairwise comparison and the first column is refers to 'sample' while the second column assigns the 'reference'.
#' Otheriwse, one may provide a character vector as \code{useComparison} where each entry referes to a new pairwise comparison. 
#' In this case it is recommended to use '--' as separator for combining the group-names to be used as 'sample' and 'reference'.
#' In case the names of groups (argument \code{grp}) do not contain any '-' the single character separator '-' may also be used in compatibility with use in bioconductor packae limma.
#' However, if \code{grp}) does contain any '-', results will be presented using '--' as separator to prevent any confusion.
#' 
#' @return This function returns a limma-type MA-object ('MArrayLM', can be handeled like a list).
#' For details 'on choice of NA-impuation procedures with arguments 'imputMethod' and 'avSd' please see  \code{\link{matrixNAneighbourImpute}}.
#' 
#' @param dat (matrix or data.frame) main data (may contain \code{NA}); if \code{dat} is list containing $quant and $annot as matrix, the element $quant will be used
#' @param gr (character or factor) replicate association; if \code{dat} contains a list-element \code{$sampleSetup$groups} or \code{$sampleSetup$lev} this may be used in case \code{gr=NULL}
#' @param useComparison (list, character or matrix) optional argument allowing to specify which pairwise comparions sould be performed, 
#'   default \code{useComparison=NULL} will run all pairwise comparisons;
#'   may be character combining two group-names (from argument \code{grp}) separated by a '-' (eg 'A-B') or matrix where the rownames design the elements to be compared as pairwise;
#'   It is also possible to give a list with $sep (separator to be used when combining) and $useComparison (regular comparisons)
#'   Note : the names of the groups may not contain any '-' to avoid confucing them with pairwise separators !
#' @param annot (matrix or data.frame) annotation (lines must match lines of data !), if \code{annot} is \code{NULL} and argument \code{dat} is a list containing both $quant and $annot, the element $annot will be used 
#' @param retnNA (logical) retain and report number of \code{NA}
#' @param avSd (numerical,length=2) population characteristics (mean and sd) for >1 \code{NA}-neighbours (per line)
#' @param avSdH deprecated, please use \code{avSd} inestad; (numerical,length=2) population characteristics 'high' (mean and sd) for >1 \code{NA}-neighbours (per line)
#' @param plotHist (logical) additional histogram of original, imputed and resultant distribution (made using \code{\link{matrixNAneighbourImpute}} )
#' @param xLab (character) custom x-axis label 
#' @param tit (character) custom title
#' @param imputMethod (character) choose the imputation method (may be 'mode2'(default), 'mode1', 'datQuant', 'modeAdopt', 'informed' or 'none', for details see \code{\link{matrixNAneighbourImpute}} )
#' @param seedNo (integer) seed-value for normal random values
#' @param multCorMeth (character) define which method(s) for correction of multipl testing should be run (for choice : 'BH','lfdr','BY','tValTab', choosing several is possible)
#' @param nLoop (integer) number of runs of independent \code{NA}-imputation
#' @param lfdrInclude (logical) deprecated, please used \code{multCorMeth} instead (include lfdr estimations, may cause warning message(s) concerning convergence if few too lines/proteins in dataset tested).
#' @param ROTSn (integer) deprecated, please used \code{multCorMeth} instead (number of repeats by \code{ROTS}, if \code{NULL} \code{ROTS} will not be called)
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) This function allows easier tracking of messages produced
#' @return This function returns a limma-type S3 object of class 'MArrayLM' (which can be accessed lika a list); multiple results of testing or multiple testing correction types may get included ('p.value','FDR','BY','lfdr' or 'ROTS.BH')
#' @seealso NA-imputation via \code{\link{matrixNAneighbourImpute}}, modereated t-test without NA-imputation \code{\link[wrMisc]{moderTest2grp}}, calculating lfdr \code{\link[wrMisc]{pVal2lfdr}}, \code{eBayes} in Bioconductor package \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{limma}, \code{\link[stats]{t.test}},\code{ROTS} of Bioconductor package \href{https://bioconductor.org/packages/release/bioc/html/ROTS.html}{ROTS}   
#' @examples
#' set.seed(2015); rand1 <- round(runif(600) +rnorm(600,1,2),3)
#' dat1 <- matrix(rand1,ncol=6) + matrix(rep((1:100)/20,6),ncol=6)
#' dat1[13:16,1:3] <- dat1[13:16,1:3] +2      # augment lines 13:16 
#' dat1[19:20,1:3] <- dat1[19:20,1:3] +3      # augment lines 19:20
#' dat1[15:18,4:6] <- dat1[15:18,4:6] +1.4    # augment lines 15:18 
#' dat1[dat1 <1] <- NA                        # mimick some NAs for low abundance
#' ## normalize data
#' boxplot(dat1, main="Data Before Normalization", las=1)
#' dat1 <- wrMisc::normalizeThis(as.matrix(dat1), meth="median")
#' ## designate replicate relationships in samples ...  
#' grp1 <- gl(2, 3, labels=LETTERS[1:2])                   
#' ## moderated t-test with repeated imputations (may take >10 sec, >60 sec if ROTSn >0 !) 
#' PLtestR1 <- testRobustToNAimputation(dat=dat1, gr=grp1, retnNA=TRUE, nLoop=20)
#' names(PLtestR1)
#' head(PLtestR1$p.value)
#' head(PLtestR1$BH)
#' head(PLtestR1$means)
#' boxplot(PLtestR1$datImp, main="Data At Normalization & Imputation", las=1)
#' 
#' ## custom selection of comparisons (incl custom orientation)
#' useComp <- c("A-B", "B-A")      # You can choose orientation sample/reference
#' PLtestR2 <- testRobustToNAimputation(dat=dat1, gr=grp1, useComparison=useComp, 
#'   retnNA=TRUE, nLoop=20)
#' head(PLtestR2$BH)
#' head(PLtestR2$means)
#' 
#' @export
testRobustToNAimputation <- function(dat, gr=NULL, useComparison=NULL, annot=NULL, retnNA=TRUE, avSd=c(0.15,0.5), avSdH=NULL, plotHist=FALSE, xLab=NULL, tit=NULL, imputMethod="mode2", 
  seedNo=NULL,  multCorMeth=c("lfdr","FDR"), nLoop=100, lfdrInclude=NULL, ROTSn=NULL, silent=FALSE, debug=FALSE, callFrom=NULL) {
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="testRobustToNAimputation")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  if(!isTRUE(plotHist)) plotHist <- FALSE
  datOK <- length(dat) >0
  msg <- grIni <- sampleSetup <- sep1 <- useComparisonNa <-NULL
  ## start testing input
  if(datOK >0 && is.list(dat)) { if(all(c("quant","annot") %in% names(dat))) {
    if(length(dim(dat$annot)) ==2 && length(annot) <1) annot <- dat$annot else if(!silent) message(fxNa,"Invalid '$annot' (or ignoring due to conflict with argument 'annot')") # recover$annot if not given separately
    if(debug) {message(fxNa,"  tRN1"); tRN1 <- list(dat=dat, gr=gr,useComparison=useComparison,annot=annot,retnNA=retnNA,avSd=avSd )}   # dat=dat,
     
    if(debug) message(fxNa," .. tRN1b")
    if("sampleSetup" %in% names(dat) && length(gr) <1) {
      gr <- if("groups" %in% names(dat$sampleSetup)) dat$sampleSetup$groups else dat$sampleSetup$level
      ## rather use names instead of index
      if(all(is.integer(gr)) && all(nchar(names(gr)) > 0)) { grIni <- gr
        gr <- names(gr)
        names(gr) <- as.character(grIni)
      }
      if("quant" %in% names(dat) && length(gr) != ncol(dat$quant) && all(c("sdrfDat","col") %in% names(dat$sampleSetup))) gr <- dat$sampleSetup$sdrfDat[,dat$sampleSetup$col] 
      if(length(names(gr)) <1) names(gr) <-  match(gr, unique(gr))
    } else grIni <- gr
    if(debug) {message(fxNa," .. tRN1c"); tRN1c <- list(dat=dat, gr=gr,grIni=grIni,useComparison=useComparison,annot=annot,retnNA=retnNA,avSd=avSd )}
    
    if("sampleSetup" %in% names(dat)) {       # retreive useComparison from dat$sampleSetup (if available)
      sampleSetup <- dat$sampleSetup
      if(length(useComparison) ==0 && "compareTo" %in% names(sampleSetup) && all(dim(sampleSetup$compareTo) > c(0,1)) ) {
        useComparison <- sampleSetup$compareTo}  
      if(length(useComparison) >0 && is.list(useComparison) && "useComparison" %in% names(useComparison))  {
        if("sep" %in% names(useComparison)) sep1 <- useComparison$sep
        if("GrpNames" %in% names(useComparison)) useComparisonNa <- useComparison$GrpNames
        useComparison <- useComparison$useComparison
      }
      if(debug) message(fxNa,"sampleSetup$compareTo provided in 'dat' - using all ",nrow(useComparison)," comparisons")
      ## recuperate groups as substitute to gr
      if(length(gr) <1  && "groups" %in% names(dat$sampleSetup)) {
        gr <- dat$sampleSetup$groups }
    } 
    
    dat <- dat$quant 
  } else { datOK <- FALSE; msg <- "Invalid 'dat' : does NOT contain both '$quant' and '$annot !"} } 
  
  if(datOK) { if(length(unique(gr))==length(gr)) { datOK <- FALSE
    msg <- "Argument 'gr' is empty or does NOT design any replicates !  (nothing to do)"} } 
  if(datOK) if(any(length(dim(dat)) !=2, dim(dat) < 1:2, na.rm=TRUE)) { datOK <- FALSE
    msg <- "'dat' must be matrix or data.frame with >1 columns"}  
  if(debug) {message(fxNa," .. tRN1d"); tRN1d <- list(dat=dat,datOK=datOK, gr=gr,grIni=grIni,useComparison=useComparison,annot=annot,retnNA=retnNA,avSd=avSd )}

  if(datOK) {
    if(is.data.frame(dat)) dat <- as.matrix(dat)
    if(length(gr) != ncol(dat)) { datOK <- FALSE
      msg <- "Number of columns in 'dat' and number of (group-)elements in 'gr' do NOT match !"} }
  if(datOK) {
    if(!is.factor(gr)) gr <- as.factor(gr)
    if(is.null(xLab)) xLab <- "values"            
    if(length(annot) <1) annot <- matrix(NA, nrow=nrow(dat), ncol=1, dimnames=list(rownames(dat),"rowNa"))
    if(length(ROTSn) >0) message(fxNa,"Argument 'ROTSn' is deprecated, please used argument 'multCorMeth' instead (like multCorMeth=c(ROTSn='10'))")
    if(length(lfdrInclude) >0) message(fxNa,"Argument 'lfdrInclude' is deprecated - IGNORING, please used argument 'multCorMeth' instead (like multCorMeth='lfdrInclude')")

    ## get compatible to old arguments lfdrInclude & ROTSn
    ROTSn <- NULL
    multCorMeth <- if(length(multCorMeth) <1) c("lfdr","FDR","means") else unique(c(multCorMeth, "means"))
    if(length(multCorMeth) ==1 && is.numeric(multCorMeth)) { 
       multCorMeth <- if(multCorMeth >1) c("lfdr", ROTSn=as.integer(multCorMeth), "means") else "lfdr"}
    
    if("ROTSn" %in% names(multCorMeth)) { ROTSn <- try(as.integer(multCorMeth["ROTSn"]), silent=TRUE)
      if(inherits(ROTSn, "try-error")) {ROTSn <- NULL; comp <- NULL }} 
  
    if("lfdr" %in% multCorMeth) { lfdrInclude <- TRUE
    } else if("lfdr" %in% names(multCorMeth)) { lfdrInclude <- try(as.logical(multCorMeth["lfdr"]), silent=TRUE)
        if(inherits(lfdrInclude, "try-error")) { lfdrInclude <- FALSE; multCorMeth <- multCorMeth[-which(names(multCorMeth)== "lfdr")] } }
    if(length(lfdrInclude) <1) lfdrInclude <- FALSE            # if for some reason whatsoever ...
        
    ## main
    isNA <- is.na(dat)
    chNA <- any(isNA)
    nNAmat <- matrix(0, nrow=nrow(dat), ncol=length(levels(gr)), dimnames=list(NULL,levels(gr)))
    seedNo <- as.integer(seedNo)[1]
    gr <- try(as.factor(gr)) 
    if(inherits(gr, "try-error")) message("+++++\n",fxNa," MAJOR PROBLEM with argument 'gr' !!  (possibly not sufficient level-names ?) \n+++++")
    
    if(length(avSdH) >1 && length(avSd) <1)  { avSd <- avSdH
      if(!silent) message(fxNa,"Using deprecated 'avSdH' as substitute of 'avSd', please adopt your code to use 'avSd' !!")
    } else if(length(avSdH) >1 && !silent) message(fxNa,"Argument 'avSdH' has been deprecated, 'avSd' is used instead")
    
    ## 1st pass
    if(debug) message(fxNa,"imputMethod= ",imputMethod,"  tRN1e")

    ## prepare useComparison & sep1 for later use
    if(length(sep1) !=1) sep1 <- wrMisc::getPWseparator(compNames=NULL, grp=gr, silent=!debug, callFrom=fxNa,debug=debug)
    if(length(useComparison) ==0) {
      useComparisonNa <- t(utils::combn(levels(gr), 2 ))
      dimnames(useComparisonNa) <- list(paste0(useComparisonNa[,1],sep1, useComparisonNa[,2]), c("samp","ref"))
      useComparison <- matrix(match(useComparisonNa, wrMisc::naOmit(unique(as.character(gr)))), ncol=2, dimnames=dimnames(useComparisonNa))   # refers to unique(gr) and NOT to levels(gr) !!
      #useCompRef <- wrMisc::indexGroupsFromPW(compNames=rownames(useComparisonNa), grp=gr)
    }
    
    if("none" %in% tolower(imputMethod)) {
      if(debug) {message(fxNa," tRN1f"); tRN1f <- list(dat=dat,gr=gr,useComparison=useComparison,imputMethod=imputMethod,annot=annot,seedNo=seedNo, retnNA=retnNA, avSd=avSd,ROTSn=ROTSn,lfdrInclude=lfdrInclude,multCorMeth=multCorMeth)}
      #nLoop <- 1
      chFin <- is.finite(dat)
      if(any(!chFin, na.rm=TRUE)) dat[which(!chFin)] <- NA   # need to replace Inf & -Inf by NA to avoid problems at limma::lmFit() 
      datI <- list(data=dat)
      datFi <- combineMultFilterNAimput(dat=dat, imputed=datI, grp=gr, useComparison=rownames(useComparison),annDat=annot, abundThr=stats::quantile(if(is.list(dat)) dat$quant else dat, 0.02,na.rm=TRUE), silent=silent, debug=debug,callFrom=fxNa)  # number of unique peptides unknown !
       
    } else {
      ## do first imputation befor applying filtering
      datI <- matrixNAneighbourImpute(dat, gr, imputMethod=imputMethod, retnNA=retnNA ,avSd=avSd, plotHist=plotHist, xLab=xLab, tit=tit, seedNo=seedNo, silent=!debug, callFrom=fxNa,debug=debug)
      if(debug) { message(fxNa,"Start combineMultFilterNAimput   tRN1g"); tRN1g <- list(dat=dat,datI=datI,gr=gr,useComparison=useComparison,imputMethod=imputMethod,annot=annot,seedNo=seedNo, retnNA=retnNA, avSd=avSd,ROTSn=ROTSn,lfdrInclude=lfdrInclude,multCorMeth=multCorMeth,plotHist=plotHist, xLab=xLab,tit=tit, seedNo=seedNo)}
      datFi <- combineMultFilterNAimput(dat=dat, imputed=datI, grp=gr, useComparison=useComparison,annDat=annot, abundThr=stats::quantile(if(is.list(dat)) dat$quant else dat, 0.02,na.rm=TRUE), silent=silent, debug=debug,callFrom=fxNa)  # number of unique peptides unknown !
      if(debug) {message(fxNa,"Done combineMultFilterNAimput   tRN1h")}     # done combineMultFilterNAimput
    }
    
    ## prepare for testing
    if(lfdrInclude) {
      chLfdr <- try(find.package("fdrtool"), silent=TRUE)
      if(inherits(chLfdr, "try-error")) { 
        message(fxNa,"Package 'fdrtool' NOT found ! Please install first from CRAN for calculating lfdr-values. Omitting (defaut) 'lfdr' option from argument 'multCorMeth' ..")
        lfdrInclude <- FALSE } }
    if(debug) {message(fxNa,"  tRN2"); tRN2 <- list(dat=dat,gr=gr,useComparison=useComparison,datFi=datFi,multCorMeth=multCorMeth,datI=datI,imputMethod=imputMethod)}
    
    if(debug) message(fxNa,"Start 1st moderTestXgrp()")
    pwComb <- NULL    # initailize

    ### ADJUST argument 'useComparison' to pwComb    (wr 29dec25)
    if(length(useComparison) >0) {
      ch1 <- any(is.na(useComparison) | useComparison=="")
      if(any(ch1)) {useComparison <- NULL; if(debug) message(fxNa,"NOTE : argument 'useComparison' contains NAs or empty fields - IGNORING !!")}
    }

    if(length(useComparison) >0) {
      ## assume useComparison is 2-column matrix OR character-vector assiging the levels of gr to be compared  (and replacing default pwComb)
      if(debug) message(fxNa,"  tRN2aa")
      if(length(dim(useComparison)) <2) {  # treat as character vector
        useComparison <- sub(" {0,3}- {0,3}","-", useComparison)    # remove spaces around separator
        aa <- strsplit(useComparison,"-")
        ch1 <- unique(unlist(aa))
        ch2 <- ch1 %in% levels(gr)
        if(all(!ch2)) { warning(fxNa,"Cannot relate these groups of 'useComparison' to 'gr' : ",wrMisc::pasteC(ch1[which(!ch2)], quoteC="'")," ..  ignoring")
          aa <- NULL
        } else {          
          if(any(!ch2)) { ch3 <- which(sapply(aa, function(x) all(x %in% unique(gr))))
            aa <- aa[ch3]   
            if(!silent) message(fxNa, " ",sum(!ch2)," elements not found in 'gr', eliminating ",length(aa)- length(ch3)," comparisons" ) 
          }
          if(length(aa) >0) {          
            chLe <- sapply(aa, length) ==2
            if(!all(chLe)) { warning(fxNa,"Having trouble understanding/relating 'useComparison' to 'gr', ignoring")         
              aa <- NULL }
          }    
        } 
        if(length(aa) >0) {
          ab <- sapply(aa, match, levels(gr))
          if(is.matrix(ab)) { comp <- t(ab); dimnames(comp) <- list(useComparison, c("samp","ref")) 
            if(debug) message(fxNa,"Successfully set custom comparisons (from character vector)") }        
          useComparison <- comp       # replace initial
        }        
      } else {
        ## check if 2-dim
        if(debug) message(fxNa," tRN2ac")
        if(ncol(useComparison) >1 && nrow(useComparison) >0) {
          useComparison <- useComparison[,1:2]  
          #chNum <- suppressWarnings(try(as.numeric(as.matrix(useComparison)), silent=TRUE))
          #if(any(is.na(chNum), inherits(chNum, "try-error"))) {   # convert to index
        }
      }
      pwComb <- useComparison           
    } else {
      if(length(pwComb) <1) { pwComb <- wrMisc::triCoord(length(levels(gr)))
        rownames(pwComb) <- paste(levels(gr)[pwComb[,1]],levels(gr)[pwComb[,2]], sep="-")   # might better first test for best separator
        if(debug) message(fxNa,"Aborting import/use of useComparison") } else if(debug) message(
        fxNa,"Sucessfully imported custom choice of comparisons from 'useComparison'  (or from  dat$sampleSetup$compareTo )")
    }   
    ## from now on rather use pwComb instead of useComparison    
    if(debug) {message(fxNa,"  tRN2b"); tRN2b <- list(dat=dat,gr=gr,datFi=datFi,multCorMeth=multCorMeth,datI=datI,pwComb=pwComb,useComparison=useComparison,sep1=sep1,chNA=chNA)}

    if(utils::packageVersion("wrMisc") < "2.0.1") { 
      pwComb2 <- if(length(dim(pwComb)) !=2) matrix(pwComb, nrow=1) else pwComb
      pwComb2 <- paste(levels(gr)[pwComb2[,1]], levels(gr)[pwComb2[,2]], sep=sep1)       
      # pwComb2 <- paste(levels(gr)[if(length(dim(pwComb))==2) pwComb[,1] else pwComb[1]], levels(gr)[if(length(dim(pwComb))==2) pwComb[,2] else pwComb[2]], sep=sep1)
      out <- suppressWarnings(wrMisc::moderTestXgrp(datFi$data, grp=gr, useComparison=pwComb2, limmaOutput=TRUE, addResults=multCorMeth, silent=silent, debug=debug, callFrom=fxNa))
    } else out <- suppressWarnings(wrMisc::moderTestXgrp(datFi$data, grp=gr, useComparison=pwComb, limmaOutput=TRUE, addResults=multCorMeth, silent=silent, debug=debug, callFrom=fxNa))   # can't do question specific filtering w/o explicit loop
    if(debug) {message(fxNa," tRN3"); tRN3 <- list(dat=dat,gr=gr,datFi=datFi,multCorMeth=multCorMeth,datI=datI,out=out,pwComb=pwComb,useComparison=useComparison,useComparisonNa=useComparisonNa,sep1=sep1,chNA=chNA,ROTSn=ROTSn,nLoop=nLoop)}
   
    chFDR <- names(out) =="FDR"
    if(any(chFDR, na.rm=TRUE)) names(out)[which(chFDR)] <- "BH"                 # rename $FDR to $BH

    if(length(useComparison) <1) rownames(pwComb) <- colnames(out$t)         # needed ??

    ## need to add $ROTS.p
    if(length(ROTSn)==1) { if(ROTSn >0 && !is.na(ROTSn)) {  
      chPa <- requireNamespace("ROTS", quietly=TRUE)
      if(!chPa) { message(fxNa,"Package 'ROTS' not found/installed (please install from Bioconductor), omitting argument 'ROTSn'")
        ROTSn <- 0 }
      } else ROTSn <- NULL }
    if(length(ROTSn)==1 && ROTSn >0) {
      ## this requires package ROTS
      if(debug) message(fxNa,"Start ROTS   n=",ROTSn)
      comp <- useComparison
      #comp <- wrMisc::triCoord(length(levels(gr)))
      #rownames(comp) <- paste(levels(gr)[comp[,1]], levels(gr)[comp[,2]],sep="-")
      tmRO <- matrix(nrow=nrow(datFi$data), ncol=nrow(comp))
      comPair <- matrix(unlist(strsplit(rownames(comp),"-")), ncol=nrow(comp))
      useCol <- apply(comPair, 2, function(x) gr %in% x)
      for(i in 1:nrow(comp)) tmRO[which(datFi$filt[,i]),i] <- ROTS::ROTS(datFi$data[which(datFi$filt[,i]), 
        which(useCol[,i])], groups=as.numeric(as.factor(gr[which(useCol[,i])]))-1, B=ROTSn)$pvalue        # K=500  
      out$ROTS.p <- tmRO
      if(length(ROTSn)==1 && ROTSn >0) {
        out$ROTS.BH <- apply(tmRO, 2, stats::p.adjust, method="BH") 
        if(lfdrInclude) out$ROTS.lfdr <- suppressWarnings(apply(tmRO, 2, wrMisc::pVal2lfdr)) }
    } 
    if(debug) {message(fxNa," tRN4"); tRN4 <- list(dat=dat,datI=datI, gr=gr,datFi=datFi,multCorMeth=multCorMeth,datI=datI,out=out,pwComb=pwComb,useComparison=useComparison,useComparisonNa=useComparisonNa,sep1=sep1,chNA=chNA,nLoop=nLoop,seedNo=seedNo,imputMethod=imputMethod,avSd=avSd,retnNA=retnNA)}
          
    ## subsequent rounds of NA-imputation  
    if(chNA && nLoop >1) { 
      if(debug) message(fxNa,"Subsequent rounds of NA-imputation   nLoop=",nLoop)
      pValTab <- tValTab <- array(NA_real_, dim=c(nrow(dat), if(length(dim(pwComb))==2) nrow(pwComb) else 1, nLoop))
      datIm <- array(NA_real_, dim=c(nrow(dat), ncol(dat), nLoop))
      datIm[,,1] <- datFi$data
      pValTab[,,1] <- out$p.value
      tValTab[,,1] <- out$t
      if(length(ROTSn)==1) if(ROTSn >0) {
        pVaRotsTab <- array(NA_real_, dim=c(nrow(dat), nrow(pwComb), min(10,nLoop)))
        pVaRotsTab[,,1] <- out$ROTS.p }
      if(debug) {message(fxNa," tRN5"); tRN5 <- list(out=out,dat=dat,datI=datI,gr=gr,datFi=datFi,multCorMeth=multCorMeth,nLoop=nLoop,imputMethod=imputMethod,seedNo=seedNo,retnNA=retnNA, avSd=avSd, NAneigLst=datI$NAneigLst,datIm=datIm,pValTab=pValTab,tValTab=tValTab)}  # 


      for(i in 2:nLoop) {
        ## the repeated NA-imputation & testing   
        if(length(seedNo)==1) seedNo <- seedNo +i 
      
        datX <- matrixNAneighbourImpute(dat, gr, imputMethod=imputMethod, seedNo=seedNo, retnNA=retnNA, avSd=avSd, NAneigLst=datI$NAneigLst, plotHist=FALSE, silent=TRUE, callFrom=fxNa)$data
        #datI <- matrixNAneighbourImpute(dat, gr, imputMethod=imputMethod, seedNo=seedNo, retnNA=retnNA ,avSd=avSd, plotHist=plotHist, xLab=xLab, tit=tit, seedNo=seedNo, silent=!debug, callFrom=fxNa,debug=debug)
        if(debug) {message(fxNa,"Passed matrixNAneighbourImpute()   in loop no ",i,"  tRN5b"); tRN5b <- list(dat=dat,gr=gr,datX=datX,out=out,datIm=datIm,pValTab=pValTab) }
        if(is.list(datX) && "data" %in% names(datX)) datX <- datX$data
  
        #1st round# datI <- matrixNAneighbourImpute(dat, gr, seedNo=seedNo, retnNA=retnNA ,avSd=avSd, plotHist=plotHist, xLab=xLab, tit=tit, silent=silent, callFrom=fxNa)
        #1st round# datFi <- combineMultFilterNAimput(dat=dat, imputed=datI, grp=gr, annDat=annot, abundThr=stats::quantile(dat,0.02,na.rm=TRUE), silent=silent, callFrom=fxNa) 
        #1st round# out <- wrMisc::moderTestXgrp(datFi$data, grp=gr, limmaOutput=TRUE, addResults="", silent=silent, callFrom=fxNa)       
        fitX <- limma::eBayes(limma::contrasts.fit(limma::lmFit(datX, out$design), contrasts=out$contrasts))
        datIm[,,i] <- datX
        pValTab[,,i] <- fitX$p.value
        tValTab[,,i] <- fitX$t
        if(debug) message(fxNa,"  tRN5c - ",i)
        if(length(ROTSn)==1) if(ROTSn >0 & i < min(10, nLoop)) {       # test using ROTS (TAKES MUCH TIME !!)
          for(i in 1:nrow(comp)) tmRO[which(datFi$filt[,i]),i] <- ROTS::ROTS(datFi$data[which(datFi$filt[,i]),which(useCol[,i])], groups=as.numeric(as.factor(gr[which(useCol[,i])])), B=ROTSn)$pvalue   # ,K=500  
          pVaRotsTab[,,i] <- tmRO } 
      }      # end loop of repetitive NA-imputation & testing
      
      ## define optimal sep (ie sep1)
      ## from wrProteo::combineMultFilterNAimput()
        sepX <- c("-","--","_","__"," - "," -- ","[[:upper:]]|[[:digit:]]vs[[:upper:]]|[[:digit:]]","-vs-","_vs_","--vs--","__vs__"," ","  ")  # candidates as sep, integrate "\." ?? (see sampNoDeMArrayLM())
        sepInGrp <- sapply(sepX, function(x) any(grepl(x, gr)))
        if(isTRUE(sepInGrp[7])) names(sepInGrp)[7] <- "vs"
        if(all(sepInGrp)) {
          sepX <- c(strrep("-", times=4:6), strrep("_", times=4:6))
          sepInGrp <- sapply(sepX, function(x) any(grepl(x, gr)))
          stop(fxNa,"Cannot find a suitable separator !")}
        if(length(useComparison) >0) {
          ## ?? moove block below (from moderTestXgrp() ) to separate fx ??  : check for suitable sep for grp, only once occuring in useComparison (if given)
          if(length(dim(useComparison))==0) {
            sepOnceInComp <- sapply(sepX, function(x) nchar(gsub(x,"", useComparison)) + nchar(x) ==nchar(useComparison))
            if(length(dim(sepOnceInComp))==2) sepOnceInComp <- colSums(sepOnceInComp) == nrow(sepOnceInComp)
            if(any(sepOnceInComp)) sep1 <- names(sepInGrp)[sepOnceInComp & !sepInGrp][1] else {
              if(!silent) message(fxNa,"Problem understanding 'useComparison', no reasonable group-separator found, ignoring 'useComparison'")
              sep1 <- useComparison <- NULL }
          } else sep1 <- names(sepInGrp)[which(!sepInGrp)[1]] 
        } else sep1 <- names(sepInGrp)[which(!sepInGrp)[1]]


      if(debug) { message(fxNa," tRN6"); tRN6 <- list(out=out,datFi=datFi,annot=annot,dat=dat,gr=gr,retnNA=retnNA, useComparison=useComparison,pwComb=pwComb,avSd=avSd,datI=datI,multCorMeth=multCorMeth, pValTab=pValTab,datIm=datIm,nLoop=nLoop,pwComb=pwComb,useComparison=useComparison,sep1=sep1 )  }
      
      ## better to combine values from interative procedure  before filtering
      if(length(dim(pValTab))==3) {
        pValTab <- apply(pValTab, c(1,2), mean, na.rm=TRUE)
        out$datImp <- as.matrix(apply(datIm, 1:2, mean, na.rm=TRUE))
      } else out$datImp <- datIm
      
      ## prepepare for adjusting to custom useComparison : adjCol gives index to use from complete combinatorics 
      refComb <- pwComb
      adjCol <- match(rownames(pwComb), colnames(datFi$filt))
      if(any(is.na(adjCol))) {
        tmp <- as.data.frame(strsplit(rownames(pwComb)[which(is.na(adjCol))], sep1))
        tmp <- paste(tmp[2,], tmp[1,],sep=sep1)
        adjCol[which(is.na(adjCol))] <- match(tmp, colnames(datFi$filt))
      } 
         
      ## propagate filtering results to p-values (disqualify to NA) : adjust pValTab
      if(any(!datFi$filt, na.rm=TRUE)) {    ## slow !?!
        ## need to adjust if not all pwComb tested
        useFi <- datFi$filt[,adjCol]
        pValTab[which(!useFi)] <- NA
      }   
      if(debug) { message(fxNa," tRN6b"); tRN6b <- list(out=out,datFi=datFi,annot=annot,dat=dat,gr=gr,retnNA=retnNA, useComparison=useComparison,pwComb=pwComb,adjCol=adjCol,avSd=avSd,datI=datI,multCorMeth=multCorMeth, pValTab=pValTab,datIm=datIm,nLoop=nLoop,pwComb=pwComb,useComparison=useComparison,sep1=sep1 )  }
         
      ## resume indiv rounds of imputation, optinal return details       
      rownames(out$datImp) <- if(is.null(rownames(dat))) rownames(annot) else rownames(dat)
      colnames(out$datImp) <- colnames(dat)
      if("tValTab" %in% multCorMeth) { out$tValArr <- tValTab
        out$pValArr <- pValTab }
      if("noLimma" %in% multCorMeth) out$simple.p.value <- out$p.value
      out$p.value <- if(length(dim(pValTab))==3) as.matrix(apply(pValTab, 1:2, stats::median, na.rm=FALSE)) else as.matrix(pValTab)
      out$t <- if(length(dim(tValTab))==3) as.matrix(apply(tValTab, 1:2, stats::median, na.rm=FALSE)) else as.matrix(tValTab)
      ###out$t <- as.matrix(apply(tValTab, 1:2, stats::median, na.rm=FALSE))
      colnames(out$p.value) <- colnames(out$t) <- rownames(pwComb)
      ## when converting t-value to p how to consider n due to nLoop ??
    } else {
      ## no iterative loop .. 
      out$datImp <- datFi$data
      filtZ <- which(!datFi$filt, arr.ind=TRUE)     
      if(any(!datFi$filt, na.rm=TRUE)) out$p.value[filtZ] <- out$t[filtZ] <- NA     # need to integrate filtering
    }
    
    
    out$annot <- annot
    out$filter <- datFi$filt[,adjCol]
    if(debug) {message(fxNa," tRN7"); tRN7 <- list(out=out,annot=annot,dat=dat,gr=gr,retnNA=retnNA, pwComb=pwComb,adjCol=adjCol,avSd=avSd,datI=datI,multCorMeth=multCorMeth,datFi=datFi, pValTab=pValTab,datIm=datIm,lfdrInclude=lfdrInclude ) }
  
    ## update dimnames of out$datImp
    dimnames(out$datImp) <- list(if(is.null(rownames(out$lods))) rownames(out$annot) else rownames(out$lods), colnames(dat))
    rownames(out$t) <- rownames(out$p.value) <- rownames(out$datImp)
    ## integrate column specific filtering
    #if(any(!datFi$filt, na.rm=TRUE)) out$p.value[which(!datFi$filt)] <- NA   # problem
    
    if(lfdrInclude) { out$lfdr <- as.matrix(suppressWarnings(apply(out$p.value, 2, function(x) suppressWarnings(wrMisc::pVal2lfdr(x, callFrom=fxNa)))))
      dimnames(out$lfdr) <- list(rownames(out$lods), if(length(dim(pwComb))==2) rownames(pwComb) else colnames(out$contrasts))
      if("noLimma" %in% multCorMeth) out$simple.lfdr <- if(ncol(out$simple.p.value) >1) as.matrix(apply(out$simple.p.value, 2, function(x) suppressWarnings(wrMisc::pVal2lfdr(x, callFrom=fxNa)))) else suppressWarnings(wrMisc::pVal2lfdr(out$simple.p.value))
    }
  
    if(any(c("FDR","BH") %in% multCorMeth, na.rm=TRUE)) { out$BH <- as.matrix(apply(out$p.value, 2, stats::p.adjust, method="BH"))
      dimnames(out$BH) <- list(rownames(out$lods), if(length(dim(pwComb))==2) rownames(pwComb) else colnames(out$contrasts))
      if("noLimma" %in% multCorMeth) out$simple.BH <- if(ncol(out$simple.p.value) >1) apply(out$simple.p.value, 2, stats::p.adjust, method="BH") else stats::p.adjust(out$simple.p.value, method="BH")
    } 
    if(debug) {message(fxNa," tRN8"); tRN8 <- list(out=out,annot=annot,dat=dat,gr=gr,retnNA=retnNA, pwComb=pwComb,avSd=avSd,datI=datI,multCorMeth=multCorMeth,datFi=datFi, pValTab=pValTab,datIm=datIm,lfdrInclude=lfdrInclude ) }
    if("BY" %in% multCorMeth) { out$BY <- as.matrix(apply(out$p.value, 2, stats::p.adjust, method="BY"))
      dimnames(out$BY) <- list(rownames(out$lods), if(length(dim(pwComb))==2) rownames(pwComb) else colnames(out$contrasts))}
    if(length(ROTSn)==1) if(ROTSn >0 && chNA && nLoop >1) {
      out$ROTS.p <- apply(pVaRotsTab, 1:2, stats::median, na.rm=TRUE)
      if(any(!datFi$filt, na.rm=TRUE)) out$ROTS.p[which(!datFi$filt)] <- NA    
      out$ROTS.BH <- as.matrix(apply(out$ROTS.p, 2, stats::p.adjust, method="BH"))
      dimnames(out$ROTS.BH) <- list(rownames(out$lods), if(length(dim(pwComb))==2) rownames(pwComb) else colnames(out$contrasts) )
      if(lfdrInclude) {out$ROTS.lfdr <- as.matrix(apply(out$ROTS.p, 2, function(x) suppressWarnings(wrMisc::pVal2lfdr(x, callFrom=fxNa))))
        dimnames(out$ROTS.lfdr) <- list(rownames(out$lods), if(length(dim(pwComb))==2) rownames(pwComb) else colnames(out$contrasts))} }
      if(debug) {message(fxNa," tRN9"); tRN9 <- list(out=out,annot=annot,dat=dat,gr=gr,retnNA=retnNA, avSd=avSd,datI=datI,multCorMeth=multCorMeth, pValTab=pValTab,datIm=datIm )    }
      out
  } else { warning(fxNa, msg)
    return(NULL) }
  }
    

