#' Volcano-plot (Statistical Test Outcome versus Relative Change)   
#'
#' This type of plot is very common in high-throughput biology, see \href{https://en.wikipedia.org/wiki/Volcano_plot_(statistics)}{Volcano-plot}.
#' Basically, this plot allows comparing the outcome of a statistical test to the differential of the group means (ie log fold-change),
#' 
#' In high-throughput biology data are typically already transformed to log2 and thus, the 'M'-value represents a relative change.
#' Besides, output from statistical testing by \code{\link[wrMisc]{moderTest2grp}} or \code{\link[wrMisc]{moderTestXgrp}} can be directly read to produce Volcano plots for diagnostic reasons.
#' Please note, that plotting a very number of points in transparency (eg >10000) may take several seconds.
#'	 
#' @param Mvalue (numeric or matrix) data to plot; M-values are typically calculated as difference of log2-abundance values and 'pValue' the mean of log2-abundance values;
#'   M-values and p-values may be given as 2 columsn of a matrix, in this case the argument \code{pValue} should remain NULL 
#' @param pValue (numeric, list or data.frame) if \code{NULL} it is assumed that 2nd column of 'Mvalue' contains the p-values to be used
#' @param useComp (integer, length=1) choice of which of multiple comparisons to present in \code{Mvalue} (if generated using \code{moderTestXgrp()})  
#' @param filtFin (matrix or logical) The data may get filtered before plotting: If \code{FALSE} no filtering will get applied; if matrix of \code{TRUE}/\code{FALSE} it will be used as optional custom filter, otherwise (if \code{Mvalue} if an \code{MArrayLM}-object eg from limma) a default filtering based on the \code{filtFin} element will be applied 
#' @param ProjNa (character) custom title
#' @param FCthrs (numeric) Fold-Change threshold (display as line) give as Fold-change and NOT log2(FC), default at 1.5, set to \code{NA} for omitting
#' @param FdrList (numeric) FDR data or name of list-element
#' @param FdrThrs (numeric) FDR threshold (display as line), default at 0.05, set to \code{NA} for omitting 
#' @param FdrType (character) FDR-type to extract if \code{Mvalue} is 'MArrayLM'-object (eg produced by from \code{moderTest2grp} etc);
#'   if \code{NULL} it will search for suitable fields/values in this order : 'FDR','BH',"lfdr" and 'BY'
#' @param subTxt (character) custom sub-title
#' @param grayIncrem (logical) if \code{TRUE}, display overlay of points as increased shades of gray
#' @param col (character) custom color(s) for points of plot (see also \code{\link[graphics]{par}})
#' @param pch (integer) type of symbol(s) to plot (default=16) (see also \code{\link[graphics]{par}}) 
#' @param compNa (character) names of groups compared
#' @param batchFig (logical) if \code{TRUE} figure title and axes legends will be kept shorter for display on fewer splace  
#' @param cexMa (numeric) font-size of title, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}})
#' @param cexLa (numeric) size of axis-labels, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}}) 
#' @param limM (numeric, length=2) range of axis M-values 
#' @param limp (numeric, length=2) range of axis FDR / p-values
#' @param annotColumn (character) column names of annotation to be extracted (only if \code{Mvalue} is \code{MArrayLM}-object containing matrix $annot).
#'   The first entry (typically 'SpecType') is used for different symbols in figure, the second (typically 'GeneName') is used as prefered text for annotating the best points (if \code{namesNBest} allows to do so.)
#' @param annColor (character or integer) colors for specific groups of annoatation (only if \code{Mvalue} is \code{MArrayLM}-object containing matrix $annot)
#' @param cexPt (numeric) size of points, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}})
#' @param cexSub (numeric) size of subtitle, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}})
#' @param cexTxLab (numeric) size of text-labels for points, as expansion factor (see also \code{cex} in \code{\link[graphics]{par}})
#' @param namesNBest (integer or character) number of best points to add names in figure; if 'passThr' all points passing FDR and FC-filtes will be selected; 
#'   if the initial object \code{Mvalue} contains a list-element called 'annot' the second of the column specified in argument \code{annotColumn} will be used as text
#' @param NbestCol (character or integer) colors for text-labels of best points
#' @param sortLeg (character) sorting of 'SpecType' annotation either ascending ('ascend') or descending ('descend'), no sorting if \code{NULL}
#' @param NaSpecTypeAsContam (logical) consider lines/proteins with \code{NA} in Mvalue$annot[,"SpecType"] as contaminants (if a 'SpecType' for contaminants already exits)
#' @param useMar (numeric,length=4) custom margings (see also \code{\link[graphics]{par}})
#' @param returnData (logical) optional returning data.frame with (ID, Mvalue, pValue, FDRvalue, passFilt) 
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @param debug (logical) additional messages for debugging 
#' @return MA-plot only
#' @seealso (for PCA) \code{\link[wrGraph]{plotPCAw}})
#' @examples
#' library(wrMisc)
#' set.seed(2005); mat <- matrix(round(runif(900),2), ncol=9)
#' rownames(mat) <- paste0(rep(letters[1:25],each=4), rep(letters[2:26],4))
#' mat[1:50,4:6] <- mat[1:50,4:6] + rep(c(-1,1)*0.1,25)
#' mat[3:7,4:9] <- mat[3:7,4:9] + 0.7
#' mat[11:15,1:6] <- mat[11:15,1:6] - 0.7
#' ## assume 2 groups with 3 samples each
#' gr3 <- gl(3,3,labels=c("C","A","B"))
#' tRes2 <- moderTest2grp(mat[,1:6], gl(2,3), addResults = c("FDR","means"))
#' # Note: due to the small number of lines only FDR chosen to calculate 
#' VolcanoPlotW2(tRes2)
#' ## Add names of points passing custom filters
#' VolcanoPlotW2(tRes2, FCth=1.3, FdrThrs=0.2, namesNBest="passThr")
#'
#' ## assume 3 groups with 3 samples each
#' tRes <- moderTestXgrp(mat, gr3, addResults = c("FDR","means"))
#' # Note: due to the small number of lines only FDR chosen to calculate 
#' VolcanoPlotW2(tRes)
#' VolcanoPlotW2(tRes, FCth=1.3, FdrThrs=0.2)
#' VolcanoPlotW2(tRes, FCth=1.3, FdrThrs=0.2, useComp=2)
#'  
#' @export
VolcanoPlotW2 <- function(Mvalue, pValue=NULL, useComp=1, filtFin=NULL, ProjNa=NULL, FCthrs=NULL, FdrList=NULL, FdrThrs=NULL, FdrType=NULL,
  subTxt=NULL, grayIncrem=TRUE, col=NULL, pch=16, compNa=NULL, batchFig=FALSE, cexMa=1.8, cexLa=1.1, limM=NULL, limp=NULL,
  annotColumn=c("SpecType","GeneName","EntryName","Accession","Species","Contam"), annColor=NULL, cexPt=NULL, cexSub=NULL, 
  cexTxLab=0.7, namesNBest=NULL, NbestCol=1, sortLeg="descend", NaSpecTypeAsContam=TRUE, useMar=c(6.2,4,4,2), returnData=FALSE, callFrom=NULL, silent=FALSE,debug=FALSE) {
  ## MA plot
  ## optional arguments for explicit title in batch-mode
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="VolcanoPlotW2")
  message("++ NOTE : OLD VERSION !! ++ please use VolcanoPlotW() from package wrGraph")
  opar <- graphics::par(no.readonly=TRUE) 
  on.exit(graphics::par(opar$mar)) 
  on.exit(graphics::par(opar$cex.main)) 
  on.exit(graphics::par(opar$cex.lab)) 
  on.exit(graphics::par(opar$las)) 
  plotByFdr <- TRUE  #
  namesIn <- c(deparse(substitute(Mvalue)), deparse(substitute(pValue)), deparse(substitute(filtFin)))
  basRGB <- c(0.3,0.3,0.3)           # grey
  fcRGB <- c(1,0,0)                  # red        for points passing  FC filt line
  splNa <- annot <- ptType <- colPass <- ptBg <- grpMeans <- pcol <- FDRvalue <- NULL      # initialize
  multiComp <- TRUE      # initialize  
  if(debug) silent <- FALSE
  if(length(Mvalue) <1) message(" nothing to do, 'Mvalue' seems to be empty !") else  {
    ## data seem valid to make MAplot
    if(length(cexTxLab) <0) cexTxLab <- 0.7
    if("MArrayLM" %in% class(Mvalue)) {
      ## try working based on MArrayLM-object (Mvalue)
      ## initial check of  useComp
      if(length(useComp) >1) { useComp <- wrMisc::naOmit(useComp)[1]
        if(!silent) message(fxNa," argument 'useComp' should be integer of length=1; using only 1st entry") }
      if(length(useComp) <1) { useComp <- 1
        if(!silent) message(fxNa," argument 'useComp' invalid, setting to 1") }
      if(length(pValue) <1) {
        ## no spearate pValues provided : extract from MArrayLM-object (Mvalue)
        pcol <- wrMisc::naOmit(match(c("p.value","pvalue","pval","p"), tolower(names(Mvalue))))
        if(length(pcol) >0) pValue <- Mvalue[[pcol[1]]] else stop("can't find suitable element for p-Values from MArrayLM-object")
        if(length(dim(pValue)) >0) if(colnames(pValue)[1]=="(Intercept)" & ncol(pValue) >1) {
          ## extract 2nd col if result from wrMisc::moderTest2grp() 
          pNa <- rownames(pValue)
          pValue <- as.numeric(pValue[,2])
          names(pValue) <- pNa
          multiComp <- FALSE 
        } else {
          ## select corresponding of multiple comparisons
          if(useComp > ncol(pValue)) { useComp <- 1
            if(!silent) message(fxNa," argument 'useComp' for pValues invalid or too high; reset to 1") }
          names(useComp) <- colnames(pValue)[useComp]
          pNa <- rownames(pValue)
          pValue <- as.numeric(pValue[,useComp])
          if(length(pNa) >0) names(pValue) <- pNa }
      }  ## .. otherwise spearate pValue was provided
      ## look for M-values (need to create if not available - using useComp checked when extracting pValue)
      Melem <- wrMisc::naOmit(match(c("mvalues","mvalue","m"), tolower(names(Mvalue))))      # which list-element
      Fcol <- wrMisc::naOmit(match(if(length(FdrType) <1) c("fdr","bh","lfdr","by","p.value") else FdrType, tolower(names(Mvalue))))    # needed for matching means to pair-wise names and for extracting FDR values
      if(length(Fcol) <1) stop("Can't find elment suitable for statistical values", if(length(FdrType)==1) c(" to '",FdrType,"'")) else Fcol <- Fcol[1]
      if(!silent & length(FdrType) <1) message(fxNa,"Using element '",names(Mvalue)[Fcol],"' as FDR-values for plot")
      if("lfdr" %in% names(Mvalue)[Fcol]) {   # switch from directly plotting FDR-values to uncorrected p-values
        plotByFdr <- FALSE
      }
      
      ## look for group-means & identify column association to current question/pairwise comparison
      if("means" %in% names(Mvalue)) {
        ## identify sammple-groups to comparsison(s) - needed lateron
        pairwCol <- .sampNoDeMArrayLM2(Mvalue, useComp, lstMeans="means", lstP=Fcol,callFrom=fxNa,silent=silent) 
        grpMeans <- cbind(mean1=Mvalue$means[,pairwCol[1]], mean2=Mvalue$means[,pairwCol[2]])  
        ## are all group-means needed (for exporting) ??
      } else warning("Could not find suitable field '$means' in '",namesIn[1],"'")    
      
      if(length(Melem) >0) {            ## M-values are available
        Mvalue$Mval <- Mvalue[[Melem]]
        if(length(dim(Mvalue$Mval)) >0) if(ncol(Mvalue$Mval) >1) {          
          if(useComp[1] > ncol(Mvalue$Mval)) { if(!silent) message(fxNa," 'useComp' is too high, set to 1"); useComp <- 1 }
          Mvalue$Mval <- Mvalue$Mval[,useComp]}
      } else {                               # need to construct M-values based on means
        if("means" %in% names(Mvalue)) {
          ## construct Mvalue based on means (only one/current pairwise comparison needed)
          Mvalue$Mval <- grpMeans[,2] - grpMeans[,1]                           
          Melem <- which(names(Mvalue)=="Mval")                # update
        } else stop("Can't construct M-values since suitable field '$means' missing in '",namesIn[1],"' !")    
      }
      ## now one can check if 'pValue' & Mvalue match
      chPM <- length(pValue) >0 & length(as.numeric(pValue)) == length(as.numeric(Mvalue[[Melem]]))
      if(!chPM) {
        if(length(pcol) >0) {            # pVal avilable in Mvalue 
          if(length(as.numeric(Mvalue[[pcol[1]]][,useComp])) ==length(as.numeric(Mvalue[[Melem]]))) {  # pVal from Mvalue seems to fit => use
             pValue <- as.numeric(Mvalue[[pcol[1]]][,useComp])          
          } else stop("'pValue' & 'Mvalue' don't match")  
        } else stop("'pValue' & 'Mvalue' don't match (no field in MArrayLM-object available)") 
      }       
      ## extract FDR from MArrayLM-object 
      if(length(FdrList) <1) {
        ## no explicit pValue, try to extract from MArrayLM-object (Mvalue)
        if(length(Fcol) >0) { FDRvalue <- Mvalue[[Fcol[1]]]
          ## extract 2nd col if result from wrMisc::moderTest2grp() 
          if(length(dim(FDRvalue)) >0) if(colnames(FDRvalue)[1]=="(Intercept)" & ncol(FDRvalue) >1) {
            pNa <- rownames(FDRvalue)
            FDRvalue <- as.numeric(FDRvalue[,2])
            names(FDRvalue) <- pNa 
            multiComp <- FALSE }
        } else {
          FDRvalue <- stats::p.adjust(pValue) 
          if(!silent) message(fxNa,"No FDR data found, generating BH-FDR")}
        ## need to find corresponding of multiple comparisons
        if(length(dim(FDRvalue)) >0) {
          pNa <- rownames(FDRvalue)
          if(ncol(FDRvalue) >1) { if(useComp > ncol(FDRvalue)) { useComp <- 1
            if(!silent) message(fxNa," argument 'useComp' for FDRvalues invalid or too high; reset to 1") }
          FDRvalue <- as.numeric(FDRvalue[,useComp])
          if(length(pNa) >0) names(FDRvalue) <- pNa }}
      }
      ## recuperate filtering - if present, but only when no custom filtering provided
      if(length(filtFin) <1 | identical(filtFin, FALSE)) {
        Fcol <- wrMisc::naOmit(match(c("filtfin","filter","filt","finfilt"), tolower(names(Mvalue))))
        filtFin <- if(length(Fcol) >0) Mvalue[[Fcol[1]]] else rep(TRUE,length(pValue))
        if(length(dim(filtFin)) >1) filtFin <- filtFin[,useComp]        
      }
      ## recuperate $annot if present and use for symbol
      if("annot" %in% names(Mvalue)) {
        useAnnCol <- match(annotColumn, colnames(Mvalue$annot))      
        if(!is.na(useAnnCol[1])) {                         # annotation (for multiple groups) exists
          ptType <- Mvalue$annot[,useAnnCol[1]]            # SpecType
          chNA <- is.na(ptType)
          ## associate NAs from 'SpecType' in ptType with conta ?
          if(NaSpecTypeAsContam) {
            chConta <- tolower(ptType) %in% c("contaminant","contam","conta","cont")
            if(any(chConta)) ptType[which(is.na(ptType))] <- unique(ptType[which(chConta)])[1]}          
          if(any(is.na(ptType))) ptType[which(chNA)] <- "NA" 
          if(length(pch) < length(pValue) & length(unique(wrMisc::naOmit(ptType))) >1) {
            if(length(pch) >1 & !silent) message(fxNa," (invalid pch) using default 'pch' oriented by $annot and starting from 15")
            pch <- 14 + as.integer(as.factor(ptType))
          } 
          useAnnCol <- wrMisc::naOmit(useAnnCol)
          annot <- Mvalue$annot[,useAnnCol] 
          if(annotColumn[1] %in% colnames(annot)) annot[,annotColumn[1]] <- ptType           
      } }
      if(length(pch)==1) pch <- rep(as.integer(pch), length(pValue))
      
      ## recuperate M values (& dismiss rest of MArrayLM-object)  
      Mvalue <- Mvalue$Mval
      if(length(dim(Mvalue)) >1) { MNa <- rownames(Mvalue)
        Mvalue <- as.numeric(Mvalue)
        if(length(MNa) >0) names(Mvalue) <- MNa }
      ## additional check for length 
      chpM <- length(Mvalue)==length(pValue)  
      if(!chpM & !silent) message(fxNa,"trouble ahead ? p- and M- values have different length !!  (M=",length(Mvalue)," vs p=",length(pValue),")")

      ## done with extracing MArrayLM-object    
      if(!silent) message(fxNa,"Successfully extracted  ",length(Mvalue)," Mvalues and  ",length(pValue)," pValues", if(length(annot) >0) c(" plus anotation"))      
    } else {
      ## thus argument 'Mvalue' is not 'MArrayLM'-object
      ## ... case of explicit pValue argument
      if(length(pValue) <1) stop(" argument 'pValue' is required (if 'Mvalue' not 'MArrayLM'-type object) !")
      if(length(dim(pValue)) >1) if(ncol(pValue) >1) {
        if(!silent) message(fxNa," Note, ",namesIn[2]," has ",ncol(pValue)," columns, using last column")
        pNa <- rownames(pValue)
        pValue <- as.numeric(pValue[,ncol(pValue)] )  
        names(pValue) <- pNa} 
      FDRvalue <- if(length(FdrList) <1) NULL else FdrList
    }
    
    ## need to introduce -log10 to pValue
    chNA <- is.na(pValue)
    if(all(chNA)) stop(fxNa," All p-values are NA, nothing to draw !")
    pValue <- -log10(pValue) 
    ## check for (same) order, adjust Mvalue & pValue according to names
    chNa <- list(MNa=if(length(dim(Mvalue)) >1) rownames(Mvalue) else names(Mvalue),
      pNa=if(length(dim(pValue)) >1) rownames(pValue) else names(pValue))
    nIni <- c(M=length(Mvalue),p=length(pValue))
    if(length(chNa$MNa) >0 & length(chNa$pNa) >0) {        # ie both have names, so one can match names
      if(!identical(chNa$MNa,chNa$pNa)) {
        matchNa <- wrMisc::naOmit(match(chNa$MNa,chNa$pNa))
        if(length(matchNa) <1) stop("Both 'Mvalue' and 'pValue' have names, but none of them match !!")
        pValue <- pValue[matchNa]
        Mvalue <- wrMisc::naOmit(Mvalue[match(names(pValue),names(Mvalue))])
      } } else {
        if(length(Mvalue) != length(pValue)) stop("p- and M- values have different length, but no names to match !!  (M=",length(Mvalue)," vs p=",length(pValue),")")
      }
    if(length(grpMeans) <1) grpMeans <- matrix(rep(NA,2*length(Mvalue)), ncol=2, dimnames=list(names(Mvalue),c("mean1","mean2")))
   
    ## start creating merged data for plot (& export)
    merg <- if(length(annot) >0) data.frame(ID=NA, grpMeans, Mvalue=Mvalue, pValue=pValue, FDR=if(length(FDRvalue) >0) FDRvalue else rep(NA,length(pValue)), 
      filtFin=filtFin, annot, pch=pch, stringsAsFactors=FALSE) else {
      data.frame(ID=NA, grpMeans, Mvalue=Mvalue, pValue=pValue, FDR=FDRvalue, filtFin=filtFin, pch=pch, stringsAsFactors=FALSE) }
    if(length(names(Mvalue)) >0) merg[,1] <- names(Mvalue) else {if(length(names(pValue)) >0) merg[,1] <- names(pValue)}
    ## replace NA in 'SpecType' by 'NA'

    if(annotColumn[1] %in% colnames(merg)) { chNa <- is.na(merg[,annotColumn[1]])     # replace NAs in col "SpecType" by "NA"
      if(any(chNa)) merg[which(chNa),annotColumn[1]] <- "NA"
    } else { merg <- cbind(merg, rep(1,nrow(merg)))            # add colum for 'SpecType'
      colnames(merg)[ncol(merg)] <- annotColumn[1] }
    
    ## adjust col & pch
    if(!any(c(1,length(Mvalue)) %in% length(pch))) {
      if(!silent) message(fxNa,"argument 'pch' should be either length=1 or correspond to length of data, reset to default=16")
      pch <- 16 }
    if(length(col) >1 & length(col) <length(Mvalue)) {
      if(!silent) message(fxNa,"argument 'col' should be either length=1 or correspond to length of data, reset to default=NULL")
      col <- NULL }

    ## prepare/integrate FILTERING
    if(length(filtFin) >0) {
      ## if filtFin is matrix use each line with min 1 instance of TRUE,
      if(length(dim(filtFin)) >1) filtFin <- as.logical(as.matrix(filtFin)[,useComp])    # use rows with >= 1 TRUE
      if(length(names(filtFin)) >0) {
        matchNa <- wrMisc::naOmit(match(rownames(merg), names(filtFin)))       
        if(length(matchNa)==nrow(merg)) merg[,"filtFin"] <- filtFin[matchNa]
      } else if(length(filtFin)==nrow(merg)) merg[,"filtFin"] <- filtFin        # no proof that order of filtFin is correct
    } else filtFin <- rep(TRUE, nrow(merg)) 
    if(debug) message(fxNa," ++ DONE extracting columns : ",wrMisc::pasteC(colnames(merg),quo="'"))
    
    ## apply filtering
    msg <- " data provided in 'Mvalue' and 'pValue' "
    if(!silent & nrow(merg) < round(length(Mvalue)/10)) message(" .. note : less than 10% of",msg," were matched") else {
      if(!silent & nrow(merg) < length(Mvalue)/2) message(" .. NOTE : less than 50% of",msg," were matched !!")}
    if(debug) message(msg," were matched to ",nrow(merg)," common entries")    
    ## apply filtering (keep all lines where at least one condition passes)
    if(length(filtFin) >0 & !identical(filtFin, FALSE)) {                         #  use filtering provided
      if(sum(filtFin) >0 & sum(filtFin) < nrow(merg)) { 
        whFilt <- which(merg$filtFin)
        if(length(pch) >1) pch <- pch[whFilt]
        if(length(col) >1) col <- col[whFilt]
        merg <- merg[whFilt,]
        if(!silent) message(fxNa," filtered (based on 'filtFin') from ",length(filtFin)," to  ",nrow(merg)," lines")
      }
    } else filtFin <- rep(TRUE, nrow(merg))
    
    ## sort merg, so that legend always gets constructed the same order, ascending ('ascend') or descending ('descend')    
    sortLeg <- if(identical(sortLeg,"ascend")) FALSE else {if(identical(sortLeg,"descend")) TRUE else NULL}
    if(length(sortLeg) ==1 & annotColumn[1] %in% colnames(merg)) merg <- merg[order(merg[,annotColumn[1]], decreasing=sortLeg),]     
     
    ## update ..
    nIDco <- sum(c("ID","nredID","uniqID") %in% colnames(merg))                   #  number of heading columns in 'merg'
    Mvalue <- as.numeric(if("Mvalue" %in% colnames(merg)) merg[,"Mvalue"] else merg[,nIDco+1])
    pValue <- as.numeric(if("pValue" %in% colnames(merg)) merg[,"pValue"] else {
      if(length(dim(Mvalue)) >0) merg[,ncol(Mvalue) +nIDco +1] else merg[,nIDco+2]})
    if("Lfdr" %in% colnames(merg)) FdrList <- merg[,"Lfdr"] else {
      if("lfdr" %in% colnames(merg)) FdrList <- merg[,"lfdr"]}
    pch <- merg[,"pch"]            # update
    ptType <- if(annotColumn[1] %in% colnames(merg)) merg[,annotColumn[1]] else rep(1,nrow(merg))     # update "SpecType"

    ## prepare for  plotting
    if(is.null(cexSub)) cexSub <- cexLa +0.05  
    xLab <- "M-value (log2 fold-change)"
    tit1 <- paste(c(if(!batchFig) c(ProjNa, if(!is.null(ProjNa)) ": ","Volcano-plot"),
      if(!is.null(compNa)) c(compNa[1]," vs ",compNa[2])), collapse=" ")    # but what title if batchFig=NULL & compNa=NULL -> only "Volcano-plot"
    if(length(FCthrs) <1) FCthrs <- 1.5 
    if(length(FdrThrs) <1) FdrThrs <- 0.05 

    ## count no of passing
    passFC <- if(length(FCthrs) ==1 & !any(is.na(FCthrs))) abs(merg[,"Mvalue"]) > log2(FCthrs) else merg[,"filtFin"]      ## convert FCthrs to log2
    passFdr <- if(length(FdrThrs) ==1 & !any(is.na(FdrThrs))) {merg[,"FDR"] <= FdrThrs} else merg[,"filtFin"]
    passAll <- merg[,"filtFin"] & passFC & passFdr
    chNA <- is.na(passAll)                              # passFdr may contain NAs
    if(any(chNA)) passAll[which(chNA)] <- FALSE 
    if(debug) message(fxNa,"  ",sum(passFC,na.rm=TRUE)," passing FCthrs ; ",sum(passFdr,na.rm=TRUE)," passing FdrThrs ; combined ",sum(passAll,na.rm=TRUE))
    ## color for points passing filter
    if(length(col) >0) if(length(col) != nrow(merg)) { col <- NULL
      if(!silent) message(fxNa," invalid entry for 'col', should be of length=",nrow(merg),", resetting to default")}
    if(length(col) <1) {
      alph <- sort(c(0.14, round(0.6/log10(length(Mvalue)),2), 0.8))[2]       # alph <- round(12/sqrt(nrow(eBayesLst$pValue)),2)
      alph2 <- sort(c(round(7/(5 +sum(passAll)^0.7),2), alph,0.9))[2]                   # for points passing thresholds
      useCol <- if(grayIncrem) grDevices::rgb(0.35,0.35,0.35,alph) else grDevices::rgb(0.7,0.7,0.7)  # basic color
      useCex <- if(length(cexPt) >0) cexPt else max(round(0.8 + 2/(1 +sum(filtFin, na.rm=TRUE))^0.28,2), 1.1)
      chCol <- unique(merg[, annotColumn[1]])      # check how many different colors may be needed
      chNaC <- is.na(chCol)
      if(any(chNaC)) chCol[which(chNaC)] <- "NA"    
      if(length(annColor) >0) {colPass <- annColor} else if(length(chCol) >4) {
        colPass <- cbind(red=c(141,72,90,171, 220,253,244,255), green=c(129,153,194,221, 216,174,109,0), blue=c(194,203,185,164, 83,97,67,0))       
        colPass <- grDevices::rgb(red=colPass[,1], green=colPass[,2], blue=colPass[,3], alph2, maxColorValue=255)
        if(length(chCol) >8) { colPass <- c(colPass, rep(colPass[8], length(chCol) -8))
          if(!silent) message(fxNa," > 8 different groups found, using 8th color after 7th group")}
      } else colPass <- grDevices::rgb(c(0.95,0.2,0,0.75), c(0.15,0.2,0.9,0.35), c(0.15,0.95,0,0.8), alph2)    # red, blue, green, purple (luminosity adjusted) 
      useCol <- rep(useCol[1], nrow(merg))         # fuse basic gray to colors for different types

      ## assign color for those passing
      if(any(passAll)) useCol[which(passAll)] <- colPass[if(length(unique(merg[which(passAll),annotColumn[1]])) >1) .levIndex(merg[which(passAll),annotColumn[1]]) else rep(1,sum(passAll))]  # assign colors for those passing 
    } else useCol <- col
    ## adjust fill color for open symbols
    chPch <- pch %in% c(21:25)
    if(any(chPch)) { ptBg <- useCol
      ptBg[which(chPch)] <- useCol[which(chPch)]    # background color for filled symbols
      useCol[which(chPch)] <- 1                     # contour as black
    }

    ## main graphic
    graphics::par(mar=c(6.5,4,4,2), cex.main=cexMa, las=1)
    ## rather directly plot FDR
    graphics::plot(Mvalue, if(plotByFdr) -1*log10(merg[,"FDR"]) else merg[,"pValue"], pch=pch, cex=useCex, main=tit1, 
      ylab=if(plotByFdr) "- log10 FDR" else "- log10 p-value (uncorrected)", col=useCol, xlab=xLab, cex.lab=cexLa, xlim=limM,ylim=limp, pt.bg=ptBg)    
    
    sTxt <- if(length(subTxt) ==1) subTxt else { if(multiComp) paste0(if(length(names(useComp)) >0) names(useComp) else c("useComp=",useComp),"; ",collapse="")}
    sTxt <- paste0(sTxt,"n=",length(Mvalue),
      if(!all(is.na(c(FCthrs,FdrThrs)))) paste(";",sum(passAll, na.rm=TRUE),"(color) points passing",
        if(!is.na(FCthrs)) paste0("(FCthr=", as.character(FCthrs),", ") else " (", paste0("FdrThrs=",as.character(FdrThrs),")")))
    graphics::mtext(sTxt,cex=0.75,line=0.2)
    if(!all(is.na(c(FCthrs,FdrThrs)))) { 
      if(debug) message(fxNa," n=",length(Mvalue),"  FCthrs=",as.character(FCthrs),"  filt.ini=", sum(filtFin, na.rm=TRUE),
        "  passAll=",sum(passAll,na.rm=TRUE)," ; range Mva ",wrMisc::pasteC(signif(range(Mvalue,na.rm=TRUE),3))," ;  alph=",alph,"  useCex=",useCex,"  alph2=",alph2)
      graphics::abline(v=c(-1,1)*(log2(FCthrs) + diff(graphics::par("usr")[1:2])/500), col=grDevices::rgb(0.87,0.72,0.72), lty=2) }
    if(sum(passFdr) >0) {
      if(plotByFdr) {
        graphics::abline(h=-1*log10(max(merg[passAll,"FDR"], na.rm=TRUE)) -diff(graphics::par("usr")[3:4])/400, col=grDevices::rgb(0.87,0.72,0.72), lty=2)
      } else { 
        graphics::mtext("Note, that FDR and p-value may not correlate perfectly, thus points may appear at good p-value but finally don't get retained",line=-1.4,cex=0.7)
        pRa <- range(merg[which(passFdr),"pValue"], na.rm=TRUE)
        graphics::abline(h=pRa[1] +diff(graphics::par("usr")[3:4])/400, col=grDevices::rgb(0.87,0.72,0.72), lty=2) }}
    
    ## add names to best points
    if(length(namesNBest) >0) { 
      if(identical(namesNBest,"passThr") | identical(namesNBest,"signif")) namesNBest <- sum(passAll) 
      if(!is.integer(namesNBest)) namesNBest <- try(as.integer(namesNBest))
      if(namesNBest >0 & any(passAll)) {      
        useLi <- if(any(!passAll)) which(passAll) else 1:nrow(merg)
        tmP <- as.numeric(merg[useLi,"pValue"])
        names(tmP) <- rownames(merg)[useLi]
        ## look for more informative names to display
        if(length(annot) >0) {
          proNa <- annot[match(names(tmP), rownames(annot)), annotColumn[2]]   # normally 'Description'
          chNa <- is.na(proNa)
          if(!all(chNa)) names(tmP)[which(!chNa)] <- proNa[which(!chNa)]
        }        
        useL2 <- order(tmP, decreasing=TRUE)[1:min(namesNBest,sum(passAll))]      
        xOffs <- signif(diff(graphics::par("usr")[1:2])/170,3)
        yOffs <- signif(diff(graphics::par("usr")[3:4])/90,3)
        noNa <- if(is.null(names(tmP[useL2]))) 1:length(tmP) else which(is.na(names(tmP)[useL2]))
        if(length(noNa) >0 & all(annotColumn %in% colnames(merg))) names(tmP)[useL2[noNa]] <- merg[useLi[useL2[noNa]], wrMisc::naOmit(match(annotColumn[-1], colnames(merg)))[1]]
        if(length(NbestCol) <1) NbestCol <- 1
        if(is.null(names(tmP[useL2]))) {if(!silent) message(fxNa," No names available for displaying names of best in plot")
        } else graphics::text(Mvalue[useLi[useL2]] +xOffs, yOffs -1*log10(merg[useLi[useL2], if(plotByFdr) "FDR" else "pValue"]),
          names(tmP)[useL2], cex=cexTxLab, col=NbestCol, adj=0) 
      }
    }              

    ## legend (if multiple symbols)
    pch[which(is.na(pch))] <- -2
    ch1 <- unique(pch)
    if(length(ch1) >1) {
      legInd <- which(!duplicated(merg[which(passAll), annotColumn[1]], fromLast=FALSE))
      legPch <- pch[which(passAll)[legInd]]
      legCol <- useCol[which(passAll)[legInd]]
      legBg <- ptBg[which(passAll)[legInd]]
      if(alph2 <1) {legCol <- substr(legCol,1,7); legBg <- substr(legBg,1,7)}  # reset to no transparency
      legLab <- merg[which(passAll)[legInd], annotColumn[1]]
      chNa <- is.na(legLab)
      if(any(chNa)) legLab[chNa] <- "NA"
      legOr <- if(length(legLab) >1) order(legLab) else 1   # not used so far 
      legLoc <- wrGraph::checkForLegLoc(cbind(Mvalue, pValue), sampleGrp=legLab, showLegend=FALSE)
      legCex <- stats::median(c(useCex,cexTxLab,1.2), na.rm=TRUE)
      graphics::legend(legLoc$loc, legend=legLab, col=legCol, text.col=1, pch=legPch, if(length(ptBg) >0) pt.bg=ptBg, cex=legCex, pt.cex=1.2*legCex, xjust=0.5, yjust=0.5)  # as points
    }

  ## export results
  if(returnData) {
    merg <- merg[,-1*c(1,ncol(merg))]        # remove col 'ID' 'redundant' & 'pch'
    annCo <- wrMisc::naOmit(match(annotColumn, colnames(merg)))
    if(length(annCo) >0) cbind(merg[,annCo],  merg[,-annCo]) else merg }
  } }
     
#' @export
.sampNoDeMArrayLM2 <- function(MArrayObj, useComp, groupSep="-",lstMeans="means",lstP="BH",silent=FALSE,callFrom=NULL) {
  ## locate sample index from index or name of pair-wise comparisons in list or MArrayLM-object
  fxNa <- wrMisc::.composeCallName(callFrom, newNa=".sampNoDeMArrayLM")
  errMsg <- c("argument 'MArrayObj' is ","empty","doesn't contain the list-element needed  ('",lstMeans,"') !")
  if(length(MArrayObj) <1) stop(errMsg[1:2])
  if(length(MArrayObj[[lstMeans]]) <1) stop(errMsg[-2])
  if(length(colnames(MArrayObj[[lstMeans]])) <1)  stop(" problem with 'MArrayObj$lstMeans' (does not contain matrix of means)")
  if(ncol(MArrayObj[[lstMeans]])==2) {             # only 2 mean-values, no other choice, don't need to try matching anything
    if(!identical(as.character(useComp),"1") & !silent) message(fxNa,"Only 2 columns of mean-values available, can't interpret properly 'useComp=",useComp,"'")
    out <- 1:2
  } else {
    if(length(lstP) <0) stop(" 'lstP' is empty !")
    if(length(colnames(MArrayObj[[lstP]])) <1)  stop(" problem with 'MArrayObj' (does not contain matrix of p-values)")
    ## convert/locate names to index 
    if(is.character(useComp) & length(grep("[[:alpha:]]",useComp)) >0) useComp <- wrMisc::naOmit(match(useComp, MArrayObj[[lstP]] ))
    if(length(useComp) <1) stop("argument 'useComp' is empty or can't locate in comparison-names")
    ## main
    out <- if(ncol(MArrayObj[[lstP]])==2 & colnames(MArrayObj[[lstP]])[1] =="(Intercept)") 1:2 else {
      wrMisc::matchSampToPairw(grpNa=colnames(MArrayObj[[lstMeans]]), pairwNa=colnames(MArrayObj[[lstP]])[useComp], sep=groupSep,silent=silent,callFrom=fxNa)}
    if(length(useComp)==1) out <- as.integer(out)}
  out }
    
.levIndex <- function(dat,asSortedLevNa=FALSE) {
  ## transform levels into index; should get integrated to wrMisc
  if(asSortedLevNa) out <- as.integer(as.factor(dat)) else {
    out <- dat
    levU <- wrMisc::naOmit(unique(out))   # levels in orig order (non-alpahbetical)
    for(i in 1:length(levU)) out[which(out==levU[i])] <- i
    out <- as.integer(out)}
  out }
   
