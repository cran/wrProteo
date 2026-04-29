#' Combine Multiple Filters On NA-imputed Data
#'
#' In most omics data-analysis one needs to employ a certain number of filtering strategies to avoid getting artifacts to the step of statistical testing.
#' This function takes both the origial data and the NA-imputed data to create several different filters (number of NAs, min abundance, ..) and finally combines them.
#' The filter-component aiming to take away the least abundant values (using the imputede data) can be fine-tuned by the argument \code{abundThr}.
#' 
#' @details
#' The main filtering step compares the means for each group and line (from argument \code{imputed}), at least one group-mean has to be superior the threshold (\code{abundThr}, 
#' based on hypothesis that if all conditions represent extrememy low measures their differential may not be determined with certainty).
#' 
#' The filter addressing the number of missing values (\code{NA}) uses the original data (argument \code{dat}), the arguments \code{colTotNa}, \code{minSpeNo} and \code{minTotNo} 
#' are used during this step. Basically, this step allows defining a minimum content of 'real' (ie non-\code{NA}) values for further considering the measurements as reliable.
#' This part uses internally \code{\link[wrMisc]{presenceFilt}} for filtering elevated content of \code{NA} per line.
#' 
#' Finally, this function combines both filters (as matrix of \code{FALSE} and \code{TRUE}) on NA-imputed and original data 
#' and retruns a vector of logical values if corresponding lines passe all filter criteria.
#'  
#' @param dat (matrix or data.frame) main data (may contain \code{NA})
#' @param imputed (character)  same as 'dat' but with all \code{NA} imputed
#' @param grp (character or factor) define groups of replicates (in columns of 'dat')
#' @param useComparison (character or matrix) argument allowing to specify which pairwise comparions sould be performed, default \code{useComparison=NULL} will run all pairwise comparisons;
#'   may be character as result of combining two group-names (from argument \code{grp}) (eg 'A-B', beware, the separator may not appear inside group-names) 
#'   or 2-column matrix of group-names (thus, without separator; names combined with separator may be shown as rownames) or 2-column matrix of indexes to sorted group-names;
#'   default setting of \code{useComparison='all'} will select all possible pairwise combinations
#' @param experSetup (list) optional esperimental setup (list with $ind, $pwGrpIndex,$sep and $pwGrpNa  as obtained using \code{getPairwiseSetup})
#' @param annDat (matrix or data.frame) annotation data (should match lines of 'dat')
#' @param abundThr (numeric) optional threshold filter for minimumn abundance
#' @param colRazNa (character) if razor peptides are used: column name for razor peptide count 
#' @param colTotNa (character) column name for total peptide count 
#' @param minSpeNo (integer) minimum number of specific peptides for maintaining proteins
#' @param minTotNo (integer) minimum total ie max razor number of peptides
#' @param maxGrpMiss (numeric) passd to \code{\link[wrMisc]{presenceFilt}} : at least 1 group has not more than this number of NAs (otherwise marke line as bad)
#' @param ratMaxNA (numeric) passd to \code{\link[wrMisc]{presenceFilt}} : at least 1 group below this content of \code{NA} values
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This function returns a list with $data, $nNA, $randParam, $NAneigLst, $seed, $filt, and $annot  
#' @seealso \code{\link[wrMisc]{presenceFilt}} 
#' @examples
#' set.seed(2013)
#' datT6 <- matrix(round(rnorm(300)+3,1), ncol=6,
#'   dimnames=list(paste0("li",1:50), letters[19:24]))
#' datT6 <- datT6 +matrix(rep(1:nrow(datT6),ncol(datT6)), ncol=ncol(datT6))
#' datT6[6:7,c(1,3,6)] <- NA
#' datT6[which(datT6 < 11 & datT6 > 10.5)] <- NA
#' datT6[which(datT6 < 6 & datT6 > 5)] <- NA
#' datT6[which(datT6 < 4.6 & datT6 > 4)] <- NA
#' datT6b <- matrixNAneighbourImpute(datT6, gr=gl(2,3))
#' datT6c <- combineMultFilterNAimput(datT6, datT6b, grp=gl(2,3), abundThr=2)
#' head(datT6c$filt)
#' ## with custom choice of comparisons :
#' gr2 <- gl(2, 3, labels=c("A","B"))
#' datT6d <- combineMultFilterNAimput(datT6, datT6b, grp=gr2, useCom="B-A", abundThr=2)
#' head(datT6d$filt)
#' @export
combineMultFilterNAimput <- function(dat, imputed, grp, useComparison="all", experSetup=NULL, annDat=NULL, abundThr=NULL, colRazNa=NULL, colTotNa=NULL, minSpeNo=1, minTotNo=2, maxGrpMiss=1, ratMaxNA=0.8, silent=FALSE,debug=FALSE,callFrom=NULL){
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="combineMultFilterNAimput")
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  if(!isTRUE(silent)) silent <- FALSE
  pwComb <- pwIndex <- pwGrpNa <- sep1 <- NULL
  ## checks
  if(length(dat)==0 || length(dim(dat)) !=2 || any(dim(dat) < c(1,2))) stop(fxNa,"Invalid input for argument 'dat' : Must be matrix or data.frame with min 2 columns !")
  #if(length(imputed) >0 && is.list(imputed) && "data" %in% names(imputed)) imputed <- imputed$data
  
  if(!is.list(imputed) || length(dim(imputed$data)) !=2) imputed <- list(data=imputed, nNA=sum(is.na(dat)))
  if(length(imputed$data)==0 || !identical(dim(imputed$data), dim(dat))) stop(fxNa,"Invalid input for argument 'imputed' : Must be matrix or data.frame of same dimensions as 'dat' !")
  if(length(grp)==0 && length(experSetup) > 0 && length(experSetup$grp) >0) grp <- experSetup$grp 
  if(!is.factor(grp)) { grp <- try(as.factor(grp), silent=TRUE)
    if(inherits(grp, "try-error") || length(grp) != ncol(dat)) stop(fxNa,"Invalid input for argument 'grp' : Must be factor with at least 2 levels where each element corresponds to one column fo 'dat' !") } 
  if(length(levels(grp)) <2) stop(fxNa,"Can't run function when 'grp' has less than 2 level ")
  if(debug) {message(fxNa,"cMF00"); cMF00 <- list(dat=dat,imputed=imputed,grp=grp,useComparison=useComparison,experSetup=experSetup,annDat=annDat,abundThr=abundThr,pwComb=pwComb,sep1=sep1)}

  if(length(experSetup) >0) {
    if(length(experSetup$sep)==1) sep1 <- experSetup$sep
    if(length(useComparison)==0 && length(sep1) ==0) sep1 <- experSetup$sep 
  }  
  if(debug) {message(fxNa,"cMF0b");  cMF0b <- list(useComparison=useComparison, grp=grp, experSetup=experSetup, sep1=sep1)}
  
  if(length(experSetup)==0) {
    experSetup <- wrMisc::convPairwiseSetup(useComparison=useComparison, grp=grp, experSetup=experSetup, sep=sep1, silent=silent, debug=debug, callFrom=fxNa)
    ## note $grp from convPairwiseSetup() not useful (is rather levels(grp) )
    if(length(sep1)==0 && "sep" %in% names(experSetup))  sep1 <- experSetup$sep
    if("pwGrpNa" %in% names(experSetup)) pwGrpNa <- experSetup$pwGrpNa
    if(length(grp)==0 && length(experSetup$grp)==ncol(dat)) grp <- experSetup$grp
    useComparison <- if("concat" %in% names(experSetup)) experSetup$concat else paste(experSetup$pwGrpNa[,1], sep=sep1, experSetup$pwGrpNa[,2])      
    
  } else {   # possibly from getPairwiseSetup()
    if("ind" %in% names(experSetup))  pwIndex <- experSetup$index
    if("pwGrpNames" %in% names(experSetup)) pwGrpNa <- experSetup$pwGrpNames else if("pwGrpNa" %in% names(experSetup)) pwGrpNa <- experSetup$pwGrpNa
    if("sep" %in% names(experSetup))  sep <- experSetup$sep
    if("ind" %in% names(experSetup))  pwIndex <- experSetup$index    
  }
  if(debug) {message(fxNa,"cMF0d");  cMF0d <- list(useComparison=useComparison, grp=grp, experSetup=experSetup, sep1=sep1)}
  grpNa <- wrMisc::naOmit(unique(sort(grp)))
  
  if(length(sep1)==0) sep1 <- wrMisc::getPWseparator(grp=grp, silent=silent, debug=debug, callFrom=fxNa) 
  if(length(useComparison)==0) useComparison <- paste(pwGrpNa[,1], pwGrpNa[,2], sep=sep1)
  if(debug) {message(fxNa,"cMF0e"); cMF0e <- list(dat=dat,imputed=imputed,grp=grp,useComparison=useComparison,experSetup=experSetup,annDat=annDat,abundThr=abundThr,pwComb=pwComb,sep1=sep1)} 
  datFi <- try(wrMisc::presenceFilt(dat, grp=grp, useComparison=useComparison, maxGrpMiss=maxGrpMiss, ratMaxNA=ratMaxNA, sep=sep1, silent=silent, debug=debug,callFrom=fxNa), silent=TRUE)

  if(inherits(datFi, "try-error")) stop(fxNa,"Failed to run presenceFilt")
  
  #pwComb <- useComparison   #
  #if(ncol(useComparison)) 
  if(debug) {message(fxNa,"cMF0f   at presenceFilt:  ",paste(colSums(datFi), collapse=" "),"  out of ",nrow(dat)); cMF0f <- list(dat=dat,imputed=imputed,grp=grp,useComparison=useComparison,colRazNa=colRazNa,datFi=datFi,pwComb=pwComb,experSetup=experSetup)}
  
  if(length(colRazNa) >0 && length(annDat) >0) {
    razFilt <- razorNoFilter(annot=annDat, totNa=colTotNa, minRazNa=colRazNa, minSpeNo=minSpeNo, minTotNo=minTotNo, silent=silent, debug=debug, callFrom=fxNa)
    datFi[which(!razFilt),] <- rep(FALSE, ncol(datFi)) 
    if(debug) message(fxNa,"cMF0g   at razorNoFilter: ",paste(colSums(datFi), collapse=" "))    
  }

  if(length(pwComb)==0 || ncol(pwComb) !=2) {
    pwComb <- experSetup$pwGrpNa
    if(length(dim(pwComb)) !=2) pwComb <- matrix(pwComb, byrow=TRUE, ncol=2)
    #pwComb <- experSetup$pwGrpIndex
  }
  
  ## filter mostly low abundance (using imputed), see also .filterMinAv
  grpMeans <- wrMisc::rowGrpMeans(imputed$data, grp)
  if(debug) {message(fxNa,"cMF1"); cMF1 <- list(dat=dat,imputed=imputed,grpMeans=grpMeans,datFi=datFi,grp=grp,useComparison=useComparison,experSetup=experSetup,annDat=annDat,abundThr=abundThr,pwComb=pwComb)}
  if(any(!(colnames(grpMeans) == colnames(imputed$nNA))) && !silent) message(fxNa,"Problem with order of columns of imputed$nNA !?")
  ## determine optimal separator
  ## ?? moove block below (from moderTestXgrp() ) to separate fx ??  : check for suitable sep for grp, only once occuring in useComparison (if given)
  # if(length(useComparison) >0 && length(dim(useComparison))==2 ) {
  #   pwComb <- useComparison
  # } else {
  #   pwComb <- wrMisc::indexGroupsFromPW(compNames=if(length(dim(useComparison))==2) rownames(useComparison) else useComparison, grp=unique(grp), includeGrp=FALSE, potSep=sep1, silent=silent, debug=debug, callFrom=fxNa)
  #   # sep1 <- pwComb$sep            # needed ?
  #   # pwCombNa <- pwComb$GrpNames   # needed ?
  #   # pwComb <- pwComb$ind          # needed ?
  # }
  #deprecated# pwComb <- wrMisc::indexPairWise(grp, useComparison=useComparison, silent=silent, debug=debug, callFrom=fxNa)
  
  
  if(debug) {message(fxNa,"cMF2"); cMF2 <- list(dat=dat,imputed=imputed,grp=grp,useComparison=useComparison,experSetup=experSetup,annDat=annDat,abundThr=abundThr,pwComb=pwComb,grpNa=grpNa,grpMeans=grpMeans)}
  
  if(length(abundThr)==1 && is.numeric(abundThr)) {
    for(i in 1:nrow(pwComb)) {                                                # loop along all pair-wise questions => (update filter) datFi
      chLi <- all(grpMeans[,if(is.numeric(pwComb)) pwComb[i,] else grpNa[match(pwComb[i,], grpNa)]] < abundThr)
      if(any(chLi)) datFi[which(chLi), i] <- FALSE }
    if(debug) message(fxNa,"cMF2b   at abundanceFilt: ",paste(colSums(datFi), collapse=" ")) }
  ## check if set of mostly imputed data higher than measured -> filter
  ## number of NAs per line & group
  nNAbyGroup <- wrMisc::rowGrpNA(dat, grp)
  if(debug) {message(fxNa,"cMF3"); cMF3 <- list(dat=dat,datFi=datFi,grpMeans=grpMeans,imputed=imputed,grp=grp,useComparison=useComparison,experSetup=experSetup,annDat=annDat,abundThr=abundThr,pwComb=pwComb,nNAbyGroup=nNAbyGroup,grpNa=grpNa)}
  
  tabGrp <- table(grp)
  for(i in 1:nrow(pwComb)) {
    critNAGrp <- tabGrp[if(is.numeric(pwComb)) pwComb[i,] else grpNa[match(pwComb[i,], grpNa)]]
    critNAGrp <- critNAGrp/2 -0.1  
      #re-check ?# potentially filter when min 50% of data NA
    #chLi <- nNAbyGroup[,pwComb[i,]] > matrix(rep(critNAGrp, each=nrow(grpMeans)), ncol=2)   # use imputed$nNA; return T when need to filter
    chLi <- nNAbyGroup[,if(is.numeric(pwComb)) pwComb[i,] else grpNa[match(pwComb[i,], grpNa)]] > matrix(rep(critNAGrp, each=nrow(grpMeans)), ncol=2)   # use imputed$nNA; return T when need to filter
    if(any(chLi)) {
      chLi2 <- cbind(chLi[,1] & grpMeans[,if(is.numeric(pwComb)) pwComb[i,] else grpNa[match(pwComb[i,], grpNa)]] > grpMeans[,if(is.numeric(pwComb)) pwComb[i,] else grpNa[match(pwComb[i,], grpNa)]], chLi[,2] & grpMeans[,if(is.numeric(pwComb)) pwComb[i,] else grpNa[match(pwComb[i,], grpNa)]] > grpMeans[,if(is.numeric(pwComb)) pwComb[i,] else grpNa[match(pwComb[i,], grpNa)]]) # is T if bad
      datFi[,i] <- datFi[,i] & !chLi[,1] & !chLi[,2] }
  }                
  if(debug) message(fxNa,"   at NA > mean:   ",wrMisc::pasteC(colSums(datFi)))    
  imputed$filt <- datFi
  if(!is.null(annDat)) imputed$annot <- annDat    
  imputed }
          
