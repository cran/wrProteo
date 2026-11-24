#' Combine Multiple Filters On NA-Imputed Data
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
#' @param experSetup (list) optional esperimental setup (list with $ind, $pwGrpIndex, $sep and $pwGrpNa  as obtained using \code{getPairwiseSetup})
#' @param annDat (matrix or data.frame) annotation data (should match lines of 'dat')
#' @param abundThr (numeric) optional threshold filter for minimumn abundance
#' @param colRazNa (character) if razor peptides are used: column name for razor peptide count 
#' @param colTotNa (character) column name for total peptide count 
#' @param minSpeNo (integer) minimum number of specific peptides for maintaining proteins
#' @param minTotNo (integer) minimum total ie max razor number of peptides
#' @param maxGrpMiss (numeric) passd to \code{\link[wrMisc]{presenceFilt}} : at least 1 group has equal or less n NAs ;; at least 1 group has not more than this number of NAs (otherwise marke line as bad)
#' @param ratMaxNA (numeric) passd to \code{\link[wrMisc]{presenceFilt}} : at least 1 group has equal or less percent NAs ;; at least 1 group below this content of \code{NA} values
#' @param ratMaxTotNA (default NULL, numeric) passd to \code{\link[wrMisc]{presenceFilt}} : max ration of NAs in all samples per protein
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This function returns a list with $data (unmodified input from dat), $nNA, $randParam, $NAneigLst, $seed, $filt (the reommended filtering : TRUE for keeping value), and $annot  
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
#' datT6[8,1:3] <- NA
#' datT6b <- matrixNAneighbourImpute(datT6, gr=gl(2,3))
#' datT6c <- combineMultFilterNAimput(datT6, datT6b, grp=gl(2,3), abundThr=2)
#' head(datT6c$filt)
#' ## with custom choice of comparisons :
#' gr2 <- gl(2, 3, labels=c("A","B"))
#' datT6d <- combineMultFilterNAimput(datT6, datT6b, grp=gr2, useCom="B-A", abundThr=2)
#' head(datT6d$filt)
#' @export
combineMultFilterNAimput <- function(dat, imputed, grp, useComparison="all", experSetup=NULL, annDat=NULL, abundThr=NULL, colRazNa=NULL, colTotNa=NULL,
  minSpeNo=1, minTotNo=2, maxGrpMiss=9, ratMaxNA=0.5, ratMaxTotNA=NULL, silent=FALSE,debug=FALSE,callFrom=NULL){
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="combineMultFilterNAimput")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) { silent <- FALSE } else { debug <- FALSE }
  pwComb <- pwIndex <- pwGrpNa <- sep1 <- datFi <- NULL
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
  
  ## extract/define pairwise setup
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
  
  if(length(sep1) ==0) sep1 <- wrMisc::getPWseparator(grp=grp, silent=silent, debug=debug, callFrom=fxNa) 
  if(length(useComparison) ==0) useComparison <- paste(pwGrpNa[,1], pwGrpNa[,2], sep=sep1)  # default compare 1st & 2nd (alphabet soted) group

  ## extract & check of pwComb
  if(length(useComparison) !=0 && length(dim(useComparison))==2 && ncol(useComparison)==2) {
    pwComb <- if(is.numeric(useComparison)) experSetup$pwGrpNa else useComparison }
  if(length(pwComb)==0 || ncol(pwComb) !=2) {
    pwComb <- experSetup$pwGrpNa
    if(all(is.na(pwComb))) warning(fxNa,"PROBLEM : All 'experSetup$pwGrpNa' are NA !!")
    if(length(dim(pwComb)) !=2) pwComb <- matrix(pwComb, byrow=TRUE, ncol=2)
    if(is.numeric(pwComb)) warning(fxNa,"PROBLEM : Expecting 'experSetup$pwGrpNa' as matrix of group-names (character) anf NOT as index !!")
    ## possibly check that pwComb appears in grp 
    #pwComb <- experSetup$pwGrpIndex
  }
  
  if(debug) {message(fxNa,"cMF0e"); cMF0e <- list(dat=dat,imputed=imputed,grp=grp,useComparison=useComparison,experSetup=experSetup,annDat=annDat,abundThr=abundThr,pwComb=pwComb,sep1=sep1,maxGrpMiss=maxGrpMiss,ratMaxNA=ratMaxNA)} 
  
  
  datFi <- try(wrMisc::presenceFilt(dat, grp=grp, useComparison=useComparison, maxGrpMiss=maxGrpMiss, ratMaxNA=ratMaxNA, sep=sep1, silent=silent, debug=debug,callFrom=fxNa), silent=TRUE)

  if(inherits(datFi, "try-error")) stop(fxNa,"Failed to run presenceFilt")
  
  #pwComb <- useComparison   #
  #if(ncol(useComparison)) 
  if(debug) {message(fxNa,"cMF0f   at presenceFilt:  ",paste(colSums(datFi), collapse=" "),"  out of ",nrow(dat)); cMF0f <- list(dat=dat,imputed=imputed,grp=grp,useComparison=useComparison,colRazNa=colRazNa,datFi=datFi,pwComb=pwComb,experSetup=experSetup)}
  
  ## F1) filter  for max no of peptides (from $annot)
  if(length(colRazNa) >0 && length(annDat) >0) {
    razFilt <- razorNoFilter(annot=annDat, totNa=colTotNa, minRazNa=colRazNa, minSpeNo=minSpeNo, minTotNo=minTotNo, silent=silent, debug=debug, callFrom=fxNa)
    datFi[which(!razFilt),] <- rep(FALSE, ncol(datFi)) 
    if(debug) message(fxNa,"cMF0g   at razorNoFilter: ",paste(colSums(datFi), collapse=" "))    
  }

  ## F2) filter mostly low abundance (using imputed), see also .filterMinAv
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
  ## F3) low abundance filter
  if(length(abundThr)==1 && !is.na(abundThr)) {
    abundbyGrp2 <- grpMeans >= abundThr      # init passing min abund
    filtAbundByGrp <- apply(pwComb, 1, function(x) {rowSums(abundbyGrp2[, which(colnames(grpMeans) %in% x), drop=FALSE]) != 0})
  } else filtAbundByGrp <- matrix(TRUE, ncol=ncol(pwComb), nrow=nrow(grpMeans))
  

  # if(length(abundThr)==1 && is.numeric(abundThr)) {
  #   for(i in 1:nrow(pwComb)) {                                                # loop along all pair-wise questions => (update filter) datFi
  #     chLi <- all(grpMeans[ , if(is.numeric(pwComb)) pwComb[i,] else grpNa[match(pwComb[i,], grpNa)]] < abundThr)
  #     if(any(chLi, na.rm=TRUE)) datFi[which(chLi), i] <- FALSE }
  #   if(debug) message(fxNa,"cMF2b   at abundanceFilt: ",paste(colSums(datFi), collapse=" ")) }

  ## F4) future option - check if set of mostly imputed data higher than measured -> filter
  ## F5a) number of NAs per line & group
  nNAbyGroup <- wrMisc::rowGrpNA(dat, grp)
  tabGrp <- table(grp)
  if(debug) {message(fxNa,"cMF3"); cMF3 <- list(dat=dat,datFi=datFi,grpMeans=grpMeans,imputed=imputed,grp=grp,useComparison=useComparison,experSetup=experSetup,annDat=annDat,abundThr=abundThr,maxGrpMiss=maxGrpMiss,ratMaxNA=ratMaxNA,pwComb=pwComb,nNAbyGroup=nNAbyGroup,grpNa=grpNa,filtAbundByGrp=filtAbundByGrp)}  

  ## F5b) min nNA & %NA by pw question
  nNAbyGrp2 <- sapply(unique(as.character(pwComb)), function(x) rowSums(is.na(dat[, which(grp %in% x), drop=FALSE])))     # no of NAs in groups used in pw questions
  percNAbyGrp2 <- nNAbyGrp2/ matrix( rep(table(grp)[colnames(nNAbyGrp2)], each=nrow(nNAbyGrp2)), nrow=nrow(nNAbyGrp2))  # initial just %NA per group

  minnNAbyGrp <- apply(pwComb, 1, function(x) {aa <- nNAbyGrp2[, which(colnames(nNAbyGrp2) %in% x), drop=FALSE]; pmin(aa[, 1], aa[, 2])})            # observed min NA count  by pw group
  minpercNAbyGrp <- apply(pwComb, 1, function(x) {aa <- percNAbyGrp2[, which(colnames(percNAbyGrp2) %in% x), drop=FALSE]; pmin(aa[, 1], aa[, 2])})   # observed min NA%  by pw group
  ## apply thrsh  (T means to)   maxGrpMiss=2;  ratMaxNA=0.5   .. at least one pw group with 50% data AND at least one group with no more than 2 missing

  ## F6) combine filtering 
  imputed$filt <- minnNAbyGrp <= maxGrpMiss & minpercNAbyGrp <= ratMaxNA 


  ## F7) Combine: add abundance filter
  if(length(filtAbundByGrp) !=0 && !all(filtAbundByGrp)) imputed$filt <- imputed$filt & filtAbundByGrp
  
  ## F8)  combine with optional total %NA (PAIRWISE) => allow removing lines with NAs
  if(length(ratMaxTotNA) ==1 && is.finite(ratMaxTotNA) && ratMaxTotNA >= 0) {
    ## need  maxnNAbyGrp
    maxnNAbyGrp <- apply(pwComb, 1, function(x) {aa <- nNAbyGrp2[, which(colnames(nNAbyGrp2) %in% x), drop=FALSE]; pmax(aa[, 1], aa[, 2])})            # observed min NA count  by pw group
    imputed$filt <- imputed$filt & maxnNAbyGrp ==ratMaxTotNA 
  }
  
  ## F8b) combine with optional total %NA  => allow removing ALL lines with ANY NAs (even for pw comparisons where all data present)
  #if(length(ratMaxTotNA) ==1 && is.finite(ratMaxTotNA) && ratMaxTotNA >= 0) imputed$filt <- imputed$filt & (rowSums(is.na(dat))/ ncol(dat) <= ratMaxTotNA)
  

  # ## filter for each pw question
  # for(i in 1:nrow(pwComb)) {
  #   critNAGrp <- tabGrp[if(is.numeric(pwComb)) pwComb[i,] else grpNa[match(pwComb[i,], grpNa)]]
    
  #   nNAbyGrp2 <- nNAbyGroup[,if(is.numeric(pwComb)) pwComb[i,] else grpNa[match(pwComb[i,], grpNa)]]  # limit to groups used
  #   percNAbyGrp2 <- nNAbyGrp2/ matrix( rep(table(grp)[colnames(nNAbyGrp2)], each=nrow(nNAbyGrp2)), nrow=nrow(nNAbyGrp2))

  #       minNNA2 <- do.call(pmin, as.data.frame(nNAbyGrp2))    # presumed faster than apply
  #       minPercNA2 <- do.call(pmin, as.data.frame(percNAbyGrp2))

  #   ## filter for min NA count (ie at least 1 group at this level)
  #   ## argumnet : maxGrpMiss .. at least 1 group has equal or less n NAs ... (numeric) passd to \code{\link[wrMisc]{presenceFilt}} : at least 1 group has not more than this number of NAs (otherwise marke line as bad)
  #   minNNA3 <- apply(pwComb, 1,  function(x) do.call(pmin, as.data.frame(nNAbyGrp2[,match(x, grpNa)])) <= maxGrpMiss )


  #   ## filter for min NA percentage (ie at least 1 group at this level)
  #   ## argumnet : ratMaxNA .. at least 1 group has equal or less percent NAs
  #   minNNA3 <- apply(pwComb, 1,  function(x) do.call(pmin, as.data.frame(percNAbyGrp2[,match(x, grpNa)])) <= ratMaxNA )   # each pw comp in new column, TRUE => use

  #   ## comnine
  #   filtMinNA <- minNNA3 & minNNA3

  #   #old# critNAGrp <- critNAGrp/2 -0.1
  #   critNAGrp <- critNAGrp/2 -0.1
  

  #     #re-check ?# potentially filter when min 50% of data NA
  #   #chLi <- nNAbyGroup[,pwComb[i,]] > matrix(rep(critNAGrp, each=nrow(grpMeans)), ncol=2)   # use imputed$nNA; return T when need to filter
  #   ## check if any line & grp is exceeding its critNAGrp
  #   chLi <- nNAbyGroup[,if(is.numeric(pwComb)) pwComb[i,] else grpNa[match(pwComb[i,], grpNa)]] > matrix(rep(critNAGrp, each=nrow(grpMeans)), ncol=2)   # use imputed$nNA; return T when need to filter
  #   if(any(chLi)) {
  #     chLi2 <- cbind(chLi[,1] & grpMeans[,if(is.numeric(pwComb)) pwComb[i,] else grpNa[match(pwComb[i,], grpNa)]] > grpMeans[,if(is.numeric(pwComb)) pwComb[i,] else grpNa[match(pwComb[i,], grpNa)]], chLi[,2] & grpMeans[,if(is.numeric(pwComb)) pwComb[i,] else grpNa[match(pwComb[i,], grpNa)]] > grpMeans[,if(is.numeric(pwComb)) pwComb[i,] else grpNa[match(pwComb[i,], grpNa)]]) # is T if bad
  #     datFi[,i] <- datFi[,i] & !chLi[,1] & !chLi[,2] }
  # }                
  # if(debug) message(fxNa,"   at NA > mean:   ",wrMisc::pasteC(colSums(datFi)))    
  # imputed$filt <- datFi
  if(!is.null(annDat)) imputed$annot <- annDat    
  imputed }
          
