#' Combine Multiple Filters On NA-imputed Data
#'
#' In most omics data-analysis one needs to employ a certain number of filtering strategies to avoid getting artifacts to the step of statistical testing.
#' \code{combineMultFilterNAimput} takes on one side the origial data and on the other side NA-imputed data to create several differnet filters and to finally combine them.
#' A filter aiming to take away the least abundant values (using the imputede data) can be fine-tuned by the argument \code{abundThr}. 
#' This step compares the means for each group and line, at least one grou-mean has to be > the threshold (based on hypothesis 
#' that if all conditions represent extrememy low measures their diffrenetial may not be determined with certainty).
#' In contrast, the filter addressing the number of missing values (\code{NA}) uses the original data, the arguments \code{colTotNa},\code{minSpeNo} and \code{minTotNo} 
#' are used at this step. Basically, this step allows defining a minimum content of 'real' (ie non-NA) values for further considering the measurements as reliable.
#' This part uses internally \code{\link[wrMisc]{presenceFilt}} for filtering elevated content of \code{NA} per line.
#' Finally, this function combines both filters (as matrix of \code{FALSE} and \code{TRUE}) on NA-imputed and original data 
#' and retruns a vector of logical values if corresponding lines passe all filter criteria.
#'  
#' @param dat (matrix or data.frame) main data (may contain \code{NA})
#' @param imputed (character)  same as 'dat' but with all \code{NA} imputed
#' @param grp (character or factor) define groups of replicates (in columns of 'dat')
#' @param annDat (matrix or data.frame) annotation data (should match lines of 'dat')
#' @param abundThr (numeric) optional threshold filter for minimumn abundance
#' @param colRazNa (character) if razor peptides are used: column name for razor peptide count 
#' @param colTotNa (character) column name for total peptide count 
#' @param minSpeNo (integer) minimum number of specific peptides for maintaining proteins
#' @param minTotNo (integer) minimum total ie max razor number of peptides
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This function returns a vector of logical values if corresponding line passes filter criteria  
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
#' 
#' @export
combineMultFilterNAimput <- function(dat, imputed, grp, annDat=NULL, abundThr=NULL, colRazNa=NULL, colTotNa=NULL, minSpeNo=1, minTotNo=2, silent=FALSE,debug=FALSE,callFrom=NULL){
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="combineMultFilterNAimput")
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  if(!isTRUE(silent)) silent <- FALSE
  datFi <- wrMisc::presenceFilt(dat, grp=grp, maxGrpM=1, ratMa=0.8, silent=silent, callFrom=fxNa)
  if(debug) message(fxNa,"   at presenceFilt:  ",paste(colSums(datFi),collapse=" "),"  out of ",nrow(dat))
  if(length(colRazNa) >0 & length(annDat) >0) {
    razFilt <- razorNoFilter(annot=annDat, totNa=colTotNa, minRazNa=colRazNa, minSpeNo=minSpeNo, minTotNo=minTotNo, silent=silent,callFrom=fxNa)
    datFi[which(!razFilt),] <- rep(FALSE,ncol(datFi)) 
    if(debug) message(fxNa,"   at razorNoFilter: ",paste(colSums(datFi),collapse=" "))    
    }
  ## filter mostly low abundance (using imputed), see also .filterMinAv
  grpMeans <- wrMisc::rowGrpMeans(imputed$data,grp)
  if(any(!(colnames(grpMeans) == colnames(imputed$nNA)))) message(fxNa," Problem with order of columns of imputed$nNA !?")
  pwComb <- wrMisc::triCoord(ncol(grpMeans))

  if(is.numeric(abundThr) & length(abundThr)==1) {
    for(i in 1:nrow(pwComb)) {                                                # loop along all pair-wise questions => (update filter) datFi
      chLi <- grpMeans[,pwComb[i,1]] < abundThr & grpMeans[,pwComb[i,2]] < abundThr
      if(any(chLi)) datFi[which(chLi),i] <- FALSE}
    if(debug) message(fxNa,"   at abundanceFilt: ",paste(colSums(datFi),collapse=" ")) }
  ## check if set of mostly imputed data higher than measured -> filter
  ## number of NAs per line & group
  nNAbyGroup <- wrMisc::rowGrpNA(dat,grp)
  
  for(i in 1:nrow(pwComb)) {
    critNAGrp <- table(grp)[colnames(grpMeans)[pwComb[i,]]]
    critNAGrp <- critNAGrp/2 -0.1  
      #re-check ?# potentially filter when min 50% of data NA
    chLi <- nNAbyGroup[,pwComb[i,]] > matrix(rep(critNAGrp, each=nrow(grpMeans)), ncol=2)   # use imputed$nNA; return T when need to filter
    if(any(chLi)) {
      chLi2 <- cbind(chLi[,1] & grpMeans[,pwComb[i,1]] > grpMeans[,pwComb[i,2]], chLi[,2] & grpMeans[,pwComb[i,2]] > grpMeans[,pwComb[i,1]]) # is T if bad
      datFi[,i] <- datFi[,i] & !chLi[,1] & !chLi[,2] }
  }                
  if(debug) message(fxNa,"   at NA> mean:   ",wrMisc::pasteC(colSums(datFi)))    
  imputed$filt <- datFi
  if(!is.null(annDat)) imputed$annot <- annDat    
  imputed }
    
