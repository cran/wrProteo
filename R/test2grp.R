#' t-test each line of 2 groups of data
#'
#' \code{test2grp} performs t-test on two groups of data using \href{https://bioconductor.org/packages/release/bioc/html/limma.html}{limma},
#' this is a custom implementation of \code{\link[wrMisc]{moderTest2grp}} for proteomics.
#' The final obkect also includes the results without moderation by \code{limma} (eg BH-FDR in \code{$nonMod.BH}). 
#' Furthermore, there is an option to make use of package ROTS (note, this will increase the time of computatins considerably).
#'  
#' @param dat (matrix or data.frame) main data (may contain NAs)
#' @param questNo (integer) specify here which question, ie comparison should be adressed
#' @param useCol (integer or character) 
#' @param grp (character or factor) 
#' @param annot (matrix or data.frame) 
#' @param ROTSn (integer) number of iterations ROTS runs (stabilization of reseults may be seen with >300) 
#' @param silent (logical) suppress messages
#' @param debug (logical) display additional messages for debugging
#' @param callFrom (character) allow easier tracking of message(s) produced
#' @return This function returns a limma-type S3 object of class 'MArrayLM' (which can be accessed like a list); multiple testing correction types or modified testing by ROTS may get included ('p.value','FDR','BY','lfdr' or 'ROTS.BH')
#' @seealso  \code{\link[wrMisc]{moderTest2grp}}, \code{\link[wrMisc]{pVal2lfdr}}, \code{\link[stats]{t.test}}, \code{ROTS} from the Bioconductor package \href{https://www.bioconductor.org/packages/release/bioc/html/ROTS.html}{ROTS}  
#' @examples
#' set.seed(2018);  datT8 <- matrix(round(rnorm(800)+3,1), nc=8, dimnames=list(paste(
#'   "li",1:100,sep=""), paste(rep(LETTERS[1:3],c(3,3,2)),letters[18:25],sep="")))
#' datT8[3:6,1:2] <- datT8[3:6,1:2] +3   # augment lines 3:6 (c-f) 
#' datT8[5:8,5:6] <- datT8[5:8,5:6] +3   # augment lines 5:8 (e-h) 
#' grp8 <- gl(3,3,labels=LETTERS[1:3],length=8)
#' datL <- list(data=datT8, filt= wrMisc::presenceFilt(datT8,grp=grp8,maxGrpM=1,ratMa=0.8))
#' testAvB0 <- wrMisc::moderTest2grp(datT8[,1:6], gl(2,3))
#' testAvB <- test2grp(datL, questNo=1)
#' @export
test2grp <- function(dat, questNo, useCol=NULL, grp=NULL, annot=NULL, ROTSn=0, silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## custom extracting data from list with $data and $filt
  ## return MA-type list with test resuls
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="test2grp")
  msg <- " 'dat' must be list containing $data and $filt, with same number of rows !!"
  if(!all(c("data","filt") %in% names(dat))) stop(msg)
  if(nrow(dat$filt) != nrow(dat$data)) stop(msg)
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE else debug <- FALSE
  
  questNa <- colnames(dat$filt)[questNo]
  questNa <- unlist(strsplit(questNa, "-"))
  if(is.null(useCol)) useCol <- lapply(questNa, grep, colnames(dat$data))
  if(is.null(grp)) {grp <- rep(questNa, sapply(useCol, length))
    if(!silent) message("Groups set automatic to ",wrMisc::pasteC(grp,quoteC="'")) }
  out <- wrMisc::moderTest2grp(dat$data[which(dat$filt[,questNo]), unlist(useCol)], gr=grp, limmaOutput=TRUE, addResults=c("lfdr","FDR","Mval","means","nonMod"), callFrom=fxNa)
  out$BH <- apply(out$p.value, 2, stats::p.adjust, method="BH")
  out$nonMod.BH <- stats::p.adjust(out$nonMod.p, method="BH")
  chLfdr <- try(find.package("fdrtools"), silent=TRUE)
  if(inherits(chLfdr, "try-error")) { 
      message(fxNa,"Package 'fdrtool' not found !  Please install for calculating lfdr-values ..") 
  } else out$nonMod.lfdr <- wrMisc::pVal2lfdr(out$nonMod.p)
  ## need to add : (non-moderated test and) ROTS
  if(length(ROTSn==1)) if(ROTSn >0 && !is.na(ROTSn)) {
    chPa <- requireNamespace("ROTS", quietly=TRUE)
    if(!chPa) { ROTSn <- NULL
      message(fxNa,"Package 'ROTS' not found ! Please install first .. setting  ROTSn=NULL") }
  } else ROTSn <- NULL      
  if(length(ROTSn)==1) if(ROTSn >0) {
    ## this part requires ROTS
    tmp <- ROTS::ROTS(dat$data[which(dat$filt[,questNo]), unlist(useCol)], groups=as.numeric(as.factor(grp)) -1, B=ROTSn)   # ,K=500  
    out$ROTS.p <- tmp$pvalue
    out$ROTS.BH <- stats::p.adjust(tmp$pvalue, method="BH") 
    if(! inherits(chLfdr, "try-error")) out$ROTS.lfdr <- wrMisc::pVal2lfdr(tmp$pvalue)
  }
  if(!is.null(annot)) out$annot <- as.matrix(annot[which(dat$filt[,questNo]),])
  out }
   
