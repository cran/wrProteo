#' Compare in-silico digested proteomes for unique and shared peptides, counts per protein or as peptides
#'   
#' Compare in-silico digested proteomes for unique and shared peptides, counts per protein or as peptides.
#' The in-silico digestion may be performed separately using the package \href{https://bioconductor.org/packages/release/bioc/html/cleaver.html}{cleaver}.
#' Note: input must be list (or multiple names lists) of proteins with their respective peptides (eg by in-silico digestion).
#'  
#' @param ... (list) multiple lists of (ini-silico) digested proteins (typically protein ID as names) with their respectice peptides (AA sequence), one entry for each species
#' @param prefix (character) optional (species-) prefix for entries in '...', will be only considered if '...' has no names
#' @param sep (character) concatenation symbol
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of message(s) produced
#' @return This function returns a list with $byPep as list of logical matrixes for each peptide (as line) and unique/shared/etc for each species; $byProt as list of matrixes with count data per proten (as line) for each species; $tab with simple summary-type count data   
#' @seealso  \code{\link{readFasta2}} and/or \code{cleave-methods} in package \href{https://bioconductor.org/packages/release/bioc/html/cleaver.html}{cleaver}    
#' @examples
#' ## The example mimics a proteomics experiment where extracts form E coli and 
#' ## Saccharomyces cerevisiae were mixed, thus not all peptdes may occur unique.  
#' (mi2 = countNoOfCommonPeptides(Ec=list(E1=letters[1:4],E2=letters[c(3:7)],
#'   E3=letters[c(4,8,13)],E4=letters[9]),Sc=list(S1=letters[c(2:3,6)], 
#'   S2=letters[10:13],S3=letters[c(5,6,11)],S4=letters[c(11)],S5="n")))
#' ##  a .. uni E, b .. inteR, c .. inteR(+intra E), d .. intra E  (no4), e .. inteR, 
#' ##  f .. inteR +intra E   (no6), g .. uni E, h .. uni E  no 8), i .. uni E, 
#' ##  j .. uni S (no10), k .. intra S  (no11), l .. uni S (no12), m .. inteR  (no13)
#' lapply(mi2$byProt,head)
#' mi2$tab
#' @export
countNoOfCommonPeptides <- function(...,prefix=c("Hs","Sc","Ec"),sep="_",silent=FALSE,callFrom=NULL) {
  ## compare in-silico digested proteomes for unique and shared peptides, counts per protein
  ## .. input must be lists of proteins wther their respective peptides
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="countNoOfCommonPeptides")
  inp <- list(...)
  chInp <- c("prefix","sep","silent","callFrom") %in% names(inp)
  if(any(chInp)) inp <- inp[which(!chInp)]
  if(length(inp) <2) stop("not sufficient input elements - nothing to do !")
  chSep <- sapply(inp,function(x) length(grep(sep,names(x))) >0)
  if(any(chSep)) message("Trouble ahead : Separator 'sep' also appears in sequence  names !!")
  ## main
  chN <- names(inp)
  if(length(names(inp)) >0) if(all(nchar(names(inp)) >0)) prefix <- names(inp)     # use preferetially names of proteins as in main input
  nBySet <- lapply(inp,function(x) sapply(x,length))
  nBySet2 <- sapply(nBySet,sum)
  seqs <- unlist(inp,recursive=TRUE,use.names=FALSE)
  .countFra <- function(x) paste(1:x,sep,rep(x,x),sep="")  
  .firstOfRep <- function(x) duplicated(x,fromLast=TRUE) & !duplicated(x,fromLast=FALSE)          # find first of replicated (mark as T)
  names(seqs) <- paste(rep(prefix,nBySet2),sep,unlist(lapply(unlist(nBySet,use.names=FALSE),.countFra)),sep, rep(unlist(lapply(inp,names)),unlist(nBySet)),sep="")
  seqAnn <- cbind(species=rep(prefix,nBySet2), pepNo=unlist(lapply(unlist(nBySet),function(x) cbind(1:x))),
    pepTot=unlist(lapply(unlist(nBySet),function(x) cbind(rep(x,x)))), protID=rep(unlist(lapply(inp,names)),unlist(nBySet)))
  ## make matrix with prot names & number if pep ? (no need for strsplit later - alternative to concatenate as names in seqs)
  staSto <- cumsum(nBySet2)
  staSto <- cbind(sta=c(1,staSto[-length(staSto)]+1),stop=staSto)
  out <- list(byPep=list(),byProt=list(),tab=list())
  matColNa <- paste(rep(c("uni","shared","red"),2), rep(c("IntrA","InteR"),each=3),sep="")
  matColNa <- c("uni0","shared0","red0",  "unique","sharedInteR","sharedIntrA","redInteR","redIntraA")
  for(i in 1:length(inp)) {
    mat <- matrix(FALSE,nrow=nBySet2[i],ncol=length(matColNa),dimnames=list(unlist(inp[[i]],use.names=FALSE),matColNa))
    frL <- duplicated(seqs[staSto[i,1]:staSto[i,2]],fromLast=TRUE)
    frB <- duplicated(seqs[staSto[i,1]:staSto[i,2]],fromLast=FALSE)
    mat[,1:3] <- matrix(c(!frL & !frB, frL & !frB, frB),ncol=3)
    ## now compose unique, shared (as 1st occurance of non-unique) and redundant (further repeated sequences) etc
    com <- unique(seqs[staSto[i,1]:staSto[i,2]])
    chSh <- com %in% unique(seqs[-1*(staSto[i,1]:staSto[i,2])])                   # check for common with any other spec
    if(any(chSh)) {
      com <- seqs[staSto[i,1]:staSto[i,2]] %in% com[which(chSh)]     # index of common (among cur spec)
      uniqCom <- unique(seqs[staSto[i,1]:staSto[i,2]][which(com)])
      fir <- .firstOfRep(c(seqs[staSto[i,1]:staSto[i,2]][which(com)],uniqCom))[1:sum(com)]    # find 1st occur in subgroup of common (needed to re-inject unique common to make sure 1st always shows up)
      mat[which(mat[,1] & !com),4] <- TRUE           # unique : unique intra and NOT in common
      mat[which(com)[which(fir)],5] <- TRUE          # shared inter : first inst of common
      mat[which(com)[which(!fir)],7] <- TRUE         # redun inter : not first inst of common      
      fir <- .firstOfRep(seqs[staSto[i,1]:staSto[i,2]][which(!com)])    # find 1st occur in subgroup of non-common
      mat[which(!com)[which(fir)],6] <- TRUE         # shared intra : first inst of non-common (+ later remove unique)
      mat[which(!com)[which(!fir)],8] <- TRUE        # redun intra : not first inst of non-common (+later remove unique)
      mat[which(mat[,4]),7:8] <- FALSE               # unique can't be redundant (or shared)
    }
    out$byPep[[i]] <- mat
    names(out$byPep)[i] <- prefix[i] 
    ## exploit by prot
    tmp <- matrix(unlist(by(mat, rep(names(inp[[i]]), nBySet[[i]]), function(x) colSums(as.matrix(x),na.rm=TRUE))),
      ncol=ncol(mat), byrow=TRUE, dimnames=list(names(inp[[i]]),colnames(mat)))
    out[["byProt"]][[i]] <- cbind(nPep=nBySet[[i]], matrix(unlist(by(mat,rep(names(inp[[i]]), nBySet[[i]]), function(x) colSums(as.matrix(x),na.rm=TRUE))),
      ncol=ncol(mat), byrow=TRUE, dimnames=list(names(inp[[i]]),colnames(mat))) )                            # cut in list of matrixes
    names(out[["byProt"]])[[i]] <- prefix[i] 
    ## summarize
    supl <- c(
      n1pep=sum(nBySet[[i]]==1), n2pep=sum(nBySet[[i]] ==2),  n3pep=sum(nBySet[[i]]==3), n4fpep=sum(nBySet[[i]] >3),               # n at 1 pep, 2 pep, 3 pep, n at >3 pep
      n1pepSpeIntra=sum(nBySet[[i]]==1 & out[["byProt"]][[i]][,"uni0"]==1),                    # n at 1 pep which is unique.intra
      n1pepSpeInter=sum(nBySet[[i]]==1 & out[["byProt"]][[i]][,"unique"]==1),                  # n at 1 pep which is unique, #1 spec pep,
      n2pep2SpeIntra=sum(nBySet[[i]] ==2 & out[["byProt"]][[i]][,"uni0"] ==2),                 # n = 2pep with  2 spec intra
      n2pep2SpeInter=sum(nBySet[[i]] ==2 & out[["byProt"]][[i]][,"unique"] ==2),               # n = 2pep with  2 spec iner
      inf2SpePepIntra=sum(out[["byProt"]][[i]][,"uni0"] <2), inf2SpePep=sum(out[["byProt"]][[i]][,"unique"] <2),                   # <2 spec pep at intra-species, <2 spec pep
      min2pep1Intra=sum(nBySet[[i]] >1 & out[["byProt"]][[i]][,"uni0"] >0),                                      # min 2 pep & min 1 specif
      min2pep1Inter=sum(nBySet[[i]] >1 & out[["byProt"]][[i]][,"unique"] >0),                                      # min 2 pep & min 1 specif
      min2speIntra=sum(out[["byProt"]][[i]][,"uni0"] >1),                                         # (min 2 pep &) min 2 specif intra
      min2speIner=sum(out[["byProt"]][[i]][,"unique"] >1),                                        # (min 2 pep &) min 2 specif inter
      nLost2pepSpec=sum(out[["byProt"]][[i]][,"uni0"] >1 & out[["byProt"]][[i]][,"unique"] <2) )   # min 2 pep as single species but when mult species combined less than 2 pep
    out$tab <- if(length(out$tab) <1) as.matrix(c(nProt=length(inp[[i]]),nTotPep=nBySet2[[i]],nPepSpec=sum(mat[,"uni0"])+sum(mat[,"shared0"]), 
      colSums(mat),supl)) else cbind(out$tab,c(nProt=length(inp[[i]]),nTotPep=nBySet2[[i]],nPepSpec=sum(mat[,"uni0"])+sum(mat[,"shared0"]),colSums(mat),supl))
    colnames(out$tab)[ncol(out$tab)] <- prefix[i]    
    }
  out } 
   
