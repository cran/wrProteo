#' Read file of protein sequences in fasta format
#'   
#' Read fasta formatted file (from \href{https://www.uniprot.org}{UniProt}) to extract (protein) sequences and name.
#' If \code{tableOut=TRUE} output may be organized as matrix for separating meta-annotation (eg GeneName, OrganismName, ProteinName) in separate columns.
#'  
#' @param filename (character) names fasta-file to be read
#' @param delim (character) delimeter at header-line
#' @param databaseSign (character) characters at beginning right afetr the '>' (typically specifying the data-base-origin), they will be excluded from the sequance-header
#' @param tableOut (logical) toggle to return named character-vector or matrix with enhaced parsing of fasta-header. The resulting matrix will contain the comumns 'database','uniqueIdentifier','entryName','proteinName','sequence' and further columns depending on argument \code{UniprSep} 
#' @param UniprSep (character) separators for further separating entry-fields if \code{tableOut=TRUE}, see also \href{https://www.uniprot.org/help/fasta-headers}{UniProt-FASTA-headers}  
#' @param cleanCols (logical) remove columns with all entries NA, if \code{tableOut=TRUE}
#' @param debug (logical) supplemental messages for debugging
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of message(s) produced
#' @return return (based on 'tableOut') simple character vector (of sequence) with Uniprot ID as name or matrix with columns: 'database','uniqueIdentifier','entryName','proteinName','sequence' and further columns depending on argument \code{UniprSep}
#' @seealso  \code{\link[base]{scan}} or \code{\link[seqinr]{read.fasta}}  
#' @examples
#' # tiny example with common contaminants 
#' path1 <- system.file('extdata',package='wrProteo')
#' fiNa <-  "conta1.fasta"
#' fasta1 <- readFasta2(file.path(path1,fiNa))
#' ## now let's read and further separate annotation-fields
#' fasta2 <- readFasta2(file.path(path1,fiNa),tableOut=TRUE)
#' str(fasta1)
#' @export
readFasta2 <- function(filename,delim="|",databaseSign=c("sp","tr","generic","gi"),tableOut=FALSE,UniprSep=c("OS=","OX=","GN=","PE=","SV="),cleanCols=TRUE,silent=FALSE,callFrom=NULL,debug=FALSE){
  ## read fasta formatted file (from Uniprot) to extract (protein) sequences and name
  ## info about Uniprot fasta https://www.uniprot.org/help/fasta-headers
  ## return (based on 'tableOut') simple character vector (of sequence) with Uniprot ID as name or matrix with cols:'ID','Sequence','EntryName','ProteinName','OS','GN'
  fxNa <- wrMisc::.composeCallName(callFrom,newNa="readFasta2")
  ## protect characters from delim
  deli2 <- if(nchar(delim==1)) paste("\\",delim,sep="") else paste(paste("\\",unlist(strsplit(delim,"")),sep=""),collapse="") #"
  ## initial test reading
  if(!file.exists(filename)) stop(" file ",filename," not existing")
  sca <- try(readLines(filename))
  if(any(class(sca) == "try-error")) stop(fxNa," file ",filename," exits but could not read !")
  ## abandon using scan due to cases of EOL during text read interfering with ...
  newLi <- grep("^>",sca)
  newLi <- if(is.list(newLi)) newLi <- sort(unlist(newLi))  else as.numeric(newLi)
  if(length(newLi) <1) stop(fxNa," no instances of 'databaseSign', ie '",paste(databaseSign,collapse=""),"' found")
  byDBsig <- sapply(databaseSign,function(x) grep(paste("^>",x,sep=""),sca[newLi]))
  names(byDBsig) <- databaseSign
  ## count occurance of prefix types
  chLe <- sapply(byDBsig,length)
  out <- NULL
  if(any(chLe >0)) {  # has prefix
    if(any(chLe <1)) byDBsig <- byDBsig[which(chLe >0)]   # eliminate prefix types not found
    dbSig <- rep(names(byDBsig)[1],length(newLi))
    if(length(byDBsig) >1) for(i in 2:length(byDBsig)) dbSig[byDBsig[[i]]] <- names(byDBsig)[i]
    id0 <- substr(sca[newLi],nchar(dbSig)+2+nchar(delim),nchar(sca[newLi]))     # head wo prefix
  } else {                                          # in case no separator found, use text available as ID and as name  
    id0 <- sca[newLi] ; dbSig <- NULL} 
  ## isolate ID : strsplit by delim
  sep1 <- strsplit(id0,delim,fixed=TRUE) 
  chLe <- sapply(sep1,length)
  id <- entryName <- sub("^ ","",sub(" $","",sapply(sep1,function(x) x[1])))    # remove heading or tailing space 
  if(any(chLe >1)) {           # use 1st as ID and last as names+further
    whCh <- which(chLe >0)
    entryName[whCh] <- sapply(sep1,function(x) x[length(x)]) } 
  if(any(chLe <2)) {  
    whCh <- which(chLe <1)
    id[whCh] <- entryName[whCh] <- id0[whCh] }                # no more separator, use text available as ID AND as name  
  entryName <- sub("^ ","",sub("\\.$","",entryName))          # remove heading space or tailing point
  useLi <- cbind(newLi+1,c(newLi[-1]-1,length(sca)))
  if(useLi[nrow(useLi),2] - useLi[nrow(useLi),1] ==0) useLi <- useLi[-nrow(useLi),]    # omit list if empty
  chLe <- useLi[,2] - useLi[,1] <0
  if(any(chLe)) {
    if(!silent) message(fxNa," found ",sum(chLe)," case(s) of headers without any sequence underneith; bizzare !")
    useLi <- useLi[which(!chLe),] }
  seqs <- apply(useLi,1,function(x) paste(sca[x[1]:x[2]],collapse=""))
  if(debug) message(fxNa," len seqs ",length(seqs))
  if(tableOut) {           # further separating name/description field 
    ## Uniprot headers : https://www.uniprot.org/help/fasta-headers
    ## >db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion
    UniprSep <- sub("^ ","",sub(" $","",UniprSep))                              # remove heading or tailing space 
    out <- matrix(NA,nrow=length(id0),ncol=length(UniprSep)+5,dimnames=list(NULL,c("database","uniqueIdentifier","entryName","proteinName","sequence",sub("=$","",sub("^ +","",UniprSep)) )))
    aftS <- paste(sapply(UniprSep,function(x) paste("\ [[:print:]]+",x,sep="")),collapse="|")
    aftS <- "[[:print:]]*[[:upper:]]+[[:digit:]]*_[[:upper:]]{2,} "                               # some CAPs + ev some digits
    ## suplID (if available)
    entryNameS <- sub(aftS,"",entryName)                                    # get everything before Uniprot like sparator  (space+2upper+"="+anyText)
    nch <- nchar(entryName)
    ncha <- nchar(entryNameS)
    suplID <- if(any(ncha >0)) substr(entryName,1,nch -1-ncha) else rep(NA,length(entryName))
    chS <- ncha == nchar(entryName)
    out[,c("database","uniqueIdentifier","entryName","proteinName","sequence")] <- cbind(dbSig,id,entryNameS,suplID,seqs)
    ## extract part after current uniprSep and the before something looking like next uniprSep
    grUni <- lapply(UniprSep,grep,entryNameS)                   # which entires/lines concerned
    chUni <- which(sapply(grUni,length) >0)              # which separators concerned
    UniprSep <- c(UniprSep,"ZYXWVUTSR=")                 # need to ad dummy sequence for last
    if(any(chUni)) for(i in chUni) {
      aftS <- paste(sapply(UniprSep[-1*(1:i)],function(x) paste("\ ",x,"[[:alnum:]]+[[:print:]]*",sep="")),collapse="|")
      curS <- paste(sapply(UniprSep[i],function(x) paste("^[[:print:]]* ",x,sep="")),collapse="|")
      out[grUni[[i]],c(3,5+i)] <- sub(aftS,"", sub(curS,"",entryNameS[grUni[[i]]]))
    }
    chNA <- colSums(!is.na(out)) <1
    if(any(chNA) & cleanCols) out <- if(sum(!chNA) >1) out[,which(!chNA)] else matrix(out[,which(!chNA)],ncol=1,dimnames=list(NULL,colnames(out)[which(!chNA)])) # remove columns with NA only
  } else {out <- seqs; names(out) <- entryName}  
  out }
       
