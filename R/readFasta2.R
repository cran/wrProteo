#' Read file of protein sequences in fasta format
#'   
#' Read fasta formatted file (from \href{https://www.uniprot.org}{UniProt}) to extract (protein) sequences and name.
#' If \code{tableOut=TRUE} output may be organized as matrix for separating meta-annotation (eg uniqueIdentifier, entryName, proteinName, GN) in separate columns.
#'  
#' @param filename (character) names fasta-file to be read
#' @param delim (character) delimeter at header-line
#' @param databaseSign (character) characters at beginning right after the '>' (typically specifying the data-base-origin), they will be excluded from the sequance-header
#' @param removeEntries (character) if \code{'empty'} allows removing entries without any sequence entries; set to \code{'duplicated'} to remove duplicate entries (same sequence and same header)
#' @param tableOut (logical) toggle to return named character-vector or matrix with enhaced parsing of fasta-header. The resulting matrix will contain the comumns 'database','uniqueIdentifier','entryName','proteinName','sequence' and further columns depending on argument \code{UniprSep} 
#' @param UniprSep (character) separators for further separating entry-fields if \code{tableOut=TRUE}, see also \href{https://www.uniprot.org/help/fasta-headers}{UniProt-FASTA-headers}  
#' @param cleanCols (logical) remove columns with all entries NA, if \code{tableOut=TRUE}
#' @param debug (logical) supplemental messages for debugging
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This function returns (depending on parameter \code{tableOut}) a) a simple character vector (of sequence) with Uniprot ID as name or b) a matrix with columns: 'database','uniqueIdentifier','entryName','proteinName','sequence' and further columns depending on argument \code{UniprSep}
#' @seealso  \code{\link{writeFasta2}} for writing as fasta, or for reading \code{\link[base]{scan}} or  \code{read.fasta} from the package \href{https://CRAN.R-project.org/package=seqinr}{seqinr} 
#' @examples
#' ## Tiny example with common contaminants 
#' path1 <- system.file('extdata',package='wrProteo')
#' fiNa <-  "conta1.fasta.gz"
#' fasta1 <- readFasta2(file.path(path1,fiNa))
#' ## now let's read and further separate annotation-fields
#' fasta2 <- readFasta2(file.path(path1,fiNa),tableOut=TRUE)
#' str(fasta1)
#' @export
readFasta2 <- function(filename, delim="|", databaseSign=c("sp","tr","generic","gi"), removeEntries=NULL, tableOut=FALSE, UniprSep=c("OS=","OX=","GN=","PE=","SV="),
  cleanCols=TRUE, silent=FALSE, callFrom=NULL, debug=FALSE){
  ## read fasta formatted file (from Uniprot) to extract (protein) sequences and name
  ## info about Uniprot fasta https://www.uniprot.org/help/fasta-headers
  ## return (based on 'tableOut') simple character vector (of sequence) with Uniprot ID as name or matrix with cols:'ID','Sequence','EntryName','ProteinName','OS','GN'
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readFasta2")
  if(isTRUE(debug)) silent <- FALSE
  if(!isTRUE(silent)) silent <- FALSE
  ## protect characters from delim
  deli2 <- if(nchar(delim==1)) paste0("\\",delim) else paste(paste0("\\",unlist(strsplit(delim,""))),collapse="") #"
  ## initial test reading
  if(!file.exists(filename)) stop(" file ",filename," not existing")
  sca <- try(readLines(filename))
  ## faster reading of file ? see https://www.r-bloggers.com/2011/08/faster-files-in-r/
  if(inherits(sca, "try-error")) stop(fxNa,"File ",filename," exits but could not read !")
  if(debug) {message(fxNa,"Successfully read file '",filename,"'")}
  ## abandon using scan due to cases of EOL during text read interfering with ...
  newLi <- grep("^>", sca)
  newLi <- if(is.list(newLi)) newLi <- sort(unlist(newLi))  else as.numeric(newLi)
  if(length(newLi) <1) stop(fxNa,"No instances of 'databaseSign', ie '",paste(databaseSign,collapse=""),"' found !  Maybe this is NOT a real fasta-file ?")
  byDBsig <- sapply(databaseSign, function(x) grep(paste0("^>",x), sca[newLi]))
  names(byDBsig) <- databaseSign
  if(debug) {message(fxNa,"Checking for database signs ",wrMisc::pasteC(databaseSign, quoteC="'"))}

  ## count occurance of prefix types
  chLe <- sapply(byDBsig,length)
  out <- NULL
  if(any(chLe >0, na.rm=TRUE)) {  # has prefix
    if(any(chLe <1, na.rm=TRUE)) byDBsig <- byDBsig[which(chLe >0)]   # eliminate prefix types not found
    dbSig <- rep(names(byDBsig)[1],length(newLi))
    if(length(byDBsig) >1) for(i in 2:length(byDBsig)) dbSig[byDBsig[[i]]] <- names(byDBsig)[i]
    id0 <- substr(sca[newLi], nchar(dbSig) +2 +nchar(delim), nchar(sca[newLi]))     # head wo prefix
  } else {                                          # in case no separator found, use text available as ID and as name  
    id0 <- sca[newLi] ; dbSig <- NULL} 
  if(debug) {message(fxNa," found ",wrMisc::pasteC(chLe)," occurances of database signs" )}

  ## check for empty sequences
  useLi <- cbind(newLi +1, c(newLi[-1] -1, length(sca)))
  ## note : if single line of sequence both values on same line have same index
  chLe <- useLi[,2] - useLi[,1] <0
  if(any(c("empty","removeempty") %in% tolower(removeEntries), na.rm=TRUE) & any(chLe, na.rm=TRUE)) {
    if(!silent) message(fxNa," found ",sum(chLe)," case(s) of entries without any sequence underneith - omitting; bizzare !")
    useLi <- useLi[which(!chLe),] 
    id0 <- id0[which(!chLe)]
    #sca <- sca[which(!chLe)]
  }
  if(debug) {message(fxNa," rf4")}

  ## isolate ID : strsplit by delim
  sep1 <- strsplit(id0, delim, fixed=TRUE) 
  chLe <- sapply(sep1, length)
  id <- sub("^ ","",sub(" $","",sapply(sep1, function(x) x[1])))    # remove heading or tailing space 
  chNa <- is.na(id)
  if(any(chNa, na.rm=TRUE)) {
   if(!silent) message(fxNa,"Note:  ",sum(chNa)," entries have no names, will be given names 'NONAME01'", if(sum(chNa)>1)" etc...")
   id[which(chNa)] <- paste0("NONAME",sprintf(paste0("%0",max(2,nchar(sum(chNa))),"d"), 1:sum(chNa)))
  }
  entryName <- id
  ## Try extracting 2nd part after ID
  if(any(chLe >1, na.rm=TRUE)) {           # use 1st as ID and last as names+further
    entryName <- unlist(sapply(sep1, function(x) if(length(x) <1) NA else x[min(2,length(x))])) }   # use 2nd after separator (or 1st if)
  chNa <- is.na(entryName)
  if(any(chNa, na.rm=TRUE)) entryName[which(chNa)] <- id[which(chNa)]      # no separator, use text available as ID and as name 

  if(debug) {message(fxNa,"Isolated ",length(id)," ids ( ",sum(chNa)," with same text as ID and sequence name)" )}
  entryName <- sub("^ ","", sub("\\.$","",entryName))          # remove heading space or tailing point
  seqs <- apply(useLi, 1, function(x) paste(sca[x[1]:x[2]], collapse=""))

  if(debug) {message(fxNa," rf5")}

  if(any(c("duplicate","duplicated") %in% tolower(removeEntries), na.rm=TRUE) ) {
    chDup <- duplicated(seqs, fromLast=FALSE) & duplicated(entryName, fromLast=FALSE)
    if(any(chDup, na.rm=TRUE)) {
      if(!silent) message(fxNa,"Removing ",sum(chDup)," duplicated entries (same sequence AND same header)")
      seqs <- seqs[which(!chDup)]
      entryName <- entryName[which(!chDup)]
  } }

  if(debug) message(fxNa," len seqs ",length(seqs))

  if(isTRUE(tableOut)) {           # further separating name/description field 
    ## Uniprot headers : https://www.uniprot.org/help/fasta-headers
    ## >db|uniqueIdentifier|entryName proteinName OS=OrganismName OX=OrganismIdentifier [GN=GeneName] PE=ProteinExistence SV=SequenceVersion
    UniprSep <- sub("^ ","",sub(" $","",UniprSep))                              # remove heading or tailing space 
    out <- matrix(NA, nrow=length(id0), ncol=length(UniprSep) +5, 
      dimnames=list(NULL,c("database","uniqueIdentifier","entryName","proteinName","sequence",sub("=$","",sub("^ +","",UniprSep)) )))
    aftS <- "[[:alpha:]][[:upper:]]*[[:digit:]]*[[:upper:]]*_[[:upper:]]+[[:digit:]]* "                               # some CAPs + ev some digits
    if(debug) {message(fxNa," rf6")}

    ## suplID (if available)
    entryNameS <- sub(aftS, "", entryName)                                    # get everything before Uniprot like sparator  (space+2upper+"="+anyText)
    nch <- nchar(entryName)
    ncha <- nchar(entryNameS)
    suplID <- if(any(ncha >0, na.rm=TRUE)) substr(entryName, 1, nch -1-ncha) else rep(NA, length(entryName))
    chS <- ncha == nchar(entryName)
    out[,c("database","uniqueIdentifier","entryName","proteinName","sequence")] <- cbind(dbSig, id, entryNameS, suplID, seqs)
    if(debug) {message(fxNa," rf7")}

    ## extract part after current uniprSep and the before something looking like next uniprSep
    grUni <- lapply(UniprSep, grep, entryNameS)                   # which entires/lines concerned
    chUni <- which(sapply(grUni, length) >0)              # which separators concerned
    UniprSep <- c(UniprSep, "ZYXWVUTSR=")                 # need to add dummy sequence for last
    if(any(chUni, na.rm=TRUE)) for(i in chUni) {
      aftS <- paste(sapply(UniprSep[-1*(1:i)], function(x) paste0("\ ",x,"[[:alnum:]]+[[:print:]]*")),collapse="|")
      curS <- paste(sapply(UniprSep[i], function(x) paste0("^[[:print:]]* ",x)),collapse="|")
      out[grUni[[i]], sub("=$","",UniprSep[i]) ] <- sub(aftS,"", sub(curS,"", entryNameS[grUni[[i]]]))           # c(3,5+i)
    }
    ## propagate NONAME to empty proteinName
    ch1 <- grepl("NONAME", out[,2]) & nchar(out[,"proteinName"]) <1
    if(any(ch1, na.rm=TRUE)) out[which(ch1),"proteinName"] <- out[which(ch1),2] 
    
    ## remove cols with all NA
    chNA <- colSums(!is.na(out)) <1
    if(any(chNA) && isTRUE(cleanCols)) out <- if(sum(!chNA) >1) out[,which(!chNA)] else matrix(out[,which(!chNA)], ncol=1, dimnames=list(NULL,colnames(out)[which(!chNA)])) # remove columns with NA only
  } else {out <- seqs; names(out) <- entryName}  
  out }
       
