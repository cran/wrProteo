#' Read File Of Protein Sequences In Fasta Format
#'
#' Read fasta formatted file (from \href{https://www.uniprot.org}{UniProt}) to extract (protein) sequences and name.
#' If \code{tableOut=TRUE} output may be organized as matrix for separating meta-annotation (eg uniqueIdentifier, entryName, proteinName, GN) in separate columns.
#'
#' @param filename (character) names fasta-file to be read
#' @param delim (character) delimeter at header-line
#' @param databaseSign (character) characters at beginning right after the '>' (typically specifying the data-base-origin), they will be excluded from the sequance-header
#' @param removeEntries (character) if \code{'empty'} allows removing entries without any sequence entries; set to \code{'duplicated'} to remove duplicate entries (same sequence and same header)
#' @param tableOut (logical) toggle to return named character-vector or matrix with enhaced parsing of fasta-header. 
#'   The resulting matrix will contain the comumns 'database','uniqueIdentifier','entryName','proteinName','sequence' and further columns depending on argument \code{UniprSep}
#' @param UniprSep (character) separators for further separating entry-fields if \code{tableOut=TRUE}, see also \href{https://www.uniprot.org/help/fasta-headers}{UniProt-FASTA-headers}
#' @param strictSpecPattern (logical or character) pattern for recognizing EntryName which is typically preceeding ProteinName (separated by ' '); if \code{TRUE} the 
#'   name (capital letters and digits) must contain in the second part '_' plus capital letters, if \code{FALSE} the second part may be absent; if not matching pattern the text will be at the beggining of the ProteinName
#' @param cleanCols (logical) remove columns with all entries NA, if \code{tableOut=TRUE}
#' @param debug (logical) supplemental messages for debugging
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This function returns (depending on argument \code{tableOut}) a simple character vector (of sequences) with (entire) Uniprot annotation as name or 
#'   b) a matrix with columns: 'database','uniqueIdentifier','entryName','proteinName','sequence' and further columns depending on argument \code{UniprSep}
#' @seealso  \code{\link{writeFasta2}} for writing as fasta; for reading \code{\link[base]{scan}} or  \code{read.fasta} from the package \href{https://CRAN.R-project.org/package=seqinr}{seqinr}
#' @examples
#' ## Tiny example with common contaminants
#' path1 <- system.file('extdata', package='wrProteo')
#' fiNa <-  "conta1.fasta.gz"
#' fasta1 <- readFasta2(file.path(path1, fiNa))
#' ## now let's read and further separate annotation-fields
#' fasta2 <- readFasta2(file.path(path1, fiNa), tableOut=TRUE)
#' str(fasta1)
#' @export
readFasta2 <- function(filename, delim="|", databaseSign=c("sp","tr","generic","gi"), removeEntries=NULL, tableOut=FALSE, UniprSep=c("OS=","OX=","GN=","PE=","SV="),
  strictSpecPattern=TRUE, cleanCols=TRUE, silent=FALSE, callFrom=NULL, debug=FALSE){
  ##  strictSpecPattern .. decide if species MUST be cited in 'entryName' as "_[[:upper]]+"  (eg 'THMS2_HUMAN')
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
  sca <- try(suppressWarnings(readLines(filename)), silent=TRUE)
  ## faster reading of file ? see https://www.r-bloggers.com/2011/08/faster-files-in-r/
  if(inherits(sca, "try-error")) stop(fxNa,"File ",filename," exits but could not read !")
  if(debug) {message(fxNa,"Successfully read file '",filename,"'")}
  ## abandon using scan due to cases of EOL during text read interfering with ...
  newLi <- grep("^>", sca)
  newLi <- if(is.list(newLi)) newLi <- sort(unlist(newLi))  else as.numeric(newLi)
  if(length(newLi) <1) stop(fxNa,"No instances of 'databaseSign', ie '",paste(databaseSign,collapse=""),"' found !  Maybe this is NOT a real fasta-file ?")
  byDBsig <- sapply(databaseSign, function(x) grep(paste0("^>",x), sca[newLi]))
  names(byDBsig) <- databaseSign
  if(debug) {message(fxNa,"Checking for database signs ",wrMisc::pasteC(databaseSign, quoteC="'"),"   rf2"); rf2 <- list(sca=sca,byDBsig=byDBsig,newLi=newLi,filename=filename,delim=delim,databaseSign=databaseSign,removeEntries=removeEntries,tableOut=tableOut,UniprSep=UniprSep) }

  ## count occurance of prefix types
  chLe <- sapply(byDBsig, length)
  out <- NULL
  if(sum(chLe) < length(newLi)) {       # non-standard format or unknown databaseSign
    id0 <- sca[newLi]
    ch1 <- nchar(gsub("\\|","", id0)) - nchar(id0)
    badFo <- ch1 > min(ch1, na.rm=TRUE)
    if(any(badFo)) {    ## inconsistent format, suppose databaseSign is missing
      if(!silent) message(fxNa,"Found ",sum(badFo)," inconsistent entries with missing databaseSign, adding 'xx|'" )
      databaseSign <- union(databaseSign,"xx")
      sca[newLi[which(badFo)]] <- paste0(">xx|", sub("^>","", sca[newLi[which(badFo)]]))   # add unknwn  db sign
      spePat <- " - [[:upper:]][[:lower:]]+ [[:lower:]]+ \\([[:upper:]][[:lower:]]+\\)$"
      chSpe <- which(grepl(spePat, sca[newLi[which(badFo)]]) & !grepl(" OS=[[:upper:]][[:lower:]]", sca[newLi[which(badFo)]]))
      if(length(chSpe) >0) {         ## found potential species designation (eg ' - Homo sapiens (Human)' while 'OS=' is missing)
        if(!silent) message(fxNa,"Found ",length(chSpe)," potential species designations, adding as OS=" )
        woSpe <- nchar(sub(spePat,"", sca[newLi[which(badFo)]]))
        newSpe <- sub(" \\([[:upper:]][[:lower:]]+\\)$","", substr(sca[newLi[which(badFo)]], woSpe +4, nchar(sca[newLi[which(badFo)]])))
        for(i in unique(newSpe)) sca[newLi[which(badFo)]] <- sub(paste0(" - ",i),"", sca[newLi[which(badFo)]])
        ## future : try converting OS to OX - need to extend  .commonSpecies()
        sca[newLi[which(badFo)]] <- paste0(sca[newLi[which(badFo)]], " OS=",newSpe)
      }
    }
  }       # finish adding missing databaseSign and trying to recover species from non-standard format

  if(any(chLe >0, na.rm=TRUE)) {  # has prefix
    if(any(chLe <1, na.rm=TRUE)) byDBsig <- byDBsig[which(chLe >0)]   # eliminate prefix types not found
    dbSig <- sub("\\|.+","", sub("^>","", sca[newLi]))     # correct database column
    if(length(byDBsig) >1) for(i in 2:length(byDBsig)) dbSig[byDBsig[[i]]] <- names(byDBsig)[i]
    id0 <- substr(sca[newLi], nchar(dbSig) +2 +nchar(delim), nchar(sca[newLi]))     # head wo prefix
  } else {                                          # in case no separator found, use text available as ID and as name
    id0 <- sca[newLi] ; dbSig <- NULL}


  if(debug) {message(fxNa," found ",wrMisc::pasteC(chLe)," occurances of database signs   rf3"); rf3 <- list(sca=sca,out=out,newLi=newLi,chLe=chLe,byDBsig=byDBsig,newLi=newLi,filename=filename,delim=delim,databaseSign=databaseSign,removeEntries=removeEntries,tableOut=tableOut,UniprSep=UniprSep) }

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
  if(debug) {message(fxNa,"rf4"); rf4 <- list(sca=sca,out=out,chLe=chLe,newLi=newLi,id0=id0,filename=filename,delim=delim,databaseSign=databaseSign,removeEntries=removeEntries,tableOut=tableOut,UniprSep=UniprSep) }

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
  if(debug) {message(fxNa,"rf5"); rf5 <- list(out=out,chLe=chLe,newLi=newLi,entryName=entryName,id=id,filename=filename,delim=delim,databaseSign=databaseSign,removeEntries=removeEntries,tableOut=tableOut,UniprSep=UniprSep) }
  ## Try extracting 2nd part after ID
  if(any(chLe >1, na.rm=TRUE)) {           # use 1st as ID and last as names+further
    entryName <- unlist(sapply(sep1, function(x) if(length(x) <1) NA else x[min(2,length(x))])) }   # use 2nd after separator (or 1st if)
  chNa <- is.na(entryName)
  if(any(chNa, na.rm=TRUE)) entryName[which(chNa)] <- id[which(chNa)]      # no separator, use text available as ID and as name

  if(debug) {message(fxNa,"Isolated ",length(id)," ids ( ",sum(chNa)," with same text as ID and sequence name)" )}
  entryName <- sub("^ ","", sub("\\.$","",entryName))          # remove heading space or tailing point
  seqs <- apply(useLi, 1, function(x) paste(sca[x[1]:x[2]], collapse=""))

  if(debug) {message(fxNa," rf5b")}

  if(any(c("duplicate","duplicated") %in% tolower(removeEntries), na.rm=TRUE) ) {
    chDup <- duplicated(seqs, fromLast=FALSE) & duplicated(entryName, fromLast=FALSE)
    if(any(chDup, na.rm=TRUE)) {
      if(!silent) message(fxNa,"Removing ",sum(chDup)," duplicated entries (same sequence AND same header)")
      seqs <- seqs[which(!chDup)]
      entryName <- entryName[which(!chDup)]
  } }

  if(debug) {message(fxNa," len seqs ",length(seqs),"  rf6");
     rf6 <- list(out=out,chLe=chLe,newLi=newLi,entryName=entryName,id=id,filename=filename,delim=delim,databaseSign=databaseSign,removeEntries=removeEntries,tableOut=tableOut,UniprSep=UniprSep,strictSpecPattern=strictSpecPattern,seqs=seqs) }
  ## Uniil here 'entryName' is entire annot wo 1st block   


  if(isTRUE(tableOut)) {           # further separating name/description field
    ## Uniprot headers : https://www.uniprot.org/help/fasta-headers
    ## >db|uniqueIdentifier|entryName proteinName OS=OrganismName OX=OrganismIdentifier [GN=GeneName] PE=ProteinExistence SV=SequenceVersion
    UniprSep <- sub("^ ","",sub(" $","", UniprSep))                              # remove heading or tailing space
    out <- matrix(NA, nrow=length(id0), ncol=length(UniprSep) +5,
      dimnames=list(NULL,c("database","uniqueIdentifier","entryName","proteinName","sequence",sub("=$","",sub("^ +","",UniprSep)) )))
    if(debug) {message(fxNa," rf7"); rf7 <- list(out=out,chLe=chLe,newLi=newLi,entryName=entryName,id=id,filename=filename,delim=delim,databaseSign=databaseSign,removeEntries=removeEntries, strictSpecPattern=strictSpecPattern,tableOut=tableOut,UniprSep=UniprSep,dbSig=dbSig,seqs=seqs)}

    ## suplID (if available)
    ## check if 1st word of entryName looks like ID
    entryNameS <- sub(" [[:upper:]].+","", entryName)     # first block before space+Upper  => should become 'entryName'
    ## define pattern how 'entryName' should look like
    paEN <- if(is.character(strictSpecPattern) && length(strictSpecPattern) ==1) strictSpecPattern else {
      if(isTRUE(strictSpecPattern)) "^([[:upper:]]+[[:digit:]]*)([[:upper:]]+[[:digit:]]*){0,1}(_[[:upper:]]+[[:digit:]]*){0,3}$" else 
                                    "^[[:upper:]]+[[:digit:]]([[:upper:]]|[[:digit:]])*(_[[:upper:]]+){0,1}"}
    chEn <- grepl(paEN, entryNameS)
    
    ## 20sep24 problem : nothing remains !
    
    if(any(!chEn)) entryNameS[which(!chEn)] <- ""
    if(any(chEn)) entryName[which(chEn)] <- substring(entryName[which(chEn)], nchar(entryNameS[which(chEn)]) +2)     # shorten entryName since UniqueIdentifier was split off
    out[,c("database","uniqueIdentifier","entryName","proteinName","sequence")] <- cbind(dbSig, id, entryNameS, entryName, seqs)
    if(debug) {message(fxNa," rf8"); rf8 <- list(out=out,chLe=chLe,newLi=newLi,entryName=entryName,entryNameS=entryNameS,id=id,filename=filename,delim=delim,databaseSign=databaseSign,removeEntries=removeEntries,tableOut=tableOut,UniprSep=UniprSep,seqs=seqs)}

    ## extract part after current uniprSep and the before something looking like next uniprSep
    grUni <- lapply(UniprSep, grep, entryName)                   # which entires/lines concerned
    chUni <- which(sapply(grUni, length) >0)              # which separators concerned
    UniprSep <- c(UniprSep, "ZYXWVUTSR=")                 # need to add dummy sequence for last
    if(any(chUni, na.rm=TRUE)) for(i in chUni) {
      aftS <- paste(sapply(UniprSep[-1*(1:i)], function(x) paste0("\ ",x,"[[:alnum:]]+[[:print:]]*")),collapse="|")
      curS <- paste(sapply(UniprSep[i], function(x) paste0("^[[:print:]]* ",x)),collapse="|")
      out[grUni[[i]], sub("=$","",UniprSep[i]) ] <- sub(aftS,"", sub(curS,"", entryName[grUni[[i]]]))           # c(3,5+i)
      out[grUni[[i]], "proteinName"] <- sub(paste0(" ",UniprSep[i],".*"),"", out[grUni[[i]], "proteinName"])                   # trim current & all behind
    }
    if(debug) {message(fxNa," rf9"); rf9 <- list(out=out,chLe=chLe,newLi=newLi,entryName=entryName,entryNameS=entryNameS,id=id,filename=filename,delim=delim,databaseSign=databaseSign,removeEntries=removeEntries,tableOut=tableOut,UniprSep=UniprSep)}
    ## propagate NONAME to empty proteinName
    ch1 <- grepl("NONAME", out[,2]) & nchar(out[,"proteinName"]) <1
    if(any(ch1, na.rm=TRUE)) out[which(ch1),"proteinName"] <- out[which(ch1),2]
  
    ## remove cols with all NA
    chNA <- colSums(!is.na(out)) <1
    if(any(chNA) && isTRUE(cleanCols)) out <- if(sum(!chNA) >1) out[,which(!chNA)] else matrix(out[,which(!chNA)], ncol=1, dimnames=list(NULL,colnames(out)[which(!chNA)])) # remove columns with NA only
  } else {out <- seqs; names(out) <- paste(id, entryName)}
  if(debug) {message(fxNa,"finish with dim out ",nrow(out)," x ",ncol(out),"   rf9b")}
  out }
  
