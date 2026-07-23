#' Read File Of Protein Sequences In Fasta Format
#'
#' Read fasta formatted file (from \href{https://www.uniprot.org}{UniProt}) to extract (protein) sequences and name.
#'  
#' @details
#' Read fasta-header -as is (ie without any parsing) : Set argument \code{tableOut=FALSE}.
#' If \code{tableOut=TRUE} the output will be organized as matrix for separating meta-annotation (eg uniqueIdentifier, entryName, proteinName, GN) in separate columns.
#' Please keep in mind that parsers wer primarily designed for the UniProt format.
#'
#' @param filename (character) names fasta-file to be read; .gz compressed files can be read, too (see examples)
#' @param delim (character) delimeter at header-line
#' @param databaseSign (character) characters at beginning right after the '>' (typically specifying the data-base-origin), they will be excluded from the sequance-header
#' @param removeEntries (character) if \code{removeEntries='empty'} allows removing entries without any sequence entries; 
#'   set to \code{removeEntries='duplicated'} to remove duplicate entries (same sequence and same header)  
#'   \code{removeEntries='allNA'} remove columns with all entries NA (if \code{tableOut=TRUE})
#' @param tableOut (logical) toggle to return named character-vector or matrix with enhaced parsing of fasta-header. 
#'   The resulting matrix will contain the comumns 'database','uniqueIdentifier','entryName','proteinName','sequence' and further columns depending on argument \code{UniprSep}
#' @param UniprSep (character) separators for further separating entry-fields if \code{tableOut=TRUE}, see also \href{https://www.uniprot.org/help/fasta-headers}{UniProt-FASTA-headers}
#' @param strictSpecPattern (logical or character) deprecated, this argument is not used any more
#' @param cleanCols (logical) deprecated, please use argument \code{removeEntries="allNA"} 
#' @param debug (logical) supplemental messages for debugging
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This function returns (depending on argument \code{tableOut}) a simple character vector (of sequences) with (entire) Uniprot annotation as name or 
#'   b) a matrix with columns: 'database','uniqueIdentifier','entryName','proteinName','sequence' and further columns depending on argument \code{UniprSep}
#' @seealso  \code{\link{writeFasta2}} for writing as fasta; for reading \code{\link[base]{readLines}} or  \code{read.fasta} from the package \href{https://CRAN.R-project.org/package=seqinr}{seqinr}
#' @examples
#' ## Tiny example with common contaminants
#' path1 <- system.file('extdata', package='wrProteo')
#' fiNa <-  "conta1.fasta.gz"
#' fasta1 <- readFasta2(file.path(path1, fiNa))
#' ## now let's read and further separate annotation-fields
#' fasta2 <- readFasta2(file.path(path1, fiNa), tableOut=TRUE)
#' str(fasta1)
#' @export
readFasta2 <- function(filename, delim="|", databaseSign=c("sp","tr","generic","conta","synt","gi"), removeEntries=NULL, tableOut=FALSE, UniprSep=c("OS=","OX=","GN=","PE=","SV="),
  strictSpecPattern=TRUE, cleanCols=TRUE, silent=FALSE, callFrom=NULL, debug=FALSE){
  ##  strictSpecPattern .. decide if species MUST be cited in 'entryName' as "_[[:upper]]+"  (eg 'THMS2_HUMAN')
  ## read fasta formatted file (from Uniprot) to extract (protein) sequences and name
  ## info about Uniprot fasta https://www.uniprot.org/help/fasta-headers
  ## return (based on 'tableOut') simple character vector (of sequence) with Uniprot ID as name or matrix with cols:'ID','Sequence','EntryName','ProteinName','OS','GN'
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="readFasta2")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) { silent <- FALSE } else { debug <- FALSE }
  ## initial test reading
  if(!file.exists(filename)) stop(" file ",filename," not existing")
  if(length(databaseSign) >0 && any(is.na(databaseSign))) databaseSign <- databaseSign[which(!is.na(databaseSign))]   # remove NAs
  out <- chLe <- dbSig <- byDBsig <- NULL

  ## start reading
  sca <- try(suppressWarnings(readLines(filename)), silent=TRUE)
  ## faster reading of file ? see https://www.r-bloggers.com/2011/08/faster-files-in-r/
  if(inherits(sca, "try-error")) stop(fxNa,"File ",filename," exits but could not read !")
  if(debug) {message(fxNa,"Successfully read file '",filename,"'")}
  ## abandon using scan due to cases of EOL during text read interfering with ...

  newLi <- grep("^>", sca)
  newLi <- if(is.list(newLi)) newLi <- sort(unlist(newLi))  else as.numeric(newLi)
  if(length(newLi) <1) stop(fxNa,"No instances of '>' found !  Maybe this is NOT a real fasta-file ?")

  if(debug) {message(fxNa,if(length(databaseSign) >0) c("Checking for database signs ",wrMisc::pasteC(databaseSign, quoteC="'"))," rf2"); 
    rf2 <- list(sca=sca,newLi=newLi,filename=filename,delim=delim,databaseSign=databaseSign,removeEntries=removeEntries,tableOut=tableOut,UniprSep=UniprSep) }
  
  id0 <- sca[newLi]
  ## combine sequence lines
  useLi <- cbind(newLi +1, c(newLi[-1] -1, length(sca)))
  ## note : if single line of sequence both values on same line have same index
  chLe <- useLi[,2] - useLi[,1] <0
  useLi <- cbind(useLi, empty=chLe)    # add 3rd column to properly consider/treat empty sequences (head present but not no sequence)
  if(any(chLe, na.rm=TRUE)) {
    if(any(c("empty","removeempty") %in% tolower(removeEntries), na.rm=TRUE) ) {
      if(!silent) message(fxNa,"Found ",sum(chLe)," case(s) of entries without any sequence underneith - omitting; bizzare !")
      useLi <- useLi[which(!chLe),]
      id0 <- id0[which(!chLe)]
    }
  }    
  seqs <- apply(useLi, 1, function(x) if(x[3]==0) paste(sca[x[1]:x[2]], collapse="") else "")
  if(debug) {message(fxNa,"rf3"); rf3 <- list(sca=sca,newLi=newLi,filename=filename,delim=delim,databaseSign=databaseSign,removeEntries=removeEntries,tableOut=tableOut,UniprSep=UniprSep,seqs=seqs,useLi=useLi,id0=id0) }
  
  if(isTRUE(tableOut)) {
    out <- .parseFastaHeader(header=sca[newLi], delim=delim, databaseSign=databaseSign, UniprSep=UniprSep, asList=TRUE, silent=silent, callFrom=fxNa, debug=debug)
    out <- cbind(out[[1]], sequence=seqs, out[[2]])
    if(any(c("allNA","emptyColumns") %in% removeEntries)) {
      chCol <- colSums(!is.na(out)) ==0
      if(all(chCol)) warning(fxNa,"All columns of file '",filename,"' seem empty !?!")
      if(any(chCol)) out <- out[, chCol, drop=FALSE]
    }
    if(debug) { message(fxNa,"finish with dim out ",nrow(out)," x ",ncol(out),"   rf10")}
  } else {         ## no mining/splitting of fasta-header
    out <- seqs
    names(out) <- sca[newLi]   #  '>' already removed
  }
  if(length(removeEntries) >0) {
    if(any(c("duplicated","dupl") %in% removeEntries)){
      chSe <- duplicated(seqs)
      if(any(chSe)) out <- if(length(dim(out))==2) out[which(!chSe), , drop=FALSE] else out[which(!chSe)] }
  }
  out }


#' Parse Fasta Header
#'
#' Parse fasta header (from \href{https://www.uniprot.org}{UniProt}) to extract different annotation fields
#'  
#' @param header (character) fasta-header
#' @param delim (character) delimeter (ie primary separator)
#' @param databaseSign (character) characters at beginning right after the '>' (typically specifying the data-base-origin), they will be excluded from the sequance-header
#' @param UniprSep (character) separators for further separating entry-fields if \code{tableOut=TRUE}; with these delimeter fields a space is assumed in addition to the separators;
#'   see also \href{https://www.uniprot.org/help/fasta-headers}{UniProt-FASTA-headers}
#' @param asList (logical) if \code{asList=TRUE},the function returns a list with two matrixes, one for primary parsing and 
#'   another matrix for further parsing (using \code{UniprSep}), otherwise all will be combined in single matrix
#' @param debug (logical) supplemental messages for debugging
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This function returns (depending on argument \code{asList}) a)  a matrix with columns: 'db','uniqueIdentifier','entryName','proteinName' and further columns depending on argument \code{UniprSep}
#'   of b) a list with matrix of primary parsing (argument \code{delim}) and matrix from further parsing (argument \code{UniprSep})
#' @seealso This function is use by  \code{\link{readFasta2}},  \code{\link{writeFasta2}} for writing as fasta; for reading \code{\link[base]{readLines}} or  \code{read.fasta} from the package \href{https://CRAN.R-project.org/package=seqinr}{seqinr}
#' @examples
#' .parseFastaHeader(">sp|P00760|TRY1_BOVIN Serine protease 1 OS=Bos taurus OX=9913 GN=PRSS1 PE=1")
#' @export
.parseFastaHeader <- function(header, delim="|", databaseSign=c("sp","tr","generic","conta","synt","gi"), UniprSep=c("OS=","OX=","GN=","PE=","SV="), asList=FALSE,
  silent=FALSE, callFrom=NULL, debug=FALSE) {
  ## mine fasta-header (see https://www.uniprot.org/help/fasta-headers)
  fxNa <- wrMisc::.composeCallName(callFrom, newNa=".parseFastaHeader")
  header <- sub("^>","", header)
  
  hasDBname <- FALSE
  mat2 <- out2 <- NULL
  useDel <- wrMisc::protectSpecChar(delim)
  ## split by delim
  spl <- strsplit(header, useDel)
  chLe <- sapply(spl, length)
  chLeUni <- sort(unique(chLe))
  if(length(chLeUni) >1) {            # multiple numbers of separators
    ## see if those that are shorter don't have databaseSign - if so assume as missing and replace by ''
    chDBsign <- sapply(spl, function(x) x[1] %in% databaseSign)
    chDBsum <- sum(chDBsign)
    if(chDBsum >0 && chDBsum < length(spl)) {
      hasDBname <- TRUE   
      spl[which(!chDBsign & chLe==chLeUni[1])] <- lapply(spl[which(!chDBsign & chLe==chLeUni[1])], function(x) c("",x))
    }
    ## update
    chLe <- sapply(spl, length)
    chLeUni <- sort(unique(chLe))
    if(length(chLeUni) >1) {  # make all of same number of splitted parts (add empty fields at end)
      if(!silent) message(fxNa,"Beware, variable number of separators in fasta-headers")
      spl <- lapply(spl, function(x) {if(length(x) < max(chLeUni)) c(x, rep("", max(chLeUni) -length(x))) else x})  
      if(debug) {message(fxNa,"rf4"); rf4 <- list(spl=spl,header=header,delim=delim,databaseSign=databaseSign,UniprSep=UniprSep) }
    }
  }    
  spl <- matrix(unlist(spl, use.names=FALSE), nrow=length(spl), byrow=TRUE)        # convert to matrix
  ## now see if 1st col contains databaseSign
  chDB <- spl[,1] %in% databaseSign
  if(any(chDB)) {
    hasDBname <- TRUE
    if(sum(!chDB) >0 && !silent) message(fxNa,"Note : ",sum(!chDB)," fasta-headers do not look like having data-base names (missing or not recognized ? check content of argument 'chDB' ?)") 
  }
  ## colnames
  defColNa <- c("database","uniqueIdentifier", "entryName","proteinName","other1")    # UniProt field names
  if(ncol(spl) > 4) defColNa <- c(defColNa, paste0("other",1+ 1:(ncol(spl) -4)))
  colnames(spl) <- defColNa[!hasDBname + 1:ncol(spl)]
  if(debug) {message(fxNa,"rf4b"); rf4b <- list(spl=spl,header=header, delim=delim,databaseSign=databaseSign,UniprSep=UniprSep) }
  
  ## Uniprot : split last part further
  ## check for 'NAME_SPECIES Name bla bla' in 2nd to 4th col
  if(ncol(spl) >2) {
    splCol <- 0
    EntryNameProteinNamePat <- "^ {0,1}[[:upper:]]+([[:digit:]]|[[:upper:]])+_[[:upper:]]+ [[:upper:]]+[[:lower:]]*"   # pattern for EntryName and ProteinName 
    chUniPr2 <- grepl(EntryNameProteinNamePat, spl[,2])   
    chUniPr3 <- grepl(EntryNameProteinNamePat, spl[,3])   
    chUniPr4 <- if(ncol(spl) >3) grepl(EntryNameProteinNamePat, spl[,4]) else rep(FALSE,nrow(spl))
    if(any(chUniPr2)) splCol <- 2 else {if(any(chUniPr3)) splCol <- 3 else {if(any(chUniPr4)) splCol <- 4 }}       
    if(debug) {message(fxNa,"rf4c"); rf4c <- list(spl=spl,header=header, delim=delim,databaseSign=databaseSign,UniprSep=UniprSep,splCol=splCol,chUniPr2=chUniPr2,chUniPr3=chUniPr3) }

    if(splCol >1) {
      colnames(spl) <- letters[1:ncol(spl)]                                 # default colnames
      useLiS <- which(get(paste0("chUniPr",splCol)) )
      splSec <- sub("^ {0,1}[[:upper:]]([[:upper:]]|[[:digit:]])*_[[:upper:]]+ ","", spl[useLiS, splCol])  # after ' ', ie keep ProteinName (wo ProteinName OS=..)
      splPri <- sub(" [[:upper:]].+","", sub("^ ","", spl[useLiS, splCol]))              # remove ProteinName OS=...  ie EntryName remains
      spl[useLiS, splCol] <- sub("^ ","", splPri)
      if(any(nchar(splSec) >0)) {
        spl <- cbind(spl, ProteinName=NA) 
        spl[useLiS, ncol(spl)] <- sub("^ ","", splSec) 
      }
      if(any(grepl("^[[:upper:]]+([[:digit:]]|[[:upper:]])+$", spl[,splCol-1]))) colnames(spl)[splCol-1] <- "UniqueIdentifier"
      if(any(grepl("^[[:upper:]]+([[:digit:]]|[[:upper:]])+_[[:upper:]]+", spl[,splCol]))) colnames(spl)[splCol] <- "EntryName"
      if(length(databaseSign >0) && any( grep(paste(paste0("(",databaseSign,")"), collapse="|"), spl[,1]) )) colnames(spl)[1] <- "db"        
      if(debug) {message(fxNa,"rf4d"); rf4d <- list(spl=spl,header=header, delim=delim,databaseSign=databaseSign,UniprSep=UniprSep,useLiS=useLiS,splPri=splPri) }
        
      ## further split (for "OS=","OX=" ...)
      if(length(UniprSep) >0) {
        UniprSep <- sub("^ ","",sub(" $","", UniprSep))                              # remove heading or tailing space
        out2 <- matrix(NA_character_, nrow=nrow(spl), ncol=length(UniprSep), dimnames=list(NULL, sub("=$","", UniprSep)))
        grUni <- lapply(UniprSep, grep, spl[,ncol(spl)])                   # which entires/lines concerned
        chUni <- which(sapply(grUni, length) >0)              # which separators concerned
        if(debug) {message(fxNa," rf6b"); rf6b <- list(spl=spl,out2=out2,grUni=grUni,chUni=chUni,header=header, delim=delim,databaseSign=databaseSign,UniprSep=UniprSep) }
       
        if(length(chUni) >0) {
          for(i in chUni) {
            tmp <- spl[, ncol(spl)] 
            aftS <- paste(sapply(UniprSep[-1*(1:i)], function(x) paste0("\ ",x,"[[:alnum:]]+[[:print:]]*")),collapse="|")   # regex pattern
            curS <- paste(sapply(UniprSep[i], function(x) paste0("^[[:print:]]* ",x)),collapse="|")                         # regex pattern
            out2[grUni[[i]],i] <- sub(aftS,"", sub(curS,"", tmp[grUni[[i]]]))           # extract ith part
            tmp[grUni[[i]]] <- sub(paste0(" ",UniprSep[i],".*"),"", spl[grUni[[i]], ncol(spl)])                  # all behind ith part
          }

          aftS <- paste0(" (",paste(sapply(UniprSep, function(x) paste0("(",x,")")),collapse="|"),").+")   # regex pattern          
          spl[,ncol(spl)] <- sub(aftS,"", spl[,ncol(spl)])
        }
        if(debug) {message(fxNa,"rf7"); rf7 <- list(spl=spl,out2=out2,grUni=grUni,chUni=chUni,header=header, delim=delim,databaseSign=databaseSign,UniprSep=UniprSep) }
        ## last attempt to complete specied based on species-part of EntryName
        if("EntryName" %in% colnames(spl) && length(out2) >0 && "OS" %in% colnames(out2)){
          noSpec <- is.na(out2[,"OS"]) & !is.na(spl[,"EntryName"])
          if(any(noSpec)) {
            secPart <- sub("^[[:alnum:]]+_","", spl[,"EntryName"])
            cheLi <- which(nchar(secPart) >1)
            if(length(cheLi) >0) { newSpec <- .commonSpecies()[match(paste0("_",secPart[cheLi]), .commonSpecies()[,1]), 2]
              if(any(!is.na(newSpec))) { out2[cheLi,"OS"] <- newSpec
                if(debug) message(fxNa,"Managed to recuperate ",sum(!is.na(newSpec))," (commonly known) species names (based in EntryName)") } 
            }
          }
        }
        if(debug) {message(fxNa,"rf8"); rf8 <- list(spl=spl,out2=out2,grUni=grUni,chUni=chUni,header=header, delim=delim,databaseSign=databaseSign,UniprSep=UniprSep) }
      
      }  
    }    
  } 
  if(debug) {message(fxNa," rf9"); rf9 <- list(spl=spl,header=header, delim=delim,databaseSign=databaseSign,UniprSep=UniprSep,mat2=mat2) }
  out <- if(isTRUE(asList)) list(spl, out2) else {if(length(out2) >0) cbind(spl, out2) else spl}
  out
}  
       
