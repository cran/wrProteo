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
    out <- .parseFastaHeader(header=sca[newLi], delim=delim, databaseSign=databaseSign, UniprSep=UniprSep, silent=silent, callFrom=fxNa, debug=debug)
    out <- cbind(out, sequence=seqs)
    if(any(c("allNA","empty","emptyColumns") %in% removeEntries)) {
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
#' @details
#' Expecting complete or partial of : >db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion
#' 
#' @param header (character) fasta-header
#' @param delim (character) delimeter (ie primary separator)
#' @param databaseSign (character) characters at beginning right after the '>' (typically specifying the data-base-origin), they will be excluded from the sequance-header
#' @param UniprSep (character) separators for further separating entry-fields if \code{tableOut=TRUE}; with these delimeter fields a space is assumed in addition to the separators;
#'   see also \href{https://www.uniprot.org/help/fasta-headers}{UniProt-FASTA-headers}
#' @param debug (logical) supplemental messages for debugging
#' @param silent (logical) suppress messages
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This function returns a matrix with columns: 'db','uniqueIdentifier','entryName','proteinName' and further columns depending on argument \code{UniprSep}
#' @seealso This function is use by  \code{\link{readFasta2}},  \code{\link{writeFasta2}} for writing as fasta; for reading \code{\link[base]{readLines}} or  \code{read.fasta} from the package \href{https://CRAN.R-project.org/package=seqinr}{seqinr}
#' @examples
#' .parseFastaHeader(">sp|P00760|TRY1_BOVIN Serine protease 1 OS=Bos taurus OX=9913 GN=PRSS1 PE=1")
#' @export
.parseFastaHeader <- function(header, delim="|", databaseSign=c("sp","tr","generic","conta","synt","gi"), UniprSep=c("OS=","OX=","GN=","PE=","SV="), 
  silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## mine fasta-header (see https://www.uniprot.org/help/fasta-headers)
  fxNa <- wrMisc::.composeCallName(callFrom, newNa=".parseFastaHeader")
  header <- sub("^>","", header)
  
  chSplit <- FALSE    # hasDBname <- 
  mat2 <- out2 <- splX <- NULL
  ## split by delim
                                                 
  con <- try(textConnection(header), silent=TRUE)  # only those with some XX=
  if(!inherits(con, "try-error")) {
    spl <- try(utils::read.table(con, sep = delim, fill = TRUE, as.is = TRUE, row.names = NULL, header=FALSE, check.names = FALSE), silent=TRUE)
    close(con) } else spl <- FALSE
  ## inspect 1st (as database-sign or UniqueIdentifier) & 2nd col (as UniqueIdentifier) 

  if(!isFALSE(spl) && !inherits(spl, "try-error")) {  
    ch1stCo <- cbind(db=if(length(databaseSign) ==0 || all(is.na(databaseSign))) grepl("^[[:lower:]]+$", spl[,1]) else spl[,1] %in% databaseSign , 
      ui=grepl("^[A-Z][A-Z0-9]+$", spl[,1]), shift=FALSE, seUI= if(ncol(spl) >1) grepl("^[A-Z][A-Z0-9]+$",spl[,2]) else FALSE)  # 2nd UniqueIdentifier
    if(ncol(spl) >1) ch1stCo[,3] <- ch1stCo[,2] & !ch1stCo[,1] & !ch1stCo[,4]        # test if 1st col is UniqueIdentifier (and not db sign) -> shift
    
    if(any(ch1stCo[,3])) {                                                            # mixed, ie some have db sign others don't -> shift row
      spl[which(ch1stCo[,3]),] <- cbind( rep("", sum(ch1stCo[,3])), spl[which(ch1stCo[,3]),-ncol(spl)])         # replace missing db-sign by ""
    }   
                                                                                                                  
    ## now see if 1st col contains databaseSign
    chDB <- if(length(databaseSign) !=0 && !all(is.na(databaseSign))) spl[,1] %in% databaseSign else grepl("[a-z]{2,9}", spl[,1])
    if(any(chDB)) {
      hasDBname <- TRUE
      if(sum(!chDB) >0 && !silent) message(fxNa,"Note : ",sum(!chDB)," fasta-headers do not look like having data-base names (missing or not recognized ? check content of argument 'chDB' ?)") 
    }
    ## colnames
    defColNa <- c("database","UniqueIdentifier", "EntryName","ProteinName","other1")    # UniProt field names
    if(ncol(spl) > 4) defColNa <- c(defColNa, paste0("other",1+ 1:(ncol(spl) -4)))
    colnames(spl) <- defColNa[!hasDBname + 1:ncol(spl)]
    if(debug) {message(fxNa,"rf4b"); rf4b <- list(spl=spl,header=header, delim=delim,databaseSign=databaseSign,UniprSep=UniprSep) }
    
    ## Uniprot : split last part further
    ## check for 'NAME_SPECIES Name bla bla' in 2nd to 4th col
    if(length(UniprSep) !=0 && !isFALSE(UniprSep) && !identical(UniprSep, NA) && ncol(spl) >2) {
      splCol <- 0
      #EntryNameProteinNamePat <- "^[A-Z0-9]+(_[A-Z][A-Z0-9_]*)? [[:upper:]][[:alnum:]]"   # pattern for EntryName AND ProteinName (open end)
      EntryNameProteinNamePat <- "^[A-Z0-9]+(_[A-Z][A-Z0-9_]*)? .+"   # pattern for EntryName AND something (open end) => still worth splitting to separate EntryName
      #GeneNamePat <- "^[[:alpha:]][[:alnum:]]+? .+"   # pattern for GN= AND something (open end) => still worth splitting to separate EntryName
      ## try to find out which column should be split (typically 3rd, sometimes 2nd or 4th)
#      chUniPr2 <- grepl(EntryNameProteinNamePat, spl[,2])   
#      chUniPr3 <- grepl(EntryNameProteinNamePat, spl[,3])   
#      chUniPr4 <- if(ncol(spl) >3) grepl(EntryNameProteinNamePat, spl[,4]) else rep(FALSE, nrow(spl))
      
      chUniPr <- cbind(c2=grepl(EntryNameProteinNamePat, spl[,2]) , c3=grepl(EntryNameProteinNamePat, spl[,3]) ,
        c4=if(ncol(spl) >3) grepl(EntryNameProteinNamePat, spl[,4]) else rep(FALSE, nrow(spl)) )
      chSplit <- colSums(chUniPr) >0   
      
      #chSplit <- c(c2=any(chUniPr2), c3=any(chUniPr3), c4=any(chUniPr4))            # which column looks like candidate for split (pick 1st possible)
      if(debug) {message(fxNa," rf5"); rf5 <- list(spl=spl,header=header, delim=delim,databaseSign=databaseSign,UniprSep=UniprSep,chSplit=chSplit) }
      if(any(chSplit)) {
        ## continue for searching other
        tmp <- as.character(spl[, 1 +which.max(chSplit)])   # this col should be further split (XX=xxx)     
        anySpl <- chUniPr[, which.max(chSplit)]

        hasXX <- if(length(UniprSep) ==0 || all(is.na(UniprSep))) rep(NA_character_, sum(anySpl)) else {
          if(isTRUE(UniprSep)) grepl("[A-Z]{2,5}=", tmp[which(anySpl)]) else grepl(paste0("(",paste(UniprSep, collapse=")|("),")" ), tmp[which(anySpl)]) }    # length adopted to those which need any split

        ## Need to split 'base-part' (EntryName ProteinName) in all cases (w & wo OS= ...)
        supl <- sub("^[A-Z0-9]+(_[A-Z][A-Z0-9_]*)? ","", tmp)
        #entrNa <- substr(tmp[which(!hasXX)], 1, nchar(tmp[which(!hasXX)]) - nchar(supl) -1 )
        splY <- cbind(EntryName=substr(tmp, 1, nchar(tmp) - nchar(supl) -1 ), ProteinName=supl)   # in some cases ProteinName need to get split further
          
        if(any(hasXX)) {
          ## split off XX=.. part
          #tmp <- as.character(spl[which(hasXX), 1+which.max(chSplit), drop=FALSE])
          tmp2 <- tmp[which(hasXX)]
          partBas <- sub(" [A-Z]{2}=.+", "", tmp2)    # all left part, still need to spearate EntryName AND ProteinName !
          #protNa <- sub("^[A-Z][A-Z_]+ ","", partBas)
          
          
          partXX <- substring(tmp2, nchar(partBas) +2 )  # only concerned XX= part  ==> SPLIT
            #cbind(partBas, partXX)
          maxXX <- max((nchar(partXX) - nchar(gsub("[A-Z]{2}=","", partXX)))/3)     # max number of XX= elements
          
          
      #    spl <- cbind(spl[which(hasXX), which.max(chSplit), drop=FALSE], protNa=sub("^[^ ]* ", "", partBas), entryName=sub(" .*", "", partBas),  matrix(NA_character_, nrow=sum(hasXX), ncol=maxXX))  # wo XX= entries
          if(debug) {message(fxNa,"rf6"); rf6 <- list(spl=spl,header=header, delim=delim,databaseSign=databaseSign,UniprSep=UniprSep,partBas=partBas,partXX=partXX,maxXX=maxXX,hasXX=hasXX, tmp=tmp,chSplit=chSplit) }           
    
          ## now need to fill : split XX= elements
          con <- textConnection(partXX)  # only those with some XX=
          splX <- try(utils::read.table(con, sep = "=", fill = TRUE, as.is = TRUE, row.names = NULL, header=FALSE, check.names = FALSE), silent=TRUE)
          close(con)
    
          ## need to extract last 2 caps of all entries (exept 1st col)
          tmp1 <- as.character(as.matrix(splX[,-1]))
          tmpEntr <- sub(" [A-Z]{2}$","", tmp1)              ##  left parts
          chNA <- is.na(tmpEntr)
          if(any(chNA)) tmpEntr[which(chNA)] <- ""
    #             matrix(tmpEntr, nr=nrow(splX))         ##  left parts
    #             matrix(substr(tmp1, nchar(tmpEntr)+2, nchar(tmp1)), nr=nrow(splX))  # right parts
          if(debug) {message(fxNa,"rf6b"); rf6b <- list(spl=spl,header=header, delim=delim,databaseSign=databaseSign,UniprSep=UniprSep,partBas=partBas,partXX=partXX,maxXX=maxXX,hasXX=hasXX, tmp=tmp,tmp1=tmp1,tmpEntr=tmpEntr,splX=splX,chSplit=chSplit) }           
    
          tmpTy <- matrix(substr(tmp1, nchar(tmpEntr) +2, nchar(tmp1)), nrow=nrow(splX))
          tmpTy <- cbind(splX[,1], tmpTy[,-ncol(tmpTy), drop=FALSE])
          typesXX <- unique(as.character(tmpTy))        
          typesXX <- typesXX[which(typesXX != "")][1:maxXX]         #  eg "OS" "GN" "PE" "SV"
          
          pos <- rep(1:nrow(tmpTy), length(typesXX)) + (match(tmpTy, typesXX) -1)* sum(hasXX)
          chNA <- is.na(pos) | tmpEntr %in% ""        
          if(any(chNA)) { pos <- pos[which(!chNA)]; tmpEntr <- tmpEntr[which(!chNA)]  }
          if(debug) {message(fxNa,"rf6c"); rf6c <- list(spl=spl,splX=splX,splY=splY, header=header, delim=delim,databaseSign=databaseSign,UniprSep=UniprSep,partBas=partBas,typesXX=typesXX,partXX=partXX,maxXX=maxXX,hasXX=hasXX,tmpTy=tmpTy,tmpEntr=tmpEntr,pos=pos, tmp=tmp,chSplit=chSplit) }           
    
          ## fill in to spl2
          spl2 <- matrix(NA_character_, nrow=sum(hasXX), ncol=length(typesXX), dimnames=list(NULL,typesXX))
          spl2 <- rep(NA_character_, prod(length(typesXX), sum(hasXX)))    
          spl2[pos] <- sub("^ ","", tmpEntr)
          spl2 <-  matrix(spl2, ncol=length(typesXX), dimnames=list(NULL, typesXX))        # OK
          splX <- cbind(EntryName=splY[which(hasXX),1], ProteinName=substring(partBas, nchar(splY[which(hasXX),1]) +2), spl2)
          colnames(splX) <- c("EntryName","ProteinName", typesXX )
        }
        if(debug) {message(fxNa,"rf7"); rf7 <- list(spl=spl,splY=splY,header=header, delim=delim,databaseSign=databaseSign,UniprSep=UniprSep,partBas=partBas,partXX=partXX,maxXX=maxXX,hasXX=hasXX,splY=splY,splX=splX,chDB=chDB, tmp=tmp,chSplit=chSplit) }           
        ## combine
        out <- cbind(spl[,1:+which.max(chSplit), drop=FALSE], splY)
        colnames(out) <- c(if(any(chDB)) "db","UniqueIdentifier","EntryName","ProteinName")[1:ncol(out)]
        dimN <- dim(out)
        if(any(hasXX)) {
          out <- cbind(out, matrix(NA_character_, nrow=nrow(spl), ncol=ncol(splX) -2, dimnames=list(out[,1], colnames(splX)[-(1:2)])) )
          out[which(hasXX), (dimN[2]):ncol(out) ] <- splX[,-1]
        }
        rownames(out) <- NULL    # duplicate rownames not allowed in matrix  
      }
    }     
    out
  }  
}        
 

