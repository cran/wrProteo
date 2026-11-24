#' Write Sequences In Fasta Format To File
#'   
#' This function writes sequences from character vector or parsed matrix to file as fasta formatted file (close to \href{https://www.uniprot.org}{UniProt}) 
#' This function also allows comparing the main vector of sequences with a reference vector \code{ref} to check if any of the sequences therein are truncated.
#'  
#' @details 
#' If a named character vector is entered as \code{prot}, names will be used as annotatin of the sequence entries
#' Sequences without any names will be given generic headers like protein01 ... etc.
#' 
#' If \code{prot} is a matrix the column named 'sequence' or 'seq' will be considered as sequence and the rest will get concatenated 
#' to form theb fasta-header using the separator given by argument \code{sep}
#'  
#'  
#' @param prot (character) vector of sequenes, names will be used for fasta-header
#' @param fileNa (character) name (and path) for file to be written
#' @param sep (character) primary separator for fasta-header  (when \code{prot} is matrix)
#' @param suplSep (character) supplemental saparator for additional UniProt fields
#' @param suplColumns (character) column names to be used for supplemental saparator (for additional UniProt fields)
#' 
#' @param ref (character) optional/additional set of (reference-) sequences (only for comparison to \code{prot}), length of proteins from \code{prot} will be checked to mark truncated proteins by '_tru'
#' @param lineLength (integer, length=1) number of sequence characters per line (default 60, should be >1 and <10000)
#' @param eol (character) the character(s) to print at the end of each line (row); for example, eol = "\\r\\n" will produce Windows' line endings on a Unix-alike OS 
#' @param truSuf (character) suffix to be added for sequences found truncated when comparing with \code{ref} 
#' @param substName (character, length=1) in case no names/elements are given to be used as fasta-headers this term will be enumerated
#' @param silent (logical) suppress messages
#' @param debug (logical) supplemental messages for debugging
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This function writes the sequences from \code{prot} as fasta formatted-file
#' @seealso \code{\link{readFasta2}} for reading fasta, \code{write.fasta} from the package \href{https://CRAN.R-project.org/package=seqinr}{seqinr} 
#' @examples
#' prots <- c(SEQU1="ABCDEFGHIJKL", SEQU2="CDEFGHIJKLMNOP")
#' writeFasta2(prots, fileNa=file.path(tempdir(),"testWrite.fasta"), lineLength=6)
#' @export
writeFasta2 <- function(prot, fileNa=NULL, sep="|", suplSep=" ", suplColumns=c("ProteinName","OS=","OX=","GN=","PE=","SV="), ref=NULL, lineLength=60, eol="\n", truSuf="_tru", substName="protein", silent=FALSE, debug=FALSE, callFrom=NULL) {
  ##
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="writeFasta2")
  if(isTRUE(debug)) silent <- FALSE
  if(!isTRUE(silent)) silent <- FALSE
  namesInp <- c(deparse(substitute(prot)), deparse(substitute(fileNa)), deparse(substitute(ref)))

  if(length(prot) <1) stop(" '",namesInp,"' seems to be empty !  Expecting named character vector with (protein) sequence")
  chN <- length(names(prot)) <1                    # check names of 'prot'
  
  ## convert prot if matrix
  if(length(dim(prot))==2 && ncol(prot)==1) { 
    pNa <- names(prot) 
    prot <- as.character(prot)
    if(length(pNa) >0) names(prot) <- pNa else message(fxNa,"Input given as matrix but no column available for fasta-header") } 
  if(debug) {message(fxNa,"wFa0"); wFa0 <- list(prot=prot,fileNa=fileNa,sep=sep,ref=ref,lineLength=lineLength,eol=eol,truSuf=truSuf,substName=substName)}
    
  if(length(dim(prot))==2) {
    seqCol <- which(tolower(colnames(prot)) %in% c("sequence","seq"))
    if(length(seqCol) >0) {
      nonSeqCol <- (1:ncol(prot))[-seqCol]
      prot1 <- prot[,nonSeqCol, drop=FALSE]
      if(debug) {message(fxNa,"wFa1"); wFa1 <- list(prot=prot,fileNa=fileNa,sep=sep,suplSep=suplSep,suplColumns=suplColumns,ref=ref, seqCol=seqCol,nonSeqCol=nonSeqCol,prot1=prot1)}
      ## start combining columns with suplSep (if given & columns apply)
      if(length(suplSep)==1 && length(suplColumns) >0) {
        chSup <- colnames(prot1) %in% sub("=$","", suplColumns)
        if(any(chSup)) { combColS <- which(chSup)
          for(i in combColS) prot1[,which(!chSup)[sum(!chSup)]] <- paste(prot1[,which(!chSup)[sum(!chSup)]], prot1[,i], sep=suplSep)
          prot1[,which(!chSup)[sum(!chSup)]] <- gsub(paste0(wrMisc::protectSpecChar(suplSep),"NA"),"", prot1[,which(!chSup)[sum(!chSup)]])    # writing long series of ' NA NA'
          prot1 <- prot1[,which(!chSup)]          
        }
      }
      if(debug) {message(fxNa,"wFa2"); wFa2 <- list(prot=prot,fileNa=fileNa,sep=sep,suplSep=suplSep,suplColumns=suplColumns,ref=ref, seqCol=seqCol,nonSeqCol=nonSeqCol,prot1=prot1)}
      if(length(nonSeqCol) >1) {for(i in nonSeqCol[-1]) prot1 <- paste0(prot1, sep, prot[,i])
        prot1 <- gsub(paste0(wrMisc::protectSpecChar(sep),"NA"),"", prot1) }          # writing long series of '|NA|NA'      
    } 
    prot <- prot[,seqCol[1]]
    names(prot) <- prot1
    chN <- length(names(prot)) <1                    # update check names of 'prot'

  }  
  if(debug) {message(fxNa,"wFa3"); wFa3 <- list(prot=prot,fileNa=fileNa,sep=sep,suplSep=suplSep,suplColumns=suplColumns,ref=ref)}

  ## compare to 'ref' -if provided
  if(length(ref)==length(prot) && length(prot) >1) {
    if(length(names(ref)) <1) { if(!silent) message(fxNa," NOTE: '",namesInp[2],"' has no names ! Assuming same order as '",namesInp[1],"' !")
    } else {
      ## check order based on names
      heP <- sub(" [[:print:]]+", "", names(ref))       # get 1st word/identifyer per line
      heR <- sub(" [[:print:]]+", "", names(prot))
      if(!all(heP==heR)) {
        message(fxNa,"...adjusting order realtif to '",namesInp[2],"'")
        newOrd <- match(heR, heP)
        if(sum(is.na(newOrd)) >0 && !silent) message(fxNa," NOTE : ",sum(is.na(newOrd))," entries of 'ref' not found in 'prot'" )
        prot <- prot[newOrd] }
    }
    ## compare proteins : count number of characters (AAs)
    chLe <- sapply(prot, nchar)
    chLeRe <- sapply(ref, nchar)
    if(!silent) message(fxNa," (out of ",length(prot)," proteins) ",sum(chLe < chLeRe)," proteins shorter, ",sum(chLe > chLeRe)," proteins longer as in '",namesInp[2],"'")
    ## propagate names from ref, adj if truncated
    if(debug) {message(fxNa,"wFa4\n"); wFa4 <- list(prot=prot,ref=ref,heP=heP,heR=heR,chLe=chLe,chLeRe=chLeRe,chN=chN)}
    ## adjust names
    if(chN && length(names(ref))==length(ref)) names(prot) <- names(ref)    
    ## look for truncated proteins
    if(any(chLe < chLeRe)) { 
      redLe <- which(chLe < chLeRe)
      if(length(truSuf) <1 | any(is.na(truSuf))) {truSuf <- ""; if(!silent) message(fxNa," Suffix 'truSuf' set as empty")}
      ch1 <- grepl(" $",truSuf[1])
      if(!ch1) truSuf <- paste0(truSuf[1]," ")
      names(prot)[redLe] <- sub(" ", truSuf, names(prot)[redLe]) }      # append suffix to truncated protein sequences 
  } else { 
    if(length(ref) >0 && !silent) message(fxNa,"NOTE : '",namesInp[2],"' does NOT match '",namesInp[1],"', ignoring ...")
    if(chN) { if(!silent) message(fxNa," Note : '",namesInp[1],"' have NO NAMES, renaming to protein01 ... proteinN for fasta headers")
      names(prot) <- paste0(substName[1],sprintf(paste("%0",nchar(length(prot)),"d",sep=""),1:length(prot)))} }

  ## check names of lines/sequences for fasta-format  (heading '>')
  ch1 <- !grepl("^>", names(prot))
  if(any(ch1)) {if(all(ch1)) names(prot) <- paste0(">",names(prot)) else names(prot)[which(ch1)] <- paste0(">",names(prot)[which(ch1)])}
  if(debug) {message(fxNa,"wFa5"); wFa5 <- list(prot=prot, ch1=ch1,ref=ref,fileNa=fileNa,lineLength=lineLength)}

  ## prepare for adding header, split sequence in fixed length blocks
  out <- rep(">", length(prot)*2)
  out[2*(1:length(prot)) -1] <- names(prot)
  if(all(length(lineLength)==1, is.numeric(lineLength), !is.na(lineLength))) {
    if(lineLength <2 || lineLength >= 1e4) { lineLength <- 60
    if(!silent) message(fxNa,"Invalid entry for 'lineLength' (setting to default=60)")}
  }
  if(all(length(lineLength)==1, is.numeric(lineLength), !is.na(lineLength))) {
    out[2*(1:length(prot))] <- strsplit(gsub(paste0("([[:alpha:]]{",lineLength,"})"), "\\1 ", as.character(prot)), " ")  # first, use pattern recall to add spaces used then with strsplit
    out <- unlist(out, use.names=FALSE)}
  if(debug) {message(fxNa,"wFa6"); wFa6 <- list(prot=prot, ch1=ch1,ref=ref,fileNa=fileNa,lineLength=lineLength, out=out)}

  ## check fileName
  if(length(fileNa) <1) fileNa <- namesInp[2]
  if(!grepl("\\.fasta$",fileNa)) fileNa <- paste0(fileNa,".fasta")
  chFi <- file.exists(fileNa)
  if(!silent) {if(chFi) message(fxNa," NOTE : file '",fileNa,"' will be overwritten !") else if(debug) message(fxNa," Ready to write file ",fileNa)}

  ## setup file connection (see https://stackoverflow.com/questions/36933590/how-to-write-files-with-unix-end-of-lines-on-r-for-windows)
  con <- try(file(fileNa), silent=TRUE)
  
  if(inherits(con, "try-error")) warning(fxNa,"Failed to call 'file()' for '",fileNa,"' (check authorization/access)") else {
    if(debug) message(" file '",fileNa,"' : conection error : ",  "try-error" %in% class(con)) }

  ## open connection
  if(!isOpen(con=con, rw="wb")) { tmp <- try(open(con, open="wb"), silent=TRUE) } else tmp <- NULL   # 'wb' means : open for writing, binary mode
  if(debug) {message(fxNa,"wF7"); wF7 <- list(prot=prot, ch1=ch1,ref=ref,fileNa=fileNa,lineLength=lineLength,con=con,out=out)}

  if(inherits(tmp, "try-error")) warning("Failed to open connection to file '",fileNa,"' (check authorization/access)") else {
    if(debug) message(fxNa,"Open connection error ", inherits(tmp, "try-error"))
    on.exit(try(close(con), silent=TRUE), add=TRUE)
    ## write to file
    tmp <- try(writeLines(as.character(out), con=con, sep=eol), silent=TRUE)    #
    if(debug) message(fxNa," writeLines() error ",  "try-error" %in% class(tmp))
    ## close connection                
    close(con)
 }
}
   
