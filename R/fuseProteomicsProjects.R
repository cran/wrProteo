#' Combine Multiple Proteomics Data-Sets
#'
#' This function allows combining up to 3 separate data-sets previously imported using wrProteo.
#' 
#' @details 
#' Some quantification software way give some identifyers multiple times, ie as multiple lines (eg for different modifictions or charge states, etc).
#' In this case this function tries first to summarize all lines with identical identifyers (using the function \code{\link[wrMisc]{combineRedundLinesInList}}
#' which used by default the median value). 
#' Thus, it is very important to know your data and to understand when lines that appear with the same identifyers should/may be fused/summarized without 
#' doing damage to the later biological interpretation ! The user may specify for each dataset the colum out of the protein/peptide-annotation to use
#' via the argument \code{columnNa}. 
#' Then, this content will be matched as identical match, so when combining data from different software special care shoud be taken !
#' 
#' Please note, that (at this point) the data from different series/objects will be joined as they are, ie without any additional normalization.
#' It is up to the user to inspect the resulting data and to decide if and which type of normalization may be suitable !
#' 
#' Please do NOT try combining protein and peptide quntification data.
#' 
#' @param x (list) First Proteomics data-set
#' @param y (list) Second Proteomics data-set
#' @param z (list) optional third Proteomics data-set
#' @param columnNa (character) column names from annotation 
#' @param NA.rm (logical) remove \code{NA}s
#' @param listNa (character) names of key list-elemnts from \code{x} to be treated; the first one is used as pattern for the format of quantitation data, 
#' ,  the last one for the annotation data
#' @param all (logical) union of intersect or merge should be performed between x, y and z
#' @param textModif (character) Additional modifications to the identifiers from argument \code{columnNa}; 
#'    so far intregrated: \code{rmPrecAA} for removing preceeding caps letters (amino-acids, eg [KR].AGVIFPVGR.[ML] => AGVIFPVGR) 
#'    or \code{rmTerminalDigit} for removing terminal digits (charge-states)
#' @param shortNa (character) for appending to output-colnames
#' @param retProtLst (logical) return list-object similar to input, otherwise a matrix of fused/aligned quantitation data
#' @param silent (logical) suppress messages
#' @param debug (logical) additional messages for debugging
#' @param callFrom (character) allow easier tracking of messages produced
#' @return This function returns a list with the same number of list-elements as  \code{$x}, ie typically this contains :
#'   \code{$raw} (initial/raw abundance values), \code{$quant} with final normalized quantitations, 
#'   \code{$annot}, optionally \code{$counts} an array with number of peptides, \code{$quantNotes} or \code{$notes}
#' @seealso \code{\link[stats]{sd}}
#' @examples
#' path1 <- system.file("extdata", package="wrProteo")
#' dataMQ <- readMaxQuantFile(path1, specPref=NULL, normalizeMeth="median")
#' MCproFi1 <- "tinyMC.RData"
#' dataMC <- readMassChroQFile(path1, file=MCproFi1, plotGraph=FALSE)
#' dataFused <- fuseProteomicsProjects(dataMQ, dataMC)
#' dim(dataMQ$quant)
#' dim(dataMC$quant)
#' dim(dataFused$quant)
#' @export
fuseProteomicsProjects <- function(x, y, z=NULL, columnNa="Accession", NA.rm=TRUE, listNa=c(quant="quant",annot="annot"), all=FALSE, textModif=NULL, shortNa=NULL, retProtLst=FALSE, silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## 
  ##
  ## 'listNa' names of list-elements containing quantitation data (1st position) and protein/line annotation (2nd position)
  ## 'columnNa' column names from annotation (equiv to argument 'by' in merge() )
  ## 'shortNa' for appending to output-colnames
  ## 'retProtLst' return list-object similar to input, otherwise a matrix of fused/aligned quantitation data
  ## set 'textModif' (character) to 'rmTerminalDigit' for treating peptides from DiaNN
  ##  note : no normalization by default
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="fuseProteomicsProjects")
  #namesXY <- c(deparse(substitute(x)), deparse(substitute(y)), if(length(z) >0) deparse(substitute(z)))
  namesXY <- c(deparse(substitute(x)), deparse(substitute(y)), deparse(substitute(z)))
  out <- NULL
  if(length(x) <1) stop("Argument 'x' seems empty !")
  if(length(y) <1) stop("Argument 'y' seems empty !")
  if(debug) message(fxNa,"  fPL0")
  ## check argument 'shortNa'
  if(length(shortNa) <1) shortNa <- wrMisc::trimRedundText(namesXY)
  chNA <- is.na(shortNa)
  if(any(chNA)) shortNa[which(chNA)] <- ""
  ## check argument 'columnNa'
  if(length(columnNa) <1) stop("Missing 'columnNa' (should designate colnames of x$listNa[2], etc)")
  if(length(columnNa)==1) columnNa <- rep(columnNa, 2 +(length(z) >0))
  
  datOK <- columnNa[1] %in% colnames(x[[listNa[length(listNa)]]])
  if(!datOK) stop("missing column '",columnNa[1],"' in ",namesXY[1])
  ch2 <- columnNa[2] %in% colnames(y[[listNa[length(listNa)]]])
  if(!ch2) stop("missing column '",columnNa[2],"' in ",namesXY[2])
  if(length(z) >0) {
    ch1 <- columnNa[3] %in% colnames(z[[listNa[length(listNa)]]])
    if(!ch1) stop("missing column '",columnNa[3],"' in ",namesXY[3])
  }
      
  annN <- list(x=x[[listNa[length(listNa)]]][,which(colnames(x[[listNa[length(listNa)]]])==columnNa[1])[1]],
    y=y[[listNa[length(listNa)]]][,which(colnames(y[[listNa[length(listNa)]]])==columnNa[2])[1]], 
    z=if(length(z) >0) z[[listNa[length(listNa)]]][,which(colnames(z[[listNa[length(listNa)]]])==columnNa[3])[1]] else NULL)
    
  ## check IDs (all identical or all NA)
  badID <- lapply(annN, function(w) if(length(w) >0) {all(is.na(w)) || sum(duplicated(w)) == length(w) -1} else NULL)
  if(debug) {message(fxNa," fPL0b ++"); fPL0b <- list()} 
  if(any(unlist(badID[1:2]))) stop(fxNa," The annotation given for ",wrMisc::pasteC(namesXY[which(unlist(badID))], quoteC="'")," with column '",columnNa,"' is all NA or all redundant")
                                                                                                             
  if("rmPrecAA" %in% textModif) {          ## remove preceeding and following AAs eg '[KR].AGVIFPVGR.[ML]' => 'AGVIFPVGR'
    annN <- lapply(annN, function(x) sub("^\\[[[:upper:]]+\\]\\.","", sub("\\.\\[[[:upper:]]+\\]$","", x)) ) }

  if("rmTerminalDigit" %in% textModif) {       ## remove terminal digit (eg from charge-state)
    annN <- lapply(annN, function(x) sub("[[:digit:]]+$","", x) ) }
  if(debug) {message(fxNa," fPL1 ++"); fPL1 <- list(x=x,y=y,z=z,columnNa=columnNa,listNa=listNa,shortNa=shortNa,annN=annN,textModif=textModif,all=all)} #xAnn=xAnn,xQua=xQua,yAnn=yAnn,yQua=yQua,  

  ## combine redundant IDs ?
  chUni <- lapply(annN, duplicated)
  ##
  if(datOK) {
    ## IDs : combine from data-sets   
    rmRed <- list(x=x, y=y, z=z)
    chLe <- sapply(rmRed, length) >0
    
    ### CHECK THIS WHEN FINISHED modifying combineRedundLinesInList()  !!
    allID3 <- lapply(rmRed, function(w) wrMisc::naOmit(w[[listNa[length(listNa)]]][,columnNa]))
    if(isTRUE(all)) allID <- unique(unlist(allID3)) else {
      k <- which(sapply(allID3, length) >0)  
      allID <- if(length(k) >1) intersect(allID3[[k[1]]], allID3[[k[2]]]) else NULL
      if(length(k) >2) for(i in k[-1:-2]) allID <- intersect(allID, allID3[[k[i]]])
    }
    if(debug) {message(fxNa," fPL1b ++"); fPL1b <- list()} #  
    if(length(allID) <1) { datOK <- FALSE
      if(debug) message(fxNa,"Found NO COMMON IDs !   Nothing to do ..") }                                                                                                               
  }
          
  if(datOK) {
    if(debug) {message(fxNa," fPL2"); fPL2 <- list(all=all,rmRed=rmRed,allID=allID,x=x,y=y,z=z,columnNa=columnNa,listNa=listNa,shortNa=shortNa,annN=annN,textModif=textModif,namesXY=namesXY)} #xAnn=xAnn,xQua=xQua,yAnn=yAnn
    ## extract column names of quantitation data (to fuse)
    colNaAnn <- lapply(rmRed, function(w) colnames(w[[listNa[length(listNa)]]]))
    colNa <- lapply(rmRed, function(w) colnames(w[[listNa[1]]]))
    colInd <- cumsum(sapply(colNa, length))
    colInd <- cbind(beg=c(1, 1 +colInd[-length(colNa)]), end=colInd)    
    
    ## fuse quantitative data (add columns)
    useLst <- which(sapply(rmRed, length) >0)       # which datasets contain data for fusing
    if(debug) {message(fxNa," fPL2c"); fPL2c <- list()} 
    j <- 1
    for(i in useLst) {
      if(j==1) { out <- rmRed[[i]]
        dim2 <- lapply(out, dim)
        arr3dim <- sapply(dim2, function(w) if(length(w) ==3) all(w[1:2]==dim(rmRed[[i]][[listNa[1]]])) else FALSE )     # all elements with 3 dims & same nrow & ncol as listNa[1]
        mat2dim <- sapply(dim2, function(w) if(length(w) ==2) all(w[1:2]==dim(rmRed[[i]][[listNa[1]]])) else FALSE )     # all elements with 2 dims & same nrow & ncol as listNa[1]
        chAnn <- names(mat2dim) %in% listNa[length(listNa)]
        if(any(chAnn, na.rm=TRUE)) mat2dim[which(chAnn)] <- FALSE
        if(debug) {message(fxNa,"Matrix-type elements found for fusing : ",wrMisc::pasteC(names(out)[which(mat2dim)], quoteC="'")) 
          message(fxNa,"3-dim arrays elements found for fusing : ",wrMisc::pasteC(names(out)[which(arr3dim)], quoteC="'"))  }
      if(debug) {message(fxNa," fPL2d  i=",i,"  j=",j); fPL2d <- list(all=all,rmRed=rmRed,allID=allID,x=x,y=y,z=z,out=out,columnNa=columnNa,listNa=listNa,shortNa=shortNa,annN=annN,textModif=textModif,namesXY=namesXY,useLst=useLst,colInd=colInd,colNa=colNa,arr3dim=arr3dim)
      }       
      ## fuse quantitative data (add columns)
      if(sum(arr3dim) >0) for(k in which(arr3dim)) out[[k]] <- array(dim=c(length(allID), sum(sapply(colNa, length)), dim(rmRed[[i]][[k]])[3]), 
        dimnames=list(allID, paste0(rep(namesXY,sapply(colNa, length)),".",unlist(colNa)), dimnames(rmRed[[i]][[k]])[[3]]))
      if("notes" %in% names(out)) out$notes <- matrix(out$notes, ncol=1, dimnames=list(names(out$notes),namesXY[i]))
      if("sampleSetup" %in% names(out)) out$sampleSetup <- list(out$sampleSetup$groups)           # very minimal fusion ... ($groups ONLY !!)
      if("quantNotes" %in% names(out)) out$quantNotes <- list(out$quantNotes) }
      ## prepare for fusing quantitative data (columns)
      if(j==1) { tmp <- matrix(NA, nrow=length(allID), ncol=sum(sapply(colNa, length)), dimnames=list(allID, paste0(rep(namesXY,sapply(colNa, length)),".",unlist(colNa))))  ## initialize
        if(sum(mat2dim) >1) for(k in which(mat2dim))  out[[k]] <- tmp
        if(sum(arr3dim) >0) {
          for(k in which(arr3dim)) tmp <- array(NA, dim=c(length(allID), sum(sapply(colNa, length)), arr3dim[[k]]), 
            dimnames=list(allID, paste0(rep(namesXY,sapply(colNa, length)),".",unlist(colNa)), dimnames(dim2[[k]])[[3]]))  ## initialize
        }       
      }      
      ## assign
      inLi <- wrMisc::naOmit(match(allID, rmRed[[i]][[listNa[length(listNa)]]][,columnNa[i]]))                # where current fits to allID
      if(debug) message(fxNa,"i=",i," extract ",length(inLi)," common out of current ",nrow(rmRed[[i]][[listNa[length(listNa)]]]),"")
      
      for(k in which(mat2dim)) out[[k]][which(allID %in% rmRed[[i]][[listNa[length(listNa)]]][,columnNa[i]]), colInd[i,1]:colInd[i,2]] <- rmRed[[i]][[k]][ inLi, ]
      if(j==1) out[[listNa[length(listNa)]]] <- out[[listNa[length(listNa)]]][inLi, ]       # annotation

      if(debug) {message(fxNa," fPL2e  i=",i,"  j=",j); fPL2e <- list()}
      if(sum(arr3dim) >0) for(k in which(arr3dim)) if(j >1 && names(out)[k] %in% names(rmRed[[i]])) {   # possibly also check if name of 3rd dim is consistent ..
        if(length(dim(rmRed[[i]][[k]])) ==3) out[[k]][which(allID %in% rmRed[[i]][[listNa[length(listNa)]]][,columnNa[i]]), colInd[i,1]:colInd[i,2],] <- rmRed[[i]][[k]][ inLi,, ] }
      

      ## also try fusing annotation data ?  
            
      ## fuse other (notes, sampleSetup, quantNotes)
      if(debug) message(fxNa," fPL3a ")
      if(j >1 && "notes" %in% names(rmRed[[i]])) out$notes <- cbind(out$notes, unlist(rmRed[[i]]$notes)[match(rownames(out$notes), names(rmRed[[i]]$notes))])
      if(j >1 && "sampleSetup" %in% names(rmRed[[i]])) out$sampleSetup[[j]] <- rmRed[[i]]$sampleSetup$groups
      if(j >1 && "quantNotes" %in% names(rmRed[[i]])) out$quantNotes[[j]] <- rmRed[[i]]$quantNotes
      j <- j +1      
    }
  }
  out
}
    
