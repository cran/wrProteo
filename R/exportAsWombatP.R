#' Export As Wombat-P Set Of Files
#'
#' This function allows exporting objects created from wrProteo to the format of Wombat-P \href{https://github.com/wombat-p}{Wombat-P}.
#'  
#' 
#' @param wrProtObj (list produced by any import-function from wrProteo) object which will be exported as Wombat-P format
#' @param path (character) the location where the data should be exorted to
#' @param combineFractions (\code{NULL} or character (length=1)) if not \code{NULL} this assigns the method how multiple farctions should be combined
#'   (at this point only the method 'mean' is implemented)
#' @param silent (logical) suppress messages
#' @param debug (logical) display additional messages for debugging
#' @param callFrom (character) allows easier tracking of messages produced
#' @return This function creates a set of files (\code{README.md}, \code{test_params.yml}), plus a sud-directory containig file(s) (\code{stand_prot_quant_method.csv}); finally the function returns  (\code{NULL}),
#' @seealso  \code{\link{readMaxQuantFile}}, \code{\link{readProteomeDiscovererFile}}; \code{\link[wrMisc]{moderTestXgrp}} or \code{\link[wrMisc]{moderTest2grp}}
#' @examples
#' path1 <- system.file("extdata", package="wrProteo")
#' fiNa <- "proteinGroupsMaxQuant1.txt.gz"
#' specPr <- c(conta="conta|CON_|LYSC_CHICK", mainSpecies="YEAST", spike="HUMAN_UPS")
#' dataMQ <- readMaxQuantFile(path1, file=fiNa, specPref=specPr, tit="tiny MaxQuant")
#' 
#' exportAsWombatP(dataMQ, path=tempdir())
#' @export
exportAsWombatP <- function(wrProtObj, path=".", combineFractions="mean", silent=FALSE, debug=FALSE, callFrom=NULL) {
  ## This function allows exporting objects created from wrProteo to the format of Wombat-P \href{https://github.com/wombat-p}{Wombat-P}.

  fxNa <- wrMisc::.composeCallName(callFrom, newNa="exportAsWombatP")
  if(isTRUE(debug)) silent <- FALSE
  if(!isTRUE(silent)) silent <- FALSE
  subDir <- "dev"
  if(identical(subDir,".")) subDir <- NULL
  ## check path
  msg1 <- "Invalid argument 'path' setting to '.' (current working path)" 
  if(!identical(path,".")) {
    if(length(path) < 1) {path <- "."; if(!silent) message(msg1)}
  } else {
    if(length(path) >1) { path <- path[1]
      if(!silent) message(fxNa, "Path should be length=1, using 1st")} 
    if(is.na(path))  { path <- path[1]
      if(!silent) message(fxNa, msg1) }
    if(path !="." && !dir.exists(path))  {
      ch1 <- try(dir.create(path), silent=TRUE)
      if(inherits(ch1, "try-error")) {  
        warning(fxNa, "Unable to create 'path' (",path,") as requested, setting to '.' (current working path) ")
        path <- "."
      } else if(debug) message(fxNa,"Sucessfully created path '",path,"'")
    } else if(debug) message(fxNa,"Path seems OK")
  }

  
  datOK <- length(wrProtObj) >1 && is.list(wrProtObj) && all(c("quant","annot") %in% names(wrProtObj)) && length(wrProtObj$quant) >1
  if(!datOK) warning(fxNa,"This does NOT look like a valid object from wrProteo, abortiong export-function")

  ## check creating folder /dev
  if(!dir.exists(file.path(path,subDir))) { ch1 <- try(dir.create(file.path(path,subDir)), silent=TRUE)
    if(inherits(ch1, "try-error")) {datOK <- FALSE; message(fxNa, "Unable to create subdir '",subDir,"', abortiong export-function")}  
  }

  ## write .md
  if(datOK) {
    txt <- c("# Processed (PRIDE) data sets",
      "Collection of data set(s) that has/have been re-analyzed, imported using [wrProteo](https://CRAN.R-project.org/package=wrProteo) and exported as [WOMBAT-P](https://github.com/wombat-p/WOMBAT-pipelines).",
      "",
      "Each data set is represented by a folder named as the ProteomeXChange accession number.",
      paste0("Using version wrProteo-",as.character(utils::packageVersion("wrProteo")))
    )
    ch1 <- try(writeLines(txt, sep="\n", con=file.path(path, "README.md")), silent=silent)
    if(inherits(ch1, "try-error")) {datOK <- FALSE; warning(fxNa,"Failed writing README file, check for permissions to write in this directory") 
    } else if(debug) message(fxNa,"Wrote README.md file succesfully ")

  }

  ## write testparams.yml
  if(datOK) {
    txt <- matrix(c("params: ","",
      "  enzyme: ",sub("^NT=","",sub(";[[:upper:]]+.+","",wrProtObj$sampleSetup$sdrfDat$comment.cleavage.agent.details.[1])),  #Trypsin/P
      "  fions: ",NA,  #b
      "  rions: ",NA,  #y
      "  isotope_error_range: ",NA,   #1
      "  add_decoys: ",NA,   #true
      "  num_hits: ",NA,   #1
      "  miscleavages: ",NA,   #1
      "  min_precursor_charge: ",NA,   #2
      "  max_precursor_charge: ",NA,   #3
      "  min_peptide_length: ",NA,   #5
      "  max_peptide_length: ",NA,   #15
      "  max_mods: ",NA,   #4
      "  ident_fdr_psm: ",NA,   #0.01
      "  ident_fdr_peptide: ",NA,   #0.01
      "  ident_fdr_protein: ",NA,   #0.01
      "  enable_match_between_runs: ",NA,   #True
      "  protein_inference: ",NA,   #unique
      "  quantification_method: ","precursor",
      "  summarization_proteins: ",NA,
      "  min_num_peptides: ",NA,
      "  summarization_psms: ",NA,   #sum_abs
      "  quant_transformation: ","log",
      "  normalization_method: ",wrProtObj$notes["normalizeMeth"],
      "  run_statistics: ","true",
      "  fdr_method: ","qvalue",
      "  fdr_threshold: ",NA,
      "rawfiles: ","None",
      "fastafile: ","None"
    ), ncol=2, byrow=TRUE)
    txt <- paste0(txt[,1],txt[,2])
    ## possible to complete based on sdrf or other info ?  

    ch1 <- try(writeLines(txt, sep="\n", con=file.path(path, "testparams.yml")), silent=silent)
    if(inherits(ch1, "try-error")) {datOK <- FALSE; warning(fxNa,"Failed writing 'testparams.yml' file, check for permissions to write in this directory")}
  }

  ## write main quantitation data to 'dev'
  if(datOK) {
    expDat <- wrProtObj$quant
    colnames(expDat) <- paste0("abundance_",colnames(expDat))
    expDat <- data.frame( protein_group=paste(wrProtObj$annot[,"Accession"], wrProtObj$annot[,"EntryName"],sep="|"), wrProtObj$annot, expDat)
    ch1 <- try(utils::write.csv(expDat, file=file.path(path, subDir, paste0("std_p",if(isTRUE(wrProtObj$notes["identType"]=="peptide")) "ep" else "rot","_quant_",  ".csv")) ,
      row.names=FALSE, quote=FALSE), silent=silent)
    if(inherits(ch1, "try-error")) {datOK <- FALSE; warning(fxNa,"Failed writing data to file, check for permissions to write in this directory")}
  }
  if(datOK && !silent) message(fxNa,"Succesfully exported data as Wombat-P type of files")
}

  