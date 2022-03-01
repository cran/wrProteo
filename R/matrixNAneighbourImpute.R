#' Imputation of NA-values based on non-NA replicates
#' 
#' It is assumed that \code{NA}-values appear in data when quantitation values are very low (as this appears eg in quantitative shotgun proteomics). 
#' Here, the concept of (technical) replicates is used to investigate what kind of values appear in the other replicates next to NA-values for the same line/protein.
#' Groups of replicate samples  are defined via argument \code{gr} which descibes the columns of \code{dat}). 
#' Then, they are inspected for each line to gather NA-neighbour values (ie those values where NAs and regular measures are observed the same time).
#' Eg, let's consider a line contains a set of 4 replicates for a given group. Now, if 2 of them are \code{NA}-values, the remaining 2 non-\code{NA}-values will be considered as NA-neighbours.
#' Ultimately, the aim is to replaces all \code{NA}-values based on values from a normal distribution ressembling theire respective NA-neighbours.
#' 
#' By default a histogram gets plotted showing the initial, imputed and final distribution to check the global hypothesis that \code{NA}-values arose 
#' from very low measurements and to appreciate the impact of the imputed values to the overall final distribution.
#' 
#' @details
#' There are a number of experimental settings where low measurements may be reported as \code{NA}. 
#' Sometimes an arbitrary defined baseline (as 'zero') may provoke those values found below being unfortunately reported as \code{NA} or as 0 (in case of MaxQuant). 
#' In quantitative proteomics (DDA-mode) the presence of numerous high-abundance peptides will lead to the fact that a number of less 
#' intense MS-peaks don't get identified properly and will then be reported as \code{NA} in the respective samples, 
#' while the same peptides may by correctly identified and quantified in other (replicate) samples.
#' So, if a given protein/peptide gets properly quantified in some replicate samples but reported as \code{NA} in other replicate samples 
#' one may thus speculate that similar values like in the successful quantifications may have occored. 
#' Thus, imputation of \code{NA}-values may be done on the basis of \code{NA}-neighbours.
#'  
#'  
#'  
#' When extracting \code{NA}-neighbours, a slightly more focussed approach gets checked, too, the 2-\code{NA}-neighbours : In case a set of replicates for a given protein 
#' contains at least 2 non-\code{NA}-values (instead of just one) it will be considered as a (min) 2-\code{NA}-neighbour as well as regular \code{NA}-neighbour. 
#' If >300 of these (min) 2-\code{NA}-neighbours get found, they will be used instead of the regular \code{NA}-neighbours.
#' For creating a collection of normal random values one may use directly the mode of the \code{NA}-neighbours (or 2-\code{NA}-neighbours, if >300 such values available).
#' To do so, the first value of argument \code{avSdH} must be set to \code{NA}. Otherwise, the first value \code{avSdH} will be used as quantile of all data to define the mean   
#' for the imputed data (ie as \code{quantile(dat, avSdH[1], na.rm=TRUE)}). The sd for generating normal random values will be taken from the sd of all  \code{NA}-neighbours (or 2-\code{NA}-neighbours)
#' multiplied by the second value in argument \code{avSdH} (or \code{avSdH}, if >300 2-\code{NA}-neighbours), since the sd of the \code{NA}-neighbours is usually quite high.
#' In extremely rare cases it may happen that no \code{NA}-neighbours are found (ie if \code{NA}s occur, all replicates are \code{NA}). 
#' Then, this function replaces \code{NA}-values based on the normal random values obtained as dscribed above.
#' 
#' @param dat (matrix or data.frame) main data (may contain \code{NA})
#' @param gr (character or factor) grouping of columns of 'dat', replicate association
#' @param imputMethod (character) choose the imputation method (may be 'mode2'(default), 'mode1', 'datQuant', 'modeAdopt' or 'informed')
#' @param retnNA (logical) decide (if =\code{TRUE}) only NA-substuted data should be returned, or if list with $data, $nNA, $NAneighbour and $randParam should be returned
#' @param avSdH (numerical,length=2) population characteristics 'high' (mean and sd) for >1 \code{NA}-neighbours (per line)
#' @param NAneigLst (list) option for repeated rounds of imputations: list of \code{NA}-neighbour values can be furnished for slightly faster processing 
#' @param plotHist (character or logical) decide if supplemental figure with histogram shoud be drawn, the details 'Hist','quant' (display quantile of originak data), 'mode' (display mode of original data) can be chosen explicitely
#' @param xLab (character) label on x-axis on plot
#' @param xLim (numeric, length=2) custom x-axis limits
#' @param yLab (character) label on y-axis on plot
#' @param yLim (numeric, length=2) custom y-axis limits
#' @param tit (character) title on plot
#' @param figImputDetail (logical) display details about data (number of NAs) and imputation in graph (min number of NA-neighbours per protein and group, quantile to model, mean and sd of imputed)
#' @param seedNo (integer) seed-value for normal random values
#' @param silent (logical) suppress messages
#' @param callFrom (character) allow easier tracking of messages produced
#' @param debug (logical) supplemental messages for debugging
#' @return This function returns a list with \code{$data} .. matrix of data where \code{NA} are replaced by imputed values, \code{$nNA} .. number of \code{NA} by group, \code{$randParam} .. parameters used for making random data 
#' @seealso this function gets used by \code{\link{testRobustToNAimputation}}; estimation of mode \code{\link[wrMisc]{stableMode}}; detection of NAs \code{\link[stats]{na.fail}}
#' @examples
#' set.seed(2013)
#' datT6 <- matrix(round(rnorm(300)+3,1), ncol=6, dimnames=list(paste("li",1:50,sep=""),
#'   letters[19:24]))
#' datT6 <- datT6 +matrix(rep(1:nrow(datT6), ncol(datT6)), ncol=ncol(datT6))
#' datT6[6:7, c(1,3,6)] <- NA
#' datT6[which(datT6 < 11 & datT6 > 10.5)] <- NA
#' datT6[which(datT6 < 6 & datT6 > 5)] <- NA
#' datT6[which(datT6 < 4.6 & datT6 > 4)] <- NA
#' datT6b <- matrixNAneighbourImpute(datT6, gr=gl(2,3))
#' head(datT6b$data)
#' @export
matrixNAneighbourImpute <- function(dat, gr, imputMethod="mode2", retnNA=TRUE, avSdH=c(0.1,0.5), NAneigLst=NULL, 
  plotHist=c("hist","mode"), xLab=NULL, xLim=NULL, yLab=NULL, yLim=NULL, tit=NULL, figImputDetail=TRUE, 
  seedNo=NULL, silent=FALSE, callFrom=NULL, debug=FALSE){
  ## replace NA values based on group neigbours (based on grouping of columns in gr), overall assumption of close to Gaussian distrib
  ## return matrix including imputed values or list of final & matrix with number of imputed by group
  ## 'batch-mode' (iterated runs) furnish NAneigLst (with $nNaNei, $charAll, $all.lm or $linMod, $NAneighbour), avSdH (postition 3+ for medMode)
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="matrixNAneighbourImpute")
  if(!isTRUE(silent)) silent <- FALSE
  if(isTRUE(debug)) silent <- FALSE
  datOK <- TRUE
  tx1 <- "; invalid entry - nothing to do"
  if(datOK) { if(length(unique(gr))==length(gr)) { datOK <- FALSE
    msg <- "Argument 'gr' does not designate any replicates !  .. ie NA-neighbours can't be determined - nothing to do"} } 
  if(length(dat) <1) { datOK <- FALSE
    msg <- c("'dat' must be list, S3-object, matrix or data.frame with min 2 rows and min 2 cols",tx1)}   
  if(datOK) if(is.list(dat)) {
    if("quant" %in% names(dat)) {dat <- dat$quant
      if(debug) message(fxNa,"'dat' is list, using element 'quant' as 'dat'")
    } else { datOK <- FALSE
      msg <- c("'dat' is list but does NOT contain element named 'quant'",tx1)}
  } else { if(any(dim(dat) < c(2,1))) { datOK <- FALSE
    msg <- c("'dat' matrix or data.frame with min 2 rows and min 2 cols",tx1)}} 
  if(datOK) { if(length(gr) != ncol(dat)) { datOK <- FALSE
    msg <- c("Number of columns of 'dat' must match length of 'gr'",tx1)}}  
        
  if(datOK) {
    if(is.data.frame(dat)) dat <- as.matrix(dat)
    if(!is.factor(gr)) gr <- as.factor(gr)
    if(is.null(xLab)) xLab <- "Values"
    if(length(wrMisc::naOmit(imputMethod)) <1) { imputMethod <- "mode2"
      if(!silent) message(fxNa,"Invalid entry for 'imputMethod' setting to default")}
    chMeth <- imputMethod %in% c("datQuant", "medmode", "mode1", "mode2", "mode3", "modeadopt","informed")
    if(!chMeth) { 
      if(!silent) message(fxNa,"Unknown method '",imputMethod,"' as entry for 'imputMethod', setting to default")
      imputMethod <- "mode2"}
    ## extract elements for 'batch-mode' 
    if(debug) message(fxNa,"Extract elements for 'batch-mode'    mn0")   
    if(is.list(NAneigLst) & length(NAneigLst) >1) { 
      nNaNei <- NAneigLst$nNaNei; charAll <- NAneigLst$charAll;
      medMod <- NAneigLst$medMod 
      all.lm <- if("all.lm" %in% names(NAneigLst)) NAneigLst$all.lm else NAneigLst$linMod
      NAneighbour <- NAneigLst$NAneighbour
    } else { NAneighbour <- NAneigLst <- nNaNei <- medMod <- all.lm <- charAll <- NULL }    # initialize
    if(length(seedNo) >1) { seedNo <- seedNo[1]
      if(!silent) message(fxNa,"Invalid entry for argument 'seedNo', it may be single integer or NULL, setting to NULL")}  
    if(length(seedNo) >0) if(any(is.na(seedNo))) seedNo <- NULL
    if(is.logical(plotHist)) { plotHist <- if(identical(TRUE,plotHist)) c("hist","quant","mode") else NULL}
    medMod <- NULL
    
    ## main
    isNA <- is.na(dat)
    chNA <- any(isNA)
    if(debug) { message(fxNa,"Starting main,  mn0");  mn0 <- list(dat=dat,gr=gr,isNA=isNA,chNA=chNA,imputMethod=imputMethod,NAneigLst=NAneigLst,NAneighbour=NAneighbour)}
  
    if(!chNA) {            
      ## no NAs, nothing to impute ...
      if(debug) message(fxNa,"No NAs, nothing to impute ...")
      if("hist" %in% plotHist) {graphics::hist(dat, br="FD", border=grDevices::grey(0.85), col=grDevices::grey(0.92), xlab=xLab, las=1, main=tit)
        graphics::mtext("No NA-replacement needed  ", adj=1, cex=0.6, line=-0.3)    
        graphics::mtext(paste("  n=",length(dat)), side=3, line=-0.3, cex=0.55, adj=0,col=grDevices::grey(0.3)) }
      return( if(isTRUE(retnNA)) list(data=dat, nNA=0, NAneighbour=NULL, randParam=NULL) else dat)
    } else {
      modNa <- NULL
      if(length(NAneighbour) <1 & (length(grep("mode", imputMethod)) >0 | any(c("informed","datQuant") %in% imputMethod))) {   ## all methods using mode need NA-neighbours ...
        if(debug) message(fxNa," mn1")
        if(length(NAneighbour) <1) NAneighbour <- isolNAneighb(dat, gr, iniCheck=FALSE, silent=silent, callFrom=fxNa)     #
        if(debug) {message(fxNa," mn2"); mn2 <- list(NAneighbour=NAneighbour)}
        chLast <- length(NAneighbour[[length(NAneighbour)]])
        if(chLast==0) NAneighbour <- NAneighbour[-1*length(NAneighbour)]     # remove last empty field if empty (as usual)

      } else { nNA <- sum(isNA)}
      if(debug) {message(fxNa," mn3"); mn3 <- list(dat=dat,gr=gr,isNA=isNA,chNA=chNA,chLast=chLast,imputMethod=imputMethod,NAneigLst=NAneigLst,NAneighbour=NAneighbour)}
      nNaNei <- sapply(NAneighbour, length)  
      datIni <- dat
        
      ## check if sufficient NAs for mode-based methods
      if(debug) message(fxNa,"Checking if sufficient NAs for mode-based methods")
      if(sum(nNaNei) <10 & length(grep("mode", imputMethod)) >0) {
        ##number of NA neighb not yet known#
        if(!silent) message(fxNa,"Only ",sum(nNaNei)," NA-neighbour values available, ie insufficient to calculate representative mode, using instead 10%ile of global distribution")
        imputMethod <- "datQuant"
        avSdH[1] <- 0.1
      }
      if(debug) {message(fxNa," mn4") }
        
      ## IMPUTATIONS
      ## choose method
      if("datQuant" %in% imputMethod) {   # quantile of data
        if(debug) message(fxNa,"Starting method 'datQuant'")
        if(length(avSdH) >2 & !is.na(avSdH[3])) {useQu <- avSdH[3]} else {   # 3rd value (if specified) may be used as custom mean for normal distrib (instead of determining as xth quantile) 
          if(is.na(avSdH[1])) { avSdH[1] <- 0.1           # check quantile value (ie where to check full data)
            if(!silent) message(fxNa," avSdH not valid, using 10%quantile instead")}
          useQu <- stats::quantile(dat, avSdH[1], na.rm=TRUE) }
        if(length(seedNo) ==1) set.seed(seedNo)
        randVa <- signif(stats::rnorm(sum(isNA), useQu, avSdH[2]),5)
        plotHist <- unique(c(plotHist,"quantile"))     # for ploting quantile-guides
        msg <- paste(signif(avSdH[1],3),"quantile, ie mean=",signif(useQu,4),"and sd=", signif(avSdH[2],4))
      }
      if("medmode" %in% tolower(imputMethod)) {        # whatever is lowest: global median or global mode of all NA neighbours     
        if(debug) message(fxNa,"Starting method 'medmode'")
        if(length(medMod) <1) medMod <- c(med=stats::median(unlist(NAneighbour)), mod=wrMisc::stableMode(unlist(NAneighbour), method="density", callFrom=fxNa, silent=silent))
        if(debug) {message(fxNa," mn5a") }
        if(length(seedNo) ==1) set.seed(seedNo)
        randVa <- signif(stats::rnorm(sum(isNA), min(medMod), avSdH[2]), 5)      
        msg <- paste("mean=",signif(min(medMod),4),"and sd=", signif(avSdH[2],4))
      }
      if(any(c("mode1","mode3","informed") %in% imputMethod)) {      # global mode of all NA neighbours
        if(debug) message(fxNa,"Prepare for methods 'mode1','mode2' and 'informed'")
        if(length(medMod) <1) medMod <- wrMisc::stableMode(unlist(NAneighbour), method="density",silent=silent,callFrom=fxNa) 
        if(length(seedNo) ==1) set.seed(seedNo)
        randVa <- signif(stats::rnorm(sum(isNA), medMod, avSdH[2]),5)
        msg <- paste("mode=",signif(medMod,4),"and sd=", signif(avSdH[2],4))  # correct ??
      }
        if(debug) {message(fxNa," mn5b") }
          
      if("mode2" %in% imputMethod) {      # selective mode of NA of '2-NA-neighbours' if n.2NA > 300
        if(debug) message(fxNa,"Starting methods 'mode2'")
        if(length(medMod) <1) medMod <- wrMisc::stableMode(if(sum(sapply(NAneighbour[-1], length)) >300) unlist(NAneighbour[-1]) else unlist(NAneighbour), method="density", callFrom=fxNa,silent=silent) 
        if(length(seedNo) ==1) set.seed(seedNo)
        randVa <- signif(stats::rnorm(sum(isNA), medMod, avSdH[2]),5)
        msg <- paste("mean=",signif(medMod,4),"and sd=", signif(avSdH[2],4))
      }
      if("informed" %in% imputMethod) {      # informed about min abundance in line with any NA, use mean of line-min and non-biased NA-neigh random value
        ## determine "biased" view : min per line with any NA 
        ## still has problem when one all repl of one grp NA and highly abundant in other group
        if(debug) message(fxNa,"Starting method 'informed'")
        NAliCo <- rowSums(isNA)   
        NAli <- which(NAliCo >0 & NAliCo < ncol(dat))             # lines with NA but not all NA
        NAmin <- rep(NA, nrow(dat))
        if(debug) {message(fxNa," mn5c") }
        
         if(length(NAli) >0) {
          ## prepare
          indC <- matrix(0, nrow=nrow(isNA), ncol=ncol(isNA))
          indC[NAli,] <- 1                              # matrix (full dim) indicating where min-informed correction can be done
          indC <- 0 + isNA + indC                       # is 2 if NA and eligible to min-informed
          isNA2 <- 0 + isNA
          isNA2[which(indC==2)] <- 2                    # matrix (full dim), 2 .. NA eligible to informed cor; 1.. NA not elig
          if(debug) {message(fxNa," mn5c2") }
          
          NAmin <- apply(dat[NAli,], 1, min, na.rm=TRUE)   # min of line with any (but not all) NA (for biased estimate of NA-replacement
          ReMa <- matrix(rep(NAmin, ncol(isNA)), ncol=ncol(isNA))    # matrix of row-min values
          ReMa[which(indC[NAli,] !=2)] <- NA                         # leave only value at position eligible to min-informed correction
          useInd <- which(isNA2[which(isNA2 >0)] >1)    # which randVa is eligible to min-informed cor
          ## compare to level of regualar NA-neighb
          corF <- wrMisc::naOmit(as.numeric(ReMa)) - medMod
          chUp <- corF >0
          if(any(chUp)) corF[which(chUp)] <- 0            # don't use if higher than mode of NA-neighb
          ##  apply min-informed cor
          randVa[useInd] <- randVa[useInd] +corF          # mean of imputed and non-biased rand Va
        }
        
      }    
      if("modeadopt" %in% tolower(imputMethod)) {   # flexible/adopt
        ## median and mode of NA-neighbours all/ by group :
        if(debug) {message(fxNa,"Starting method 'modeadopt'","  mn5d") }
        if(length(charAll) <1) { charAll <- c(mean=mean(unlist(NAneighbour)), med=stats::median(unlist(NAneighbour)), 
            mode=as.numeric(wrMisc::stableMode(unlist(NAneighbour), method="density", silent=TRUE)))
          charAll <- rbind(charAll,cbind(sapply(NAneighbour, mean), sapply(NAneighbour, stats::median), sapply(NAneighbour, wrMisc::stableMode, method="density", silent=TRUE, callFrom=fxNa)))}
        modNa <- charAll[1,3]
        ## model for dynamic NA-neighbour mean  (hypoth for further interpolating)      
        nRa <- as.integer(sub("n","", names(NAneighbour) ))     
        if(length(all.lm) <1) all.lm <- stats::lm(a~n, data=data.frame(n=rep(nRa,nNaNei), a=unlist(NAneighbour)))
        pVaSlo <- summary(stats::aov(all.lm))[[1]][1,"Pr(>F)"]                  # not  < 0.05
        if(pVaSlo < 0.1 & stats::coef(all.lm)[2] <0) {
          ## dynamic NA-substit depending on number of NAs per prot&repl
          prLev <- stats::predict(all.lm, new=data.frame(n=1:(max(nRa) +1)))
          ranOff <- prLev -charAll[1,1]                  # chose to subtract mode of all NA-neighbours (as reference)
          nNaGrp <- wrMisc::rowGrpNA(dat, gr)            # check : sum(nNaGrp) == sum(isNA)    # here 2845        
          ranOff <- charAll[1,1] - rep(ranOff[nNaGrp], nNaGrp[which(nNaGrp >0)])        # final offset in order for is.na(dat)
          msg <- paste("mean=",wrMisc::pasteC(signif(prLev,4)),"(for",wrMisc::pasteC(1:(1+max(nRa))),"NAs) and sd=", signif(avSdH[2],4))
          if(!silent) message(fxNa,"Substituting dynamically based on mean per number of NAs")
        } else {
          ## constant NA-substit, use min of mean, median and mode
          ranOff <- charAll[1,1] - rep(min(charAll[1,]), sum(isNA))
          msg <- paste("mean=",signif(min(charAll[1,]),4),"and sd=", signif(avSdH[2],4))
          if(!silent) message(fxNa,"Substituting based on ",c("mean","median","mode")[which.min(charAll[1,])]," of all ",sum(nNaNei)," NA-neighbours")
        }
        if(length(seedNo) ==1) set.seed(seedNo)
        randVa <- signif(stats::rnorm(sum(isNA), charAll[1,1], avSdH[2]) -ranOff, 5)    # initial global value (median for all NA-neighb) + offset/correction
      }
      
      ## replace NAs
      if(debug) {message(fxNa," mn6") }
      dat[which(isNA)] <- randVa
      chImp <- stats::quantile(datIni, c(0.05,0.15), na.rm=TRUE)
      chImp <- c(mean(randVa) < chImp[1], mean(randVa) > chImp[2])
  
      msg <- list(li1=c(" n.woNA=",sum(!isNA),", n.NA =",sum(isNA)),
        li2=c("impute based on",paste0("'",imputMethod,"'"),"using",msg),
        li3=if(any(chImp)) c("note mean for impuation is ",if(chImp[1]) "below 0.05 " else "above 0.15", "quantile !!") )
      if(!silent) message(fxNa, paste(sapply(msg, paste, collapse=" "), collapse="\n    "))
  
      ## FIGURE  
      if("hist" %in% tolower(plotHist)) {              
        hi1 <- graphics::hist(as.numeric(dat), breaks="FD", col=grDevices::grey(0.9), border=grDevices::grey(0.8),
          xlab=xLab, ylab=yLab, las=1, ylim=yLim, main=paste(tit,"at NA-replacement"))     #xlim=xLim, #  grey cols (final distr)      
        colPanel <- c(grDevices::grey(0.6), grDevices::rgb(0,0.7,0,0.6), grDevices::rgb(0.15,0.15,0.7,0.7), grDevices::rgb(0.7,0.5,0.2,0.6), grDevices::rgb(0.8,0.2,0.7,0.7))
        graphics::hist(datIni, breaks=hi1$breaks, border=grDevices::grey(0.75), col=grDevices::rgb(0.1,1,0.1,0.15), add=TRUE)                  # orig data in green
        if(length(randVa) >5) graphics::hist(randVa, br=hi1$breaks, border=grDevices::grey(0.75), col=grDevices::rgb(0,0,0.7,0.2), add=TRUE)           # add purple hist to
        nextLi <- -1.7
        if(any(c("quant","quantile") %in% plotHist)) {
          graphics::abline(v=stats::quantile(datIni, c(0.05,0.1,0.15), na.rm=TRUE), lty=2, col=c(colPanel[4:5],"tomato4"))
          nextLi <- nextLi -c(0, 0.5, 1.1)
          graphics::mtext(paste(" - - ",c(0.05,0.1,0.15),"quantile (initial data)"), col=c(colPanel[4:5],"tomato4"), cex=0.7, adj=0, line=nextLi, side=3)
          nextLi <- min(nextLi) -0.6 }
  
        if(length(modNa) >0) {    # display mode
          yLim <- signif(graphics::par("usr")[3:4], 3)             # current y-limits
          if(any(c("mode") %in% plotHist)) {
            graphics::mtext(paste(" (arrow) mode of", NULL, " NA-neighbours :", signif(modNa,4)), col="sienna2", cex=0.7, adj=0, line=nextLi, side=3)
            graphics::arrows(modNa, yLim[1] - (yLim[2] -yLim[1])/18, modNa, 0, length=0.1, col="sienna2",lwd=2) }
        }    
        
        if(isTRUE(figImputDetail)) graphics::mtext(paste(sapply(msg[1:2],paste,collapse=" "),collapse="\n "), side=3, line=-1.2, cex=0.75, adj=0,col=grDevices::grey(0.3))
        graphics::legend("topright",c("final","initial","imputed"), col=colPanel, text.col=colPanel, cex=0.9, seg.len=0.3, lwd=4)
      }
      # prev return NAneighbour= 
      return(if(isTRUE(retnNA)) list(data=dat, nNA=sum(isNA) , randParam=imputMethod, NAneigLst=list(NAneighbour=NAneighbour, nNaNei=nNaNei, medMod=medMod, charAll=charAll, linMod=all.lm)) else dat)
    } 
  } else { if(!silent) message(fxNa,msg)}
  dat }

  
  #' @export  
.imputeNA <- function(dat,gr,impParam,exclNeg=TRUE,inclLowValMod=TRUE) {
  ## basic NA imputation for 'dat' using 'impParam' (mean, sd, (3rd not used) and optional seed (as 4th))
  ## 'impParam' .. (numeric) 1st for mean; 2nd for sd; 3rd for seed
  ## used (so far) in subsequent loops of testRobustToNAimputation
  isNa <- is.na(dat)
  if(length(impParam) >3) set.seed(as.integer(impParam[4]))
  impDat <- stats::rnorm(round(1.5*sum(isNa)), mean=impParam[1], sd=impParam[2])
  if(exclNeg) impDat <- impDat[which(impDat >0)]
  dat[which(isNa)] <- impDat[1:sum(isNa)]
  if(inclLowValMod) list(data=dat, datImp=impDat)
  dat }
       
