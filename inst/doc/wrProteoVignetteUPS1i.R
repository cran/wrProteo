## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment = "#>")

## ----install, echo=TRUE, eval=FALSE-------------------------------------------
#  # If not already installed, you'll have to install this package and wrMisc first.
#  install.packages("wrMisc")
#  install.packages("wrProteo")
#  
#  # The package wrGraph is recommended for better graphics
#  install.packages("wrGraph")
#  
#  # You cat start the vignettes for this package by typing :
#  browseVignettes("wrProteo")    #  ... and the select the html output

## ----setup, echo=TRUE, messages=FALSE, warnings=FALSE-------------------------
library(knitr)
library(wrMisc)
library(wrGraph)
library(wrProteo)

# Version number for wrProteo :
packageVersion("wrProteo")

## ----functions1, echo=TRUE----------------------------------------------------
## A few functions we'll need lateron
mergeVectors <- function(...,callFrom=NULL,silent=FALSE) {
  ## merge for simple named vectors (each element needs to be named)
  namesXYZ <- c(deparse(substitute(...)))
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="mergeVectors")
  inpL <- list(...)
  chNa <- sapply(inpL, function(x) length(unique(names(x)))==length(x))
   cat("yy\n"); yy <<- list(inpL=inpL,chNa=chNa,namesXYZ=namesXYZ)
  if(any(!chNa)) {if(!silent) message(fxNa," Vectors must have names on each element for merging; omit vectors ")
    inpL <- inpL }
  if(length(names(inpL)) <1) { names(inpL) <- 1:length(inpL)}
  if(length(inpL) >0) {
    spe <- sort(unique(unlist(lapply(inpL,names))))
    ta3 <- matrix(0, nrow=length(inpL), ncol=length(spe), dimnames=list(names(inpL),spe))
    for(i in 1:length(inpL)) ta3[i, match(names(inpL[[i]]),spe)] <- inpL[[i]]
    ta3 
  } else NULL }

replSpecType <- function(x, annCol="SpecType", replBy=cbind(old=c("mainSpe","species2"), new=c("Yeast","UPS1")), silent=TRUE) {
  ## rename $annot[,"SpecType"] to more specific names
  chCol <- annCol[1] %in% colnames(x$annot)
  if(chCol) { chCol <- which(colnames(x$annot)==annCol[1])
    chIt <- replBy[,1] %in% unique(x$annot[,chCol])    # check items to replace if present
    if(any(chIt)) for(i in which(chIt)) {useLi <- which(x$annot[,chCol] %in% replBy[i,1]); cat("useLi",head(useLi),"\n"); x$annot[useLi,chCol] <- replBy[i,2]}
  } else if(!silent) message(" replSpecType: 'annCol' not found in x$annot !")
  x }
  
replNAProtNames <- function(x,annCol=c("EntryName","Accession","SpecType"), silent=FALSE) {
  ## replace in $annot missing EntryNames by concatenating Accession + SpecType (ie 2nd & 3rd of annCol)
  chCol <- annCol %in% colnames(x$annot)
  if(all(chCol)) {
    chNA <- is.na(x$annot[,annCol[1]])
    if(any(chNA)) { if(!silent) message(" ..replacing ",sum(chNA)," entry-names")
      x$annot[which(chNA),annCol[1]] <- paste(x$annot[which(chNA),annCol[2]],x$annot[which(chNA),annCol[3]],sep="_")}
  } else message(" replNAProtNames: some of the columnnames from 'annCol' found in x$annot !")
  x }

standEntMatr <- function(mat,useCol=NULL) {
  ## standardize selected content of entire matrix (ie not column-wise, relative differences in row get thus conserved), return selected columns only 
  if(is.null(useCol)) useCol <- 1:ncol(mat)
  out <- as.matrix(mat[,useCol])
  (out - mean(out,na.rm=TRUE))/sd(out,na.rm=TRUE)
}

reorgByCluNo <- function(mat, cluNo, useMeth=1:3,useCol="sco", cluCol="cluNo", retList=FALSE, silent=FALSE, callFrom=NULL) {
  ## reorganize input matrix as sorted by cluster numbers (and geometric mean) according to vector with cluster names, and index for sorting per cluster and per geometric mean
  ## mat (matrix or data.frame) main input
  ## cluNo (integer) initial cluster numbers for each line of 'mat' (obtained by separate clustering or other segmentation)
  ## useMeth (character or integer) the columns to use from mat
  ## useCol (character or integer) the column to use for buidling geometrix mean
  ## retList (character or integer) the culumn of 'mat' which will be used for sorting the clusters
  ## retList (logical) decide if return list of clusters with data from 'cluNo' or matrix in order of input 'mat' with index,cluNo,geoMean
  fxNa <- wrMisc::.composeCallName(callFrom, newNa="reorgByCluNo")
  iniCla <- class(mat)
  chDim <- dim(mat)
  dataOK <- FALSE
  if(length(chDim) >1) if(all(chDim >1)) dataOK <- TRUE
  if(!silent & !dataOK) message(fxNa," invalid input ... (returning entry)")
  ## main
  if(dataOK){
    mat1 <- if(length(chDim) ==3) mat[,useMeth,useCol] else mat[,useMeth]
    nClu <- length(unique(wrMisc::naOmit(cluNo)))
    ## construct geometric mean for sorting
    mat1 <- cbind(mat1, index=1:nrow(mat), geoMean=apply(mat1,1, function(x) prod(x,na.rm=TRUE)^(1/sum(!is.na(x)))))    
    ## 1: split in list, determine clu median, overall score & sort clusters  
    cluL <- by(mat1, cluNo, as.matrix)
    clTab <- table(cluNo)                # already sorted by cluNo
    if(length(clTab) < max(cluNo) & !silent) message(fxNa," Note: Some cluster-names seem to be absent (no-consecutive numbers for names) !")
    ## need to correct when single occurance
    if(any(clTab ==1)) for(i in which(clTab ==1)) cluL[[i]] <- matrix(cluL[[i]],nrow=1,dimnames=list(rownames(mat1)[which(cluNo==i)] ,colnames(mat1)))
    ## sort intra
    cluL <- lapply(cluL, function(x) if(nrow(x) >1) x[order(x[,ncol(x)], decreasing=TRUE),] else x)
    ## sort inter
    cluL <- cluL[order(sapply(cluL, function(x) median(x[,ncol(x)],na.rm=TRUE)),  decreasing=TRUE)]        
    names(cluL) <- 1:length(cluL)
    nByClu  <- sapply(cluL,nrow)
    if(retList) { out <- cluL
      for(i in 1:nClu) out[[i]] <- cbind(out[[i]], cluNo=1)
       if("data.frame" %in% iniCla) out <- lapply(out,wrMisc::convMatr2df, addIniNa=FALSE, silent=silent,callFrom=fxNa) 
    } else {
      out <- cbind(wrMisc::lrbind(cluL), cluNo=rep(1:nClu,nByClu)) }   # in order of input
  } else {out <- NULL; if(!silent) message(fxNa," invalid input, return NULL")}
  out }


methNa <- c("ProteomeDiscoverer","MaxQuant","Proline" )

## The accession numbers for the UPS1 proteins
UPS1ac <- c("P00915", "P00918", "P01031", "P69905", "P68871", "P41159", "P02768", "P62988",
  "P04040", "P00167", "P01133", "P02144", "P15559", "P62937", "Q06830", "P63165",
  "P00709", "P06732", "P12081", "P61626", "Q15843", "P02753", "P16083", "P63279",
  "P01008", "P61769", "P55957", "O76070", "P08263", "P01344", "P01127", "P10599",
  "P99999", "P06396", "P09211", "P01112", "P01579", "P02787", "O00762", "P51965",
  "P08758", "P02741", "P05413", "P10145", "P02788", "P10636-8", "P00441", "P01375" ) 

## ----readMaxQuant, fig.height=8, fig.width=9.5, fig.align="center", echo=TRUE----
path1 <- system.file("extdata", package="wrProteo")
fiNaMQ <- "proteinGroups.txt.gz"
specPrefMQ <- list(conta="CON_|LYSC_CHICK", mainSpecies="OS=Saccharomyces cerevisiae", spike=UPS1ac)

dataMQ <- readMaxQuantFile(path1, file=fiNaMQ, specPref=specPrefMQ, refLi="mainSpe")

## ----readMaxQuant2, fig.height=8, fig.width=9.5, fig.align="center", echo=TRUE----
## the number of lines and colums
dim(dataMQ$quant)
## a summary of the quantitation data
summary(dataMQ$quant[,1:8])                # the first 8 cols

## ----readProteomeDiscoverer1, fig.height=8, fig.width=9.5, fig.align="center", echo=TRUE----
path1 <- system.file("extdata", package="wrProteo")
fiNaPd <- "pxd001819_PD24_Proteins.txt.gz"

## Note: data exported from ProteomeDiscoverer do not have proper column-names, providing names here
UPSconc <- c(50,125,250,500,2500,5000,12500,25000,50000)  
sampNa <- paste(rep(UPSconc, each=3),"amol_",rep(1:3,length(UPSconc)),sep="") 
specPrefPD <- list(conta="Bos tauris|Gallus", mainSpecies="OS=Saccharomyces cerevisiae", spike=UPS1ac)

dataPD <- readPDExport(file=fiNaPd, path=path1, sampleNames=sampNa, refLi="mainSpe", specPref=specPrefPD)

## ----readProteomeDiscoverer2, fig.height=8, fig.width=9.5, fig.align="center", echo=TRUE----
## the number of lines and colums
dim(dataPD$quant)
## a summary of the quantitation data
summary(dataPD$quant[,1:8])        # the first 8 cols

# there are proteins where the 'OS='-tag won't be visible as Species (orig Fasta-header and Protein-name not accessible) :
which(is.na(dataPD$annot[,"Species"]) & dataPD$annot[,"SpecType"]=="species2")    # is NA as Species

## ----readProline, fig.height=8, fig.width=9.5, fig.align="center", echo=TRUE----
## shifted for not printing
path1 <- system.file("extdata", package="wrProteo")
fiNaPl <- "pxd001819_PL.xlsx"

specPrefPL <- c(conta="_conta", mainSpecies="OS=Saccharomyces cerevisiae", spike="_ups")  
dataPL <- readProlineFile(file.path(path1,fiNaPl), specPref=specPrefPL, normalizeMeth="median", refLi="mainSpe")

## ----postTreatmPL, echo=TRUE--------------------------------------------------
head(colnames(dataPL$raw),8)
sub1 <- cbind(ini=paste0("-",c("100a","250a","500a","1f","5f","10f","25f","50f","100f"),"mol-"),
  paste0(fin=c("50a","125a","250a","500a","2500a","5000a","12500a","25000a","50000a"),"mol_"))
dataPL <- cleanListCoNames(dataPL, rem=c("abundance_","Levure2ug+ UPS1"), subst=sub1)

## ----readProlineInfo, fig.height=8, fig.width=9.5, fig.align="center", echo=TRUE----
## the number of lines and colums
dim(dataPL$quant)
## a summary of the quantitation data
summary(dataPL$quant[,1:8])        # the first 8 cols

## ----rearrange1, echo=TRUE----------------------------------------------------
# get all results (MaxQuant,ProteomeDiscoverer, ...) in same order
grp9 <- paste0(rep(UPSconc,each=3),"amol") 

## it is more convenient to re-order columns this way in each project
dataPD <- corColumnOrder(dataPD,sampNames=sampNa)          # already in good order
dataMQ <- corColumnOrder(dataMQ,sampNames=sampNa) 
dataPL <- corColumnOrder(dataPL,sampNames=sampNa) 

## ----postTreatm1, echo=TRUE---------------------------------------------------
## Need to rename $annot[,"SpecType"]  
dataPD <- replSpecType(dataPD, replBy=cbind(old=c("mainSpe","species2"), new=c("Yeast","UPS1")))
dataMQ <- replSpecType(dataMQ, replBy=cbind(old=c("mainSpe","species2"), new=c("Yeast","UPS1")))
dataPL <- replSpecType(dataPL, replBy=cbind(old=c("mainSpe","species2"), new=c("Yeast","UPS1")))

## Need to addres missing ProteinNames (UPS1) due to missing tags in Fasta
dataPD <- replNAProtNames(dataPD) 
dataMQ <- replNAProtNames(dataMQ) 
dataPL <- replNAProtNames(dataPL) 

## ----postTreatmCheck, echo=TRUE-----------------------------------------------
## extract names of quantified UPS1-proteins
NamesUpsPD <- dataPD$annot[which(dataPD$annot[,"SpecType"]=="UPS1"),"Accession"]
NamesUpsMQ <- dataMQ$annot[which(dataMQ$annot[,"SpecType"]=="UPS1"),"Accession"]
NamesUpsPL <- dataPL$annot[which(dataPL$annot[,"SpecType"]=="UPS1"),"Accession"]

## ----postTreatmTables, echo=TRUE----------------------------------------------
tabS <- mergeVectors(PD=table(dataPD$annot[,"SpecType"]), MQ=table(dataMQ$annot[,"SpecType"]), PL=table(dataPL$annot[,"SpecType"]))  
knitr::kable(tabS, caption="Number of proteins identified, by custom tags and software")
tabT <- mergeVectors(PD=table(dataPD$annot[,"Species"]), MQ=table(dataMQ$annot[,"Species"]), PL=table(dataPL$annot[,"Species"]))  
knitr::kable(tabT, caption="Number of proteins identified, by species and software")

## ----NA_ProteomeDiscoverer, echo=TRUE-----------------------------------------
## Let's inspect NA values as graphic
matrixNAinspect(dataPD$quant, gr=grp9, tit="ProteomeDiscoverer")  

## ----NA_MaxQuant, echo=TRUE---------------------------------------------------
## Let's inspect NA values as graphic
matrixNAinspect(dataMQ$quant, gr=gl(length(UPSconc),3), tit="MaxQuant") 

## ----NA_Proline, echo=TRUE----------------------------------------------------
## Let's inspect NA values as graphic
matrixNAinspect(dataPL$quant, gr=as.factor(substr(colnames(dataPL$quant),1,1)), tit="Proline") 

## ----nNA1, echo=TRUE----------------------------------------------------------
## Let's look at the number of NAs. Is there an accumulated number in lower UPS1 samples ?
tabSumNA <- rbind(PD=sumNAperGroup(dataPD$raw, grp9), MQ=sumNAperGroup(dataMQ$raw, grp9), PL=sumNAperGroup(dataPL$raw, grp9) )
knitr::kable(tabSumNA, caption="Number of NAs per group of samples", align="r")

## ----testProteomeDiscoverer, echo=TRUE----------------------------------------
testPD <- testRobustToNAimputation(dataPD, gr=grp9, lfdrInclude=TRUE)     # ProteomeDiscoverer

## ----testMaxQuant, echo=TRUE--------------------------------------------------
testMQ <- testRobustToNAimputation(dataMQ, gr=grp9, lfdrInclude=TRUE)      # MaxQuant

## ----testProline, echo=TRUE---------------------------------------------------
testPL <- testRobustToNAimputation(dataPL, gr=grp9, lfdrInclude=TRUE)      # Proline

## ----testReorganize1, echo=TRUE-----------------------------------------------
dataPD$datImp <- testPD$datImp       # recuperate imputeded data to main data-object
dataMQ$datImp <- testMQ$datImp
dataPL$datImp <- testPL$datImp

## ----PCA2, fig.height=12, fig.width=9.5, fig.align="center", echo=TRUE--------
# limit to UPS1 
plotPCAw(testPD$datImp[which(testPD$annot[,"SpecType"]=="UPS1"),], sampleGrp=grp9, tit="PCA on ProteomeDiscoverer, UPS1 only (NAs imputed)",rowTyName="proteins", useSymb2=0)
plotPCAw(testMQ$datImp[which(testMQ$annot[,"SpecType"]=="UPS1"),], sampleGrp=grp9, tit="PCA on MaxQuant, UPS1 only (NAs imputed)",rowTyName="proteins", useSymb2=0)
plotPCAw(testPL$datImp[which(testPL$annot[,"SpecType"]=="UPS1"),], sampleGrp=grp9, tit="PCA on Proline, UPS1 only (NAs imputed)",rowTyName="proteins", useSymb2=0)

## ----pairWise2, fig.height=4.5, fig.width=9.5, fig.align="center", echo=TRUE----
## The number of differentially abundant proteins passing 5% FDR (ProteomeDiscoverer and MaxQuant) 
signCount <- cbind( sig.PD.BH=colSums(testPD$BH < 0.05, na.rm=TRUE), sig.PD.lfdr=if("lfdr" %in% names(testPD)) colSums(testPD$lfdr < 0.05, na.rm=TRUE),
  sig.MQ.BH=colSums(testMQ$BH < 0.05, na.rm=TRUE), sig.MQ.lfdr=if("lfdr" %in% names(testMQ)) colSums(testMQ$lfdr < 0.05, na.rm=TRUE),
  sig.PL.BH=colSums(testPL$BH < 0.05, na.rm=TRUE), sig.PL.lfdr=if("lfdr" %in% names(testPL)) colSums(testPL$lfdr < 0.05, na.rm=TRUE)  )

table1 <- numPairDeColNames(testPD$BH, stripTxt="amol", sortByAbsRatio=TRUE)
table1 <- cbind(table1, signCount[table1[,1],])
knitr::kable(table1, caption="All pairwise comparisons and number of significant proteins", align="c")

## ----pairWise3, fig.height=4.5, fig.width=9.5, fig.align="center", echo=TRUE----
par(mar=c(6.2, 4.7, 4, 1))   
imageW(table1[,c("sig.PD.BH","sig.MQ.BH","sig.PL.BH" )], tit="Number of BH.FDR signif proteins by the quantification approaches")
mtext("red for high number signif proteins", cex=0.7)

## ----pairWiseSelect2, echo=TRUE-----------------------------------------------
## Selection in Ramus paper 
knitr::kable(table1[which(rownames(table1) %in% colnames(testPD$BH)[c(2,21,27)]),], caption="Selected pairwise comparisons (as in Ramus et al)", align="c")

## ----ROC_main1, echo=TRUE-----------------------------------------------------
layout(1)
rocPD <- lapply(table1[,1], function(x) summarizeForROC(testPD, annotCol="SpecType", spec=c("Yeast","UPS1"), columnTest=x, tyThr="BH", plotROC=FALSE))
rocMQ <- lapply(table1[,1], function(x) summarizeForROC(testMQ, annotCol="SpecType", spec=c("Yeast","UPS1"), columnTest=x, tyThr="BH", plotROC=FALSE))
rocPL <- lapply(table1[,1], function(x) summarizeForROC(testPL, annotCol="SpecType", spec=c("Yeast","UPS1"), columnTest=x, tyThr="BH", plotROC=FALSE))

names(rocPD) <- colnames(testPD$BH)
names(rocMQ) <- colnames(testMQ$BH)
names(rocPL) <- colnames(testPL$BH)

## calulate  AUC for each ROC 
AucAll <- cbind(ind=table1[match(names(rocPD),rownames(table1)),"index"], comb=NA, clu=NA, 
  PD=sapply(rocPD,AucROC), MQ=sapply(rocMQ,AucROC), PL=sapply(rocPL,AucROC) )

## ----ROC_segm, fig.height=9, fig.width=9.5, fig.align="center", echo=TRUE-----
## number of groups for clustering
nGr <- 5
## K-Means clustering
kMAx <- stats::kmeans(standEntMatr(AucAll[,c("PD","MQ","PL")]), nGr)$cluster  
   table(kMAx)
AucAll[,"clu"] <- kMAx

## ----ROC_segm2, echo=TRUE-----------------------------------------------------
AucAll <- reorgByCluNo(AucAll,cluNo=kMAx,useMeth=c("PD","MQ","PL"))
AucAll <- cbind(AucAll, iniInd=table1[match(rownames(AucAll), rownames(table1)), "index"])
colnames(AucAll)[1:(which(colnames(AucAll)=="index")-1)] <- paste("Auc",colnames(AucAll)[1:(which(colnames(AucAll)=="index")-1)], sep=".")

kMAx <- AucAll[,"cluNo"]   # update
  table(AucAll[,"cluNo"])
 ## note : column 'index' is relative to table1, iniInd to ordering inside objects from clustering 

## ----ROC_segm3, echo=TRUE-----------------------------------------------------
layout(1)
plot(AucAll[,"geoMean"], type="l", col=grey(0.7), lty=2, las=1, ylab="AUC",
  main="Pairwise Comparisons as Clustered AUC from ROC Curves", xlab="Comparison number")
abline(v=cumsum(table(AucAll[,"cluNo"])[-length(unique(AucAll[,"cluNo"]))])+0.5 ,lty=4,col=grey(0.8))      # clu borders
points(1:nrow(AucAll), AucAll[,"Auc.PD"], pch=AucAll[,"cluNo"], col=1)
points(1:nrow(AucAll), AucAll[,"Auc.MQ"], pch=AucAll[,"cluNo"], col=2)
points(1:nrow(AucAll), AucAll[,"Auc.PL"], pch=AucAll[,"cluNo"], col=3)
legend("bottomleft",c("PD","MQ","PL","geomMean"), text.col=c(1:3,1),pch=c(8,8,8,NA),lty=c(NA,NA,NA,2),lwd=c(NA,NA,NA,1.5),col=c(1:3,1),cex=0.85,xjust=0.5,yjust=0.5)

## ----ROC_segm4, echo=TRUE-----------------------------------------------------
AucRep <- table(AucAll[,"cluNo"])[rank(unique(AucAll[,"cluNo"]))]   # representative for each cluster
AucRep <- round(cumsum(AucRep) -AucRep/2 +0.1) 

## select representative for each cluster
knitr::kable(round(AucAll[AucRep,c("Auc.PD","Auc.MQ","Auc.PL","cluNo")],3), caption="Selected representative for each cluster ", align="c")
  

## ----ROC_segm5, fig.height=9, fig.width=9.5, fig.align="center", echo=TRUE----
## appendix : what characterizes a good or bad auc ?
biplot(prcomp(AucAll[,1:3]), cex=0.7, main="PCA of AUC from ROC Curves")   

## ----ROC_grp1, fig.height=10, fig.width=9.5, fig.align="center", echo=TRUE----
gr <- 1 
colPanel <- 2:5
layout(1)

j <- match(rownames(AucAll)[AucRep[gr]],names(rocPD)) 

## table of all of best cluster
useLi <- which(AucAll[,"cluNo"]==gr)
knitr::kable(cbind(round(AucAll[useLi,c("cluNo","Auc.PD","Auc.MQ","Auc.PL")],3), 
  table1[match(names(which(AucAll[,"cluNo"]==gr)),rownames(table1)),c(2,5,7,9)]), caption="AUC details for best pairwise-comparisons ", align="c")  
## frequent concentrations :
layout(matrix(1:2), heights=c(1,2.5)) 
tbl <- table(table1[match(names(which(AucAll[,"cluNo"]==gr)), rownames(table1)),c(3:4)])  # with(mydata, table(Species, Depth))
barplot(tbl, las=1, beside = TRUE, main=paste("Frequency of UPS-1 Concentrations Appearing in Cluster",gr))
    
## representative ROC    
plotROC(rocPD[[j]],rocMQ[[j]],rocPL[[j]], col=colPanel, methNames=methNa, pointSi=0.8, xlim=c(0,0.45),
  txtLoc=c(0.12,0.1,0.03), tit=paste("Cluster",gr," Example: ",rownames(AucAll)[AucRep[gr]]), legCex=1)

## ----VolcanoClu1, fig.height=10, fig.width=9.5, fig.align="center", echo=TRUE----
layout(matrix(1:4,ncol=2))
VolcanoPlotW(testPD, useComp=j, FCthrs=2, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[1], silent=TRUE)
VolcanoPlotW(testMQ, useComp=j, FCthrs=2, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[2], silent=TRUE)
VolcanoPlotW(testPL, useComp=j, FCthrs=2, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[3], silent=TRUE)

## ----ROC_grp2, fig.height=10, fig.width=9.5, fig.align="center", echo=TRUE----
gr <- 2
j <- match(rownames(AucAll)[AucRep[gr]],names(rocPD)) 

## table of all of best cluster
useLi <- which(AucAll[,"cluNo"]==gr)
knitr::kable(cbind(round(AucAll[useLi,c("cluNo","Auc.PD","Auc.MQ","Auc.PL")],3), 
  table1[match(names(which(AucAll[,"cluNo"]==gr)),rownames(table1)),c(2,5,7,9)]), caption="AUC details for 2nd best pairwise-comparisons ", align="c")  
## frequent concentrations :
layout(matrix(1:2), heights=c(1,3)) 
tbl <- table(table1[match(names(which(AucAll[,"cluNo"]==gr)),rownames(table1)),c(3:4)])  # with(mydata, table(Species, Depth))
barplot(tbl, las=1, beside = TRUE, main=paste("Frequency of UPS-1 Concentrations Appearing in Cluster",gr))
    
## repreentative ROC    
plotROC(rocPD[[j]],rocMQ[[j]],rocPL[[j]], col=colPanel, methNames=methNa, pointSi=0.8, xlim=c(0,0.45),
  txtLoc=c(0.12,0.1,0.03), tit=paste("Cluster",gr," Example: ",rownames(AucAll)[AucRep[gr]]), legCex=1)

## ----VolcanoClu2, fig.height=10, fig.width=9.5, fig.align="center", echo=TRUE----
layout(matrix(1:4, ncol=2)) 
VolcanoPlotW(testPD, useComp=j, FCthrs=2, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[1], silent=TRUE)
VolcanoPlotW(testMQ, useComp=j, FCthrs=2, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[2], silent=TRUE)
VolcanoPlotW(testPL, useComp=j, FCthrs=2, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[3], silent=TRUE)

## ----ROC_grp3, fig.height=10, fig.width=9.5, fig.align="center", echo=TRUE----
gr <- 3 
j <- match(rownames(AucAll)[AucRep[gr]],names(rocPD)) 

## table of all of best cluster
useLi <- which(AucAll[,"cluNo"]==gr)
knitr::kable(cbind(round(AucAll[useLi,c("cluNo","Auc.PD","Auc.MQ","Auc.PL")],3), 
  table1[match(names(which(AucAll[,"cluNo"]==gr)),rownames(table1)),c(2,5,7,9)]), caption="AUC details for 3rd best pairwise-comparisons ", align="c")  
## frequent concentrations :
layout(matrix(1:2), heights=c(1,3)) 
tbl <- table(table1[match(names(which(AucAll[,"cluNo"]==gr)),rownames(table1)),c(3:4)])  # with(mydata, table(Species, Depth))
barplot(tbl, las=1, beside = TRUE, main=paste("Frequency of UPS-1 Concentrations Appearing in Cluster",gr))
    
## representative ROC    
plotROC(rocPD[[j]],rocMQ[[j]],rocPL[[j]], col=colPanel, methNames=methNa, pointSi=0.8, xlim=c(0,0.45),
  txtLoc=c(0.12,0.1,0.03), tit=paste("Cluster",gr," Example: ",rownames(AucAll)[AucRep[gr]]), legCex=1)

## ----VolcanoClu3, fig.height=10, fig.width=9.5, fig.align="center", echo=TRUE----
layout(matrix(1:4, ncol=2)) 
VolcanoPlotW(testPD, useComp=j, FCthrs=2, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[1], silent=TRUE)
VolcanoPlotW(testMQ, useComp=j, FCthrs=2, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[2], silent=TRUE)
VolcanoPlotW(testPL, useComp=j, FCthrs=2, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[3], silent=TRUE)

## ----ROC_grp4, fig.height=10, fig.width=9.5, fig.align="center", echo=TRUE----
gr <- 4 
j <- match(rownames(AucAll)[AucRep[gr]],names(rocPD)) 

## table of all of best cluster
useLi <- which(AucAll[,"cluNo"]==gr)
knitr::kable(cbind(round(AucAll[useLi,c("cluNo","Auc.PD","Auc.MQ","Auc.PL")],3), 
  table1[match(names(which(AucAll[,"cluNo"]==gr)),rownames(table1)),c(2,5,7,9)]), caption="AUC details for 4th pairwise-comparisons ", align="c")  
## frequent concentrations :
layout(matrix(1:2), heights=c(1,3)) 
tbl <- table(table1[match(names(which(AucAll[,"cluNo"]==gr)),rownames(table1)),c(3:4)])  # with(mydata, table(Species, Depth))
barplot(tbl, las=1, beside = TRUE, main=paste("Frequency of UPS-1 Concentrations Appearing in Cluster",gr))
    
## representative ROC    
plotROC(rocPD[[j]],rocMQ[[j]],rocPL[[j]], col=colPanel, methNames=methNa, pointSi=0.8, xlim=c(0,0.45),
  txtLoc=c(0.12,0.1,0.03), tit=paste("Cluster",gr," Example: ",rownames(AucAll)[AucRep[gr]]), legCex=1)

## ----VolcanoClu4, fig.height=10, fig.width=9.5, fig.align="center", echo=TRUE----
layout(matrix(1:4, ncol=2)) 
VolcanoPlotW(testPD, useComp=j, FCthrs=2, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[1], silent=TRUE)
VolcanoPlotW(testMQ, useComp=j, FCthrs=2, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[2], silent=TRUE)
VolcanoPlotW(testPL, useComp=j, FCthrs=2, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[3], silent=TRUE)

## ----ROC_grp5, fig.height=10, fig.width=9.5, fig.align="center", echo=TRUE----
gr <- 5 
j <- match(rownames(AucAll)[AucRep[gr]],names(rocPD)) 

## table of all of best cluster
useLi <- which(AucAll[,"cluNo"]==gr)
knitr::kable(cbind(round(AucAll[useLi,c("cluNo","Auc.PD","Auc.MQ","Auc.PL")],3), 
  table1[match(names(which(AucAll[,"cluNo"]==gr)),rownames(table1)),c(2,5,7,9)]), caption="AUC details for 5th pairwise-comparisons ", align="c")  
## frequent concentrations :
layout(matrix(1:2), heights=c(1,3)) 
tbl <- table(table1[match(names(which(AucAll[,"cluNo"]==gr)),rownames(table1)),c(3:4)])  # with(mydata, table(Species, Depth))
barplot(tbl, las=1, beside = TRUE, main=paste("Frequency of UPS-1 Concentrations Appearing in Cluster",gr))
    
## representative ROC    
plotROC(rocPD[[j]],rocMQ[[j]],rocPL[[j]], col=colPanel, methNames=methNa, pointSi=0.8, xlim=c(0,0.45),
  txtLoc=c(0.12,0.1,0.03), tit=paste("Cluster",gr," Example: ",rownames(AucAll)[AucRep[gr]]), legCex=1)

## ----VolcanoClu5, fig.height=10, fig.width=9.5, fig.align="center", echo=TRUE----
layout(matrix(1:4, ncol=2)) 
VolcanoPlotW(testPD, useComp=j, FCthrs=2, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[1], silent=TRUE)
VolcanoPlotW(testMQ, useComp=j, FCthrs=2, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[2], silent=TRUE)
VolcanoPlotW(testPL, useComp=j, FCthrs=2, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[3], silent=TRUE) 

## ----nNA2, echo=TRUE----------------------------------------------------------
tab1 <- rbind(PD=sumNAperGroup(dataPD$raw[which(dataPD$annot[,"SpecType"]=="UPS1"),], grp9),
  MQ=sumNAperGroup(dataMQ$raw[which(dataMQ$annot[,"SpecType"]=="UPS1"),], grp9),
  PL= sumNAperGroup(dataPL$raw[which(dataPL$annot[,"SpecType"]=="UPS1"),], grp9)  ) 
knitr::kable(tab1, caption="The number of NAs in the UPS1 proteins", align="c")

## ----nNAfig1, fig.height=3.5, fig.width=9.5, fig.align="center", echo=TRUE----
countRawNA <- function(dat, newOrd=UPS1ac, relative=FALSE) {  # count number of NAs per UPS protein and order as UPS
  out <- rowSums(is.na(dat$raw[match(newOrd,rownames(dat$raw)),])) 
  if(relative) out/nrow(dat$raw) else out }

sumNAperMeth <- cbind(PD=countRawNA(dataPD), MQ=countRawNA(dataMQ), PL=countRawNA(dataPL) )
UPS1na <- sub("_UPS","",dataPL$annot[UPS1ac,"EntryName"])
par(mar=c(6.8, 3.5, 4, 1))   
imageW(sumNAperMeth, rowNa=UPS1na, tit="Number of NAs in UPS proteins")
mtext("red for high number of NAs",cex=0.7)

## ----intraReplicCV1, fig.height=10, fig.width=12, fig.align="center", echo=TRUE----
## combined plot : all data (left), Ups1 (right)
layout(1:3)
sumNAinPD <- list(length=18)
sumNAinPD[2*(1:length(unique(grp9))) -1] <- as.list(as.data.frame(log2(rowGrpCV(testPD$datImp, grp9))))
sumNAinPD[2*(1:length(unique(grp9))) ] <- as.list(as.data.frame(log2(rowGrpCV(testPD$datImp[which(testPD$annot[,"SpecType"]=="UPS1"),], grp9))))
names(sumNAinPD)[2*(1:length(unique(grp9))) -1] <-  sub("amol","",unique(grp9))
names(sumNAinPD)[2*(1:length(unique(grp9))) ] <- paste(sub("amol","",unique(grp9)),"Ups",sep=".")
vioplotW(sumNAinPD, halfViolin="pairwise", tit="CV Intra Replicate, ProteomeDiscoverer", cexNameSer=0.6) 
mtext("left part : all data\nright part: UPS1",adj=0,cex=0.8)

sumNAinMQ <- list(length=18)
sumNAinMQ[2*(1:length(unique(grp9))) -1] <- as.list(as.data.frame(log2(rowGrpCV(testMQ$datImp, grp9))))
sumNAinMQ[2*(1:length(unique(grp9))) ] <- as.list(as.data.frame(log2(rowGrpCV(testMQ$datImp[which(testMQ$annot[,"SpecType"]=="UPS1"),], grp9))))
names(sumNAinMQ)[2*(1:length(unique(grp9))) -1] <- sub("amol","",unique(grp9))                        # paste(unique(grp9),"all",sep=".")
names(sumNAinMQ)[2*(1:length(unique(grp9))) ] <- paste(sub("amol","",unique(grp9)),"Ups",sep=".")      #paste(unique(grp9),"Ups1",sep=".")
vioplotW(sumNAinMQ, halfViolin="pairwise", tit="CV intra replicate, MaxQuant",cexNameSer=0.6) 
mtext("left part : all data\nright part: UPS1",adj=0,cex=0.8)

sumNAinPL <- list(length=18)
sumNAinPL[2*(1:length(unique(grp9))) -1] <- as.list(as.data.frame(log2(rowGrpCV(testPL$datImp, grp9))))
sumNAinPL[2*(1:length(unique(grp9))) ] <- as.list(as.data.frame(log2(rowGrpCV(testPL$datImp[which(testPL$annot[,"SpecType"]=="UPS1"),], grp9))))
names(sumNAinPL)[2*(1:length(unique(grp9))) -1] <-  sub("amol","",unique(grp9))
names(sumNAinPL)[2*(1:length(unique(grp9))) ] <- paste(sub("amol","",unique(grp9)),"Ups",sep=".")
vioplotW(sumNAinPL, halfViolin="pairwise", tit="CV Intra Replicate, Proline", cexNameSer=0.6) 
mtext("left part : all data\nright part: UPS1",adj=0,cex=0.8)

## ----linModel0, echo=TRUE-----------------------------------------------------
## prepare object for storing all results
datUPS1 <- array(NA, dim=c(length(UPS1ac),length(methNa),7), dimnames=list(UPS1ac,c("PD","MQ","PL"),
  c("sco","nPep","medAbund", "logp","slope","startFr","cluNo")))

## ----linModelPD, fig.height=17, fig.width=9.5, fig.align="center", echo=TRUE----
lmPD <- list(length=length(NamesUpsPD))
doPl <- FALSE
lmPD[1:length(NamesUpsPD)] <- lapply(NamesUpsPD[1:length(NamesUpsPD)], linModelSelect, dat=dataPD, 
  expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=doPl, silent=TRUE)
names(lmPD) <- NamesUpsPD

## ----linModelPD2, echo=TRUE---------------------------------------------------
## We make a little summary of regression-results (ProteomeDiscoverer)
linIn <- match(names(lmPD), UPS1ac)
datUPS1[linIn,1,c("logp","slope","startFr")] <- cbind(log10(sapply(lmPD, function(x) x$coef[2,4])), 
  sapply(lmPD, function(x) x$coef[2,1]), sapply(lmPD, function(x) x$startLev) )
## need correct rownames in  dataPD$datImp !!! 
datUPS1[,1,"medAbund"] <- apply(wrMisc::.scale01(dataPD$datImp)[match(UPS1ac,rownames(dataPD$datImp)),],1,median,na.rm=TRUE)

## ----linModelMQ, fig.height=17, fig.width=9.5, fig.align="center", echo=TRUE----
lmMQ <- list(length=length(NamesUpsMQ))
lmMQ[1:length(NamesUpsMQ)] <- lapply(NamesUpsMQ[1:length(NamesUpsMQ)], linModelSelect, dat=dataMQ, 
  expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=doPl, silent=TRUE)
names(lmMQ) <- NamesUpsMQ

## ----linModelMQ2, fig.height=17, fig.width=9.5, fig.align="center", echo=TRUE----
## We make a little summary of regression-results (MaxQuant)
linIn <- match(names(lmMQ), UPS1ac)
datUPS1[linIn,2,c("logp","slope","startFr")] <- cbind( log10(sapply(lmMQ, function(x) x$coef[2,4])), 
  sapply(lmMQ, function(x) x$coef[2,1]), sapply(lmMQ, function(x) x$startLev) )
datUPS1[,2,"medAbund"] <- apply(wrMisc::.scale01(dataMQ$datImp)[match(UPS1ac,rownames(dataMQ$datImp)),],1,median,na.rm=TRUE)

## ----linModelPL, fig.height=17, fig.width=9.5, fig.align="center", echo=TRUE----
lmPL <- list(length=length(NamesUpsPL))
lmPL[1:length(NamesUpsPL)] <- lapply(NamesUpsPL[1:length(NamesUpsPL)], linModelSelect, dat=dataPL, 
  expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=doPl, silent=TRUE)
names(lmPL) <- NamesUpsPL   

## ----linModelPLsum, fig.height=17, fig.width=9.5, fig.align="center", echo=TRUE----
## We make a little summary of regression-results (Proline including Percolator)
linIn <- match(names(lmPL), UPS1ac)
datUPS1[linIn,3,c("logp","slope","startFr")] <- cbind(log10(sapply(lmPL, function(x) x$coef[2,4])), 
  sapply(lmPL, function(x) x$coef[2,1]), sapply(lmPL, function(x) x$startLev) )
datUPS1[,3,"medAbund"] <- apply(wrMisc::.scale01(dataPL$datImp)[match(UPS1ac,rownames(dataPL$datImp)),],1,median,na.rm=TRUE)

## ----linModelStartStat,  echo=TRUE--------------------------------------------
## at which concentration of UPS1 did the best regression start ?
stTab <- sapply(1:5, function(x) apply(datUPS1[,,"startFr"],2,function(y) sum(x==y)))
colnames(stTab) <- paste("lev",1:5,sep="_")
knitr::kable(stTab, caption = "Frequency of starting levels for regression")

## ----linModelPlotAll, fig.height=12, fig.width=9.5, fig.align="center", echo=TRUE----
layout(matrix(1:4,ncol=2))
subTi <- "fill according to median abundance (violet=low - green - red=high)"

tit <- "ProteomeDiscoverer UPS1, p-value vs slope"
useCol <- colorAccording2(datUPS1[,1,"medAbund"], gradTy="rainbow", revCol=TRUE, nEndOmit=14)
plot(datUPS1[,1,c("logp","slope")], main=tit, type="n",xlim=c(-25,-1),ylim=c(0.1,1.6))   #col=1, bg.col=useCol, pch=20+lmPDsum[,"startFr"],
points(datUPS1[,1,c("logp","slope")], col=1, bg=useCol, pch=20+datUPS1[,1,"startFr"],)
legend("topright",paste("best starting from ",1:5), text.col=1, pch=21:25, col=1, pt.bg="white", cex=0.9, xjust=0.5, yjust=0.5)
mtext(subTi,cex=0.9)
  abline(v=c(-12,-10),lty=2,col="grey") ; abline(h=c(0.7,0.75),lty=2,col="grey")
hi1 <- hist(datUPS1[,1,"medAbund"], plot=FALSE)
legendHist(sort(datUPS1[,1,"medAbund"]), colRamp=useCol[order(datUPS1[,1,"medAbund"])][cumsum(hi1$counts)], location="bottomleft", legTit="median raw abundance")  #
##
##
tit <- "MaxQuant UPS1, p-value vs slope"
useCol <- colorAccording2(datUPS1[,2,"medAbund"], gradTy="rainbow", revCol=TRUE, nEndOmit=14)
plot(datUPS1[,2,c("logp","slope")], main=tit, type="n",xlim=c(-25,-1),ylim=c(0.1,1.6))   
points(datUPS1[,2,c("logp","slope")], col=1, bg=useCol, pch=20+datUPS1[,2,"startFr"],)
legend("topright",paste("best starting from ",1:5), text.col=1, pch=21:25, col=1, pt.bg="white", cex=0.9, xjust=0.5, yjust=0.5)
mtext(subTi,cex=0.9)
  abline(v=c(-12,-10),lty=2,col="grey") ; abline(h=c(0.7,0.75),lty=2,col="grey")
hi1 <- hist(datUPS1[,2,"medAbund"], plot=FALSE)
legendHist(sort(datUPS1[,2,"medAbund"]), colRamp=useCol[order(datUPS1[,2,"medAbund"])][cumsum(hi1$counts)], location="bottomleft", legTit="median raw abundance")  #
##
##
tit <- "Proline UPS1, p-value vs slope"
useCol <- colorAccording2(datUPS1[,3,"medAbund"], gradTy="rainbow", revCol=TRUE, nEndOmit=14)
plot(datUPS1[,3,c("logp","slope")], main=tit, type="n",xlim=c(-25,-1),ylim=c(0.1,1.6))   #
points(datUPS1[,3,c("logp","slope")], col=1, bg=useCol, pch=20+datUPS1[,3,"startFr"],)
legend("topright",paste("best starting from ",1:5), text.col=1, pch=21:25, col=1, pt.bg="white", cex=0.9, xjust=0.5, yjust=0.5)
mtext(subTi,cex=0.9)
  abline(v=c(-12,-10),lty=2,col="grey") ; abline(h=c(0.7,0.75),lty=2,col="grey")
hi1 <- hist(datUPS1[,3,"medAbund"], plot=FALSE)
legendHist(sort(datUPS1[,3,"medAbund"]), colRamp=useCol[order(datUPS1[,3,"medAbund"])][cumsum(hi1$counts)], location="bottomleft", legTit="median raw abundance")  #

## ----combRegrScore1, echo=TRUE------------------------------------------------
for(i in 1:(dim(datUPS1)[2])) datUPS1[,i,"sco"] <- -datUPS1[,i,"logp"] - (datUPS1[,i,"slope"] -1)^2    # cut at > 8

## ----combRegrScore2, echo=TRUE------------------------------------------------
datUPS1[,1,2] <- rowSums(dataPD$count[match(UPS1ac,dataPD$annot[,1]),,1], na.rm=TRUE)
datUPS1[,2,2] <- rowSums(dataMQ$count[match(UPS1ac,dataMQ$annot[,1]),,2], na.rm=TRUE)
datUPS1[,3,2] <- rowSums(dataPL$count[match(UPS1ac,dataPL$annot[,1]),,2], na.rm=TRUE)

## ----combRegrScore3, fig.height=6, fig.width=9.5, fig.align="center", echo=TRUE----
layout(matrix(1:4,ncol=2))
par(mar=c(5.5, 3, 4, 0.4))  
imageW(datUPS1[,,1],col=heat.colors(15), tit="Linear regression score")
mtext("red for bad score", cex=0.8)

imageW(log(datUPS1[,,2]),col=rev(heat.colors(15)),tit="number of peptides")
mtext("red for high number of peptides", cex=0.8)

## ratio : regression score vs no of peptides
imageW(datUPS1[,,1]/log(datUPS1[,,2]),col=rev(heat.colors(15)),tit="Regression score / Number of peptides")
mtext("red for high (good) lmScore/peptide ratio)", cex=0.8)

## score vs abundance
imageW(datUPS1[,,1]/datUPS1[,,3], col=rev(heat.colors(15)),tit="Regression score / median Abundance")
mtext("red for high (good) lmScore/abundance ratio)", cex=0.8)

## ----combScore1, echo=TRUE----------------------------------------------------
## number of groups for clustering
nGr <- 5

## clustering using kMeans
kMx <- stats::kmeans(standEntMatr(datUPS1[,1:3,"sco"] ,useCol=c("PD","MQ","PL")), nGr)$cluster  
datUPS1[,,"cluNo"] <- matrix(rep(kMx,dim(datUPS1)[2]), nrow=length(kMx))

datUPS1clu <- reorgByCluNo(datUPS1,cluNo=kMx,useMeth=1:3)
## bring datUPS1 in order of clustering
datUPS1 <- datUPS1[match(rownames(datUPS1clu),rownames(datUPS1)),,]
datUPS1[,,"cluNo"] <- rep(datUPS1clu[,"cluNo"], ncol(datUPS1))                                # as renamed clusters

## graphcal inspection of clustering of UPS proteins based on regression 
## this graph could be improved
    cluNum <- datUPS1[,1,"cluNo"]
    plot(datUPS1[,"MQ","sco"], type="p",pch=cluNum,col=2,ylab="Regression score",xlab="UPS1 proteins",main="Clustered Regression Score for UPS1 Proteins",las=1)
    abline(v=cumsum(table(cluNum)[-length(unique(cluNum))])+0.5 ,lty=2,col=grey(0.8))      # clu borders
    points(1:nrow(datUPS1),datUPS1[,"PD","sco"], pch=cluNum, col=1)
    points(1:nrow(datUPS1),datUPS1[,"PL","sco"], pch=cluNum, col=3)

    lines(1:nrow(datUPS1),datUPS1clu[,"geoMean"],col=grey(0.5),lty=2,lwd=1.6)   # mean
    legend("bottomleft",c("PD","MQ","PL","geomMean"), text.col=c(1:3,1),pch=c(8,8,8,NA),lty=c(NA,NA,NA,2),lwd=c(NA,NA,NA,1.5),col=c(1:3,1),cex=0.85,xjust=0.5,yjust=0.5)        # 1st as point, 2nd as text

## ----combScore3, echo=TRUE----------------------------------------------------
## representative for each cluster  (median position inside cluster)
UPSrep <- table(datUPS1[,1,"cluNo"])[rank(unique(datUPS1[,1,"cluNo"]))]
UPSrep <- round(cumsum(UPSrep) -UPSrep/2 +0.1) 

## prepare annotation of UPS proteins
annUPS1 <- dataPL$annot[match(rownames(datUPS1),dataPL$annot[,1]), c(1,3)]
annUPS1[,2] <- substr(sub("_UPS","",sub("generic_ups\\|[[:alnum:]]+-{0,1}[[:digit:]]\\|","",annUPS1[,2])),1,42)

## ----combRegrClu1Tab, echo=TRUE-----------------------------------------------
gr <- 1
useLi <- which(datUPS1[,1,"cluNo"]==gr)
colNa <- c("Protein",paste(colnames(datUPS1),rep(c("slope","logp"),each=ncol(datUPS1)),sep=" "))
knitr::kable(cbind(annUPS1[useLi,2], signif(datUPS1[useLi,,"slope"],3), signif(datUPS1[useLi,,"logp"],3)), 
  caption="Regression details for cluster of best UPS1 proteins ", col.names=colNa, align="l")

## ----regrPlot5star, fig.height=9, fig.width=9.5, fig.align="center", echo=TRUE----
## Plotting the best regressions
layout(matrix(1:4, ncol=2))
tit <- paste0(methNa,", ",annUPS1[UPSrep[gr],1])
tm <- linModelSelect(annUPS1[UPSrep[gr],1],dat=dataPD, tit=tit[1], expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE)
tm <- linModelSelect(annUPS1[UPSrep[gr],1],dat=dataMQ, tit=tit[2],expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE)
tm <- linModelSelect(annUPS1[UPSrep[gr],1],dat=dataPL, tit=tit[3], expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE)

## ----combRegrClu2Tab, echo=TRUE-----------------------------------------------
gr <- 2
useLi <- which(datUPS1[,1,"cluNo"]==gr)
knitr::kable(cbind(annUPS1[useLi,2], signif(datUPS1[useLi,,"slope"],3), signif(datUPS1[useLi,,"logp"],3)), 
  caption="Regression details for cluster of 2nd best UPS1 proteins ", col.names=colNa, align="l")

## ----regrPlot4star, fig.height=9, fig.width=9.5, fig.align="center", echo=TRUE----
layout(matrix(1:4,ncol=2))
tit <- paste0(methNa,", ",annUPS1[UPSrep[gr],1])
tm <- linModelSelect(annUPS1[UPSrep[gr],1],dat=dataPD, tit=tit[1], expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE)
tm <- linModelSelect(annUPS1[UPSrep[gr],1],dat=dataMQ, tit=tit[2],expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE)
tm <- linModelSelect(annUPS1[UPSrep[gr],1],dat=dataPL, tit=tit[3], expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE)

## ----combRegrClu3Tab, echo=TRUE-----------------------------------------------
gr <- 3
useLi <- which(datUPS1[,1,"cluNo"]==gr)
knitr::kable(cbind(annUPS1[useLi,2], signif(datUPS1[useLi,,"slope"],3), signif(datUPS1[useLi,,"logp"],3)), 
  caption="Regression details for 3rd cluster UPS1 proteins ", col.names=colNa, align="l")

## ----regrPlot3star, fig.height=9, fig.width=9.5, fig.align="center", echo=TRUE----
layout(matrix(1:4, ncol=2))
tit <- paste0(methNa,", ",annUPS1[UPSrep[gr],1])
tm <- linModelSelect(annUPS1[UPSrep[gr],1],dat=dataPD, tit=tit[1], expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE)
tm <- linModelSelect(annUPS1[UPSrep[gr],1],dat=dataMQ, tit=tit[2],expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE)
tm <- linModelSelect(annUPS1[UPSrep[gr],1],dat=dataPL, tit=tit[3], expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE)

## ----regrPlot2star, fig.height=9, fig.width=9.5, fig.align="center", echo=TRUE----
gr <- 4
useLi <- which(datUPS1[,1,"cluNo"]==gr)
tit <- paste0(methNa,", ",annUPS1[UPSrep[gr],1])
layout(matrix(1:4, ncol=2))
tm <- linModelSelect(annUPS1[UPSrep[gr],1],dat=dataPD, tit=tit[1], expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE)
tm <- linModelSelect(annUPS1[UPSrep[gr],1],dat=dataMQ, tit=tit[2],expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE)
tm <- linModelSelect(annUPS1[UPSrep[gr],1],dat=dataPL, tit=tit[3], expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE)

## ----regrPlot1star, fig.height=9, fig.width=9.5, fig.align="center", echo=TRUE----
gr <- 5
useLi <- which(datUPS1[,1,"cluNo"]==gr)
knitr::kable(cbind(annUPS1[useLi,2], signif(datUPS1[useLi,,"slope"],3), signif(datUPS1[useLi,,"logp"],3)), 
  caption="Regression details for 5th cluster UPS1 proteins ", col.names=colNa, align="l")
tit <- paste0(methNa,", ",annUPS1[UPSrep[gr],1])
layout(matrix(1:4, ncol=2))
tm <- linModelSelect(annUPS1[UPSrep[gr],1],dat=dataPD, tit=tit[1], expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE)
tm <- linModelSelect(annUPS1[UPSrep[gr],1],dat=dataMQ, tit=tit[2],expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE)
tm <- linModelSelect(annUPS1[UPSrep[gr],1],dat=dataPL, tit=tit[3], expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

