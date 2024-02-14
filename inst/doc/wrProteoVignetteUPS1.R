## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment = "#>")

## ----install, echo=TRUE, eval=FALSE-------------------------------------------
#  ## This is R code, you can run this to redo all analysis presented here.
#  ## If not already installed, you'll have to install wrMisc and wrProteo first.
#  install.packages("wrMisc")
#  install.packages("wrProteo")
#  ## These packages are used for the graphics
#  install.packages("wrGraph")
#  install.packages("RColorBrewer")
#  
#  ## Installation of limma from Bioconductor
#  if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
#  BiocManager::install("limma")
#  
#  ## You cat also see all vignettes for this package by typing :
#  browseVignettes("wrProteo")    #  ... and the select the html output

## ----setup, echo=TRUE, messages=FALSE, warnings=FALSE-------------------------
## Let's assume this is a fresh R-session
library(knitr)
library(wrMisc)
library(wrGraph)
library(wrProteo)

# Version number for wrProteo :
packageVersion("wrProteo")

## ----metaData1, echo=TRUE-----------------------------------------------------
## Read meta-data from  github.com/bigbio/proteomics-metadata-standard/
pxd001819meta <- readSdrf("PXD001819")

## The concentration of the UPS1 spike-in proteins in the samples
if(length(pxd001819meta) >0) {
  UPSconc <- sort(unique(as.numeric(wrMisc::trimRedundText(pxd001819meta$characteristics.spiked.compound.))))  # trim to get to 'essential' info
} else {
  UPSconc <- c(50, 125, 250, 500, 2500, 5000, 12500, 25000, 50000)       # in case access to github failed
}

## ----functions1, echo=TRUE----------------------------------------------------
## A few elements and functions we'll need lateron
methNa <- c("ProteomeDiscoverer","MaxQuant","Proline")
names(methNa) <- c("PD","MQ","PL")

## The accession numbers for the UPS1 proteins
UPS1 <- data.frame(ac=c("P00915", "P00918", "P01031", "P69905", "P68871", "P41159", "P02768", "P62988",
  "P04040", "P00167", "P01133", "P02144", "P15559", "P62937", "Q06830", "P63165",
  "P00709", "P06732", "P12081", "P61626", "Q15843", "P02753", "P16083", "P63279",
  "P01008", "P61769", "P55957", "O76070", "P08263", "P01344", "P01127", "P10599",
  "P99999", "P06396", "P09211", "P01112", "P01579", "P02787", "O00762", "P51965",
  "P08758", "P02741", "P05413", "P10145", "P02788", "P10636-8", "P00441", "P01375"),
  species=rep("Homo sapiens", 48),
  name=NA)

## ----functions2, echo=TRUE----------------------------------------------------
## additional functions
replSpecType <- function(x, annCol="SpecType", replBy=cbind(old=c("mainSpe","species2"), new=c("Yeast","UPS1")), silent=TRUE) {
  ## rename $annot[,"SpecType"] to more specific names
  chCol <- annCol[1] %in% colnames(x$annot)
  if(chCol) { chCol <- which(colnames(x$annot)==annCol[1])
    chIt <- replBy[,1] %in% unique(x$annot[,chCol])    # check items to replace if present
    if(any(chIt)) for(i in which(chIt)) {useLi <- which(x$annot[,chCol] %in% replBy[i,1]); cat("useLi",head(useLi),"\n"); x$annot[useLi,chCol] <- replBy[i,2]}
  } else if(!silent) message(" replSpecType: 'annCol' not found in x$annot !")
  x }

plotConcHist <- function(mat, ref, refColumn=3:4, matCluNa="cluNo", lev=NULL, ylab=NULL, tit=NULL) {
  ## plot histogram like counts of UPS1 concentrations
  if(is.null(tit)) tit <- "Frequency of UPS1 Concentrations Appearing in Cluster"
  gr <- unique(mat[,matCluNa])
  ref <- ref[,refColumn]
  if(length(lev) <2) lev <- sort(unique(as.numeric(as.matrix(ref))))
  if(length(ylab) !=1) ylab <- "Frequency"
  tbl <- table(factor( as.numeric(ref[which(rownames(ref) %in% rownames(mat)),]), levels=lev))
  graphics::barplot(tbl, las=1, beside=TRUE, main=paste(tit,gr), col=grDevices::gray(0.8), ylab=ylab)
}

plotMultRegrPar <- function(dat, methInd, tit=NULL, useColumn=c("logp","slope","medAbund","startFr"), lineGuide=list(v=c(-12,-10),h=c(0.7,0.75),col="grey"), xlim=NULL,ylim=NULL,subTit=NULL) {
  ## scatter plot logp (x) vs slope (y) for all UPS proteins, symbol by useColumn[4], color by hist of useColumn[3]
  ## dat (array) UPS1 data
  ## useColumn (character) 1st as 'logp', 2nd as 'slope', 3rd as median abundance, 4th as starting best regression from this point
  fxNa <- "plotMultRegrPar"
   #fxNa <- wrMisc::.composeCallName(callFrom,newNa="plotMultRegrPar")
  if(length(dim(dat)) !=3) stop("invalid input, expecting as 'dat' array with 3 dimensions (proteins,Softw,regrPar)")
  if(any(length(methInd) >1, methInd > dim(dat)[2], !is.numeric(methInd))) stop("invalid 'methInd'")
  chCol <- useColumn %in% dimnames(dat)[[3]]
  if(any(!chCol)) stop("argument 'useColumn' does not fit to 3rd dim dimnames of 'dat'")
  useCol <- colorAccording2(dat[,methInd,useColumn[3]], gradTy="rainbow", revCol=TRUE, nEndOmit=14)
  graphics::plot(dat[,methInd,useColumn[1:2]], main=tit, type="n",xlim=xlim,ylim=ylim)   #col=1, bg.col=useCol, pch=20+lmPDsum[,"startFr"],
  graphics::points(dat[,methInd,useColumn[1:2]], col=1, bg=useCol, pch=20+dat[,methInd,useColumn[4]],)
  graphics::legend("topright",paste("best starting from ",1:5), text.col=1, pch=21:25, col=1, pt.bg="white", cex=0.9, xjust=0.5, yjust=0.5)
  if(length(subTit)==1) graphics::mtext(subTit,cex=0.9)
  if(is.list(lineGuide) & length(lineGuide) >0) {if(length(lineGuide$v) >0) graphics::abline(v=lineGuide$v,lty=2,col=lineGuide$col)
    if(length(lineGuide$h) >0) graphics::abline(h=lineGuide$h,lty=2,col=lineGuide$col)}
  hi1 <- graphics::hist(dat[,methInd,useColumn[3]], plot=FALSE)
  wrGraph::legendHist(sort(dat[,methInd,useColumn[3]]), colRamp=useCol[order(dat[,methInd,useColumn[3]])][cumsum(hi1$counts)],
    cex=0.5, location="bottomleft", legTit="median raw abundance")  #
}

## ----readMaxQuant, fig.height=8, fig.width=9.5, fig.align="center", echo=TRUE----
path1 <- system.file("extdata", package="wrProteo")
fiNaMQ <- "proteinGroups.txt.gz"

## We need to define the setup of species
specPrefMQ <- list(conta="CON_|LYSC_CHICK", mainSpecies="OS=Saccharomyces cerevisiae", spike=UPS1$ac)
dataMQ <- readMaxQuantFile(path1, file=fiNaMQ, specPref=specPrefMQ, refLi="mainSpe",
  sdrf=c("PXD001819","max"), suplAnnotFile=TRUE, plotGraph=FALSE)

## ----readMaxQuant2, fig.height=8, fig.width=9.5, fig.align="center", echo=TRUE----
## The number of lines and colums
dim(dataMQ$quant)
## A quick summary of some columns of quantitation data
summary(dataMQ$quant[,1:7])                # the first 8 cols
table(dataMQ$annot[,"SpecType"], useNA="always")

## ----readProteomeDiscoverer1, fig.height=8, fig.width=9.5, fig.align="center", echo=TRUE----
path1 <- system.file("extdata", package="wrProteo")
fiNaPd <- "pxd001819_PD24_Proteins.txt.gz"
## Next, we define the setup of species
specPrefPD <- list(conta="Bos tauris|Gallus", mainSpecies="Saccharomyces cerevisiae", spike=UPS1$ac)
dataPD <- readProteomeDiscovererFile(file=fiNaPd, path=path1, refLi="mainSpe", specPref=specPrefPD,
  sdrf=c("PXD001819","max"), plotGraph=FALSE)

## ----readProteomeDiscoverer2, fig.height=8, fig.width=9.5, fig.align="center", echo=TRUE----
## The number of lines and colums
dim(dataPD$quant)
## A quick summary of some columns of quantitation data
summary(dataPD$quant[,1:7])        # the first 8 cols
table(dataPD$annot[,"SpecType"], useNA="always")

## ----readProline, fig.height=8, fig.width=9.5, fig.align="center", echo=TRUE----
path1 <- system.file("extdata", package="wrProteo")
fiNaPl <- "pxd001819_PL.xlsx"

specPrefPL <- list(conta="_conta", mainSpecies="Saccharomyces cerevisiae", spike=UPS1$ac)
dataPL <- readProlineFile(fiNaPl, path=path1, specPref=specPrefPL, normalizeMeth="median", refLi="mainSpe",
  sdrf=c("PXD001819","max"), plotGraph=FALSE)

## ----postTreatmPL, echo=TRUE--------------------------------------------------
head(colnames(dataPL$raw), 7)
dataPL <- cleanListCoNames(dataPL, rem=c("Levure2ug+ UPS1-"), subst=cbind(c("fmol","mol-"), c("000amol","mol_R")), mathOper="/2")

## let's check the result
head(colnames(dataPL$raw),8)

## ----readProlineInfo, fig.height=8, fig.width=9.5, fig.align="center", echo=TRUE----
## The number of lines and colums
dim(dataPL$quant)
## A quick summary of some columns of quantitation data
summary(dataPL$quant[,1:8])        # the first 8 cols
table(dataPL$annot[,"SpecType"], useNA="always")

## ----rearrange1, echo=TRUE----------------------------------------------------
## bring all results (MaxQuant,ProteomeDiscoverer, ...) in same ascending order
## as reference will use the order from ProteomeDiscoverer, it's output is already in a convenient order
sampNa <- colnames(dataPD$quant)

## it is more convenient to re-order columns this way in each project
dataPD <- corColumnOrder(dataPD, sampNames=sampNa)          # already in good order
dataMQ <- corColumnOrder(dataMQ, replNames=paste0("UPS1_",sub("amol_", "amol_R", colnames(dataMQ$quant))), sampNames=sampNa)  # incl canged names
dataPL <- corColumnOrder(dataPL, replNames=paste0("UPS1_",colnames(dataPL$quant)), sampNames=sampNa)               # incl canged names

## ----postTreatm1, echo=TRUE---------------------------------------------------
## Need to rename $annot[,"SpecType"]
dataPD <- replSpecType(dataPD, replBy=cbind(old=c("mainSpe","species2"), new=c("Yeast","UPS1")))
dataMQ <- replSpecType(dataMQ, replBy=cbind(old=c("mainSpe","species2"), new=c("Yeast","UPS1")))
dataPL <- replSpecType(dataPL, replBy=cbind(old=c("mainSpe","species2"), new=c("Yeast","UPS1")))

## Need to address missing ProteinNames (UPS1) due to missing tags in Fasta
dataPD <- replMissingProtNames(dataPD)
dataMQ <- replMissingProtNames(dataMQ)
dataPL <- replMissingProtNames(dataPL)
    table(dataPD$annot[,"SpecType"])

## synchronize order of groups
(grp9 <- dataMQ$sampleSetup$level)
names(grp9) <- rep(paste0(UPSconc,"amol"), each=3)
dataPL$sampleSetup$groups <- dataMQ$sampleSetup$groups <- dataPD$sampleSetup$groups <- grp9  # synchronize order of groups

## ----postTreatmCheck, echo=TRUE-----------------------------------------------
## extract names of quantified UPS1-proteins
NamesUpsPD <- dataPD$annot[which(dataPD$annot[,"SpecType"]=="spike"), "Accession"]
NamesUpsMQ <- dataMQ$annot[which(dataMQ$annot[,"SpecType"]=="spike"), "Accession"]
NamesUpsPL <- dataPL$annot[which(dataPL$annot[,"SpecType"]=="spike"), "Accession"]

## ----postTreatmTables, echo=TRUE----------------------------------------------
tabS <- mergeVectors(PD=table(dataPD$annot[,"SpecType"]), MQ=table(dataMQ$annot[,"SpecType"]), PL=table(dataPL$annot[,"SpecType"]))
tabT <- mergeVectors(PD=table(dataPD$annot[,"Species"]), MQ=table(dataMQ$annot[,"Species"]), PL=table(dataPL$annot[,"Species"]))
tabS[which(is.na(tabS))] <- 0
tabT[which(is.na(tabT))] <- 0
kable(cbind(tabS[,2:1], tabT), caption="Number of proteins identified, by custom tags, species and software")

## ----metaData2, echo=TRUE-----------------------------------------------------
kable(cbind(dataMQ$sampleSetup$sdrf[,c(23,7,19,22)], groups=dataMQ$sampleSetup$groups))

## ----NA_ProteomeDiscoverer, echo=TRUE-----------------------------------------
## Let's inspect NA values from ProteomeDiscoverer as graphic
matrixNAinspect(dataPD$quant, gr=grp9, tit="ProteomeDiscoverer")

## ----NA_MaxQuant, echo=TRUE---------------------------------------------------
## Let's inspect NA values from MaxQuant as graphic
matrixNAinspect(dataMQ$quant, gr=grp9, tit="MaxQuant")

## ----NA_Proline, echo=TRUE----------------------------------------------------
## Let's inspect NA values from Proline as graphic
matrixNAinspect(dataPL$quant, gr=grp9, tit="Proline")

## ----nNA1, echo=TRUE----------------------------------------------------------
## Let's look at the number of NAs. Is there an accumulated number in lower UPS1 samples ?
tabSumNA <- rbind(PD=sumNAperGroup(dataPD$raw, grp9), MQ=sumNAperGroup(dataMQ$raw, grp9), PL=sumNAperGroup(dataPL$raw, grp9) )
kable(tabSumNA, caption="Number of NAs per group of samples", align="r")

## ----testProteomeDiscoverer, echo=TRUE----------------------------------------
testPD <- testRobustToNAimputation(dataPD, imputMethod="informed")     # ProteomeDiscoverer

## ----testMaxQuant, echo=TRUE--------------------------------------------------
testMQ <- testRobustToNAimputation(dataMQ, imputMethod="informed")      # MaxQuant , ok

## ----testProline, echo=TRUE---------------------------------------------------
testPL <- testRobustToNAimputation(dataPL, imputMethod="informed")      # Proline

## ----testReorganize1, echo=TRUE-----------------------------------------------
dataPD$datImp <- testPD$datImp       # recuperate imputeded data to main data-object
dataMQ$datImp <- testMQ$datImp
dataPL$datImp <- testPL$datImp

## ----pairWise2, echo=TRUE-----------------------------------------------------
## The number of differentially abundant proteins passing 5% FDR (ProteomeDiscoverer and MaxQuant)
signCount <- cbind( sig.PD.BH=colSums(testPD$BH < 0.05, na.rm=TRUE), sig.PD.lfdr=if("lfdr" %in% names(testPD)) colSums(testPD$lfdr < 0.05, na.rm=TRUE),
  sig.MQ.BH=colSums(testMQ$BH < 0.05, na.rm=TRUE), sig.MQ.lfdr=if("lfdr" %in% names(testMQ)) colSums(testMQ$lfdr < 0.05, na.rm=TRUE),
  sig.PL.BH=colSums(testPL$BH < 0.05, na.rm=TRUE), sig.PL.lfdr=if("lfdr" %in% names(testPL)) colSums(testPL$lfdr < 0.05, na.rm=TRUE)  )

table1 <- numPairDeColNames(testPD$BH, stripTxt="amol", sortByAbsRatio=TRUE)
table1 <- cbind(table1, signCount[table1[,1],])
rownames(table1) <- colnames(testMQ$BH)[table1[,1]]

kable(table1, caption="All pairwise comparisons and number of significant proteins", align="c")

## ----check2, echo=TRUE--------------------------------------------------------
resMQ1 <- extractTestingResults(testMQ, compNo=1, thrsh=0.05, FCthrs=2)
resPD1 <- extractTestingResults(testPD, compNo=1, thrsh=0.05, FCthrs=2)
resPL1 <- extractTestingResults(testPL, compNo=1, thrsh=0.05, FCthrs=2)

## ----pairWise3, fig.height=4.5, fig.width=9.5, fig.align="center", echo=TRUE----
par(mar=c(5.5, 4.7, 4, 1))
imageW(table1[,c("sig.PD.BH","sig.MQ.BH","sig.PL.BH" )], col=rev(RColorBrewer::brewer.pal(9,"YlOrRd")),
  transp=FALSE, tit="Number of BH.FDR passing proteins by the quantification approaches")
mtext("Dark red for high number signif proteins", cex=0.75)

## ----pairWiseSelect2, echo=TRUE-----------------------------------------------
## Selection in Ramus paper
kable(table1[which(rownames(table1) %in% colnames(testPD$BH)[c(2,21,27)]),], caption="Selected pairwise comparisons (as in Ramus et al)", align="c")

## ----ROC_main1, echo=TRUE-----------------------------------------------------
## calulate  AUC for each ROC
layout(1)
rocPD <- lapply(table1[,1], function(x) summarizeForROC(testPD, useComp=x, annotCol="SpecType", spec=c("mainSpecies","spike"), tyThr="BH", plotROC=FALSE,silent=TRUE))
rocMQ <- lapply(table1[,1], function(x) summarizeForROC(testMQ, useComp=x, annotCol="SpecType", spec=c("mainSpecies","spike"), tyThr="BH", plotROC=FALSE,silent=TRUE))
rocPL <- lapply(table1[,1], function(x) summarizeForROC(testPL, useComp=x, annotCol="SpecType", spec=c("mainSpecies","spike"), tyThr="BH", plotROC=FALSE,silent=TRUE))

# we still need to add the names for the pair-wise groups:
names(rocPD) <- names(rocMQ) <- names(rocPL) <- rownames(table1)

## ----ROC_main2, echo=TRUE-----------------------------------------------------
AucAll <- cbind(ind=table1[match(names(rocPD), rownames(table1)),"index"], clu=NA,
  PD=sapply(rocPD, AucROC), MQ=sapply(rocMQ, AucROC), PL=sapply(rocPL, AucROC) )

## ----ROC_biplot, fig.height=9, fig.width=9.5, fig.align="center", echo=TRUE----
try(biplot(prcomp(AucAll[,names(methNa)]), cex=0.7, main="PCA of AUC from ROC Curves"))

## ----ROC_segm, fig.height=9, fig.width=9.5, fig.align="center", echo=TRUE-----
## number of groups for clustering
nGr <- 5
## K-Means clustering
kMAx <- stats::kmeans(standardW(AucAll[,c("PD","MQ","PL")]), nGr)$cluster
   table(kMAx)
AucAll[,"clu"] <- kMAx

## ----ROC_segm2, echo=TRUE-----------------------------------------------------
AucAll <- reorgByCluNo(AucAll, cluNo=kMAx, useColumn=c("PD","MQ","PL"))
AucAll <- cbind(AucAll, iniInd=table1[match(rownames(AucAll), rownames(table1)), "index"])
colnames(AucAll)[1:(which(colnames(AucAll)=="index")-1)] <- paste("Auc",colnames(AucAll)[1:(which(colnames(AucAll)=="index")-1)], sep=".")
AucAll[,"cluNo"] <- rep(nGr:1, table(AucAll[,"cluNo"]))        # make cluNo descending

kMAx <- AucAll[,"cluNo"]      # update
  table(AucAll[,"cluNo"])
 ## note : column 'index' is relative to table1, iniInd to ordering inside objects from clustering

## ----ROC_profFig, echo=TRUE---------------------------------------------------
try(profileAsClu(AucAll[,c(1:length(methNa),(length(methNa)+2:3))], clu="cluNo", meanD="geoMean", tit="Pairwise Comparisons as Clustered AUC from ROC Curves",
  xlab="Comparison number", ylab="AUC", meLty=1, meLwd=3))

## ----ROC_segmTable, echo=TRUE-------------------------------------------------
AucRep <- table(AucAll[,"cluNo"])[rank(unique(AucAll[,"cluNo"]))]   # representative for each cluster
AucRep <- round(cumsum(AucRep) -AucRep/2 +0.1)

## select representative for each cluster
kable(round(AucAll[AucRep,c("Auc.PD","Auc.MQ","Auc.PL","cluNo")],3), caption="Selected representative for each cluster ", align="c")

## ----freqOfFCperClu, echo=TRUE------------------------------------------------
ratTab <- sapply(5:1, function(x) { y <- table1[match(rownames(AucAll),rownames(table1)),]
  table(factor(signif(y[which(AucAll[,"cluNo"]==x),"log2rat"],1), levels=unique(signif(table1[,"log2rat"],1))) )})
colnames(ratTab) <- paste0("\nclu",5:1,"\nn=",rev(table(kMAx)))
layout(1)
imageW(ratTab, tit="Frequency of rounded log2FC in the 5 clusters", xLab="log2FC (rounded)", col=RColorBrewer::brewer.pal(9,"YlOrRd"),las=1)
mtext("Dark red for enrichment of given pair-wise ratio", cex=0.7)

## ----ROC_grp5tab, echo=TRUE---------------------------------------------------
colPanel <- 2:5
gr <- 5
j <- match(rownames(AucAll)[AucRep[6-gr]], colnames(testPD$t))

## table of all proteins in cluster
useLi <- which(AucAll[,"cluNo"]==gr)
tmp <- cbind(round(as.data.frame(AucAll)[useLi,c("cluNo","Auc.PD","Auc.MQ","Auc.PL")],3),
  as.data.frame(table1)[match(names(useLi),rownames(table1)), c(2,5,7,9)])
kable(tmp, caption="AUC details for best pairwise-comparisons ", align="c")

## ----ROC_grp5fig, fig.height=9, fig.width=9.5, fig.align="center", echo=TRUE----
## frequent concentrations :
layout(matrix(1:2), heights=c(1,2.5))
plotConcHist(mat=tmp, ref=table1)

## representative ROC
jR <- match(rownames(AucAll)[AucRep[6-gr]], names(rocPD))
plotROC(rocPD[[jR]], rocMQ[[jR]], rocPL[[jR]], col=colPanel, methNames=methNa, pointSi=0.8, xlim=c(0,0.45),
  txtLoc=c(0.12,0.1,0.033), tit=paste("Cluster",gr," Example: ",names(rocPD)[jR]), legCex=1)

## ----VolcanoClu5, fig.height=10, fig.width=9.5, fig.align="center", echo=TRUE----
## This required package 'wrGraph' at version 1.2.5 (or higher)
if(packageVersion("wrGraph")  >= "1.2.5") {
  layout(matrix(1:4,ncol=2))
  try(VolcanoPlotW(testPD, useComp=j, FCthrs=1.5, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[1], expFCarrow=TRUE, silent=TRUE),silent=TRUE)
  try(VolcanoPlotW(testMQ, useComp=j, FCthrs=1.5, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[2], expFCarrow=TRUE, silent=TRUE),silent=TRUE)
  try(VolcanoPlotW(testPL, useComp=j, FCthrs=1.5, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[3], expFCarrow=TRUE, silent=TRUE),silent=TRUE)}

## ----ROC_grp4tab, echo=TRUE---------------------------------------------------
gr <- 4
j <- match(rownames(AucAll)[AucRep[6-gr]], colnames(testPD$t))

## table of all proteins in cluster
useLi <- which(AucAll[,"cluNo"]==gr)
tmp <- cbind(round(as.data.frame(AucAll)[useLi,c("cluNo","Auc.PD","Auc.MQ","Auc.PL")],3),
  as.data.frame(table1)[match(names(useLi),rownames(table1)), c(2,5,7,9)])
kable(tmp, caption="AUC details for cluster '++++' pairwise-comparisons ", align="c")

## ----ROC_grp4fig, fig.height=10, fig.width=9.5, fig.align="center", echo=TRUE----
## frequent concentrations :
layout(matrix(1:2), heights=c(1,2.5))
plotConcHist(mat=tmp, ref=table1)

## representative ROC
jR <- match(rownames(AucAll)[AucRep[6-gr]], names(rocPD))
plotROC(rocPD[[jR]], rocMQ[[jR]], rocPL[[jR]], col=colPanel, methNames=methNa, pointSi=0.8, xlim=c(0,0.45),
  txtLoc=c(0.12,0.1,0.033), tit=paste("Cluster",gr," Example: ",names(rocPD)[jR]), legCex=1)

## ----VolcanoClu4, fig.height=10, fig.width=9.5, fig.align="center", echo=TRUE----
if(packageVersion("wrGraph")  >= "1.2.5"){
  layout(matrix(1:4,ncol=2))
  try(VolcanoPlotW(testPD, useComp=j, FCthrs=1.5, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[1], expFCarrow=TRUE, silent=TRUE),silent=TRUE)
  try(VolcanoPlotW(testMQ, useComp=j, FCthrs=1.5, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[2], expFCarrow=TRUE, silent=TRUE),silent=TRUE)
  try(VolcanoPlotW(testPL, useComp=j, FCthrs=1.5, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[3], expFCarrow=TRUE, silent=TRUE),silent=TRUE)}

## ----ROC_grp3tab, echo=TRUE---------------------------------------------------
gr <- 3
j <- match(rownames(AucAll)[AucRep[6-gr]], colnames(testPD$t))

## table of all proteins in cluster
useLi <- which(AucAll[,"cluNo"]==gr)
tmp <- cbind(round(as.data.frame(AucAll)[useLi,c("cluNo","Auc.PD","Auc.MQ","Auc.PL")],3),
  as.data.frame(table1)[match(names(useLi),rownames(table1)), c(2,5,7,9)])
kable(tmp, caption="AUC details for cluster '+++' pairwise-comparisons ", align="c")

## ----ROC_grp3fig, fig.height=10, fig.width=9.5, fig.align="center", echo=TRUE----
## frequent concentrations :
layout(matrix(1:2), heights=c(1,2.5))
plotConcHist(mat=tmp, ref=table1)

## representative ROC
jR <- match(rownames(AucAll)[AucRep[6-gr]], names(rocPD))
plotROC(rocPD[[jR]],rocMQ[[jR]],rocPL[[jR]], col=colPanel, methNames=methNa, pointSi=0.8, xlim=c(0,0.45),
  txtLoc=c(0.12,0.1,0.033), tit=paste("Cluster",gr," Example: ",names(rocPD)[jR]), legCex=1)

## ----VolcanoClu3, fig.height=10, fig.width=9.5, fig.align="center", echo=TRUE----
if(packageVersion("wrGraph")  >= "1.2.5"){
  layout(matrix(1:4,ncol=2))
  try(VolcanoPlotW(testPD, useComp=j, FCthrs=1.5, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[1], expFCarrow=TRUE, silent=TRUE),silent=TRUE)
  try(VolcanoPlotW(testMQ, useComp=j, FCthrs=1.5, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[2], expFCarrow=TRUE, silent=TRUE),silent=TRUE)
  try(VolcanoPlotW(testPL, useComp=j, FCthrs=1.5, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[3], expFCarrow=TRUE, silent=TRUE),silent=TRUE)}

## ----ROC_grp2tab, echo=TRUE---------------------------------------------------
gr <- 2
j <- match(rownames(AucAll)[AucRep[6-gr]], colnames(testPD$t))

## table of all proteins in cluster
useLi <- which(AucAll[,"cluNo"]==gr)
tmp <- cbind(round(as.data.frame(AucAll)[useLi,c("cluNo","Auc.PD","Auc.MQ","Auc.PL")],3),
  as.data.frame(table1)[match(names(useLi),rownames(table1)), c(2,5,7,9)])
kable(tmp, caption="AUC details for cluster '++' pairwise-comparisons ", align="c")

## ----ROC_grp2fig, fig.height=10, fig.width=9.5, fig.align="center", echo=TRUE----
## frequent concentrations :
layout(matrix(1:2), heights=c(1,2.5))
plotConcHist(mat=tmp, ref=table1)

## representative ROC
jR <- match(rownames(AucAll)[AucRep[6-gr]], names(rocPD))
plotROC(rocPD[[jR]], rocMQ[[jR]], rocPL[[jR]], col=colPanel, methNames=methNa, pointSi=0.8, xlim=c(0,0.45),
  txtLoc=c(0.12,0.1,0.033), tit=paste("Cluster",gr," Example: ",names(rocPD)[jR]), legCex=1)

## ----VolcanoClu2, fig.height=10, fig.width=9.5, fig.align="center", echo=TRUE----
if(packageVersion("wrGraph")  >= "1.2.5"){
  layout(matrix(1:4,ncol=2))
  try(VolcanoPlotW(testPD, useComp=j, FCthrs=1.5, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[1], expFCarrow=TRUE, silent=TRUE),silent=TRUE)
  try(VolcanoPlotW(testMQ, useComp=j, FCthrs=1.5, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[2], expFCarrow=TRUE, silent=TRUE),silent=TRUE)
  try(VolcanoPlotW(testPL, useComp=j, FCthrs=1.5, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[3], expFCarrow=TRUE, silent=TRUE),silent=TRUE)}

## ----ROC_grp1tab, echo=TRUE---------------------------------------------------
gr <- 1
j <- match(rownames(AucAll)[AucRep[6-gr]], colnames(testPD$t))

## table of all proteins in cluster
useLi <- which(AucAll[,"cluNo"]==gr)
tmp <- cbind(round(as.data.frame(AucAll)[useLi,c("cluNo","Auc.PD","Auc.MQ","Auc.PL")],3),
  as.data.frame(table1)[match(names(useLi),rownames(table1)), c(2,5,7,9)])
kable(tmp, caption="AUC details for cluster '+' pairwise-comparisons ", align="c")

## ----ROC_grp1fig, fig.height=10, fig.width=9.5, fig.align="center", echo=TRUE----
## frequent concentrations :
layout(matrix(1:2, ncol=1), heights=c(1,2.5))
plotConcHist(mat=tmp, ref=table1)

## representative ROC
jR <- match(rownames(AucAll)[AucRep[6-gr]], names(rocPD))
plotROC(rocPD[[jR]], rocMQ[[jR]], rocPL[[jR]], col=colPanel, methNames=methNa, pointSi=0.8, xlim=c(0,0.45),
  txtLoc=c(0.12,0.1,0.033), tit=paste("Cluster",gr," Example: ",names(rocPD)[jR]), legCex=1)

## ----VolcanoClu1, fig.height=10, fig.width=9.5, fig.align="center", echo=TRUE----
if(packageVersion("wrGraph")  >= "1.2.5"){
  layout(matrix(1:4,ncol=2))
  try(VolcanoPlotW(testPD, useComp=j, FCthrs=1.5, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[1], expFCarrow=TRUE, silent=TRUE),silent=TRUE)
  try(VolcanoPlotW(testMQ, useComp=j, FCthrs=1.5, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[2], expFCarrow=TRUE, silent=TRUE),silent=TRUE)
  try(VolcanoPlotW(testPL, useComp=j, FCthrs=1.5, FdrThrs=0.05, annColor=c(4,2,3), ProjNa=methNa[3], expFCarrow=TRUE, silent=TRUE),silent=TRUE)}

## ----nNA2, echo=TRUE----------------------------------------------------------
tab1 <- rbind(PD=sumNAperGroup(dataPD$raw[which(dataPD$annot[,"SpecType"]=="UPS1"),], grp9),
  MQ=sumNAperGroup(dataMQ$raw[which(dataMQ$annot[,"SpecType"]=="UPS1"),], grp9),
  PL= sumNAperGroup(dataPL$raw[which(dataPL$annot[,"SpecType"]=="UPS1"),], grp9)  )
kable(tab1, caption="The number of NAs in the UPS1 proteins", align="c")

## ----nNAfig1, fig.height=3.5, fig.width=9.5, fig.align="center", echo=TRUE----
countRawNA <- function(dat, newOrd=UPS1$ac, relative=FALSE) {  # count number of NAs per UPS protein and order as UPS
  out <- rowSums(is.na(dat$raw[match(newOrd,rownames(dat$raw)),]))
  if(relative) out/nrow(dat$raw) else out }

sumNAperMeth <- cbind(PD=countRawNA(dataPD), MQ=countRawNA(dataMQ), PL=countRawNA(dataPL) )
UPS1na <- sub("_UPS","",dataPL$annot[UPS1$ac,"EntryName"])
par(mar=c(6.8, 3.5, 4, 1))
imageW(sumNAperMeth, rowNa=UPS1na, tit="Number of NAs in UPS proteins", xLab="", yLab="",
  transp=FALSE, col=rev(RColorBrewer::brewer.pal(9,"YlOrRd")))
mtext("Dark red for high number of NAs",cex=0.7)

## ----PCA2PD, fig.height=12, fig.width=9.5, fig.align="center", echo=TRUE------
try(plotPCAw(testPD$datImp[which(testPD$annot[,"SpecType"]=="spike"),], sampleGrp=grp9, tit="PCA on ProteomeDiscoverer, UPS1 only (NAs imputed)", rowTyName="proteins", useSymb2=0, silent=TRUE), silent=TRUE)

## ----PCA2MQ, fig.height=12, fig.width=9.5, fig.align="center", echo=TRUE------
try(plotPCAw(testMQ$datImp[which(testMQ$annot[,"SpecType"]=="spike"),], sampleGrp=grp9, tit="PCA on MaxQuant, UPS1 only (NAs imputed)", rowTyName="proteins", useSymb2=0, silent=TRUE), silent=TRUE)

## ----PCA2PL, fig.height=12, fig.width=9.5, fig.align="center", echo=TRUE------
try(plotPCAw(testPL$datImp[which(testPL$annot[,"SpecType"]=="spike"),], sampleGrp=grp9, tit="PCA on Proline, UPS1 only (NAs imputed)", rowTyName="proteins", useSymb2=0, silent=TRUE), silent=TRUE)

## ----intraReplicCV1, fig.height=10, fig.width=12, fig.align="center", echo=TRUE----
## combined plot : all data (left), Ups1 (right)
layout(1:3)
sumNAinPD <- list(length=18)
sumNAinPD[2*(1:length(unique(grp9))) -1] <- as.list(as.data.frame(log2(rowGrpCV(testPD$datImp, grp9))))
sumNAinPD[2*(1:length(unique(grp9))) ] <- as.list(as.data.frame(log2(rowGrpCV(testPD$datImp[which(testPD$annot[,"SpecType"]=="spike"),], grp9))))
names(sumNAinPD)[2*(1:length(unique(grp9))) -1] <-  sub("amol","",unique(grp9))
names(sumNAinPD)[2*(1:length(unique(grp9))) ] <- paste(sub("amol","",unique(grp9)),"Ups",sep=".")
try(vioplotW(sumNAinPD, halfViolin="pairwise", tit="CV Intra Replicate, ProteomeDiscoverer", cexNameSer=0.6))
mtext("left part : all data\nright part: UPS1",adj=0,cex=0.8)

sumNAinMQ <- list(length=18)
sumNAinMQ[2*(1:length(unique(grp9))) -1] <- as.list(as.data.frame(log2(rowGrpCV(testMQ$datImp, grp9))))
sumNAinMQ[2*(1:length(unique(grp9))) ] <- as.list(as.data.frame(log2(rowGrpCV(testMQ$datImp[which(testMQ$annot[,"SpecType"]=="spike"),], grp9))))
names(sumNAinMQ)[2*(1:length(unique(grp9))) -1] <- sub("amol","",unique(grp9))                        # paste(unique(grp9),"all",sep=".")
names(sumNAinMQ)[2*(1:length(unique(grp9))) ] <- paste(sub("amol","",unique(grp9)),"Ups",sep=".")      #paste(unique(grp9),"Ups1",sep=".")
try(vioplotW(sumNAinMQ, halfViolin="pairwise", tit="CV intra replicate, MaxQuant",cexNameSer=0.6))
mtext("left part : all data\nright part: UPS1",adj=0,cex=0.8)

sumNAinPL <- list(length=18)
sumNAinPL[2*(1:length(unique(grp9))) -1] <- as.list(as.data.frame(log2(rowGrpCV(testPL$datImp, grp9))))
sumNAinPL[2*(1:length(unique(grp9))) ] <- as.list(as.data.frame(log2(rowGrpCV(testPL$datImp[which(testPL$annot[,"SpecType"]=="spike"),], grp9))))
names(sumNAinPL)[2*(1:length(unique(grp9))) -1] <-  sub("amol","",unique(grp9))
names(sumNAinPL)[2*(1:length(unique(grp9))) ] <- paste(sub("amol","",unique(grp9)),"Ups",sep=".")
try(vioplotW(sumNAinPL, halfViolin="pairwise", tit="CV Intra Replicate, Proline", cexNameSer=0.6))
mtext("left part : all data\nright part: UPS1",adj=0,cex=0.8)

## ----linModel0, echo=TRUE-----------------------------------------------------
## prepare object for storing all results
datUPS1 <- array(NA, dim=c(length(UPS1$ac),length(methNa),7), dimnames=list(UPS1$ac,c("PD","MQ","PL"),
  c("sco","nPep","medAbund", "logp","slope","startFr","cluNo")))

## ----linModelPD, fig.height=17, fig.width=9.5, fig.align="center", echo=TRUE----
lmPD <- list(length=length(NamesUpsPD))
doPl <- FALSE
lmPD[1:length(NamesUpsPD)] <- lapply(NamesUpsPD[1:length(NamesUpsPD)], linModelSelect, dat=dataPD,
  expect=names(grp9), startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=doPl, silent=TRUE)
names(lmPD) <- NamesUpsPD

## ----linModelPD2, echo=TRUE---------------------------------------------------
## We make a little summary of regression-results (ProteomeDiscoverer)
tmp <- cbind(log10(sapply(lmPD, function(x) x$coef[2,4])), sapply(lmPD, function(x) x$coef[2,1]), sapply(lmPD, function(x) x$startLev))
datUPS1[,1,c("logp","slope","startFr")] <- tmp[match(rownames(datUPS1), names(lmPD)), ]
datUPS1[,1,"medAbund"] <- apply(wrMisc::.scale01(dataPD$datImp)[match(UPS1$ac,rownames(dataPD$datImp)),],1,median,na.rm=TRUE)

## ----linModelMQ, echo=TRUE----------------------------------------------------
lmMQ <- list(length=length(NamesUpsMQ))
lmMQ[1:length(NamesUpsMQ)] <- lapply(NamesUpsMQ[1:length(NamesUpsMQ)], linModelSelect, dat=dataMQ,
  expect=names(grp9), startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=doPl, silent=TRUE)
names(lmMQ) <- NamesUpsMQ

## ----linModelMQ2, fig.height=17, fig.width=9.5, fig.align="center", echo=TRUE----
## We make a little summary of regression-results (MaxQuant)
tmp <- cbind(log10(sapply(lmMQ, function(x) x$coef[2,4])), sapply(lmMQ, function(x) x$coef[2,1]), sapply(lmMQ, function(x) x$startLev))
datUPS1[,2,c("logp","slope","startFr")] <- tmp[match(rownames(datUPS1), names(lmMQ)), ]
datUPS1[,2,"medAbund"] <- apply(wrMisc::.scale01(dataMQ$datImp)[match(UPS1$ac,rownames(dataMQ$datImp)),],1,median,na.rm=TRUE)

## ----linModelPL, echo=TRUE----------------------------------------------------
lmPL <- list(length=length(NamesUpsPL))
lmPL[1:length(NamesUpsPL)] <- lapply(NamesUpsPL[1:length(NamesUpsPL)], linModelSelect, dat=dataPL,
  expect=names(grp9), startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=doPl, silent=TRUE)
names(lmPL) <- NamesUpsPL

## ----linModelPLsum, fig.height=17, fig.width=9.5, fig.align="center", echo=TRUE----
tmp <- cbind(log10(sapply(lmPL, function(x) x$coef[2,4])), sapply(lmPL, function(x) x$coef[2,1]), sapply(lmPL, function(x) x$startLev))
datUPS1[,3,c("logp","slope","startFr")] <- tmp[match(rownames(datUPS1), names(lmPL)), ]
datUPS1[,3,"medAbund"] <- apply(wrMisc::.scale01(dataPL$datImp)[match(UPS1$ac,rownames(dataPL$datImp)),],1,median,na.rm=TRUE)

## ----linModelStartStat,  echo=TRUE--------------------------------------------
## at which concentration of UPS1 did the best regression start ?
stTab <- sapply(1:5, function(x) apply(datUPS1[,,"startFr"],2,function(y) sum(x==y)))
colnames(stTab) <- paste("lev",1:5,sep="_")
kable(stTab, caption = "Frequency of starting levels for regression")

## ----linModelPlotAll, fig.height=12, fig.width=9.5, fig.align="center", echo=TRUE----
layout(matrix(1:4,ncol=2))
subTi <- "fill according to median abundance (blue=low - green - red=high)"
xyRa <- apply(datUPS1[,,4:5], 3, range, na.rm=T)

plotMultRegrPar(datUPS1, 1, xlim=xyRa[,1], ylim=xyRa[,2],tit="ProteomeDiscoverer UPS1, p-value vs slope",subTit=subTi)    # adj wr 9jan23
plotMultRegrPar(datUPS1, 2, xlim=xyRa[,1], ylim=xyRa[,2],tit="MaxQuant UPS1, p-value vs slope",subTit=subTi)
plotMultRegrPar(datUPS1, 3, xlim=xyRa[,1], ylim=xyRa[,2],tit="Proline UPS1, p-value vs slope",subTit=subTi)

## ----combRegrScore1, echo=TRUE------------------------------------------------
for(i in 1:(dim(datUPS1)[2])) datUPS1[,i,"sco"] <- -datUPS1[,i,"logp"] - (datUPS1[,i,"slope"] -1)^2    # cut at > 8

## ----combRegrScore2, echo=TRUE------------------------------------------------
datUPS1[,1,2] <- rowSums(dataPD$count[match(UPS1$ac,dataPD$annot[,1]),,"NoOfPeptides"], na.rm=TRUE)
datUPS1[,2,2] <- rowSums(dataMQ$count[match(UPS1$ac,dataMQ$annot[,1]),,1], na.rm=TRUE)
datUPS1[,3,2] <- rowSums(dataPL$count[match(UPS1$ac,dataPL$annot[,1]),,"NoOfPeptides"], na.rm=TRUE)

## ----combRegrScore3, fig.height=6, fig.width=9.5, fig.align="center", echo=TRUE----
layout(matrix(1:4, ncol=2))
par(mar=c(5.5, 2.2, 4, 0.4))
col1 <- RColorBrewer::brewer.pal(9,"YlOrRd")
imageW(datUPS1[,,1], col=col1, tit="Linear regression score", xLab="",yLab="",transp=FALSE)
mtext("red for bad score", cex=0.75)

imageW(log(datUPS1[,,2]), tit="Number of peptides", xLab="",yLab="", col=col1, transp=FALSE)
mtext("dark red for high number of peptides", cex=0.75)

## ratio : regression score vs no of peptides
imageW(datUPS1[,,1]/log(datUPS1[,,2]), col=rev(col1), tit="Regression score / Number of peptides", xLab="",yLab="", transp=FALSE)
mtext("dark red for high (good) lmScore/peptide ratio)", cex=0.75)

## score vs abundance
imageW(datUPS1[,,1]/datUPS1[,,3], col=rev(col1), tit="Regression score / median Abundance", xLab="",yLab="", transp=FALSE)
mtext("dark red for high (good) lmScore/abundance ratio)", cex=0.75)

## ----combScore1, echo=TRUE----------------------------------------------------
## number of groups for clustering
nGr <- 5
chFin <- is.finite(datUPS1[,,"sco"])
if(any(!chFin)) datUPS1[,,"sco"][which(!chFin)] <- -1      # just in case..


## clustering using kMeans
kMx <- stats::kmeans(standardW(datUPS1[,,"sco"], byColumn=FALSE), nGr)$cluster
datUPS1[,,"cluNo"] <- matrix(rep(kMx, dim(datUPS1)[2]), nrow=length(kMx))

geoM <- apply(datUPS1[,,"sco"], 1, function(x) prod(x)^(1/length(x)))        # geometric mean across analysis soft
geoM2 <- lrbind(by(cbind(geoM,datUPS1[,,"sco"], clu=kMx), kMx, function(x) x[order(x[,1],decreasing=TRUE),]))  # organize by clusters
tmp <- tapply(geoM2[,"geoM"], geoM2[,"clu"], median)
geoM2[,"clu"] <- rep(rank(tmp, ties.method="first"), table(kMx))
geoM2 <- geoM2[order(geoM2[,"clu"],geoM2[,"geoM"],decreasing=TRUE),]         # order as decreasing median.per.cluster
geoM2[,"clu"] <- rep(1:max(kMx), table(geoM2[,"clu"])[rank(unique(geoM2[,"clu"]))])    # replace cluster-names to increasing

try(profileAsClu(geoM2[,2:4], geoM2[,"clu"], tit="Clustered Regression Results for UPS1 Proteins", ylab="Linear regression score"))

## ----combScore2, echo=TRUE----------------------------------------------------
datUPS1 <- datUPS1[match(rownames(geoM2), rownames(datUPS1)),,]               # bring in new order
datUPS1[,,"cluNo"] <- geoM2[,"clu"]                                          # update cluster-names

### prepare annotation of UPS proteins
annUPS1 <- dataPL$annot[match(rownames(datUPS1), dataPL$annot[,1]), c(1,3)]
annUPS1[,2] <- substr(sub("_UPS","",sub("generic_ups\\|[[:alnum:]]+-{0,1}[[:digit:]]\\|","",annUPS1[,2])),1,42)

## ----combScore3, echo=TRUE----------------------------------------------------
## index of representative for each cluster  (median position inside cluster)
UPSrep <- tapply(geoM2[,"geoM"], geoM2[,"clu"], function(x) ceiling(length(x)/2)) + c(0, cumsum(table(geoM2[,"clu"]))[-nGr])

## ----regr5star, echo=TRUE-----------------------------------------------------
gr <- 1
useLi <- which(datUPS1[,1,"cluNo"]==gr)
colNa <- c("Protein",paste(colnames(datUPS1), rep(c("slope","logp"), each=ncol(datUPS1)), sep=" "))
try(kable(cbind(annUPS1[useLi,2], signif(datUPS1[useLi,,"slope"],3), signif(datUPS1[useLi,,"logp"],3)),
  caption=paste("Regression details for cluster of the",length(useLi),"best UPS1 proteins "), col.names=colNa, align="l"),silent=TRUE)

## ----regrPlot5star, fig.height=9, fig.width=9.5, fig.align="center", echo=TRUE----
## Plotting the best regressions, this required package wrGraph version 1.2.5 (or higher)
if(packageVersion("wrGraph")  >= "1.2.5"){
  layout(matrix(1:4, ncol=2))
  tit <- paste0(methNa,", ",annUPS1[UPSrep[gr],1])
  try(tm <- linModelSelect(annUPS1[UPSrep[gr],1], dat=dataPD, tit=tit[1], expect=names(grp9), startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE),silent=TRUE)
  try(tm <- linModelSelect(annUPS1[UPSrep[gr],1], dat=dataMQ, tit=tit[2], expect=names(grp9), startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE),silent=TRUE)
  try(tm <- linModelSelect(annUPS1[UPSrep[gr],1], dat=dataPL, tit=tit[3], expect=names(grp9), startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE),silent=TRUE) }

## ----regr4star, echo=TRUE-----------------------------------------------------
gr <- 2
useLi <- which(datUPS1[,1,"cluNo"]==gr)
try(kable(cbind(annUPS1[useLi,2], signif(datUPS1[useLi,,"slope"],3), signif(datUPS1[useLi,,"logp"],3)),
  caption=paste("Regression details for cluster of the",length(useLi),"2nd best UPS1 proteins "), col.names=colNa, align="l"),silent=TRUE)

## ----regrPlot4star, fig.height=9, fig.width=9.5, fig.align="center", echo=TRUE----
if(packageVersion("wrGraph")  >= "1.2.5"){
  layout(matrix(1:4, ncol=2))
  tit <- paste0(methNa,", ",annUPS1[UPSrep[gr],1])
  try(tm <- linModelSelect(annUPS1[UPSrep[gr],1], dat=dataPD, tit=tit[1], expect=names(grp9), startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE),silent=TRUE)
  try(tm <- linModelSelect(annUPS1[UPSrep[gr],1], dat=dataMQ, tit=tit[2], expect=names(grp9), startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE),silent=TRUE)
  try(tm <- linModelSelect(annUPS1[UPSrep[gr],1], dat=dataPL, tit=tit[3], expect=names(grp9), startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE),silent=TRUE) }

## ----regr3star, echo=TRUE-----------------------------------------------------
gr <- 3
useLi <- which(datUPS1[,1,"cluNo"]==gr)
try(kable(cbind(annUPS1[useLi,2], signif(datUPS1[useLi,,"slope"],3), signif(datUPS1[useLi,,"logp"],3)),
  caption="Regression details for 3rd cluster UPS1 proteins ", col.names=colNa, align="l"),silent=TRUE)

## ----regrPlot3star, fig.height=9, fig.width=9.5, fig.align="center", echo=TRUE----
if(packageVersion("wrGraph")  >= "1.2.5"){
  layout(matrix(1:4, ncol=2))
  tit <- paste0(methNa,", ",annUPS1[UPSrep[gr],1])
  try(tm <- linModelSelect(annUPS1[UPSrep[gr],1], dat=dataPD, tit=tit[1], expect=names(grp9), startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE),silent=TRUE)
  try(tm <- linModelSelect(annUPS1[UPSrep[gr],1], dat=dataMQ, tit=tit[2], expect=names(grp9), startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE),silent=TRUE)
  try(tm <- linModelSelect(annUPS1[UPSrep[gr],1], dat=dataPL, tit=tit[3], expect=names(grp9), startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE),silent=TRUE) }

## ----regrPlot2star, fig.height=9, fig.width=9.5, fig.align="center", echo=TRUE----
gr <- 4
useLi <- which(datUPS1[,1,"cluNo"]==gr)
try(kable(cbind(annUPS1[useLi,2], signif(datUPS1[useLi,,"slope"],3), signif(datUPS1[useLi,,"logp"],3)),
  caption="Regression details for 3rd cluster UPS1 proteins ", col.names=colNa, align="l"),silent=TRUE)

if(packageVersion("wrGraph")  >= "1.2.5"){
  layout(matrix(1:4, ncol=2))
  tit <- paste0(methNa,", ",annUPS1[UPSrep[gr],1])
  try(tm <- linModelSelect(annUPS1[UPSrep[gr],1], dat=dataPD, tit=tit[1], expect=names(grp9), startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE),silent=TRUE)
  try(tm <- linModelSelect(annUPS1[UPSrep[gr],1], dat=dataMQ, tit=tit[2], expect=names(grp9), startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE),silent=TRUE)
  try(tm <- linModelSelect(annUPS1[UPSrep[gr],1], dat=dataPL, tit=tit[3], expect=names(grp9), startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE),silent=TRUE) }

## ----regrPlot1star, fig.height=9, fig.width=9.5, fig.align="center", echo=TRUE----
gr <- 5
useLi <- which(datUPS1[,1,"cluNo"]==gr)
try(kable(cbind(annUPS1[useLi,2], signif(datUPS1[useLi,,"slope"],3), signif(datUPS1[useLi,,"logp"],3)),
  caption="Regression details for 5th cluster UPS1 proteins ", col.names=colNa, align="l"),silent=TRUE)
if(packageVersion("wrGraph")  >= "1.2.5"){
  layout(matrix(1:4, ncol=2))
  tit <- paste0(methNa,", ",annUPS1[UPSrep[gr],1])
  try(tm <- linModelSelect(annUPS1[UPSrep[gr],1], dat=dataPD, tit=tit[1], expect=names(grp9), startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE),silent=TRUE)
  try(tm <- linModelSelect(annUPS1[UPSrep[gr],1], dat=dataMQ, tit=tit[2], expect=names(grp9), startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE),silent=TRUE)
  try(tm <- linModelSelect(annUPS1[UPSrep[gr],1], dat=dataPL, tit=tit[3], expect=names(grp9), startLev=1:5, cexXAxis=0.7, logExpect=TRUE, plotGraph=TRUE, silent=TRUE),silent=TRUE) }

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

