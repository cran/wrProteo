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
library(wrMisc)
library(wrProteo)
library(wrGraph)

# Version number for wrProteo :
packageVersion("wrProteo")

## ----functions1, echo=TRUE----------------------------------------------------
## Two small functions we'll need lateron

replSpecType <- function(x, annCol="SpecType", replBy=cbind(old=c("mainSpe","species2"), new=c("Yeast","UPS1"))) {
  ## rename $annot[,"SpecType"] to more specific names
  chCol <- annCol[1] %in% colnames(x$annot)
  if(chCol) { chCol <- which(colnames(x$annot)==annCol[1])
    chIt <- replBy[,1] %in% unique(x$annot[,chCol])    # check items to replace if present
    if(any(chIt)) for(i in which(chIt)) {useLi <- which(x$annot[,chCol] %in% replBy[i,1]); cat("useLi",head(useLi),"\n"); x$annot[useLi,chCol] <- replBy[i,2]}
  } else message(" replSpecType: 'annCol' not found in x$annot !")
  x }
  
replNAProtNames <- function(x,annCol=c("ProteinName","Accession","SpecType")) {
  ## replace in $annot missing ProteinNames by concatenating Accession + SpecType (ie 2nd & 3rd of annCol)
  chCol <- annCol %in% colnames(x$annot)
  if(all(chCol)) {
    chNA <- is.na(x$annot[,annCol[1]])
    if(any(chNA)) x$annot[which(chNA),annCol[1]] <- paste(x$annot[which(chNA),annCol[2]],x$annot[which(chNA),annCol[3]],sep="_")
  } else message(" replNAProtNames: none of the columnnames 'annCol' found in x$annot !")
  x }

## ----readMaxQuant, fig.height=8, fig.width=9.5, fig.align="center", echo=TRUE----
path1 <- system.file("extdata", package="wrProteo")
fiNaMa <- "proteinGroups.txt.gz"
specPrefMQ <- c(conta="CON_|LYSC_CHICK", mainSpecies="OS=Saccharomyces cerevisiae", spike="HUMAN_UPS")

dataMQ <- readMaxQuantFile(path1, file=fiNaMa, specPref=specPrefMQ, refLi="mainSpe")

## ----readMaxQuant2, fig.height=8, fig.width=9.5, fig.align="center", echo=TRUE----
## a summary of the quantitation data
dim(dataMQ$quant)
summary(dataMQ$quant[,1:8])       # the first 8 cols
colnames(dataMQ$annot)[1:12]
table(dataMQ$annot[,"Species"])
table(dataMQ$annot[,"SpecType"])

## ----readProteomeDiscoverer1, fig.height=8, fig.width=9.5, fig.align="center", echo=TRUE----
path1 <- system.file("extdata", package="wrProteo")
fiNaPd <- "pxd001819_PD2.4_Proteins.txt.gz"
 file.exists(file.path(path1,fiNaPd))
## Note: data exported from ProteomeDiscoverer does not have proper column-names 
sampNa <- paste(rep(c(50,125,250,500,2500,5000,12500,25000,50000), each=3),"amol_R",rep(1:3,9),sep="") 
specPrefPD <- c(conta="Bos tauris|Gallus", mainSpecies="OS=Saccharomyces cerevisiae", spike="OS=Homo sapiens")

dataPD <- readPDExport(file=fiNaPd, path=path1, sampleNames=sampNa, refLi="mainSpe", specPref=specPrefPD)

## ----readProteomeDiscoverer2, fig.height=8, fig.width=9.5, fig.align="center", echo=TRUE----
## a summary of the quantitation data
summary(dataPD$quant[,1:8])        # the first 8 cols
dim(dataPD$quant)
colnames(dataPD$annot)[]
#head(dataPD$annot)
table(dataPD$annot[,"Species"])
table(dataPD$annot[,"SpecType"])

## ----ProteomeDiscoverer2, echo=FALSE------------------------------------------
# get ProteomeDiscoverer2 in same order

## ----readProline, echo=FALSE--------------------------------------------------
## shifted for not printing
path1 <- system.file("extdata", package="wrProteo")
#fiNaPl <- "xxProlineABC.csv"
#dataPL <- readProlineFile(file.path(path1,fiNaPl))
## a summary of the quantitation data
#summary(dataPL$quant)

## ----rearrange1, echo=TRUE----------------------------------------------------
# get all results (MaxQuant,ProteomeDiscoverer, ...) in same order
sampNa <- paste0(rep(c(50,125,250,500,2500,5000,12500,25000,50000),each=3),"amol_R",rep(1:3,9))
grp9 <- paste0(rep(c(50,125,250,500,2500,5000,12500,25000,50000),each=3),"amol") 

## it is more convenient to re-order columns this way in each project
dataPD <- corColumnOrder(dataPD,sampNames=sampNa)          # already in good order
dataMQ <- corColumnOrder(dataMQ,sampNames=sampNa) 
#dataPL <- corColumnOrder(dataPL,sampNames=sampNa) 

## ----postTreatm1, echo=FALSE--------------------------------------------------
## The automatic separation by custom (species) search terms filled in just generic category-names
table(dataPD$annot[,"SpecType"])

## Need to rename $annot[,"SpecType"]  
    dataPDx <- dataPD   # backup
    dataMQx <- dataMQ   # backup
dataPD <- replSpecType(dataPD, replBy=cbind(old=c("mainSpe","species2"), new=c("Yeast","UPS1")))
dataMQ <- replSpecType(dataMQ, replBy=cbind(old=c("mainSpe","species2"), new=c("Yeast","UPS1")))

## Need to addres missing ProteinNames (UPS1) due to missing tags in Fasta
dataPD <- replNAProtNames(dataPD) 
dataMQ <- replNAProtNames(dataMQ) 

## ----NA_ProteomeDiscoverer, echo=TRUE-----------------------------------------
## Let's inspect NA values as graphic
matrixNAinspect(dataPD$quant, gr=grp9, tit="ProteomeDiscoverer")  # gl(9,3)

## ----NA_MaxQuant, echo=TRUE---------------------------------------------------
## Let's inspect NA values as graphic
matrixNAinspect(dataMQ$quant, gr=gl(9,3), tit="MaxQuant") 
## why only 24 columns => reprocess ?

## ----NA_Proline, echo=FALSE---------------------------------------------------
## shifted for not printing
## Let's inspect NA values as graphic
statusPL <- NULL
#matrixNAinspect(dataPL$quant, gr=as.factor(substr(colnames(dataPL$quant),1,1)), tit="Tiny example data from Proline") 

## ----nNA1, echo=TRUE----------------------------------------------------------
## Let's look at the number of NAs. Is there an accumulated number in lower UPS1 semples ?
sumNAperGroup(dataPD$raw, grp9) 
sumNAperGroup(dataMQ$raw, grp9) 


## ----testProteomeDiscoverer, echo=TRUE----------------------------------------
## Let's run pairwise-testing for ProteomeDiscoverer
testPD <- testRobustToNAimputation(dataPD$quant, gr=grp9, lfdrInclude=TRUE, annot=dataPD$annot)       # gl(9,3

## ----testMaxQuant, echo=TRUE--------------------------------------------------
## Let's run pairwise-testing for MaxQuant
testMQ <- testRobustToNAimputation(dataMQ$quant, gr=grp9, lfdrInclude=TRUE, annot=dataMQ$annot) 

## ----testProline, echo=FALSE--------------------------------------------------
## shifted for not printing
## Let's run pairwise-testing for Proline
statusPL <- NULL
#testPL <- testRobustToNAimputation(dataPL$quant, gr=gl(3,4), lfdrInclude=FALSE, annot=dataPL$annot) 

## ----testReorganize, echo=TRUE------------------------------------------------
## recuperate imputeded data to main data-object
dataPD$datImp <- testPD$datImp
dataMQ$datImp <- testMQ$datImp

## ----PCA1, fig.height=12, fig.width=9.5, fig.align="center", echo=TRUE--------
plotPCAw(testPD$datImp, sampleGrp=grp9, tit="PCA on ProteomeDiscoverer (NAs imputed)", rowTyName="proteins", useSymb2=0)
plotPCAw(testMQ$datImp ,sampleGrp=grp9, tit="PCA on MaxQuant (NAs imputed)", rowTyName="proteins", useSymb2=1:9)

## ----PCA2, fig.height=12, fig.width=9.5, fig.align="center", echo=TRUE--------
# limit to UPS1 
plotPCAw(testPD$datImp[which(testPD$annot[,"SpecType"]=="UPS1"),], sampleGrp=grp9, tit="PCA on ProteomeDiscoverer, UPS1 only (NAs imputed)",rowTyName="proteins", useSymb2=0)

plotPCAw(testMQ$datImp[which(testMQ$annot[,"SpecType"]=="UPS1"),], sampleGrp=grp9, tit="PCA on MaxQuant, UPS1 only (NAs imputed)",rowTyName="proteins", useSymb2=1:9)


## ----pairWise1, echo=TRUE-----------------------------------------------------
## The names of all the pair-wise comparisons possible
colnames(testPD$BH)

## ----pairWise2, echo=TRUE-----------------------------------------------------
## The number of differentially abundant proteins passing 5% FDR (ProteomeDiscoverer and MaxQuant) 
signCount <- cbind( sig.PD.BH=colSums(testPD$BH < 0.05, na.rm=TRUE), sig.PD.lfdr=if("lfdr" %in% names(testPD)) colSums(testPD$lfdr < 0.05, na.rm=TRUE),
  sig.MQ.BH=colSums(testMQ$BH < 0.05, na.rm=TRUE), sig.MQ.lfdr=if("lfdr" %in% names(testMQ)) colSums(testMQ$lfdr < 0.05, na.rm=TRUE) )

table1 <- numPairDeColNames(testPD$BH, stripTxt="amol", sortByAbsRatio=TRUE)
table1 <- cbind(table1, signCount[table1[,1],])
knitr::kable(table1, caption="All pairwise comparisons (extended from Ramus et al)", align="c")

## ----pairWiseSelect2, echo=TRUE-----------------------------------------------
## In Ramus paper selection
colnames(testPD$BH)[c(2,21,27)]   

## ----pairWiseSelect3, echo=TRUE-----------------------------------------------
## extended selection
useCompNo <- c(2,21,27, 14,15)
colnames(testPD$BH)[useCompNo]

## Let's extract the concentration part to numeric
numNamePart <- numPairDeColNames(testPD$BH, selComp=useCompNo, stripTxt="amol", sortByAbsRatio=TRUE)
head(numNamePart)

## table with concentrations in selected comparisons
table2 <- cbind(numNamePart, signCount[numNamePart[,1],])
knitr::kable(table2, caption="Selected pairwise comparisons (extended from Ramus et al)", align="c")

## ----Volcano0, echo=TRUE------------------------------------------------------
## the selected comparisons to check
cbind(no=useCompNo, name=colnames(testPD$t)[useCompNo])

## ----Volcano1, fig.height=16, fig.width=9.5, fig.align="center", echo=TRUE----
## check presence and good version of package wrGraph
doVolc <- requireNamespace("wrGraph", quietly=TRUE)
if(doVolc) doVolc <- packageVersion("wrGraph") >= "1.0.6"

## ProteomeDiscoverer
layout(matrix(1:6, ncol=2)) 
if(doVolc) {
  for(i in useCompNo) VolcanoPlotW(testPD, useComp=i, FCthrs=2, FdrThrs=0.05, annColor=c(4,2,3),silent=TRUE)}

## ----Volcano2, fig.height=16, fig.width=9.5, fig.align="center", echo=TRUE----
## MaxQuant
layout(matrix(1:6, ncol=2))
if(doVolc) {
  for(i in useCompNo) VolcanoPlotW(testMQ, useComp=i, FCthrs=2, FdrThrs=0.05, annColor=c(4,2,3),silent=TRUE)}

## ----nNA2, echo=TRUE----------------------------------------------------------
## The number of NAs, just the UPS1 proteins (in ProteomeDiscoverer):
sumNAperGroup(dataPD$raw[which(dataPD$annot[,"SpecType"]=="species2"),], grp9) 
sumNAperGroup(dataMQ$raw[which(dataMQ$annot[,"SpecType"]=="species2"),], grp9) 


## ----intraReplicCV1, fig.height=12, fig.width=10, fig.align="center", echo=TRUE----
## combined plot : all data (left), Ups1 (right)
layout(1:2)
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

## decent compromise based on CV : focus on 250 amol  vs 50000 amol
##    ... or for low no of NAs:           2500 amol   vs 50000 amol

 ## for linear modeling rather  500 amol  vs 50000 amol (last of 'high' NA counts)
 
## Ramus:  500  vs  50000    (PD 28/0 NA, MQ 97/1 NA)
##        5000  vs  50000    (PD  0/0 NA, MQ  8/1 NA)
##       12500  vs  25000

## ----linModel0, fig.height=17, fig.width=9.5, fig.align="center", echo=TRUE----
## the quantified UPS1 names
table(dataPD$annot[,"SpecType"])              # 46
table(dataMQ$annot[,"SpecType"])              # 48

## extract names of quantified UPS1-proteins
NamesUpsPD <- dataPD$annot[which(dataPD$annot[,"SpecType"]=="UPS1"),"Accession"]
NamesUpsMQ <- dataMQ$annot[which(dataMQ$annot[,"SpecType"]=="UPS1"),"Accession"]


## ----linModelPD, fig.height=17, fig.width=9.5, fig.align="center", echo=TRUE----
## ProteomeDiscoverer
lmPD <- list(length=length(NamesUpsPD))

layout(matrix(1:12, ncol=2))
lmPD[1:12] <- lapply(NamesUpsPD[1:12], linModelSelect, dat=dataPD, expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, silent=TRUE)
lmPD[13:24] <- lapply(NamesUpsPD[13:24], linModelSelect, dat=dataPD, expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, silent=TRUE)
lmPD[25:36] <- lapply(NamesUpsPD[25:36], linModelSelect, dat=dataPD, expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, silent=TRUE)
lmPD[37:46] <- lapply(NamesUpsPD[37:46], linModelSelect, dat=dataPD, expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, silent=TRUE)
names(lmPD) <- NamesUpsPD

## ----linModelPD2, fig.height=17, fig.width=9.5, fig.align="center", echo=TRUE----
## We make a little summary of regression-results (ProteomeDiscoverer)
lmPDsum <- cbind(pVal=sapply(lmPD,function(x) x$coef[2,4]),logp=NA,slope=sapply(lmPD,function(x) x$coef[2,1]), startFr=sapply(lmPD,function(x) x$startLev), medRawAbund=apply(log2(dataPD$raw[NamesUpsPD,]),1,median,na.rm=TRUE),good=0)

lmPDsum[,"logp"] <- log10(lmPDsum[,"pVal"])
lmPDsum[which(lmPDsum[,"logp"] < -12 & lmPDsum[,"slope"] >0.75),"good"] <- 1
lmPDsum[which(lmPDsum[,"logp"] < -10 & lmPDsum[,"slope"] >0.7),"good"] <- lmPDsum[which(lmPDsum[,"logp"] < -10 & lmPDsum[,"slope"] >0.7),"good"]+ 1

## now we can check the number of high-confidence quantifications (0 means bad linear model) 
table(lmPDsum[,"good"])           # 24 good quantifications

## at which concentration of UPS1 did one et the best regression results ?
table(lmPDsum[,"startFr"])        # most starting at 1

## a brief summary/overview of regression-results
summary(lmPDsum)

## ----linModelMQ, fig.height=17, fig.width=9.5, fig.align="center", echo=TRUE----
## Now for MaxQuant
lmMQ <- list(length=length(NamesUpsMQ))

layout(matrix(1:12, ncol=2))
lmMQ[1:12] <- lapply(NamesUpsMQ[1:12], linModelSelect, dat=dataMQ, expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, silent=TRUE)
lmMQ[13:24] <- lapply(NamesUpsMQ[13:24], linModelSelect, dat=dataMQ, expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, silent=TRUE)
lmMQ[25:36] <- lapply(NamesUpsMQ[25:36], linModelSelect, dat=dataMQ, expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, silent=TRUE)
lmMQ[37:48] <- lapply(NamesUpsMQ[37:48], linModelSelect, dat=dataMQ, expect=grp9, startLev=1:5, cexXAxis=0.7, logExpect=TRUE, silent=TRUE)
names(lmMQ) <- NamesUpsMQ

## ----linModelMQ2, fig.height=17, fig.width=9.5, fig.align="center", echo=TRUE----
## We make a little summary of regression-results (MaxQuant)
## Regressions with bad slope and/or p-value will be marked as O
lmMQsum <- cbind(pVal=sapply(lmMQ,function(x) x$coef[2,4]),logp=NA,slope=sapply(lmMQ,function(x) x$coef[2,1]), startFr=sapply(lmMQ,function(x) x$startLev), medRawAbund=apply(log2(dataMQ$raw[NamesUpsMQ,]),1,median,na.rm=TRUE),good=0)
lmMQsum[,"logp"] <- log10(lmMQsum[,"pVal"])
lmMQsum[which(lmMQsum[,"logp"] < -12 & lmMQsum[,"slope"] >0.75),"good"] <- 1
lmMQsum[which(lmMQsum[,"logp"] < -10 & lmMQsum[,"slope"] >0.7),"good"] <- lmMQsum[which(lmMQsum[,"logp"] < -10 & lmMQsum[,"slope"] >0.7),"good"]+ 1

## now we can check the number of high-confidence quantifications (0 means bad linear model) 
table(lmMQsum[,"good"])           # 26 good quantifications

## at which concentration of UPS1 did one et the best regression results ?
table(lmMQsum[,"startFr"])        # most starting at 5 !

## a brief summary/overview of regression-results
summary(lmMQsum)

## ----linModelPlot1, fig.height=12, fig.width=9.5, fig.align="center", echo=TRUE----
## summary graphics on all indiv protein regressions for ProteomeDiscoverer
layout(matrix(c(1:3,3), ncol=2, byrow=TRUE))
hist(log10(sapply(lmPD,function(x) x$coef[2,4])), br=15,las=1, main="PD: hist of regr p-values",xlab="log10 p-values")     # good p < 1e-12
hist( sapply(lmPD,function(x) x$coef[2,1]), br=15,las=1, main="PD: hist of regr slopes",xlab="slope")     # good 

tit <- "ProteomeDiscoverer, UPS1 regressions :  p-value vs slope"
useCol <- colorAccording2(lmPDsum[,"medRawAbund"], gradTy="rainbow", revCol=TRUE, nEndOmit=14)
plot(lmPDsum[,c(2,3)], main=tit, type="n")   #col=1, bg.col=useCol, pch=20+lmPDsum[,"startFr"],
points(lmPDsum[,c(2,3)], col=1, bg=useCol, pch=20+lmPDsum[,"startFr"],)
legend("topright",paste("best starting from ",1:5), text.col=1, pch=21:25, col=1, pt.bg="white", cex=0.9, xjust=0.5, yjust=0.5)
mtext("fill color according to median (raw) abundance (violet/blue/low -> green -> red/high)",cex=0.9)
  abline(v=c(-12,-10),lty=2,col="grey") ; abline(h=c(0.7,0.75),lty=2,col="grey")

hi1 <- hist(lmPDsum[,"medRawAbund"], plot=FALSE)
legendHist(sort(lmPDsum[,5]), colRamp=useCol[order(lmPDsum[,"medRawAbund"])][cumsum(hi1$counts)], location="bottomleft", legTit="median raw abundance")  #


## ----linModelPlot2, fig.height=12, fig.width=9.5, fig.align="center", echo=TRUE----
## now for MaxQuant
layout(matrix(c(1:3,3), ncol=2, byrow=TRUE))
hist(log10(sapply(lmMQ,function(x) x$coef[2,4])), br=15,las=1, main="MQ: hist of regr p-values",xlab="log10 p-values")     # good p < 1e-12
hist( sapply(lmMQ,function(x) x$coef[2,1]), br=15,las=1, main="MQ: hist of regr slopes",xlab="slope")     # good 

tit <- "MaxQuant, UPS1 regressions :  p-value vs slope"
useCol <- colorAccording2(lmMQsum[,"medRawAbund"], gradTy="rainbow", revCol=TRUE, nEndOmit=14)
plot(lmMQsum[,c(2,3)], main=tit, type="n")   #col=1, bg.col=useCol, pch=20+lmMQsum[,"startFr"],
points(lmMQsum[,c(2,3)], col=1, bg=useCol, pch=20+lmMQsum[,"startFr"],)
legend("topright",paste("best starting from ",1:5), text.col=1, pch=21:25, col=1, pt.bg="white", cex=0.9, xjust=0.5, yjust=0.5)
mtext("fill color according to median (raw) abundance (red/high -> blue/low)",cex=0.9)
  abline(v=c(-12,-10),lty=2,col="grey") ; abline(h=c(0.7,0.75),lty=2,col="grey") 

hi1 <- hist(lmMQsum[,"medRawAbund"], plot=FALSE)
legendHist(sort(lmMQsum[,5]), colRamp=useCol[order(lmMQsum[,"medRawAbund"])][cumsum(hi1$counts)], location="bottomleft", legTit="median raw abundance")  


## ----ROC_single1, echo=TRUE---------------------------------------------------
## single comparison data for ROC
rocPD.2 <- summarizeForROC(testPD, annotCol="SpecType", spec=c("Yeast","UPS1"), columnTest=2, tyThr="BH",overl=F,color=5)       # 12500amol-25000amol
tail(signif(rocPD.2,3))

## ----ROC_main1, echo=TRUE-----------------------------------------------------

layout(1)
rocPD <- lapply(table2[,1],function(x) summarizeForROC(testPD, annotCol="SpecType", spec=c("Yeast","UPS1"), columnTest=x, tyThr="BH", plotROC=FALSE))
rocMQ <- lapply(table2[,1],function(x) summarizeForROC(testMQ, annotCol="SpecType", spec=c("Yeast","UPS1"), columnTest=x, tyThr="BH", plotROC=FALSE))

names(rocPD) <- colnames(testPD$BH)[useCompNo] 
names(rocMQ) <- colnames(testMQ$BH)[useCompNo] 

## ----ROC_PD1, fig.height=9, fig.width=9.5, fig.align="center", echo=TRUE------
layout(1)
colPanel <- 2:6                                              #c(grey(0.4),2:5)
methNa <- paste(table2[,1],", ie",table2[,3],"-",table2[,4])
methNa <- paste0(rep(c("PD","MQ"), each=length(useCompNo)), methNa)
plotROC(rocPD[[1]],rocPD[[2]],rocPD[[3]],rocPD[[4]],rocPD[[5]], col=colPanel, methNames=methNa[1:5], pointSi=0.8, tit="ProteomeDiscoverer at 5 ratios",legCex=1)

## ----ROC_MQ1, fig.height=9, fig.width=9.5, fig.align="center", echo=TRUE------
plotROC(rocMQ[[1]],rocMQ[[2]],rocMQ[[3]],rocMQ[[4]],rocMQ[[5]], col=colPanel, methNames=methNa[6:10], pointSi=0.8, xlim=c(0,0.27),txtLoc=c(0.09,0.3,0.03), tit="MaxQuant selected ratios",legCex=1)

## ----ROC_PD+MQ, fig.height=9, fig.width=9.5, fig.align="center", echo=TRUE----

colPan10 <- rainbow(13)[c(-3,-5,-13)]
plotROC(rocPD[[1]],rocPD[[2]],rocPD[[3]],rocPD[[4]],rocPD[[5]], rocMQ[[1]],rocMQ[[2]],rocMQ[[3]],rocMQ[[4]],rocMQ[[5]], col=colPan10, methNames=methNa, pointSi=0.8, tit="PD and MQ at selected ratios",legCex=1)

## ----ROC_PL1, echo=FALSE------------------------------------------------------
## shifted for not printing
## Proline
statusPL <- NULL

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

