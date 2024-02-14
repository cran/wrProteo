## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment = "#>")

## ----install, echo=TRUE, eval=FALSE-------------------------------------------
#  ## if you need to install the packages 'wrMisc','wrProteo' and 'wrGraph' from CRAN :
#  install.packages("wrMisc")
#  install.packages("wrProteo")
#  ## The package 'wrGraph' is not obligatory, but it allows making better graphs
#  install.packages("wrGraph")
#  
#  ## Installation of limma from Bioconductor
#  if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
#  BiocManager::install("limma")

## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE------------------------
suppressPackageStartupMessages({
    library(wrMisc)
    library(wrProteo)
    library(wrGraph)
    library(knitr)
    library(rmarkdown) 
}) 

## ----setup2-------------------------------------------------------------------
## Let's assume this is a fresh R-session
## Get started by loading the packages
library("knitr")
library("wrMisc")
library("wrProteo")
library("wrGraph")
# This is wrProteo version no :
packageVersion("wrProteo")

## ----Vigenttes1, echo=TRUE, eval=FALSE----------------------------------------
#  browseVignettes("wrProteo")

## ----ChemFormMolMass1, echo=TRUE----------------------------------------------
massDeFormula(c("12H12O", "HO", " 2H 1 Se, 6C 2N", "HSeCN", " ", "e"))

# Ignore empty/invalid entries
massDeFormula(c("12H12O", "HO", " 2H 1 Se, 6C 2N", "HSeCN"), rmEmpty=TRUE)


## ----ChemFormMolMass2, echo=TRUE----------------------------------------------
massDeFormula(c("12H12O", "HO", " 2H 1 Se, 6C 2N", "HSeCN"), massTy="aver")

## ----AAseqMolMass, echo=TRUE--------------------------------------------------
AAmass()

## ----AAseqMolMass2, echo=TRUE-------------------------------------------------
## mass of peptide (or protein)
pep1 <- c(aa="AAAA",de="DEFDEF")
convAASeq2mass(pep1, seqN=FALSE)

## ----readFasta, echo=TRUE-----------------------------------------------------
path1 <- system.file('extdata', package='wrProteo')
fiNa <- "conta1.fasta.gz"
## basic reading of Fasta
fasta1 <- readFasta2(file.path(path1, fiNa))
str(fasta1)

## now let's read and further separate details in annotation-fields
fasta1b <- readFasta2(file.path(path1, fiNa), tableOut=TRUE)
str(fasta1b)

## ----treatFasta2, echo=TRUE---------------------------------------------------
dupEntry <- duplicated(fasta1)
table(dupEntry)

## ----treatFasta3, echo=TRUE---------------------------------------------------
fasta3 <- fasta1[which(!dupEntry)]

length(fasta3)

## ----writeFasta1, echo=TRUE, eval=FALSE---------------------------------------
#  writeFasta2(fasta3, fileNa="testWrite.fasta")

## ----readMaxQuant1, fig.height=8, fig.width=9.5, fig.align="center", echo=TRUE----
path1 <- system.file("extdata", package="wrProteo")
dataMQ <- readMaxQuantFile(path1, specPref=NULL, normalizeMeth="median")
## number of lines and columns of quantitation data
dim(dataMQ$quant)

## ----readMaxQuant2, fig.height=8, fig.width=9.5, fig.align="center", echo=TRUE----
## The grouping of replicates
grp9 <- rep(1:9,each=3)
head(grp9)

## special group of proteins (we want to differentiate/ highlight lateron)
UPS1ac <- c("P00915", "P00918", "P01031", "P69905", "P68871", "P41159", "P02768", "P62988",
  "P04040", "P00167", "P01133", "P02144", "P15559", "P62937", "Q06830", "P63165", "P00709", "P06732",
  "P12081", "P61626", "Q15843", "P02753", "P16083", "P63279", "P01008", "P61769", "P55957", "O76070",
  "P08263", "P01344", "P01127", "P10599", "P99999", "P06396", "P09211", "P01112", "P01579", "P02787",
  "O00762", "P51965", "P08758", "P02741", "P05413", "P10145", "P02788", "P10636-8", "P00441", "P01375")

specPrefMQ <- list(conta="CON_|LYSC_CHICK", mainSpecies="OS=Saccharomyces cerevisiae", spike=UPS1ac)

dataMQ <- readMaxQuantFile(path1, specPref=specPrefMQ, suplAnnotFile=TRUE, groupPref=list(lowNumberOfGroups=FALSE), gr=grp9, plotGraph=FALSE)

## the quantifiation data is the same as before
dim(dataMQ$quant)

## ----readMaxQuant3,  echo=TRUE------------------------------------------------
## count of tags based on argument specPref
table(dataMQ$annot[,"SpecType"])

## ----readMaxQuant4,  echo=TRUE------------------------------------------------
dataMQ <- readMaxQuantFile(path1, specPref=specPrefMQ, sdrf="PXD001819", suplAnnotFile=TRUE, groupPref=list(lowNumberOfGroups=FALSE), plotGraph=FALSE)

## ----exportSdrfDraftMaxQuant5,  echo=TRUE-------------------------------------
path1 <- system.file("extdata", package="wrProteo")
fiNaMQ <- "proteinGroups.txt.gz"
dataMQ2 <- readMaxQuantFile(path1, file=fiNaMQ, refLi="mainSpe", sdrf=FALSE, suplAnnotFile=TRUE)
## Here we'll write simply in the current temporary directory of this R-session
exportSdrfDraft(dataMQ2, file.path(tempdir(),"testSdrf.tsv"))

## ----readMaxQuantPeptides,  echo=TRUE-----------------------------------------
MQpepFi1 <- "peptides_tinyMQ.txt.gz"
path1 <- system.file("extdata", package="wrProteo")
specPref1 <- c(conta="conta|CON_|LYSC_CHICK", mainSpecies="YEAST", spec2="HUMAN")
dataMQpep <- readMaxQuantPeptides(path1, file=MQpepFi1, specPref=specPref1, tit="Tiny MaxQuant Peptides")
summary(dataMQpep$quant)

## ----readProteomeDiscovererProt1,  echo=TRUE----------------------------------
fiNa <- "tinyPD_allProteins.txt.gz"
dataPD <- readProteomeDiscovererFile(file=fiNa, path=path1, suplAnnotFile=FALSE, plotGraph=FALSE)
summary(dataPD$quant)

## ----readDiaNN1, fig.height=8, fig.width=9.5, fig.align="center", echo=TRUE----
diaNNFi1 <- "tinyDiaNN1.tsv.gz"
## This file contains much less identifications than one may usually obtain
path1 <- system.file("extdata", package="wrProteo")
## let's define the main species and allow tagging some contaminants
specPref1 <- c(conta="conta|CON_|LYSC_CHICK", mainSpecies="HUMAN")
dataNN <- readDiaNNFile(path1, file=diaNNFi1, specPref=specPref1, tit="Tiny DIA-NN Data", plotGraph=FALSE)
summary(dataNN$quant)

## ----readProlineProt1,  echo=TRUE---------------------------------------------
path1 <- system.file("extdata", package="wrProteo")
fiNa <- "exampleProlineABC.csv.gz"                  # gz compressed data can be read, too
dataPL <- readProlineFile(file=fiNa, path=path1, plotGraph=FALSE)
summary(dataPL$quant[,1:8])

## ----readFragpipe1,  echo=TRUE------------------------------------------------
FPproFi1 <- "tinyFragpipe1.tsv.gz"
## let's define the main species and allow tagging some contaminants
specPref1 <- c(conta="conta|CON_|LYSC_CHICK", mainSpecies="MOUSE")
dataFP <- readFragpipeFile(path1, file=FPproFi1, specPref=specPref1, tit="Tiny Fragpipe Example", plotGraph=FALSE)
summary(dataFP$quant)

## ----readMassChroq1,  echo=TRUE-----------------------------------------------
MCproFi1 <- "tinyMC.RData"
dataMC <- readMassChroQFile(path1, file=MCproFi1, tit="Tiny MassChroq Example", plotGraph=FALSE)
summary(dataMC$quant)

## ----readAlphaPept1,  echo=TRUE-----------------------------------------------
APproFi1 <- "tinyAlpaPeptide.csv.gz"
## let's define the main species and allow tagging some contaminants
specPref1 <- c(conta="conta|CON_|LYSC_CHICK")
dataAP <- readAlphaPeptFile(path1, file=APproFi1, specPref=specPref1, tit="Tiny AlphaPept Example", plotGraph=FALSE)
summary(dataAP$quant)

## ----readWombarP1,  echo=TRUE-------------------------------------------------
WBproFi1 <- "tinyWombCompo1.csv.gz"
## let's define the main species and allow tagging some contaminants
specPref1 <- c(conta="conta|CON_|LYSC_CHICK", mainSpecies="YEAST")
dataWB <- readWombatNormFile(path1, file=WBproFi1, specPref=specPref1, tit="Tiny Wombat-P Example", plotGraph=FALSE)
summary(dataWB$quant)

## ----readSampleMetaData2,  echo=TRUE------------------------------------------
MQsdrf001819Setup <- readSampleMetaData(quantMeth="MQ", sdrf="PXD001819", path=path1, suplAnnotFile="summary.txt.gz", abund=dataMQ$quant)
str(MQsdrf001819Setup)

## ----fuseProteomicsProjects1,  echo=TRUE--------------------------------------
path1 <- system.file("extdata", package="wrProteo")
dataMQ <- readMaxQuantFile(path1, specPref=NULL, normalizeMeth="median")
dataMC <- readMassChroQFile(path1, file="tinyMC.RData", tit="Tiny MassChroq Example", plotGraph=FALSE)
dataFused <- fuseProteomicsProjects(dataMQ, dataMC)
str(dataFused$quant)

## ----NA_MaxQuant, echo=TRUE---------------------------------------------------
## Let's inspect NA values as graphic
matrixNAinspect(dataMQ$quant, gr=grp9, tit="Histogram of Protein Abundances and NA-Neighbours")

## ----NArepl_MaxQuant, echo=TRUE-----------------------------------------------
## MaxQuant simple NA-imputation (single round)
dataMQimp <- matrixNAneighbourImpute(dataMQ$quant, gr=grp9, tit="Histogram of Imputed and Final Data")

## ----testRobustToNAimputation_MQ1, echo=TRUE----------------------------------
## Impute NA-values repeatedly and run statistical testing after each round of imputations
testMQ <- testRobustToNAimputation(dataMQ, gr=grp9)

## Example of the data after repeated NA-imputation
head(testMQ$datImp[,1:6])

## ----PCA1MQ, fig.height=12, fig.width=9.5, fig.align="center", echo=TRUE------
# limit to UPS1
plotPCAw(testMQ$datImp, sampleGrp=grp9, tit="PCA on Protein Abundances (MaxQuant,NAs imputed)", rowTyName="proteins", useSymb2=0)

## ----MAplot1, fig.height=6.5, fig.width=9.5, fig.align="center", echo=TRUE----
# By default this plots at the first of all pairwise questions
MAplotW(testMQ)

## ----MAplot2, fig.height=6.5, fig.width=9.5, fig.align="center", echo=TRUE----
res1 <- NULL
MAplotW(testMQ, useComp=2, namesNBest="passFC")

## ----VolcanoPlot1MQ, fig.height=6.5, fig.width=9.5, fig.align="center", echo=TRUE----
## by default the first pairwise comparison is taken
## using the argument 'namesNBest' we can add names from the annotation
VolcanoPlotW(testMQ, useComp=2, namesNBest="passFDR")

## ----results1, echo=TRUE------------------------------------------------------
res1 <- extractTestingResults(testMQ, compNo=1, thrsh=0.05, FCthrs=2)

## ----results2, echo=TRUE------------------------------------------------------
knitr::kable(res1[,-1], caption="5%-FDR (BH) Significant results for 1st pairwise set", align="c")

## ----readUCSC1, echo=TRUE-----------------------------------------------------
path1 <- system.file("extdata", package="wrProteo")
gtfFi <- file.path(path1, "UCSC_hg38_chr11extr.gtf.gz")
UcscAnnot1 <- readUCSCtable(gtfFi)

# The Ensemble transcript identifyers and their chromosomal locations :
head(UcscAnnot1)

## ----readUCSC2, echo=TRUE-----------------------------------------------------
# Here we'll redo reading the UCSC table, plus immediatley write the file for UniProt conversion
#  (in this vignette we write to tempdir() to keep things tidy)
expFi <- file.path(tempdir(),"deUcscForUniProt2.txt")
UcscAnnot1 <- readUCSCtable(gtfFi, exportFileNa=expFi)

## ----readUniProt1, echo=TRUE--------------------------------------------------
deUniProtFi <- file.path(path1, "deUniProt_hg38chr11extr.tab")
deUniPr1 <- readUniProtExport(UniP=deUniProtFi, deUcsc=UcscAnnot1, targRegion="chr11:1-135,086,622")
str(deUniPr1)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

