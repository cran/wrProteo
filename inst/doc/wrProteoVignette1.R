## ---- include = FALSE---------------------------------------------------------
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
fiNa <-  "conta1.fasta"
## basic reading of Fasta
fasta1 <- readFasta2(file.path(path1,fiNa))
str(fasta1)

## now let's read and further separate details in annotation-fields
fasta1det <- readFasta2(file.path(path1,fiNa), tableOut=TRUE)
str(fasta1det)

## ----readMaxQuant1, fig.height=8, fig.width=9.5, fig.align="center", echo=TRUE----
path1 <- system.file("extdata", package="wrProteo")
dataMQ <- readMaxQuantFile(path1, specPref=NULL, normalizeMeth="median")

## ----readMaxQuant2,  echo=TRUE------------------------------------------------
## the number of lines and colums
dim(dataMQ$quant)
## a summary of the quantitation data
summary(dataMQ$quant[,1:8])        # the first 8 cols

## ----sampNa1,  echo=TRUE------------------------------------------------------
UPSconc <- c(50,125,250,500,2500,5000,12500,25000,50000)  
sampNa <- paste0(rep(UPSconc, each=3),"amol_",rep(1:3,length(UPSconc))) 
grp9 <- paste0(rep(UPSconc,each=3),"amol") 
head(grp9)

## ----metaData1, echo=TRUE-----------------------------------------------------
## Read meta-data from  github.com/bigbio/proteomics-metadata-standard/
pxd001819meta <- readSdrf("PXD001819")

str(pxd001819meta)

## ----NA_MaxQuant, echo=TRUE---------------------------------------------------
## Let's inspect NA values as graphic
matrixNAinspect(dataMQ$quant, gr=grp9, tit="Histogram of Protein Abundances and NA-Neighbours") 

## ----NArepl_MaxQuant, echo=TRUE-----------------------------------------------
## MaxQuant simple NA-imputation (single round)
dataMQimp <- matrixNAneighbourImpute(dataMQ$quant, gr=grp9, tit="Histogram of Imputed and Final Data") 

## ----testRobustToNAimputation_MQ1, echo=TRUE----------------------------------
## Impute NA-values repeatedly and run statistical testing after each round of imputations
testMQ <- testRobustToNAimputation(dataMQ, gr=grp9) 

## the data after repeated NA-imputation
head(testMQ$datImp[,1:8])

## ----PCA1MQ, fig.height=12, fig.width=9.5, fig.align="center", echo=TRUE------
# limit to UPS1 
plotPCAw(testMQ$datImp, sampleGrp=grp9, tit="PCA on Protein Abundances (MaxQuant,NAs imputed)", rowTyName="proteins", useSymb2=0)

## ----MAplot1, fig.height=6.5, fig.width=9.5, fig.align="center", echo=TRUE----
# By default this plots at the first of all pairwise questions
MAplotW(testMQ)

## ----MAplot2, fig.height=6.5, fig.width=9.5, fig.align="center", echo=TRUE----
MAplotW(testMQ, useComp=2, namesNBest="passFC") 

## ----VolcanoPlot1MQ, fig.height=6.5, fig.width=9.5, fig.align="center", echo=TRUE----
## by default the first pairwise comparison is taken
## using the argument 'namesNBest' we can add names from the annotation
VolcanoPlotW(testMQ, useComp=2, namesNBest="passFDR")

## ----results1, echo=TRUE------------------------------------------------------
res1 <- extractTestingResults(testMQ, compNo=1, thrsh=0.05, FCthrs=2)

## ----results2, echo=TRUE------------------------------------------------------
kable(res1[,-1], caption="5%-FDR (BH) Significant results for 1st pairwise set", align="c")

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

