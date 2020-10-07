## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment = "#>")

## ----install, echo=TRUE, eval=FALSE-------------------------------------------
#  ## if you need to install the packages 'wrMisc','wrProteo' and 'wrGraph' from CRAN :
#  install.packages("wrMisc")
#  install.packages("wrProteo")
#  ## The package 'wrGraph' is not obligatory, but it allows making better graphs
#  install.packages("wrGraph")

## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE------------------------
suppressPackageStartupMessages({
    library(wrMisc)
    library(wrProteo)
    library(knitr)
    library(rmarkdown) 
}) 

## ----setup2-------------------------------------------------------------------
library("wrMisc")
library("wrProteo")
# This is wrProteo version no :
packageVersion("wrProteo")

## ----ChemFormMolMass1, echo=TRUE----------------------------------------------
massDeFormula(c("12H12O", "HO", " 2H 1 Se, 6C 2N", "HSeCN", " ", "e"))
# Note, that empty/invalid entries will be returned as a mass of 0.0 .

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

## ----readUCSC1, echo=TRUE-----------------------------------------------------
path1 <- system.file("extdata", package="wrProteo")
gtfFi <- file.path(path1, "UCSC_hg38_chr11extr.gtf")
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

## ----readProteomeDiscoverer, echo=TRUE----------------------------------------
path1 <- system.file("extdata", package="wrProteo")
fiNaPd <- "exampleProtDiscov1.txt"
dataPD <- readPDExport(file=fiNaPd, path=path1, normalizeMeth="median")
## a summary of the quantitation data
summary(dataPD$quant)

## ----readMaxQuant, echo=TRUE--------------------------------------------------
path1 <- system.file("extdata", package="wrProteo")
fiNaMa <- "proteinGroupsMaxQuantUps1.txt"
specPref1 <- c(conta="conta|CON_|LYSC_CHICK", mainSpecies="YEAST", spike="HUMAN_UPS")
dataMQ <- readMaxQuantFile(path1, file=fiNaMa, specPref=specPref1, normalizeMeth="median")
## a summary of the quantitation data
summary(dataMQ$quant)

## ----readProline, echo=TRUE---------------------------------------------------
normalizeMeth <- NULL
path1 <- system.file("extdata", package="wrProteo")
fiNaPl <- "exampleProlineABC.csv"
dataPL <- readProlineFile(file.path(path1,fiNaPl) )
## a summary of the quantitation data
summary(dataPL$quant)

## ----NA_ProteomeDiscoverer, echo=TRUE-----------------------------------------
## Let's inspect NA values as graphic
matrixNAinspect(dataPD$quant, gr=gl(2,3), tit="Tiny example data from ProteomeDiscoverer") 

## ----NA_MaxQuant, echo=TRUE---------------------------------------------------
## Let's inspect NA values as graphic
matrixNAinspect(dataMQ$quant, gr=gl(3,3), tit="Tiny example data from MaxQuant") 

## ----NA_Proline, echo=TRUE----------------------------------------------------
## Let's inspect NA values as graphic
matrixNAinspect(dataPL$quant, gr=as.factor(substr(colnames(dataPL$quant),1,1)),
  tit="Tiny example data from Proline") 

## ----NArepl_ProteomeDiscoverer, echo=TRUE-------------------------------------
## ProteomeDiscoverer
dataPD$NArepl <- matrixNAneighbourImpute(dataPD$quant,gr=gl(2,3),
  tit="Tiny example from ProteomeDiscoverer") 

## ----NArepl_MaxQuant, echo=TRUE-----------------------------------------------
## MaxQuant
dataMQ$NArepl <- matrixNAneighbourImpute(dataMQ$quant,gr=gl(3,3),tit="Tiny example from MaxQuant") 

## ----NArepl_Proline, echo=TRUE------------------------------------------------
## Proline
dataPL$NArepl <- matrixNAneighbourImpute(dataPL$quant,gr=gl(3,4),tit="Tiny example from Proline") 

## ----testRobustToNAimputation_PL1, echo=TRUE----------------------------------
## Impute NA-values repeatedly and run statistical testing after each round of imputations
testPL <- testRobustToNAimputation(dataPL$quant, gr=gl(3,4), lfdrInclude=FALSE) 
## Note, this example is too small for reliable lfdr estimation (would give warning)

## test results: classical BH FDR
head(testPL$BH)
sum(testPL$BH[,1] < 0.05, na.rm=TRUE)

## the data after repeated NA-imputation
head(testPL$datImp)


## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

