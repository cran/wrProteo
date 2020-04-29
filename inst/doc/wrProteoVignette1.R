## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment = "#>")

## ----setup, echo=FALSE, messages=FALSE, warnings=FALSE------------------------
suppressPackageStartupMessages({
    library(wrMisc)
    library(wrProteo)
})

## ----ChemFormMolMass, echo=TRUE-----------------------------------------------
massDeFormula(c("12H12O","HO"," 2H 1 Se, 6C 2N","HSeCN"," ","e"))

## ----AAseqMolMass, echo=TRUE--------------------------------------------------
AAmass()

## ----AAseqMolMass2, echo=TRUE-------------------------------------------------
## mass of peptide (or protein)
pep1 <- c(aa="AAAA",de="DEFDEF")
convAASeq2mass(pep1,seqN=FALSE)

## ----readFasta, echo=TRUE-----------------------------------------------------
path1 <- system.file('extdata',package='wrProteo')
fiNa <-  "conta1.fasta"
## basic reading of Fasta
fasta1 <- readFasta2(file.path(path1,fiNa))
str(fasta1)

## now let's read and further separate annotation-fields
fasta2 <- readFasta2(file.path(path1,fiNa),tableOut=TRUE)
str(fasta2)

## ----readProteomeDiscoverer, echo=TRUE----------------------------------------
path1 <- system.file("extdata",package="wrProteo")
fiNa <- "exampleProtDiscov1.txt"
dataPD <- readPDExport(file=fiNa,path=path1)
## a summary of the quantitation data
summary(dataPD$quant)

## ----readMaxQuant, echo=TRUE--------------------------------------------------
fiNa <- "proteinGroupsMaxQuantUps1.txt"
specPref1=c(conta="conta|CON_|LYSC_CHICK",mainSpecies="YEAST",spike="HUMAN_UPS")
dataMQ <- readMaxQuantFile(path1,file=fiNa,specPref=specPref1)
## a summary of the quantitation data
summary(dataMQ$quant)

## ----readProline, echo=TRUE---------------------------------------------------
path1 <- system.file("extdata",package="wrProteo")
fiNa <- "exampleProlineABC.csv"
dataPL <- readProlineFile(file.path(path1,fiNa))
## a summary of the quantitation data
summary(dataPL$quant)

## ----NA_ProteomeDiscoverer, echo=TRUE-----------------------------------------
## Let's inspect NA values as graphic
matrixNAinspect(dataPD$quant,gr=gl(2,3),tit="Tiny example data from ProteomeDiscoverer") 

## ----NA_MaxQuant, echo=TRUE---------------------------------------------------
## Let's inspect NA values as graphic
matrixNAinspect(dataMQ$quant,gr=gl(3,3),tit="Tiny example data from MaxQuant") 

## ----NA_Proline, echo=TRUE----------------------------------------------------
## Let's inspect NA values as graphic
matrixNAinspect(dataPL$quant,gr=as.factor(substr(colnames(dataPL$quant),1,1)),
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

## ----testRobustToNAimputation_ProteomeDiscoverer, echo=TRUE-------------------
## Impute NA-values repeatedly and run statistical testing after each round of imputations
testPL <- testRobustToNAimputation(dataPL$quant,gr=gl(3,4),lfdrInclude=FALSE) 
## Note, this example is too small for reliable lfdr estimation (would give warning)

## test results: classical BH FDR
head(testPL$BH)
sum(testPL$BH[,1] < 0.05,na.rm=TRUE)
## the data after repeated NA-imputation
head(testPL$datImp)


## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

