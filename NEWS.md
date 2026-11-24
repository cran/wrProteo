# wrProteo 2.2.0

## new
* findOnOff() .. find proteins all NA/absent in one condition and present (ie measured in other) 
* resumeOnOff()

## Bugfixes and Improvements
* All import-functions for proteomics quantitation data have been revised
* combineMultFilterNAimput()
* extractTestingResults() .. enhanced mode 'slim4'
* fuseProteomicsProjects()
* matrixNAinspect()
* matrixNAneighbourImpute()
* readFasta2() .. more flexible when parsing additional non-conventional annotation
* writeFasta2() .. allows adding add'l information
* readSampleMetaData() .. better matching of raw-file names and names form setup, integration reading pairwise setup from 'compreTo'-columns
* matrixNAinspect()
* testRobustToNAimputation()


# wrProteo 2.1.1

## Bugfixes and Improvements
* Improved `extractTestingResults()`     better cretion of sample short-names
* Improved `testRobustToNAimputation()`  better integration of sample short-names
* Improved `readSampleMetaData()`        better integration of sample short-names 
* Improved `combineMultFilterNAimput()`  better integration of sample short-names
* Example sdrf for vignette also provided in package and used in vignette if corresponding page on github not responding (or no internet connection)


# wrProteo 2.1.0

## Bugfixes and Improvements
* Improved `readSampleMetaData()` to better recognize and interpret
    experimental setup information from sdrf-like files.
* Improved `readSdrf()` to eliminate annotation column-names repeated inside the main annotation
* The function `extractTestingResults()` has been re-written to several problems

# wrProteo 2.0.2

## Bugfixes and Improvements	
* From this version on it is possible to custom choose pairwise groups for statistical testing 
    via an extended sdrf-file which will be mined by the function `readSampleMetaData()`
* Partial rewriting of `testRobustToNAimputation()` to accomodate custom choice of pairwise combinations from `readSampleMetaData()`	
* Bugfixes in the function `readProteomeDiscovererFile()` and `extractTestingResults()`
* The functions `readFasta2()`, `combineMultFilterNAimput()` have been rewritten
	