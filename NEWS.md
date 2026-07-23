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
	