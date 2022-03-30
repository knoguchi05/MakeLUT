# MakeLUT
1. Creating a correlation plot between Beta and Pt.
### setup
- Set up the input files in the BetaStudy.h file.
- Set up the following parameters in the BetaStudy.C file.
----------------------------------
| parameters | description |
|:------------:|:------------:|
| IsManyDisplay | |
| IsChargeSidePlus | true(Charge *Side > 0)| 
| outputFileName | | 
### run
```sh
$./BetaBatch.sh
```
- The following files will be created.
----------------------------------
| parameters | description |
|:------------:|:------------:|
| outputFileName.pdf | correlation plot between Beta and Pt |
| outputFileName.root | Input to FittingLUT.C | 

2. Fitting Beta and Pt correlations, creating LUTs
### setup

## run
```sh
$./FitBatch.sh
```
