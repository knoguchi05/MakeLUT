# MakeLUT
## Creating a correlation plot between Beta and Pt.
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
### output 
- The following files will be created.
----------------------------------
| File Name | description |
|:------------:|:------------:|
| ./PDF/outputFileName_plus(minus).pdf | correlation plot between Beta and Pt |
| ./root/outputFileName_plus(minus).root | Input to FittingLUT.C | 

## Fitting Beta and Pt correlations, creating LUTs
### setup
- Set up the following parameters in the FittingLUT.C file.
----------------------------------
| parameters | description |
|:------------:|:------------:|
| IsChargeSidePlus | true(Charge *Side > 0)| 
| inputRoot | input root (outputFileName_plus(minus).root)|
| outputFileName | |

## run
```sh
$./FitBatch.sh
```
### output 
- The following files will be created.
----------------------------------
| File Name | description |
|:------------:|:------------:|
| ./LUT/outputFileName_newLUT_plus(minus).lut | New LUT |
| ./PDF/outputFileName_fit_plus(minus).pdf | Fitting Result | 
