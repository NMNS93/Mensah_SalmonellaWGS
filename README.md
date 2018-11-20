
# AMR in S. Typhimurium through WGS

**Nana Mensah - 1st July 2018**

This repository contains analysis scripts for the publication "Mensah et al., 2018. Determining antimicrobial resistance in Salmonella typhimurium by Whole Genome Sequencing: A comparison against multiple phenotypic susceptibility testing methods".

This research benchmarks WGS antimicrobial resistance predictions, using ARIBA and the Resfinder database (v.3), with the gold standard broth microdilution method. Additionally, WGS predictions were evaluated in tandem with legacy methods (disk diffusion, agar dilution) that are routinely used due to their utility for AMR surveillance and clinical diagnostics.

## Directories

|Directory|Description|
|---------|---------|
|README.md| High-level information|
|Makefile|GNU make file to run analysis scripts in order|
|data/|Genotype and phenotype data for analysis. Raw data can be found in 'data/genotype' and 'data/phenotype'. Additionally, 'data/clean' contains processed analysis-ready data|
|files/|Files containing parameters passed to scripts e.g. breakpoints, genotype-conversions|
|results/|Outputs from analysis scripts|
|scripts/|Analysis scripts|
|plots/|Plots generated for figures|

## Analysis

The Makefile included contains instructions in the GNU make utility to run the analysis scripts in the correct order. With R and GNU make installed, call `make` from this directory on the command line to run the analysis (Linux tested only).

**Warnings** - As R performs type conversions throughout the scripts, many warnings are printed. An error from make as displayed below is the indicator that the script has not completed:
```
Makefile:17: recipe for target 'output' failed
make: *** [output] Error 1
```

## scripts/

- preprocessing/* - Scripts to create files for analysis. Output to data/clean/
- functions.R - All scripts import functions defined in this file
- tables.R    - Generates data for tables in the manuscript
- create_db.R - Generates a csv that can be filtered to query individual isolate-test concordance results
- analysis.R  - Performs statistical analyses on cleaned input data

## Input data

- genotype.csv ; per-isolate gene presence/absence
- broth_ecoff.csv ; per-isolate broth microdilution results interpreted with ECOFFs
- broth_clinbp.csv ; per-isolate broth microdilution results interpreted with CBPs
- disk_clinbp.csv ; per-isolate disk diffusion results interpreted with CBPs
- agar_clinbp.csv ; per-isolate agar diliution results interpreted with CBPs

## Outputs

All outputs of this analysis can be found in the 'results/' and 'plots/' directories.

** Tested with R version 3.5.1 and RStudio 1.1.463 **
