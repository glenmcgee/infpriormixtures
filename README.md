# infpriormixtures
Code for reproducing "Integrating Biological Knowledge in Kernel-Based Analyses of Environmental Mixtures and Health"

## Install bsmim2 for functions used to fit a BMIM
```
library(devtools)
devtools::install_github("glenmcgee/bsmim2") 
library(bsmim2)
```

## Download NHANES data
- Download "studypop.csv" at https://github.com/lizzyagibson/SHARP.Mixtures.Workshop/blob/1f2da3a14bb096d99b2c45a69d11053b0ef60088/Data/studypop.csv
- Save csv file in the "NHANES_analysis" folder

## Simulations
Code to replicate simulation scenarios 
- NHANES_infprior_sim.R: code to replicate simulation A--C
- summarize_NHANES_infprior_sim_results.R: code to report results of NHANES analysis (generates tables and figures)
If submitting job on compute cluster via slurm:
- submitJob_NHANES_infprior_sim.sh: submit batch job to run NHANES simulations via slurm
- combineResults.R: function to combine results of simulations

## Analysis
Code to replicate NHANES data analysis.
- NHANES_infprior_analysis.R: code to run NHANES analysis
- report_NHANES_infprior_analysis.Rmd: code to report results of NHANES analysis (generates tables and figures)
If submitting job on compute cluster via slurm:
- submitJob_NHANES_infprior_analysis.sh: submit job to run NHANES analysis via slurm

