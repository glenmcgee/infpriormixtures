#!/bin/bash
#SBATCH -c 1    # Request one core
#SBATCH -t 6-12:00  # Runtime in D-HH:MM format max 12:00
#SBATCH --mem=8000                          # Memory total in MB (for all cores)
#SBATCH -o NHANES_infprior_analysis-%A_%a.out # File to which STDOUT will be written, including job ID
#SBATCH -e NHANES_infprior_analysis-%A_%a.err # File to which STDERR will be written, including job ID
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=gmcgee@uwaterloo.ca # Email to which notifications will be sent

module load gcc/7.3.0 R/4.0.2

Rscript NHANES_infprior_analysis.R $1 $2
