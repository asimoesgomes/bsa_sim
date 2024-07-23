#!/bin/bash

#SBATCH --job-name=run_all_cases.R
#SBATCH --error=./sbatch_logs/run-all-cases.error
#SBATCH --account=pi-rglennerster
#SBATCH --partition=caslake
#SBATCH --time=05:0:00
#SBATCH --mem=4G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=squaade@uchicago.edu

module load R/4.2.0
module load gsl/2.7

Rscript run-all-cases.R > sbatch_logs/run-all-cases.Rout

