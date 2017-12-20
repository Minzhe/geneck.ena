#!/bin/bash
#SBATCH --job-name=ggm.middle.nSample200.sigma0.2.sh
#SBATCH --partition=128G
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --output=./sbatch_output_%j
#SBATCH --error=./sbatch_error_%j



cd /project/bioinformatics/Xiao_lab/s418336/projects/geneck.ena/code/replicate
Rscript auc.replicate.R 10 ggm middle 200 0.2 -p -v
