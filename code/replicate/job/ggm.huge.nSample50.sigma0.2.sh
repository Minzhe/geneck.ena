#!/bin/bash
#SBATCH --job-name=ggm.huge.nSample50.sigma0.2.sh
#SBATCH --partition=256G
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --output=./sbatch_output_%j
#SBATCH --error=./sbatch_error_%j



cd /project/bioinformatics/Xiao_lab/s418336/projects/geneck.ena/code/replicate
Rscript auc.replicate.R 10 ggm huge 50 0.2 -p -v
