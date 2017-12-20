#!/bin/bash
#SBATCH --job-name=simulation.small.nSample50.sigma1.sh
#SBATCH --partition=super
#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --output=./sbatch_output_%j
#SBATCH --error=./sbatch_error_%j



cd /project/bioinformatics/Xiao_lab/s418336/projects/geneck.ena/code/replicate
Rscript auc.replicate.R 10 simulation small 50 1 -p -v
