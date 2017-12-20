###############################################################
###                        job.py                           ###
###############################################################
# This python script writes shell scripts and submits them.

import os

cur_dir = os.path.dirname(os.path.realpath(__file__))

############   functions   #############
def writeHead(job_name, node, file):
    print('#!/bin/bash', file=file)
    print('#SBATCH --job-name={}'.format(job_name), file=file)
    print('#SBATCH --partition={}'.format(node), file=file)
    print('#SBATCH --nodes=1', file=file)
    print('#SBATCH --time=10-00:00:00', file=file)
    print('#SBATCH --output=./sbatch.output.{}'.format(job_name), file=file)
    print('#SBATCH --error=./sbatch_error.{}'.format(job_name), file=file)
    print('\n\n', file=file)





####################  write job shell script  #########################
simus = ['ggm', 'simulation']
sizes = ['tiny', 'small', 'moderate', 'middle', 'large', 'huge']
sigmas = [0.2, 0.5, 1]
n_samples = [50, 100, 200]
nodes = {'tiny':'super', 'small':'super', 'moderate':'128G', 'middle':'128G', 'large':'256G', 'huge':'256G'}

for simu in simus:
    for size in sizes:
        for sigma in sigmas:
            for n_sample in n_samples:
                job_name = simu + '.' + size + '.' + 'nSample' + str(n_sample) + '.' + 'sigma' + str(sigma) + '.sh'

                with open(job_name, 'w') as f:
                    ### write header
                    writeHead(job_name=job_name, node=nodes[size], file=f)
				
                    ### wd
                    wd_cmd = 'cd /project/bioinformatics/Xiao_lab/s418336/projects/geneck.ena/code/replicate'
                    print(wd_cmd, file=f)

                    ### call R script
                    R_cmd = 'Rscript auc.replicate.R {} {} {} {} {} -p -v'.format(10, simu, size, n_sample, sigma)
                    print(R_cmd, file=f)
