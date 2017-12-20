###############################################################
###                        job.py                           ###
###############################################################
# This python script writes shell scripts and submits them.

import os

cur_dir = os.path.dirname(os.path.realpath(__file__))

############   functions   #############
def writeHead(name, file):
    print('#PBS -S /bin/bash', file=file)
    print('#PBS -q qbrc', file=file)
    print('#PBS -l nodes=1:ppn=8', file=file)
    print('#PBS -e ./{}.error.qsub'.format(name), file=file)
    print('#PBS -o ./{}.log.qsub'.format(name), file=file)
    print('\n\n', file=file)





####################  write job shell script  #########################
sizes = ['tiny', 'small', 'moderate', 'middle', 'large', 'huge']
sigma2s = [0.2, 0.5, 1]
n_samples = [50, 100, 200, 500]

for size in sizes:
    for sigma2 in sigma2s:
        for n_sample in n_samples:
            job_name = size + '.' + 'nSample' + str(n_sample) + '.' + 'sigma' + str(sigma2) + '.sh'

            with open(job_name, 'w') as f:
                ### write header
                writeHead(name=job_name, file=f)
				
				### wd
                wd_cmd = 'cd /qbrc/home/mzhang/projects/geneck.ena/code/analysis.on.simulation'
                print(wd_cmd, file=f)

                ### call R script
                R_cmd = 'Rscript Compare.ROC.R --size {} --sigma2 {} -n {}'.format(size, sigma2, n_sample)
                print(R_cmd, file=f)
