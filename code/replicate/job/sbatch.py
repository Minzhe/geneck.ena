#################################################################
###                        sbatch.py                          ###
#################################################################
# This python script is to submit job.sh in jobs folder

import os
import glob

cur_dir = os.path.dirname(os.path.realpath(__file__))
all_jobs = glob.glob('ggm.tiny*.sh')

for job in all_jobs:
    os.system('sbatch {}'.format(job))
