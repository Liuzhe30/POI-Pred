import os

for i in range(92):
	with open('cmds/test_' + str(i+1) + '.slurm','w+') as w:
		w.write('#!/bin/bash\n\n#SBATCH --job-name=testev\n#SBATCH --partition=64c512g\n#SBATCH -N 1\n#SBATCH --ntasks-per-node=4\n#SBATCH --output=evjob.out\n#SBATCH --error=evjob.err\n\n')
		w.write('python SNP_check.py -i ' + str(i)) 

with open('cmds/cmds.txt','w+') as w:
	for i in range(92):
		w.write('sbatch cmds/test_' + str(i+1) + '.slurm\n')