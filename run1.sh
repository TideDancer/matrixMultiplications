#!/bin/bash
#SBATCH -n 4
#SBATCH -o outfile.txt
#SBATCH -e errfile.txt
#SBATCH --partition=HaswellPriority
#SBATCH --account=rajasek

srun -l --multi-prog run1.conf
