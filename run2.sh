#!/bin/bash
#SBATCH -n 1
#SBATCH -o outfile.txt
#SBATCH -e errfile.txt

srun -l --multi-prog run2.conf
