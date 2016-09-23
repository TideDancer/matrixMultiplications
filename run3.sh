#!/bin/bash
#SBATCH -n 6
#SBATCH -o outfile.txt
#SBATCH -e errfile.txt

srun -l --multi-prog run3.conf
