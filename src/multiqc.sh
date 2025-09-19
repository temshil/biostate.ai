#!/bin/bash
#SBATCH -p normal
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 10G
#SBATCH --mail-type ALL
#SBATCH -D /home/n-z/ts16/biostate_ai
#SBATCH --mail-user ts16@illinois.edu
#SBATCH -J multiqc
#SBATCH -o /home/n-z/ts16/biostate_ai/slurm/slurm-multiqc.txt

module load MultiQC/1.28-IGB-gcc-8.2.0-Python-3.10.1

multiqc /home/n-z/ts16/biostate_ai