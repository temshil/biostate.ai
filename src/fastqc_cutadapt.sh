#!/bin/bash
#SBATCH -p normal
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 10G
#SBATCH --mail-type ALL
#SBATCH -D /home/n-z/ts16/biostate_ai
#SBATCH --mail-user ts16@illinois.edu
#SBATCH -J FastQC_cutadapt
#SBATCH -o /home/n-z/ts16/biostate_ai/slurm/slurm-FastQC_cutadapt.txt

module load FastQC/0.12.1-Java-15.0.1

fastqc cutadapt/*.fastq.gz -o fastqc_cutadapt