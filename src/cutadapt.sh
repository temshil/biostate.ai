#!/bin/bash
#SBATCH -p normal
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 10G
#SBATCH --mail-type ALL
#SBATCH -D /home/n-z/ts16/biostate_ai
#SBATCH --mail-user ts16@illinois.edu
#SBATCH -J cutadapt
#SBATCH -o /home/n-z/ts16/biostate_ai/slurm/slurm-cutadapt.txt
#SBATCH --array 1-24

line=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" src/file_list.txt)

module load cutadapt/3.7-IGB-gcc-8.2.0-Python-3.7.2

cutadapt \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o cutadapt/${line}_1_trimmed.fastq.gz -p cutadapt/${line}_2_trimmed.fastq.gz \
fastq/${line}_1.fastq.gz fastq/${line}_2.fastq.gz