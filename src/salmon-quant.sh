#!/bin/bash
#SBATCH -p normal
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --mem 50G
#SBATCH --mail-type ALL
#SBATCH -D /home/n-z/ts16/biostate_ai
#SBATCH --mail-user ts16@illinois.edu
#SBATCH -J salmon-quant
#SBATCH -o /home/n-z/ts16/biostate_ai/slurm/slurm-salmon-quant.txt
#SBATCH --array 1-24

line=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" src/file_list.txt)

module load Salmon/1.10.0-IGB-gcc-8.2.0

salmon quant -i salmon/salmon_index -l A \
-1 cutadapt/${line}_1_trimmed.fastq.gz \
-2 cutadapt/${line}_2_trimmed.fastq.gz \
-p 4 --seqBias --gcBias --recoverOrphans --validateMappings -o results/${line}
