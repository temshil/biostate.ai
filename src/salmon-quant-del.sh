#!/bin/bash
#SBATCH -p normal
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --mem 50G
#SBATCH --mail-type ALL
#SBATCH -D /home/n-z/ts16
#SBATCH --mail-user ts16@illinois.edu
#SBATCH -J salmon-quant-del
#SBATCH -o /home/n-z/ts16/biostate_ai/slurm/slurm-salmon-quant-del.txt
#SBATCH --array 1-6

line=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" biostate_ai/src/file_list-del.txt)

module load Salmon/1.10.0-IGB-gcc-8.2.0

salmon quant -i biostate_ai/salmon/salmon_index -l A \
-1 DeletionWT_BulkRNAseq_210225/data/${line}_1.fastq.gz \
-2 DeletionWT_BulkRNAseq_210225/data/${line}_2.fastq.gz \
-p 4 --seqBias --gcBias --recoverOrphans --validateMappings -o biostate_ai/results-del/${line}