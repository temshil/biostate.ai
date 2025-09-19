#!/bin/bash
#SBATCH -p normal
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 50G
#SBATCH --mail-type ALL
#SBATCH -D /home/n-z/ts16/biostate_ai
#SBATCH --mail-user ts16@illinois.edu
#SBATCH -J salmon-index
#SBATCH -o /home/n-z/ts16/biostate_ai/slurm/slurm-salmon-index.txt

grep "^>" <(gunzip -c salmon/GRCm39.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > salmon/decoys.txt

sed -i.bak -e 's/>//g' salmon/decoys.txt

cat salmon/gencode.vM37.transcripts.fa.gz salmon/GRCm39.primary_assembly.genome.fa.gz > salmon/gentrome.fa.gz

module load Salmon/1.10.0-IGB-gcc-8.2.0

salmon index -t salmon/gentrome.fa.gz -d salmon/decoys.txt -p 12 -i salmon/salmon_index --gencode