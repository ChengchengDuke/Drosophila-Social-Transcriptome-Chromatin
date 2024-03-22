#!/bin/sh
#SBATCH -J fastqc
#SBATCH --mem=2G
#SBATCH -c 1
#SBATCH -a 1-1

module load FastQC/0.11.7

names_file=fastq_list.txt
fastq=$(awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i' $names_file)

fastqc $fastq

l