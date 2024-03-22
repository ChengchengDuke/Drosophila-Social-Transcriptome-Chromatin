#!/bin/sh
#SBATCH -J seqtk_sample
#SBATCH --mem=4G
#SBATCH -a 1-16

names_file=fastq_list.txt

fastq1=$(awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i' $names_file | cut -f1)
fastq2=$(awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i' $names_file | cut -f2)

out_dir=./

fastq1_sub=${fastq1%.*}_sub10M.fastq.gz
fastq1_sub=$out_dir${fastq1_sub##*/}
fastq2_sub=${fastq2%.*}_sub10M.fastq.gz
fastq2_sub=$out_dir${fastq2_sub##*/}

zcat $fastq1 | seqtk sample -s100 - 10000000 | gzip > $fastq1_sub
zcat $fastq2 | seqtk sample -s100 - 10000000 | gzip > $fastq2_sub
