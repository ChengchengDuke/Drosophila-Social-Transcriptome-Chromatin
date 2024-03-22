#!/bin/sh
#SBATCH -J bedgraphtobigwig
#SBATCH --mem=4G
#SBATCH -a 1-16
genome_file=/work/jes157/cutnrun_socialexperience/genome_files/dmel-all-chromosome-r6.45.chrom.sizes

bedgraph_files=bedgraph_files.txt
bedgraph=$(awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i' $bedgraph_files)

bigwig=${bedgraph%.*}.bigWig
bedgraph_sort=${bedgraph%.*}_sorted.bedGraph

sort -k1,1 -k2,2n $bedgraph > $bedgraph_sort
bedGraphToBigWig $bedgraph_sort $genome_file $bigwig