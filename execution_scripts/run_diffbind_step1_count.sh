#!/bin/sh
#SBATCH -J diffbind
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -a 1-3

module load R/4.1.1-rhel8

samples=(RNAPolII H3K27ac H3K4me3)
sample=${samples[${SLURM_ARRAY_TASK_ID}-1]}

Rscript ./scripts/diffbind_step1_count.R \
  /work/jes157/cutnrun_socialexperience/differential_analysis/peaks_dir/over_igg_peaks/${sample}/ \
  /work/jes157/cutnrun_socialexperience/genome_files/dm6-blacklist.v2.bed \
  /work/jes157/cutnrun_socialexperience/differential_analysis/bams_dir/${sample} \
  ./data_output/cutnrun_${sample}_seacrigg_peaks_union.bed \
  ./data_output/cutnrun_${sample}_seacrigg_peaks_union_deseq_counts.Rds
  
Rscript ./scripts/diffbind_step1_count.R \
  /work/jes157/cutnrun_socialexperience/differential_analysis/peaks_dir/over_igg_peaks/${sample}/ \
  /work/jes157/cutnrun_socialexperience/genome_files/dm6-blacklist.v2.bed \
  /work/jes157/cutnrun_socialexperience/differential_analysis/bams_dir/IgG \
  ./data_output/cutnrun_${sample}_IgG_seacrigg_peaks_union.bed \
  ./data_output/cutnrun_${sample}_IgG_seacrigg_peaks_union_deseq_counts.Rds
  
