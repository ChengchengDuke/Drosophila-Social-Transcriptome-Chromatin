#!/bin/sh
#SBATCH -J diffbind
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -a 1-4

module load R/4.1.1-rhel8

samples=(RNAPolII H3K27ac H3K4me3 IgG)
sample=${samples[${SLURM_ARRAY_TASK_ID}-1]}

Rscript ./scripts/diffbind_step1b_count_promoter.R \
  /work/jes157/cutnrun_socialexperience/differential_analysis/bams_dir/${sample} \
  ./data_output/cutnrun_${sample}_promoter1kb_deseq_counts.Rds
  
