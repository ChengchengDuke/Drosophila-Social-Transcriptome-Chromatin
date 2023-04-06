#!/bin/sh
#SBATCH -J genomecov
#SBATCH --mem=4G
#SBATCH -a 1-6
source activate deeptools

bigwig_files=bigwig_gh.txt
bigwig=$(awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i' $bigwig_files)
igg_ctrl=/work/jes157/cutnrun_socialexperience/seacr_igg/IgG_CS_GH_male_brain_R1.target.dedup.sorted.fragments.bigWig

out_dir=/work/jes157/cutnrun_socialexperience/seacr_igg/

sample=${bigwig%.*}
sample=$out_dir${sample##*/}

bigwig_out=${sample}l.bigWig

# Run seacr
bigwigCompare --binSize 10 --bigwig1 $bigwig --bigwig2 $igg_ctrl -o $bigwig_out

