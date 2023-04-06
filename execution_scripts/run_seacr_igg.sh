#!/bin/sh
#SBATCH -J genomecov
#SBATCH --mem=4G
#SBATCH -a 1-16
module load R

bedgraph_files=bedgraph_gh.txt
bedgraph=$(awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i' $bedgraph_files)
igg_ctrl=/work/jes157/cutnrun_socialexperience/seacr_igg/IgG_CS_GH_male_brain_R1.target.dedup.sorted.fragments.bedgraph

out_dir=/work/jes157/cutnrun_socialexperience/seacr_igg/

sample=${bedgraph%.*}
sample=$out_dir${sample##*/}

sample_igg=${sample}_igg
sample_noigg=${sample}_cpm

# Run seacr
SEACR_1.3.sh $bedgraph $igg_ctrl norm stringent $sample_igg

SEACR_1.3.sh $bedgraph 0.01 non stringent $sample_noigg
