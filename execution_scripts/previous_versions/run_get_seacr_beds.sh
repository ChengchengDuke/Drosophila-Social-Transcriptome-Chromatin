#!/bin/sh
#SBATCH -J seacr_bed
#SBATCH --mem=16G
#SBATCH -c 8
#SBATCH -a 1-16

module load samtools

names_file=bam_list.txt
bam=$(awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i' $names_file)

genome_file=/work/jes157/cutnrun_socialexperience/genome_files/dmel-all-chromosome-r6.45.chrom.sizes

out_dir=/work/jes157/cutnrun_socialexperience/seacr_igg/

sample=${bam%.*}
sample=$out_dir${sample##*/}

sample_bam=$out_dir${sample##*/}.name_sorted.bam


scaling_factor=$(samtools idxstats $bam | awk '{sum+=$3}END{print sum}')
scaling_factor=$(echo "1000000 / $scaling_factor" | bc -l)

samtools sort -n --threads 8 -o $sample_bam $bam

# bedgraph from bam (from seacr's rescommendations)
bedtools bamtobed -bedpe -i $sample_bam > ${sample}.bed
awk '$1==$4 && $6-$2 < 1000 {print $0}' ${sample}.bed > ${sample}.clean.bed
cut -f 1,2,6 ${sample}.clean.bed | sort -k1,1 -k2,2n -k3,3n > ${sample}.fragments.bed
bedtools genomecov -scale $scaling_factor -bg -i ${sample}.fragments.bed -g $genome_file > ${sample}.fragments.bedgraph

rm ${sample}.bed
rm ${sample}.clean.bed
rm ${sample}.fragments.bed


