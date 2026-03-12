#!/bin/sh
#SBATCH -J CutAndRun
#SBATCH -c 16
#SBATCH --mem=128G
module load Java/11.0.8

# run bowtie2-build before: bowtie2-build dmel-all-chromosome-r6.45.fasta dm6_bwt2idx https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#getting-started-with-bowtie-2-lambda-phage-example

export NXF_OPTS="-Xms2g -Xmx20g"
# https://nf-co.re/cutandrun/3.2.2
# Genome size:  https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
nextflow run nf-core/cutandrun -r 3.0 --input /work/jes157/20231003_cutandrun_nfrun/fastq_list_cutnrun_all_samples_oct23_v2.csv --outdir /work/jes157/20231003_cutandrun_nfrun/nfrun_cutnrun_all_samples_oct23 --bowtie2 /work/jes157/20231003_cutandrun_nfrun/genome_files/dm6_bwt2idx --gtf /work/jes157/20231003_cutandrun_nfrun/genome_files/dmel-all-r6.45.gtf --fasta /work/jes157/20231003_cutandrun_nfrun/genome_files/dmel-all-chromosome-r6.45.fasta --blacklist /work/jes157/20231003_cutandrun_nfrun/genome_files/dm6-blacklist.v2.bed --macs_gsize 142573017 --normalisation_mode CPM -profile singularity --normalisation_binsize 1 --use_control false --peakcaller seacr,macs2 --macs2_narrow_peak false -c config.txt -resume
