#!/bin/sh
#SBATCH -J nfcore-rnaseq
#SBATCH -c 8
#SBATCH --mem=64G
module load Java/11.0.8

nextflow run nf-core/rnaseq --save_reference --input 20231002_nfrnaseq_sample_sheet_full_mutants.csv --outdir /work/jes157/20231002_rnaseq_nfrun/nfrnaseq_fly-social-evo_bdg6 --fasta /work/jes157/20231002_rnaseq_nfrun/genomes/Drosophila_melanogaster.BDGP6.46.dna_sm.toplevel.fa.gz --gtf /work/jes157/20231002_rnaseq_nfrun/genomes/Drosophila_melanogaster.BDGP6.46.110.gtf.gz -profile singularity -c config.txt --extra_salmon_quant_args '--seqBias --gcBias' -resume
