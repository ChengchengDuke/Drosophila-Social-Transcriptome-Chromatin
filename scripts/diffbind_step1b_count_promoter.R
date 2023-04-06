## load cutnrun peaks as granges objects
library(rtracklayer)
library(tidyverse)
library(csaw)
library(DESeq2)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(GenomicFeatures)
inputs <- commandArgs(trailingOnly = TRUE)

# Inputs
bam_dir <- inputs[1]
bam_files <- list.files(bam_dir,pattern = ".bam$",full.names = T)
output_deseq <- inputs[2]


## Get TSSs
dm6_promoters <- promoters(transcripts(TxDb.Dmelanogaster.UCSC.dm6.ensGene),
          upstream = 500,
          downstream = 500)
seqlevels(dm6_promoters) <- str_remove(seqlevels(dm6_promoters) ,"chr")

##### Count bam files into counts matrix
# Get counts matrix

counts <- csaw::regionCounts(bam.files = bam_files,
                             regions = dm6_promoters,
                             param = csaw::readParam(dedup = FALSE, minq = 20, pe = "both"))


colnames(counts) <- basename(bam_files)


# Get DESEq object

deseq_obj <- DESeqDataSet(counts, design = ~ 1)

saveRDS(deseq_obj, output_deseq)
