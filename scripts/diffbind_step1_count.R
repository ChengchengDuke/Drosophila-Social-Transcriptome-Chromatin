## load cutnrun peaks as granges objects
library(rtracklayer)
library(tidyverse)
library(csaw)
library(DESeq2)

# Inputs
cutnrun_peaks_dir <- "/work/jes157/cutnrun_socialexperience/differential_analysis/peaks_dir/RNAPolII/"
blacklist_bed <- "/work/jes157/cutnrun_socialexperience/genome_files/dm6-blacklist.v2.bed"
bam_dir <- "/work/jes157/cutnrun_socialexperience/differential_analysis/bams_dir/RNAPolII"
bam_files <- list.files(bam_dir,pattern = ".bam$",full.names = T)
output_consensus <- "./data_output/cutnrun_RNAPolII_macs2_peaks_union.bed"
output_deseq <- "./data_output/cutnrun_RNAPolII_macs2_peaks_union_deseq_counts.Rds"

cutnrun_peaks_files <- list.files(cutnrun_peaks_dir,pattern = ".bed")
cutnrun_peaks_full <- paste0(cutnrun_peaks_dir,cutnrun_peaks_files)


cutnrun_peaks_granges <- lapply(cutnrun_peaks_full, import.bed)

names(cutnrun_peaks_granges) <- basename(cutnrun_peaks_files)


## get union
cutnrun_peaks_union_granges <- GenomicRanges::union(unlist(GRangesList(cutnrun_peaks_granges)),
                                                    unlist(GRangesList(cutnrun_peaks_granges)))
## remove peaks that overlap blacklist

blacklist_gr <- import(blacklist_bed)

seqlevels(blacklist_gr) <- str_remove(seqlevels(blacklist_gr),"chr")

cutnrun_peaks_union_granges <- cutnrun_peaks_union_granges[!overlapsAny(cutnrun_peaks_union_granges, blacklist_gr)]

# Export union
export.bed(cutnrun_peaks_union_granges,output_consensus)


##### Count bam files into counts matrix
# Get counts matrix

counts <- csaw::regionCounts(bam.files = bam_files,
                             regions = cutnrun_peaks_union_granges,
                             param = csaw::readParam(dedup = FALSE, minq = 20, pe = "both"))


colnames(counts) <- basename(bam_files)


# Get DESEq object

deseq_obj <- DESeqDataSet(counts, design = ~ 1)

saveRDS(deseq_obj, output_deseq)

