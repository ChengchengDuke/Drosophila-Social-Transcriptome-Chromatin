library(DESeq2)
library(tidyverse)

# Inputs
inputs <- commandArgs(trailingOnly = TRUE)


contrasts_list <- list(c(inputs[1],inputs[2],inputs[3]))

res_names <- inputs[4]
out_dir <- inputs[5]
deseq_obj_file <- inputs[6]

# Load deseq obj
deseq_obj <- readRDS(deseq_obj_file)
colnames(deseq_obj) <- str_remove(basename(colData(deseq_obj)$bam.files),".target.dedup.sorted.bam")
colData(deseq_obj)$housing <- as.factor(str_split(colnames(deseq_obj),pattern = "_",n = 6) %>% map_chr(3))

deseq_obj@design <- ~ housing

# Prefiltering
keep <- rowSums(counts(deseq_obj)) >= 20
deseq_obj <- deseq_obj[keep,]

# Build model
deseq_obj_contrasts <- DESeq(deseq_obj)


# Test
res_list <- map(contrasts_list, function(x) results(contrast = x, 
                                                    object = deseq_obj_contrasts,
                                                    alpha = 0.05) )


# Shrinkage
res_ashr_list <- map(res_list, function(x) lfcShrink(dds = deseq_obj_contrasts,
                                                     res = x,
                                                     type = "ashr"))

res_ashr_list <- res_list
# Get diff ranges

get_ranges_up <- function(res, deseq, padj = 0.05,lfc = 0){
  granges(rowRanges(deseq)[tidyr::replace_na(res$padj,1) < padj & tidyr::replace_na(res$log2FoldChange,0) > lfc])
}
get_ranges_down <- function(res, deseq, padj = 0.05,lfc = 0){
  granges(rowRanges(deseq)[tidyr::replace_na(res$padj,1) < padj & tidyr::replace_na(res$log2FoldChange,0) < lfc])
}

res_up_list <- map(res_ashr_list, get_ranges_up, deseq = deseq_obj_contrasts)
res_down_list <- map(res_ashr_list, get_ranges_down, deseq = deseq_obj_contrasts)

names(res_ashr_list) <- res_names

# Write results objects

names(res_ashr_list) <- paste0(out_dir,res_names,"_ashr.Rds")

iwalk(res_ashr_list, ~ saveRDS(.x, .y))

# Write results tables
deseq_obj_contrasts_ranges <- as.data.frame(granges(rowRanges(deseq_obj_contrasts)))[,1:3]

names(res_ashr_list) <- paste0(out_dir,res_names,"_ashr.tsv")
iwalk(res_ashr_list, ~ saveRDS(cbind(deseq_obj_contrasts_ranges,as.data.frame(.x)),
                               .y))

# Write deseq fitted objects
saveRDS(deseq_obj_contrasts, paste0(out_dir,res_names,"_deseq_obj_fitted.Rds"))

# DiffBind bedfiles
names(res_up_list) <- paste0(out_dir, res_names,"_up.bed")
iwalk(res_up_list, ~ rtracklayer::export.bed(.x,.y))

names(res_down_list) <- paste0(out_dir, res_names,"_down.bed")
iwalk(res_down_list, ~ rtracklayer::export.bed(.x,.y))

