library(DESeq2)
library(tidyverse)


# Inputs
inputs <- commandArgs(trailingOnly = TRUE)


res_names <- inputs[1]
out_dir <- inputs[2]
deseq_obj_file <- inputs[3]

full_model <- inputs[4]
reduced_model <- inputs[5]


# Load deseq obj
deseq_obj <- readRDS(deseq_obj_file)
colnames(deseq_obj) <- str_remove(basename(colData(deseq_obj)$bam.files),".target.dedup.sorted.bam")
colData(deseq_obj)$housing <- as.factor(str_split(colnames(deseq_obj),pattern = "_",n = 6) %>% map_chr(3))
colData(deseq_obj)$genotype <- as.factor(str_split(colnames(deseq_obj),pattern = "_",n = 6) %>% map_chr(2))

deseq_obj@design <- as.formula(full_model)

# Prefiltering
keep <- rowSums(counts(deseq_obj)) >= 10
deseq_obj <- deseq_obj[keep,]

# Build model
deseq_obj_contrasts <- DESeq(deseq_obj,test = "LRT",reduced =as.formula(reduced_model))

# Test
res <- results(object = deseq_obj_contrasts, alpha = 0.05)

# Get diff ranges
get_ranges_diff <- function(res, deseq, padj = 0.05){
  granges(rowRanges(deseq)[tidyr::replace_na(res$padj,1) < padj])
}

res_granges_diff <- get_ranges_diff(res, deseq = deseq_obj_contrasts)


# Write results objects
saveRDS(res, paste0(out_dir,res_names,"_lrt_res.Rds"))

# Write results tables
deseq_obj_contrasts_ranges <- as.data.frame(granges(rowRanges(deseq_obj_contrasts)))[,1:3]

write_tsv(cbind(deseq_obj_contrasts_ranges,as.data.frame(res)),paste0(out_dir,res_names,"_lrt_res.tsv"))

# Write deseq fitted objects
saveRDS(deseq_obj_contrasts, paste0(out_dir,res_names,"_lrt_deseq_obj_fitted.Rds"))

# DiffBind bedfiles
rtracklayer::export.bed(res_granges_diff,paste0(out_dir, res_names,"_lrt_diff.bed"))



