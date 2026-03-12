library(DESeq2)
library(tidyverse)
# Functions
pre_filtering <- function(deseq_obj, threhshold = 20){
  keep <- rowSums(counts(deseq_obj) >= threhshold) >= 3
  
  deseq_obj <- deseq_obj[keep,]
  return(deseq_obj)
}


# Inputs
inputs <- commandArgs(trailingOnly = TRUE)


res_names <- inputs[1]
out_dir <- inputs[2]
deseq_obj_file <- inputs[3]


# Load deseq obj
deseq_obj <- readRDS(deseq_obj_file)
colnames(deseq_obj) <- str_remove(basename(colData(deseq_obj)$bam.files),".target.dedup.sorted.bam")
colData(deseq_obj)$housing <- as.factor(str_split(colnames(deseq_obj),pattern = "_",n = 6) %>% map_chr(3))
colData(deseq_obj)$genotype <- as.factor(str_split(colnames(deseq_obj),pattern = "_",n = 6) %>% map_chr(2))
colData(deseq_obj)$condition <- as.factor(paste0(colData(deseq_obj)$genotype, "_",colData(deseq_obj)$housing))
deseq_obj@design <- ~ condition


# Prefiltering
deseq_obj <- pre_filtering(deseq_obj) # remove really low counts
deseq_obj <- deseq_obj[width(ranges(deseq_obj) ) < quantile(width(ranges(deseq_obj) ), probs = c(0.95)),] # remove really large peaks



# Build model
deseq_obj <- DESeq(deseq_obj)

# Test

contrast_list <- data.frame(contrast = "condition",
                                  level1 = levels(colData(deseq_obj)$condition)[-c(1:2)],
                                  level2 = rep(levels(colData(deseq_obj)$condition)[c(1:2)]))
contrast_list_names <- contrast_list$level1

contrast_list <- split(contrast_list, 1:nrow(contrast_list))
names(contrast_list) <- contrast_list_names
contrast_list <- lapply(contrast_list, unlist)

res_list <- map(contrast_list, ~ results(contrast = .x, 
                                         object = deseq_obj,
                                         alpha = 0.01) )


# Shrinkage
#res_ashr_list <- map( res_list,~ lfcShrink(dds = deseq_obj,
#                                                    res = .x,
#                                                    type = "ashr"))
# Try no shrinkage
res_ashr_list <- res_list

# Get diff ranges

get_ranges_up <- function(res, deseq, padj = 0.01,lfc = 0.1){
  granges(rowRanges(deseq)[tidyr::replace_na(res$padj,1) < padj & tidyr::replace_na(res$log2FoldChange,0) > lfc])
}
get_ranges_down <- function(res, deseq, padj = 0.01,lfc = 0.1){
  granges(rowRanges(deseq)[tidyr::replace_na(res$padj,1) < padj & tidyr::replace_na(res$log2FoldChange,0) < lfc])
}

res_up_list <- map(res_ashr_list, ~ get_ranges_up(res = .x, deseq = deseq_obj))
res_down_list <- map(res_ashr_list, ~ get_ranges_down(res = .x, deseq = deseq_obj))

names(res_ashr_list) <- paste0(res_names,"_",names(res_ashr_list) )

# Write results objects

names(res_ashr_list) <- paste0(out_dir,names(res_ashr_list) ,"_genotype_ashr.Rds")

iwalk(res_ashr_list, ~ saveRDS(.x, .y))

# Write results tables
deseq_obj_contrasts_ranges <- as.data.frame(granges(rowRanges(deseq_obj)))[,1:3]

res_ashr_gr_df_list <- map(res_ashr_list, ~ cbind(deseq_obj_contrasts_ranges,as.data.frame(.x)))
names(res_ashr_gr_df_list) <- str_replace(names(res_ashr_list),".Rds",".tsv")

iwalk(res_ashr_gr_df_list, ~ write_tsv(.x, .y))

# Write deseq fitted objects
saveRDS(object = deseq_obj, paste0(out_dir,res_names,"_deseq_obj_genotype_fitted.Rds"))

# DiffBind bedfiles
names(res_up_list) <- str_replace(names(res_ashr_list),".Rds","_up.bed")
iwalk(res_up_list, ~ rtracklayer::export.bed(.x,.y))

names(res_down_list) <- str_replace(names(res_ashr_list),".Rds","_down.bed")
iwalk(res_down_list, ~ rtracklayer::export.bed(.x,.y))

