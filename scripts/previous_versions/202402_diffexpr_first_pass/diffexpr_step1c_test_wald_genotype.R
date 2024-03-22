library(DESeq2)
library(tidyverse)
# Functions

# Inputs
inputs <- commandArgs(trailingOnly = TRUE)


res_names <- inputs[1]
out_dir <- inputs[2]
deseq_obj_file <- inputs[3]

# Load deseq obj
deseq_obj <- readRDS(deseq_obj_file)

deseq_obj@design <- ~ condition


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
res_ashr_list <- map( res_list,~ lfcShrink(dds = deseq_obj,
                                                   res = .x,
                                                    type = "ashr"))
# Try no shrinkage
#res_ashr_list <- res_list

# Write results objects

names(res_ashr_list) <- paste0(out_dir,names(res_ashr_list) ,"_genotype_ashr.Rds")

iwalk(res_ashr_list, ~ saveRDS(.x, .y))

# Write results tables
res_ashr_gr_df_list <- map2(res_ashr_list,res_list, ~ cbind(data.frame("gene_name"=rownames(deseq_obj),
                                                                       "stat" = .y$stat,
                                                                       "rawlog2FoldChange"=.y$log2FoldChange,
                                                                       "rawpvalue"=.y$pvalue,
                                                                       "rawpadj"=.y$padj),as.data.frame(.x)))
names(res_ashr_gr_df_list) <- str_replace(names(res_ashr_list),".Rds",".tsv")

iwalk(res_ashr_gr_df_list, ~ write_tsv(.x, .y))

# Write deseq fitted objects
saveRDS(object = deseq_obj, paste0(out_dir,res_names,"_deseq_obj_genotype_fitted.Rds"))

