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

deseq_obj@design <- ~ condition # design 


# Build model
deseq_obj <- DESeq(deseq_obj) # Default is pairwise comparison 

# Test
## Creating a list that looks like this: # list(c("contrast" = "condition", "CS_SH" ,"CS_GH"), c("contrast" = "condition", "DSX_SH" ,"DSX_GH"), ...)
contrast_list <- split(data.frame(contrast = "condition", 
                                  level1 = paste0(levels(colData(deseq_obj)$genotype),"_","SH"),
                                  level2 = paste0(levels(colData(deseq_obj)$genotype),"_","GH")), 1:length(levels(colData(deseq_obj)$genotype)))
names(contrast_list) <- levels(colData(deseq_obj)$genotype)
contrast_list <- lapply(contrast_list, unlist) # After splitting data frame to get the vector storing the contrast (re)

#Iterate over contrasts list to apply the results function to each contrast
# This function does the actual statistical comparison and outputs the results: pvalue, fold change, etc
res_list <- map(contrast_list, ~ DESeq2::results(contrast = .x, 
                                         object = deseq_obj,
                                         lfcThreshold = 0.1, # log2fc threshold
                                         altHypothesis = "greaterAbs",
                                         alpha = 0.05) ) # p value cutoff


# Shrinkage 
res_ashr_list <- map( res_list,~ lfcShrink(dds = deseq_obj,
                                                   res = .x,
                                                    type = "ashr"))

# Write results objects

names(res_ashr_list) <- paste0(out_dir,names(res_ashr_list) ,"_housing_ashr.Rds")

iwalk(res_ashr_list, ~ saveRDS(.x, .y)) 

# Write results tables
res_ashr_gr_df_list <- map2(res_ashr_list, # results objects with shrunk log2fc estimates
                            res_list, # results objects with raw log2fc estimates
                            ~ cbind(data.frame("gene_name"=rownames(deseq_obj),
                                                                       "stat" = .y$stat,
                                                                       "rawlog2FoldChange"=.y$log2FoldChange,
                                                                       "rawpvalue"=.y$pvalue,
                                                                       "rawpadj"=.y$padj),as.data.frame(.x))) # paste together both results
names(res_ashr_gr_df_list) <- str_replace(names(res_ashr_list),".Rds",".tsv")

# Save each dataframe of results to a separate tsv file
iwalk(res_ashr_gr_df_list, ~ write_tsv(.x, .y))

# Write deseq fitted objects
saveRDS(object = deseq_obj, paste0(out_dir,res_names,"_deseq_obj_housing_fitted.Rds"))

