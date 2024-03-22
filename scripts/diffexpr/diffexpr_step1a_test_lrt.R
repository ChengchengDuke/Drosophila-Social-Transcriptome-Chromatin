library(DESeq2)
library(tidyverse)

# Inputs
inputs <- commandArgs(trailingOnly = TRUE)


res_names <- inputs[1]
out_dir <- inputs[2]
deseq_obj_file <- inputs[3]

full_model <- inputs[4] # Full model includes all the terms in your design (including the one you are testing)
reduced_model <- inputs[5] # Is the null model, it does not include the term you want to test

# Load deseq obj
deseq_obj <- readRDS(deseq_obj_file)

deseq_obj@design <- as.formula(full_model) # Adding the full model as design 

# Build model
deseq_obj_contrasts <- DESeq(deseq_obj,test = "LRT",reduced =as.formula(reduced_model))

# Test
res <- results(object = deseq_obj_contrasts, alpha = 0.05) # Alpha is pval cutoff

# Write results objects
saveRDS(res, paste0(out_dir,res_names,"_lrt_res.Rds")) 

# Write results tables
write_tsv(cbind(data.frame("gene_name"=rownames(deseq_obj_contrasts)),as.data.frame(res)),paste0(out_dir,res_names,"_lrt_res.tsv"))

# Write deseq fitted objects
saveRDS(deseq_obj_contrasts, paste0(out_dir,res_names,"_lrt_deseq_obj_fitted.Rds"))



