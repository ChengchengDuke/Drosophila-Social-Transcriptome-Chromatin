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

deseq_obj@design <- as.formula(full_model)

# Build model
deseq_obj_contrasts <- DESeq(deseq_obj,test = "LRT",reduced =as.formula(reduced_model))

# Test
res <- results(object = deseq_obj_contrasts, alpha = 0.01)

# Write results objects
saveRDS(res, paste0(out_dir,res_names,"_lrt_res.Rds"))

# Write results tables
write_tsv(cbind(data.frame("gene_name"=rownames(deseq_obj_contrasts)),as.data.frame(res)),paste0(out_dir,res_names,"_lrt_res.tsv"))

# Write deseq fitted objects
saveRDS(deseq_obj_contrasts, paste0(out_dir,res_names,"_lrt_deseq_obj_fitted.Rds"))



