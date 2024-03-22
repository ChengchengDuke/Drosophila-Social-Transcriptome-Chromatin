library(DESeq2) 

############
# Read data
############

data_dir <- "data_output/2024-1_diffexpr_input/"
data_out <- "data_output/2024-1_diffexpr_output/"
# DESeq files raw

deseq_files <- list.files(data_dir,pattern = "[^$]+rds")
deseq_list <- lapply(paste0(data_dir,deseq_files), readRDS)
names(deseq_list) <- deseq_files

# VST for PCA
vsd_combat <- vst(deseq_list$social_mutants.salmon.merged.gene_counts.rds)

saveRDS(vsd_combat, "data_output/2024-1_diffexpr_output/vsd_combat.Rds")