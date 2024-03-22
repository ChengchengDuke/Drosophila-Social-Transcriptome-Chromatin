library(DESeq2)
# Read data
data_dir <- "data_output/2024-01_diffbind_cutnrun/"

deseq_files <- list.files(data_dir,pattern = "[^$]+_deseq_obj_fitted.Rds")
deseq_list <- lapply(paste0(data_dir,deseq_files), readRDS)
names(deseq_list) <- deseq_files

# VST for PCA
vsd_h3k4me3 <- vst(deseq_list$cutnrun_H3K4me3_seacrigg_peaks_union_deseq_counts_genotype_lrt_deseq_obj_fitted.Rds)
vsd_rnapolii <- vst(deseq_list$cutnrun_RNAPolII_seacrigg_peaks_union_deseq_counts_genotype_lrt_deseq_obj_fitted.Rds)

vsd_h3k4me3$condition <- as.factor(paste0(vsd_h3k4me3$genotype,"_", vsd_h3k4me3$housing))

vsd_rnapolii$condition <- as.factor(paste0(vsd_rnapolii$genotype,"_", vsd_rnapolii$housing))

saveRDS(vsd_h3k4me3, "data_output/2024-01_diffbind_cutnrun/vsd_h3k4me3.Rds")
saveRDS(vsd_rnapolii, "data_output/2024-01_diffbind_cutnrun/vsd_rnapolii.Rds")


