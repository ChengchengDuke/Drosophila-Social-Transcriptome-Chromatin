# helper_filter_out_low_count_samples.R
library(DESeq2)
library(tidyverse)
deseq_obj <- readRDS("data_input/2024-01_diff_cutandrun_data_output/cutnrun_H3K4me3_seacrigg_peaks_union_deseq_counts.Rds")
deseq_obj <- deseq_obj[,!str_detect(colnames(deseq_obj), "H3K4me3_CS_SH_male_brain_R2")]
deseq_obj <- deseq_obj[,!str_detect(colnames(deseq_obj), "H3K4me3_dsxmut_GH_male_brain_R4")]
deseq_obj <- deseq_obj[,!str_detect(colnames(deseq_obj), "H3K4me3_CS_GH_male_brain_R1")]
deseq_obj <- deseq_obj[,!str_detect(colnames(deseq_obj), "H3K4me3_CS_SH_male_brain_R1")]
deseq_obj <- deseq_obj[,!str_detect(colnames(deseq_obj), "H3K4me3_CS_GH_male_brain_R2")]
deseq_obj <- deseq_obj[,!str_detect(colnames(deseq_obj), "H3K4me3_CS_SH_male_brain_R2")]

saveRDS(deseq_obj,"data_output/2024-01_diff_cutandrun_data_output/cutnrun_H3K4me3_seacrigg_peaks_union_deseq_counts.Rds")

deseq_obj <- readRDS("data_input/2024-01_diff_cutandrun_data_output/cutnrun_RNAPolII_seacrigg_peaks_union_deseq_counts.Rds")

deseq_obj <- deseq_obj[,!str_detect(colnames(deseq_obj), "RNAPolII_CS_GH_male_brain_R1")]
deseq_obj <- deseq_obj[,!str_detect(colnames(deseq_obj), "RNAPolII_CS_SH_male_brain_R1")]
deseq_obj <- deseq_obj[,!str_detect(colnames(deseq_obj), "RNAPolII_CS_GH_male_brain_R2")]
deseq_obj <- deseq_obj[,!str_detect(colnames(deseq_obj), "RNAPolII_CS_SH_male_brain_R2")]
saveRDS(deseq_obj,"data_output/2024-01_diff_cutandrun_data_output/cutnrun_RNAPolII_seacrigg_peaks_union_deseq_counts.Rds")



