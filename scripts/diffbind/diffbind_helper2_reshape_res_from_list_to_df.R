library(tidyverse)

# Functions

# Rename columns of data frame with by appending .y except for the first 3 columns
rename_cols_after3 <- function(.x,.y){
  colnames(.x) <- c(colnames(.x[,c(1:3)]), paste0(colnames(.x[,-c(1:3)]),"_",.y))
  return(.x)
}

#Find common set of ranges
find_common_ranges <- function(res_list){
  res_df <- map(res_list,~ .x[,c("seqnames","start","end")])
  res_df <- do.call(rbind, res_df)
  res_df <- unique(res_df)
  return(res_df)
}

# Change names of columns, and left_join recursively 
reshape_res_list_to_df <- function(res_list) {
  res_df <- find_common_ranges(res_list)
  res_df <- reduce(c(list(res_df),imap(res_list, rename_cols_after3)),left_join, by = c("seqnames","start","end"))
  return(res_df)
}

# Main

# Read data
data_dir <- "data_output/2024-01_diffbind_cutnrun/"

#Res files
res_files <- list.files(data_dir,pattern = "[^$]+[res|ashr].tsv")
res_list <- lapply(paste0(data_dir,res_files), read_tsv)
names(res_list) <- res_files

#Filter by chromatin mark 
res_list_h3k4me3 <- str_detect(names(res_list),"H3K4me3")
res_list_h3k27ac <- str_detect(names(res_list),"H3K27ac")
res_list_rnapolii <- str_detect(names(res_list),"RNAPolII")


# Data reshaping
# Clean list names
names(res_list) <- str_remove(names(res_list),"cutnrun_")
names(res_list) <- str_remove(names(res_list),"_seacrigg_peaks_union_deseq_counts")
names(res_list) <- str_remove(names(res_list), ".tsv")
names(res_list) <- str_remove(names(res_list), "_res")

# Reshape results from list to DF  
res_df_h3k4me3 <- reshape_res_list_to_df(res_list[res_list_h3k4me3])
res_df_h3k27ac <- reshape_res_list_to_df(res_list[res_list_h3k27ac])
res_df_rnapolii <- reshape_res_list_to_df(res_list[res_list_rnapolii])

# Add sign combinations

# Split  based on log2fc sign combinations

res_df_h3k4me3 <- res_df_h3k4me3 %>%
  mutate(across(starts_with("log2FoldChange"),sign,.names = "sign_{.col}"))

res_df_h3k4me3 <- res_df_h3k4me3 %>%
  mutate(sign_paste_genotype_gh = paste0(
    "dsx-gh:",sign_log2FoldChange_H3K4me3_dsxmut_GH_genotype_ashr,
    "_fru-gh:",sign_log2FoldChange_H3K4me3_frumut_GH_genotype_ashr,
    "_Or47b-gh:",sign_log2FoldChange_H3K4me3_Or47bmut_GH_genotype_ashr,
    "_Or67d-gh:",sign_log2FoldChange_H3K4me3_Or67dmut_GH_genotype_ashr),
    sign_paste_genotype_sh = paste0(
      "dsx-sh:",sign_log2FoldChange_H3K4me3_dsxmut_SH_genotype_ashr,
      "_fru-sh:",sign_log2FoldChange_H3K4me3_frumut_SH_genotype_ashr,
      "_Or47b-sh:",sign_log2FoldChange_H3K4me3_Or47bmut_SH_genotype_ashr,
      "_Or67d-sh:",sign_log2FoldChange_H3K4me3_Or67dmut_SH_genotype_ashr),
    sign_paste_genotype = paste0(
      "dsx-gh:",sign_log2FoldChange_H3K4me3_dsxmut_GH_genotype_ashr,
      "dsx-sh:",sign_log2FoldChange_H3K4me3_dsxmut_SH_genotype_ashr,
      "_fru-gh:",sign_log2FoldChange_H3K4me3_frumut_GH_genotype_ashr,
      "_fru-sh:",sign_log2FoldChange_H3K4me3_frumut_SH_genotype_ashr,
      "_Or47b-gh:",sign_log2FoldChange_H3K4me3_Or47bmut_GH_genotype_ashr,
      "_Or47b-sh:",sign_log2FoldChange_H3K4me3_Or47bmut_SH_genotype_ashr,
      "_Or67d-gh:",sign_log2FoldChange_H3K4me3_Or67dmut_GH_genotype_ashr,
      "_Or67d-sh:",sign_log2FoldChange_H3K4me3_Or67dmut_SH_genotype_ashr))



# Split  based on log2fc sign combinations

res_df_rnapolii <- res_df_rnapolii %>%
  mutate(across(starts_with("log2FoldChange"),sign,.names = "sign_{.col}"))

res_df_rnapolii <- res_df_rnapolii %>%
  mutate(sign_paste_genotype_gh = paste0(
    "dsx-gh:",sign_log2FoldChange_RNAPolII_dsxmut_GH_genotype_ashr,
    "_fru-gh:",sign_log2FoldChange_RNAPolII_frumut_GH_genotype_ashr,
    "_Or47b-gh:",sign_log2FoldChange_RNAPolII_Or47bmut_GH_genotype_ashr,
    "_Or67d-gh:",sign_log2FoldChange_RNAPolII_Or67dmut_GH_genotype_ashr),
    sign_paste_genotype_sh = paste0(
      "dsx-sh:",sign_log2FoldChange_RNAPolII_dsxmut_SH_genotype_ashr,
      "_fru-sh:",sign_log2FoldChange_RNAPolII_frumut_SH_genotype_ashr,
      "_Or47b-sh:",sign_log2FoldChange_RNAPolII_Or47bmut_SH_genotype_ashr,
      "_Or67d-sh:",sign_log2FoldChange_RNAPolII_Or67dmut_SH_genotype_ashr),
    sign_paste_genotype = paste0(
      "dsx-gh:",sign_log2FoldChange_RNAPolII_dsxmut_GH_genotype_ashr,
      "dsx-sh:",sign_log2FoldChange_RNAPolII_dsxmut_SH_genotype_ashr,
      "_fru-gh:",sign_log2FoldChange_RNAPolII_frumut_GH_genotype_ashr,
      "_fru-sh:",sign_log2FoldChange_RNAPolII_frumut_SH_genotype_ashr,
      "_Or47b-gh:",sign_log2FoldChange_RNAPolII_Or47bmut_GH_genotype_ashr,
      "_Or47b-sh:",sign_log2FoldChange_RNAPolII_Or47bmut_SH_genotype_ashr,
      "_Or67d-gh:",sign_log2FoldChange_RNAPolII_Or67dmut_GH_genotype_ashr,
      "_Or67d-sh:",sign_log2FoldChange_RNAPolII_Or67dmut_SH_genotype_ashr))



# Write tsv
write_tsv(res_df_h3k4me3,file = "data_output/2024-01_diffbind_cutnrun/res_df_h3k4me3_all_join.tsv")
write_tsv(res_df_h3k27ac,file = "data_output/2024-01_diffbind_cutnrun/res_df_h3k27ac_all_join.tsv")
write_tsv(res_df_rnapolii,file = "data_output/2024-01_diffbind_cutnrun/res_df_rnapolii_all_join.tsv")
