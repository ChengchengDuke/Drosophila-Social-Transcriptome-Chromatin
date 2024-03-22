library(tidyverse)
library(clusterProfiler) # Go terms
library(org.Dm.eg.db) # GO annotation

# Functions

# Rename columns of data frame with by appending .y except for the first column
rename_cols_after1 <- function(.x,.y){
  # .x is the data.frame of results
  # .y is the comparison name
  colnames(.x) <- c(colnames(.x[,c(1:1)]), paste0(colnames(.x[,-c(1:1)]),"_",.y))
  return(.x)
}

# Find genes that are in any results tables (union of all genes )
find_common_genes <- function(res_list){
  res_df <- map(res_list,~ .x[,c("gene_name")])
  res_df <- do.call(rbind, res_df)
  res_df <- unique(res_df)
  # Returns a data.frame where there is only one column of gene_name
  return(res_df)
}

# Change names of columns, and left_join recursively 
reshape_res_list_to_df <- function(res_list) {
  # first find common genes
  res_df <- find_common_genes(res_list)
  # purrr::reduce function to join tables with left_join
  res_df <- purrr::reduce(c(list(res_df),imap(res_list, rename_cols_after1)),left_join, by = c("gene_name"))
  return(res_df)
}

# Add sign column to data.frames by a padding cutoff
my_sign <- function(x, cutoff = 0) {
  ifelse(between(x, -cutoff, cutoff),sign(0),sign(x))
}
# Main

# Read data
data_dir <- "data_output/2024-1_diffexpr_output/"

#Res files
res_files <- list.files(data_dir,pattern = "[^$]+[res|ashr].tsv")
res_list <- lapply(paste0(data_dir,res_files), read_tsv)
names(res_list) <- res_files


# Data reshaping
# Clean list names

names(res_list) <- str_remove(names(res_list), ".tsv")
names(res_list) <- str_remove(names(res_list), "_res")
names(res_list) <- str_remove(names(res_list), "_ashr")

# Reshape results from list to DF  
res_df <- reshape_res_list_to_df(res_list)

############
## Add entrezid column
############

symbol_etrez_map_df <- clusterProfiler::bitr(geneID = res_df$gene_name, fromType = "SYMBOL",toType ="ENTREZID",OrgDb = org.Dm.eg.db )
symbol_etrez_map <- symbol_etrez_map_df$ENTREZ 
names(symbol_etrez_map) <-  symbol_etrez_map_df$SYMBOL

res_df$entrezid <- symbol_etrez_map[res_df$gene_name]

############
# Adding sign column to data.frames
############

res_df <- res_df %>%
  # Apply the my_sign function to every column that starts with log2foldchange
  mutate(across(starts_with("log2FoldChange"),my_sign,.names = "sign_{.col}"))

res_df <- res_df %>%
  mutate(sign_paste_genotype_gh = paste0(
    "dsx-gh:",sign_log2FoldChange_DSX_GH_genotype,
    "_fru-gh:",sign_log2FoldChange_FRU_GH_genotype,
    "_Or47b-gh:",sign_log2FoldChange_Or47b_GH_genotype,
    "_Or67d-gh:",sign_log2FoldChange_Or67d_GH_genotype),
    sign_paste_genotype_sh = paste0(
      "dsx-sh:",sign_log2FoldChange_DSX_SH_genotype,
      "_fru-sh:",sign_log2FoldChange_FRU_SH_genotype,
      "_Or47b-sh:",sign_log2FoldChange_Or47b_SH_genotype,
      "_Or67d-sh:",sign_log2FoldChange_Or67d_SH_genotype),
    sign_paste_housing = paste0(
      "cs:",sign_log2FoldChange_CS_housing,
      "_dsx:",sign_log2FoldChange_DSX_housing,
      "_fru:",sign_log2FoldChange_FRU_housing,
      "_Or47b:",sign_log2FoldChange_Or47b_housing,
      "_Or67d:",sign_log2FoldChange_Or67d_housing))


# Write tsv
write_tsv(res_df,file = "data_output/2024-1_diffexpr_output/res_df_all_join.tsv")


