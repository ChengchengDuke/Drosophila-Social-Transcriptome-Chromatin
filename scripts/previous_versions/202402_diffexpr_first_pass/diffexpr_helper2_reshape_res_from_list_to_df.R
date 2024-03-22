library(tidyverse)

# Functions

# Rename columns of data frame with by appending .y except for the first 3 columns
rename_cols_after1 <- function(.x,.y){
  colnames(.x) <- c(colnames(.x[,c(1:1)]), paste0(colnames(.x[,-c(1:1)]),"_",.y))
  return(.x)
}

#Find common set of ranges
find_common_genes <- function(res_list){
  res_df <- map(res_list,~ .x[,c("gene_name")])
  res_df <- do.call(rbind, res_df)
  res_df <- unique(res_df)
  return(res_df)
}

# Change names of columns, and left_join recursively 
reshape_res_list_to_df <- function(res_list) {
  res_df <- find_common_genes(res_list)
  res_df <- reduce(c(list(res_df),imap(res_list, rename_cols_after1)),left_join, by = c("gene_name"))
  return(res_df)
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

# Write tsv
write_tsv(res_df,file = "data_output/2024-1_diffexpr_output/res_df_all_join.tsv")


