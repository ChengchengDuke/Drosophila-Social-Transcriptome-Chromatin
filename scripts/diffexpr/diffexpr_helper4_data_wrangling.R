library(DESeq2) 
library(clusterProfiler) # Go terms
library(org.Dm.eg.db) # GO annotation
library(tidyverse)

############
# Read data
############
# DESEQ Object
deseq_obj_fitted <- readRDS("data_output/2024-1_diffexpr_output/social_mutants_housing_lrt_deseq_obj_fitted.Rds")

# Results joined files
res_df <- read_tsv(file = "data_output/2024-1_diffexpr_output/res_df_all_join.tsv")

############
### Data transform for volcano plots
############

## Pivot Longer
log2fc_long <- res_df %>%
  dplyr::select(gene_name,entrezid, starts_with("log2FoldChange")) %>%
  pivot_longer(c(-gene_name,-entrezid), names_to = "comparison", values_to = "log2FoldChange") %>%
  mutate(comparison = str_remove(comparison,"log2FoldChange_" ))

rawlog2fc_long <- res_df %>%
  dplyr::select(gene_name,entrezid, starts_with("rawlog2FoldChange")) %>%
  pivot_longer(c(-gene_name,-entrezid), names_to = "comparison", values_to = "rawlog2FoldChange") %>%
  mutate(comparison = str_remove(comparison,"rawlog2FoldChange_" ))

padj_long <- res_df %>%
  dplyr::select(gene_name,entrezid, starts_with("padj")) %>%
  pivot_longer(c(-gene_name,-entrezid), names_to = "comparison", values_to = "padj") %>%
  mutate(comparison = str_remove(comparison,"padj_" ))

pval_long <- res_df %>%
  dplyr::select(gene_name,entrezid, starts_with("pvalue")) %>%
  pivot_longer(c(-gene_name,-entrezid), names_to = "comparison", values_to = "pvalue") %>%
  mutate(comparison = str_remove(comparison,"pvalue_" ))

basemean_long <- res_df %>%
  dplyr::select(gene_name,entrezid, starts_with("baseMean")) %>%
  pivot_longer(c(-gene_name,-entrezid), names_to = "comparison", values_to = "baseMean") %>%
  mutate(comparison = str_remove(comparison,"baseMean_" ))

res_long_list <- list(basemean_long,rawlog2fc_long, log2fc_long, padj_long, pval_long)
# Join every column
res_long_df <- purrr::reduce(res_long_list,left_join,by =c("gene_name","entrezid","comparison") )

############
## Make lists of differentially expressed genes for different combinations of comparisons
############

# Get all genes that are diff expressed in at least one comparison
res_long_df_diff_any <- res_long_df %>%
  drop_na(padj) %>%
  filter(padj < 0.05, abs(log2FoldChange) > 0.1)

## Housing within the same genotype

comp_list <- list() # Initialize empty list

### All
comp_list$housing_all_symbol <- res_long_df_diff_any %>%
  filter(str_detect(comparison,"housing") & ! str_detect(comparison,"lrt")) %>%
  dplyr::select(gene_name, comparison) %>%
  mutate(comparison = str_remove(comparison, "_housing")) %>%
  split(.$comparison) %>%
  map("gene_name")

comp_list$housing_all_entrezid <- res_long_df_diff_any %>%
  filter(str_detect(comparison,"housing") & ! str_detect(comparison,"lrt")) %>%
  dplyr::select(entrezid, comparison) %>%
  mutate(comparison = str_remove(comparison, "_housing")) %>%
  split(.$comparison) %>%
  map("entrezid")

### Up

comp_list$housing_up_symbol <- res_long_df_diff_any %>%
  filter(str_detect(comparison,"housing") & ! str_detect(comparison,"lrt")) %>%
  filter(log2FoldChange > 0) %>% # Separate by direction of change
  dplyr::select(gene_name, comparison) %>%
  mutate(comparison = str_remove(comparison, "_housing")) %>%
  split(.$comparison) %>%
  map("gene_name")

comp_list$housing_up_entrezid <- res_long_df_diff_any %>%
  filter(str_detect(comparison,"housing") & ! str_detect(comparison,"lrt")) %>%
  filter(log2FoldChange > 0) %>%
  dplyr::select(entrezid, comparison) %>%
  mutate(comparison = str_remove(comparison, "_housing")) %>%
  split(.$comparison) %>%
  map("entrezid")

### Down

comp_list$housing_down_symbol <- res_long_df_diff_any %>%
  filter(str_detect(comparison,"housing") & ! str_detect(comparison,"lrt")) %>%
  filter(log2FoldChange < 0) %>%
  dplyr::select(gene_name, comparison) %>%
  mutate(comparison = str_remove(comparison, "_housing")) %>%
  split(.$comparison) %>%
  map("gene_name")

comp_list$housing_down_entrezid <- res_long_df_diff_any %>%
  filter(str_detect(comparison,"housing") & ! str_detect(comparison,"lrt")) %>%
  filter(log2FoldChange < 0) %>%
  dplyr::select(entrezid, comparison) %>%
  mutate(comparison = str_remove(comparison, "_housing")) %>%
  split(.$comparison) %>%
  map("entrezid")

## GH: all to CS

### All

comp_list$genotype_gh_all_symbol <- res_long_df_diff_any %>%
  filter(str_detect(comparison,"_GH") & str_detect(comparison,"genotype") & ! str_detect(comparison,"lrt")) %>%
  dplyr::select(gene_name, comparison) %>%
  mutate(comparison = str_remove(comparison, "_genotype")) %>%
  split(.$comparison) %>%
  map("gene_name")

comp_list$genotype_gh_all_entrezid <- res_long_df_diff_any %>%
  filter(str_detect(comparison,"_GH") & str_detect(comparison,"genotype") & ! str_detect(comparison,"lrt")) %>%
  dplyr::select(entrezid, comparison) %>%
  mutate(comparison = str_remove(comparison, "_genotype")) %>%
  split(.$comparison) %>%
  map("entrezid")

### Up

comp_list$genotype_gh_up_symbol <- res_long_df_diff_any %>%
  filter(str_detect(comparison,"_GH") & str_detect(comparison,"genotype") & ! str_detect(comparison,"lrt")) %>%
  filter(log2FoldChange > 0) %>%
  dplyr::select(gene_name, comparison) %>%
  mutate(comparison = str_remove(comparison, "_genotype")) %>%
  split(.$comparison) %>%
  map("gene_name")

comp_list$genotype_gh_up_entrezid <- res_long_df_diff_any %>%
  filter(str_detect(comparison,"_GH") & str_detect(comparison,"genotype") & ! str_detect(comparison,"lrt")) %>%
  filter(log2FoldChange > 0) %>%
  dplyr::select(entrezid, comparison) %>%
  mutate(comparison = str_remove(comparison, "_genotype")) %>%
  split(.$comparison) %>%
  map("entrezid")

### Down

comp_list$genotype_gh_down_symbol <- res_long_df_diff_any %>%
  filter(str_detect(comparison,"_GH") & str_detect(comparison,"genotype") & ! str_detect(comparison,"lrt")) %>%
  filter(log2FoldChange < 0) %>%
  dplyr::select(gene_name, comparison) %>%
  mutate(comparison = str_remove(comparison, "_genotype")) %>%
  split(.$comparison) %>%
  map("gene_name")

comp_list$genotype_gh_down_entrezid <- res_long_df_diff_any %>%
  filter(str_detect(comparison,"_GH") & str_detect(comparison,"genotype") & ! str_detect(comparison,"lrt")) %>%
  filter(log2FoldChange < 0) %>%
  dplyr::select(entrezid, comparison) %>%
  mutate(comparison = str_remove(comparison, "_genotype")) %>%
  split(.$comparison) %>%
  map("entrezid")

## SH: all to CS

### All

comp_list$genotype_sh_all_symbol <- res_long_df_diff_any %>%
  filter(str_detect(comparison,"_SH") & str_detect(comparison,"genotype") & ! str_detect(comparison,"lrt")) %>%
  dplyr::select(gene_name, comparison) %>%
  mutate(comparison = str_remove(comparison, "_genotype")) %>%
  split(.$comparison) %>%
  map("gene_name")

comp_list$genotype_sh_all_entrezid <- res_long_df_diff_any %>%
  filter(str_detect(comparison,"_SH") & str_detect(comparison,"genotype") & ! str_detect(comparison,"lrt")) %>%
  dplyr::select(entrezid, comparison) %>%
  mutate(comparison = str_remove(comparison, "_genotype")) %>%
  split(.$comparison) %>%
  map("entrezid")

### Up

comp_list$genotype_sh_up_symbol <- res_long_df_diff_any %>%
  filter(str_detect(comparison,"_SH") & str_detect(comparison,"genotype") & ! str_detect(comparison,"lrt")) %>%
  filter(log2FoldChange > 0) %>%
  dplyr::select(gene_name, comparison) %>%
  mutate(comparison = str_remove(comparison, "_genotype")) %>%
  split(.$comparison) %>%
  map("gene_name")

comp_list$genotype_sh_up_entrezid <- res_long_df_diff_any %>%
  filter(str_detect(comparison,"_SH") & str_detect(comparison,"genotype") & ! str_detect(comparison,"lrt")) %>%
  filter(log2FoldChange > 0) %>%
  dplyr::select(entrezid, comparison) %>%
  mutate(comparison = str_remove(comparison, "_genotype")) %>%
  split(.$comparison) %>%
  map("entrezid")

### Down

comp_list$genotype_sh_down_symbol <- res_long_df_diff_any %>%
  filter(str_detect(comparison,"_SH") & str_detect(comparison,"genotype") & ! str_detect(comparison,"lrt")) %>%
  filter(log2FoldChange < 0) %>%
  dplyr::select(gene_name, comparison) %>%
  mutate(comparison = str_remove(comparison, "_genotype")) %>%
  split(.$comparison) %>%
  map("gene_name")

comp_list$genotype_sh_down_entrezid <- res_long_df_diff_any %>%
  filter(str_detect(comparison,"_SH") & str_detect(comparison,"genotype") & ! str_detect(comparison,"lrt")) %>%
  filter(log2FoldChange < 0) %>%
  dplyr::select(entrezid, comparison) %>%
  mutate(comparison = str_remove(comparison, "_genotype")) %>%
  split(.$comparison) %>%
  map("entrezid")

# Define the desired order (prefixes)
desired_order <- c("CS", "Or47b", "Or67d", "FRU", "DSX")

# Apply function to each sublist within the big list
comp_list <- lapply(comp_list, function(sublist) {
  # Ensure we match elements that start with the desired order using grep
  ordered_names <- unlist(lapply(desired_order, function(pattern) {
    grep(paste0("^", pattern), names(sublist), value = TRUE)
  }))
  
  # Reorder the sublist
  sublist <- sublist[ordered_names]
  
  # Return the reordered sublist while keeping original names
  return(sublist)
})


############
# Data transform for heatmaps
############

### Build norm counts z score matrix

norm_cts <- counts(deseq_obj_fitted,norm =TRUE) # Get normalized counts from DESEQ2 object

de <- list()
de$cs_genes_id <- res_long_df %>%
  dplyr::select(gene_name, comparison, padj, log2FoldChange) %>%
  filter(comparison == "CS_housing") %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0.1) %>%
  arrange(log2FoldChange) %>%
  pull(gene_name) %>%
  unique()

de$genotype_gh_genes_id <- res_long_df %>%
  dplyr::select(gene_name, comparison, padj, log2FoldChange) %>%
  filter(str_detect(comparison, "GH_genotype")) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0.1) %>%
  pull(gene_name) %>%
  unique()

de$genotype_sh_genes_id <- res_long_df %>%
  dplyr::select(gene_name, comparison, padj, log2FoldChange) %>%
  filter(str_detect(comparison, "SH_genotype")) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0.1) %>%
  pull(gene_name) %>%
  unique()

de$genotype_Or_gh_genes_id <- res_long_df %>%
  dplyr::select(gene_name, comparison, padj, log2FoldChange) %>%
  filter(str_detect(comparison, "GH_genotype")) %>%
  filter(str_detect(comparison, "Or47b|Or67d")) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0.1) %>%
  pull(gene_name) %>%
  unique()

de$genotype_Or_sh_genes_id <- res_long_df %>%
  dplyr::select(gene_name, comparison, padj, log2FoldChange) %>%
  filter(str_detect(comparison, "SH_genotype")) %>%
  filter(str_detect(comparison, "Or47b|Or67d")) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0.1) %>%
  pull(gene_name) %>%
  unique()

de$genotype_TF_gh_genes_id <- res_long_df %>%
  dplyr::select(gene_name, comparison, padj, log2FoldChange) %>%
  filter(str_detect(comparison, "GH_genotype")) %>%
  filter(str_detect(comparison, "FRU|DSX")) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0.1) %>%
  pull(gene_name) %>%
  unique()

de$genotype_TF_sh_genes_id <- res_long_df %>%
  dplyr::select(gene_name, comparison, padj, log2FoldChange) %>%
  filter(str_detect(comparison, "SH_genotype")) %>%
  filter(str_detect(comparison, "FRU|DSX")) %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0.1) %>%
  pull(gene_name) %>%
  unique()

# Norm counts median per group

# Remove rows with zero counts
norm_cts_filt <- norm_cts[rowSums(norm_cts) > 0,]

# Get mean matrix
norm_cts_median <- as.data.frame(norm_cts_filt) %>% # convert matrix to data frame, makes it compatible with dplyr
  rownames_to_column(var = "gene_symbol") %>% # move rownames to a new column, makes it compatible with dplyr 
  pivot_longer(cols = -gene_symbol, # make data.frame in tidy format
               names_to = "sample_name",
               values_to = "normalizedCounts") %>% 
  mutate(condition = str_remove(sample_name, "_R\\d$")) %>% # Make condition variables 
  group_by(gene_symbol, condition) %>% # Group by experimental condition 
  summarise(normalizedCounts_median = median(normalizedCounts)) %>% # Calculate mean per groups 
  pivot_wider(names_from = condition, # Return from tidy to matrix format 
              values_from = normalizedCounts_median)
# process matrix further as a numeric matrix format with rownames, compatible with ComplexHeatmap
norm_cts_median <- as.data.frame(norm_cts_median) 
rownames(norm_cts_median) <- norm_cts_median$gene_symbol
norm_cts_median <- norm_cts_median %>%
  dplyr::select(- gene_symbol)
norm_cts_median <- as.matrix(norm_cts_median)

### Build log2Fold change matrix
# Select only log2FoldChane variables

log2fc_mat_housing <- res_df %>%
  drop_na(padj_CS_housing) %>%
  filter(padj_CS_housing < 0.05 & abs(log2FoldChange_CS_housing) > 0.1) %>% # Get diff genes
  dplyr::select(gene_name, starts_with("log2FoldChange") & contains("housing") & ! contains("lrt")) %>% # Select foldchange columns
  arrange(log2FoldChange_CS_housing) %>% # Sort by fold change 
  column_to_rownames(var = "gene_name") %>% # Move gene names to rownames
  as.matrix() 

log2fc_mat_genotype_gh<- res_df %>%
  filter(padj_DSX_GH_genotype < 0.05 | 
           padj_FRU_GH_genotype < 0.05 |
           padj_Or67d_GH_genotype < 0.05 |
           padj_Or47b_GH_genotype < 0.05) %>% # Get any gene that is significant in at least one mutant
  dplyr::select(gene_name, starts_with("log2FoldChange") & contains("GH_genotype") & ! contains("lrt")) %>%
  column_to_rownames(var = "gene_name") %>%
  as.matrix()

log2fc_mat_genotype_sh<- res_df %>%
  filter(padj_DSX_SH_genotype < 0.05 | 
           padj_FRU_SH_genotype < 0.05 |
           padj_Or67d_SH_genotype < 0.05 |
           padj_Or47b_SH_genotype < 0.05) %>%
  dplyr::select(gene_name, starts_with("log2FoldChange") & contains("SH_genotype") & ! contains("lrt")) %>%
  column_to_rownames(var = "gene_name") %>%
  as.matrix()

# Removing extra words to make each sample name easier to read in the heatmap
colnames(log2fc_mat_housing) <- str_remove(colnames(log2fc_mat_housing) ,"log2FoldChange_")
colnames(log2fc_mat_housing) <- str_remove(colnames(log2fc_mat_housing) ,"_housing")

colnames(log2fc_mat_genotype_gh) <- str_remove(colnames(log2fc_mat_genotype_gh) ,"log2FoldChange_")
colnames(log2fc_mat_genotype_gh) <- str_remove(colnames(log2fc_mat_genotype_gh) ,"_genotype")

colnames(log2fc_mat_genotype_sh) <- str_remove(colnames(log2fc_mat_genotype_sh) ,"log2FoldChange_")
colnames(log2fc_mat_genotype_sh) <- str_remove(colnames(log2fc_mat_genotype_sh) ,"_genotype")



##########
# significant ashr genotype GH regions
##########

# Filter results dataframe to get significant ashr genotype regions
# Select relevant columns for pairwise LogFC
pairwise_log2fc_genotype_df_ashr_GH <- res_df %>%
  dplyr::filter(padj_DSX_GH_genotype < 0.05 | padj_FRU_GH_genotype < 0.05 | padj_Or47b_GH_genotype < 0.05 | padj_Or67d_GH_genotype < 0.05) %>% # get diff genes by ashr genotype 
  dplyr::select(gene_name, entrezid,sign_paste_genotype_gh,starts_with("log2FoldChange")) %>% # selecting relevant columns
  dplyr::select(!contains("lrt")) %>% # REmove lrt columns
  dplyr::select(!contains("housing")) %>% # remove housing pairwise comparison 
  dplyr::select(!contains("_SH_")) # remove SH columns 

# Filter by absolute Fold change
# Create a matrix with only fold change 
pairwise_log2fc_genotype_mat <- pairwise_log2fc_genotype_df_ashr_GH %>%
  dplyr::select(starts_with("log2FoldChange"))
# Make a logical filterl
logfc_genotype_filter <- rowSums(abs(pairwise_log2fc_genotype_mat) > 0.1) >= 1

# Filter out the dataframe
pairwise_log2fc_genotype_df_ashr_GH <- pairwise_log2fc_genotype_df_ashr_GH %>%
  dplyr::filter(logfc_genotype_filter)

# Change the column names for the data.frame to keep only relevant info 
colnames(pairwise_log2fc_genotype_df_ashr_GH) <- str_remove(colnames(pairwise_log2fc_genotype_df_ashr_GH), "log2FoldChange_")
colnames(pairwise_log2fc_genotype_df_ashr_GH) <- str_remove(colnames(pairwise_log2fc_genotype_df_ashr_GH), "_GH_genotype")

# Change log2fc columns to a single log2c column for compatibility with ggplot2
pairwise_log2fc_genotype_df_ashr_GH <- pairwise_log2fc_genotype_df_ashr_GH %>%
  pivot_longer(cols = - c(gene_name,entrezid,sign_paste_genotype_gh), names_to = "genotype",values_to = "log2FoldChange")





##########
# significant ashr genotype SH regions
##########

# Filter results dataframe to get significant ashr genotype regions
# Select relevant columns for pairwise LogFC
pairwise_log2fc_genotype_df_ashr_SH <- res_df %>%
  dplyr::filter(padj_DSX_SH_genotype < 0.05 | padj_FRU_SH_genotype < 0.05 | padj_Or47b_SH_genotype < 0.05 | padj_Or67d_SH_genotype < 0.05) %>% # get diff genes by ashr genotype 
  dplyr::select(gene_name, entrezid,sign_paste_genotype_sh,starts_with("log2FoldChange")) %>% # selecting relevant columns
  dplyr::select(!contains("lrt")) %>% # REmove lrt columns
  dplyr::select(!contains("housing")) %>% # remove housing pairwise comparison 
  dplyr::select(!contains("_GH_")) # remove GH columns 

# Filter by absolute Fold change
# Create a matrix with only fold change 
pairwise_log2fc_genotype_mat <- pairwise_log2fc_genotype_df_ashr_SH %>%
  dplyr::select(starts_with("log2FoldChange"))
# Make a logical filterl
logfc_genotype_filter <- rowSums(abs(pairwise_log2fc_genotype_mat) > 0.1) >= 1

# Filter out the dataframe
pairwise_log2fc_genotype_df_ashr_SH <- pairwise_log2fc_genotype_df_ashr_SH %>%
  dplyr::filter(logfc_genotype_filter)

# Change the column names for the data.frame to keep only relevant info 
colnames(pairwise_log2fc_genotype_df_ashr_SH) <- str_remove(colnames(pairwise_log2fc_genotype_df_ashr_SH), "log2FoldChange_")
colnames(pairwise_log2fc_genotype_df_ashr_SH) <- str_remove(colnames(pairwise_log2fc_genotype_df_ashr_SH), "_SH_genotype")

# Change log2fc columns to a single log2c column for compatibility with ggplot2
pairwise_log2fc_genotype_df_ashr_SH <- pairwise_log2fc_genotype_df_ashr_SH %>%
  pivot_longer(cols = - c(gene_name,entrezid,sign_paste_genotype_sh), names_to = "genotype",values_to = "log2FoldChange")





##########
# significant lrt genotype regions
##########

# Filter results dataframe to get significant lrt genotype regions
# Select relevant columns for pairwise LogFC
pairwise_log2fc_genotype_df <- res_df %>%
  dplyr::filter(padj_social_mutants_genotype_lrt < 0.05) %>% # get diff genes by lrt genotype 
  dplyr::select(gene_name, entrezid,sign_paste_genotype_gh,starts_with("log2FoldChange")) %>% # selecting relevant columns
  dplyr::select(!contains("lrt")) %>% # REmove lrt columns
  dplyr::select(!contains("housing")) %>% # remove housing pairwise comparison 
  dplyr::select(!contains("_SH_")) # remove SH columns 

# Filter by absolute Fold change
# Create a matrix with only fold change 
pairwise_log2fc_genotype_mat <- pairwise_log2fc_genotype_df %>%
  dplyr::select(starts_with("log2FoldChange"))
# Make a logical filterl
logfc_genotype_filter <- rowSums(abs(pairwise_log2fc_genotype_mat) > 0.1) >= 1

# Filter out the dataframe
pairwise_log2fc_genotype_df <- pairwise_log2fc_genotype_df %>%
  dplyr::filter(logfc_genotype_filter)

# Change the column names for the data.frame to keeep only relevant info 
colnames(pairwise_log2fc_genotype_df) <- str_remove(colnames(pairwise_log2fc_genotype_df), "log2FoldChange_")
colnames(pairwise_log2fc_genotype_df) <- str_remove(colnames(pairwise_log2fc_genotype_df), "_GH_genotype")

# Change log2fc columns to a single log2c column for compatibility with ggplot2
pairwise_log2fc_genotype_df <- pairwise_log2fc_genotype_df %>%
  pivot_longer(cols = - c(gene_name,entrezid,sign_paste_genotype_gh), names_to = "genotype",values_to = "log2FoldChange")

##########
# significant lrt housing genes
##########

# Filter results dataframe to get significant lrt genotype regions
# Select relevant columns for pairwise LogFC
pairwise_log2fc_housing_df <- res_df %>%
  dplyr::filter(padj_social_mutants_housing_lrt < 0.05) %>%
  dplyr::select(gene_name, entrezid,sign_paste_housing,starts_with("log2FoldChange")) %>%
  dplyr::select(!contains("lrt")) %>%
  dplyr::select(!contains("genotype")) 

# Filter by absolute Fold change
# Create a matrix with only fold change 
pairwise_log2fc_housing_mat <- pairwise_log2fc_housing_df %>%
  dplyr::select(starts_with("log2FoldChange"))
# Make a logical filterl
logfc_housing_filter <- rowSums(abs(pairwise_log2fc_housing_mat) > 0.1) >= 1

# Filter out the dataframe
pairwise_log2fc_housing_df <- pairwise_log2fc_housing_df %>%
  dplyr::filter(logfc_housing_filter)

# Change the column names for the data.frame to keeep only relevant info 
colnames(pairwise_log2fc_housing_df) <- str_remove(colnames(pairwise_log2fc_housing_df), "log2FoldChange_")
colnames(pairwise_log2fc_housing_df) <- str_remove(colnames(pairwise_log2fc_housing_df), "_housing")

pairwise_log2fc_housing_df <- pairwise_log2fc_housing_df %>%
  pivot_longer(cols = - c(gene_name,entrezid,sign_paste), names_to = "genotype",values_to = "log2FoldChange")


##########
# significant lrt interaction regions
##########

# Filter results dataframe to get significant lrt genotype regions
# Select relevant columns for pairwise LogFC
pairwise_log2fc_interaction_df <- res_df %>%
  dplyr::filter(padj_social_mutants_interaction_lrt < 0.05) %>%
  dplyr::select(gene_name, entrezid, sign_paste_housing,starts_with("log2FoldChange")) %>%
  dplyr::select(!contains("lrt")) %>%
  dplyr::select(!contains("genotype")) 

# Filter by absolute Fold change
# Create a matrix with only fold change 
pairwise_log2fc_interaction_mat <- pairwise_log2fc_interaction_df %>%
  dplyr::select(starts_with("log2FoldChange"))
# Make a logical filterl
logfc_interaction_filter <- rowSums(abs(pairwise_log2fc_interaction_mat) > 0.1) >= 1

# Filter out the dataframe
pairwise_log2fc_interaction_df <- pairwise_log2fc_interaction_df %>%
  dplyr::filter(logfc_interaction_filter)

# Change the column names for the data.frame to keeep only relevant info 
colnames(pairwise_log2fc_interaction_df) <- str_remove(colnames(pairwise_log2fc_interaction_df), "log2FoldChange_")
colnames(pairwise_log2fc_interaction_df) <- str_remove(colnames(pairwise_log2fc_interaction_df), "_housing")

pairwise_log2fc_interaction_df <- pairwise_log2fc_interaction_df %>%
  pivot_longer(cols = - c(gene_name,entrezid,sign_paste), names_to = "genotype",values_to = "log2FoldChange")

##########
# Annotate deseq obj with behavior data
##########
# Read behavior data
behavior_data <- read_tsv("data_input/behavior_data/MF_SP_w1118_CI.txt",col_names = FALSE)

# Transform behavior data
behavior_data <- t(behavior_data) # transpose matrix

colnames(behavior_data) <- c("genotype","GH","SH") # add column names
behavior_data <- behavior_data[-1,] # remove extra row that indicated column names
behavior_data <- as.data.frame(behavior_data) 
rownames(behavior_data) <- NULL # remove previous rownames
# Convert courtship data columns to numeric 
behavior_data$GH <- as.numeric(behavior_data$GH)
behavior_data$SH <- as.numeric(behavior_data$SH)
# Subset into two tables, one for SH, one for GH
behavior_data_gh <- behavior_data[,c("genotype","GH")]
behavior_data_sh <- behavior_data[,c("genotype","SH")]
colnames(behavior_data_gh) <- c("genotype","CI")
colnames(behavior_data_sh) <- c("genotype","CI")
# Create the housing variable 
behavior_data_gh$housing <- "GH"
behavior_data_sh$housing <- "SH"
# Merge single house and group house together
behavior_data <- rbind(behavior_data_gh, behavior_data_sh)
behavior_data <- drop_na(behavior_data,CI)

# Summarise by mean
behavior_data <- behavior_data %>%
  group_by(genotype, housing) %>%
  summarise(CI = mean(CI)) %>%
  ungroup()
# Change names to match format for genotype names in other tables
behavior_data$genotype <- c("CS","CS","Or47b","Or47b","Or67d","Or67d","FRU","FRU")
# Join two tables by genotype and housing
colData(deseq_obj_fitted) <- DataFrame(left_join(as.data.frame(colData(deseq_obj_fitted)),behavior_data, by = c("genotype", "housing")))

##########
# Save data
##########
# res_long_df
saveRDS(deseq_obj_fitted,file = "data_output/2024-1_diffexpr_output/deseq_obj_fitted_annotated.Rds")

# res_long_df
saveRDS(res_long_df,file = "data_output/2024-1_diffexpr_output/res_long_df.Rds")

# Comparison lists
saveRDS(comp_list,file = "data_output/2024-1_diffexpr_output/comp_list.Rds")

# Counts matrix for heatmaps
saveRDS(norm_cts,file = "data_output/2024-1_diffexpr_output/norm_cts.Rds")
saveRDS(norm_cts_median,file = "data_output/2024-1_diffexpr_output/norm_cts_median.Rds")
saveRDS(log2fc_mat_housing,file = "data_output/2024-1_diffexpr_output/log2fc_mat_housing.Rds")
saveRDS(log2fc_mat_genotype_gh,file = "data_output/2024-1_diffexpr_output/log2fc_mat_genotype_gh.Rds")
saveRDS(log2fc_mat_genotype_sh,file = "data_output/2024-1_diffexpr_output/log2fc_mat_genotype_sh.Rds")

# Pairwise data for ashr diff genes
saveRDS(pairwise_log2fc_genotype_df_ashr_GH,file = "data_output/2024-1_diffexpr_output/pairwise_log2fc_genotype_df_ashr_GH.Rds")
saveRDS(pairwise_log2fc_genotype_df_ashr_SH,file = "data_output/2024-1_diffexpr_output/pairwise_log2fc_genotype_df_ashr_SH.Rds")


# Pairwise data for lrt diff genes
saveRDS(pairwise_log2fc_interaction_df,file = "data_output/2024-1_diffexpr_output/pairwise_log2fc_interaction_df.Rds")
saveRDS(pairwise_log2fc_housing_df,file = "data_output/2024-1_diffexpr_output/pairwise_log2fc_housing_df.Rds")
saveRDS(pairwise_log2fc_genotype_df,file = "data_output/2024-1_diffexpr_output/pairwise_log2fc_genotype_df.Rds")

# DE genes list
saveRDS(de,file = "data_output/2024-1_diffexpr_output/de.Rds")












