######
### Libraries
######

library(tidyverse)
library(GenomicRanges)
library(plyranges)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(org.Dm.eg.db)
library(GenomicFeatures)
library(clusterProfiler)

######
### Functions
######

# Calculate the median value of Normalized counts median per group (summarize replicates into one column)
get_median_per_group <- function(norm_cts) {
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
  return(norm_cts_median)
}

# Filter and Rename matrices: convert from peak-based matrices to 
rownames_filt_and_convert_to_gene <- function(mat,gene_map) {
  # Filter and Rename mats
  # Keep only rows that exist in gene map
  mat <- mat[rownames(mat) %in% names(gene_map),]
  # Change rownames from id to gene name
  rownames(mat) <- gene_map[rownames(mat)]
  # Remove rows that failed to map
  mat <- mat[!is.na(rownames(mat)),]
  # Remove duplicated rows
  mat <- mat[!duplicated(rownames(mat)),]
  return(mat)
}
######
### Read data
######

# Results data
res_df_h3k4me3 <- read_tsv(file = "data_output/2024-01_diffbind_cutnrun/res_df_h3k4me3_all_join.tsv")
res_df_h3k27ac <- read_tsv(file = "data_output/2024-01_diffbind_cutnrun/res_df_h3k27ac_all_join.tsv")
res_df_rnapolii <- read_tsv(file = "data_output/2024-01_diffbind_cutnrun/res_df_rnapolii_all_join.tsv")

# DESEQ Objects
data_dir <- "data_output/2024-01_diffbind_cutnrun/"
deseq_files <- list.files(data_dir,pattern = "[^$]+_deseq_obj_fitted.Rds")
deseq_list <- lapply(paste0(data_dir,deseq_files), readRDS)
names(deseq_list) <- deseq_files

# Data re-shaping

deseq_list_h3k4me3 <- str_detect(names(deseq_list),"H3K4me3")
deseq_list_h3k27ac <- str_detect(names(deseq_list),"H3K27ac")
deseq_list_rnapolii <- str_detect(names(deseq_list),"RNAPolII")

names(deseq_list) <- str_remove(names(deseq_list),"cutnrun_")



#######
### Volcano plots
#######

# Make long dataframe versions of log2fc, padj, etc
log2fc_long <- res_df_h3k4me3 %>%
  dplyr::select(seqnames,start,end, starts_with("log2FoldChange")) %>%
  pivot_longer(c(-seqnames,-start,-end), names_to = "comparison", values_to = "log2FoldChange") %>%
  mutate(comparison = str_remove(comparison,"log2FoldChange_" ))

padj_long <- res_df_h3k4me3 %>%
  dplyr::select(seqnames,start,end, starts_with("padj")) %>%
  pivot_longer(c(-seqnames,-start,-end), names_to = "comparison", values_to = "padj") %>%
  mutate(comparison = str_remove(comparison,"padj_" ))

pval_long <- res_df_h3k4me3 %>%
  dplyr::select(seqnames,start,end, starts_with("pvalue")) %>%
  pivot_longer(c(-seqnames,-start,-end), names_to = "comparison", values_to = "pvalue") %>%
  mutate(comparison = str_remove(comparison,"pvalue_" ))

basemean_long <- res_df_h3k4me3 %>%
  dplyr::select(seqnames,start,end, starts_with("baseMean")) %>%
  pivot_longer(c(-seqnames,-start,-end), names_to = "comparison", values_to = "baseMean") %>%
  mutate(comparison = str_remove(comparison,"baseMean_" ))

# Join all tables to make a single long data frame
res_long_list <- list(basemean_long, log2fc_long, padj_long, pval_long)
res_long_h3k4me3_df <- purrr::reduce(res_long_list,left_join,by =c("seqnames","start","end","comparison") )

# Make long dataframe versions of log2fc, padj, etc
log2fc_long <- res_df_rnapolii %>%
  dplyr::select(seqnames,start,end, starts_with("log2FoldChange")) %>%
  pivot_longer(c(-seqnames,-start,-end), names_to = "comparison", values_to = "log2FoldChange") %>%
  mutate(comparison = str_remove(comparison,"log2FoldChange_" ))

padj_long <- res_df_rnapolii %>%
  dplyr::select(seqnames,start,end, starts_with("padj")) %>%
  pivot_longer(c(-seqnames,-start,-end), names_to = "comparison", values_to = "padj") %>%
  mutate(comparison = str_remove(comparison,"padj_" ))

pval_long <- res_df_rnapolii %>%
  dplyr::select(seqnames,start,end, starts_with("pvalue")) %>%
  pivot_longer(c(-seqnames,-start,-end), names_to = "comparison", values_to = "pvalue") %>%
  mutate(comparison = str_remove(comparison,"pvalue_" ))

basemean_long <- res_df_rnapolii %>%
  dplyr::select(seqnames,start,end, starts_with("baseMean")) %>%
  pivot_longer(c(-seqnames,-start,-end), names_to = "comparison", values_to = "baseMean") %>%
  mutate(comparison = str_remove(comparison,"baseMean_" ))


res_long_list <- list(basemean_long, log2fc_long, padj_long, pval_long)
res_long_rnapolii_df <- purrr::reduce(res_long_list,left_join,by =c("seqnames","start","end","comparison") )

saveRDS(res_long_h3k4me3_df, "data_output/2024-01_diffbind_cutnrun/res_long_h3k4me3_df.Rds")
saveRDS(res_long_rnapolii_df, "data_output/2024-01_diffbind_cutnrun/res_long_rnapolii_df.Rds")

######
### Heatmaps
######

### Build norm counts z score matrix
norm_cts_h3k4me3 <- counts(deseq_list$H3K4me3_seacrigg_peaks_union_deseq_counts_genotype_lrt_deseq_obj_fitted.Rds,norm =TRUE)
norm_cts_rnapolii <- counts(deseq_list$RNAPolII_seacrigg_peaks_union_deseq_counts_genotype_lrt_deseq_obj_fitted.Rds,norm =TRUE)

rownames(norm_cts_h3k4me3) <- as.data.frame(granges(rowRanges(deseq_list$H3K4me3_seacrigg_peaks_union_deseq_counts_genotype_lrt_deseq_obj_fitted.Rds))) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  pull(range_id)

rownames(norm_cts_rnapolii) <- as.data.frame(granges(rowRanges(deseq_list$RNAPolII_seacrigg_peaks_union_deseq_counts_genotype_lrt_deseq_obj_fitted.Rds))) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  pull(range_id)

saveRDS(norm_cts_h3k4me3, "data_output/2024-01_diffbind_cutnrun/norm_cts_h3k4me3.Rds")
saveRDS(norm_cts_rnapolii, "data_output/2024-01_diffbind_cutnrun/norm_cts_rnapolii.Rds")


# DE genes
de_genotype <- list()
de_genotype$gh_h3k4me3_id <- res_long_h3k4me3_df %>%
  dplyr::select(seqnames, start, end, comparison, padj, log2FoldChange) %>%
  filter(str_detect(comparison, "GH_genotype")) %>%
  filter(padj < 0.01) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  pull(range_id) %>%
  unique()

de_genotype$sh_h3k4me3_id <- res_long_h3k4me3_df %>%
  dplyr::select(seqnames, start, end, comparison, padj, log2FoldChange) %>%
  filter(str_detect(comparison, "SH_genotype")) %>%
  filter(padj < 0.01) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  pull(range_id) %>%
  unique()


de_genotype$gh_rnapolii_id <- res_long_rnapolii_df %>%
  dplyr::select(seqnames, start, end, comparison, padj, log2FoldChange) %>%
  filter(str_detect(comparison, "GH_genotype")) %>%
  filter(padj < 0.01) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  pull(range_id) %>%
  unique()

de_genotype$sh_rnapolii_id <- res_long_rnapolii_df %>%
  dplyr::select(seqnames, start, end, comparison, padj, log2FoldChange) %>%
  filter(str_detect(comparison, "SH_genotype")) %>%
  filter(padj < 0.01) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  pull(range_id) %>%
  unique()

norm_cts_h3k4me3_median <- get_median_per_group(norm_cts_h3k4me3)
norm_cts_rnapolii_median <- get_median_per_group(norm_cts_rnapolii)

saveRDS(de_genotype, "data_output/2024-01_diffbind_cutnrun/de_genotype.Rds")
saveRDS(norm_cts_h3k4me3_median, "data_output/2024-01_diffbind_cutnrun/norm_cts_h3k4me3_median.Rds")
saveRDS(norm_cts_rnapolii_median, "data_output/2024-01_diffbind_cutnrun/norm_cts_rnapolii_median.Rds")

### Build log2Fold change matrix
# Select only log2FoldChange variables

log2fc_mat_h3k4me3_genotype_gh <- res_df_h3k4me3 %>%
  filter(padj_H3K4me3_dsxmut_GH_genotype_ashr < 0.01 | 
           padj_H3K4me3_frumut_GH_genotype_ashr < 0.01 |
           padj_H3K4me3_Or67dmut_GH_genotype_ashr < 0.01 |
           padj_H3K4me3_Or47bmut_GH_genotype_ashr < 0.01) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  dplyr::select(range_id, starts_with("log2FoldChange") & contains("GH_genotype") & ! contains("lrt")) %>%
  column_to_rownames(var = "range_id") %>%
  as.matrix()

log2fc_mat_h3k4me3_genotype_sh <- res_df_h3k4me3 %>%
  filter(padj_H3K4me3_dsxmut_SH_genotype_ashr < 0.01 | 
           padj_H3K4me3_frumut_SH_genotype_ashr < 0.01 |
           padj_H3K4me3_Or67dmut_SH_genotype_ashr < 0.01 |
           padj_H3K4me3_Or47bmut_SH_genotype_ashr < 0.01) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  dplyr::select(range_id, starts_with("log2FoldChange") & contains("SH_genotype") & ! contains("lrt")) %>%
  column_to_rownames(var = "range_id") %>%
  as.matrix()


colnames(log2fc_mat_h3k4me3_genotype_gh) <- str_remove(colnames(log2fc_mat_h3k4me3_genotype_gh) ,"log2FoldChange_")
colnames(log2fc_mat_h3k4me3_genotype_gh) <- str_remove(colnames(log2fc_mat_h3k4me3_genotype_gh) ,"_genotype_ashr")

colnames(log2fc_mat_h3k4me3_genotype_sh) <- str_remove(colnames(log2fc_mat_h3k4me3_genotype_sh) ,"log2FoldChange_")
colnames(log2fc_mat_h3k4me3_genotype_sh) <- str_remove(colnames(log2fc_mat_h3k4me3_genotype_sh) ,"_genotype_ashr")

log2fc_mat_rnapolii_genotype_gh <- res_df_rnapolii %>%
  filter(padj_RNAPolII_dsxmut_GH_genotype_ashr < 0.01 | 
           padj_RNAPolII_frumut_GH_genotype_ashr < 0.01 |
           padj_RNAPolII_Or67dmut_GH_genotype_ashr < 0.01 |
           padj_RNAPolII_Or47bmut_GH_genotype_ashr < 0.01) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  dplyr::select(range_id, starts_with("log2FoldChange") & contains("GH_genotype") & ! contains("lrt")) %>%
  column_to_rownames(var = "range_id") %>%
  as.matrix()

log2fc_mat_rnapolii_genotype_sh <- res_df_rnapolii %>%
  filter(padj_RNAPolII_dsxmut_SH_genotype_ashr < 0.01 | 
           padj_RNAPolII_frumut_SH_genotype_ashr < 0.01 |
           padj_RNAPolII_Or67dmut_SH_genotype_ashr < 0.01 |
           padj_RNAPolII_Or47bmut_SH_genotype_ashr < 0.01) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  dplyr::select(range_id, starts_with("log2FoldChange") & contains("SH_genotype") & ! contains("lrt")) %>%
  column_to_rownames(var = "range_id") %>%
  as.matrix()


colnames(log2fc_mat_rnapolii_genotype_gh) <- str_remove(colnames(log2fc_mat_rnapolii_genotype_gh) ,"log2FoldChange_")
colnames(log2fc_mat_rnapolii_genotype_gh) <- str_remove(colnames(log2fc_mat_rnapolii_genotype_gh) ,"_genotype_ashr")

colnames(log2fc_mat_rnapolii_genotype_sh) <- str_remove(colnames(log2fc_mat_rnapolii_genotype_sh) ,"log2FoldChange_")
colnames(log2fc_mat_rnapolii_genotype_sh) <- str_remove(colnames(log2fc_mat_rnapolii_genotype_sh) ,"_genotype_ashr")

saveRDS(log2fc_mat_h3k4me3_genotype_gh, "data_output/2024-01_diffbind_cutnrun/log2fc_mat_h3k4me3_genotype_gh.Rds")
saveRDS(log2fc_mat_h3k4me3_genotype_sh, "data_output/2024-01_diffbind_cutnrun/log2fc_mat_h3k4me3_genotype_sh.Rds")
saveRDS(log2fc_mat_rnapolii_genotype_gh, "data_output/2024-01_diffbind_cutnrun/log2fc_mat_rnapolii_genotype_gh.Rds")
saveRDS(log2fc_mat_rnapolii_genotype_sh, "data_output/2024-01_diffbind_cutnrun/log2fc_mat_rnapolii_genotype_sh.Rds")

######
# Connect to genes
######

# Make ranges for both peaks objects
res_gr_h3k4me3 <- makeGRangesFromDataFrame( res_df_h3k4me3,keep.extra.columns = TRUE)
res_gr_rnapolii <- makeGRangesFromDataFrame( res_df_rnapolii,keep.extra.columns = TRUE)

# Get promoters info
promoters_dm6 <- promoters(TxDb.Dmelanogaster.UCSC.dm6.ensGene, columns = "gene_id",upstream=0,downstream=1)
seqlevels(promoters_dm6) <- str_remove(seqlevels(promoters_dm6),"chr")

# Map gene names
promoter_map_df <- clusterProfiler::bitr(geneID = unlist(promoters_dm6$gene_id), fromType = "ENSEMBL",toType ="SYMBOL",OrgDb = org.Dm.eg.db )
promoter_map <- promoter_map_df$SYMBOL 
names(promoter_map) <-  promoter_map_df$ENSEMBL
promoters_dm6$gene_name <- promoter_map[unlist(promoters_dm6$gene_id)]

# Map entrez id
promoter_map_df_entrez <- clusterProfiler::bitr(geneID = unlist(promoters_dm6$gene_id), fromType = "ENSEMBL",toType ="ENTREZID",OrgDb = org.Dm.eg.db )
promoter_map_entrez <- promoter_map_df_entrez$ENTREZ 
names(promoter_map_entrez) <-  promoter_map_df_entrez$ENSEMBL
promoters_dm6$entrez_id <- promoter_map_entrez[unlist(promoters_dm6$gene_id)]

# Join H3K4me3 peaks and promoters data by the nearest TSS (keep distance info)
res_gr_h3k4me3_genes <- plyranges::join_nearest_left(x = res_gr_h3k4me3,y = promoters_dm6, distance = TRUE) 

res_gr_h3k4me3_genes$range_id <- as.data.frame(granges(res_gr_h3k4me3_genes)) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  pull(range_id)

names(res_gr_h3k4me3_genes) <- as.data.frame(elementMetadata(res_gr_h3k4me3_genes)) %>%
  pull(range_id)

# Join RNAPolII peaks and promoters data by the nearest (keep distance info)
res_gr_rnapolii_genes <- plyranges::join_nearest_left(x = res_gr_rnapolii,y = promoters_dm6, distance = TRUE)
res_gr_rnapolii_genes$range_id <- as.data.frame(granges(res_gr_rnapolii_genes)) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  pull(range_id)

names(res_gr_rnapolii_genes) <- as.data.frame(elementMetadata(res_gr_rnapolii_genes)) %>%
  pull(range_id)

saveRDS(res_gr_h3k4me3_genes, "data_output/2024-01_diffbind_cutnrun/res_gr_h3k4me3_genes.Rds")
saveRDS(res_gr_rnapolii_genes, "data_output/2024-01_diffbind_cutnrun/res_gr_rnapolii_genes.Rds")

# Get differential peaks ranges
diff_ranges <- list()
diff_ranges$h3k4me3_gh <- as.data.frame(elementMetadata(res_gr_h3k4me3_genes)) %>%
  filter(padj_H3K4me3_dsxmut_GH_genotype_ashr < 0.01 | 
           padj_H3K4me3_frumut_GH_genotype_ashr < 0.01 |
           padj_H3K4me3_Or67dmut_GH_genotype_ashr < 0.01 |
           padj_H3K4me3_Or47bmut_GH_genotype_ashr < 0.01) %>%
  pull(range_id)


diff_ranges$h3k4me3_sh <- as.data.frame(elementMetadata(res_gr_h3k4me3_genes)) %>%
  filter(padj_H3K4me3_dsxmut_SH_genotype_ashr < 0.01 | 
           padj_H3K4me3_frumut_SH_genotype_ashr < 0.01 |
           padj_H3K4me3_Or67dmut_SH_genotype_ashr < 0.01 |
           padj_H3K4me3_Or47bmut_SH_genotype_ashr < 0.01) %>%
  pull(range_id)


diff_ranges$rnapolii_gh <- as.data.frame(elementMetadata(res_gr_rnapolii_genes)) %>%
  filter(padj_RNAPolII_dsxmut_GH_genotype_ashr < 0.01 | 
           padj_RNAPolII_frumut_GH_genotype_ashr < 0.01 |
           padj_RNAPolII_Or67dmut_GH_genotype_ashr < 0.01 |
           padj_RNAPolII_Or47bmut_GH_genotype_ashr < 0.01) %>%
  pull(range_id)


diff_ranges$rnapolii_sh <- as.data.frame(elementMetadata(res_gr_rnapolii_genes)) %>%
  filter(padj_RNAPolII_dsxmut_SH_genotype_ashr < 0.01 | 
           padj_RNAPolII_frumut_SH_genotype_ashr < 0.01 |
           padj_RNAPolII_Or67dmut_SH_genotype_ashr < 0.01 |
           padj_RNAPolII_Or47bmut_SH_genotype_ashr < 0.01) %>%
  pull(range_id)

diff_ranges$h3k4me3_gh_genes <- as.data.frame(elementMetadata(res_gr_h3k4me3_genes)) %>%
  filter(padj_H3K4me3_dsxmut_GH_genotype_ashr < 0.01 | 
           padj_H3K4me3_frumut_GH_genotype_ashr < 0.01 |
           padj_H3K4me3_Or67dmut_GH_genotype_ashr < 0.01 |
           padj_H3K4me3_Or47bmut_GH_genotype_ashr < 0.01) %>%
  pull(gene_name) %>%
  unique() %>%
  na.omit()


diff_ranges$h3k4me3_sh_genes <- as.data.frame(elementMetadata(res_gr_h3k4me3_genes)) %>%
  filter(padj_H3K4me3_dsxmut_SH_genotype_ashr < 0.01 | 
           padj_H3K4me3_frumut_SH_genotype_ashr < 0.01 |
           padj_H3K4me3_Or67dmut_SH_genotype_ashr < 0.01 |
           padj_H3K4me3_Or47bmut_SH_genotype_ashr < 0.01) %>%
  pull(gene_name)%>%
  unique() %>%
  na.omit()


diff_ranges$rnapolii_gh_genes <- as.data.frame(elementMetadata(res_gr_rnapolii_genes)) %>%
  filter(padj_RNAPolII_dsxmut_GH_genotype_ashr < 0.01 | 
           padj_RNAPolII_frumut_GH_genotype_ashr < 0.01 |
           padj_RNAPolII_Or67dmut_GH_genotype_ashr < 0.01 |
           padj_RNAPolII_Or47bmut_GH_genotype_ashr < 0.01) %>%
  pull(gene_name)%>%
  unique() %>%
  na.omit()


diff_ranges$rnapolii_sh_genes <- as.data.frame(elementMetadata(res_gr_rnapolii_genes)) %>%
  filter(padj_RNAPolII_dsxmut_SH_genotype_ashr < 0.01 | 
           padj_RNAPolII_frumut_SH_genotype_ashr < 0.01 |
           padj_RNAPolII_Or67dmut_SH_genotype_ashr < 0.01 |
           padj_RNAPolII_Or47bmut_SH_genotype_ashr < 0.01) %>%
  pull(gene_name)%>%
  unique() %>%
  na.omit()

saveRDS(diff_ranges, "data_output/2024-01_diffbind_cutnrun/diff_ranges.Rds")


### Read target genes

#### Fru

fru_ma_genes <- read_csv("data_input/gene_annotations/fru_targets/Sheet 1-fruma_targets.csv")[,-2]
fru_mc_genes <- read_csv("data_input/gene_annotations/fru_targets/Sheet 2-frumc_targets.csv")[,-2]
fru_mb_genes <- read_csv("data_input/gene_annotations/fru_targets/Sheet 3-frumb_targets.csv")[,-2]


fru_ma_genes <- fru_ma_genes$FruMA
fru_mc_genes <- fru_mc_genes$FruMC
fru_mb_genes <- fru_mb_genes$FruMB

fru_target_genes <- unique(c(fru_ma_genes, fru_mc_genes, fru_mb_genes))

de_genotype$h3k4me3_gh_genes_id_fru <- res_gr_h3k4me3_genes[diff_ranges$h3k4me3_gh]$range_id[res_gr_h3k4me3_genes[diff_ranges$h3k4me3_gh]$gene_name %in% fru_target_genes]

de_genotype$h3k4me3_sh_genes_id_fru <- res_gr_h3k4me3_genes[diff_ranges$h3k4me3_sh]$range_id[res_gr_h3k4me3_genes[diff_ranges$h3k4me3_sh]$gene_name %in% fru_target_genes]


de_genotype$rnapolii_gh_genes_id_fru <- res_gr_rnapolii_genes[diff_ranges$rnapolii_gh]$range_id[res_gr_rnapolii_genes[diff_ranges$rnapolii_gh]$gene_name %in% fru_target_genes]

de_genotype$rnapolii_sh_genes_id_fru <- res_gr_rnapolii_genes[diff_ranges$rnapolii_sh]$range_id[res_gr_rnapolii_genes[diff_ranges$rnapolii_sh]$gene_name %in% fru_target_genes]


#### dsx
dsx_target_genes <- readxl::read_excel("data_input/gene_annotations/dsx_target_genes.xlsx",sheet = 3)
dsx_target_genes <- dsx_target_genes$`Male_Fatbody_DamIDseq_Gene_Level_Occupancy> 0.9`


de_genotype_h3k4me3_gh_genes_id_dsx <- res_gr_h3k4me3_genes[diff_ranges$h3k4me3_gh]$range_id[res_gr_h3k4me3_genes[diff_ranges$h3k4me3_gh]$gene_name %in% dsx_target_genes]

de_genotype_h3k4me3_sh_genes_id_dsx <- res_gr_h3k4me3_genes[diff_ranges$h3k4me3_sh]$range_id[res_gr_h3k4me3_genes[diff_ranges$h3k4me3_sh]$gene_name %in% dsx_target_genes]


de_genotype_rnapolii_gh_genes_id_dsx <- res_gr_rnapolii_genes[diff_ranges$rnapolii_gh]$range_id[res_gr_rnapolii_genes[diff_ranges$rnapolii_gh]$gene_name %in% dsx_target_genes]

de_genotype_rnapolii_sh_genes_id_dsx <- res_gr_rnapolii_genes[diff_ranges$rnapolii_sh]$range_id[res_gr_rnapolii_genes[diff_ranges$rnapolii_sh]$gene_name %in% dsx_target_genes]


# Make gene name to range map 
res_gr_h3k4me3_genes <- unique(res_gr_h3k4me3_genes)
res_gr_rnapolii_genes <- unique(res_gr_rnapolii_genes)

gene_map_h3k4me3 <- res_gr_h3k4me3_genes$gene_name
names(gene_map_h3k4me3) <- res_gr_h3k4me3_genes$range_id

gene_map_rnapolii <- res_gr_rnapolii_genes$gene_name
names(gene_map_rnapolii) <- res_gr_rnapolii_genes$range_id

# Filter matrices and rename to have them gene-based instead of peak-based 
norm_cts_h3k4me3_gene <- rownames_filt_and_convert_to_gene(norm_cts_h3k4me3,gene_map_h3k4me3)
norm_cts_h3k4me3_median_gene <- rownames_filt_and_convert_to_gene(norm_cts_h3k4me3_median,gene_map_h3k4me3)
norm_cts_rnapolii_gene <- rownames_filt_and_convert_to_gene(norm_cts_rnapolii,gene_map_rnapolii)
norm_cts_rnapolii_median_gene <- rownames_filt_and_convert_to_gene(norm_cts_rnapolii_median,gene_map_rnapolii)

log2fc_mat_h3k4me3_genotype_gh_gene <- rownames_filt_and_convert_to_gene(log2fc_mat_h3k4me3_genotype_gh,gene_map_h3k4me3)
log2fc_mat_rnapolii_genotype_gh_gene <- rownames_filt_and_convert_to_gene(log2fc_mat_rnapolii_genotype_gh,gene_map_rnapolii)

log2fc_mat_h3k4me3_genotype_sh_gene <- rownames_filt_and_convert_to_gene(log2fc_mat_h3k4me3_genotype_sh,gene_map_h3k4me3)
log2fc_mat_rnapolii_genotype_sh_gene <- rownames_filt_and_convert_to_gene(log2fc_mat_rnapolii_genotype_sh,gene_map_rnapolii)


saveRDS(norm_cts_h3k4me3_gene, "data_output/2024-01_diffbind_cutnrun/norm_cts_h3k4me3_gene.Rds")
saveRDS(norm_cts_h3k4me3_median_gene, "data_output/2024-01_diffbind_cutnrun/norm_cts_h3k4me3_median_gene.Rds")
saveRDS(norm_cts_rnapolii_gene, "data_output/2024-01_diffbind_cutnrun/norm_cts_rnapolii_gene.Rds")
saveRDS(norm_cts_rnapolii_median_gene, "data_output/2024-01_diffbind_cutnrun/norm_cts_rnapolii_median_gene.Rds")
saveRDS(log2fc_mat_h3k4me3_genotype_gh_gene, "data_output/2024-01_diffbind_cutnrun/log2fc_mat_h3k4me3_genotype_gh_gene.Rds")
saveRDS(log2fc_mat_rnapolii_genotype_gh_gene, "data_output/2024-01_diffbind_cutnrun/log2fc_mat_rnapolii_genotype_gh_gene.Rds")
saveRDS(log2fc_mat_h3k4me3_genotype_sh_gene, "data_output/2024-01_diffbind_cutnrun/log2fc_mat_h3k4me3_genotype_sh_gene.Rds")
saveRDS(log2fc_mat_rnapolii_genotype_sh_gene, "data_output/2024-01_diffbind_cutnrun/log2fc_mat_rnapolii_genotype_sh_gene.Rds")



## Make lists for different combinations of comparisons
res_long_h3k4me3_df_any <- res_long_h3k4me3_df %>%
  drop_na(padj) %>%
  filter(padj < 0.01)

res_long_rnapolii_df_any <- res_long_rnapolii_df %>%
  drop_na(padj) %>%
  filter(padj < 0.01)


## GH: all to CS

### All
comp_list <- list()

comp_list$genotype_gh_all_h3k4me3_rangeid <- res_long_h3k4me3_df_any %>%
  filter(str_detect(comparison,"GH_genotype_ashr")) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  dplyr::select(range_id, comparison) %>%
  mutate(comparison = str_remove(comparison, "_GH_genotype_ashr")) %>%
  split(.$comparison) %>%
  map("range_id")

comp_list$genotype_gh_all_rnapolii_rangeid <- res_long_rnapolii_df_any %>%
  filter(str_detect(comparison,"GH_genotype_ashr")) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  dplyr::select(range_id, comparison) %>%
  mutate(comparison = str_remove(comparison, "_GH_genotype_ashr")) %>%
  split(.$comparison) %>%
  map("range_id")

### Up
comp_list$genotype_gh_up_h3k4me3_rangeid <- res_long_h3k4me3_df_any %>%
  filter(str_detect(comparison,"GH_genotype_ashr")) %>%
  filter(log2FoldChange > 0) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  dplyr::select(range_id, comparison) %>%
  mutate(comparison = str_remove(comparison, "_GH_genotype_ashr")) %>%
  split(.$comparison) %>%
  map("range_id")

comp_list$genotype_gh_up_rnapolii_rangeid <- res_long_rnapolii_df_any %>%
  filter(str_detect(comparison,"GH_genotype_ashr")) %>%
  filter(log2FoldChange > 0) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  dplyr::select(range_id, comparison) %>%
  mutate(comparison = str_remove(comparison, "_GH_genotype_ashr")) %>%
  split(.$comparison) %>%
  map("range_id")

### Down
comp_list$genotype_gh_down_h3k4me3_rangeid <- res_long_h3k4me3_df_any %>%
  filter(str_detect(comparison,"GH_genotype_ashr")) %>%
  filter(log2FoldChange < 0) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  dplyr::select(range_id, comparison) %>%
  mutate(comparison = str_remove(comparison, "_GH_genotype_ashr")) %>%
  split(.$comparison) %>%
  map("range_id")

comp_list$genotype_gh_down_rnapolii_rangeid <- res_long_rnapolii_df_any %>%
  filter(str_detect(comparison,"GH_genotype_ashr")) %>%
  filter(log2FoldChange < 0) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  dplyr::select(range_id, comparison) %>%
  mutate(comparison = str_remove(comparison, "_GH_genotype_ashr")) %>%
  split(.$comparison) %>%
  map("range_id")

## SH: all to CS

### All
comp_list$genotype_sh_all_h3k4me3_rangeid <- res_long_h3k4me3_df_any %>%
  filter(str_detect(comparison,"SH_genotype_ashr")) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  dplyr::select(range_id, comparison) %>%
  mutate(comparison = str_remove(comparison, "_SH_genotype_ashr")) %>%
  split(.$comparison) %>%
  map("range_id")

comp_list$genotype_sh_all_rnapolii_rangeid <- res_long_rnapolii_df_any %>%
  filter(str_detect(comparison,"SH_genotype_ashr")) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  dplyr::select(range_id, comparison) %>%
  mutate(comparison = str_remove(comparison, "_SH_genotype_ashr")) %>%
  split(.$comparison) %>%
  map("range_id")


### Up
comp_list$genotype_sh_up_h3k4me3_rangeid <- res_long_h3k4me3_df_any %>%
  filter(str_detect(comparison,"SH_genotype_ashr")) %>%
  filter(log2FoldChange > 0) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  dplyr::select(range_id, comparison) %>%
  mutate(comparison = str_remove(comparison, "_SH_genotype_ashr")) %>%
  split(.$comparison) %>%
  map("range_id")

comp_list$genotype_sh_up_rnapolii_rangeid <- res_long_rnapolii_df_any %>%
  filter(str_detect(comparison,"SH_genotype_ashr")) %>%
  filter(log2FoldChange > 0) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  dplyr::select(range_id, comparison) %>%
  mutate(comparison = str_remove(comparison, "_SH_genotype_ashr")) %>%
  split(.$comparison) %>%
  map("range_id")

### Down
comp_list$genotype_sh_down_h3k4me3_rangeid <- res_long_h3k4me3_df_any %>%
  filter(str_detect(comparison,"SH_genotype_ashr")) %>%
  filter(log2FoldChange < 0) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  dplyr::select(range_id, comparison) %>%
  mutate(comparison = str_remove(comparison, "_SH_genotype_ashr")) %>%
  split(.$comparison) %>%
  map("range_id")

comp_list$genotype_sh_down_rnapolii_rangeid <- res_long_rnapolii_df_any %>%
  filter(str_detect(comparison,"SH_genotype_ashr")) %>%
  filter(log2FoldChange < 0) %>%
  mutate(range_id = paste0(seqnames,":",start,"-",end)) %>%
  dplyr::select(range_id, comparison) %>%
  mutate(comparison = str_remove(comparison, "_SH_genotype_ashr")) %>%
  split(.$comparison) %>%
  map("range_id")

### Export to file
saveRDS(comp_list, "data_output/2024-01_diffbind_cutnrun/comp_list.Rds")


######
### Significant LRT genotype regions
######

# Filter results dataframe to get significant lrt genotype regions
# Select relevant columns for pairwise LogFC
pairwise_log2fc_h3k4me3_df <- res_df_h3k4me3 %>%
  filter(padj_H3K4me3_genotype_lrt < 0.01) %>%
  select(seqnames, start, end, sign_paste_genotype_gh,starts_with("log2FoldChange")) %>%
  select(!contains("lrt")) %>%
  select(!contains("housing")) %>%
  select(!contains("_SH_"))

# Filter by absolute Fold change
# Create a matrix with only fold change 
pairwise_log2fc_h3k4me3_mat <- pairwise_log2fc_h3k4me3_df %>%
  select(starts_with("log2FoldChange"))
# Make a logical filterl
logfc_filter <- rowSums(abs(pairwise_log2fc_h3k4me3_mat) > 1) >= 1

# Filter out the dataframe
pairwise_log2fc_h3k4me3_df <- pairwise_log2fc_h3k4me3_df %>%
  filter(logfc_filter)

# Change the column names for the data.frame to keeep only relevant info 
colnames(pairwise_log2fc_h3k4me3_df) <- str_remove(colnames(pairwise_log2fc_h3k4me3_df), "log2FoldChange_")
colnames(pairwise_log2fc_h3k4me3_df) <- str_remove(colnames(pairwise_log2fc_h3k4me3_df), "_genotype_ashr")
colnames(pairwise_log2fc_h3k4me3_df) <- str_remove(colnames(pairwise_log2fc_h3k4me3_df), "_GH")
colnames(pairwise_log2fc_h3k4me3_df) <- str_remove(colnames(pairwise_log2fc_h3k4me3_df), "H3K4me3_") 
colnames(pairwise_log2fc_h3k4me3_df) <- str_remove(colnames(pairwise_log2fc_h3k4me3_df), "mut") 


pairwise_log2fc_h3k4me3_df <- pairwise_log2fc_h3k4me3_df %>%
  pivot_longer(cols = - c(seqnames, start, end,sign_paste_genotype_gh), names_to = "genotype",values_to = "log2FoldChange") %>%
  mutate(gene_name = paste0(seqnames,":",start,"-",end))


# Filter results dataframe to get significant lrt genotype regions
# Select relevant columns for pairwise LogFC
pairwise_log2fc_rnapolii_df <- res_df_rnapolii %>%
  filter(padj_RNAPolII_genotype_lrt < 0.01) %>%
  select(seqnames, start, end, sign_paste_genotype_gh,starts_with("log2FoldChange")) %>%
  select(!contains("lrt")) %>%
  select(!contains("housing")) %>%
  select(!contains("_SH_"))

# Filter by absolute Fold change
# Create a matrix with only fold change 
pairwise_log2fc_rnapolii_mat <- pairwise_log2fc_rnapolii_df %>%
  select(starts_with("log2FoldChange"))
# Make a logical filterl
logfc_filter <- rowSums(abs(pairwise_log2fc_rnapolii_mat) > 1) >= 1

# Filter out the dataframe
pairwise_log2fc_rnapolii_df <- pairwise_log2fc_rnapolii_df %>%
  filter(logfc_filter)

# Change the column names for the data.frame to keeep only relevant info 
colnames(pairwise_log2fc_rnapolii_df) <- str_remove(colnames(pairwise_log2fc_rnapolii_df), "log2FoldChange_")
colnames(pairwise_log2fc_rnapolii_df) <- str_remove(colnames(pairwise_log2fc_rnapolii_df), "_genotype_ashr")
colnames(pairwise_log2fc_rnapolii_df) <- str_remove(colnames(pairwise_log2fc_rnapolii_df), "_GH")
colnames(pairwise_log2fc_rnapolii_df) <- str_remove(colnames(pairwise_log2fc_rnapolii_df), "RNAPolII_") 
colnames(pairwise_log2fc_rnapolii_df) <- str_remove(colnames(pairwise_log2fc_rnapolii_df), "mut") 


pairwise_log2fc_rnapolii_df <- pairwise_log2fc_rnapolii_df %>%
  pivot_longer(cols = - c(seqnames, start, end,sign_paste_genotype_gh), names_to = "genotype",values_to = "log2FoldChange") %>%
  mutate(gene_name = paste0(seqnames,":",start,"-",end))

# Save data
saveRDS(pairwise_log2fc_h3k4me3_df, "data_output/2024-01_diffbind_cutnrun/pairwise_log2fc_h3k4me3_df.Rds")
saveRDS(pairwise_log2fc_rnapolii_df, "data_output/2024-01_diffbind_cutnrun/pairwise_log2fc_rnapolii_df.Rds")


