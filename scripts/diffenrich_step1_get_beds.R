library(tidyverse)
library(rtracklayer)
library(org.Dm.eg.db)
library(GenomicFeatures)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

# Inputs
diff_expressed_table_file <- "data_output/2024-1_diffexpr_output/res_df_all_join.tsv"
res_df <- read_tsv(diff_expressed_table_file)



# Get TSS
promoters_dm6 <- promoters(TxDb.Dmelanogaster.UCSC.dm6.ensGene,use.names = TRUE,upstream = 0,downstream = 1, columns = "gene_id")

# Map ids
gene_id_map_df <- clusterProfiler::bitr(fromType = "ENSEMBL",toType = "SYMBOL",geneID = promoters_dm6$gene_id,OrgDb = org.Dm.eg.db)
gene_id_map <- gene_id_map_df$SYMBOL
names(gene_id_map) <- gene_id_map_df$ENSEMBL

promoters_dm6 <- promoters_dm6[unlist(promoters_dm6$gene_id %in% gene_id_map_df$ENSEMBL)]
promoters_dm6$gene_name <- gene_id_map[unlist(promoters_dm6$gene_id)]
names(promoters_dm6) <- promoters_dm6$gene_name

seqlevels(promoters_dm6) <- str_remove(seqlevels(promoters_dm6),"chr")

res_df <- res_df[res_df$gene_name %in% names(promoters_dm6),]
promoters_dm6 <- promoters_dm6[names(promoters_dm6) %in% res_df$gene_name]


# All promoters
# Random selection of promoters
# Remove the sample part if you want all genes 
export.bed(promoters_dm6[sample(1:length(promoters_dm6), 1000)], con = "data_output/2024-1_diffexpr_output/all_genes.bed")


# Selection of promoters of specific comparisons 

# CS housing SH / GH up
CS_housing_up <- res_df %>%
  filter(padj_CS_housing < 0.05, log2FoldChange_CS_housing > 0.1) %>%
  pull(gene_name)

# Export as bed file
export.bed(promoters_dm6[CS_housing_up], con = "data_output/2024-1_diffexpr_output/CS_housing_up.bed")

# CS housing SH / GH down
CS_housing_down <- res_df %>%
  filter(padj_CS_housing < 0.05, log2FoldChange_CS_housing < -0.1) %>%
  pull(gene_name)
export.bed(promoters_dm6[CS_housing_down], con = "data_output/2024-1_diffexpr_output/CS_housing_down.bed")

# Mutant GH comparisons

# mutant GH / CantonS GH
DSX_GH_genotype_up <- res_df %>%
  filter(padj_DSX_GH_genotype < 0.05, log2FoldChange_DSX_GH_genotype > 0.1) %>%
  pull(gene_name)

export.bed(promoters_dm6[DSX_GH_genotype_up], con = "data_output/2024-1_diffexpr_output/DSX_GH_genotype_up.bed")

DSX_GH_genotype_down <- res_df %>%
  filter(padj_DSX_GH_genotype < 0.05, log2FoldChange_DSX_GH_genotype < -0.1) %>%
  pull(gene_name)

export.bed(promoters_dm6[DSX_GH_genotype_down], con = "data_output/2024-1_diffexpr_output/DSX_GH_genotype_down.bed")

FRU_GH_genotype_up <- res_df %>%
  filter(padj_FRU_GH_genotype < 0.05, log2FoldChange_FRU_GH_genotype > 0.1) %>%
  pull(gene_name)

export.bed(promoters_dm6[FRU_GH_genotype_up], con = "data_output/2024-1_diffexpr_output/FRU_GH_genotype_up.bed")

FRU_GH_genotype_down <- res_df %>%
  filter(padj_FRU_GH_genotype < 0.05, log2FoldChange_FRU_GH_genotype < -0.1) %>%
  pull(gene_name)

export.bed(promoters_dm6[FRU_GH_genotype_down], con = "data_output/2024-1_diffexpr_output/FRU_GH_genotype_down.bed")

Or47b_GH_genotype_up <- res_df %>%
  filter(padj_Or47b_GH_genotype < 0.05, log2FoldChange_Or47b_GH_genotype > 0.1) %>%
  pull(gene_name)

export.bed(promoters_dm6[Or47b_GH_genotype_up], con = "data_output/2024-1_diffexpr_output/Or47b_GH_genotype_up.bed")

Or47b_GH_genotype_down <- res_df %>%
  filter(padj_Or47b_GH_genotype < 0.05, log2FoldChange_Or47b_GH_genotype < -0.1) %>%
  pull(gene_name)

export.bed(promoters_dm6[Or47b_GH_genotype_down], con = "data_output/2024-1_diffexpr_output/Or47b_GH_genotype_down.bed")

Or67d_GH_genotype_up <- res_df %>%
  filter(padj_Or67d_GH_genotype < 0.05, log2FoldChange_Or67d_GH_genotype > 0.1) %>%
  pull(gene_name)

export.bed(promoters_dm6[Or67d_GH_genotype_up], con = "data_output/2024-1_diffexpr_output/Or67d_GH_genotype_up.bed")

Or67d_GH_genotype_down <- res_df %>%
  filter(padj_Or67d_GH_genotype < 0.05, log2FoldChange_Or67d_GH_genotype < -0.1) %>%
  pull(gene_name)

export.bed(promoters_dm6[Or67d_GH_genotype_down], con = "data_output/2024-1_diffexpr_output/Or67d_GH_genotype_down.bed")

# mutant SH / CantonS SH
DSX_SH_genotype_up <- res_df %>%
  filter(padj_DSX_SH_genotype < 0.05, log2FoldChange_DSX_SH_genotype > 0.1) %>%
  pull(gene_name)

export.bed(promoters_dm6[DSX_SH_genotype_up], con = "data_output/2024-1_diffexpr_output/DSX_SH_genotype_up.bed")

DSX_SH_genotype_down <- res_df %>%
  filter(padj_DSX_SH_genotype < 0.05, log2FoldChange_DSX_SH_genotype < -0.1) %>%
  pull(gene_name)

export.bed(promoters_dm6[DSX_SH_genotype_down], con = "data_output/2024-1_diffexpr_output/DSX_SH_genotype_down.bed")

FRU_SH_genotype_up <- res_df %>%
  filter(padj_FRU_SH_genotype < 0.05, log2FoldChange_FRU_SH_genotype > 0.1) %>%
  pull(gene_name)

export.bed(promoters_dm6[FRU_SH_genotype_up], con = "data_output/2024-1_diffexpr_output/FRU_SH_genotype_up.bed")

FRU_SH_genotype_down <- res_df %>%
  filter(padj_FRU_SH_genotype < 0.05, log2FoldChange_FRU_SH_genotype < -0.1) %>%
  pull(gene_name)

export.bed(promoters_dm6[FRU_SH_genotype_down], con = "data_output/2024-1_diffexpr_output/FRU_SH_genotype_down.bed")

Or47b_SH_genotype_up <- res_df %>%
  filter(padj_Or47b_SH_genotype < 0.05, log2FoldChange_Or47b_SH_genotype > 0.1) %>%
  pull(gene_name)

export.bed(promoters_dm6[Or47b_SH_genotype_up], con = "data_output/2024-1_diffexpr_output/Or47b_SH_genotype_up.bed")

Or47b_SH_genotype_down <- res_df %>%
  filter(padj_Or47b_SH_genotype < 0.05, log2FoldChange_Or47b_SH_genotype < -0.1) %>%
  pull(gene_name)

export.bed(promoters_dm6[Or47b_SH_genotype_down], con = "data_output/2024-1_diffexpr_output/Or47b_SH_genotype_down.bed")

Or67d_SH_genotype_up <- res_df %>%
  filter(padj_Or67d_SH_genotype < 0.05, log2FoldChange_Or67d_SH_genotype > 0.1) %>%
  pull(gene_name)

export.bed(promoters_dm6[Or67d_SH_genotype_up], con = "data_output/2024-1_diffexpr_output/Or67d_SH_genotype_up.bed")

Or67d_SH_genotype_down <- res_df %>%
  filter(padj_Or67d_SH_genotype < 0.05, log2FoldChange_Or67d_SH_genotype < -0.1) %>%
  pull(gene_name)

export.bed(promoters_dm6[Or67d_SH_genotype_down], con = "data_output/2024-1_diffexpr_output/Or67d_SH_genotype_down.bed")

