library(tidyverse)
library(rtracklayer)
library(org.Dm.eg.db)
library(GenomicFeatures)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

# Inputs
diff_expressed_table_file = "~/repos/rnaseq_social_isolation/data_output/combat/res_housing_cs_combat.tsv"
bed_file_up <- "data_output/res_housing_cs_combat_up.bed"
bed_file_down <- "data_output/res_housing_cs_combat_down.bed"
bed_file_nochange <- "data_output/res_housing_cs_combat_nochange.bed"
bed_file_sort <- "data_output/res_housing_cs_combat_sort.bed"
# Read RNA
res_housing_cs_combat <- read_tsv("~/repos/rnaseq_social_isolation/data_output/combat/res_housing_cs_combat.tsv")

# Map ids
gene_id_map <- clusterProfiler::bitr(fromType = "ENSEMBL",toType = "SYMBOL",geneID = res_housing_cs_combat$geneid,OrgDb = org.Dm.eg.db)
res_housing_cs_combat <- res_housing_cs_combat[res_housing_cs_combat$geneid %in% gene_id_map$ENSEMBL,]
res_housing_cs_combat <- left_join(res_housing_cs_combat,gene_id_map,by = c("geneid" = "SYMBOL"))


# Get TSS
promoters_dm6 <- promoters(genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene),upstream = 1,downstream = 1)


res_housing_cs_combat <- res_housing_cs_combat[res_housing_cs_combat$geneid %in% names(promoters_dm6),]

# Genes up
genes_up <- res_housing_cs_combat %>%
  filter(padj < 0.05, log2FoldChange > 0.5)
# Genes down
genes_down <- res_housing_cs_combat %>%
  filter(padj < 0.05, log2FoldChange < - 0.5)

# Genes no change
genes_nochange <- res_housing_cs_combat %>%
  filter(padj > 0.5, abs(log2FoldChange) < 0.1) %>%
  slice_sample(n = 50)

seqlevels(promoters_dm6) <- str_remove(seqlevels(promoters_dm6),"chr")

# Split diff expressed genes tss to bed files

tss_up <- promoters_dm6[genes_up$geneid]
tss_down <- promoters_dm6[genes_down$geneid]
tss_nochange <- promoters_dm6[genes_nochange$geneid]


export.bed(tss_up,con = bed_file_up)
export.bed(tss_down, con = bed_file_down)
export.bed(tss_nochange, con = bed_file_nochange)

# Sort all genes by lfc
gene_id_sort_lfc <- res_housing_cs_combat %>%
  arrange(log2FoldChange) %>%
  pull(geneid)

tss_sort <- promoters_dm6[gene_id_sort_lfc]
export.bed(tss_sort, con = bed_file_sort)


