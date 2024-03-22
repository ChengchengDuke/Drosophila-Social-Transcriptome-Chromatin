######
### Libraries
######
library(tidyverse)
library(clusterProfiler)
library(org.Dm.eg.db)

######
### Read data
######

res_gr_h3k4me3_genes <- readRDS( "data_output/2024-01_diffbind_cutnrun/res_gr_h3k4me3_genes.Rds")
res_gr_rnapolii_genes <- readRDS( "data_output/2024-01_diffbind_cutnrun/res_gr_rnapolii_genes.Rds")


######
### GO terms
######

res_gr_h3k4me3_genes_list <- as.data.frame(unique(res_gr_h3k4me3_genes)) %>%
  dplyr::select(entrez_id, contains("genotype") & starts_with("padj") & ! contains("lrt")) %>%
  pivot_longer(-entrez_id, names_to ="comparison",values_to = "padj") %>%
  filter(padj < 0.01) %>%
  mutate(comparison = str_remove(comparison, "padj_H3K4me3_")) %>%
  mutate(comparison = str_remove(comparison, "_genotype_ashr") ) %>%
  drop_na() %>%
  split(f = .$comparison) %>%
  map(~unique(.x$entrez_id))

res_gr_h3k4me3_genes_go <- compareCluster(geneClusters = res_gr_h3k4me3_genes_list,
                                          OrgDb = org.Dm.eg.db,
                                          fun = enrichGO,
                                          qvalueCutoff = 0.05,
                                          readable = TRUE, ont = "ALL")


res_gr_rnapolii_genes_list <- as.data.frame(unique(res_gr_rnapolii_genes)) %>%
  dplyr::select(entrez_id, contains("genotype")  & starts_with("padj") & ! contains("lrt")) %>%
  pivot_longer(-entrez_id, names_to ="comparison",values_to = "padj") %>%
  filter(padj < 0.01) %>%
  mutate(comparison = str_remove(comparison, "padj_RNAPolII_")) %>%
  mutate(comparison = str_remove(comparison, "_genotype_ashr") ) %>%
  drop_na() %>%
  split(f = .$comparison) %>%
  map(~unique(.x$entrez_id))


res_gr_rnapolii_genes_go <- compareCluster(geneClusters = res_gr_rnapolii_genes_list,
                                           OrgDb = org.Dm.eg.db,
                                           fun = enrichGO,
                                           qvalueCutoff = 0.05,
                                           readable = TRUE, ont = "ALL")


saveRDS(res_gr_rnapolii_genes_go,file = "data_output/2024-1_diffexpr_output/res_gr_rnapolii_genes_go.Rds")
saveRDS(res_gr_h3k4me3_genes_go,file = "data_output/2024-1_diffexpr_output/res_gr_h3k4me3_genes_go.Rds")

