######
### Libraries
######
library(tidyverse)
library(clusterProfiler)
library(org.Dm.eg.db)
library(ReactomePA)

######
### Read data
######

res_gr_h3k4me3_genes <- readRDS( "data_output/2024-01_diffbind_cutnrun/res_gr_h3k4me3_genes.Rds")
res_gr_rnapolii_genes <- readRDS( "data_output/2024-01_diffbind_cutnrun/res_gr_rnapolii_genes.Rds")


######
### reactome terms
######

res_gr_h3k4me3_genes_list <- as.data.frame(unique(res_gr_h3k4me3_genes)) %>%
  dplyr::select(ENTREZID, contains("genotype") & starts_with("padj") & ! contains("lrt")) %>%
  pivot_longer(-ENTREZID, names_to ="comparison",values_to = "padj") %>%
  filter(padj < 0.01) %>%
  mutate(comparison = str_remove(comparison, "padj_H3K4me3_")) %>%
  mutate(comparison = str_remove(comparison, "_genotype_ashr") ) %>%
  drop_na() %>%
  split(f = .$comparison) %>%
  map(~unique(.x$ENTREZID))

res_gr_h3k4me3_genes_reactome <- compareCluster(geneClusters = res_gr_h3k4me3_genes_list,
                                                organism = "fly",
                                                fun = enrichPathway,
                                                qvalueCutoff = 0.05,
                                                readable = TRUE)




res_gr_rnapolii_genes_list <- as.data.frame(unique(res_gr_rnapolii_genes)) %>%
  dplyr::select(ENTREZID, contains("genotype")  & starts_with("padj") & ! contains("lrt")) %>%
  pivot_longer(-ENTREZID, names_to ="comparison",values_to = "padj") %>%
  filter(padj < 0.01) %>%
  mutate(comparison = str_remove(comparison, "padj_RNAPolII_")) %>%
  mutate(comparison = str_remove(comparison, "_genotype_ashr") ) %>%
  drop_na() %>%
  split(f = .$comparison) %>%
  map(~unique(.x$ENTREZID))


res_gr_rnapolii_genes_reactome <- compareCluster(geneClusters = res_gr_rnapolii_genes_list,
                                                 organism = "fly",
                                                 fun = enrichPathway,
                                                 qvalueCutoff = 0.05,
                                                 readable = TRUE)



saveRDS(res_gr_h3k4me3_genes_reactome,file = "data_output/2024-01_diffbind_cutnrun/res_gr_h3k4me3_genes_reactome.Rds")
saveRDS(res_gr_rnapolii_genes_reactome,file = "data_output/2024-01_diffbind_cutnrun/res_gr_rnapolii_genes_reactome.Rds")

