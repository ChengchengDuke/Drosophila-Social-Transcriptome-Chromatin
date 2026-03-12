######
### Libraries
######
library(tidyverse)
library(clusterProfiler)
library(org.Dm.eg.db)

######
### Read data
######

res_gr_h3k4me3_genes <- readRDS( "./data_output/2024-01_diffbind_cutnrun/res_gr_h3k4me3_genes.Rds")
res_gr_rnapolii_genes <- readRDS( "./data_output/2024-01_diffbind_cutnrun/res_gr_rnapolii_genes.Rds")


######
### KEGG terms
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

res_gr_h3k4me3_genes_KEGG <- compareCluster(geneClusters = res_gr_h3k4me3_genes_list,
                                            organism = "dme",
                                            keyType = "ncbi-geneid",
                                            fun = enrichKEGG,
                                          qvalueCutoff = 0.05)




res_gr_rnapolii_genes_list <- as.data.frame(unique(res_gr_rnapolii_genes)) %>%
  dplyr::select(ENTREZID, contains("genotype")  & starts_with("padj") & ! contains("lrt")) %>%
  pivot_longer(-ENTREZID, names_to ="comparison",values_to = "padj") %>%
  filter(padj < 0.01) %>%
  mutate(comparison = str_remove(comparison, "padj_RNAPolII_")) %>%
  mutate(comparison = str_remove(comparison, "_genotype_ashr") ) %>%
  drop_na() %>%
  split(f = .$comparison) %>%
  map(~unique(.x$ENTREZID))


res_gr_rnapolii_genes_KEGG <- compareCluster(geneClusters = res_gr_rnapolii_genes_list,
                                             organism = "dme",
                                             keyType = "ncbi-geneid",
                                             fun = enrichKEGG,
                                           qvalueCutoff = 0.05)


#change format
df_results <- res_gr_rnapolii_genes_KEGG@compareClusterResult

# Function to convert Entrez IDs to Gene Symbols while handling multiple IDs
convert_entrez_to_symbol <- function(entrez_ids) {
  genes <- unlist(strsplit(entrez_ids, "/"))  # Split multiple gene IDs
  gene_symbols <- mapIds(
    org.Dm.eg.db,
    keys = genes,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"  # Get first match if multiple exist
  )
  return(paste(na.omit(gene_symbols), collapse = "/"))  # Recombine into string
}

# Apply conversion function to each row
df_results$geneSymbol <- sapply(df_results$geneID, convert_entrez_to_symbol)
df_results$geneID <- df_results$geneSymbol
df_results$Description <- gsub("- Drosophila melanogaster (fruit fly)", "", df_results$Description, fixed = TRUE)

res_gr_rnapolii_genes_KEGG@compareClusterResult <- df_results




#saveRDS(res_gr_h3k4me3_genes_KEGG,file = "./data_output/2024-01_diffbind_cutnrun/res_gr_h3k4me3_genes_KEGG.Rds")
saveRDS(res_gr_rnapolii_genes_KEGG,file = "./data_output/2024-01_diffbind_cutnrun/res_gr_rnapolii_genes_KEGG.Rds")

