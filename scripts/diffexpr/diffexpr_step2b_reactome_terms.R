library(tidyverse)
library(clusterProfiler)
library(org.Dm.eg.db)
library(ReactomePA)
# Read data
comp_list <- readRDS(file = "data_output/2024-1_diffexpr_output/comp_list.Rds")

# Pairwise data for lrt diff genes
pairwise_log2fc_interaction_df <- readRDS(file = "data_output/2024-1_diffexpr_output/pairwise_log2fc_interaction_df.Rds")
pairwise_log2fc_housing_df <- readRDS(file = "data_output/2024-1_diffexpr_output/pairwise_log2fc_housing_df.Rds")
pairwise_log2fc_genotype_df <- readRDS(file = "data_output/2024-1_diffexpr_output/pairwise_log2fc_genotype_df.Rds")


# Go terms for pairwise comparisons
## Housing
go_housing <- list()

go_housing$all_entrezid <- compareCluster(geneClusters = comp_list$housing_all_entrezid,
                                          organism = "fly",
                                          fun = enrichPathway, qvalueCutoff = 0.05, readable = TRUE)

go_housing$up_entrezid <- compareCluster(geneClusters = comp_list$housing_up_entrezid,
                                         organism = "fly",
                                         fun = enrichPathway, qvalueCutoff = 0.05, readable = TRUE)

go_housing$down_entrezid <- compareCluster(geneClusters = comp_list$housing_down_entrezid,
                                           organism = "fly",
                                           fun = enrichPathway, qvalueCutoff = 0.05, readable = TRUE)

## Genotype
go_genotype <- list()

go_genotype$gh_all_entrezid <- compareCluster(geneClusters = comp_list$genotype_gh_all_entrezid,
                                              organism = "fly",
                                              fun = enrichPathway, qvalueCutoff = 0.05, readable = TRUE)

go_genotype$gh_up_entrezid <- compareCluster(geneClusters = comp_list$genotype_gh_up_entrezid,
                                             organism = "fly",
                                             fun = enrichPathway, qvalueCutoff = 0.05, readable = TRUE)

go_genotype$gh_down_entrezid <- compareCluster(geneClusters = comp_list$genotype_gh_down_entrezid,
                                               organism = "fly",
                                               fun = enrichPathway, qvalueCutoff = 0.05, readable = TRUE)

go_genotype$sh_all_entrezid <- compareCluster(geneClusters = comp_list$genotype_sh_all_entrezid,
                                              organism = "fly",
                                              fun = enrichPathway, qvalueCutoff = 0.05, readable = TRUE)

go_genotype$sh_up_entrezid <- compareCluster(geneClusters = comp_list$genotype_sh_up_entrezid,
                                             organism = "fly",
                                             fun = enrichPathway, qvalueCutoff = 0.05, readable = TRUE)

go_genotype$sh_down_entrezid <- compareCluster(geneClusters = comp_list$genotype_sh_down_entrezid,
                                               organism = "fly",
                                               fun = enrichPathway, qvalueCutoff = 0.05, readable = TRUE)

## Interaction

go_lrt <- list()

comp_list_housing_entrezid <- pairwise_log2fc_housing_df %>%
  dplyr::select(entrezid, sign_paste) %>%
  split(.$sign_paste) %>%
  map("entrezid")

go_lrt$housing <- compareCluster(geneClusters = comp_list_housing_entrezid,
                                       organism = "fly",
                                       fun = enrichPathway, qvalueCutoff = 0.05, readable = TRUE)


comp_list_genotype_entrezid <- pairwise_log2fc_genotype_df %>%
  dplyr::select(entrezid, sign_paste_genotype_gh) %>%
  split(.$sign_paste_genotype_gh) %>%
  map("entrezid")

go_lrt$genotype <- compareCluster(geneClusters = comp_list_genotype_entrezid,
                                          organism = "fly",
                                          fun = enrichPathway, qvalueCutoff = 0.05, readable = TRUE)


comp_list_interaction_entrezid <- pairwise_log2fc_interaction_df %>%
  dplyr::select(entrezid, sign_paste) %>%
  split(.$sign_paste) %>%
  map("entrezid")

go_lrt$interaction <- compareCluster(geneClusters = comp_list_interaction_entrezid,
                                          organism = "fly",
                                          fun = enrichPathway, qvalueCutoff = 0.05, readable = TRUE)
## Write go objects
saveRDS(go_housing,file = "data_output/2024-1_diffexpr_output/reactome_housing.Rds")
saveRDS(go_genotype,file = "data_output/2024-1_diffexpr_output/reactome_genotype.Rds")
saveRDS(go_lrt,file = "data_output/2024-1_diffexpr_output/reactome_lrt.Rds")


