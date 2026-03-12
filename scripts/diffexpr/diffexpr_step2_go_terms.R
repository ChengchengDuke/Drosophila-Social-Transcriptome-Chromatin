library(tidyverse)
library(clusterProfiler)
library(org.Dm.eg.db)

# Read data Ashr results
comp_list <- readRDS(file = "data_output/2024-1_diffexpr_output/comp_list.Rds")

# Read data Ashr genotype results - by sign
pairwise_log2fc_genotype_df_ashr_GH <- readRDS("./data_output/2024-1_diffexpr_output/pairwise_log2fc_genotype_df_ashr_GH.Rds")

pairwise_log2fc_genotype_df_ashr_SH <- readRDS("./data_output/2024-1_diffexpr_output/pairwise_log2fc_genotype_df_ashr_SH.Rds")


# Pairwise data for lrt diff genes
pairwise_log2fc_interaction_df <- readRDS(file = "data_output/2024-1_diffexpr_output/pairwise_log2fc_interaction_df.Rds")
pairwise_log2fc_housing_df <- readRDS(file = "data_output/2024-1_diffexpr_output/pairwise_log2fc_housing_df.Rds")
pairwise_log2fc_genotype_df <- readRDS(file = "data_output/2024-1_diffexpr_output/pairwise_log2fc_genotype_df.Rds")


# Go terms for pairwise comparisons
## Housing
go_housing <- list()

go_housing$all_entrezid <- compareCluster(geneClusters = comp_list$housing_all_entrezid,
                                          OrgDb = org.Dm.eg.db,
                                          fun = enrichGO, qvalueCutoff = 0.05, readable = TRUE, ont = "ALL")

go_housing$up_entrezid <- compareCluster(geneClusters = comp_list$housing_up_entrezid,
                                         OrgDb = org.Dm.eg.db,
                                         fun = enrichGO, qvalueCutoff = 0.05, readable = TRUE, ont = "ALL")

go_housing$down_entrezid <- compareCluster(geneClusters = comp_list$housing_down_entrezid,
                                           OrgDb = org.Dm.eg.db,
                                           fun = enrichGO, qvalueCutoff = 0.05, readable = TRUE, ont = "ALL")

go_housing_all_simplify <- clusterProfiler::simplify(x = go_housing$all_entrezid, cutoff = 0.6)



####################
CS_list <- list()
CS_list$CS_all <- comp_list$housing_all_entrezid$CS
CS_list$CS_SH_up <- comp_list$housing_up_entrezid$CS
CS_list$CS_SH_down <- comp_list$housing_down_entrezid$CS

go_housing_CS <- list()

go_housing_CS <- compareCluster(geneClusters = CS_list,
                                OrgDb = org.Dm.eg.db,
                                fun = enrichGO, qvalueCutoff = 0.05, readable = TRUE, ont = "ALL")

go_housing_CS_simplify <- clusterProfiler::simplify(x = go_housing_CS, cutoff = 0.6)


## Genotype
go_genotype <- list()

go_genotype$gh_all_entrezid <- compareCluster(geneClusters = comp_list$genotype_gh_all_entrezid,
                                              OrgDb = org.Dm.eg.db,
                                              fun = enrichGO, qvalueCutoff = 0.05, readable = TRUE, ont = "ALL")

go_genotype$gh_up_entrezid <- compareCluster(geneClusters = comp_list$genotype_gh_up_entrezid,
                                             OrgDb = org.Dm.eg.db,
                                             fun = enrichGO, qvalueCutoff = 0.05, readable = TRUE, ont = "ALL")

go_genotype$gh_down_entrezid <- compareCluster(geneClusters = comp_list$genotype_gh_down_entrezid,
                                               OrgDb = org.Dm.eg.db,
                                               fun = enrichGO, qvalueCutoff = 0.05, readable = TRUE, ont = "ALL")

go_genotype_gh_all_simplify <- clusterProfiler::simplify(x = go_genotype$gh_all_entrezid, cutoff = 0.6)


go_genotype$sh_all_entrezid <- compareCluster(geneClusters = comp_list$genotype_sh_all_entrezid,
                                              OrgDb = org.Dm.eg.db,
                                              fun = enrichGO, qvalueCutoff = 0.05, readable = TRUE, ont = "ALL")

go_genotype$sh_up_entrezid <- compareCluster(geneClusters = comp_list$genotype_sh_up_entrezid,
                                             OrgDb = org.Dm.eg.db,
                                             fun = enrichGO, qvalueCutoff = 0.05, readable = TRUE, ont = "ALL")

go_genotype$sh_down_entrezid <- compareCluster(geneClusters = comp_list$genotype_sh_down_entrezid,
                                               OrgDb = org.Dm.eg.db,
                                               fun = enrichGO, qvalueCutoff = 0.05, readable = TRUE, ont = "ALL")

go_genotype_sh_all_simplify <- clusterProfiler::simplify(x = go_genotype$sh_all_entrezid, cutoff = 0.6)



## Ashr by sign 

go_Ashr_sign <- list()

Ashr_sign_genotype_GH_entrezid <- pairwise_log2fc_genotype_df_ashr_GH %>%
  dplyr::select(entrezid, sign_paste_genotype_gh) %>%
  split(.$sign_paste_genotype_gh) %>%
  map("entrezid") %>%
  map(unique)

go_Ashr_sign$genotype_GH <- compareCluster(geneClusters = Ashr_sign_genotype_GH_entrezid,
                                  OrgDb = org.Dm.eg.db,
                                  fun = enrichGO, qvalueCutoff = 0.05, readable = TRUE, ont = "ALL")


Ashr_sign_genotype_SH_entrezid <- pairwise_log2fc_genotype_df_ashr_SH %>%
  dplyr::select(entrezid, sign_paste_genotype_sh) %>%
  split(.$sign_paste_genotype_sh) %>%
  map("entrezid") %>%
  map(unique)

go_Ashr_sign$genotype_SH <- compareCluster(geneClusters = Ashr_sign_genotype_SH_entrezid,
                                           OrgDb = org.Dm.eg.db,
                                           fun = enrichGO, qvalueCutoff = 0.05, readable = TRUE, ont = "ALL")

go_Ashr_sign_GH_simplify <- clusterProfiler::simplify(go_Ashr_sign$genotype_GH, cutoff = 0.6)
go_Ashr_sign_SH_simplify <- clusterProfiler::simplify(go_Ashr_sign$genotype_SH, cutoff = 0.6)





## LRT Interaction

go_lrt <- list()

comp_list_housing_entrezid <- pairwise_log2fc_housing_df %>%
  dplyr::select(entrezid, sign_paste) %>%
  split(.$sign_paste) %>% # splits the data frame by sign paste
  map("entrezid")

go_lrt$housing <- compareCluster(geneClusters = comp_list_housing_entrezid,
                                       OrgDb = org.Dm.eg.db,
                                       fun = enrichGO, qvalueCutoff = 0.05, readable = TRUE, ont = "ALL")


comp_list_genotype_entrezid <- pairwise_log2fc_genotype_df %>%
  dplyr::select(entrezid, sign_paste_genotype_gh) %>%
  split(.$sign_paste_genotype_gh) %>%
  map("entrezid")

go_lrt$genotype <- compareCluster(geneClusters = comp_list_genotype_entrezid,
                                          OrgDb = org.Dm.eg.db,
                                          fun = enrichGO, qvalueCutoff = 0.05, readable = TRUE, ont = "ALL")


comp_list_interaction_entrezid <- pairwise_log2fc_interaction_df %>%
  dplyr::select(entrezid, sign_paste) %>%
  split(.$sign_paste) %>%
  map("entrezid")

go_lrt$interaction <- compareCluster(geneClusters = comp_list_interaction_entrezid,
                                          OrgDb = org.Dm.eg.db,
                                          fun = enrichGO, qvalueCutoff = 0.05, readable = TRUE, ont = "ALL")

go_lrt_interaction_simplify <- clusterProfiler::simplify(go_lrt$interaction, cutoff = 0.6)


## Write go objects
saveRDS(go_housing,file = "data_output/2024-1_diffexpr_output/go_housing.Rds")
saveRDS(go_housing_all_simplify,file = "data_output/2024-1_diffexpr_output/go_housing_all_simplify.Rds")
saveRDS(go_housing_CS_simplify,file = "data_output/2024-1_diffexpr_output/go_housing_CS_simplify.Rds")

saveRDS(go_genotype,file = "data_output/2024-1_diffexpr_output/go_genotype.Rds")
saveRDS(go_genotype_gh_all_simplify,file = "data_output/2024-1_diffexpr_output/go_genotype_gh_all_simplify.Rds")
saveRDS(go_genotype_sh_all_simplify,file = "data_output/2024-1_diffexpr_output/go_genotype_sh_all_simplify.Rds")

saveRDS(go_Ashr_sign,file = "data_output/2024-1_diffexpr_output/go_Ashr_sign.Rds")
saveRDS(go_Ashr_sign_GH_simplify,file = "data_output/2024-1_diffexpr_output/go_Ashr_sign_GH_simplify.Rds")
saveRDS(go_Ashr_sign_SH_simplify,file = "data_output/2024-1_diffexpr_output/go_Ashr_sign_SH_simplify.Rds")

saveRDS(go_lrt,file = "data_output/2024-1_diffexpr_output/go_lrt.Rds")
saveRDS(go_lrt_interaction_simplify,file = "data_output/2024-1_diffexpr_output/go_lrt_interaction_simplify.Rds")


