library(DESeq2)
library(tidyverse)
# Functions
pre_filtering <- function(deseq_obj, threhshold = 10){
  keep <- rowSums(counts(deseq_obj)) >= threhshold
  deseq_obj <- deseq_obj[keep,]
}


# Inputs
inputs <- commandArgs(trailingOnly = TRUE)


res_names <- inputs[1]
out_dir <- inputs[2]
deseq_obj_file <- inputs[3]



# Load deseq obj
deseq_obj <- readRDS(deseq_obj_file)
colnames(deseq_obj) <- str_remove(basename(colData(deseq_obj)$bam.files),".target.dedup.sorted.bam")
colData(deseq_obj)$housing <- as.factor(str_split(colnames(deseq_obj),pattern = "_",n = 6) %>% map_chr(3))
colData(deseq_obj)$genotype <- as.factor(str_split(colnames(deseq_obj),pattern = "_",n = 6) %>% map_chr(2))
colData(deseq_obj)$condition <- as.factor(paste(colData(deseq_obj)$genotype,colData(deseq_obj)$housing, sep = "_"))
deseq_obj@design <- ~ condition

genotype_list <- levels(colData(deseq_obj)$genotype)[-1]

contrast_list <- lapply(split(data.frame("contrast"= "condition",
"level1" = levels(colData(deseq_obj)$condition)[-(1:2)],
"level2" = rep(c("CS_GH","CS_SH"),4)),f = 1:(2*(length(genotype_list)))), unlist)

genotype_deseq_list <- map(contrast_list,function(.x, deseq_obj){
  
  res <- deseq_obj[,colData(deseq_obj)$condition %in% .x]
  res$condition <- droplevels(res$condition)
  return(res)
}, deseq_obj = deseq_obj)

# Prefiltering
genotype_deseq_list <- map(genotype_deseq_list,pre_filtering)


# Build model
genotype_deseq_list <- map(genotype_deseq_list,DESeq)


# Test
res_list <- map2(contrast_list, genotype_deseq_list,function(.x,.y) results(contrast = .x, 
                                                    object = .y,
                                                    alpha = 0.05) )


# Shrinkage
res_ashr_list <- map2(genotype_deseq_list, res_list,function(x, y) lfcShrink(dds = x,
                                                     res = y,
                                                     type = "ashr"))


# Get diff ranges

get_ranges_up <- function(res, deseq, padj = 0.05,lfc = 0){
  granges(rowRanges(deseq)[tidyr::replace_na(res$padj,1) < padj & tidyr::replace_na(res$log2FoldChange,0) > lfc])
}
get_ranges_down <- function(res, deseq, padj = 0.05,lfc = 0){
  granges(rowRanges(deseq)[tidyr::replace_na(res$padj,1) < padj & tidyr::replace_na(res$log2FoldChange,0) < lfc])
}

res_up_list <- map2(res_ashr_list,genotype_deseq_list, function(x,y)get_ranges_up(res = x, deseq = y))
res_down_list <- map2(res_ashr_list,genotype_deseq_list, function(x,y)get_ranges_down(res = x, deseq = y))

names(res_ashr_list) <- paste0(res_names,"_",unlist(contrast_list %>% map(2)))

# Write results objects

names(res_ashr_list) <- paste0(out_dir,names(res_ashr_list) ,"_genotype_ashr.Rds")

iwalk(res_ashr_list, ~ saveRDS(.x, .y))

# Write results tables
deseq_obj_contrasts_ranges <- map(genotype_deseq_list,~as.data.frame(granges(rowRanges(.x)))[,1:3])

res_ashr_gr_df_list <- map2(res_ashr_list,deseq_obj_contrasts_ranges, ~ cbind(.y,as.data.frame(.x)))
names(res_ashr_gr_df_list) <- str_replace(names(res_ashr_list),".Rds",".tsv")

iwalk(res_ashr_gr_df_list, ~ write_tsv(.x, .y))

# Write deseq fitted objects
names(deseq_obj_contrasts_ranges) <- str_replace(names(res_ashr_list),".Rds","_deseq_obj_fitted.Rds")
iwalk(deseq_obj_contrasts_ranges, saveRDS)

# DiffBind bedfiles
names(res_up_list) <- str_replace(names(res_ashr_list),".Rds","_up.bed")
iwalk(res_up_list, ~ rtracklayer::export.bed(.x,.y))

names(res_down_list) <- str_replace(names(res_ashr_list),".Rds","_down.bed")
iwalk(res_down_list, ~ rtracklayer::export.bed(.x,.y))

