library(DESeq2)
library(tidyverse)
library(sva)
# Read matrix
salmon_merged_gene_counts <- read_tsv("data_input/2024-01_diff_rnaseq_data/salmon.merged.gene_counts.tsv")
# Remove id cols, convert to matrix, add back gene names as rownmames
deseq_mat <- salmon_merged_gene_counts %>%
                         select(-gene_id, -gene_name) %>%
  mutate(across(everything(), as.integer)) %>%
  as.matrix()
rownames(deseq_mat) <- salmon_merged_gene_counts$gene_name

# Build deseq obj
deseq_obj <- DESeqDataSetFromMatrix(countData = deseq_mat,colData = DataFrame("sample_name"=colnames(deseq_mat)),design = ~1)
# Remove not relevant columns
deseq_obj <- deseq_obj[,!str_detect(colnames(deseq_obj), "B04")]
deseq_obj <- deseq_obj[,!str_detect(colnames(deseq_obj), "sub")]
deseq_obj <- deseq_obj[,!str_detect(colnames(deseq_obj), "FEM")]
deseq_obj <- deseq_obj[,!str_detect(colnames(deseq_obj), "Z144")]
deseq_obj <- deseq_obj[,!str_detect(colnames(deseq_obj), ".HIGH")]
deseq_obj <- deseq_obj[,!str_detect(colnames(deseq_obj), ".LOW")]
deseq_obj <- deseq_obj[,!str_detect(colnames(deseq_obj), "_M_")]

# Add batch info
rnaseq_batches <- read_csv("data_input/2024-01_diff_rnaseq_data/rnaseq_batches.csv")
rnaseq_batches$rownames <- str_remove(rnaseq_batches$rownames,"_V_MAL")
rnaseq_batches$rownames <- str_replace(rnaseq_batches$rownames,"7B","7b")
rnaseq_batches$rownames <- str_replace(rnaseq_batches$rownames,"7D","7d")

rnaseq_batches <- rnaseq_batches %>%
  select(rownames,dissection)

deseq_obj_coldata <- left_join(as.data.frame(colData(deseq_obj)),rnaseq_batches, by = c("sample_name" = "rownames"))
deseq_obj_coldata <- deseq_obj_coldata %>%
  mutate(dissection = ifelse(is.na(dissection),"cd",dissection)) %>%
  mutate(dissection = ifelse(dissection == "ac","ac","notac"))
deseq_obj_coldata <- DataFrame(deseq_obj_coldata)
rownames(deseq_obj_coldata) <- deseq_obj_coldata$sample_name
colData(deseq_obj) <- deseq_obj_coldata

######
colData(deseq_obj)$housing <- as.factor(str_split(colnames(deseq_obj),pattern = "_",n = 3) %>% map_chr(2))
colData(deseq_obj)$genotype <- as.factor(str_split(colnames(deseq_obj),pattern = "_",n = 3) %>% map_chr(1))
colData(deseq_obj)$condition <- as.factor(paste0(colData(deseq_obj)$genotype,"_",colData(deseq_obj)$housing))

### Filter out non mRNA genes
filt_out_rna_rownames <- rownames(deseq_obj) %>%
  str_split(":") %>% map_chr(1) 

filt_out_rna <- rownames(deseq_obj)[str_detect(rownames(deseq_obj),":")] %>%
  str_split(":") %>% map_chr(1) 
filt_out_rna <- filt_out_rna %>%
  unique()
filt_out_rna <- filt_out_rna[!filt_out_rna %in% c("lncRNA","asRNA")]
deseq_obj <- deseq_obj[!filt_out_rna_rownames %in% filt_out_rna,]

### Filter out mobile elements
deseq_obj <- deseq_obj[!str_detect(rownames(deseq_obj),"\\{")]

# Prefilter low counts

pre_filtering <- function(deseq_obj, threhshold = 10){
  
  keep <- rowSums(counts(deseq_obj) >= threhshold) >= 3
  
  deseq_obj <- deseq_obj[keep,]
  return(deseq_obj)
}

# Prefiltering
deseq_obj <- pre_filtering(deseq_obj) # remove really low counts

###
saveRDS(deseq_obj,"data_output/2024-1_diffexpr_input/social_mutants_raw_counts.salmon.merged.gene_counts.rds")

# Batch correction 

###### Batch correction
# Batch correction with ComBat 
combat_cts <- ComBat_seq(counts(deseq_obj),batch =deseq_obj$dissection) # apply combat function to raw counts, removing the variation due to dissection 
combat_int <- apply(combat_cts,2,as.integer) # Convert each column to integer
rownames(combat_int) <- rownames(combat_cts)

# Export batch corrected counts to a  text table
write.table(combat_cts, "data_output/2024-1_diffexpr_input/count_matrix_combat.tsv", quote = FALSE, sep = "\t",row.names = T, col.names = T)

# Change the raw counts to combat adujsted in our DESeq2 object
counts(deseq_obj) <- combat_int

saveRDS(deseq_obj,"data_output/2024-1_diffexpr_input/social_mutants.salmon.merged.gene_counts.rds")






