counts_path=./data_output/2024-1_diffexpr_input/
rds=${counts_path}social_mutants.salmon.merged.gene_counts.rds
basename=social_mutants

# Build deseq obj
Rscript ./scripts/diffexpr/diffexpr_helper1_build_deseqobj.R

# LRT tests (anova-like)  
# Rscript script.R output_name output_dir/ deseq_obj_file.rds full_model reduced_model
Rscript ./scripts/diffexpr/diffexpr_step1a_test_lrt.R ${basename}_interaction data_output/2024-1_diffexpr_output/ $rds "~ housing + genotype + housing:genotype" "~housing + genotype";
Rscript ./scripts/diffexpr/diffexpr_step1a_test_lrt.R ${basename}_housing data_output/2024-1_diffexpr_output/ $rds "~  genotype + housing" "~genotype";
Rscript ./scripts/diffexpr/diffexpr_step1a_test_lrt.R ${basename}_genotype data_output/2024-1_diffexpr_output/ $rds "~ housing + genotype" "~housing";

# Wald tests (pairwise)  
# Rscript script.R output_name output_dir/ deseq_obj_file.rds
Rscript ./scripts/diffexpr/diffexpr_step1b_test_wald_housing.R $basename data_output/2024-1_diffexpr_output/ $rds;
Rscript ./scripts/diffexpr/diffexpr_step1c_test_wald_genotype.R $basename data_output/2024-1_diffexpr_output/ $rds;

# Join together results in single by gene results data frame
Rscript ./scripts/diffexpr/diffexpr_helper2_reshape_res_from_list_to_df.R

# Get vst for PCA
Rscript ./scripts/diffexpr/diffexpr_helper3_counts_to_vst.R

# Data wrangling to get matrices and gene data for plots
Rscript ./scripts/diffexpr/diffexpr_helper4_data_wrangling.R

# Go term enrichment 
Rscript ./scripts/diffexpr/diffexpr_step2_go_terms.R