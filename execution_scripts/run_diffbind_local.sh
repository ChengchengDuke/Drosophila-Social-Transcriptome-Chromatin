counts_path=./data_output/2024-01_diff_cutandrun_data_output/

# Filter out low count samples
Rscript ./scripts/diffbind/diffbind_helper1_filter_out_low_count_samples.R

# Pairwise Wald test housing
for rds in ${counts_path}*seacrigg*Rds;
  do 
    echo $rds;
    basename=${rds##*/}
    basename=${basename%.*}
    Rscript ./scripts/diffbind/diffbind_step2a_test_wald_housing.R $basename data_output/2024-01_diffbind_cutnrun/ $rds;
  done

# Pairwise Wald test genotype
for rds in ${counts_path}*seacrigg*Rds;
  do 
    echo $rds;
    basename=${rds##*/}
    basename=${basename%.*}
    Rscript ./scripts/diffbind/diffbind_step2b_test_wald_genotype.R $basename data_output/2024-01_diffbind_cutnrun/ $rds;
  done

# LRT interaction
for rds in ${counts_path}*seacrigg*Rds;
  do 
    echo $rds;
    basename=${rds##*/}
    basename=${basename%.*}
    Rscript ./scripts/diffbind/diffbind_step2c_test_lrt.R ${basename}_interaction data_output/2024-01_diffbind_cutnrun/ $rds "~ housing + genotype + housing:genotype" "~housing + genotype";
  done
  
# LRT housing
for rds in ${counts_path}*seacrigg*Rds;
  do 
    echo $rds;
    basename=${rds##*/}
    basename=${basename%.*}
    Rscript ./scripts/diffbind/diffbind_step2c_test_lrt.R ${basename}_housing data_output/2024-01_diffbind_cutnrun/ $rds "~  genotype + housing" "~genotype";
  done

# LRT genotype
for rds in ${counts_path}*seacrigg*Rds;
  do 
    echo $rds;
    basename=${rds##*/}
    basename=${basename%.*}
    Rscript ./scripts/diffbind/diffbind_step2c_test_lrt.R ${basename}_genotype data_output/2024-01_diffbind_cutnrun/ $rds "~ housing + genotype" "~housing";
  done

# Data re-formatting to make a single data frame with all results
Rscript ./scripts/diffbind/diffbind_helper2_reshape_res_from_list_to_df.R

# Get vst for PCA
Rscript ./scripts/diffbind/diffbind_helper3_counts_to_vst.R

# Data wrangling to get matrices and gene data for plots
Rscript ./scripts/diffbind/diffbind_helper4_data_wrangling.R

# Go term enrichment 
Rscript ./scripts/diffbind/diffbind_step3_go_terms.R

