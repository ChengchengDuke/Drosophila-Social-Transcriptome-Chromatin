counts_path=./data_input/2024-01_diff_cutandrun_data_output/

for rds in ${counts_path}*seacrigg*Rds;
  do 
    echo $rds;
    basename=${rds##*/}
    basename=${basename%.*}
    Rscript ./scripts/diffbind_step2c_test_wald.R $basename data_output/2024-01_diffbind_cutnrun/ $rds;
  done

for rds in ${counts_path}*seacrigg*Rds;
  do 
    echo $rds;
    basename=${rds##*/}
    basename=${basename%.*}
    Rscript ./scripts/diffbind_step2d_test_wald_genotype.R $basename data_output/2024-01_diffbind_cutnrun/ $rds;
  done



for rds in ${counts_path}*seacrigg*Rds;
  do 
    echo $rds;
    basename=${rds##*/}
    basename=${basename%.*}
    Rscript ./scripts/diffbind_step2b_test_lrt.R ${basename}_interaction data_output/2024-01_diffbind_cutnrun/ $rds "~ housing + genotype + housing:genotype" "~housing + genotype";
  done
  
for rds in ${counts_path}*seacrigg*Rds;
  do 
    echo $rds;
    basename=${rds##*/}
    basename=${basename%.*}
    Rscript ./scripts/diffbind_step2b_test_lrt.R ${basename}_housing data_output/2024-01_diffbind_cutnrun/ $rds "~  genotype + housing" "~genotype";
  done
  
for rds in ${counts_path}*seacrigg*Rds;
  do 
    echo $rds;
    basename=${rds##*/}
    basename=${basename%.*}
    Rscript ./scripts/diffbind_step2b_test_lrt.R ${basename}_genotype data_output/2024-01_diffbind_cutnrun/ $rds "~ housing + genotype" "~housing";
  done


