counts_path=/Volumes/ExtremePro/Emi/social_experience/cutandrun_dec22/differential_analysis/data_output/

for rds in ${counts_path}*seacrigg*Rds;
  do 
    echo $rds;
    basename=${rds##*/}
    basename=${basename%.*}
    Rscript ./scripts/diffbind_step2_test.R housing SH GH $basename data_output/diffbind_cutnrun/ $rds;
  done

for rds in ${counts_path}*promoter*Rds;
  do 
    echo $rds;
    basename=${rds##*/}
    basename=${basename%.*}
    Rscript ./scripts/diffbind_step2_test.R housing SH GH $basename data_output/diffbind_cutnrun/ $rds;
  done

for rds in ${counts_path}*IgG_seacrigg*Rds;
  do 
    echo $rds;
    basename=${rds##*/}
    basename=${basename%.*}
    Rscript ./scripts/diffbind_step2_test.R housing SH GH $basename data_output/diffbind_cutnrun/ $rds;
  done
