# CUT&RUN data analysis

Contact:
jesus.sotelo.fonseca@duke.edu
jemilianosf@gmail.com

Volkan lab 

# Social experience project
## CUT&RUN data analysis


changes from last version
-RNAseq
--diffexpr_step1b_test_wald_housing.R SH/GH changed to GH/SH
--step2 go changed simpilified
--step2b reactome one mistake was corrected: geneotype sh down was wrongly typed as GH, changed to SH 
--step2c KEGG added change format, and one mistake was corrected: geneotype sh down was wrongly typed as GH, changed to SH 
--DEGs numbers changed a little becuase of update of packages and R
--order of genotypes to display changed  

-CUT&RUN
--diffbind step2a test_wald_housing SH/GH changed to GH/SH
--diffbind step2b test_wald_genotype saveRDS(object = deseq_obj, paste0(out_dir,res_names,"_deseq_obj_housing_fitted.Rds")) changed to saveRDS(object = deseq_obj, paste0(out_dir,res_names,"_deseq_obj_genotype_fitted.Rds"))
--diffbind helper4 saved res_gr_h3k4me3 and res_gr_rnapolii, changed peaks to genes: link all peaks to genes as possible
--added reactome and KEGG


-RNA VS CUTNRUN
--updated comparision: expanded to all peaks across the whole gene
--GO changed to RNA overlap with at least one chromatin mark concordant groups
--updated heatmconcordantion heatmaps with more peaks and color coded



