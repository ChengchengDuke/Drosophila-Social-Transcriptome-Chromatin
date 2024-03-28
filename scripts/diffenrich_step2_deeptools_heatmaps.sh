conda activate deeptools

## ## ## 
## Heatmap of all samples, all genes
## ## ## 
antibody="H3K4me3"
bigwig_path="data_output/03_bed_to_bigwig/"
# Change bw_signals, put: signal1.bigWig signal2.bigWig signal3.bigWig ...
bw_signals="${bigwig_path}${antibody}_CS_GH_male_brain_R3.bigWig ${bigwig_path}${antibody}_CS_GH_male_brain_R4.bigWig ${bigwig_path}${antibody}_CS_GH_male_brain_R5.bigWig ${bigwig_path}${antibody}_CS_SH_male_brain_R3.bigWig ${bigwig_path}${antibody}_CS_SH_male_brain_R4.bigWig ${bigwig_path}${antibody}_CS_SH_male_brain_R5.bigWig ${bigwig_path}${antibody}_Or47bmut_GH_male_brain_R1.bigWig ${bigwig_path}${antibody}_Or47bmut_GH_male_brain_R2.bigWig ${bigwig_path}${antibody}_Or47bmut_GH_male_brain_R3.bigWig ${bigwig_path}${antibody}_Or47bmut_SH_male_brain_R1.bigWig ${bigwig_path}${antibody}_Or47bmut_SH_male_brain_R2.bigWig ${bigwig_path}${antibody}_Or47bmut_SH_male_brain_R3.bigWig ${bigwig_path}${antibody}_Or67dmut_GH_male_brain_R1.bigWig ${bigwig_path}${antibody}_Or67dmut_GH_male_brain_R2.bigWig ${bigwig_path}${antibody}_Or67dmut_GH_male_brain_R3.bigWig ${bigwig_path}${antibody}_Or67dmut_SH_male_brain_R1.bigWig ${bigwig_path}${antibody}_Or67dmut_SH_male_brain_R2.bigWig ${bigwig_path}${antibody}_Or67dmut_SH_male_brain_R3.bigWig ${bigwig_path}${antibody}_dsxmut_GH_male_brain_R1.bigWig ${bigwig_path}${antibody}_dsxmut_GH_male_brain_R2.bigWig ${bigwig_path}${antibody}_dsxmut_GH_male_brain_R3.bigWig ${bigwig_path}${antibody}_dsxmut_GH_male_brain_R4.bigWig ${bigwig_path}${antibody}_dsxmut_SH_male_brain_R1.bigWig ${bigwig_path}${antibody}_dsxmut_SH_male_brain_R2.bigWig ${bigwig_path}${antibody}_dsxmut_SH_male_brain_R3.bigWig ${bigwig_path}${antibody}_frumut_GH_male_brain_R1.bigWig ${bigwig_path}${antibody}_frumut_GH_male_brain_R2.bigWig ${bigwig_path}${antibody}_frumut_GH_male_brain_R3.bigWig ${bigwig_path}${antibody}_frumut_SH_male_brain_R1.bigWig ${bigwig_path}${antibody}_frumut_SH_male_brain_R2.bigWig ${bigwig_path}${antibody}_frumut_SH_male_brain_R3.bigWig"
# Change bed_regions, put: my_genes1.bed my_genes2.bed ...
bed_regions="data_output/2024-1_diffexpr_output/all_genes.bed"
# Change the names of sets of regions (row groups): my_genes1 my_genes2 ... 
regions_label="Random_genes"
# Change the names of the columns: signal1 signal2 signal3
samples_label="${antibody}_CS_GH_R3 ${antibody}_CS_GH_R4 ${antibody}_CS_GH_R5 ${antibody}_CS_SH_R3 ${antibody}_CS_SH_R4 ${antibody}_CS_SH_R5 ${antibody}_Or47bmut_GH_R1 ${antibody}_Or47bmut_GH_R2 ${antibody}_Or47bmut_GH_R3 ${antibody}_Or47bmut_SH_R1 ${antibody}_Or47bmut_SH_R2 ${antibody}_Or47bmut_SH_R3 ${antibody}_Or67dmut_GH_R1 ${antibody}_Or67dmut_GH_R2 ${antibody}_Or67dmut_GH_R3 ${antibody}_Or67dmut_SH_R1 ${antibody}_Or67dmut_SH_R2 ${antibody}_Or67dmut_SH_R3 ${antibody}_dsxmut_GH_R1 ${antibody}_dsxmut_GH_R2 ${antibody}_dsxmut_GH_R3 ${antibody}_dsxmut_GH_R4 ${antibody}_dsxmut_SH_R1 ${antibody}_dsxmut_SH_R2 ${antibody}_dsxmut_SH_R3 ${antibody}_frumut_GH_R1 ${antibody}_frumut_GH_R2 ${antibody}_frumut_GH_R3 ${antibody}_frumut_SH_R1 ${antibody}_frumut_SH_R2 ${antibody}_frumut_SH_R3"
output_mat="data_output/deeptools_heatmaps/"${antibody}_all_samples_all_genes.mat.gz
output_png="data_output/deeptools_heatmaps/"${antibody}_all_samples_all_genes.png

# Matrix: 
computeMatrix reference-point --referencePoint center -S $bw_signals -R $bed_regions --skipZeros -b 3000 -a 3000 -o $output_mat

# Heatmap: 
plotHeatmap -m $output_mat \
--averageTypeSummaryPlot mean \
--boxAroundHeatmaps no \
--missingDataColor grey \
--regionsLabel $regions_label \
--samplesLabel $samples_label \
-out $output_png

## ## ## 
## Housing diff genes
## ## ## 

antibody="H3K4me3"
bigwig_path="data_output/03_bed_to_bigwig/"
bw_signals="${bigwig_path}${antibody}_CS_GH_male_brain_R3.bigWig ${bigwig_path}${antibody}_CS_GH_male_brain_R4.bigWig ${bigwig_path}${antibody}_CS_GH_male_brain_R5.bigWig ${bigwig_path}${antibody}_CS_SH_male_brain_R3.bigWig ${bigwig_path}${antibody}_CS_SH_male_brain_R4.bigWig ${bigwig_path}${antibody}_CS_SH_male_brain_R5.bigWig ${bigwig_path}${antibody}_Or47bmut_GH_male_brain_R1.bigWig ${bigwig_path}${antibody}_Or47bmut_GH_male_brain_R2.bigWig ${bigwig_path}${antibody}_Or47bmut_GH_male_brain_R3.bigWig ${bigwig_path}${antibody}_Or47bmut_SH_male_brain_R1.bigWig ${bigwig_path}${antibody}_Or47bmut_SH_male_brain_R2.bigWig ${bigwig_path}${antibody}_Or47bmut_SH_male_brain_R3.bigWig ${bigwig_path}${antibody}_Or67dmut_GH_male_brain_R1.bigWig ${bigwig_path}${antibody}_Or67dmut_GH_male_brain_R2.bigWig ${bigwig_path}${antibody}_Or67dmut_GH_male_brain_R3.bigWig ${bigwig_path}${antibody}_Or67dmut_SH_male_brain_R1.bigWig ${bigwig_path}${antibody}_Or67dmut_SH_male_brain_R2.bigWig ${bigwig_path}${antibody}_Or67dmut_SH_male_brain_R3.bigWig ${bigwig_path}${antibody}_dsxmut_GH_male_brain_R1.bigWig ${bigwig_path}${antibody}_dsxmut_GH_male_brain_R2.bigWig ${bigwig_path}${antibody}_dsxmut_GH_male_brain_R3.bigWig ${bigwig_path}${antibody}_dsxmut_GH_male_brain_R4.bigWig ${bigwig_path}${antibody}_dsxmut_SH_male_brain_R1.bigWig ${bigwig_path}${antibody}_dsxmut_SH_male_brain_R2.bigWig ${bigwig_path}${antibody}_dsxmut_SH_male_brain_R3.bigWig ${bigwig_path}${antibody}_frumut_GH_male_brain_R1.bigWig ${bigwig_path}${antibody}_frumut_GH_male_brain_R2.bigWig ${bigwig_path}${antibody}_frumut_GH_male_brain_R3.bigWig ${bigwig_path}${antibody}_frumut_SH_male_brain_R1.bigWig ${bigwig_path}${antibody}_frumut_SH_male_brain_R2.bigWig ${bigwig_path}${antibody}_frumut_SH_male_brain_R3.bigWig"
bed_regions="data_output/2024-1_diffexpr_output/CS_housing_up.bed data_output/2024-1_diffexpr_output/CS_housing_down.bed"
regions_label="CS_GH CS_SH"
samples_label="${antibody}_CS_GH_R3 ${antibody}_CS_GH_R4 ${antibody}_CS_GH_R5 ${antibody}_CS_SH_R3 ${antibody}_CS_SH_R4 ${antibody}_CS_SH_R5 ${antibody}_Or47bmut_GH_R1 ${antibody}_Or47bmut_GH_R2 ${antibody}_Or47bmut_GH_R3 ${antibody}_Or47bmut_SH_R1 ${antibody}_Or47bmut_SH_R2 ${antibody}_Or47bmut_SH_R3 ${antibody}_Or67dmut_GH_R1 ${antibody}_Or67dmut_GH_R2 ${antibody}_Or67dmut_GH_R3 ${antibody}_Or67dmut_SH_R1 ${antibody}_Or67dmut_SH_R2 ${antibody}_Or67dmut_SH_R3 ${antibody}_dsxmut_GH_R1 ${antibody}_dsxmut_GH_R2 ${antibody}_dsxmut_GH_R3 ${antibody}_dsxmut_GH_R4 ${antibody}_dsxmut_SH_R1 ${antibody}_dsxmut_SH_R2 ${antibody}_dsxmut_SH_R3 ${antibody}_frumut_GH_R1 ${antibody}_frumut_GH_R2 ${antibody}_frumut_GH_R3 ${antibody}_frumut_SH_R1 ${antibody}_frumut_SH_R2 ${antibody}_frumut_SH_R3"
output_mat="data_output/deeptools_heatmaps/"${antibody}_all_samples_diff_housing_cs_genes.mat.gz
output_png="data_output/deeptools_heatmaps/"${antibody}_all_samples_diff_housing_cs_genes.png

# Matrix: 
computeMatrix reference-point --referencePoint center -S $bw_signals -R $bed_regions --skipZeros -b 3000 -a 3000 -o $output_mat

# Heatmap: 
plotHeatmap -m $output_mat \
--averageTypeSummaryPlot mean \
--boxAroundHeatmaps no \
--missingDataColor grey \
--regionsLabel $regions_label \
--samplesLabel $samples_label \
-out $output_png

## ## ## 
## Genotype GH diff genes
## ## ## 

antibody="H3K4me3"
bigwig_path="data_output/03_bed_to_bigwig/"
bw_signals="${bigwig_path}${antibody}_CS_GH_male_brain_R3.bigWig ${bigwig_path}${antibody}_CS_GH_male_brain_R4.bigWig ${bigwig_path}${antibody}_CS_GH_male_brain_R5.bigWig ${bigwig_path}${antibody}_CS_SH_male_brain_R3.bigWig ${bigwig_path}${antibody}_CS_SH_male_brain_R4.bigWig ${bigwig_path}${antibody}_CS_SH_male_brain_R5.bigWig ${bigwig_path}${antibody}_Or47bmut_GH_male_brain_R1.bigWig ${bigwig_path}${antibody}_Or47bmut_GH_male_brain_R2.bigWig ${bigwig_path}${antibody}_Or47bmut_GH_male_brain_R3.bigWig ${bigwig_path}${antibody}_Or47bmut_SH_male_brain_R1.bigWig ${bigwig_path}${antibody}_Or47bmut_SH_male_brain_R2.bigWig ${bigwig_path}${antibody}_Or47bmut_SH_male_brain_R3.bigWig ${bigwig_path}${antibody}_Or67dmut_GH_male_brain_R1.bigWig ${bigwig_path}${antibody}_Or67dmut_GH_male_brain_R2.bigWig ${bigwig_path}${antibody}_Or67dmut_GH_male_brain_R3.bigWig ${bigwig_path}${antibody}_Or67dmut_SH_male_brain_R1.bigWig ${bigwig_path}${antibody}_Or67dmut_SH_male_brain_R2.bigWig ${bigwig_path}${antibody}_Or67dmut_SH_male_brain_R3.bigWig ${bigwig_path}${antibody}_dsxmut_GH_male_brain_R1.bigWig ${bigwig_path}${antibody}_dsxmut_GH_male_brain_R2.bigWig ${bigwig_path}${antibody}_dsxmut_GH_male_brain_R3.bigWig ${bigwig_path}${antibody}_dsxmut_GH_male_brain_R4.bigWig ${bigwig_path}${antibody}_dsxmut_SH_male_brain_R1.bigWig ${bigwig_path}${antibody}_dsxmut_SH_male_brain_R2.bigWig ${bigwig_path}${antibody}_dsxmut_SH_male_brain_R3.bigWig ${bigwig_path}${antibody}_frumut_GH_male_brain_R1.bigWig ${bigwig_path}${antibody}_frumut_GH_male_brain_R2.bigWig ${bigwig_path}${antibody}_frumut_GH_male_brain_R3.bigWig ${bigwig_path}${antibody}_frumut_SH_male_brain_R1.bigWig ${bigwig_path}${antibody}_frumut_SH_male_brain_R2.bigWig ${bigwig_path}${antibody}_frumut_SH_male_brain_R3.bigWig"
bed_regions="data_output/2024-1_diffexpr_output/DSX_GH_genotype_down.bed data_output/2024-1_diffexpr_output/DSX_GH_genotype_up.bed data_output/2024-1_diffexpr_output/FRU_GH_genotype_down.bed data_output/2024-1_diffexpr_output/FRU_GH_genotype_up.bed data_output/2024-1_diffexpr_output/Or47b_GH_genotype_down.bed data_output/2024-1_diffexpr_output/Or47b_GH_genotype_up.bed data_output/2024-1_diffexpr_output/Or67d_GH_genotype_down.bed data_output/2024-1_diffexpr_output/Or67d_GH_genotype_up.bed"
regions_label="DSX_GH_UP DSX_GH_DOWN FRU_GH_UP FRU_GH_DOWN Or47b_GH_UP Or47b_GH_DOWN Or67d_GH_UP Or67d_GH_DOWN"
samples_label="${antibody}_CS_GH_R3 ${antibody}_CS_GH_R4 ${antibody}_CS_GH_R5 ${antibody}_CS_SH_R3 ${antibody}_CS_SH_R4 ${antibody}_CS_SH_R5 ${antibody}_Or47bmut_GH_R1 ${antibody}_Or47bmut_GH_R2 ${antibody}_Or47bmut_GH_R3 ${antibody}_Or47bmut_SH_R1 ${antibody}_Or47bmut_SH_R2 ${antibody}_Or47bmut_SH_R3 ${antibody}_Or67dmut_GH_R1 ${antibody}_Or67dmut_GH_R2 ${antibody}_Or67dmut_GH_R3 ${antibody}_Or67dmut_SH_R1 ${antibody}_Or67dmut_SH_R2 ${antibody}_Or67dmut_SH_R3 ${antibody}_dsxmut_GH_R1 ${antibody}_dsxmut_GH_R2 ${antibody}_dsxmut_GH_R3 ${antibody}_dsxmut_GH_R4 ${antibody}_dsxmut_SH_R1 ${antibody}_dsxmut_SH_R2 ${antibody}_dsxmut_SH_R3 ${antibody}_frumut_GH_R1 ${antibody}_frumut_GH_R2 ${antibody}_frumut_GH_R3 ${antibody}_frumut_SH_R1 ${antibody}_frumut_SH_R2 ${antibody}_frumut_SH_R3"
output_mat="data_output/deeptools_heatmaps/"${antibody}_all_samples_diff_genotype_gh_genes.mat.gz
output_png="data_output/deeptools_heatmaps/"${antibody}_all_samples_diff_genotype_gh_genes.png


# Matrix: 
computeMatrix reference-point --referencePoint center -S $bw_signals -R $bed_regions --skipZeros -b 2000 -a 2000 -o $output_mat

# Heatmap: 
plotHeatmap -m $output_mat \
--averageTypeSummaryPlot mean \
--boxAroundHeatmaps no \
--missingDataColor grey \
--regionsLabel $regions_label \
--samplesLabel $samples_label \
-out $output_png

## ## ## 
## Genotype SH diff genes
## ## ## 

antibody="H3K4me3"
bigwig_path="data_output/03_bed_to_bigwig/"
bw_signals="${bigwig_path}${antibody}_CS_GH_male_brain_R3.bigWig ${bigwig_path}${antibody}_CS_GH_male_brain_R4.bigWig ${bigwig_path}${antibody}_CS_GH_male_brain_R5.bigWig ${bigwig_path}${antibody}_CS_SH_male_brain_R3.bigWig ${bigwig_path}${antibody}_CS_SH_male_brain_R4.bigWig ${bigwig_path}${antibody}_CS_SH_male_brain_R5.bigWig ${bigwig_path}${antibody}_Or47bmut_GH_male_brain_R1.bigWig ${bigwig_path}${antibody}_Or47bmut_GH_male_brain_R2.bigWig ${bigwig_path}${antibody}_Or47bmut_GH_male_brain_R3.bigWig ${bigwig_path}${antibody}_Or47bmut_SH_male_brain_R1.bigWig ${bigwig_path}${antibody}_Or47bmut_SH_male_brain_R2.bigWig ${bigwig_path}${antibody}_Or47bmut_SH_male_brain_R3.bigWig ${bigwig_path}${antibody}_Or67dmut_GH_male_brain_R1.bigWig ${bigwig_path}${antibody}_Or67dmut_GH_male_brain_R2.bigWig ${bigwig_path}${antibody}_Or67dmut_GH_male_brain_R3.bigWig ${bigwig_path}${antibody}_Or67dmut_SH_male_brain_R1.bigWig ${bigwig_path}${antibody}_Or67dmut_SH_male_brain_R2.bigWig ${bigwig_path}${antibody}_Or67dmut_SH_male_brain_R3.bigWig ${bigwig_path}${antibody}_dsxmut_GH_male_brain_R1.bigWig ${bigwig_path}${antibody}_dsxmut_GH_male_brain_R2.bigWig ${bigwig_path}${antibody}_dsxmut_GH_male_brain_R3.bigWig ${bigwig_path}${antibody}_dsxmut_GH_male_brain_R4.bigWig ${bigwig_path}${antibody}_dsxmut_SH_male_brain_R1.bigWig ${bigwig_path}${antibody}_dsxmut_SH_male_brain_R2.bigWig ${bigwig_path}${antibody}_dsxmut_SH_male_brain_R3.bigWig ${bigwig_path}${antibody}_frumut_GH_male_brain_R1.bigWig ${bigwig_path}${antibody}_frumut_GH_male_brain_R2.bigWig ${bigwig_path}${antibody}_frumut_GH_male_brain_R3.bigWig ${bigwig_path}${antibody}_frumut_SH_male_brain_R1.bigWig ${bigwig_path}${antibody}_frumut_SH_male_brain_R2.bigWig ${bigwig_path}${antibody}_frumut_SH_male_brain_R3.bigWig"
bed_regions="data_output/2024-1_diffexpr_output/DSX_SH_genotype_down.bed data_output/2024-1_diffexpr_output/DSX_SH_genotype_up.bed data_output/2024-1_diffexpr_output/FRU_SH_genotype_down.bed data_output/2024-1_diffexpr_output/FRU_SH_genotype_up.bed data_output/2024-1_diffexpr_output/Or47b_SH_genotype_down.bed data_output/2024-1_diffexpr_output/Or47b_SH_genotype_up.bed data_output/2024-1_diffexpr_output/Or67d_SH_genotype_down.bed data_output/2024-1_diffexpr_output/Or67d_SH_genotype_up.bed"
regions_label="DSX_SH_UP DSX_SH_DOWN FRU_SH_UP FRU_SH_DOWN Or47b_SH_UP Or47b_SH_DOWN Or67d_SH_UP Or67d_SH_DOWN"
samples_label="${antibody}_CS_GH_R3 ${antibody}_CS_GH_R4 ${antibody}_CS_GH_R5 ${antibody}_CS_SH_R3 ${antibody}_CS_SH_R4 ${antibody}_CS_SH_R5 ${antibody}_Or47bmut_GH_R1 ${antibody}_Or47bmut_GH_R2 ${antibody}_Or47bmut_GH_R3 ${antibody}_Or47bmut_SH_R1 ${antibody}_Or47bmut_SH_R2 ${antibody}_Or47bmut_SH_R3 ${antibody}_Or67dmut_GH_R1 ${antibody}_Or67dmut_GH_R2 ${antibody}_Or67dmut_GH_R3 ${antibody}_Or67dmut_SH_R1 ${antibody}_Or67dmut_SH_R2 ${antibody}_Or67dmut_SH_R3 ${antibody}_dsxmut_GH_R1 ${antibody}_dsxmut_GH_R2 ${antibody}_dsxmut_GH_R3 ${antibody}_dsxmut_GH_R4 ${antibody}_dsxmut_SH_R1 ${antibody}_dsxmut_SH_R2 ${antibody}_dsxmut_SH_R3 ${antibody}_frumut_GH_R1 ${antibody}_frumut_GH_R2 ${antibody}_frumut_GH_R3 ${antibody}_frumut_SH_R1 ${antibody}_frumut_SH_R2 ${antibody}_frumut_SH_R3"
output_mat="data_output/deeptools_heatmaps/"${antibody}_all_samples_diff_genotype_sh_genes.mat.gz
output_png="data_output/deeptools_heatmaps/"${antibody}_all_samples_diff_genotype_sh_genes.png


# Matrix: 
computeMatrix reference-point --referencePoint center -S $bw_signals -R $bed_regions --skipZeros -b 2000 -a 2000 -o $output_mat

# Heatmap: 
plotHeatmap -m $output_mat \
--averageTypeSummaryPlot mean \
--boxAroundHeatmaps no \
--missingDataColor grey \
--regionsLabel $regions_label \
--samplesLabel $samples_label \
-out $output_png




