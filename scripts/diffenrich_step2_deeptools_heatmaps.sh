conda activate deeptools

## Split by diff exp genes heatmap
antibody="H3K4me3"
bigwig_path="/Volumes/ExtremePro/Emi/social_experience/cutandrun_dec22/cutandrun_nfrun_dec22/03_peak_calling/03_bed_to_bigwig/"
bw_signals="${bigwig_path}${antibody}_CS_GH_male_brain_R1.bigWig ${bigwig_path}${antibody}_CS_GH_male_brain_R2.bigWig ${bigwig_path}${antibody}_CS_SH_male_brain_R1.bigWig"
bed_regions="data_output/bed/res_housing_cs_combat_down.bed data_output/bed/res_housing_cs_combat_nochange.bed data_output/bed/res_housing_cs_combat_up.bed"
regions_label="GH NoChange SH"
samples_label="${antibody}_CS_GH_R1 ${antibody}_CS_GH_R2 ${antibody}_CS_SH_R1"
output_mat="data_output/deeptools_heatmaps/"${antibody}.mat.gz
output_png="data_output/deeptools_heatmaps/"${antibody}.png

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

## Sorted by logFC heatmap

antibody="H3K27ac"
bigwig_path="/Volumes/ExtremePro/Emi/social_experience/cutandrun_dec22/cutandrun_nfrun_dec22/03_peak_calling/03_bed_to_bigwig/"
bw_signals="${bigwig_path}${antibody}_CS_GH_male_brain_R1.bigWig ${bigwig_path}${antibody}_CS_GH_male_brain_R2.bigWig ${bigwig_path}${antibody}_CS_SH_male_brain_R1.bigWig ${bigwig_path}${antibody}_CS_SH_male_brain_R2.bigWig"
bed_regions="data_output/bed/res_housing_cs_combat_sort.bed"
regions_label="Genes_sorted_by_log2fc_top:_GH_bottom:_SH"
samples_label="${antibody}_CS_GH_R1 ${antibody}_CS_GH_R2 ${antibody}_CS_SH_R1 ${antibody}_CS_SH_R2"
output_mat="data_output/deeptools_heatmaps/"${antibody}_sort.mat.gz
output_png="data_output/deeptools_heatmaps/"${antibody}_sort.png

# Matrix: 
computeMatrix reference-point --referencePoint center -S $bw_signals -R $bed_regions --sortRegions no --skipZeros -b 2000 -a 2000 -o $output_mat

# Heatmap: 
plotHeatmap -m $output_mat \
--averageTypeSummaryPlot mean \
--sortRegions no \
--boxAroundHeatmaps no \
--missingDataColor grey \
--regionsLabel $regions_label \
--samplesLabel $samples_label \
-out $output_png

## Split by diff exp genes heatmap, CPM tracks
antibody="H3K4me3"
bigwig_path="/Volumes/ExtremePro/Emi/social_experience/cutandrun_dec22/norm_bigwigs_feb23/"
bw_signals="${bigwig_path}${antibody}_CS_GH_male_brain_R1.target.dedup.sorted.fragments_igg_subtract.bigWig ${bigwig_path}${antibody}_CS_GH_male_brain_R2.target.dedup.sorted.fragments_igg_subtract.bigWig ${bigwig_path}${antibody}_CS_SH_male_brain_R1.target.dedup.sorted.fragments_igg_subtract.bigWig"
bed_regions="data_output/bed/res_housing_cs_combat_down.bed data_output/bed/res_housing_cs_combat_nochange.bed data_output/bed/res_housing_cs_combat_up.bed"
regions_label="GH NoChange SH"
samples_label="${antibody}_CS_GH_R1 ${antibody}_CS_GH_R2 ${antibody}_CS_SH_R1"
output_mat="data_output/deeptools_heatmaps/"${antibody}_iggsubtract.mat.gz
output_png="data_output/deeptools_heatmaps/"${antibody}_iggsubtract.png

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



