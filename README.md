# Project Description
Maternal deficiency of vitamin B12 (B12) is associated with neural tube defects (NTDs), fetal growth restriction, and future risk of non-communicable disease in the offspring. Little is known about the molecular basis of these associations. We investigated the association of vitamin B12 levels and other micronutrient levels in the cord-blood with gene expression in the cord blood mononuclear cells.
A Weighted Gene Co-expression Network Analysis (WGCNA) was performed on cord blood transcriptome of babies born in a pre-conceptional trial of vitamin B12 and multi-micronutrients (MMN). The modules that showed strong correlation with cord blood vitamin B12 levels were subjected for gene enrichment analysis to understand their function.
# Requirements for Analysis 
 1. Load the R script in R studio (v4.1.2 suitable for WGCNA) and download the  transcriptomic data from NCBI SRA (PRJNA756634).
 2. Before Running the script, install packages provided in the script for the analysis.
# Script description
 1. For "1_Alignment.sh" - the reference genome used for alignemnt - GRCh38_p12 and the tool used was HISAT2.
 2.  For " 2_Assignment.R" - the count matrix was generated using featureCounts function of Rsubread package, the count matrix ws normalised using variancestabilizingtransformation function of DESeq2 package in R.
 3.  For "3_WGCNA_CMC_transcriptomic vs Micronutrients.R" - Weighted Gene Co-expression Network Analysis was performed using WGCNA as the primary package in R.
 4.  For "4_Correlogram_Clin.R"- the corelogram was constructed using metan package in R.
* Please Run the script according to the chronology that is provided
