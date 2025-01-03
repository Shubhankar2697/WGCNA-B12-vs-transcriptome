## Author: Satyajeet Khare, Asst Prof, Symbiosis International (Deemed University), India. Contact: satyajeetkhare@gmail.com

# set working directory
setwd("/B12_WGCNA/nutshell/")

library(dplyr)

clin_var <- read.csv("clin_var.csv")
bam_files <- read.delim("./bam_files.txt", header = TRUE, sep = "\t")
bam_files_clin_var <- merge(x = clin_var, y = bam_files, by = "child_no", all = F)

library(Rsubread)
# Read assignment 
Counts_fc <- featureCounts(bam_files_clin_var$bam_file, 
                           annot.ext = "GRCh38_p12/gencode.v31.annotation.gtf", 
                           isGTFAnnotationFile = TRUE, GTF.attrType = "gene_name", isPairedEnd=TRUE, nthreads = 30)

library(DESeq2)
# sample names
colnames(Counts_fc$counts) <- bam_files_clin_var$child_no
# Create ColData object
samples<-subset(bam_files_clin_var, select=c("treat_gr"))
rownames(samples) <- bam_files_clin_var[,1]
# Create mycols object
mycols = data.frame(row.names = factor(bam_files_clin_var$child_no))
# Check the column names of count matrix and row names of "mycols" are same
all(mycols %in% colnames(Counts_fc$counts))

# Create DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = Counts_fc$counts, colData = samples, 
                              design = ~ treat_gr)

# Keep genes where there are more than 20 counts
idx <- rowSums(counts(dds)) >= 20
dds <- dds[idx,]

# vst normalization
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
mat <- assay(vsd)
write.table(mat, file = "vst_no_batch_correction.csv", quote = FALSE, sep = ",", row.names = TRUE)

