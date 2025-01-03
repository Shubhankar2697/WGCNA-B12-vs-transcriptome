#!/bin/bash
folder="/B12_WGCNA/"
mkdir "${folder}/hisat2"

## Author: Satyajeet Khare, Asst Prof, Symbiosis International (Deemed University), India. Contact: satyajeetkhare@gmail.com

# Human genome index
genome="/GRCh38_p12"

# Extract sample names
for i in $(ls ${folder}/fastq_files/*.gz | xargs -n 1 basename | sed 's/\(.*\)_.*/\1/' | sort -u)

do
echo ${i}

# Alignment of samples that have passed the QC (FastqQC)
hisat2 -p 30 -5 15 -3 25 --dta --known-splicesite-infile ${genome}/gencode.v31.annotation_splicesites.txt -x ${genome}/GRCh38.p12.genome -1 ${folder}/fastq_files/${i}_1.fastq.gz -2 ${folder}/fastq_files/${i}_2.fastq.gz | samtools view -@ 30 -bS - | samtools sort -@ 30 - -o ${folder}/hisat2/${i}_sorted.bam

# Index sorted bam file for browser
samtools index -@ 30 ${folder}/hisat2/${i}_sorted.bam

# For count matrix generation
echo "${folder}/hisat2/${i}_sorted.bam">>${folder}/nutshell/bam_files.txt

# Create second column with Child_no
# Create col names for both columns as 'bam_files' and 'child_no'

done
