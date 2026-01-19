# ==============================================================================
# RNA-Seq Read Counting and Normalization Script (RPKM/Log2)
# 
# Copyright (C) 2021 Danny Arends
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see https://www.gnu.org/licenses/.
#
# Authors: 
#   Original code by Danny Arends; 2021 
#   Modified version by Florian Krause (2022-2025)
#
# Dependencies:
#   This script uses Bioconductor packages licensed under their respective terms:
#     GenomicAlignments (GPL-3), GenomicFeatures (GPL-3), Rsamtools (GPL-3),
#     preprocessCore (GPL-2 or later, GPL-3 compatible), BiocParallel (GPL-3),
#     biomaRt (GPL-3)
#
# Note: All package licenses are specified in their respective documentation.
# ==============================================================================

library("GenomicAlignments")
library("GenomicFeatures")
library("Rsamtools")
library("preprocessCore")
library("BiocParallel")
library("biomaRt")

# ==========================================
# --- USER CONFIGURATION (PLEASE EDIT) ---
# ==========================================

# 1. DIRECTORIES
# Path where your BAM files are located (Output of the Bash script)
work_dir <- "/path/to/your/project_directory/Processed_BAMs/"

# Path where the GTF reference file is located
ref_dir  <- "/path/to/your/project_directory/References/"

# 2. FILE NAMES
# GTF Filename
gtf_filename <- "Mus_musculus.GRCm38.102.gtf"

# 3. THREADS
# Adjust based on your CPU cores
cores_to_use <- 4

# ==========================================
# --- END OF CONFIGURATION ---
# ==========================================

setwd(work_dir)
register(MulticoreParam(workers = cores_to_use))

# --- 1. PREPARE GENE ANNOTATION ---
print("Loading GTF annotation...")
gtf_file <- file.path(ref_dir, gtf_filename)
txdb <- makeTxDbFromGFF(gtf_file, format="gtf", organism="Mus musculus")
exons_by_gene <- exonsBy(txdb, by="gene")

# --- 2. LOAD BAM FILES ---
# Defining the two specific samples: B6N4 and BFMI4
# These names must match the output of the Bash pipeline
bam_files <- c("B6N4_final.bam", "BFMI4_final.bam")

# Verify files exist
if(!all(file.exists(bam_files))) {
  stop("Error: One or more BAM files are missing in the working directory.")
}

bam_list  <- BamFileList(bam_files, yieldSize=2000000)

# --- 3. COUNT READS (summarizeOverlaps) ---
print("Counting reads (this may take a while)...")
se <- summarizeOverlaps(features=exons_by_gene, reads=bam_list,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE, 
                        fragments=TRUE)

raw_counts <- assay(se)
# Explicitly naming columns for clarity
colnames(raw_counts) <- c("B6N4", "BFMI4")

# Save Raw Counts
write.table(raw_counts, "Raw_Reads_Per_Gene.txt", sep="\t", quote=FALSE)

# --- 4. RPKM CALCULATION & NORMALIZATION ---
print("Calculating RPKM and performing Quantile Normalization...")

# Gene lengths (sum of exons)
gene_lengths <- sum(width(reduce(exons_by_gene)))
common_genes <- intersect(rownames(raw_counts), names(gene_lengths))
raw_counts   <- raw_counts[common_genes,]
gene_lengths <- gene_lengths[common_genes]

# RPKM Calculation
total_mapped <- colSums(raw_counts)
rpkm <- raw_counts 
for(i in 1:ncol(raw_counts)){
  rpkm[,i] <- (10^9 * raw_counts[,i]) / (total_mapped[i] * gene_lengths)
}

# Quantile Normalization
# Note: With only 2 samples, this forces the distributions to be identical.
rpkm_norm <- normalize.quantiles(as.matrix(rpkm))
colnames(rpkm_norm) <- colnames(rpkm)
rownames(rpkm_norm) <- rownames(rpkm)

# --- 5. LOG2 TRANSFORMATION ---
# Adding 0.1 to avoid log(0)
log2_rpkm <- log2(rpkm_norm + 0.1)

# Save Log2 Matrix
write.table(round(log2_rpkm, 2), "RPKM_Log2_Matrix.txt", sep="\t", quote=FALSE)

# --- 6. CALCULATE FOLD CHANGES ---
# Comparison: BFMI4 vs B6N4
fc_bfmi_vs_b6n <- log2_rpkm[,"BFMI4"] - log2_rpkm[,"B6N4"]

# --- 7. ANNOTATION (mm10) ---
print("Retrieving gene info from Ensembl (Archive Nov 2020)...")
# Using the specific archive to match the GRCm38 (mm10) reference used in alignment
ensembl <- useMart(biomart="ensembl", 
                   dataset="mmusculus_gene_ensembl", 
                   host="https://nov2020.archive.ensembl.org")

gene_info <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "description"),
                   filters="ensembl_gene_id",
                   values=rownames(log2_rpkm),
                   mart=ensembl)

# Remove duplicates if any
gene_info <- gene_info[!duplicated(gene_info$ensembl_gene_id),]
rownames(gene_info) <- gene_info$ensembl_gene_id

# --- 8. FINAL OVERVIEW ---
print("Generating final results table...")

final_df <- data.frame(
  EnsemblID = rownames(log2_rpkm),
  GeneName  = gene_info[rownames(log2_rpkm), "external_gene_name"],
  
  # Log2 Values for individual samples
  Log2_B6N4  = round(log2_rpkm[,"B6N4"], 2),
  Log2_BFMI4 = round(log2_rpkm[,"BFMI4"], 2),
  
  # Fold Change (BFMI4 vs B6N4)
  Log2FC_BFMI4_vs_B6N4 = round(fc_bfmi_vs_b6n, 2),
  
  Description = gene_info[rownames(log2_rpkm), "description"]
)

# Sort by Fold Change (Descending magnitude)
final_df <- final_df[order(final_df$Log2FC_BFMI4_vs_B6N4, decreasing=TRUE),]

write.table(final_df, "Final_Results_Overview.txt", sep="\t", row.names=FALSE, quote=FALSE)

print("Analysis finished.")
print("Outputs generated: Raw_Reads_Per_Gene.txt, RPKM_Log2_Matrix.txt, Final_Results_Overview.txt")