# ==============================================================================
# Hi-C Interaction Analysis and P-Value Visualization
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
# Authors: 
#   Original code by Danny Arends; 2021 
#   Modified version by Florian Krause (2021-2025)
#
# Dependencies:
#   This script uses Bioconductor and CRAN packages:
#     strawr (MIT), RColorBrewer (Apache 2.0), biomaRt (GPL-3),
#     GenomicRanges (Artistic-2.0), AnnotationHub (Artistic-2.0), CTCF
# ==============================================================================

library(strawr)
library(RColorBrewer)
library(biomaRt)
library(GenomicRanges)
library(AnnotationHub)
library(CTCF)

# ==========================================
# --- USER CONFIGURATION (PLEASE EDIT) ---
# ==========================================

# 1. DIRECTORIES
# Path to your Hi-C valid pairs and .hic files
setwd("/path/to/your/project_directory/Hi-C_valid_pairs")

# 2. ANALYSIS PARAMETERS
chromosome <- "3"
bin_size   <- 25000  # Resolution: 25.000, 50.000, or 100.000 for p-values
rd2        <- bin_size / 2

# Range for data extraction
val_start  <- 30000000
val_end    <- 40000000

# Range for plotting (Focus area)
pS <- 37500000
pE <- 38000000

# 3. FILE NAMES (.hic files)
# These files should be located in the directory set above
files_bfmi <- c("FKR02_1bfmi_L001.allValidPairs.hic", 
                "FKR02_3bfmi_L001.allValidPairs.hic", 
                "FKR02_5bfmi_L001.allValidPairs.hic")

files_b6n  <- c("FKR02_2b6n_L001.allValidPairs.hic", 
                "FKR02_4b6n_L001.allValidPairs.hic", 
                "FKR02_6b6n_L001.allValidPairs.hic")

# ==========================================
# --- END OF CONFIGURATION ---
# ==========================================

# --- 1. GENE ANNOTATION & CTCF DATA ---
print("Fetching gene annotation from Ensembl...")
bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl", host="https://nov2020.archive.ensembl.org")
res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", 
                                    "end_position", "strand", "external_gene_name", 
                                    "mgi_id", "mgi_symbol", "mgi_description"), 
                     filters = c("chromosomal_region", "biotype"), 
                     values = list("3:35000000:40000000", "protein_coding"), mart = bio.mart)

print("Loading CTCF binding sites...")
ah <- AnnotationHub()
query_data <- query(ah, "mm10.CTCF")
CTCF_mm10 <- query_data[["AH104756"]]
CTCF_mm10 <- CTCF_mm10[CTCF_mm10$q.value < 1]
ii <- ranges(CTCF_mm10[which(seqnames(CTCF_mm10) == "chr3")])

# --- 2. DATA LOADING (STRAW) ---
print("Extracting Hi-C contact matrices...")
load_hic <- function(file) {
  as.matrix(straw("KR", file, chromosome, chromosome, unit = "BP", bin_size))
}

bed.bfmi1 <- load_hic(files_bfmi[1])
bed.bfmi2 <- load_hic(files_bfmi[2])
bed.bfmi3 <- load_hic(files_bfmi[3])
bed.b6n1  <- load_hic(files_b6n[1])
bed.b6n2  <- load_hic(files_b6n[2])
bed.b6n3  <- load_hic(files_b6n[3])

val <- as.character(seq(val_start, val_end, bin_size))

# Filter function for defined range
filter_bed <- function(bed) {
  bed[which(bed[,1] %in% val & bed[,2] %in% val & !is.na(bed[,3])),]
}

bed.bfmi1 <- filter_bed(bed.bfmi1)
bed.bfmi2 <- filter_bed(bed.bfmi2)
bed.bfmi3 <- filter_bed(bed.bfmi3)
bed.b6n1  <- filter_bed(bed.b6n1)
bed.b6n2  <- filter_bed(bed.b6n2)
bed.b6n3  <- filter_bed(bed.b6n3)

# --- 3. CALCULATE DIFFERENTIAL CONTACTS ---
print("Calculating p-values for interactions...")

bed.b6n <- c()
bed.bfmi <- c()
pvals <- c()

for(i in seq(pS, pE, bin_size)){
  for(j in seq(pS, pE, bin_size)){
    # B6N Replicates
    b6n1 <- bed.b6n1[bed.b6n1[,1] == i & bed.b6n1[,2] == j, 3]
    b6n2 <- bed.b6n2[bed.b6n2[,1] == i & bed.b6n2[,2] == j, 3]
    b6n3 <- bed.b6n3[bed.b6n3[,1] == i & bed.b6n3[,2] == j, 3]
    bed.b6n <- rbind(bed.b6n, as.numeric(c(i, j, sum(b6n1, b6n2, b6n3, na.rm=TRUE))))

    # BFMI Replicates
    bfmi1 <- bed.bfmi1[bed.bfmi1[,1] == i & bed.bfmi1[,2] == j, 3]
    bfmi2 <- bed.bfmi2[bed.bfmi2[,1] == i & bed.bfmi2[,2] == j, 3]
    bfmi3 <- bed.bfmi3[bed.bfmi3[,1] == i & bed.bfmi3[,2] == j, 3]
    bed.bfmi <- rbind(bed.bfmi, as.numeric(c(i, j, sum(bfmi1, bfmi2, bfmi3, na.rm=TRUE))))

    # Statistics
    p <- tryCatch(t.test(c(b6n1, b6n2, b6n3), c(bfmi1, bfmi2, bfmi3))$p.value, 
                  error = function(e){return(1)})
    
    pvals <- rbind(pvals, as.numeric(c(i, j, p, mean(c(b6n1, b6n2, b6n3), na.rm=TRUE), 
                                      mean(c(bfmi1, bfmi2, bfmi3), na.rm=TRUE))))
  }
}

colnames(pvals) <- c("X", "Y", "P", "mean(B6)", "mean(BFMI)")

# --- 4. PLOTTING FUNCTIONS ---
colz <- colorRampPalette(c("white", "#fde0dd", "#fa9fb5", "#c51b8a"))(30)

addGene <- function(name, col = "black"){
  sP <- res.biomart[which(res.biomart[, "mgi_symbol"] == name), "start_position"]
  eP <- res.biomart[which(res.biomart[, "mgi_symbol"] == name), "end_position"]
  lines(y = c(sP, eP), x = c(sP, eP), col = "white", lwd=6)
  lines(y = c(sP, eP), x = c(sP, eP), col = col, lwd=3)
}

createPlot <- function(bed.file, main = "Sample"){
  plot(c(pS, pE), c(pS, pE), t = 'n', xaxs='i', yaxs='i', xaxt='n', yaxt='n', 
       xlab="Position (Mb)", ylab="Position (Mb)", main = main)

  for(x in 1:nrow(bed.file)){
    v <- round(-log10(pvals[x,3]) * 10, 0)
    if(v > length(colz)) v <- length(colz)
    if(v > 0){
      rect(bed.file[x,1] - rd2, bed.file[x,2] - rd2, bed.file[x,1] + rd2, bed.file[x,2] + rd2, col = colz[v], border = NA)
      rect(bed.file[x,2] - rd2, bed.file[x,1] - rd2, bed.file[x,2] + rd2, bed.file[x,1] + rd2, col = colz[v], border = NA)
    }
  }

  abline(a = 0, b = 1)

  # Annotate Specific Genes
  addGene("Bbs7", "hotpink")
  addGene("Exosc9", "purple")
  addGene("Ccna2", "chartreuse")
  addGene("Trpc3", "darkgoldenrod1")
  addGene("4932438A13Rik", "burlywood")

  # Highlight CTCF / Bbs7 region
  bbs7CTCF <- (36599801+36600200)/2
  points(bbs7CTCF, 36854135, cex=4, col = "red")
  points(bbs7CTCF, bbs7CTCF, pch=19, cex=1.5, col = "red")
  abline(v = bbs7CTCF, col="red", lwd=2)
  
  axis(1, at = seq(pS, pE, (pE - pS) / 5), round(seq(pS, pE, (pE - pS) / 5) / 1000000, 1), cex.axis=1.5)
  axis(2, at = seq(pS, pE, (pE - pS) / 5), round(seq(pS, pE, (pE - pS) / 5) / 1000000, 1), las=2, cex.axis=1.5)
  box()
}

# --- 5. EXECUTION ---
print("Generating final plot...")
createPlot(bed.b6n, "Interaction P-values (B6N vs BFMI)")

print("Analysis finished.")