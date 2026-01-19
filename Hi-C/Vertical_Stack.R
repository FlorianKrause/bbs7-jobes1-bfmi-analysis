# ==============================================================================
# Hi-C Locus-Specific Heatmap Comparison (Vertical Stack)
# 
# Copyright (C) 2021 Danny Arends
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Authors: 
#   Original code by Danny Arends; 2021 
#   Modified version by Florian Krause (2025-2026)
#
# License: GNU General Public License v3.0 (GPL-3.0)
#
# Dependencies:
#   This script uses 'strawr' (MIT License) to extract Hi-C data.
# ==============================================================================

library(strawr)

# ==========================================
# --- USER CONFIGURATION (PLEASE EDIT) ---
# ==========================================

# 1. DIRECTORIES
# Set to the directory containing your .hic files
project_dir <- "/path/to/your/hic_data/"
setwd(project_dir)

# 2. ANALYSIS PARAMETERS
chromosome     <- "3"
resolution     <- 10000        # 10kb bins
scaling_factor <- 15000000     # Normalization constant (15M)

# 3. GENOMIC REGION (Focus: Bbs7 locus)
val_start <- 36250000
val_end   <- 37000000
sBBs      <- 36573142          # Bbs7 start position
eBBs      <- 36613477          # Bbs7 end position

# 4. INPUT FILES (Simplified Sample Names)
# This script compares three specific samples in a vertical layout
files <- c("BFMI1.hic", "B6N1.hic", "BFMI2.hic")

# ==========================================
# --- END OF CONFIGURATION ---
# ==========================================

# --- 1. DATA LOADING & FILTERING ---
print("Reading .hic files and calculating total library sizes...")

load_and_filter <- function(file, region_vals) {
  if(!file.exists(file)) stop(paste("File not found:", file))
  
  # Load raw counts (NONE normalization for manual scaling)
  data <- as.matrix(straw("NONE", file, chromosome, chromosome, unit = "BP", resolution))
  total_reads <- sum(data[,3])
  
  # Filter for the region of interest
  filtered <- data[which(data[,1] %in% region_vals & data[,2] %in% region_vals),]
  return(list(data = filtered, total = total_reads))
}

val <- as.character(seq(val_start, val_end, resolution))

res1 <- load_and_filter(files[1], val)
res2 <- load_and_filter(files[2], val)
res3 <- load_and_filter(files[3], val)

# --- 2. MATRIX GENERATION & NORMALIZATION ---
print("Generating normalized contact matrices...")

fill_matrix <- function(bed_data, total_reads, region_val) {
  mm <- matrix(0, length(region_val), length(region_val), dimnames = list(region_val, region_val))
  for(x in 1:nrow(bed_data)){
    pos1 <- as.character(bed_data[x,1])
    pos2 <- as.character(bed_data[x,2])
    # Normalization: (count / total_reads) * scaling_factor
    norm_val <- (bed_data[x,3] / total_reads) * scaling_factor
    mm[pos1, pos2] <- norm_val
    mm[pos2, pos1] <- norm_val
  }
  return(mm)
}

mm1 <- fill_matrix(res1$data, res1$total, val)
mm2 <- fill_matrix(res2$data, res2$total, val)
mm3 <- fill_matrix(res3$data, res3$total, val)

# --- 3. VISUALIZATION (Vertical Stack) ---
print("Plotting heatmaps...")

xrange <- as.numeric(val) / 1000000
yrange <- as.numeric(val) / 1000000
col_palette <- c("gray", "dark gray", "black")
breaks_val  <- c(1, 5, 10, 10000)

# Layout: 3 rows, 1 column
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))

plot_locus_heatmap <- function(mm, title) {
  image(xrange, yrange, mm, breaks = breaks_val, col = col_palette, 
        main = title, xlab = "Position (Mb)", ylab = "Position (Mb)")
  # Vertical lines highlighting the gene boundaries
  abline(v = c(sBBs, eBBs) / 1000000, col = "red", lty = 2, lwd = 1.5)
}

plot_locus_heatmap(mm1, paste("Sample:", files[1]))
plot_locus_heatmap(mm2, paste("Sample:", files[2]))
plot_locus_heatmap(mm3, paste("Sample:", files[3]))

print("Analysis finished successfully.")