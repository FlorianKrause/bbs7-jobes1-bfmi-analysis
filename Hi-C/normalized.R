# ==============================================================================
# Hi-C Matrix Normalization and Heatmap Visualization
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
#   This script uses the 'strawr' package (MIT License) to read .hic files.
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
chromosome <- "3"
resolution <- 10000        # 10kb bins
scaling_factor <- 15000000 # Normalization constant

# 3. GENOMIC REGION (Focus: Bbs7 locus)
val_start <- 36250000
val_end   <- 37000000
sBBs      <- 36573142      # Bbs7 start position
eBBs      <- 36613477      # Bbs7 end position

# 4. INPUT FILES (Simplified Sample Names)
# Ensure files are named BFMI1.hic, etc., or update extensions accordingly
bfmi_files <- c("BFMI1.hic", "BFMI2.hic", "BFMI3.hic")
b6n_files  <- c("B6N1.hic", "B6N2.hic", "B6N3.hic")

# ==========================================
# --- END OF CONFIGURATION ---
# ==========================================

# --- 1. DATA LOADING & FILTERING ---
print("Reading .hic files and calculating total library sizes...")

load_and_filter <- function(file, region_vals) {
  if(!file.exists(file)) stop(paste("File not found:", file))
  
  # Load raw counts ("NONE" normalization to calculate library size manually)
  data <- as.matrix(straw("NONE", file, chromosome, chromosome, unit = "BP", resolution))
  total_reads <- sum(data[,3])
  
  # Filter for the region of interest
  filtered <- data[which(data[,1] %in% region_vals & data[,2] %in% region_vals),]
  return(list(data = filtered, total = total_reads))
}

val <- as.character(seq(val_start, val_end, resolution))

# Process all replicates
res_bfmi1 <- load_and_filter(bfmi_files[1], val)
res_bfmi2 <- load_and_filter(bfmi_files[2], val)
res_bfmi3 <- load_and_filter(bfmi_files[3], val)
res_b6n1  <- load_and_filter(b6n_files[1], val)
res_b6n2  <- load_and_filter(b6n_files[2], val)
res_b6n3  <- load_and_filter(b6n_files[3], val)

# --- 2. MATRIX GENERATION & NORMALIZATION ---
print("Generating normalized contact matrices...")

fill_matrix <- function(bed_data, total_reads, region_val) {
  mm <- matrix(0, length(region_val), length(region_val), dimnames = list(region_val, region_val))
  for(x in 1:nrow(bed_data)){
    pos1 <- as.character(bed_data[x,1])
    pos2 <- as.character(bed_data[x,2])
    # Normalize by library size and apply scaling factor
    norm_val <- (bed_data[x,3] / total_reads) * scaling_factor
    mm[pos1, pos2] <- norm_val
    mm[pos2, pos1] <- norm_val
  }
  return(mm)
}

mm.bfmi1 <- fill_matrix(res_bfmi1$data, res_bfmi1$total, val)
mm.bfmi2 <- fill_matrix(res_bfmi2$data, res_bfmi2$total, val)
mm.bfmi3 <- fill_matrix(res_bfmi3$data, res_bfmi3$total, val)
mm.b6n1  <- fill_matrix(res_b6n1$data, res_b6n1$total, val)
mm.b6n2  <- fill_matrix(res_b6n2$data, res_b6n2$total, val)
mm.b6n3  <- fill_matrix(res_b6n3$data, res_b6n3$total, val)

# --- 3. VISUALIZATION ---
print("Plotting heatmaps...")

xrange <- as.numeric(val) / 1000000
yrange <- as.numeric(val) / 1000000
col_palette <- c("gray", "dark gray", "black")
breaks_val  <- c(1, 5, 10, 10000)

# Configure plot layout (3 rows, 2 columns)
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))

plot_hic_matrix <- function(mm, title) {
  image(xrange, yrange, mm, breaks = breaks_val, col = col_palette, 
        main = title, xlab = "Mb", ylab = "Mb")
  # Highlight Bbs7 locus in red
  abline(v = c(sBBs, eBBs) / 1000000, col = "red", lty = 2)
  abline(h = c(sBBs, eBBs) / 1000000, col = "red", lty = 2)
}

# Generate plots for all replicates
plot_hic_matrix(mm.bfmi1, "BFMI Replicate 1")
plot_hic_matrix(mm.b6n1,  "B6N Replicate 1")
plot_hic_matrix(mm.bfmi2, "BFMI Replicate 2")
plot_hic_matrix(mm.b6n2,  "B6N Replicate 2")
plot_hic_matrix(mm.bfmi3, "BFMI Replicate 3")
plot_hic_matrix(mm.b6n3,  "B6N Replicate 3")

print("Analysis finished successfully.")