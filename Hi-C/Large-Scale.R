# ==============================================================================
# Hi-C Large-Scale Interaction Heatmap (100kb Resolution)
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
# ==============================================================================

library(strawr)

# ==========================================
# --- USER CONFIGURATION (PLEASE EDIT) ---
# ==========================================

# 1. DIRECTORIES
project_dir <- "/path/to/your/hic_data/"
setwd(project_dir)

# 2. ANALYSIS PARAMETERS
chromosome <- "3"
resolution <- 100000  # 100kb resolution
range_start <- 30000000
range_end   <- 40000000

# 3. SAMPLES AND FILENAMES
# Mapping: Directory/File -> Sample Label
samples <- list(
  "BFMI1" = "1BFMI/hic/FKR02_1bfmi_L001.allValidPairs.hic",
  "B6N1"  = "2B6N/hic/FKR02_2b6n_L001.allValidPairs.hic",
  "BFMI2" = "3BFMI/hic/FKR02_3bfmi_L001.allValidPairs.hic",
  "B6N2"  = "4B6N/hic/FKR02_4b6n_L001.allValidPairs.hic",
  "BFMI3" = "5BFMI/hic/FKR02_5bfmi_L001.allValidPairs.hic",
  "B6N3"  = "6B6N/hic/FKR02_6b6n_L001.allValidPairs.hic"
)

# ==========================================
# --- END OF CONFIGURATION ---
# ==========================================

# Preparation of the index range for plotting
iix <- as.character(seq(range_start, range_end, resolution))

# Set up visual layout (e.g., 2 columns, 3 rows)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))

# --- MAIN PROCESSING LOOP ---
for (label in names(samples)) {
  file_path <- samples[[label]]
  
  if(!file.exists(file_path)) {
    warning(paste("File missing, skipping:", file_path))
    next
  }

  print(paste("Processing Sample:", label))
  
  # 1. Load Data
  bed <- as.matrix(straw("NONE", file_path, chromosome, chromosome, unit = "BP", resolution))
  
  # 2. Build Matrix
  u_pos <- unique(c(bed[,1], bed[,2]))
  mm <- matrix(0, length(u_pos), length(u_pos), dimnames = list(u_pos, u_pos))
  
  for(x in 1:nrow(bed)){
    p1 <- as.character(bed[x,1])
    p2 <- as.character(bed[x,2])
    mm[p1, p2] <- bed[x,3]
    mm[p2, p1] <- bed[x,3]
  }
  
  # 3. Visualization
  # Subsetting the matrix to the focus range
  # Using discrete breaks to emphasize structure: 0-100 (white), 100-1000 (gray), 1000+ (black)
  image(as.numeric(iix)/1e6, as.numeric(iix)/1e6, mm[iix, iix], 
        breaks = c(0, 100, 1000, 1e6), 
        col = c("white", "gray", "black"),
        main = label, xlab = "Position (Mb)", ylab = "Position (Mb)")
}

print("Large-scale Hi-C inspection finished.")