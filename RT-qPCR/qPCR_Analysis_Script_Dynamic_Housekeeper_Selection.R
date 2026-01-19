# ==============================================================================
# qPCR Analysis Script with Dynamic Housekeeper Selection
# 
# Version: 2.2 (Fixed Y-Axis at 12)
#
# This script performs relative gene expression analysis from qPCR data.
# It features an automated stability check to select the most reliable 
# housekeeper genes for normalization and generates comparative barplots.
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
#   Modified version by Florian Krause (2021-2025)
#
# Dependencies:
#   Base R (stats, graphics)
# ==============================================================================

# ==========================================
# --- USER CONFIGURATION (PLEASE EDIT) ---
# ==========================================

# 1. DIRECTORIES
# Path to your qPCR raw data files (B6N and BFMI filtered text files)
project_dir <- "/path/to/your/project_directory/qPCR_data/"
setwd(project_dir)

# 2. ANALYSIS TARGETS
# Define your candidate housekeepers, genes of interest, and tissue types
candidate_housekeepers <- c("GAPDH", "Rpl13a", "rps29")
target_genes           <- c("CX3RC1", "Aldh1l1", "Tubb3", "MBP", "BBS7")
tissue_types           <- c("AC", "MG", "FT")

# 3. FILENAMES
# Ensure these files are present in your project_dir
file_list <- list(
  "B6N 1"  = "B6N_1_HiC_filtered.txt", "B6N 2"  = "B6N_2_HiC_filtered.txt", "B6N 3"  = "B6N_3_HiC_filtered.txt",
  "BFMI 1" = "BFMI_1_HiC_filtered.txt", "BFMI 2" = "BFMI_2_HiC_filtered.txt", "BFMI 3" = "BFMI_3_HiC_filtered.txt"
)

# ==========================================
# --- END OF CONFIGURATION ---
# ==========================================

# --- 1. DATA LOADING ---
load_data <- function(files) {
  lapply(files, function(f) {
    read.table(f, sep = "\t", header = TRUE, fileEncoding = "UTF-8", check.names = FALSE)
  })
}
datasets_list <- load_data(file_list)

# --- 2. HOUSEKEEPER SELECTION LOGIC ---
# Selects stable housekeepers based on the lowest Standard Deviation
selectStableHousekeepers <- function(dataset, candidates, num_to_select = 2) {
  stability_scores <- data.frame(Housekeeper = character(), StDev = numeric())
  for (hk in candidates) {
    hk_cts <- as.numeric(as.character(dataset[dataset$`Target Name` == hk, "Cт Mean"]))
    hk_cts <- hk_cts[is.finite(hk_cts)]
    if (length(hk_cts) > 1) {
      stability_scores <- rbind(stability_scores, data.frame(Housekeeper = hk, StDev = sd(hk_cts, na.rm = TRUE)))
    }
  }
  stable_hks <- stability_scores[order(stability_scores$StDev), ]
  return(as.character(stable_hks$Housekeeper[1:num_to_select]))
}

# --- 3. ANALYSIS FUNCTIONS ---
# Performs dCt calculation
doAnalysis <- function(dataset, selected_housekeepers){
  res <- c()
  for(gene in target_genes){
    for(housekeeper in selected_housekeepers){
      for(tissue in tissue_types){
        iiG <- which(dataset[, "Sample Name"] == tissue & dataset[, "Target Name"] == gene)
        iiH <- which(dataset[, "Sample Name"] == tissue & dataset[, "Target Name"] == housekeeper)
        if(length(iiG) > 0 && length(iiH) > 0) {
            dCT <- dataset[iiG[1], "Cт Mean"] - dataset[iiH[1], "Cт Mean"]
            res <- rbind(res, c(gene, housekeeper, tissue, dCT))
        }
      }
    }
  }
  return(res)
}

# Normalizes dCt values for plotting
calculatePlotValues <- function(analysis_data, dataset_name) {
  results_df <- data.frame(Quelle = character(), Gen = character(), Gewebe = character(), 
                           PlottedValue = numeric(), Standardabweichung_dCt = numeric())
  for (gene in target_genes) {
    # Extract tissue specific values
    extract_vals <- function(t) as.numeric(analysis_data[analysis_data[, 1] == gene & analysis_data[, 3] == t, 4])
    valsAC <- extract_vals("AC"); valsMG <- extract_vals("MG"); valsFT <- extract_vals("FT")
    
    # Correction Factor calculation
    means <- c(mean(valsAC, na.rm=T), mean(valsMG, na.rm=T), mean(valsFT, na.rm=T))
    cF <- -min(-means, na.rm = TRUE) + 1
    
    results_df <- rbind(results_df,
      data.frame(Quelle = dataset_name, Gen = gene, Gewebe = "AC", PlottedValue = cF - mean(valsAC, na.rm=T), Standardabweichung_dCt = sd(valsAC, na.rm=T)),
      data.frame(Quelle = dataset_name, Gen = gene, Gewebe = "MG", PlottedValue = cF - mean(valsMG, na.rm=T), Standardabweichung_dCt = sd(valsMG, na.rm=T)),
      data.frame(Quelle = dataset_name, Gen = gene, Gewebe = "FT", PlottedValue = cF - mean(valsFT, na.rm=T), Standardabweichung_dCt = sd(valsFT, na.rm=T))
    )
  }
  return(results_df)
}

# --- 4. EXECUTION ---
all_final_results <- list()
cat("--- Starting Dynamic Housekeeper Analysis ---\n")

for(dataset_name in names(datasets_list)) {
  current_dataset <- datasets_list[[dataset_name]]
  stable_housekeepers <- selectStableHousekeepers(current_dataset, candidate_housekeepers, num_to_select = 2)
  cat(paste0("Dataset '", dataset_name, "': Selected HKs: ", paste(stable_housekeepers, collapse = ", "), "\n"))
  
  analysis_results <- doAnalysis(current_dataset, stable_housekeepers)
  all_final_results[[dataset_name]] <- calculatePlotValues(analysis_results, dataset_name)
}
final_table <- do.call(rbind, all_final_results)

# --- 5. VISUALIZATION (Aldh1l1 Example) ---
plotExpression <- function(final_data, gene_name) {
  plot_data <- final_data[final_data$Gen == gene_name, ]
  if (nrow(plot_data) == 0) return(print("No data found for plotting."))

  samples <- unique(plot_data$Quelle)
  mat_values <- matrix(NA, nrow = 3, ncol = length(samples), dimnames = list(tissue_types, samples))
  mat_sds    <- mat_values
  
  for (s in samples) {
    for (t in tissue_types) {
      idx <- plot_data$Quelle == s & plot_data$Gewebe == t
      if (any(idx)) {
        mat_values[t, s] <- plot_data$PlottedValue[idx]
        mat_sds[t, s]    <- plot_data$Standardabweichung_dCt[idx]
      }
    }
  }
  
  bp <- barplot(height = mat_values, beside = TRUE, col = c("red", "orange", "blue"),
                main = paste("Relative Expression of", gene_name), ylim = c(0, 12),
                xlab = "Sample", ylab = "Relative Expression (normalized)", legend.text = TRUE)
  arrows(x0 = bp, y0 = mat_values, x1 = bp, y1 = mat_values + mat_sds, angle = 90, code = 2, length = 0.04)
}

par(mfrow = c(1, 1))
plotExpression(final_table, "Aldh1l1")

# --- 6. EXPORT ---
write.csv(final_table, "relative_expression_values_dynamic_selection_v2.csv", row.names = FALSE)
cat("Analysis finished. Results saved to 'relative_expression_values_dynamic_selection_v2.csv'.\n")