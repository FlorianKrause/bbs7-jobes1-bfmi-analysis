# ==============================================================================
# Hi-C Significance Heatmap Visualization (ggplot2)
# 
# Copyright (C) 2025 Florian Krause
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
# Author: 
#   Florian Krause; 2025
#
# Dependencies:
#   This script uses CRAN packages:
#     ggplot2 (MIT), dplyr (MIT)
# ==============================================================================

# --- 1. PREPARE ENVIRONMENT ---
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(ggplot2)
library(dplyr)

# ==========================================
# --- USER CONFIGURATION (PLEASE EDIT) ---
# ==========================================

# 1. PLOT REGION (Mb)
plot_min <- 35.5
plot_max <- 38.0

# 2. VISUAL PARAMETERS
highlight_color <- "blue"
heatmap_high    <- "red3"
heatmap_low     <- "grey90"

# 3. SIGNIFICANCE THRESHOLDS (Used for filtering key structures)
# These p-values identify the central bins of the vertical structure
key_p_values <- c(0.04759065, 0.01473298, 0.04936459)

# ==========================================
# --- END OF CONFIGURATION ---
# ==========================================

# --- 2. INPUT DATA DEFINITION ---
# Significant bins identified from differential interaction analysis
significant_bins <- data.frame(
  Pos_X = c(36350000, 35625000, 36500000, 36500000, 36550000, 36600000, 36825000, 37125000, 37275000, 37425000, 37725000),
  Pos_Y = c(36475000, 35825000, 36575000, 36925000, 36900000, 36850000, 37000000, 37450000, 37500000, 37500000, 38000000),
  Reads_B6N = c(21.52526, 11.90252, 30.05698, 10.83007, 21.48913, 38.16311, 20.12248, 7.92363, 15.77218, 27.08969, 14.78992),
  Reads_BFMI = c(16.51918, 6.44643, 20.59095, 6.62364, 12.78576, 19.64166, 13.92431, 3.31917, 8.83577, 17.38521, 11.63829),
  p_Wert = c(0.02489024, 0.01533588, 0.01879832, 0.04759065, 0.01473298, 0.04936459, 0.03345998, 0.02348913, 0.01388127, 0.03138809, 0.04865973)
)

# --- 3. DATA TRANSFORMATION ---
print("Transforming coordinates and p-values...")

plot_data <- significant_bins %>%
  mutate(
    neg_log10_p = -log10(p_Wert),
    Pos_X_Mb = Pos_X / 1e6,
    Pos_Y_Mb = Pos_Y / 1e6
  )

# Extract key structure bins for highlighting
key_structure_plot_data <- plot_data %>%
  filter(p_Wert %in% key_p_values)

# --- 4. VISUALIZATION (ggplot2) ---
print("Generating significance heatmap...")

p_value_heatmap <- ggplot(plot_data, aes(x = Pos_X_Mb, y = Pos_Y_Mb, fill = neg_log10_p)) +
  # Main Heatmap tiles
  geom_tile(color = "white") + 
  
  # Symmetry axis (diagonal)
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  
  # Highlight key structural bins with a border
  geom_tile(data = key_structure_plot_data, color = highlight_color, linewidth = 1, fill = NA) +
  
  # Color scale configuration
  scale_fill_gradient(low = heatmap_low, high = heatmap_high, name = "-log10(p-value)") +
  
  # Labels and Titles
  labs(
    title = "Significance of Differential Interactions (BFMI vs. B6N)",
    subtitle = paste0("Region Chr3: ", plot_min, "-", plot_max, " Mb, Resolution 25 kb"),
    x = "Genomic Position (Mb)",
    y = "Genomic Position (Mb)"
  ) +
  
  # Axis limits and fixed aspect ratio for square plot
  xlim(plot_min, plot_max) +
  ylim(plot_min, plot_max) +
  coord_fixed() +
  
  # Theme and styling
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "right"
  )

# --- 5. OUTPUT ---
print(p_value_heatmap)
print("Plot generation complete.")